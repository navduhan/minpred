#!/usr/bin/env bash
#SBATCH --job-name=minpred
#SBATCH --partition=mahaguru
#SBATCH --nodelist=chela-g01
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
umask 077

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 input.fasta Phase1|Phase2|Phase3|Phase4 enzyme-class output-directory protein|nucleotide" >&2
  exit 64
fi

input_fasta=$1
level=$2
enzyme_class=$3
output_dir=$4
sequence_type=$5
app_dir="${MINPRED_APP_DIR:-$HOME/naveen_tools/minpred}"

[[ -f "$input_fasta" && -s "$input_fasta" ]] || { echo "Input FASTA is missing or empty: $input_fasta" >&2; exit 66; }
[[ -f "$app_dir/minpred.py" ]] || { echo "MINpred is missing: $app_dir" >&2; exit 69; }
[[ "$level" =~ ^Phase[1-4]$ ]] || { echo "Invalid prediction level: $level" >&2; exit 64; }
[[ "$sequence_type" == protein || "$sequence_type" == nucleotide ]] || { echo "Invalid sequence type: $sequence_type" >&2; exit 64; }
[[ "$enzyme_class" =~ ^(all|amidohydrolases|aminopeptidases|aspartic|cysteine|dipeptidases|dipeptidyl|metalloendopeptidases|metallopeptidases|omega|serine)$ ]] || { echo "Invalid Phase IV enzyme class: $enzyme_class" >&2; exit 64; }

module_name="${MINPRED_MODULE:-dl-gpu}"
if command -v module >/dev/null 2>&1; then
  module load "$module_name"
elif [[ -x /usr/bin/modulecmd ]]; then
  eval "$(/usr/bin/modulecmd bash load "$module_name")"
else
  echo "The cluster module command is unavailable; cannot load $module_name." >&2
  exit 69
fi

# Keep a fallback for dl-gpu installations whose public TransDecoder symlinks
# do not resolve. The modulefile should normally expose this directory itself.
module_python="$(command -v python)"
module_prefix="${CONDA_PREFIX:-$(dirname "$(dirname "$module_python")")}"
transdecoder_util="$module_prefix/opt/transdecoder/util"
if ! command -v TransDecoder.LongOrfs >/dev/null 2>&1 && [[ -d "$transdecoder_util" ]]; then
  export PATH="$transdecoder_util:$PATH"
fi
if [[ "$sequence_type" == nucleotide ]] && ! command -v TransDecoder.LongOrfs >/dev/null 2>&1; then
  echo "TransDecoder.LongOrfs is unavailable in the $module_name module." >&2
  exit 69
fi

if [[ -x "$app_dir/.venv/bin/python" ]]; then
  python_bin="$app_dir/.venv/bin/python"
else
  python_bin="$(command -v python || true)"
fi
[[ -n "$python_bin" && -x "$python_bin" ]] || { echo "The $module_name Python interpreter was not found." >&2; exit 69; }

mkdir -p "$output_dir"
export MINPRED_MODEL_DIR="${MINPRED_MODEL_DIR:-$app_dir/models}"
export PYTHONPATH="$app_dir${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH="${TF_FORCE_GPU_ALLOW_GROWTH:-true}"
export TF_CPP_MIN_LOG_LEVEL="${TF_CPP_MIN_LOG_LEVEL:-2}"

"$python_bin" -c 'import Bio, numpy, pandas, platformdirs, tensorflow'
"$python_bin" "$app_dir/minpred.py" \
  -i "$input_fasta" \
  -od "$output_dir" \
  -o minpred_predictions.tsv \
  -l "$level" \
  -n "$enzyme_class" \
  --sequence-type "$sequence_type"

for result in "$output_dir"/Phase_*_log.txt; do
  [[ -f "$result" ]] || continue
  cp -- "$result" "${result%.txt}.tsv"
done

prediction_file=$(find "$output_dir" -maxdepth 1 -type f -name '*.tsv' -size +0c -print -quit)
if [[ -z "$prediction_file" ]]; then
  echo "MINpred finished without a non-empty prediction table." >&2
  exit 70
fi

echo "MINpred results: $output_dir"
