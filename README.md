# MINpred

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/usubioinfo/minpred/releases/tag/v1.0.0)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Models: Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.21897103.svg)](https://doi.org/10.5281/zenodo.21897103)

MINpred is an alignment-free, four-phase deep-learning tool for identifying
and classifying nitrogen mineralization-related enzymes. It first distinguishes
enzymes from non-enzymes, identifies nitrogen mineralization enzymes, assigns
them to ten functional classes, and finally predicts class-specific EC numbers.
The standalone command-line interface accepts protein FASTA files directly and
can use TransDecoder to predict proteins from nucleotide or transcript FASTA
files. Lightweight TFLite inference and separately archived, checksum-verified
models keep the source repository small and make model retrieval reproducible.

## Key Features

**Hierarchical four-phase pipeline**

- **Phase I — Enzyme vs. non-enzyme:** binary CNN using the 640-dimensional
  DPC+NMBroto representation (average 10-fold CV accuracy: **95.70%**,
  MCC: **0.914**; independent-test accuracy: **93.43%**, MCC: **0.868**).
- **Phase II — Nitrogen mineralization vs. non-mineralization enzyme:** binary
  CNN using DPC+NMBroto (average 10-fold CV accuracy: **95.23%**, MCC:
  **0.9042**; independent-test accuracy: **92.43%**, MCC: **0.8487**).
- **Phase III — Nitrogen mineralization enzyme class:** multiclass CNN covering
  ten enzyme classes using DPC+NMBroto (average 10-fold CV overall accuracy:
  **96.23%**, MCC: **0.8706**; independent-test overall accuracy: **93.83%**,
  MCC: **0.7640**).
- **Phase IV — EC-number assignment:** ten class-specific CNNs using the
  2,400-dimensional CKSAAP representation, covering **69 specific EC numbers**
  plus wildcard and `Others` outputs. Independent-test overall accuracy across
  the reported classifiers ranges from **88.43% to 98.24%**.

**Modern packaging and portable TFLite inference**

- Managed through **PEP 621** metadata in `pyproject.toml`, with a reproducible
  Python dependency environment in `uv.lock`.
- Uses **TensorFlow Lite (`.tflite`)** models for portable CPU inference; model
  files are downloaded from Zenodo only when required and verified using file
  size and SHA-256 checksums.
- Supports protein FASTA directly and nucleotide/transcript FASTA through the
  external **TransDecoder** dependency.

## Installation

MINpred requires Python 3.9 or newer. Install the package and TensorFlow Lite
runtime support with:

```bash
pip install '.[tensorflow]'
```

For a reproducible development environment using the committed `uv.lock`:

```bash
uv sync --extra tensorflow
uv run minpred models-status
```

The trained models are distributed separately through Zenodo. A required model
is downloaded automatically on first use and validated against its published
file size and SHA-256 checksum. To download models in advance:

```bash
minpred download-models all
minpred models-status
```

By default, models are stored in the operating system's user cache directory.
Set `MINPRED_MODEL_DIR` to use a specific local model directory.

## Prediction

```bash
minpred -i example/test.fasta -l Phase1 -od MINpred_results
minpred -i example/test.fasta -l Phase4 -n amidohydrolases -od MINpred_results
```

MINpred accepts protein FASTA input by default. For transcript or other
nucleotide FASTA input, install TransDecoder and select nucleotide mode:

```bash
conda install -c conda-forge -c bioconda transdecoder
minpred -i transcripts.fasta --sequence-type nucleotide -l Phase4 \
  -n amidohydrolases -od MINpred_results
```

In nucleotide mode, MINpred runs `TransDecoder.LongOrfs`, preserves its log and
the predicted peptide FASTA inside the output directory, removes terminal stop
symbols, writes `translated_proteins.fasta`, and then applies the selected
MINpred phases to those peptides. TransDecoder is an external Conda dependency and is therefore not
included in `uv.lock`, which covers the Python environment only. The supplied
`minpred_environment.yml` installs both TransDecoder and the Python runtime.

Phases I–III use the 640-dimensional DPC+NMBroto hybrid representation. Phase
IV uses the 2,400-dimensional CKSAAP representation.

## Model files

The trained MINpred models were converted to TensorFlow Lite (`.tflite`) for
efficient and portable CPU inference. They are downloaded automatically from
Zenodo when required and verified using the file sizes and SHA-256 checksums
recorded in `minpred/model_manifest.json`.

The model files are archived in Zenodo record
[`21897103`](https://doi.org/10.5281/zenodo.21897103) under the Creative Commons
Attribution 4.0 International license (CC BY 4.0).

## License

The MINpred source code is distributed under the [MIT License](LICENSE).
The separately archived TFLite model files are distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite the
Zenodo model record when using the trained models.

## Citation

Until the associated journal article is published, cite the model release as:

> Duhan, N., Norton, J. M., & Kaundal, R. (2026). *MINpred: TFLite models for
> nitrogen mineralization-related enzyme classification* (Version 1.0.0).
> Zenodo. https://doi.org/10.5281/zenodo.21897103

---

Developed by **Naveen Duhan**
*Kaundal Artificial Intelligence & Advanced Bioinformatics Lab (KAABiL)*
Utah State University, Logan, UT

**Technical Queries / Bugs:** Naveen Duhan
([naveen.duhan@usu.edu](mailto:naveen.duhan@usu.edu))

**Scientific Queries:** Dr. Rakesh Kaundal
([rkaundal@usu.edu](mailto:rkaundal@usu.edu))

**Lab Website:** [https://kaabil.net](https://kaabil.net)

Released under the terms of the **[MIT License](LICENSE)**.
