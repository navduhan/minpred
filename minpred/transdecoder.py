"""Translate transcript FASTA input with the external TransDecoder program.

Author: Naveen Duhan
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(
            f"{name} was not found on PATH. Install TransDecoder with "
            "`conda install -c conda-forge -c bioconda transdecoder`, then "
            "activate that environment before running nucleotide input."
        )
    return executable


def translate_nucleotide_fasta(input_fasta: str, output_dir: str) -> str:
    """Predict peptide sequences from transcripts and return the peptide FASTA."""
    source = Path(input_fasta).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Nucleotide FASTA not found: {source}")

    long_orfs = _executable("TransDecoder.LongOrfs")
    predict = _executable("TransDecoder.Predict")
    run_dir = Path(output_dir).expanduser().resolve() / (
        "transdecoder_" + _file_digest(source)[:12]
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    local_input = run_dir / "transcripts.fasta"
    peptide_fasta = run_dir / "transcripts.fasta.transdecoder.pep"

    if peptide_fasta.is_file() and peptide_fasta.stat().st_size > 0:
        return str(peptide_fasta)

    shutil.copy2(source, local_input)
    commands = [
        (long_orfs, "TransDecoder.LongOrfs.log"),
        (predict, "TransDecoder.Predict.log"),
    ]
    for executable, log_name in commands:
        log_path = run_dir / log_name
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [executable, "-t", local_input.name],
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if completed.returncode != 0:
            raise RuntimeError(
                f"{Path(executable).name} failed with exit code "
                f"{completed.returncode}. See {log_path}."
            )

    if not peptide_fasta.is_file() or peptide_fasta.stat().st_size == 0:
        raise RuntimeError(
            "TransDecoder completed without producing a non-empty peptide "
            f"FASTA at {peptide_fasta}."
        )
    return str(peptide_fasta)
