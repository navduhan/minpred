"""Translate transcript FASTA input with the external TransDecoder program.

Author: Naveen Duhan
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            packaged = Path(conda_prefix) / "opt" / "transdecoder" / "util" / name
            if packaged.is_file() and os.access(packaged, os.X_OK):
                executable = str(packaged)
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
    run_dir = Path(output_dir).expanduser().resolve() / (
        "transdecoder_" + _file_digest(source)[:12]
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    local_input = run_dir / "transcripts.fasta"
    transdecoder_output = run_dir / "transdecoder_out"
    peptide_fasta = transdecoder_output / f"{local_input.name}.transdecoder_dir" / "longest_orfs.pep"
    legacy_peptide_fasta = transdecoder_output / "longest_orfs.pep"
    stable_peptide_fasta = Path(output_dir).expanduser().resolve() / "translated_proteins.fasta"

    for candidate in (peptide_fasta, legacy_peptide_fasta):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return _write_clean_peptides(candidate, stable_peptide_fasta)

    shutil.copy2(source, local_input)
    log_path = run_dir / "TransDecoder.LongOrfs.log"
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [long_orfs, "-t", local_input.name, "--output_dir", str(transdecoder_output)],
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{Path(long_orfs).name} failed with exit code "
            f"{completed.returncode}. See {log_path}."
        )

    if (not peptide_fasta.is_file() or peptide_fasta.stat().st_size == 0) and (
        not legacy_peptide_fasta.is_file() or legacy_peptide_fasta.stat().st_size == 0
    ):
        raise RuntimeError(
            "TransDecoder completed without producing a non-empty peptide "
            f"FASTA at {peptide_fasta}."
        )
    source_peptides = peptide_fasta if peptide_fasta.is_file() else legacy_peptide_fasta
    return _write_clean_peptides(source_peptides, stable_peptide_fasta)


def _write_clean_peptides(source: Path, destination: Path) -> str:
    """Write TransDecoder peptides to a stable, web-compatible FASTA path."""
    records = list(SeqIO.parse(source, "fasta"))
    if not records:
        raise RuntimeError("TransDecoder found no open reading frames.")
    for record in records:
        record.seq = Seq(str(record.seq).rstrip("*"))
    SeqIO.write(records, destination, "fasta")
    if not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError(f"Unable to write translated peptide FASTA at {destination}.")
    return str(destination)
