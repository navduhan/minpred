"""Download and validate MINpred TFLite model files.

Author: Naveen Duhan
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import urllib.request
from importlib.resources import files
from pathlib import Path

from platformdirs import user_cache_dir


def load_manifest() -> dict:
    resource = files("minpred").joinpath("model_manifest.json")
    return json.loads(resource.read_text(encoding="utf-8"))


def model_directory() -> Path:
    override = os.environ.get("MINPRED_MODEL_DIR")
    if override:
        return Path(override).expanduser().resolve()
    return Path(user_cache_dir("MINpred", "KAABiL")) / "models"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_valid(path: Path, metadata: dict) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == metadata["bytes"]
        and _sha256(path) == metadata["sha256"]
    )


def _download_url(manifest: dict, filename: str) -> str:
    base_url = os.environ.get("MINPRED_MODEL_BASE_URL")
    if base_url:
        return f"{base_url.rstrip('/')}/{filename}"
    record_id = str(manifest.get("zenodo_record_id", ""))
    if not record_id or record_id == "REPLACE_WITH_RECORD_ID":
        raise RuntimeError(
            "The Zenodo model record has not been configured yet. Set "
            "MINPRED_MODEL_DIR to a directory containing the TFLite files, or "
            "update zenodo_record_id in model_manifest.json after publication."
        )
    return f"https://zenodo.org/records/{record_id}/files/{filename}?download=1"


def get_model(model_key: str, download: bool = True) -> Path:
    manifest = load_manifest()
    try:
        metadata = manifest["models"][model_key]
    except KeyError as exc:
        raise KeyError(f"Unknown MINpred model: {model_key}") from exc

    directory = model_directory()
    destination = directory / metadata["filename"]
    if is_valid(destination, metadata):
        return destination
    if not download:
        raise FileNotFoundError(destination)

    directory.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    url = _download_url(manifest, metadata["filename"])
    print(f"Downloading {metadata['filename']} to {directory}")
    try:
        with urllib.request.urlopen(url) as response, partial.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        if not is_valid(partial, metadata):
            raise RuntimeError(
                f"Integrity check failed for {metadata['filename']}; "
                "the incomplete download was removed."
            )
        partial.replace(destination)
    finally:
        if partial.exists():
            partial.unlink()
    return destination


def select_models(selection: str) -> list[str]:
    models = load_manifest()["models"]
    normalized = selection.lower()
    if normalized == "all":
        return list(models)
    if normalized in {"phase1", "phase2", "phase3"}:
        return [normalized]
    if normalized == "phase4":
        return [key for key in models if key.startswith("phase4_")]
    if normalized in models:
        return [normalized]
    raise ValueError(f"Unknown model selection: {selection}")


def download_models(selection: str = "all") -> list[Path]:
    return [get_model(key) for key in select_models(selection)]


def model_status() -> list[tuple[str, Path, bool]]:
    manifest = load_manifest()
    directory = model_directory()
    return [
        (key, directory / metadata["filename"], is_valid(directory / metadata["filename"], metadata))
        for key, metadata in manifest["models"].items()
    ]
