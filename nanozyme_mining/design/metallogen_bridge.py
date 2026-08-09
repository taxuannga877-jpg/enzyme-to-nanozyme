"""Repo-local bridge for optional vendored MetalloGen geometry data."""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Dict


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class MetalloGenConfig:
    """Paths and settings for the optional vendored MetalloGen backend."""

    metallogen_root: Path = PROJECT_ROOT / "linker" / "metallogen"
    working_directory: Path = PROJECT_ROOT / "outputs" / "metallogen_work"
    save_directory: Path = PROJECT_ROOT / "outputs" / "metallogen_out"
    calculator: str = "orca"


def ensure_metallogen_on_path(root: Path) -> bool:
    """Return whether a repo-local MetalloGen checkout is present."""

    root = Path(root)
    return (root / "MetalloGen").is_dir()


def load_known_geometry_vectors(root: Path | None = None) -> Dict[str, object]:
    """Load vendored MetalloGen geometry vectors without mutating sys.path."""

    base = Path(root) if root is not None else MetalloGenConfig().metallogen_root
    globalvars_path = base / "MetalloGen" / "globalvars.py"
    if not globalvars_path.exists():
        return {}

    spec = importlib.util.spec_from_file_location(
        "e2n_vendored_metallogen_globalvars",
        globalvars_path,
    )
    if spec is None or spec.loader is None:
        return {}
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception:
        return {}
    vectors = getattr(module, "known_geometries_vector_dict", {})
    return dict(vectors) if isinstance(vectors, dict) else {}
