from __future__ import annotations

import sys
from pathlib import Path


def ensure_repo_imports(root: Path) -> None:
    """Make repo packages and sibling script helpers importable for CLI runners."""
    root = Path(root).resolve()
    for path in (root, root / "scripts"):
        value = str(path)
        if value not in sys.path:
            sys.path.insert(0, value)
