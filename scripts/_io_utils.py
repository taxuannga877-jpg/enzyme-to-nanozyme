from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any


def esc(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except (TypeError, ValueError):
        pass
    return html.escape(str(value), quote=True)


def as_int(value: Any, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError, OverflowError):
        return default


def as_float(value: Any, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def pct(numer: float, denom: float, digits: int = 1) -> str:
    if not denom:
        return f"{0:.{digits}f}%"
    return f"{100.0 * numer / denom:.{digits}f}%"


def pct_ratio(numer: float, denom: float) -> float:
    return 0.0 if not denom else 100.0 * numer / denom


def slug(value: Any) -> str:
    return "_".join(
        part
        for part in "".join(
            char.lower() if char.isalnum() else " " for char in str(value)
        ).split()
        if part
    )


def json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if hasattr(value, "to_dict"):
        return value.to_dict()
    return str(value)


def json_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_safe) + "\n",
        encoding="utf-8",
    )
    temp.replace(path)


def load_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default
