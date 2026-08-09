"""Bounded configuration parsing shared by web and compute entry points."""

from __future__ import annotations

import logging
import math
import os
from typing import Mapping


log = logging.getLogger("e2n.config")

_TRUE_VALUES = {"1", "true", "yes", "on"}
_FALSE_VALUES = {"0", "false", "no", "off", ""}


def clamp_int(value, *, default: int, min_value: int, max_value: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError, OverflowError):
        parsed = default
    return max(min_value, min(max_value, parsed))


def clamp_float(
    value,
    *,
    default: float,
    min_value: float,
    max_value: float,
) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        parsed = default
    if not math.isfinite(parsed):
        parsed = default
    return max(min_value, min(max_value, parsed))


def parse_bool(value, *, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            return default
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in _TRUE_VALUES:
        return True
    if normalized in _FALSE_VALUES:
        return False
    return default


def env_int(
    name: str,
    default: int,
    *,
    min_value: int | None = None,
    max_value: int | None = None,
    environ: Mapping[str, str] | None = None,
) -> int:
    source = os.environ if environ is None else environ
    raw = source.get(name)
    if raw is None or raw == "":
        parsed = default
    else:
        try:
            parsed = int(raw)
        except (TypeError, ValueError, OverflowError):
            log.warning(
                "env var %s=%r is not an integer; using default %d",
                name,
                raw,
                default,
            )
            parsed = default

    bounded = parsed
    if min_value is not None:
        bounded = max(min_value, bounded)
    if max_value is not None:
        bounded = min(max_value, bounded)
    if bounded != parsed:
        log.warning(
            "env var %s=%r is outside allowed bounds; using %d",
            name,
            raw,
            bounded,
        )
    return bounded


def env_float(
    name: str,
    default: float,
    *,
    min_value: float | None = None,
    max_value: float | None = None,
    environ: Mapping[str, str] | None = None,
) -> float:
    source = os.environ if environ is None else environ
    raw = source.get(name)
    if raw is None or raw == "":
        parsed = default
    else:
        try:
            parsed = float(raw)
        except (TypeError, ValueError, OverflowError):
            log.warning(
                "env var %s=%r is not a number; using default %s",
                name,
                raw,
                default,
            )
            parsed = default
    if not math.isfinite(parsed):
        log.warning(
            "env var %s=%r is not finite; using default %s",
            name,
            raw,
            default,
        )
        parsed = default

    bounded = parsed
    if min_value is not None:
        bounded = max(min_value, bounded)
    if max_value is not None:
        bounded = min(max_value, bounded)
    if bounded != parsed:
        log.warning(
            "env var %s=%r is outside allowed bounds; using %s",
            name,
            raw,
            bounded,
        )
    return bounded
