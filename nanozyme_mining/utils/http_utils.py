"""
PR1-3 (v4 audit): shared HTTP + filesystem helpers for external data fetchers.

Centralizes:
- HTTP session with a polite User-Agent (RCSB / UniProt expect identification)
- automatic retry with exponential backoff on 5xx / connection errors
- minimum interval throttle so a tight loop can't DOS the upstream
- atomic file writes (tempfile + os.replace) so SIGINT / disk-full doesn't
  leave a half-written file behind

Used by scripts/build_experimental_pdb_library.py and
nanozyme_mining/database/uniprot_fetcher.py.
"""
from __future__ import annotations

import contextlib
import json
import os
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Optional

import requests
from requests.adapters import HTTPAdapter

try:
    # urllib3 >= 1.26
    from urllib3.util.retry import Retry  # type: ignore
except Exception:  # pragma: no cover
    from urllib3.util import Retry  # type: ignore


DEFAULT_USER_AGENT = (
    "E2N/0.2.0 (nanozyme-mining; "
    "https://github.com/taxuannga877-jpg/enzyme-to-nanozyme) "
    "python-requests"
)


def make_session(
    *,
    user_agent: Optional[str] = None,
    total_retries: int = 3,
    backoff_factor: float = 0.6,
    pool_size: int = 10,
) -> requests.Session:
    """
    Build a polite, retrying requests.Session.

    - User-Agent identifies the caller (RCSB/UniProt log this and may throttle
      anonymous bots first).
    - Retry covers 408/425/429/500/502/503/504 with exponential backoff.
    - HEAD/GET only (we never write to these APIs).
    """
    s = requests.Session()
    s.headers["User-Agent"] = user_agent or DEFAULT_USER_AGENT
    retry = Retry(
        total=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=(408, 425, 429, 500, 502, 503, 504),
        allowed_methods=("GET", "HEAD"),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=pool_size,
        pool_maxsize=pool_size,
    )
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


class RateLimiter:
    """
    Minimum-interval throttle (token-less, monotonic-time based).

    Use one instance per host. `wait()` blocks until at least `min_interval`
    seconds have passed since the previous call.

    Thread-safe (Flask debug server is threaded).
    """

    __slots__ = ("min_interval", "_last", "_lock")

    def __init__(self, min_interval: float = 0.5):
        self.min_interval = float(min_interval)
        self._last = 0.0
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            delta = now - self._last
            if delta < self.min_interval:
                time.sleep(self.min_interval - delta)
            self._last = time.monotonic()


def atomic_write_text(path: Path | str, content: str, encoding: str = "utf-8") -> None:
    """
    Write `content` to `path` atomically.

    Pattern: write to a temp file in the same directory, then os.replace().
    SIGINT or disk-full during the write leaves the temp file behind but never
    a half-written target. The final replace is atomic on POSIX and Win32.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=f".{p.name}.",
        suffix=".tmp",
        dir=str(p.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding=encoding) as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, p)
    except Exception:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def atomic_write_json(path: Path | str, data: Any, *, indent: int = 2) -> None:
    """Convenience: atomic JSON write with stable formatting."""
    atomic_write_text(path,
                       json.dumps(data, ensure_ascii=False, indent=indent, default=str))
