"""
Flask security helpers — PR0-2 (v4 audit fixes)

All security-sensitive functions used by enzyme_viewer/app.py live here so the
fixes can be unit-tested in isolation.

Covered fixes from final audit:
- C3: debug/host from environment, not hardcoded
- C4: unified error handler (no traceback to client; rid for log correlation)
- C5: safe_join + filename/format whitelist (path traversal)
- N-C1: pdb_id whitelist (NEW-1 glob injection + accidental file enumeration)
- N-C2: subprocess stdout/stderr never returned to client
- N-H4: CORS limited to local origins by default; configurable via env
"""
from __future__ import annotations

import logging
import hmac
import math
import os
import re
import secrets
import uuid
from functools import wraps
from pathlib import Path
from typing import Iterable, Mapping, Optional, Tuple
from urllib.parse import urlparse

from flask import Flask, Response, jsonify, request
from nanozyme_mining.utils.config import parse_bool

# ---------------------------------------------------------------------------
# Configuration — driven by environment, with safe defaults
# ---------------------------------------------------------------------------

# C3: debug/host from environment (defaults to "off" + loopback)
DEBUG = os.environ.get("FLASK_DEBUG", "0").strip().lower() in {"1", "true", "yes", "on"}
HOST  = os.environ.get("FLASK_HOST", "127.0.0.1").strip() or "127.0.0.1"
PORT  = None  # populated below via env_int


# helper defined first so PORT can use it
def env_int(
    name: str,
    default: int,
    *,
    min_value: int | None = None,
    max_value: int | None = None,
    environ: Mapping[str, str] | None = None,
) -> int:
    """
    PR4-1 (M28 fix): parse an integer environment variable with safe fallback.

    int(os.environ.get(name, default)) raises ValueError on malformed input
    (e.g. FLASK_PORT='8000abc'), which would crash module import — the worst
    place for a config bug to manifest. This helper logs a warning and falls
    back to the documented default instead.
    """
    source = os.environ if environ is None else environ
    raw = source.get(name)
    if raw is None or raw == "":
        parsed = default
    else:
        try:
            parsed = int(raw)
        except (ValueError, TypeError, OverflowError):
            log.warning(
                "env var %s=%r is not an integer; using default %d",
                name,
                raw,
                default,
            )
            parsed = default
    if min_value is not None:
        parsed = max(min_value, parsed)
    if max_value is not None:
        parsed = min(max_value, parsed)
    return parsed


import logging  # late import — already imported up top, this is a no-op but explicit
log = logging.getLogger("e2n.security")
PORT = env_int("FLASK_PORT", 5000, min_value=1, max_value=65535)

# N-H4: CORS origins from environment (comma-separated). Default to localhost only.
def _parse_cors_origins() -> list:
    raw = os.environ.get("FLASK_CORS_ORIGINS", "").strip()
    if not raw:
        return [f"http://127.0.0.1:{PORT}", f"http://localhost:{PORT}"]
    if raw == "*":
        return "*"  # explicit opt-in to wildcard
    return [o.strip() for o in raw.split(",") if o.strip()]

CORS_ORIGINS = _parse_cors_origins()

# C4: error verbosity (off by default, on for dev)
EXPOSE_TRACEBACK = os.environ.get("FLASK_EXPOSE_TRACEBACK", "0").strip() == "1"

DISABLE_CSRF_CHECK = os.environ.get("E2N_DISABLE_CSRF_CHECK", "0").strip().lower() in {
    "1", "true", "yes", "on"
}
AUTH_USERNAME = os.environ.get("E2N_AUTH_USERNAME", "")
AUTH_PASSWORD = os.environ.get("E2N_AUTH_PASSWORD", "")
ALLOW_UNAUTHENTICATED_REMOTE = os.environ.get(
    "E2N_ALLOW_UNAUTHENTICATED_REMOTE",
    "0",
).strip().lower() in {"1", "true", "yes", "on"}

log = logging.getLogger("e2n.security")


def _is_loopback_host(host: str) -> bool:
    normalized = str(host or "").strip().lower().strip("[]")
    return normalized in {"127.0.0.1", "localhost", "::1"}


def install_basic_auth(
    app: Flask,
    *,
    host: str = HOST,
    username: str = AUTH_USERNAME,
    password: str = AUTH_PASSWORD,
    allow_unauthenticated_remote: bool = ALLOW_UNAUTHENTICATED_REMOTE,
) -> bool:
    """Protect remote deployments with optional HTTP Basic authentication."""

    has_username = bool(username)
    has_password = bool(password)
    if has_username != has_password:
        raise RuntimeError(
            "E2N_AUTH_USERNAME and E2N_AUTH_PASSWORD must be configured together"
        )
    if not has_username:
        if not _is_loopback_host(host) and not allow_unauthenticated_remote:
            raise RuntimeError(
                "remote FLASK_HOST requires E2N_AUTH_USERNAME/E2N_AUTH_PASSWORD "
                "or explicit E2N_ALLOW_UNAUTHENTICATED_REMOTE=1"
            )
        return False

    @app.before_request
    def _require_basic_auth():
        auth = request.authorization
        valid = (
            auth is not None
            and hmac.compare_digest(auth.username or "", username)
            and hmac.compare_digest(auth.password or "", password)
        )
        if valid:
            return None
        return Response(
            "Authentication required",
            401,
            {"WWW-Authenticate": 'Basic realm="E2N", charset="UTF-8"'},
        )

    return True


def load_secret_key(app: Flask, key_file: Optional[Path] = None) -> str:
    """Configure a stable Flask SECRET_KEY from env or a local runtime file."""
    secret = os.environ.get("E2N_SECRET_KEY", "").strip()
    if not secret:
        raw_path = os.environ.get("E2N_SECRET_KEY_FILE", "").strip()
        path = Path(raw_path) if raw_path else Path(key_file or Path(app.instance_path) / "secret_key")
        if path.exists():
            secret = path.read_text(encoding="utf-8").strip()
        if not secret:
            path.parent.mkdir(parents=True, exist_ok=True)
            secret = secrets.token_urlsafe(48)
            path.write_text(secret + "\n", encoding="utf-8")
            try:
                os.chmod(path, 0o600)
            except OSError:
                pass
    app.config["SECRET_KEY"] = secret
    app.secret_key = secret
    return secret


def _origin_from_url(value: str) -> Optional[str]:
    parsed = urlparse(value)
    if not parsed.scheme or not parsed.netloc:
        return None
    return f"{parsed.scheme}://{parsed.netloc}"


def _origin_allowed(origin: str) -> bool:
    if CORS_ORIGINS == "*":
        return True
    normalized = _origin_from_url(origin) or origin.rstrip("/")
    allowed = {str(item).rstrip("/") for item in CORS_ORIGINS}
    allowed.add(request.host_url.rstrip("/"))
    return normalized in allowed


def require_json_csrf(view_func):
    """Require a browser AJAX marker and same-origin metadata for write APIs."""
    @wraps(view_func)
    def _wrapped(*args, **kwargs):
        if not DISABLE_CSRF_CHECK:
            if request.headers.get("X-Requested-With") != "XMLHttpRequest":
                return jsonify({"status": "error", "error": "missing X-Requested-With header"}), 403
            origin = request.headers.get("Origin") or request.headers.get("Referer")
            if origin and not _origin_allowed(origin):
                return jsonify({"status": "error", "error": "origin not allowed"}), 403
        return view_func(*args, **kwargs)

    return _wrapped


def clamp_int(value, *, default: int, min_value: int, max_value: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError, OverflowError):
        parsed = default
    return max(min_value, min(max_value, parsed))


def clamp_float(value, *, default: float, min_value: float, max_value: float) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        parsed = default
    if not math.isfinite(parsed):
        parsed = default
    return max(min_value, min(max_value, parsed))


# ---------------------------------------------------------------------------
# Validators — used by API endpoints to sanitize untrusted client input
# ---------------------------------------------------------------------------

# PDB IDs are 4 alphanumeric chars per RCSB spec. Reject anything else early.
_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")

# EC numbers look like "1.11.1.7" — digits + dots, optional 'n' for unassigned subclass.
_EC_NUMBER_RE = re.compile(r"^[0-9n]+(\.[0-9n]+){0,3}$")

# Download format whitelist (matches structure_exporter outputs)
DOWNLOAD_FORMATS = frozenset({"pdb", "xyz", "sdf", "cif"})

# Motif IDs include EC numbers embedded with dots, e.g.
# 5HPW_1.11.1.7_Peroxidase.  Slash/path separators remain disallowed and
# safe_join still enforces containment before any file is opened.
_MOTIF_ID_RE = re.compile(r"^[A-Za-z0-9_.\-]{1,128}$")


def is_valid_pdb_id(pdb_id: str) -> bool:
    """N-C1 / NEW-1: strict 4-char PDB ID whitelist (rejects glob metachars * ? [)."""
    return bool(pdb_id) and bool(_PDB_ID_RE.match(pdb_id))


def is_valid_ec_number(ec: str) -> bool:
    return bool(ec) and bool(_EC_NUMBER_RE.match(ec))


def is_valid_download_format(fmt: str) -> bool:
    return fmt in DOWNLOAD_FORMATS


def is_valid_motif_id(motif_id: str) -> bool:
    return bool(motif_id) and bool(_MOTIF_ID_RE.match(motif_id))


# ---------------------------------------------------------------------------
# Path safety
# ---------------------------------------------------------------------------

def safe_join(base: Path, *parts: str) -> Path:
    """
    C5 fix: resolve `base/*parts` and assert it stays inside `base`.

    Raises ValueError on any traversal attempt (`..`, absolute path, symlink
    that escapes base after resolution).

    Example:
        safe_join(PDB_LIBRARY_DIR, "1_11_1_7", "1A2B.pdb")  # ok
        safe_join(PDB_LIBRARY_DIR, "../etc/passwd")          # raises ValueError
    """
    base = Path(base).resolve()
    if not parts:
        return base
    # Reject parts that obviously try to break out.
    for p in parts:
        s = str(p)
        if s.startswith("/") or s.startswith("\\") or s.startswith("~"):
            raise ValueError(f"absolute path component rejected: {p!r}")
        if "\x00" in s:
            raise ValueError("null byte in path component")
    target = base.joinpath(*parts).resolve()
    # Containment check using path components, not string prefix
    # (avoids /tmp vs /tmp2 confusion).
    try:
        target.relative_to(base)
    except ValueError as e:
        raise ValueError(
            f"path traversal: target {target!r} not under base {base!r}"
        ) from e
    return target


# ---------------------------------------------------------------------------
# Unified error response — C4 + N-C2
# ---------------------------------------------------------------------------

def error_response(message: str, *, status: int = 500,
                   exc: Optional[BaseException] = None,
                   extra: Optional[dict] = None):
    """
    Build a sanitized JSON error response.

    - Always emits a `request_id` so logs can be correlated to client reports.
    - Logs the full exception server-side (with traceback) at ERROR level.
    - NEVER includes the traceback in the response body unless
      FLASK_EXPOSE_TRACEBACK=1 (dev-only, opt-in).
    - N-C2: refuses to embed stdout/stderr blobs (caller passes them via
      `extra` but they're dropped here unless EXPOSE_TRACEBACK is on).
    """
    rid = uuid.uuid4().hex[:12]
    if exc is not None:
        log.exception("request_id=%s message=%r", rid, message)
    else:
        log.error("request_id=%s message=%r", rid, message)
    client_message = message
    if exc is not None and status >= 500 and not EXPOSE_TRACEBACK:
        client_message = "internal server error"

    body = {
        "status": "error",
        "error": client_message,
        "request_id": rid,
    }
    if EXPOSE_TRACEBACK:
        import traceback
        body["traceback"] = traceback.format_exc() if exc else None
        if extra:
            body["debug_extra"] = extra
    return jsonify(body), status


def install_global_error_handler(app: Flask) -> None:
    """C4 fix: catch every uncaught Exception and emit error_response()."""
    @app.errorhandler(Exception)
    def _on_error(e):
        # Don't double-log HTTPException — they're already routed by Flask.
        from werkzeug.exceptions import HTTPException
        if isinstance(e, HTTPException):
            return e
        return error_response("internal server error", status=500, exc=e)


# ---------------------------------------------------------------------------
# subprocess output sanitizer — N-C2
# ---------------------------------------------------------------------------

def sanitize_subprocess_result(parsed: dict, proc) -> dict:
    """
    N-C2: ensure subprocess stdout/stderr never leaks to client.

    The raw assembler/screen scripts may return JSON containing a `stderr` key
    populated from `proc.stderr[-4000:]`. We strip that and log it instead.
    """
    suspicious_keys = ("stderr", "stdout", "stderr_tail", "stdout_tail")
    leaked = {}
    for k in suspicious_keys:
        if k in parsed:
            leaked[k] = parsed.pop(k)
    if leaked and not EXPOSE_TRACEBACK:
        rid = uuid.uuid4().hex[:12]
        for k, v in leaked.items():
            log.error("subprocess_output request_id=%s key=%s body=%r",
                      rid, k, str(v)[-1000:])
        parsed["request_id"] = rid
    if EXPOSE_TRACEBACK and leaked:
        parsed["_dev_subprocess"] = leaked
    return parsed
