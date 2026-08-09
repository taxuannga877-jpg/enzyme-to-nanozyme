"""Subprocess orchestration for design assembly and catalysis screening."""

import json
import logging
import subprocess
import sys
import uuid
from pathlib import Path

from enzyme_viewer.design_serialization import _json_safe, _loads_subprocess_json
from enzyme_viewer.security import clamp_float, clamp_int, env_int

BASE_DIR = Path(__file__).parent.parent


def _stderr_tail(exc: subprocess.TimeoutExpired):
    stderr = exc.stderr or b""
    if isinstance(stderr, (bytes, bytearray)):
        return stderr[-1000:]
    return str(stderr)[-1000:]


def _timeout_payload(kind: str, message_prefix: str, timeout: int, exc: subprocess.TimeoutExpired):
    rid = uuid.uuid4().hex[:12]
    logging.getLogger("e2n.security").warning(
        "%s subprocess timed out request_id=%s after %ds; partial stderr=%r",
        kind,
        rid,
        timeout,
        _stderr_tail(exc),
    )
    return {
        "status": "error",
        "error": f"{message_prefix} subprocess timed out after {timeout}s",
        "request_id": rid,
    }


def _run_assemble_subprocess(spec, mode="variants"):
    payload = {"design_spec": spec.to_dict(), "mode": mode}
    script = BASE_DIR / "scripts" / "assemble_design_payload.py"
    timeout = env_int(
        "E2N_ASSEMBLY_TIMEOUT",
        600,
        min_value=1,
        max_value=24 * 3600,
    )
    try:
        proc = subprocess.run(
            [sys.executable, str(script)],
            input=json.dumps(payload, ensure_ascii=False, default=_json_safe),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=str(BASE_DIR),
        )
    except subprocess.TimeoutExpired as exc:
        return _timeout_payload("assemble", "assembly", timeout, exc)
    parsed = _loads_subprocess_json(proc)
    if parsed.get("status") == "success":
        return parsed
    if proc.returncode != 0:
        parsed.update(
            {
                "status": "error",
                "error": parsed.get("error")
                or f"assembly subprocess failed with exit code {proc.returncode}",
            }
        )
    return parsed


def _run_catalysis_screen_subprocess(result, data):
    payload = {
        "job_id": result.job_id,
        "atoms": result.atoms,
        "cores": result.cores,
        "second_shell_atoms": result.second_shell_atoms,
        "bond_graph": result.bond_graph,
        "formal_charge": result.formal_charge,
        "spin_multiplicities": result.spin_multiplicities,
        "chemistry_warnings": result.chemistry_warnings,
        "design_spec": result.design_spec.to_dict(),
        "xyz": result.xyz,
        "label": result.label,
        "nanozyme_type": data.get("nanozyme_type") or result.design_spec.nanozyme_type,
        "max_adsorption_poses": clamp_int(
            data.get("max_adsorption_poses", 6),
            default=6,
            min_value=1,
            max_value=20,
        ),
        "run_neb": bool(data.get("run_neb", False)),
        "neb_steps": clamp_int(
            data.get("neb_steps", 80),
            default=80,
            min_value=5,
            max_value=500,
        ),
        "neb_fmax": clamp_float(
            data.get("neb_fmax", 0.08),
            default=0.08,
            min_value=0.01,
            max_value=1.0,
        ),
        "run_reaction_scan": bool(data.get("run_reaction_scan", False)),
        "reaction_scan_points": clamp_int(
            data.get("reaction_scan_points", 5),
            default=5,
            min_value=3,
            max_value=25,
        ),
    }
    script = BASE_DIR / "scripts" / "screen_catalysis_payload.py"
    timeout = env_int(
        "E2N_CATALYSIS_TIMEOUT",
        600,
        min_value=1,
        max_value=24 * 3600,
    )
    try:
        proc = subprocess.run(
            [sys.executable, str(script)],
            input=json.dumps(payload, ensure_ascii=False, default=_json_safe),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=str(BASE_DIR),
        )
    except subprocess.TimeoutExpired as exc:
        return _timeout_payload("catalysis", "catalysis", timeout, exc)
    parsed = _loads_subprocess_json(proc)
    if parsed.get("status") == "success":
        return parsed
    if proc.returncode != 0:
        parsed.update(
            {
                "status": "error",
                "error": parsed.get("error")
                or f"catalysis subprocess failed with exit code {proc.returncode}",
            }
        )
    return parsed
