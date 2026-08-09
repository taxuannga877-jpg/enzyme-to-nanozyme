"""Serialization helpers for nanozyme design API payloads."""

import datetime as _dt
import json
import logging
import uuid
from pathlib import Path

from enzyme_viewer.security import sanitize_subprocess_result
from nanozyme_mining.design.constraint_scorer import ConstraintScore
from nanozyme_mining.design.design_spec import DesignSpec
from nanozyme_mining.design.nanozyme_assembler import AssemblyResult


def _json_safe(obj):
    """JSON fallback for numpy, stdlib, and design dataclass payloads."""
    try:
        import numpy as np

        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
    except ImportError:
        pass
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (_dt.datetime, _dt.date)):
        return obj.isoformat()
    if isinstance(obj, (set, frozenset)):
        return list(obj)
    if isinstance(obj, bytes):
        try:
            return obj.decode("utf-8")
        except UnicodeDecodeError:
            return obj.hex()
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _score_from_dict(data):
    allowed = ConstraintScore.__dataclass_fields__.keys()
    return ConstraintScore(**{key: value for key, value in (data or {}).items() if key in allowed})


def _score_from_persisted_payload(data):
    data = data or {}
    return ConstraintScore(
        geometry_score=float(data.get("geometry_score", data.get("geometry", 0.0)) or 0.0),
        energy_score=float(data.get("energy_score", data.get("energy", 0.0)) or 0.0),
        coordination_score=float(
            data.get("coordination_score", data.get("coordination", 0.0)) or 0.0
        ),
        steric_score=float(data.get("steric_score", data.get("steric", 0.0)) or 0.0),
        stability_score=float(
            data.get("stability_score", data.get("stability", data.get("energy", 0.0))) or 0.0
        ),
        total_score=float(data.get("total_score", data.get("total", 0.0)) or 0.0),
        passed_hard_constraints=bool(data.get("passed_hard_constraints", True)),
        errors=list(data.get("errors") or []),
        warnings=list(data.get("warnings") or []),
        method=data.get("method", "persisted"),
        backend=data.get("backend", "disk"),
        raw_energy_ev=data.get("raw_energy_ev"),
        energy_per_atom_ev=data.get("energy_per_atom_ev"),
        relaxed_energy_ev=data.get("relaxed_energy_ev"),
        max_force_ev_per_a=data.get("max_force_ev_per_a"),
        details=dict(data.get("details") or {}),
    )


def _assembly_from_dict(data):
    design_spec = DesignSpec.from_dict(data["design_spec"])
    return AssemblyResult(
        job_id=data["job_id"],
        atoms=data["atoms"],
        cores=data.get("cores", []),
        second_shell_atoms=data.get("second_shell_atoms", []),
        bond_graph=data.get("bond_graph", []),
        formal_charge=int(data.get("formal_charge", 0)),
        spin_multiplicities=data.get("spin_multiplicities", [1]),
        chemistry_warnings=data.get("chemistry_warnings", []),
        score=_score_from_dict(data.get("score")),
        design_spec=design_spec,
        smiles=data.get("smiles", ""),
        xyz=data.get("xyz", ""),
        label=data.get("label", ""),
        error=data.get("error"),
    )


def _loads_subprocess_json(proc):
    try:
        parsed = json.loads(proc.stdout)
        return sanitize_subprocess_result(parsed, proc)
    except json.JSONDecodeError:
        rid = uuid.uuid4().hex[:12]
        logging.getLogger("e2n.security").error(
            "subprocess non-JSON output request_id=%s stderr=%r stdout=%r",
            rid,
            proc.stderr[-2000:],
            proc.stdout[-2000:],
        )
        return {
            "status": "error",
            "error": "subprocess returned non-JSON output",
            "request_id": rid,
        }
