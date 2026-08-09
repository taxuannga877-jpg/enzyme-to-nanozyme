"""Persistence and lookup helpers for nanozyme design results."""

import json
import os
from pathlib import Path
from typing import Optional

from enzyme_viewer.design_io import (
    _atoms_to_xyz,
    _parse_pdb_atoms,
    _parse_xyz_atoms,
    _reconstruct_cores_for_loaded_design,
)
from enzyme_viewer.design_serialization import _json_safe, _score_from_persisted_payload
from enzyme_viewer.security import safe_join
from nanozyme_mining.design import structure_exporter
from nanozyme_mining.design.chemical_system import infer_bond_graph
from nanozyme_mining.design.design_spec import DesignSpec
from nanozyme_mining.design.nanozyme_assembler import AssemblyResult
from nanozyme_mining.design.potential_evaluator import (
    PotentialEvaluationConfig,
    evaluate_assembly,
)
from nanozyme_mining.utils.constants import CATALYTIC_METAL_ELEMENTS

BASE_DIR = Path(__file__).parent.parent
_CONFIG = None
_DESIGN_RESULTS = None


def configure_design_store(*, config, design_results) -> None:
    """Bind store helpers to the app-owned config and in-memory cache."""
    global _CONFIG, _DESIGN_RESULTS
    _CONFIG = config
    _DESIGN_RESULTS = design_results


def _config():
    if _CONFIG is None:
        raise RuntimeError("design store has not been configured")
    return _CONFIG


def _design_cache():
    if _DESIGN_RESULTS is None:
        raise RuntimeError("design store has not been configured")
    return _DESIGN_RESULTS


def _activity_validation_task_dir(task_id: str) -> Path:
    return safe_join(_config()["ACTIVITY_VALIDATION_OUTPUT_DIR"], task_id)


def _design_result_dir(job_id: str) -> Path:
    return safe_join(_config()["DESIGN_OUTPUT_DIR"], job_id)


def _score_payload(score):
    return {
        "geometry": round(score.geometry_score, 3),
        "energy": round(score.energy_score, 3),
        "coordination": round(score.coordination_score, 3),
        "steric": round(getattr(score, "steric_score", 0.0), 3),
        "stability": round(getattr(score, "stability_score", score.energy_score), 3),
        "total": round(score.total_score, 3),
        "method": getattr(score, "method", "unknown"),
        "backend": getattr(score, "backend", "unknown"),
        "raw_energy_ev": getattr(score, "raw_energy_ev", None),
        "energy_per_atom_ev": getattr(score, "energy_per_atom_ev", None),
        "relaxed_energy_ev": getattr(score, "relaxed_energy_ev", None),
        "max_force_ev_per_a": getattr(score, "max_force_ev_per_a", None),
        "warnings": getattr(score, "warnings", []),
        "details": getattr(score, "details", {}),
    }


def _assembly_result_payload(result: AssemblyResult) -> dict:
    return {
        "job_id": result.job_id,
        "label": result.label,
        "atom_count": len(result.atoms),
        "atoms": result.atoms,
        "cores": result.cores,
        "second_shell_atoms": result.second_shell_atoms,
        "bond_graph": result.bond_graph,
        "smiles": result.smiles,
        "xyz": result.xyz,
        "formal_charge": result.formal_charge,
        "spin_multiplicities": result.spin_multiplicities,
        "chemistry_warnings": result.chemistry_warnings,
        "score": _score_payload(result.score),
        "multi_metal_mode": result.design_spec.multi_metal_mode,
    }


def _read_json_file(path: Path) -> Optional[dict]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _find_validation_dir_for_job(job_id: str) -> Optional[Path]:
    root = _config()["ACTIVITY_VALIDATION_OUTPUT_DIR"]
    for report_path in sorted(root.glob("*/activity_validation_report.json"), reverse=True):
        report = _read_json_file(report_path)
        if report and report.get("job_id") == job_id:
            return report_path.parent
    return None


def _load_design_result_from_disk(job_id: str) -> Optional[AssemblyResult]:
    try:
        design_dir = _design_result_dir(job_id)
    except ValueError:
        return None

    source_dir = design_dir
    assembly_payload = _read_json_file(source_dir / "assembly_result.json")
    spec_payload = _read_json_file(source_dir / "design_spec.json")
    validation_payload = None

    if not spec_payload:
        source_dir = _find_validation_dir_for_job(job_id)
        if not source_dir:
            return None
        spec_payload = _read_json_file(source_dir / "design_spec.json")
        validation_payload = _read_json_file(source_dir / "activity_validation_report.json")
        assembly_payload = _read_json_file(source_dir / "assembly_result.json") or {}
    if not spec_payload:
        return None

    try:
        spec = DesignSpec.from_dict(spec_payload)
    except Exception:
        return None

    persisted = assembly_payload or validation_payload or {}
    pdb_path = source_dir / f"{job_id}.pdb"
    xyz_path = source_dir / f"{job_id}.xyz"
    atoms = persisted.get("atoms") or _parse_pdb_atoms(pdb_path)
    xyz_text = persisted.get("xyz") or ""
    if xyz_path.exists():
        xyz_atoms, disk_xyz_text = _parse_xyz_atoms(xyz_path)
        if not xyz_text:
            xyz_text = disk_xyz_text
        if not atoms:
            atoms = xyz_atoms
    if not atoms:
        return None
    if not xyz_text:
        xyz_text = _atoms_to_xyz(atoms)

    cores = persisted.get("cores") or _reconstruct_cores_for_loaded_design(atoms, spec)
    if assembly_payload and assembly_payload.get("score"):
        score = _score_from_persisted_payload(assembly_payload.get("score") or {})
    else:
        try:
            score = evaluate_assembly(
                {
                    "atoms": atoms,
                    "cores": cores,
                    "linker_atoms": [],
                    "mode": spec.multi_metal_mode,
                    "second_shell": [],
                },
                spec,
                config=PotentialEvaluationConfig(backend="geometry_proxy", relax=False),
            )
            score.warnings.append("assembly_score_recomputed_from_persisted_structure")
        except Exception as exc:
            score = _score_from_persisted_payload({})
            score.warnings.append(f"assembly_score_recompute_failed: {exc}")
    if not assembly_payload:
        score.method = "persisted_validation_context"
        score.backend = "disk"
    result = AssemblyResult(
        job_id=job_id,
        atoms=atoms,
        cores=cores,
        second_shell_atoms=persisted.get("second_shell_atoms") or [],
        bond_graph=persisted.get("bond_graph") or infer_bond_graph(atoms),
        formal_charge=int(
            persisted.get("formal_charge", sum(m.oxidation_state for m in spec.metals))
        ),
        spin_multiplicities=persisted.get("spin_multiplicities") or [1],
        chemistry_warnings=persisted.get("chemistry_warnings") or [],
        score=score,
        design_spec=spec,
        smiles=persisted.get("smiles", ""),
        xyz=xyz_text,
        label=persisted.get("label", job_id),
    )
    _design_cache()[job_id] = result
    return result


def _get_design_result(job_id: str) -> Optional[AssemblyResult]:
    result = _design_cache().get(job_id)
    if result:
        return result
    return _load_design_result_from_disk(job_id)


def _persist_design_result(result: AssemblyResult) -> Path:
    output_dir = _design_result_dir(result.job_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    structure_exporter.export(result, str(output_dir))
    (output_dir / "design_spec.json").write_text(
        json.dumps(
            result.design_spec.to_dict(),
            indent=2,
            ensure_ascii=False,
            default=_json_safe,
        )
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "assembly_result.json").write_text(
        json.dumps(
            _assembly_result_payload(result),
            indent=2,
            ensure_ascii=False,
            default=_json_safe,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def _persisted_validation_artifacts(task_dir: Path) -> list:
    artifacts = []
    figure_specs = [
        (
            "overview",
            "Activity validation overview",
            "activity_validation_overview.png",
            "activity_validation_overview.svg",
        ),
        (
            "reaction_profile",
            "Reaction-coordinate screening profiles",
            "reaction_coordinate_profiles.png",
            "reaction_coordinate_profiles.svg",
        ),
        (
            "mechanism_panel",
            "Mechanism-specific coordinate panels",
            "mechanism_coordinate_panels.png",
            "mechanism_coordinate_panels.svg",
        ),
        (
            "redox_state_map",
            "Redox charge/spin state map",
            "redox_state_heatmap.png",
            "redox_state_heatmap.svg",
        ),
        (
            "adsorption_volcano",
            "Adsorption-activation volcano",
            "adsorption_activation_volcano.png",
            "adsorption_activation_volcano.svg",
        ),
    ]
    for kind, label, png, svg in figure_specs:
        if (task_dir / png).exists():
            item = {"kind": kind, "label": label, "png": png}
            if (task_dir / svg).exists():
                item["svg"] = svg
            artifacts.append(item)

    for suffix, label in (
        (".pdb", "Validated design structure (PDB)"),
        (".xyz", "Validated design structure (XYZ)"),
        (".sdf", "Validated design structure (SDF)"),
    ):
        files = sorted(task_dir.glob(f"*{suffix}"))
        if files:
            key = suffix.lstrip(".")
            artifacts.append({"kind": "structure", "label": label, key: files[0].name})

    if (task_dir / "design_spec.json").exists():
        artifacts.append(
            {"kind": "design_spec", "label": "Design specification", "json": "design_spec.json"}
        )
    if (task_dir / "assembly_result.json").exists():
        artifacts.append(
            {
                "kind": "assembly_result",
                "label": "Assembly result metadata",
                "json": "assembly_result.json",
            }
        )
    if (task_dir / "activity_validation_report.json").exists():
        artifacts.append(
            {
                "kind": "json",
                "label": "Machine-readable validation report",
                "json": "activity_validation_report.json",
            }
        )
    return artifacts


def _validation_snapshot_from_disk(task_id: str) -> Optional[dict]:
    try:
        task_dir = _activity_validation_task_dir(task_id)
    except ValueError:
        return None
    report_path = task_dir / "activity_validation_report.json"
    if not report_path.exists():
        return None
    try:
        summary = json.loads(report_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    updated_at = report_path.stat().st_mtime
    activities = [
        row.get("activity")
        for row in summary.get("activity_results", [])
        if row.get("activity")
    ]
    return {
        "task_id": task_id,
        "job_id": summary.get("job_id"),
        "activities": activities,
        "status": "complete",
        "stage": "complete",
        "progress": 100.0,
        "error": None,
        "result": summary,
        "partial_results": summary.get("activity_results", []),
        "artifacts": _persisted_validation_artifacts(task_dir),
        "output_dir": str(task_dir),
        "created_at": updated_at,
        "updated_at": updated_at,
        "events": [
            {
                "time": updated_at,
                "stage": "complete",
                "message": "loaded completed validation result from disk",
            }
        ],
    }


def _activity_validation_runtime_context(result: AssemblyResult) -> dict:
    backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
    backend = backend or "geometry_proxy"
    model = os.environ.get("E2N_MLP_MODEL") or ""
    mace_model = BASE_DIR / "models" / "mace-mh-1.model"
    try:
        import tblite  # noqa: F401

        tblite_available = True
    except Exception:
        tblite_available = False
    try:
        import mace  # noqa: F401

        mace_available = True
    except Exception:
        mace_available = False
    return {
        "configured_backend": backend,
        "configured_model": model,
        "mace_available": mace_available,
        "mace_model_available": mace_model.exists(),
        "tblite_available": tblite_available,
        "recommended_flow": "MACE relaxation -> tblite / GFN2-xTB reaction screening",
    }


def _activity_validation_structure_diagnostics(result: AssemblyResult) -> dict:
    metal_count = sum(
        1
        for atom in result.atoms
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    )
    activities = list(getattr(result.design_spec, "activities", []) or [])
    atom_count = len(result.atoms)
    representative_dual = atom_count >= 100 and metal_count >= 2 and len(activities) >= 2
    warnings = []
    if atom_count < 100:
        warnings.append(
            "preview_scale_scaffold: fewer than 100 atoms; not representative of the paper-style dual-metal carbon support"
        )
    if metal_count < 2 and len(activities) >= 2:
        warnings.append("dual_activity_requested_but_less_than_two_catalytic_metals")
    if getattr(result.score, "backend", "") == "geometry_proxy":
        warnings.append("stored_score_uses_geometry_proxy; rerun with MACE/tblite for calculation figures")
    return {
        "atom_count": atom_count,
        "metal_count": metal_count,
        "activity_count": len(activities),
        "representative_dual_activity_scaffold": representative_dual,
        "warnings": warnings,
    }


def _activity_validation_reference_figures() -> list:
    figures = []
    for key, spec in (_config().get("ACTIVITY_VALIDATION_REFERENCE_FIGURES") or {}).items():
        path = Path(spec.get("path") or "")
        if path.exists():
            figures.append(
                {
                    "key": key,
                    "label": spec.get("label") or key,
                    "url": f"/api/design/activity_validation/reference/{key}",
                }
            )
    return figures
