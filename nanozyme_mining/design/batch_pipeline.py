"""Resumable physicochemical batch screening and legacy gallery audit."""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
from datetime import datetime
import itertools
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .catalysis_screening import screen_catalysis
from .design_spec import CoordAtomSpec, DesignSpec, MetalSpec
from .nanozyme_assembler import AssemblyResult, NanozymeAssembler
from .physchem_knowledge import evaluate_constructibility, get_screening_proxy_policy, knowledge_version
from .potential_evaluator import PotentialEvaluationConfig
from . import structure_exporter
from .substrate_catalog import get_reaction_task
from .validation import ValidationReport, validate_assembly
from ..utils.constants import CATALYTIC_METAL_ELEMENTS


RESULT_STATUSES = (
    "relaxed",
    "mace_failed",
    "adsorption_screened",
    "reaction_scanned",
    "redox_state_scanned",
    "not_applicable",
    "failed_with_reason",
)

PIPELINE_STAGES = (
    "requested",
    "constructible",
    "preflight_pass",
    "MACE_pass",
    "tblite_light",
    "full_reaction_scan",
    "experimental_candidate",
)


@dataclass(frozen=True)
class PhyschemScreeningConfig:
    """Public V2 configuration for MACE relaxation + tblite reaction screening."""

    activities: Tuple[str, ...] = ()
    candidates_per_combo: int = 0
    mace_model_path: str = "models/mace-mh-1.model"
    tblite_method: str = "GFN2-xTB"
    charge_scan: Tuple[int, ...] = (-1, 0, 1)
    spin_multiplicities: Tuple[int, ...] = (1, 2, 3, 4, 5)
    max_adsorption_poses: int = 8
    top_n_mace_per_combo: int = 6
    top_n_tblite_per_combo: int = 3
    run_reaction_scan_fraction: float = 0.25
    run_full_reaction_scan: bool = False
    reaction_scan_points: int = 5
    render_figures: bool = True

    def to_dict(self) -> Dict[str, Any]:
        payload = asdict(self)
        for key in ("activities", "charge_scan", "spin_multiplicities"):
            payload[key] = list(payload[key])
        return payload


@dataclass
class BatchRecord:
    record_id: str
    combo: str
    activities: Tuple[str, str]
    mode: str
    stage: str
    status: str
    calculation_status: str = "not_applicable"
    label: str = ""
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    score: Dict[str, Any] = field(default_factory=dict)
    diagnostics: Dict[str, Any] = field(default_factory=dict)
    output_files: Dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "record_id": self.record_id,
            "combo": self.combo,
            "activities": list(self.activities),
            "mode": self.mode,
            "stage": self.stage,
            "status": self.status,
            "calculation_status": self.calculation_status,
            "label": self.label,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
            "score": dict(self.score),
            "diagnostics": dict(self.diagnostics),
            "output_files": dict(self.output_files),
        }


def run_physchem_screening_batch(
    *,
    activity_templates_path: Optional[Path | str] = None,
    output_root: Path | str = "outputs/physchem_screening_batches",
    run_tblite_light: bool = False,
    run_reaction_scan_fraction: float = 0.25,
    max_adsorption_poses: int = 8,
    resume_dir: Optional[Path | str] = None,
    activities: Optional[Sequence[str]] = None,
    candidates_per_combo: int = 0,
    run_full_reaction_scan: bool = False,
    reaction_scan_points: int = 5,
    mace_model_path: str = "models/mace-mh-1.model",
    tblite_method: str = "GFN2-xTB",
    top_n_mace_per_combo: int = 6,
    top_n_tblite_per_combo: int = 3,
    render_figures: bool = True,
) -> Dict[str, Any]:
    """Generate a resumable MACE + tblite V2 screening batch manifest."""
    templates = load_activity_templates(activity_templates_path)
    if activities:
        selected = {str(activity) for activity in activities}
        templates = {key: value for key, value in templates.items() if key in selected}
    run_dir = Path(resume_dir) if resume_dir else _new_run_dir(Path(output_root))
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / "manifest.json"
    manifest = _load_manifest(manifest_path) or _new_manifest(run_dir, templates)
    batch_config = PhyschemScreeningConfig(
        activities=tuple(sorted(templates)),
        candidates_per_combo=max(0, int(candidates_per_combo)),
        mace_model_path=str(mace_model_path),
        tblite_method=tblite_method,
        max_adsorption_poses=max_adsorption_poses,
        top_n_mace_per_combo=top_n_mace_per_combo,
        top_n_tblite_per_combo=top_n_tblite_per_combo,
        run_reaction_scan_fraction=run_reaction_scan_fraction,
        run_full_reaction_scan=run_full_reaction_scan,
        reaction_scan_points=reaction_scan_points,
        render_figures=render_figures,
    )
    manifest["calculation_config"] = batch_config.to_dict()
    manifest["result_statuses"] = list(RESULT_STATUSES)
    completed = {record["record_id"] for record in manifest.get("records", [])}

    assembler = NanozymeAssembler()
    for activity_a, activity_b in itertools.combinations(sorted(templates), 2):
        combo = f"{activity_a} + {activity_b}"
        spec = design_spec_for_combo(activity_a, activity_b, templates)
        decision = evaluate_constructibility(spec)
        if not decision.constructible:
            record_id = f"{_slug(combo)}:not_constructible"
            if record_id not in completed:
                manifest["records"].append(
                    BatchRecord(
                        record_id=record_id,
                        combo=combo,
                        activities=(activity_a, activity_b),
                        mode="not_constructible",
                        stage="constructible",
                        status="not_constructible",
                        errors=list(decision.reason_codes),
                        warnings=list(decision.warnings),
                        diagnostics={"constructibility": decision.to_dict()},
                    ).to_dict()
                )
                _write_manifest(manifest_path, manifest)
            continue

        per_combo_limit = candidates_per_combo or top_n_mace_per_combo
        for idx, result in enumerate(assembler.assemble_batch(spec)):
            if per_combo_limit and idx >= per_combo_limit:
                break
            mode = _result_mode(result)
            record_id = f"{_slug(combo)}:{mode}"
            if record_id in completed:
                continue
            record = _record_from_result(record_id, combo, (activity_a, activity_b), mode, result)
            run_record_scan = run_full_reaction_scan and idx < max(1, top_n_tblite_per_combo)
            if (run_tblite_light or run_record_scan) and result.score.passed_hard_constraints:
                gate = _mace_tblite_gate(result, require_converged=run_record_scan)
                record.diagnostics["mace_tblite_gate"] = gate
                if gate["passed"]:
                    record = _attach_tblite_light(
                        record,
                        result,
                        max_adsorption_poses,
                        run_reaction_scan=run_record_scan,
                        reaction_scan_points=reaction_scan_points,
                        tblite_method=tblite_method,
                    )
                else:
                    record.status = "rejected"
                    record.stage = "MACE_pass"
                    record.calculation_status = "mace_failed"
                    record.diagnostics["calculation_status"] = record.calculation_status
                    record.errors.append(gate["reason"])
            record = _persist_batch_candidate(run_dir, record, result)
            manifest["records"].append(record.to_dict())
            completed.add(record_id)
            _write_manifest(manifest_path, manifest)

    manifest["rankings"] = rank_experimental_candidates(
        manifest["records"],
        run_reaction_scan_fraction=run_reaction_scan_fraction,
    )
    manifest["artifacts"] = (
        render_physchem_batch_figures(manifest, run_dir) if render_figures else []
    )
    _write_manifest(manifest_path, manifest)
    return manifest


def audit_legacy_bimetallic_gallery(
    gallery_dir: Path | str = "outputs/bimetallic_structure_gallery/latest",
    *,
    output_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    """Audit old gallery PDBs with the strict validation report."""
    root = Path(gallery_dir)
    templates = load_activity_templates(root / "activity_templates.json")
    records = []
    index_path = root / "index.csv"
    with index_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            pdb_path = Path(row["pdb"])
            atoms, conect = parse_pdb_atoms_and_bonds(pdb_path)
            metals = [
                atom for atom in atoms
                if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
            ]
            spec = design_spec_for_combo(row["activity_a"], row["activity_b"], templates)
            cores = _cores_from_legacy_pdb(atoms, conect, spec)
            report = validate_assembly(
                {"atoms": atoms, "cores": cores},
                spec,
                stage="legacy_gallery_audit",
            )
            records.append(
                {
                    "combo": row["combo"],
                    "activity_a": row["activity_a"],
                    "activity_b": row["activity_b"],
                    "mode": row["mode"],
                    "pdb": str(pdb_path),
                    "metal_count": len(metals),
                    "passed": report.passed,
                    "reason_codes": report.reason_codes,
                    "centers": [center.to_dict() for center in report.centers],
                }
            )
    summary = {
        "schema_version": "1.0.0",
        "knowledge_version": knowledge_version(),
        "gallery_dir": str(root),
        "record_count": len(records),
        "passed_count": sum(1 for record in records if record["passed"]),
        "failed_count": sum(1 for record in records if not record["passed"]),
        "cu_failed_count": sum(
            1
            for record in records
            if any(center["metal"] == "CU" and not center["passed"] for center in record["centers"])
        ),
        "reason_code_counts": _reason_code_counts(records),
        "records": records,
    }
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def load_activity_templates(path: Optional[Path | str] = None) -> Dict[str, Dict[str, Any]]:
    if path is None:
        path = Path("outputs/bimetallic_structure_gallery/latest/activity_templates.json")
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_pdb_atoms_and_bonds(path: Path | str) -> Tuple[List[Dict[str, Any]], Dict[int, List[int]]]:
    atoms: List[Dict[str, Any]] = []
    serial_to_index: Dict[int, int] = {}
    bonds: Dict[int, List[int]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                element = line[76:78].strip().upper() or _element_from_atom_name(atom_name)
                atom = {
                    "element": element,
                    "atom_name": atom_name,
                    "residue_name": residue_name,
                    "coords": [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    "formal_charge": 0,
                    "pdb_serial": serial,
                }
                serial_to_index[serial] = len(atoms)
                atoms.append(atom)
            elif line.startswith("CONECT"):
                values = [int(value) for value in line.split()[1:]]
                if not values:
                    continue
                left_serial = values[0]
                if left_serial not in serial_to_index:
                    continue
                left = serial_to_index[left_serial]
                for right_serial in values[1:]:
                    if right_serial in serial_to_index:
                        right = serial_to_index[right_serial]
                        bonds.setdefault(left, [])
                        bonds.setdefault(right, [])
                        if right not in bonds[left]:
                            bonds[left].append(right)
                        if left not in bonds[right]:
                            bonds[right].append(left)
    return atoms, bonds


def rank_experimental_candidates(
    records: Iterable[Dict[str, Any]],
    *,
    run_reaction_scan_fraction: float = 0.25,
    candidates_per_activity: int = 5,
) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for record in records:
        if record.get("status") != "passed":
            continue
        for activity in record.get("activities", []):
            grouped.setdefault(activity, []).append(record)
    rankings = {}
    for activity, rows in grouped.items():
        front = _pareto_front(rows)
        front = sorted(front, key=_ranking_key, reverse=True)
        cutoff = max(3, int(round(len(rows) * run_reaction_scan_fraction))) if rows else 0
        selected_ids = {row["record_id"] for row in front[:max(candidates_per_activity, cutoff)]}
        rankings[activity] = [
            {
                "record_id": row["record_id"],
                "combo": row["combo"],
                "mode": row["mode"],
                "rank": rank + 1,
                "selection": "experimental_candidate"
                if row["record_id"] in selected_ids and rank < candidates_per_activity
                else "full_reaction_scan_candidate",
                "score": row.get("score", {}),
                "diagnostic_summary": {
                    "validation_codes": row.get("diagnostics", {})
                    .get("physchem_validation", {})
                    .get("reason_codes", []),
                    "relaxation_status": row.get("diagnostics", {})
                    .get("structure_relaxation", {})
                    .get("relaxation_status"),
                },
            }
            for rank, row in enumerate(front)
        ]
    return rankings


def design_spec_for_combo(
    activity_a: str,
    activity_b: str,
    templates: Dict[str, Dict[str, Any]],
) -> DesignSpec:
    return DesignSpec(
        nanozyme_type=f"{activity_a} + {activity_b}",
        ec_numbers=[],
        activities=[activity_a, activity_b],
        metals=[
            _metal_spec_from_template(activity_a, templates[activity_a]),
            _metal_spec_from_template(activity_b, templates[activity_b]),
        ],
    )


_design_spec_for_combo = design_spec_for_combo


def _metal_spec_from_template(activity: str, template: Dict[str, Any]) -> MetalSpec:
    return MetalSpec(
        metal_type=template["metal_type"],
        oxidation_state=int(template["oxidation_state"]),
        coordination_geometry=template["coordination_geometry"],
        coordination_number=int(template["coordination_number"]),
        coord_atoms=[CoordAtomSpec(**atom) for atom in template["coord_atoms"]],
        activity_type=activity,
        prototype_id=template.get("prototype_id"),
        condition_id=template.get("condition_id"),
        microstate_id=(template.get("microstates") or [None])[0],
    )


def _record_from_result(
    record_id: str,
    combo: str,
    activities: Tuple[str, str],
    mode: str,
    result: AssemblyResult,
) -> BatchRecord:
    relaxation = result.score.details.get("structure_relaxation", {})
    stage = "preflight_pass" if result.score.passed_hard_constraints else "constructible"
    if relaxation.get("status") == "success" and result.score.passed_hard_constraints:
        stage = "MACE_pass"
    status = "passed" if result.score.passed_hard_constraints else "rejected"
    calculation_status = _relaxation_calculation_status(relaxation)
    return BatchRecord(
        record_id=record_id,
        combo=combo,
        activities=activities,
        mode=mode,
        stage=stage,
        status=status,
        calculation_status=calculation_status,
        label=result.label,
        errors=list(result.score.errors),
        warnings=list(result.score.warnings),
        score={
            "total": result.score.total_score,
            "geometry": result.score.geometry_score,
            "coordination": result.score.coordination_score,
            "distance": getattr(result.score, "distance_score", result.score.energy_score),
            "steric": result.score.steric_score,
        },
        diagnostics={
            "constructibility": result.score.details.get("constructibility"),
            "physchem_validation": result.score.details.get("physchem_validation"),
            "preflight_validation": result.score.details.get("preflight_validation"),
            "structure_relaxation": relaxation,
            "bimetallic": result.score.details.get("bimetallic"),
            "calculation_status": calculation_status,
        },
    )


def _attach_tblite_light(
    record: BatchRecord,
    result: AssemblyResult,
    max_adsorption_poses: int,
    *,
    run_reaction_scan: bool = False,
    reaction_scan_points: int = 5,
    tblite_method: str = "GFN2-xTB",
) -> BatchRecord:
    activity_payloads = []
    tblite_config = PotentialEvaluationConfig(
        backend="tblite",
        model=tblite_method,
        relax=False,
    )
    for activity in record.activities:
        task = get_reaction_task(activity)
        if task is None:
            record.warnings.append(f"No reaction task for {activity}")
            continue
        try:
            payload = screen_catalysis(
                result,
                task=task,
                config=tblite_config,
                max_adsorption_poses=max_adsorption_poses,
                run_reaction_scan=run_reaction_scan,
                reaction_scan_points=reaction_scan_points,
            )
        except Exception as exc:
            payload = {
                "status": "error",
                "calculation_status": "failed_with_reason",
                "error": str(exc),
                "task": task.to_dict(),
            }
        activity_payloads.append({"activity": activity, "payload": payload})

    record.stage = "full_reaction_scan" if run_reaction_scan else "tblite_light"
    summaries = [_activity_calculation_summary(item) for item in activity_payloads]
    record.diagnostics["activity_calculations"] = summaries
    record.diagnostics["tblite_light"] = {
        "status": "success" if summaries else "not_applicable",
        "activity_count": len(summaries),
        "best_adsorption_energy_ev": _best_activity_metric(
            summaries, "best_adsorption_energy_ev", prefer_min=True
        ),
        "max_activation_metric_ev": _best_activity_metric(
            summaries, "activation_metric_ev", prefer_min=False
        ),
    }
    status_priority = [
        "reaction_scanned",
        "redox_state_scanned",
        "adsorption_screened",
        "failed_with_reason",
        "not_applicable",
    ]
    seen = {summary.get("calculation_status") for summary in summaries}
    for status in status_priority:
        if status in seen:
            record.calculation_status = status
            break
    record.diagnostics["calculation_status"] = record.calculation_status
    if record.calculation_status == "failed_with_reason":
        record.warnings.append("tblite reaction screening failed for at least one activity")
    return record


def _mace_tblite_gate(
    result: AssemblyResult,
    *,
    require_converged: bool,
) -> Dict[str, Any]:
    relaxation = result.score.details.get("structure_relaxation", {})
    validation = result.score.details.get("physchem_validation", {})
    if str(relaxation.get("backend") or "").lower() != "mace":
        return {
            "passed": False,
            "reason": "tblite screening requires a completed MACE structure relaxation",
        }
    if relaxation.get("status") != "success":
        return {
            "passed": False,
            "reason": f"MACE relaxation status is {relaxation.get('status', 'missing')}",
        }
    relaxation_status = relaxation.get("relaxation_status")
    allowed = (
        {"converged", "converged_constrained"}
        if require_converged
        else {
            "converged",
            "converged_constrained",
            "usable_not_converged",
            "usable_constrained",
        }
    )
    if relaxation_status not in allowed:
        expected = "converged" if require_converged else "converged or usable_not_converged"
        return {
            "passed": False,
            "reason": (
                f"MACE geometry is {relaxation_status or 'unclassified'}; "
                f"tblite stage requires {expected}"
            ),
        }
    if validation and not validation.get("passed", False):
        return {
            "passed": False,
            "reason": "post-MACE physicochemical validation failed",
        }
    return {
        "passed": True,
        "reason": None,
        "relaxation_status": relaxation_status,
        "max_force_ev_per_a": relaxation.get("raw_mace_max_force_ev_per_a"),
    }


def _relaxation_calculation_status(relaxation: Dict[str, Any]) -> str:
    status = relaxation.get("status")
    backend = str(relaxation.get("backend") or "").lower()
    if status == "success":
        return "relaxed"
    if backend == "mace" and status == "failed":
        return "mace_failed"
    if status == "failed":
        return "failed_with_reason"
    return "not_applicable"


def _activity_calculation_summary(entry: Dict[str, Any]) -> Dict[str, Any]:
    activity = entry.get("activity")
    payload = entry.get("payload") or {}
    task = payload.get("task") or {}
    calculation = task.get("calculation") or {}
    redox = payload.get("redox_state_profile") or {}
    reaction = payload.get("reaction_profile") or {}
    activation = reaction.get("proxy_barrier_ev")
    if activation is None:
        activation = redox.get("redox_activation_energy_ev")
    best_adsorption = _best_adsorption_energy(payload)
    return {
        "activity": activity,
        "status": payload.get("status"),
        "calculation_status": payload.get("calculation_status"),
        "task_id": task.get("task_id"),
        "mechanism_family": calculation.get("mechanism_family"),
        "barrier_method": calculation.get("barrier_method"),
        "ml_backend": payload.get("ml_backend"),
        "active_center": payload.get("active_center"),
        "best_adsorption_energy_ev": best_adsorption,
        "activation_metric_ev": activation,
        "reaction_profile_status": reaction.get("status"),
        "redox_state_profile_status": redox.get("status"),
        "electronic_cluster": payload.get("electronic_cluster"),
        "mechanism_visualization": payload.get("mechanism_visualization"),
        "error": payload.get("error"),
    }


def _best_activity_metric(
    summaries: Sequence[Dict[str, Any]],
    key: str,
    *,
    prefer_min: bool,
) -> Optional[float]:
    values = [
        float(summary[key])
        for summary in summaries
        if summary.get(key) is not None and np.isfinite(float(summary[key]))
    ]
    if not values:
        return None
    return min(values) if prefer_min else max(values)


def _persist_batch_candidate(
    run_dir: Path,
    record: BatchRecord,
    result: AssemblyResult,
) -> BatchRecord:
    candidate_dir = run_dir / "candidates" / _slug(record.record_id)
    candidate_dir.mkdir(parents=True, exist_ok=True)
    structure_paths = structure_exporter.export(result, str(candidate_dir))
    payload = {
        "record": record.to_dict(),
        "assembly": {
            "job_id": result.job_id,
            "label": result.label,
            "atom_count": len(result.atoms),
            "atoms": result.atoms,
            "cores": result.cores,
            "formal_charge": result.formal_charge,
            "spin_multiplicities": list(result.spin_multiplicities),
            "chemistry_warnings": list(result.chemistry_warnings),
            "design_spec": result.design_spec.to_dict(),
            "score": {
                "total": result.score.total_score,
                "geometry": result.score.geometry_score,
                "coordination": result.score.coordination_score,
                "energy": result.score.energy_score,
                "steric": result.score.steric_score,
                "backend": result.score.backend,
                "method": result.score.method,
            },
        },
    }
    candidate_json = candidate_dir / "candidate.json"
    candidate_json.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=_json_safe) + "\n",
        encoding="utf-8",
    )
    record.output_files.update(
        {
            "candidate_json": _relative_to(candidate_json, run_dir),
            **{
                fmt: _relative_to(Path(path), run_dir)
                for fmt, path in structure_paths.items()
                if path
            },
        }
    )
    return record


def render_physchem_batch_figures(manifest: Dict[str, Any], run_dir: Path | str) -> List[Dict[str, str]]:
    output = Path(run_dir) / "figures"
    output.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np_local
    except Exception as exc:  # pragma: no cover
        return [{"kind": "error", "label": "batch figure rendering unavailable", "error": str(exc)}]

    _set_batch_figure_style(plt)
    artifacts = []
    for kind, label, renderer in (
        ("screening_funnel", "Screening funnel", _render_batch_funnel),
        ("combo_heatmap", "Dual-activity combination heatmap", _render_batch_heatmap),
        ("adsorption_volcano", "Batch adsorption-activation volcano", _render_batch_volcano),
        ("candidate_ranking", "Candidate ranking", _render_batch_rankings),
    ):
        png, svg = renderer(manifest, output, plt, np_local)
        artifacts.append(
            {
                "kind": kind,
                "label": label,
                "png": _relative_to(Path(png), Path(run_dir)),
                "svg": _relative_to(Path(svg), Path(run_dir)),
            }
        )
    return artifacts


def _render_batch_funnel(manifest: Dict[str, Any], output: Path, plt, np_local) -> Tuple[str, str]:
    records = manifest.get("records", [])
    stages = list(PIPELINE_STAGES)
    counts = [sum(1 for record in records if record.get("stage") == stage) for stage in stages]
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    fig.patch.set_facecolor("#fbfaf6")
    y = np_local.arange(len(stages))
    ax.barh(y, counts, color="#315f8d", edgecolor="white", height=0.62)
    for yi, count in zip(y, counts):
        ax.text(count + 0.05, yi, str(count), va="center", fontsize=8, color="#3a4256")
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.invert_yaxis()
    ax.set_xlabel("candidate count")
    ax.set_title("MACE + tblite Screening Funnel", loc="left")
    ax.set_xlim(0, max(counts + [1]) * 1.25)
    fig.tight_layout()
    png = output / "screening_funnel.png"
    svg = output / "screening_funnel.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_batch_heatmap(manifest: Dict[str, Any], output: Path, plt, np_local) -> Tuple[str, str]:
    records = manifest.get("records", [])
    activities = sorted({activity for record in records for activity in record.get("activities", [])})
    matrix = np_local.full((max(len(activities), 1), max(len(activities), 1)), np_local.nan)
    index = {activity: idx for idx, activity in enumerate(activities)}
    for record in records:
        acts = record.get("activities", [])
        if len(acts) < 2 or acts[0] not in index or acts[1] not in index:
            continue
        summaries = record.get("diagnostics", {}).get("activity_calculations", [])
        values = [
            summary.get("best_adsorption_energy_ev")
            for summary in summaries
            if summary.get("best_adsorption_energy_ev") is not None
        ]
        value = min(values) if values else record.get("score", {}).get("total")
        left, right = index[acts[0]], index[acts[1]]
        matrix[left, right] = value
        matrix[right, left] = value
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    fig.patch.set_facecolor("#fbfaf6")
    if activities:
        im = ax.imshow(matrix, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(activities)))
        ax.set_yticks(range(len(activities)))
        ax.set_xticklabels(activities, rotation=35, ha="right")
        ax.set_yticklabels(activities)
        ax.set_title("Dual-Activity Combination Heatmap", loc="left")
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="best adsorption or score")
    else:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No combinations available", ha="center", va="center")
    fig.tight_layout()
    png = output / "combo_activity_heatmap.png"
    svg = output / "combo_activity_heatmap.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_batch_volcano(manifest: Dict[str, Any], output: Path, plt, np_local) -> Tuple[str, str]:
    groups: Dict[str, List[Tuple[Dict, Dict, float, float]]] = {}
    for record in manifest.get("records", []):
        gate = record.get("diagnostics", {}).get("mace_tblite_gate") or {}
        if gate and not gate.get("passed", False):
            continue
        for summary in record.get("diagnostics", {}).get("activity_calculations", []):
            adsorption = summary.get("best_adsorption_energy_ev")
            activation = summary.get("activation_metric_ev")
            task_id = str(summary.get("task_id") or "")
            if task_id and adsorption is not None and activation is not None:
                groups.setdefault(task_id, []).append(
                    (record, summary, float(adsorption), float(activation))
                )
    points = [
        point
        for group in groups.values()
        if len(group) >= 5
        for point in group
    ]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    fig.patch.set_facecolor("#fbfaf6")
    if points:
        for record, summary, adsorption, activation in points:
            ax.scatter(
                [adsorption],
                [activation],
                s=52,
                color="#2a9d8f" if summary.get("mechanism_family") == "hydrolysis" else "#e76f51",
                edgecolor="white",
                linewidth=0.7,
                alpha=0.9,
            )
            ax.text(adsorption, activation, f" {summary.get('activity')}", fontsize=6.5)
        adsorption_optimum = get_screening_proxy_policy()["sabatier_adsorption_optimum_ev"]
        ax.axvline(adsorption_optimum, color="#264653", lw=1.0, ls="--", alpha=0.55)
        ax.set_xlabel("best adsorption energy (eV)")
        ax.set_ylabel("activation metric (eV)")
        ax.set_title("Adsorption-Activation Volcano", loc="left")
    else:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.55,
            "No comparable volcano cohort yet",
            ha="center",
            va="center",
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.42,
            "At least 5 converged candidates from the same reaction task are required.",
            ha="center",
            va="center",
            fontsize=8,
        )
    fig.tight_layout()
    png = output / "batch_adsorption_activation_volcano.png"
    svg = output / "batch_adsorption_activation_volcano.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_batch_rankings(manifest: Dict[str, Any], output: Path, plt, np_local) -> Tuple[str, str]:
    rows = sorted(
        [record for record in manifest.get("records", []) if record.get("status") == "passed"],
        key=_ranking_key,
        reverse=True,
    )[:12]
    labels = [row.get("combo", "")[:36] for row in rows] or ["No candidates"]
    scores = [float(row.get("score", {}).get("total") or 0.0) for row in rows] or [0.0]
    y = np_local.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, max(3.2, 0.38 * len(labels) + 1.2)))
    fig.patch.set_facecolor("#fbfaf6")
    ax.barh(y, scores, color="#c7792d", edgecolor="white", height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("assembly score")
    ax.set_title("Top Candidate Ranking", loc="left")
    ax.set_xlim(0, max(scores + [1.0]) * 1.12)
    fig.tight_layout()
    png = output / "candidate_ranking.png"
    svg = output / "candidate_ranking.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _set_batch_figure_style(plt) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.16,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

def _cores_from_legacy_pdb(
    atoms: List[Dict[str, Any]],
    bonds: Dict[int, List[int]],
    spec: DesignSpec,
) -> List[Dict[str, Any]]:
    metal_indices = [
        index
        for index, atom in enumerate(atoms)
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    metal_indices = sorted(metal_indices, key=lambda idx: atoms[idx]["coords"][0])
    cores = []
    for site_index, metal_index in enumerate(metal_indices[:len(spec.metals)]):
        metal_spec = spec.metals[site_index]
        site_id = f"M{site_index}"
        metal = atoms[metal_index]
        metal["site_id"] = site_id
        coord_indices = [
            neighbor
            for neighbor in bonds.get(metal_index, [])
            if str(atoms[neighbor].get("element", "")).upper() in {"N", "O", "S"}
        ]
        if not coord_indices:
            coord_indices = _distance_coordination_candidates(atoms, metal_index, metal_spec)
        coord_atoms = []
        for neighbor in coord_indices:
            atom = atoms[neighbor]
            atom["site_id"] = site_id
            atom["is_coord_atom"] = True
            coord_atoms.append(atom)
        cores.append(
            {
                "metal": metal,
                "coord_atoms": coord_atoms,
                "geometry": metal_spec.coordination_geometry,
                "metal_type": metal_spec.metal_type,
                "oxidation_state": metal_spec.oxidation_state,
                "site_id": site_id,
                "activity_type": metal_spec.activity_type,
            }
        )
    return cores


def _distance_coordination_candidates(
    atoms: List[Dict[str, Any]],
    metal_index: int,
    metal_spec: MetalSpec,
) -> List[int]:
    metal_pos = np.asarray(atoms[metal_index]["coords"], dtype=float)
    candidates = []
    for index, atom in enumerate(atoms):
        if index == metal_index:
            continue
        if str(atom.get("element", "")).upper() not in {"N", "O", "S"}:
            continue
        distance = float(np.linalg.norm(np.asarray(atom["coords"], dtype=float) - metal_pos))
        if distance <= 2.65:
            candidates.append((distance, index))
    return [index for _distance, index in sorted(candidates)[:metal_spec.coordination_number]]


def _new_run_dir(output_root: Path) -> Path:
    stamp = datetime.now().strftime("batch_%Y%m%d_%H%M%S")
    return output_root / stamp


def _new_manifest(run_dir: Path, templates: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "schema_version": "2.0.0",
        "knowledge_version": knowledge_version(),
        "run_dir": str(run_dir),
        "stages": list(PIPELINE_STAGES),
        "activity_count": len(templates),
        "records": [],
        "rankings": {},
        "artifacts": [],
    }


def _load_manifest(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_manifest(path: Path, manifest: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, default=_json_safe),
        encoding="utf-8",
    )


def _result_mode(result: AssemblyResult) -> str:
    return (
        result.score.details.get("bimetallic", {}).get("relation")
        or result.label
        or "single"
    )


cores_from_legacy_pdb = _cores_from_legacy_pdb
result_mode = _result_mode


def _best_adsorption_energy(payload: Dict[str, Any]) -> Optional[float]:
    energies = [
        item.get("adsorption_energy_ev")
        for item in payload.get("adsorption_candidates", [])
        if item.get("adsorption_energy_ev") is not None
    ]
    return min(energies) if energies else None


def _pareto_front(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    front = []
    for row in rows:
        if not any(_dominates(other, row) for other in rows if other is not row):
            front.append(row)
    return front


def _dominates(left: Dict[str, Any], right: Dict[str, Any]) -> bool:
    left_scores = _score_vector(left)
    right_scores = _score_vector(right)
    return all(l >= r for l, r in zip(left_scores, right_scores)) and any(
        l > r for l, r in zip(left_scores, right_scores)
    )


def _score_vector(row: Dict[str, Any]) -> Tuple[float, float, float, float]:
    score = row.get("score", {})
    return (
        float(score.get("geometry") or 0.0),
        float(score.get("coordination") or 0.0),
        float(score.get("distance") or 0.0),
        float(score.get("steric") or 0.0),
    )


def _ranking_key(row: Dict[str, Any]) -> Tuple[float, float, float, float, str]:
    return (*_score_vector(row), str(row.get("record_id")))


def _reason_code_counts(records: List[Dict[str, Any]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for record in records:
        for code in record.get("reason_codes", []):
            counts[code] = counts.get(code, 0) + 1
    return counts


def _element_from_atom_name(atom_name: str) -> str:
    letters = "".join(ch for ch in atom_name if ch.isalpha()).upper()
    if len(letters) >= 2 and letters[:2] in CATALYTIC_METAL_ELEMENTS:
        return letters[:2]
    return letters[:1] or "X"


def _slug(value: str) -> str:
    return (
        str(value)
        .replace(" + ", "__")
        .replace(":", "__")
        .replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
    )


def _relative_to(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, tuple):
        return list(value)
    if hasattr(value, "to_dict"):
        return value.to_dict()
    return str(value)
