"""Substrate adsorption and transition-state screening for designed nanozymes."""
from __future__ import annotations

import copy
import math
import uuid
from dataclasses import replace
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .nanozyme_assembler import AssemblyResult
from .carbon_scaffold import passivate_graphene_edges
from .chemical_system import (
    filter_spin_multiplicities_by_electron_parity,
    infer_spin_multiplicities,
)
from .potential_evaluator import (
    PotentialEvaluationConfig,
    atoms_to_ase,
    calculator_with_harmonic_restraints,
    calculator_capabilities,
    evaluate_atoms_energy,
    get_ase_calculator,
    relaxation_plan_for_atoms,
)
from .substrate_catalog import (
    ReactionTaskSpec,
    SubstrateSpec,
    expanded_substrates,
    get_reaction_task,
)
from .scientific_audit import audit_charge_microstates, derive_scan_metrics
from .validation import validate_assembly
from .physchem_knowledge import get_screening_proxy_policy
from ..utils.constants import ALL_METAL_ELEMENTS, CATALYTIC_METAL_ELEMENTS  # PR4-1 (M12/M13)


_DIRECTIONS = np.array(
    [
        [0.0, 0.0, 1.0],
        [0.90, 0.0, 0.55],
        [-0.90, 0.0, 0.55],
        [0.0, 0.90, 0.55],
        [0.0, -0.90, 0.55],
        [0.65, 0.65, 0.62],
        [-0.65, 0.65, 0.62],
        [0.65, -0.65, 0.62],
        [-0.65, -0.65, 0.62],
    ],
    dtype=float,
)
_DIRECTIONS /= np.linalg.norm(_DIRECTIONS, axis=1)[:, None]


def screen_catalysis(
    assembly: AssemblyResult,
    task: Optional[ReactionTaskSpec] = None,
    config: Optional[PotentialEvaluationConfig] = None,
    max_adsorption_poses: int = 8,
    run_neb: bool = False,
    neb_steps: int = 80,
    neb_fmax: float = 0.08,
    run_reaction_scan: bool = False,
    reaction_scan_points: int = 5,
    reaction_scan_adsorption_window_ev: Optional[Tuple[float, float]] = None,
    progress_callback: Optional[Callable[[Dict], None]] = None,
) -> Dict:
    """Generate substrate-bound candidates and mechanism-aware screening data."""
    task = task or get_reaction_task(assembly.design_spec.nanozyme_type)
    if task is None:
        return {
            "status": "error",
            "error": f"No substrate task for {assembly.design_spec.nanozyme_type}",
        }

    base_config = config or PotentialEvaluationConfig.from_env()
    task_config, method_decision = resolve_task_config(
        task,
        replace(base_config, task=task.ml_task),
        assembly=assembly,
    )
    _emit_progress(
        progress_callback,
        stage="protocol",
        message=f"resolved {task.nanozyme_type} protocol with {task_config.backend}",
        progress=8,
        backend=task_config.backend,
        ml_task=task_config.task,
    )
    active_core = select_active_core(assembly, task)
    _emit_progress(
        progress_callback,
        stage="active_center",
        message=(
            f"selected {active_core.get('metal_type') or active_core.get('metal', {}).get('element')} "
            f"center for {task.nanozyme_type}"
        ),
        progress=14,
        site_id=active_core.get("site_id") or active_core.get("metal", {}).get("site_id"),
    )
    candidates = generate_adsorption_candidates(
        assembly,
        task,
        max_poses=max_adsorption_poses,
    )
    _emit_progress(
        progress_callback,
        stage="adsorption",
        message=f"generated {len(candidates)} substrate-bound poses",
        progress=22,
        pose_count=len(candidates),
    )

    screened = []
    for idx, candidate in enumerate(candidates):
        item = _candidate_payload(candidate)
        item.update(evaluate_adsorption_energy(candidate, assembly, task, task_config))
        screened.append(item)
        _emit_progress(
            progress_callback,
            stage="adsorption",
            message=f"evaluated adsorption pose {idx + 1}/{len(candidates)}",
            progress=22 + 42 * ((idx + 1) / max(len(candidates), 1)),
            candidate_id=item.get("candidate_id"),
            ml_status=item.get("ml_status"),
            adsorption_energy_ev=item.get("adsorption_energy_ev"),
        )

    screened.sort(key=_screen_sort_key, reverse=True)
    optimized_clusters = {}
    if task_config.backend == "tblite":
        optimized_clusters = _optimize_top_tblite_poses(
            screened,
            candidates,
            assembly,
            task,
            task_config,
            limit=2,
        )
        screened.sort(key=_screen_sort_key, reverse=True)
    best = candidates[0] if not screened else next(
        c for c in candidates if c["candidate_id"] == screened[0]["candidate_id"]
    )
    _emit_progress(
        progress_callback,
        stage="ranking",
        message="ranked adsorption poses and selected the screening geometry",
        progress=68,
        best_candidate_id=screened[0].get("candidate_id") if screened else None,
    )
    ts_plan = build_transition_state_guess(best, task)
    _emit_progress(
        progress_callback,
        stage="reaction_coordinate",
        message=(
            f"built {task.transition_state.label} coordinate"
            if ts_plan.get("status") == "ready"
            else ts_plan.get("reason", "reaction coordinate unavailable")
        ),
        progress=76,
        transition_state_status=ts_plan.get("status"),
    )
    if run_neb and ts_plan.get("status") == "ready":
        _emit_progress(
            progress_callback,
            stage="neb",
            message=f"starting NEB with {task.transition_state.images} images",
            progress=80,
        )
        ts_plan["neb"] = run_neb_barrier(
            ts_plan["initial_atoms"],
            ts_plan["final_atoms"],
            task_config,
            images=task.transition_state.images,
            steps=neb_steps,
            fmax=neb_fmax,
        )
        _emit_progress(
            progress_callback,
            stage="neb",
            message=f"NEB finished with status {ts_plan['neb'].get('status')}",
            progress=87,
            neb_status=ts_plan["neb"].get("status"),
        )
    else:
        ts_plan["neb"] = {
            "status": "not_run",
            "reason": "set run_neb=true and configure a real ML backend to calculate barriers",
        }
    scan_gate = _reaction_scan_gate(
        screened,
        reaction_scan_adsorption_window_ev,
    )
    electronic_cluster = None
    reaction_atoms = best["atoms"]
    if task_config.backend == "tblite" and scan_gate["status"] != "rejected":
        electronic_cluster = optimized_clusters.get(best["candidate_id"]) or extract_electronic_cluster(
            assembly, best, task
        )
        reaction_atoms = electronic_cluster["complex_atoms"]
        _emit_progress(
            progress_callback,
            stage="electronic_cluster",
            message=(
                f"extracted local electronic cluster "
                f"({electronic_cluster['diagnostics'].get('local_complex_atom_count')} atoms)"
            ),
            progress=84,
            diagnostics=electronic_cluster["diagnostics"],
        )
    reaction_plan = build_reaction_coordinate_plan({"atoms": reaction_atoms}, task)
    redox_plan = build_redox_activation_plan({"atoms": reaction_atoms}, task)
    if run_reaction_scan and scan_gate["status"] == "rejected":
        reaction_profile = {
            "status": "screening_rejected",
            "method": reaction_plan.get("method"),
            "reason": scan_gate["reason"],
            "validity": task.calculation.validation_level,
        }
        redox_profile = {
            "status": "screening_rejected",
            "method": redox_plan.get("method"),
            "reason": scan_gate["reason"],
            "validity": task.calculation.validation_level,
        }
        _emit_progress(
            progress_callback,
            stage="reaction_scan",
            message=f"reaction scan skipped: {scan_gate['reason']}",
            progress=96,
            reaction_profile_status="screening_rejected",
            adsorption_energy_ev=scan_gate.get("adsorption_energy_ev"),
        )
    elif run_reaction_scan and reaction_plan.get("status") == "ready":
        _emit_progress(
            progress_callback,
            stage="reaction_scan",
            message=f"running {reaction_scan_points}-point protocol-specific reaction scan",
            progress=88,
            points=reaction_scan_points,
        )
        reaction_profile = run_reaction_coordinate_scan(
            reaction_plan,
            reaction_atoms,
            assembly,
            task,
            task_config,
            points=reaction_scan_points,
            electronic_cluster=electronic_cluster,
        )
        _emit_progress(
            progress_callback,
            stage="reaction_scan",
            message=f"reaction scan finished with status {reaction_profile.get('status')}",
            progress=96,
            reaction_profile_status=reaction_profile.get("status"),
            proxy_barrier_ev=reaction_profile.get("proxy_barrier_ev"),
        )
        redox_profile = {
            "status": "not_applicable",
            "method": redox_plan.get("method"),
            "reason": "hydrolysis coordinate scan used for this mechanism",
        }
    elif run_reaction_scan and redox_plan.get("status") == "ready":
        _emit_progress(
            progress_callback,
            stage="redox_state_scan",
            message=f"running {reaction_scan_points}-point redox activation/state scan",
            progress=88,
            points=reaction_scan_points,
        )
        redox_profile = run_redox_state_scan(
            redox_plan,
            reaction_atoms,
            assembly,
            task,
            task_config,
            points=reaction_scan_points,
            electronic_cluster=electronic_cluster,
        )
        _emit_progress(
            progress_callback,
            stage="redox_state_scan",
            message=f"redox state scan finished with status {redox_profile.get('status')}",
            progress=96,
            redox_profile_status=redox_profile.get("status"),
            redox_activation_energy_ev=redox_profile.get("redox_activation_energy_ev"),
        )
        reaction_profile = {
            "status": "not_applicable",
            "method": reaction_plan.get("method") or task.calculation.barrier_method,
            "reason": "redox tasks are represented by redox_state_profile, not a generic TS path",
            "validity": task.calculation.validation_level,
        }
    elif run_reaction_scan:
        reaction_profile = {
            "status": reaction_plan.get("status", "not_ready"),
            "method": reaction_plan.get("method"),
            "reason": reaction_plan.get("reason", "coordinate plan unavailable"),
            "validity": reaction_plan.get("validity", task.calculation.validation_level),
        }
        redox_profile = {
            "status": redox_plan.get("status", "not_ready"),
            "method": redox_plan.get("method"),
            "reason": redox_plan.get("reason", "redox activation plan unavailable"),
            "validity": redox_plan.get("validity", task.calculation.validation_level),
        }
        _emit_progress(
            progress_callback,
            stage="reaction_scan",
            message=f"reaction scan {reaction_profile['status']}: {reaction_profile['reason']}",
            progress=92,
            reaction_profile_status=reaction_profile["status"],
        )
    else:
        reaction_profile = {
            "status": "not_run",
            "reason": "set run_reaction_scan=true to calculate the protocol-specific screening profile",
        }
        redox_profile = {
            "status": "not_run",
            "reason": "set run_reaction_scan=true to calculate the mechanism-specific redox state profile",
        }
        _emit_progress(
            progress_callback,
            stage="reaction_scan",
            message="reaction scan skipped by request",
            progress=92,
        )
    mechanism_visualization = build_mechanism_visualization(
        task,
        reaction_plan,
        reaction_profile,
        redox_plan,
        redox_profile,
    )
    calculation_status = _calculation_status(screened, reaction_profile, redox_profile)

    result = {
        "status": "success",
        "calculation_status": calculation_status,
        "calculation_module": "mace_tblite_v2",
        "task": task.to_dict(),
        "ml_backend": task_config.backend,
        "ml_task": task_config.task,
        "method_decision": method_decision,
        "active_center": {
            "metal_type": active_core.get("metal_type") or active_core.get("metal", {}).get("element"),
            "activity_type": active_core.get("activity_type") or task.nanozyme_type,
            "site_id": active_core.get("site_id") or active_core.get("metal", {}).get("site_id"),
            "oxidation_state": active_core.get("oxidation_state"),
        },
        "structure_electronic_state": {
            "formal_charge": assembly.formal_charge,
            "spin_multiplicities": list(assembly.spin_multiplicities),
            "warnings": list(assembly.chemistry_warnings),
        },
        "assembly_job_id": assembly.job_id,
        "assembly_label": assembly.label,
        "adsorption_candidates": screened,
        "best_adsorption_structure": {
            "candidate_id": best.get("candidate_id"),
            "components": best.get("components") or [],
            "atoms": best.get("atoms") or [],
        },
        "transition_state": _strip_atoms_from_ts_plan(ts_plan),
        "reaction_coordinate_plan": reaction_plan,
        "reaction_profile": reaction_profile,
        "redox_state_plan": redox_plan,
        "redox_state_profile": redox_profile,
        "mechanism_visualization": mechanism_visualization,
        "reaction_scan_gate": scan_gate,
        "electronic_cluster": (
            electronic_cluster["diagnostics"] if electronic_cluster else {"mode": "full_system"}
        ),
    }
    result["charge_microstate_audit"] = audit_charge_microstates(result)
    return result


def _emit_progress(
    callback: Optional[Callable[[Dict], None]],
    *,
    stage: str,
    message: str,
    progress: Optional[float] = None,
    **data,
) -> None:
    if callback is None:
        return
    event = {"stage": stage, "message": message}
    if progress is not None:
        event["progress"] = float(progress)
    event.update(data)
    try:
        callback(event)
    except Exception:
        # Progress hooks are observational; never let UI bookkeeping change the
        # scientific result path.
        pass


def _reaction_scan_gate(
    screened: Sequence[Dict],
    adsorption_window_ev: Optional[Tuple[float, float]],
) -> Dict:
    if adsorption_window_ev is None:
        return {"status": "not_configured"}
    lower, upper = sorted(float(value) for value in adsorption_window_ev)
    adsorption = (
        screened[0].get("adsorption_energy_ev")
        if screened
        else None
    )
    if adsorption is None or not math.isfinite(float(adsorption)):
        return {
            "status": "rejected",
            "adsorption_energy_ev": adsorption,
            "window_ev": [lower, upper],
            "reason": "best adsorption energy is unavailable or non-finite",
        }
    adsorption = float(adsorption)
    if adsorption < lower or adsorption > upper:
        return {
            "status": "rejected",
            "adsorption_energy_ev": adsorption,
            "window_ev": [lower, upper],
            "reason": (
                f"best adsorption energy {adsorption:.3f} eV is outside "
                f"the configured scan window [{lower:.3f}, {upper:.3f}] eV"
            ),
        }
    return {
        "status": "passed",
        "adsorption_energy_ev": adsorption,
        "window_ev": [lower, upper],
    }


def generate_adsorption_candidates(
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
    max_poses: int = 8,
) -> List[Dict]:
    surface_atoms = [copy.deepcopy(a) for a in assembly.atoms]
    metal_center = np.array(select_active_core(assembly, task)["metal"]["coords"], dtype=float)
    candidates = []

    for pose_idx in range(max(max_poses, 1)):
        atoms = [copy.deepcopy(a) for a in surface_atoms]
        components = []
        for copy_idx, substrate in enumerate(expanded_substrates(task)):
            sub_atoms, anchor_idx = build_substrate_atoms(substrate, pose_idx + copy_idx)
            distance = substrate.target_distance + 0.35 * copy_idx
            direction_idx = (pose_idx + 2 * copy_idx) % len(_DIRECTIONS)
            direction_order = np.concatenate(
                (_DIRECTIONS[direction_idx:], _DIRECTIONS[:direction_idx]), axis=0
            )
            placed = None
            placed_clearance = -float("inf")
            for direction in direction_order:
                trial = _place_substrate(
                    sub_atoms,
                    anchor_idx,
                    metal_center + direction * distance,
                    pose_idx,
                    copy_idx,
                )
                trial = _resolve_substrate_clashes(atoms, trial, direction)
                clearance = _min_substrate_surface_distance(
                    atoms + trial, len(atoms)
                )
                if clearance is not None and clearance > placed_clearance:
                    placed = trial
                    placed_clearance = clearance
                if clearance is None or clearance >= 1.55:
                    break
            placed = placed or trial
            placed = _offset_bond_indices(placed, len(atoms), f"{substrate.name}:{copy_idx}")
            atoms.extend(placed)
            components.append(
                {
                    "name": substrate.name,
                    "role": substrate.role,
                    "copy_index": copy_idx,
                    "anchor_element": placed[anchor_idx]["element"],
                    "target_distance": distance,
                }
            )

        candidates.append(
            {
                "candidate_id": str(uuid.uuid4())[:8],
                "label": f"{task.task_id} pose {pose_idx + 1}",
                "task_id": task.task_id,
                "pose_index": pose_idx,
                "atoms": atoms,
                "surface_atom_count": len(surface_atoms),
                "components": components,
                "min_substrate_surface_distance": _min_substrate_surface_distance(
                    atoms, len(surface_atoms)
                ),
            }
        )
    return candidates


def extract_electronic_cluster(
    assembly: AssemblyResult,
    candidate: Dict,
    task: ReactionTaskSpec,
    support_radius: float = 5.5,
    boundary_inner_radius: float = 4.5,
) -> Dict:
    """Extract and re-passivate the task-local cluster used by tblite/GFN2-xTB."""
    core = select_active_core(assembly, task)
    center = np.array(core["metal"]["coords"], dtype=float)
    site_id = core.get("site_id") or core["metal"].get("site_id")
    bridged = _uses_bridged_dual_metal_cluster(assembly)
    included_cores = [core]
    if bridged:
        included_cores.extend(
            candidate_core
            for candidate_core in assembly.cores
            if candidate_core is not core
        )
    included_site_ids = {
        candidate_core.get("site_id") or candidate_core["metal"].get("site_id")
        for candidate_core in included_cores
    }
    included_site_ids.discard(None)
    cluster_centers = [
        np.asarray(candidate_core["metal"]["coords"], dtype=float)
        for candidate_core in included_cores
    ]
    selected = set()
    fragment_ids = set()

    for idx, atom in enumerate(assembly.atoms):
        element = str(atom.get("element", "")).upper()
        residue = str(atom.get("residue_name", "")).upper()
        atom_site = atom.get("site_id")
        position = np.asarray(atom["coords"], dtype=float)
        distance = min(
            float(np.linalg.norm(position - cluster_center))
            for cluster_center in cluster_centers
        )
        is_support = residue in {"GRA", "NDP", "SDP"}
        if atom.get("support_parent"):
            continue
        if atom_site in included_site_ids:
            selected.add(idx)
        elif atom_site and atom_site not in included_site_ids:
            continue
        elif is_support and element != "H" and distance <= support_radius:
            selected.add(idx)
        elif not is_support and element not in CATALYTIC_METAL_ELEMENTS and distance <= support_radius + 1.2:
            selected.add(idx)
        if idx in selected and atom.get("fragment_id"):
            fragment_ids.add(atom["fragment_id"])

    for idx, atom in enumerate(assembly.atoms):
        if atom.get("fragment_id") in fragment_ids:
            selected.add(idx)

    surface_atoms = []
    for idx in sorted(selected):
        atom = copy.deepcopy(assembly.atoms[idx])
        atom["source_atom_index"] = idx
        atom["distance_to_active_metal_a"] = min(
            float(
                np.linalg.norm(
                    np.asarray(atom["coords"], dtype=float) - cluster_center
                )
            )
            for cluster_center in cluster_centers
        )
        surface_atoms.append(atom)
    surface_atoms = passivate_graphene_edges(surface_atoms)
    boundary_indices = []
    cap_indices = []
    for idx, atom in enumerate(surface_atoms):
        is_cap = bool(atom.get("support_parent"))
        distance = atom.get("distance_to_active_metal_a")
        is_boundary = is_cap or (distance is not None and float(distance) >= boundary_inner_radius)
        atom["is_cap_atom"] = is_cap
        atom["is_boundary_atom"] = is_boundary
        atom["freeze_in_electronic_cluster"] = is_boundary
        if is_cap:
            cap_indices.append(idx)
        if is_boundary:
            boundary_indices.append(idx)
    surface_charge = int(sum(int(atom.get("formal_charge", 0)) for atom in surface_atoms))

    metal_specs = []
    for included_site_id in included_site_ids:
        if (
            isinstance(included_site_id, str)
            and included_site_id.startswith("M")
            and included_site_id[1:].isdigit()
        ):
            metal_idx = int(included_site_id[1:])
            if metal_idx < len(assembly.design_spec.metals):
                metal_specs.append(assembly.design_spec.metals[metal_idx])
    if not metal_specs:
        metal_specs = [
            spec
            for spec in assembly.design_spec.metals
            if (spec.activity_type or "").lower() == task.nanozyme_type.lower()
        ]
    spin_multiplicities = infer_spin_multiplicities(metal_specs)
    spin_multiplicities = filter_spin_multiplicities_by_electron_parity(
        surface_atoms,
        surface_charge,
        spin_multiplicities,
    )

    substrate_start = int(candidate.get("surface_atom_count", len(assembly.atoms)))
    substrate_indices = list(range(substrate_start, len(candidate["atoms"])))
    old_to_new = {
        old_idx: len(surface_atoms) + position
        for position, old_idx in enumerate(substrate_indices)
    }
    substrate_atoms = []
    for old_idx in substrate_indices:
        atom = copy.deepcopy(candidate["atoms"][old_idx])
        atom["bonded_atom_indices"] = [
            old_to_new[neighbor]
            for neighbor in atom.get("bonded_atom_indices", [])
            if neighbor in old_to_new
        ]
        atom["bond_orders"] = {
            str(old_to_new[int(neighbor)]): order
            for neighbor, order in atom.get("bond_orders", {}).items()
            if int(neighbor) in old_to_new
        }
        substrate_atoms.append(atom)
    complex_atoms = [copy.deepcopy(atom) for atom in surface_atoms] + substrate_atoms

    return {
        "surface_atoms": surface_atoms,
        "complex_atoms": complex_atoms,
        "surface_charge": surface_charge,
        "spin_multiplicities": spin_multiplicities,
        "diagnostics": {
            "mode": (
                "bridged_dual_metal_cluster"
                if len(included_cores) > 1
                else "task_local_cluster"
            ),
            "site_id": site_id,
            "included_site_ids": sorted(included_site_ids),
            "included_metal_count": len(included_cores),
            "included_metals": [
                str(candidate_core.get("metal_type") or candidate_core["metal"].get("element", "")).upper()
                for candidate_core in included_cores
            ],
            "support_radius_a": support_radius,
            "boundary_inner_radius_a": boundary_inner_radius,
            "full_surface_atom_count": len(assembly.atoms),
            "local_surface_atom_count": len(surface_atoms),
            "local_complex_atom_count": len(complex_atoms),
            "surface_charge": surface_charge,
            "spin_multiplicities": spin_multiplicities,
            "boundary_atom_indices": boundary_indices,
            "cap_atom_indices": cap_indices,
            "solvent": task.calculation.solvent,
            "condition_id": task.calculation.condition_id,
            "ph_range": list(task.calculation.ph_range),
            "microstates": list(task.calculation.microstates),
            "validation_level": "screening_proxy",
        },
    }


def _uses_bridged_dual_metal_cluster(assembly: AssemblyResult) -> bool:
    mode = str(getattr(assembly.design_spec, "multi_metal_mode", "") or "").lower()
    label = str(getattr(assembly, "label", "") or "").lower()
    return (
        mode == "bridged"
        or "bridged dual-metal" in label
        or any(bool(atom.get("is_bridge_atom")) for atom in assembly.atoms)
    )


def _optimize_top_tblite_poses(
    screened: List[Dict],
    candidates: List[Dict],
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
    config: PotentialEvaluationConfig,
    *,
    limit: int,
) -> Dict[str, Dict]:
    by_id = {candidate["candidate_id"]: candidate for candidate in candidates}
    optimized_clusters: Dict[str, Dict] = {}
    for item in screened[:max(0, int(limit))]:
        if item.get("ml_status") != "success":
            continue
        candidate = by_id[item["candidate_id"]]
        cluster = extract_electronic_cluster(assembly, candidate, task)
        electronic = item.get("electronic_state_scan") or {}
        selected = electronic.get("selected_complex_state") or {}
        charge = int(electronic.get("complex_charge", 0))
        multiplicity = int(selected.get("multiplicity", 1))
        optimized = _relax_tblite_cluster(cluster, config, charge, multiplicity)
        item["local_optimization"] = optimized["diagnostics"]
        if optimized["diagnostics"].get("status") != "success":
            continue
        cluster["complex_atoms"] = optimized["atoms"]
        local_core = _local_core_for_validation(cluster, assembly, task)
        if local_core is not None:
            validation = validate_assembly(
                {"atoms": cluster["complex_atoms"], "cores": [local_core]},
                _single_core_design_spec(assembly, local_core.get("site_id")),
                stage="post_tblite",
            )
            optimized["diagnostics"]["physchem_validation"] = validation.to_dict()
            if not validation.passed:
                optimized["diagnostics"]["status"] = "rejected_post_validation"
                continue
        cluster["diagnostics"]["local_optimization"] = optimized["diagnostics"]
        optimized_clusters[item["candidate_id"]] = cluster
        multiplicities = list(electronic.get("selected_complex_state", {}).get("multiplicity") for _ in [0])
        multiplicities = [int(value) for value in multiplicities if value]
        states = _evaluate_electronic_states(
            cluster["complex_atoms"],
            config,
            charge,
            multiplicities or [multiplicity],
        )
        best_state = _lowest_successful_state(states)
        item["vertical_adsorption_energy_ev"] = item.get("adsorption_energy_ev")
        item["optimized_complex_energy_ev"] = best_state["energy_ev"]
        item["reference_inconsistent_relaxation_energy_ev"] = (
            best_state["energy_ev"]
            - float(item["surface_energy_ev"])
            - float(item["substrate_energy_sum_ev"])
        )
        item["adsorption_energy_reference"] = (
            "vertical single-point energy; optimized complex is reported separately "
            "until surface and substrate references use the same optimization protocol"
        )
        electronic["optimized_complex_states"] = states
        electronic["selected_optimized_complex_state"] = best_state
    return optimized_clusters


def _local_core_for_validation(
    cluster: Dict,
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
) -> Optional[Dict]:
    source_core = select_active_core(assembly, task)
    site_id = source_core.get("site_id") or source_core.get("metal", {}).get("site_id")
    atoms = cluster["complex_atoms"]
    metal = next(
        (
            atom for atom in atoms
            if atom.get("site_id") == site_id
            and str(atom.get("element", "")).upper()
            == str(source_core.get("metal_type") or source_core.get("metal", {}).get("element", "")).upper()
        ),
        None,
    )
    if metal is None:
        return None
    expected_names = {atom.get("atom_name") for atom in source_core.get("coord_atoms", ())}
    coord_atoms = [atom for atom in atoms if atom.get("atom_name") in expected_names]
    return {
        **source_core,
        "metal": metal,
        "coord_atoms": coord_atoms,
        "site_id": site_id,
    }


def _single_core_design_spec(assembly: AssemblyResult, site_id: Optional[str]):
    spec = copy.deepcopy(assembly.design_spec)
    if isinstance(site_id, str) and site_id.startswith("M") and site_id[1:].isdigit():
        index = int(site_id[1:])
        if index < len(spec.metals):
            spec.metals = [spec.metals[index]]
    else:
        spec.metals = spec.metals[:1]
    return spec


def _relax_tblite_cluster(
    cluster: Dict,
    config: PotentialEvaluationConfig,
    charge: int,
    multiplicity: int,
) -> Dict:
    atoms = [copy.deepcopy(atom) for atom in cluster["complex_atoms"]]
    local_config = replace(
        config,
        backend="tblite",
        model=config.model or "GFN2-xTB",
        charge=charge,
        spin=multiplicity,
        solvent=config.solvent or "water",
    )
    ase_atoms = atoms_to_ase(atoms, local_config)
    frozen = [
        index
        for index, atom in enumerate(atoms)
        if atom.get("freeze_in_electronic_cluster")
    ]
    if frozen:
        from ase.constraints import FixAtoms

        ase_atoms.set_constraint(FixAtoms(indices=frozen))
    topology_config = replace(local_config, freeze_support=False)
    topology_plan = relaxation_plan_for_atoms(atoms, topology_config)
    ase_atoms.calc = calculator_with_harmonic_restraints(
        get_ase_calculator(local_config),
        topology_plan,
        topology_config,
        restraint_scale=1.0,
    )
    try:
        initial_energy = float(ase_atoms.get_potential_energy())
        from ase.optimize import FIRE

        optimizer = FIRE(ase_atoms, logfile=None)
        converged = bool(
            optimizer.run(
                fmax=local_config.tblite_local_fmax,
                steps=local_config.tblite_local_steps,
            )
        )
        retry_optimizer = None
        if not converged:
            from ase.optimize import BFGS

            retry_optimizer = BFGS(ase_atoms, logfile=None)
            converged = bool(
                retry_optimizer.run(
                    fmax=local_config.tblite_local_fmax,
                    steps=max(60, local_config.tblite_local_steps // 2),
                )
            )
        final_energy = float(ase_atoms.get_potential_energy())
        final_forces = np.asarray(ase_atoms.get_forces(), dtype=float)
    except Exception as exc:
        return {
            "atoms": atoms,
            "diagnostics": {
                "status": "failed",
                "error": str(exc),
                "method": f"{local_config.model} constrained local optimization",
                "validation_level": "screening_proxy",
            },
        }
    initial_positions = np.asarray([atom["coords"] for atom in atoms], dtype=float)
    final_positions = np.asarray(ase_atoms.get_positions(), dtype=float)
    for atom, position in zip(atoms, final_positions):
        atom["coords"] = [float(value) for value in position]
    frozen_displacement = (
        float(np.max(np.linalg.norm(final_positions[frozen] - initial_positions[frozen], axis=1)))
        if frozen else 0.0
    )
    return {
        "atoms": atoms,
        "diagnostics": {
            "status": "success",
            "method": f"{local_config.model} constrained local optimization",
            "validation_level": "screening_proxy",
            "solvent": local_config.solvent,
            "solvation_model": local_config.tblite_solvation_model,
            "charge": charge,
            "multiplicity": multiplicity,
            "steps_run": (
                int(getattr(optimizer, "nsteps", local_config.tblite_local_steps))
                + (
                    int(
                        getattr(
                            retry_optimizer,
                            "nsteps",
                            max(60, local_config.tblite_local_steps // 2),
                        )
                    )
                    if retry_optimizer is not None
                    else 0
                )
            ),
            "retry_optimizer": "BFGS" if retry_optimizer is not None else None,
            "converged": converged,
            "initial_energy_ev": initial_energy,
            "final_energy_ev": final_energy,
            "max_force_ev_per_a": float(np.max(np.linalg.norm(final_forces, axis=1))),
            "frozen_atom_indices": frozen,
            "max_frozen_displacement_a": frozen_displacement,
            "topology_restraints": topology_plan,
        },
    }


def build_substrate_atoms(substrate: SubstrateSpec, seed_offset: int = 0) -> Tuple[List[Dict], int]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise RuntimeError("RDKit is required to build substrate structures") from exc

    mol = Chem.MolFromSmiles(substrate.smiles)
    if mol is None:
        raise ValueError(f"Invalid substrate SMILES for {substrate.name}: {substrate.smiles}")
    mol = Chem.AddHs(mol)
    seed = 61453 + int(seed_offset)
    embed_status = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True)
    if embed_status != 0:
        AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True, maxAttempts=200)
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    conf = mol.GetConformer()
    adjacency = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    bond_orders = {atom.GetIdx(): {} for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        order = float(bond.GetBondTypeAsDouble())
        adjacency[left].append(right)
        adjacency[right].append(left)
        bond_orders[left][right] = order
        bond_orders[right][left] = order
    atoms = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        element = atom.GetSymbol().upper()
        atoms.append(
            {
                "element": element,
                "coords": [float(pos.x), float(pos.y), float(pos.z)],
                "residue_name": "SUB",
                "atom_name": f"{element}{idx + 1}",
                "substrate_name": substrate.name,
                "substrate_role": substrate.role,
                "substrate_charge": substrate.charge,
                "substrate_spin": substrate.spin,
                "formal_charge": int(atom.GetFormalCharge()),
                "rdkit_atom_index": idx,
                "bonded_atom_indices": sorted(adjacency[idx]),
                "bond_orders": {str(k): v for k, v in bond_orders[idx].items()},
            }
        )
    return atoms, _pick_anchor_atom(mol, substrate)


def evaluate_adsorption_energy(
    candidate: Dict,
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
    config: PotentialEvaluationConfig,
) -> Dict:
    min_dist = candidate.get("min_substrate_surface_distance")
    payload = {
        "min_substrate_surface_distance": min_dist,
        "distance_score": _distance_score(min_dist),
    }
    if config.backend not in {"fairchem", "mace", "tblite"}:
        payload.update(
            {
                "ml_status": "unavailable",
                "adsorption_energy_ev": None,
                "warning": "set E2N_MLP_BACKEND=fairchem, mace, or tblite to compute adsorption energies",
            }
        )
        return payload

    try:
        cluster = extract_electronic_cluster(assembly, candidate, task) if config.backend == "tblite" else None
        complex_atoms = cluster["complex_atoms"] if cluster else candidate["atoms"]
        surface_atoms = cluster["surface_atoms"] if cluster else assembly.atoms
        surface_charge = cluster["surface_charge"] if cluster else assembly.formal_charge
        surface_multiplicities = (
            cluster["spin_multiplicities"] if cluster else assembly.spin_multiplicities
        )
        substrate_charge = sum(s.charge * max(s.copies, 1) for s in task.substrates)
        substrate_mults = _combined_multiplicities(
            [substrate.spin for substrate in expanded_substrates(task)]
        )
        complex_charge = surface_charge + substrate_charge
        complex_mults = _couple_multiplicity_sets(
            surface_multiplicities,
            substrate_mults,
        )
        complex_states = _evaluate_electronic_states(
            complex_atoms, config, complex_charge, complex_mults
        )
        surface_states = _evaluate_electronic_states(
            surface_atoms,
            config,
            surface_charge,
            surface_multiplicities,
        )
        complex_energy = _lowest_successful_state(complex_states)
        surface_energy = _lowest_successful_state(surface_states)
        substrate_energy_sum = 0.0
        substrate_details = []
        for idx, substrate in enumerate(expanded_substrates(task)):
            sub_atoms, _ = build_substrate_atoms(substrate, seed_offset=idx)
            sub_states = _evaluate_electronic_states(
                sub_atoms,
                config,
                substrate.charge,
                [substrate.spin],
            )
            sub_energy = _lowest_successful_state(sub_states)
            substrate_energy_sum += sub_energy["energy_ev"]
            substrate_details.append(
                {
                    "name": substrate.name,
                    "energy_ev": sub_energy["energy_ev"],
                    "charge": substrate.charge,
                    "multiplicity": substrate.spin,
                    "electronic_states": sub_states,
                }
            )
    except Exception as exc:
        payload.update(
            {
                "ml_status": "failed",
                "adsorption_energy_ev": None,
                "warning": str(exc),
            }
        )
        return payload

    adsorption_energy = (
        complex_energy["energy_ev"] - surface_energy["energy_ev"] - substrate_energy_sum
    )
    payload.update(
        {
            "ml_status": "success",
            "adsorption_energy_ev": adsorption_energy,
            "complex_energy_ev": complex_energy["energy_ev"],
            "surface_energy_ev": surface_energy["energy_ev"],
            "substrate_energy_sum_ev": substrate_energy_sum,
            "max_force_ev_per_a": complex_energy.get("max_force_ev_per_a"),
            "substrate_references": substrate_details,
            "electronic_state_scan": {
                "complex_charge": complex_charge,
                "complex_states": complex_states,
                "selected_complex_state": complex_energy,
                "surface_charge": surface_charge,
                "surface_states": surface_states,
                "selected_surface_state": surface_energy,
                "cluster_diagnostics": cluster["diagnostics"] if cluster else {"mode": "full_system"},
            },
            "warning": (
                "Use adsorption energies only within the same backend/task/reference protocol; "
                "tblite/GFN2-xTB spin ordering is a screening signal, not a final kinetic claim."
            ),
        }
    )
    return payload


def build_transition_state_guess(candidate: Dict, task: ReactionTaskSpec) -> Dict:
    ts = task.transition_state
    if not task.calculation.neb_allowed:
        reason_code = (
            "electronic_state_change"
            if task.calculation.requires_spin or task.calculation.requires_charge
            else "explicit_reaction_coordinate_required"
        )
        return {
            "status": "not_applicable",
            "label": ts.label,
            "kind": ts.kind,
            "coordinate": ts.coordinate,
            "reason_code": reason_code,
            "reason": task.calculation.rationale,
            "recommended_method": task.calculation.barrier_method,
        }
    initial_atoms = [copy.deepcopy(a) for a in candidate["atoms"]]
    final_atoms = [copy.deepcopy(a) for a in candidate["atoms"]]
    pair = _find_reactive_pair(final_atoms, ts.substrate_names, ts.reactive_bond_elements)
    if pair is None:
        return {
            "status": "not_ready",
            "label": ts.label,
            "coordinate": ts.coordinate,
            "reason": "reactive atom pair not found in generated substrate geometry",
        }

    i, j = pair
    p_i = np.array(final_atoms[i]["coords"], dtype=float)
    p_j = np.array(final_atoms[j]["coords"], dtype=float)
    vec = p_j - p_i
    norm = float(np.linalg.norm(vec))
    if norm < 1e-8:
        return {
            "status": "not_ready",
            "label": ts.label,
            "coordinate": ts.coordinate,
            "reason": "reactive atom pair has zero separation",
        }

    target = ts.final_bond_distance or (norm + 0.6)
    final_atoms[j]["coords"] = (p_i + vec / norm * target).tolist()
    return {
        "status": "ready",
        "label": ts.label,
        "kind": ts.kind,
        "coordinate": ts.coordinate,
        "initial_atoms": initial_atoms,
        "final_atoms": final_atoms,
        "reactive_pair": [
            _atom_label(initial_atoms[i], i),
            _atom_label(initial_atoms[j], j),
        ],
        "initial_distance_a": norm,
        "final_distance_a": target,
        "images": ts.images,
        "description": ts.description,
    }


def build_reaction_coordinate_plan(candidate: Dict, task: ReactionTaskSpec) -> Dict:
    """Create an explicit screening coordinate without claiming a true TS."""
    atoms = candidate.get("atoms", [])
    calculation = task.calculation
    operations: List[Dict] = []

    if calculation.barrier_method != "coordinate_scan":
        return {
            "status": "not_applicable",
            "method": calculation.barrier_method,
            "mechanism_family": calculation.mechanism_family,
            "reason": calculation.rationale or "activity requires electronic-state screening",
            "validity": calculation.validation_level,
        }

    if calculation.mechanism_family == "hydrolysis":
        broken = _find_bonded_reactive_pair(
            atoms,
            task.transition_state.substrate_names,
            task.transition_state.reactive_bond_elements,
        )
        if broken:
            operations.append(
                _coordinate_operation(
                    atoms,
                    "break",
                    broken,
                    task.transition_state.final_bond_distance or 2.3,
                    move_fragment=True,
                )
            )
            electrophile = _electrophile_index(atoms, broken)
            water_o = _nearest_substrate_atom(atoms, "H2O", "O", electrophile)
            if water_o is not None:
                target = 1.85 if atoms[electrophile]["element"] == "P" else 1.50
                operations.append(
                    _coordinate_operation(
                        atoms,
                        "form",
                        (electrophile, water_o),
                        target,
                        move_component=True,
                    )
                )
    else:
        pair = _find_bonded_reactive_pair(
            atoms,
            task.transition_state.substrate_names,
            task.transition_state.reactive_bond_elements,
        )
        if pair is None and len(task.transition_state.substrate_names) >= 2:
            left_name, right_name = task.transition_state.substrate_names[:2]
            left_element, right_element = task.transition_state.reactive_bond_elements
            pair = _find_inter_substrate_pair(
                atoms, left_name, left_element, right_name, right_element
            )
        if task.nanozyme_type == "Superoxide Dismutase":
            metal_superoxide = _nearest_metal_substrate_pair(atoms, "superoxide", "O")
            if metal_superoxide:
                operations.append(
                    _coordinate_operation(
                        atoms,
                        "bind",
                        metal_superoxide,
                        2.10,
                        move_component=True,
                    )
                )
        elif task.nanozyme_type == "Glucose Oxidase" and pair:
            operations.append(_coordinate_operation(atoms, "break", pair, 1.60))
            hydrogen = pair[1]
            oxygen = _nearest_substrate_atom(atoms, "O2", "O", hydrogen)
            if oxygen is not None:
                operations.append(
                    _coordinate_operation(atoms, "form", (oxygen, hydrogen), 1.10)
                )
        elif task.nanozyme_type == "Glutathione Peroxidase" and pair:
            peroxide = _find_bonded_reactive_pair(atoms, ("H2O2",), ("O", "O"))
            if peroxide:
                operations.append(_coordinate_operation(atoms, "break", peroxide, 2.20))
            operations.append(_coordinate_operation(atoms, "form", pair, 2.00, move_component=True))
        elif pair:
            action = "form" if task.transition_state.kind in {"proton_electron_transfer"} else "stretch"
            operations.append(
                _coordinate_operation(
                    atoms,
                    action,
                    pair,
                    task.transition_state.final_bond_distance
                    or _distance_between_indices(atoms, pair),
                    move_component=action == "form",
                )
            )

    if not operations:
        return {
            "status": "not_ready",
            "method": calculation.barrier_method,
            "reason": "protocol-specific reactive atoms could not be identified",
            "validity": calculation.validation_level,
        }
    return {
        "status": "ready",
        "method": calculation.barrier_method,
        "mechanism_family": calculation.mechanism_family,
        "operations": operations,
        "validity": calculation.validation_level,
        "warning": (
            "This is a constrained MACE/tblite screening coordinate. It is useful for "
            "candidate prioritization, not as a standalone kinetic free-energy claim."
        ),
    }


def build_redox_activation_plan(candidate: Dict, task: ReactionTaskSpec) -> Dict:
    """Create a redox-specific activation/state scan plan without claiming NEB/TS semantics."""
    atoms = candidate.get("atoms", [])
    calculation = task.calculation
    if calculation.barrier_method != "electronic_state_scan":
        return {
            "status": "not_applicable",
            "method": calculation.barrier_method,
            "mechanism_family": calculation.mechanism_family,
            "reason": "non-redox tasks use the reaction coordinate plan",
            "validity": calculation.validation_level,
        }

    operations = _redox_coordinate_operations(atoms, task)
    if not operations:
        return {
            "status": "not_ready",
            "method": calculation.barrier_method,
            "mechanism_family": calculation.mechanism_family,
            "reason": "redox activation atoms could not be identified",
            "validity": calculation.validation_level,
        }
    return {
        "status": "ready",
        "method": calculation.barrier_method,
        "mechanism_family": calculation.mechanism_family,
        "operations": operations,
        "validity": calculation.validation_level,
        "warning": (
            "This is a tblite/GFN2-xTB electronic-state activation scan. It visualizes "
            "binding/activation coordinates and spin/charge ordering, not a universal TS path."
        ),
    }


def run_reaction_coordinate_scan(
    plan: Dict,
    atoms: List[Dict],
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
    config: PotentialEvaluationConfig,
    points: int = 5,
    electronic_cluster: Optional[Dict] = None,
) -> Dict:
    if plan.get("status") != "ready":
        return {
            "status": plan.get("status", "not_ready"),
            "method": plan.get("method"),
            "reason": plan.get("reason", "coordinate plan unavailable"),
        }
    if config.backend != "tblite":
        return {
            "status": "unavailable",
            "reason": "reaction coordinate screening requires tblite / GFN2-xTB",
        }
    if int(points) < 5:
        return {
            "status": "insufficient_sampling",
            "reason": "reaction coordinate profiles require at least 5 optimized scan points",
            "requested_points": int(points),
        }

    surface_charge = (
        electronic_cluster["surface_charge"] if electronic_cluster else assembly.formal_charge
    )
    surface_multiplicities = (
        electronic_cluster["spin_multiplicities"]
        if electronic_cluster
        else assembly.spin_multiplicities
    )
    charge = surface_charge + sum(
        substrate.charge * max(substrate.copies, 1) for substrate in task.substrates
    )
    substrate_mults = _combined_multiplicities(
        [substrate.spin for substrate in expanded_substrates(task)]
    )
    complex_mults = _couple_multiplicity_sets(surface_multiplicities, substrate_mults)
    fractions = list(np.linspace(0.0, 1.0, max(points, 2)))
    is_two_dimensional = (
        task.calculation.coordinate_dimensions == 2
        and len(plan["operations"]) >= 2
    )
    frames = []
    if is_two_dimensional:
        for break_index, break_fraction in enumerate(fractions):
            for form_index, form_fraction in enumerate(fractions):
                frame = [copy.deepcopy(atom) for atom in atoms]
                _apply_coordinate_operation(frame, plan["operations"][0], float(break_fraction))
                _apply_coordinate_operation(frame, plan["operations"][1], float(form_fraction))
                for operation in plan["operations"][2:]:
                    _apply_coordinate_operation(
                        frame,
                        operation,
                        float(max(break_fraction, form_fraction)),
                    )
                frame, optimization = _optimize_coordinate_frame(
                    frame,
                    plan["operations"],
                    config,
                    charge,
                    complex_mults,
                )
                states = _evaluate_electronic_states(frame, config, charge, complex_mults)
                frames.append(
                    {
                        "grid_index": [break_index, form_index],
                        "fractions": [float(break_fraction), float(form_fraction)],
                        "coordinates": [
                            {
                                "action": op["action"],
                                "distance_a": _distance_between_indices(frame, tuple(op["pair"])),
                            }
                            for op in plan["operations"]
                        ],
                        "electronic_states": states,
                        "local_optimization": optimization,
                        "optimized_atoms": frame,
                    }
                )
    else:
        previous_frame = [copy.deepcopy(atom) for atom in atoms]
        for step, fraction in enumerate(fractions):
            frame = [copy.deepcopy(atom) for atom in previous_frame]
            for operation in plan["operations"]:
                _apply_coordinate_operation(frame, operation, float(fraction))
            frame, optimization = _optimize_coordinate_frame(
                frame,
                plan["operations"],
                config,
                charge,
                complex_mults,
            )
            states = _evaluate_electronic_states(frame, config, charge, complex_mults)
            previous_frame = frame
            frames.append(
                {
                    "step": step,
                    "fraction": float(fraction),
                    "coordinates": [
                        {
                            "action": op["action"],
                            "distance_a": _distance_between_indices(frame, tuple(op["pair"])),
                        }
                        for op in plan["operations"]
                    ],
                    "electronic_states": states,
                    "local_optimization": optimization,
                    "optimized_atoms": frame,
                }
            )

    successful = [
        state
        for frame in frames
        for state in frame["electronic_states"]
        if state.get("status") == "success"
    ]
    if not successful:
        return {"status": "failed", "frames": frames, "reason": "all electronic-state calculations failed"}
    initial = min(
        state["energy_ev"]
        for state in frames[0]["electronic_states"]
        if state.get("status") == "success"
    )
    profile = []
    for frame in frames:
        energies = [
            state["energy_ev"]
            for state in frame["electronic_states"]
            if state.get("status") == "success"
        ]
        profile.append(min(energies) - initial if energies else None)
    finite = [energy for energy in profile if energy is not None and np.isfinite(energy)]
    energy_grid = None
    proxy_barrier = max(finite) if finite else None
    if is_two_dimensional:
        grid_size = len(fractions)
        energy_grid = [profile[row * grid_size:(row + 1) * grid_size] for row in range(grid_size)]
        proxy_barrier = _minimum_monotonic_grid_barrier(energy_grid)
        profile = [energy_grid[index][index] for index in range(grid_size)]
    scan_quality = _coordinate_scan_quality(frames)
    status = "success" if scan_quality["all_frames_converged"] else "incomplete"
    metrics = derive_scan_metrics(profile)
    return {
        "status": status,
        "backend": config.backend,
        "profile_type": "hydrolysis_coordinate_path",
        "charge": charge,
        "spin_multiplicities": complex_mults,
        "relative_energies_ev": profile,
        "proxy_barrier_ev": proxy_barrier,
        **metrics,
        "coordinate_dimensions": 2 if is_two_dimensional else 1,
        "energy_grid_ev": energy_grid,
        "frames": frames,
        "scan_quality": scan_quality,
        "reason": (
            None
            if status == "success"
            else "one or more constrained tblite coordinate optimizations did not converge"
        ),
        "validity": plan.get("validity", "screening_proxy"),
    }


def run_redox_state_scan(
    plan: Dict,
    atoms: List[Dict],
    assembly: AssemblyResult,
    task: ReactionTaskSpec,
    config: PotentialEvaluationConfig,
    points: int = 5,
    electronic_cluster: Optional[Dict] = None,
) -> Dict:
    if plan.get("status") != "ready":
        return {
            "status": plan.get("status", "not_ready"),
            "method": plan.get("method"),
            "reason": plan.get("reason", "redox activation plan unavailable"),
        }
    if config.backend != "tblite":
        return {
            "status": "unavailable",
            "method": plan.get("method"),
            "reason": "redox charge/spin scans require the tblite / GFN2-xTB backend",
        }
    if int(points) < 5:
        return {
            "status": "insufficient_sampling",
            "method": plan.get("method"),
            "reason": "redox activation profiles require at least 5 optimized scan points",
            "requested_points": int(points),
        }

    surface_charge = (
        electronic_cluster["surface_charge"] if electronic_cluster else assembly.formal_charge
    )
    surface_multiplicities = (
        electronic_cluster["spin_multiplicities"]
        if electronic_cluster
        else assembly.spin_multiplicities
    )
    charge = surface_charge + sum(
        substrate.charge * max(substrate.copies, 1) for substrate in task.substrates
    )
    substrate_mults = _combined_multiplicities(
        [substrate.spin for substrate in expanded_substrates(task)]
    )
    complex_mults = _couple_multiplicity_sets(surface_multiplicities, substrate_mults)
    fractions = list(np.linspace(0.0, 1.0, max(points, 2)))
    frames = []
    previous_frame = [copy.deepcopy(atom) for atom in atoms]
    for step, fraction in enumerate(fractions):
        frame = [copy.deepcopy(atom) for atom in previous_frame]
        for operation in plan["operations"]:
            _apply_coordinate_operation(frame, operation, float(fraction))
        frame, optimization = _optimize_coordinate_frame(
            frame,
            plan["operations"],
            config,
            charge,
            complex_mults,
        )
        states = _evaluate_electronic_states(frame, config, charge, complex_mults)
        selected = _lowest_successful_state_or_none(states)
        previous_frame = frame
        frames.append(
            {
                "step": step,
                "fraction": float(fraction),
                "coordinates": [
                    {
                        "action": op["action"],
                        "distance_a": _distance_between_indices(frame, tuple(op["pair"])),
                    }
                    for op in plan["operations"]
                ],
                "electronic_states": states,
                "selected_state": selected,
                "local_optimization": optimization,
                "optimized_atoms": frame,
            }
        )

    if not frames or frames[0].get("selected_state") is None:
        return {
            "status": "failed",
            "method": plan.get("method"),
            "frames": frames,
            "reason": "initial redox electronic-state calculation failed",
        }
    initial = float(frames[0]["selected_state"]["energy_ev"])
    profile = [
        (
            float(frame["selected_state"]["energy_ev"]) - initial
            if frame.get("selected_state")
            else None
        )
        for frame in frames
    ]
    finite = [value for value in profile if value is not None and np.isfinite(value)]
    if not finite:
        return {
            "status": "failed",
            "method": plan.get("method"),
            "frames": frames,
            "reason": "all redox electronic-state calculations failed",
        }
    multiplicity_columns = sorted(
        {
            int(state["multiplicity"])
            for frame in frames
            for state in frame.get("electronic_states", [])
            if "multiplicity" in state
        }
    )
    state_matrix = []
    for frame in frames:
        row = []
        for multiplicity in multiplicity_columns:
            state = next(
                (
                    item
                    for item in frame.get("electronic_states", [])
                    if int(item.get("multiplicity", -1)) == multiplicity
                    and item.get("status") == "success"
                ),
                None,
            )
            row.append(
                float(state["energy_ev"]) - initial
                if state and state.get("energy_ev") is not None
                else None
            )
        state_matrix.append(row)
    selected_spin_path = [
        (
            int(frame["selected_state"]["multiplicity"])
            if frame.get("selected_state")
            else None
        )
        for frame in frames
    ]
    scan_quality = _coordinate_scan_quality(frames)
    status = "success" if scan_quality["all_frames_converged"] else "incomplete"
    metrics = derive_scan_metrics(profile)
    energy_span = max(finite) - min(finite)
    return {
        "status": status,
        "backend": config.backend,
        "profile_type": "redox_electronic_state_path",
        "charge": charge,
        "spin_multiplicities": complex_mults,
        "relative_energies_ev": profile,
        **metrics,
        "redox_activation_energy_ev": metrics.get("forward_scan_peak_ev"),
        "redox_energy_span_ev": energy_span,
        "redox_reaction_energy_ev": finite[-1] - finite[0] if len(finite) >= 2 else None,
        "multiplicity_columns": multiplicity_columns,
        "relative_state_energy_matrix_ev": state_matrix,
        "selected_spin_path": selected_spin_path,
        "frames": frames,
        "scan_quality": scan_quality,
        "reason": (
            None
            if status == "success"
            else "one or more constrained tblite coordinate optimizations did not converge"
        ),
        "validity": plan.get("validity", "screening_proxy"),
    }


def _optimize_coordinate_frame(
    atoms: List[Dict],
    operations: Sequence[Dict],
    config: PotentialEvaluationConfig,
    charge: int,
    multiplicities: Sequence[int],
) -> Tuple[List[Dict], Dict]:
    initial_states = _evaluate_electronic_states(atoms, config, charge, multiplicities)
    reference_state = _lowest_successful_state_or_none(initial_states)
    if reference_state is None:
        return atoms, {
            "status": "failed",
            "reason": "no electronic state was available for constrained optimization",
            "converged": False,
        }

    local_config = replace(
        config,
        backend="tblite",
        model=config.model or "GFN2-xTB",
        charge=charge,
        spin=int(reference_state["multiplicity"]),
        relax=False,
    )
    ase_atoms = atoms_to_ase(atoms, local_config)
    constraints = []
    frozen = [
        index
        for index, atom in enumerate(atoms)
        if atom.get("freeze_in_electronic_cluster")
    ]
    if frozen:
        from ase.constraints import FixAtoms

        constraints.append(FixAtoms(indices=frozen))
    pairs = sorted(
        {
            tuple(int(value) for value in operation["pair"])
            for operation in operations
            if operation.get("pair") and len(operation["pair"]) == 2
        }
    )
    if pairs:
        from ase.constraints import FixBondLengths

        constraints.append(FixBondLengths(pairs))
    if constraints:
        ase_atoms.set_constraint(constraints)
    topology_config = replace(local_config, freeze_support=False)
    topology_plan = relaxation_plan_for_atoms(atoms, topology_config)
    ase_atoms.calc = calculator_with_harmonic_restraints(
        get_ase_calculator(local_config),
        topology_plan,
        topology_config,
        restraint_scale=1.0,
    )
    try:
        initial_energy = float(ase_atoms.get_potential_energy())
        from ase.optimize import FIRE

        optimizer = FIRE(ase_atoms, logfile=None)
        converged = bool(
            optimizer.run(
                fmax=local_config.tblite_local_fmax,
                steps=local_config.tblite_local_steps,
            )
        )
        retry_optimizer = None
        if not converged:
            from ase.optimize import BFGS

            retry_optimizer = BFGS(ase_atoms, logfile=None)
            converged = bool(
                retry_optimizer.run(
                    fmax=local_config.tblite_local_fmax,
                    steps=max(60, local_config.tblite_local_steps // 2),
                )
            )
        final_energy = float(ase_atoms.get_potential_energy())
        final_forces = np.asarray(ase_atoms.get_forces(), dtype=float)
    except Exception as exc:
        return atoms, {
            "status": "failed",
            "reason": str(exc),
            "converged": False,
            "multiplicity": int(reference_state["multiplicity"]),
        }

    optimized = [copy.deepcopy(atom) for atom in atoms]
    for atom, position in zip(optimized, np.asarray(ase_atoms.get_positions(), dtype=float)):
        atom["coords"] = [float(value) for value in position]
    max_force = float(np.max(np.linalg.norm(final_forces, axis=1)))
    return optimized, {
        "status": "success",
        "method": f"{local_config.model} constrained coordinate optimization",
        "converged": converged,
        "multiplicity": int(reference_state["multiplicity"]),
        "steps_run": (
            int(getattr(optimizer, "nsteps", local_config.tblite_local_steps))
            + (
                int(
                    getattr(
                        retry_optimizer,
                        "nsteps",
                        max(60, local_config.tblite_local_steps // 2),
                    )
                )
                if retry_optimizer is not None
                else 0
            )
        ),
        "retry_optimizer": "BFGS" if retry_optimizer is not None else None,
        "fmax_target_ev_per_a": float(local_config.tblite_local_fmax),
        "max_force_ev_per_a": max_force,
        "initial_energy_ev": initial_energy,
        "final_energy_ev": final_energy,
        "fixed_coordinate_pairs": [list(pair) for pair in pairs],
        "frozen_atom_count": len(frozen),
        "topology_restraints": topology_plan,
    }


def _coordinate_scan_quality(frames: Sequence[Dict]) -> Dict[str, Any]:
    optimizations = [frame.get("local_optimization") or {} for frame in frames]
    completed = sum(1 for item in optimizations if item.get("status") == "success")
    converged = sum(
        1
        for item in optimizations
        if item.get("status") == "success" and item.get("converged") is True
    )
    return {
        "frame_count": len(frames),
        "optimized_frame_count": completed,
        "converged_frame_count": converged,
        "all_frames_optimized": bool(frames) and completed == len(frames),
        "all_frames_converged": bool(frames) and converged == len(frames),
    }


def _minimum_monotonic_grid_barrier(grid: List[List[Optional[float]]]) -> Optional[float]:
    if not grid or not grid[0]:
        return None
    rows, cols = len(grid), len(grid[0])
    costs = [[float("inf")] * cols for _ in range(rows)]
    for row in range(rows):
        for col in range(cols):
            value = grid[row][col]
            if value is None or not np.isfinite(value):
                continue
            if row == 0 and col == 0:
                costs[row][col] = float(value)
                continue
            predecessors = []
            if row:
                predecessors.append(costs[row - 1][col])
            if col:
                predecessors.append(costs[row][col - 1])
            finite_predecessors = [item for item in predecessors if np.isfinite(item)]
            if finite_predecessors:
                costs[row][col] = min(max(previous, float(value)) for previous in finite_predecessors)
    final = costs[-1][-1]
    return float(final) if np.isfinite(final) else None


optimize_coordinate_frame = _optimize_coordinate_frame
coordinate_scan_quality = _coordinate_scan_quality
minimum_monotonic_grid_barrier = _minimum_monotonic_grid_barrier


def build_mechanism_visualization(
    task: ReactionTaskSpec,
    reaction_plan: Dict,
    reaction_profile: Dict,
    redox_plan: Dict,
    redox_profile: Dict,
) -> Dict[str, Any]:
    """Return compact, mechanism-aware plotting data for reports and the web UI."""
    calculation = task.calculation
    if calculation.barrier_method == "coordinate_scan":
        status = (
            "reaction_scanned"
            if reaction_profile.get("status") == "success"
            else reaction_profile.get("status", "not_ready")
        )
        return {
            "status": status,
            "kind": "hydrolysis_coordinate_path",
            "mechanism_family": calculation.mechanism_family,
            "label": task.transition_state.label,
            "method": reaction_plan.get("method"),
            "operations": reaction_plan.get("operations", []),
            "frames": _profile_frame_summaries(reaction_profile),
            "relative_energies_ev": [
                _finite_or_none(value)
                for value in reaction_profile.get("relative_energies_ev", [])
            ],
            "activation_metric_ev": _finite_or_none(
                reaction_profile.get("proxy_barrier_ev")
            ),
            "activation_metric_label": "hydrolysis coordinate proxy barrier",
        }

    status = (
        "redox_state_scanned"
        if redox_profile.get("status") == "success"
        else redox_profile.get("status", "not_ready")
    )
    return {
        "status": status,
        "kind": "redox_electronic_state_path",
        "mechanism_family": calculation.mechanism_family,
        "label": task.transition_state.label,
        "method": redox_plan.get("method"),
        "operations": redox_plan.get("operations", []),
        "frames": _profile_frame_summaries(redox_profile),
        "relative_energies_ev": [
            _finite_or_none(value)
            for value in redox_profile.get("relative_energies_ev", [])
        ],
        "activation_metric_ev": _finite_or_none(
            redox_profile.get("redox_activation_energy_ev")
        ),
        "activation_metric_label": "redox forward scan peak",
        "multiplicity_columns": redox_profile.get("multiplicity_columns", []),
        "relative_state_energy_matrix_ev": redox_profile.get(
            "relative_state_energy_matrix_ev", []
        ),
        "selected_spin_path": redox_profile.get("selected_spin_path", []),
    }


def _profile_frame_summaries(profile: Dict) -> List[Dict]:
    frames = []
    relative = profile.get("relative_energies_ev") or []
    for idx, frame in enumerate(profile.get("frames") or []):
        selected = frame.get("selected_state") or _lowest_successful_state_or_none(
            frame.get("electronic_states") or []
        )
        frames.append(
            {
                "step": frame.get("step", idx),
                "fraction": _finite_or_none(frame.get("fraction")),
                "grid_index": frame.get("grid_index"),
                "fractions": frame.get("fractions"),
                "relative_energy_ev": (
                    _finite_or_none(relative[idx]) if idx < len(relative) else None
                ),
                "coordinates": [
                    {
                        "action": coord.get("action"),
                        "distance_a": _finite_or_none(coord.get("distance_a")),
                    }
                    for coord in frame.get("coordinates", [])
                ],
                "selected_multiplicity": (
                    selected.get("multiplicity") if selected else None
                ),
                "atoms": _compact_plot_atoms(frame.get("optimized_atoms") or []),
            }
        )
    return frames


def _compact_plot_atoms(atoms: Sequence[Dict]) -> List[Dict]:
    keys = (
        "element",
        "atom_name",
        "residue_name",
        "site_id",
        "substrate_name",
        "molecule_id",
        "is_coord_atom",
        "is_embedded_metal",
    )
    compact = []
    for atom in atoms:
        row = {key: atom.get(key) for key in keys if atom.get(key) is not None}
        row["coords"] = [float(value) for value in atom.get("coords", (0.0, 0.0, 0.0))]
        compact.append(row)
    return compact


def _calculation_status(
    screened: Sequence[Dict],
    reaction_profile: Dict,
    redox_profile: Dict,
) -> str:
    if reaction_profile.get("status") == "success":
        return "reaction_scanned"
    if redox_profile.get("status") == "success":
        return "redox_state_scanned"
    if any(item.get("ml_status") == "success" for item in screened):
        return "adsorption_screened"
    if any(item.get("ml_status") == "failed" for item in screened):
        return "failed_with_reason"
    if reaction_profile.get("status") == "not_applicable" or redox_profile.get("status") == "not_applicable":
        return "not_applicable"
    return "not_applicable"


def _finite_or_none(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def run_neb_barrier(
    initial_atoms: List[Dict],
    final_atoms: List[Dict],
    config: PotentialEvaluationConfig,
    images: int = 7,
    steps: int = 80,
    fmax: float = 0.08,
) -> Dict:
    if config.backend not in {"fairchem", "mace", "tblite"}:
        return {
            "status": "unavailable",
            "reason": "NEB requires E2N_MLP_BACKEND=fairchem, mace, or tblite",
        }
    if len(initial_atoms) != len(final_atoms):
        return {"status": "error", "reason": "initial/final atom counts differ"}

    try:
        from ase.mep import NEB
        from ase.optimize import FIRE

        initial = atoms_to_ase(initial_atoms, config)
        final = atoms_to_ase(final_atoms, config)
        ase_images = [initial]
        for _ in range(max(images - 2, 1)):
            ase_images.append(initial.copy())
        ase_images.append(final)
        neb = NEB(ase_images)
        # PR2-3 (M24 fix): IDPP (image-dependent pair potential) interpolation
        # produces a much better starting band than naive linear interpolation,
        # especially when initial and final states have very different geometries
        # (bonds making/breaking). Linear interpolation can place images inside
        # other atoms, producing unphysical configurations that FIRE cannot
        # recover from. Fall back to linear only if IDPP fails (e.g. ASE
        # version without idpp_interpolate).
        try:
            neb.idpp_interpolate(traj=None, log=None, mic=False, steps=100)
        except (AttributeError, TypeError, Exception) as _idpp_err:
            import logging as _log
            _log.getLogger("e2n.catalysis_screening").warning(
                "IDPP interpolation unavailable (%s); falling back to linear",
                _idpp_err,
            )
            neb.interpolate()
        for image in ase_images:
            image.calc = get_ase_calculator(config)
        opt = FIRE(neb, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=steps))
        energies = [float(image.get_potential_energy()) for image in ase_images]
    except Exception as exc:
        return {"status": "failed", "reason": str(exc)}

    return {
        "status": "success",
        "converged": converged,
        "energies_ev": energies,
        "barrier_ev": max(energies) - energies[0],
        "reaction_energy_ev": energies[-1] - energies[0],
        "steps": steps,
        "fmax": fmax,
    }


def _candidate_payload(candidate: Dict) -> Dict:
    return {
        "candidate_id": candidate["candidate_id"],
        "label": candidate["label"],
        "pose_index": candidate["pose_index"],
        "atom_count": len(candidate["atoms"]),
        "surface_atom_count": candidate["surface_atom_count"],
        "components": candidate["components"],
    }


def _strip_atoms_from_ts_plan(ts_plan: Dict) -> Dict:
    return {k: v for k, v in ts_plan.items() if k not in {"initial_atoms", "final_atoms"}}


def _screen_sort_key(item: Dict) -> Tuple[float, float]:
    adsorption = item.get("adsorption_energy_ev")
    if adsorption is not None and np.isfinite(adsorption):
        # Mildly exergonic adsorption is desirable; very strong binding can poison.
        policy = get_screening_proxy_policy()
        optimum = policy["sabatier_adsorption_optimum_ev"]
        width = policy["sabatier_adsorption_width_ev"]
        adsorption_score = math.exp(-((adsorption - optimum) / width) ** 2)
    else:
        adsorption_score = 0.0
    return adsorption_score, float(item.get("distance_score") or 0.0)


def _pick_anchor_atom(mol, substrate: SubstrateSpec) -> int:
    anchors = {e.upper() for e in substrate.anchor_elements}
    for atom in mol.GetAtoms():
        if atom.GetSymbol().upper() in anchors and atom.GetAtomicNum() > 1:
            return atom.GetIdx()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            return atom.GetIdx()
    return 0


def _place_substrate(
    atoms: List[Dict],
    anchor_idx: int,
    target: np.ndarray,
    pose_idx: int,
    copy_idx: int,
) -> List[Dict]:
    coords = np.array([a["coords"] for a in atoms], dtype=float)
    anchor = coords[anchor_idx]
    centered = coords - anchor
    theta = 2.0 * math.pi * ((pose_idx + 1) * (copy_idx + 1)) / 7.0
    rot = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0],
            [math.sin(theta), math.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    transformed = centered @ rot.T + target
    placed = []
    for idx, atom in enumerate(atoms):
        new_atom = copy.deepcopy(atom)
        new_atom["coords"] = transformed[idx].tolist()
        new_atom["substrate_copy"] = copy_idx
        placed.append(new_atom)
    return placed


def _resolve_substrate_clashes(
    existing_atoms: List[Dict],
    placed_atoms: List[Dict],
    direction: np.ndarray,
    min_distance: float = 1.55,
    step: float = 0.45,
    max_steps: int = 12,
) -> List[Dict]:
    """
    PR2-3 (M25 fix): if the primary direction can't clear the clash within
    max_steps, fall back to perpendicular directions (cross-products) before
    giving up. The single-direction loop used to silently return the still-
    clashing structure once max_steps was exhausted, which downstream evaluators
    then scored as a normal pose.

    Strategy:
      1. Try the original direction (most relevant for substrate normal vector).
      2. If still clashing, try two perpendicular vectors.
      3. Return whichever attempt resolved the clash, else the best (largest
         minimum distance) attempt for downstream scoring.
    """
    if not existing_atoms or not placed_atoms:
        return placed_atoms

    direction = np.array(direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm < 1e-8:
        return placed_atoms
    direction = direction / norm

    # PR2-3 (M25): build candidate directions — primary + two perpendiculars.
    # Perpendiculars are derived deterministically so retries are reproducible.
    helper = np.array([1.0, 0.0, 0.0]) if abs(direction[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp1 = np.cross(direction, helper)
    n1 = float(np.linalg.norm(perp1))
    if n1 > 1e-8:
        perp1 = perp1 / n1
    else:
        perp1 = np.array([1.0, 0.0, 0.0])
    perp2 = np.cross(direction, perp1)
    n2 = float(np.linalg.norm(perp2))
    if n2 > 1e-8:
        perp2 = perp2 / n2
    else:
        perp2 = np.array([0.0, 1.0, 0.0])

    best_atoms = placed_atoms
    best_min_dist = -float("inf")
    surface_count = len(existing_atoms)

    for try_dir in (direction, perp1, perp2):
        shifted = [copy.deepcopy(a) for a in placed_atoms]
        for _ in range(max_steps + 1):
            trial_atoms = existing_atoms + shifted
            dist = _min_substrate_surface_distance(trial_atoms, surface_count)
            if dist is None or dist >= min_distance:
                return shifted
            if dist is not None and dist > best_min_dist:
                best_min_dist = dist
                best_atoms = [copy.deepcopy(a) for a in shifted]
            for atom in shifted:
                atom["coords"] = (
                    np.array(atom["coords"], dtype=float) + try_dir * step
                ).tolist()
    return best_atoms


def select_active_core(assembly: AssemblyResult, task: ReactionTaskSpec) -> Dict:
    """Select the metal center assigned to the requested activity."""
    target = task.nanozyme_type.strip().lower()
    for core in assembly.cores:
        activity = str(core.get("activity_type") or "").strip().lower()
        if activity == target:
            return core
    if assembly.cores:
        return assembly.cores[0]
    for atom in assembly.atoms:
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS:
            return {
                "metal": atom,
                "metal_type": str(atom.get("element", "")).upper(),
                "activity_type": task.nanozyme_type,
                "coord_atoms": [],
            }
    raise ValueError(f"No catalytic metal center found for {task.nanozyme_type}")


def resolve_task_config(
    task: ReactionTaskSpec,
    config: PotentialEvaluationConfig,
    assembly: Optional[AssemblyResult] = None,
) -> Tuple[PotentialEvaluationConfig, Dict]:
    """Route reaction screening to tblite/GFN2-xTB when possible."""
    capabilities = calculator_capabilities(config)
    requires_charge = task.calculation.requires_charge or bool(
        assembly is not None and assembly.formal_charge != 0
    )
    requires_spin = task.calculation.requires_spin or bool(
        assembly is not None and any(mult != 1 for mult in assembly.spin_multiplicities)
    )
    if config.backend in {"mace", "fairchem"}:
        routed = replace(config, backend="tblite", model="GFN2-xTB", mace_head=None)
        return routed, {
            "requested_backend": config.backend,
            "selected_backend": "tblite",
            "requested_backend_capable": False,
            "reason": (
                "charge_or_spin_sensitive_task"
                if requires_charge or requires_spin
                else "reaction_module_uses_tblite_gfn2_xtb_after_structure_relaxation"
            ),
        }

    capable = (
        (not requires_charge or capabilities["charge"])
        and (not requires_spin or capabilities["spin"])
    )
    if capable:
        return config, {
            "requested_backend": config.backend,
            "selected_backend": config.backend,
            "requested_backend_capable": True,
            "reason": "requested_backend_satisfies_task",
        }

    if config.backend not in {"mace", "fairchem", "tblite"}:
        return config, {
            "requested_backend": config.backend,
            "selected_backend": config.backend,
            "requested_backend_capable": False,
            "reason": "no_electronic_backend_requested",
        }

    try:
        import tblite  # noqa: F401
    except ImportError:
        return config, {
            "requested_backend": config.backend,
            "selected_backend": config.backend,
            "requested_backend_capable": False,
            "reason": "charge_or_spin_sensitive_task_but_tblite_unavailable",
        }

    routed = replace(config, backend="tblite", model="GFN2-xTB", mace_head=None)
    return routed, {
        "requested_backend": config.backend,
        "selected_backend": "tblite",
        "requested_backend_capable": False,
        "reason": "charge_or_spin_sensitive_task",
    }


def _offset_bond_indices(atoms: List[Dict], offset: int, molecule_id: str) -> List[Dict]:
    updated = []
    for atom in atoms:
        item = copy.deepcopy(atom)
        item["molecule_id"] = molecule_id
        item["bonded_atom_indices"] = [offset + int(i) for i in item.get("bonded_atom_indices", [])]
        item["bond_orders"] = {
            str(offset + int(i)): order for i, order in item.get("bond_orders", {}).items()
        }
        updated.append(item)
    return updated


def _coordinate_operation(
    atoms: List[Dict],
    action: str,
    pair: Tuple[int, int],
    target: float,
    move_component: bool = False,
    move_fragment: bool = False,
) -> Dict:
    return {
        "action": action,
        "pair": [int(pair[0]), int(pair[1])],
        "atom_labels": [_atom_label(atoms[pair[0]], pair[0]), _atom_label(atoms[pair[1]], pair[1])],
        "initial_distance_a": _distance_between_indices(atoms, pair),
        "target_distance_a": float(target),
        "move_component": bool(move_component),
        "move_fragment": bool(move_fragment),
    }


def _distance_between_indices(atoms: Sequence[Dict], pair: Tuple[int, int]) -> float:
    return float(
        np.linalg.norm(
            np.array(atoms[pair[0]]["coords"], dtype=float)
            - np.array(atoms[pair[1]]["coords"], dtype=float)
        )
    )


def _find_bonded_reactive_pair(
    atoms: List[Dict],
    substrate_names: Tuple[str, ...],
    elements: Tuple[str, str],
) -> Optional[Tuple[int, int]]:
    names = set(substrate_names)
    elem_a, elem_b = elements[0].upper(), elements[1].upper()
    candidates = []
    for i, atom in enumerate(atoms):
        if atom.get("substrate_name") not in names:
            continue
        if str(atom.get("element", "")).upper() != elem_a:
            continue
        for j in atom.get("bonded_atom_indices", []):
            if j < 0 or j >= len(atoms):
                continue
            other = atoms[j]
            if other.get("substrate_name") not in names:
                continue
            if str(other.get("element", "")).upper() != elem_b:
                continue
            score = 0
            if {elem_a, elem_b} == {"P", "O"}:
                oxygen_idx = j if elem_b == "O" else i
                oxygen_neighbors = atoms[oxygen_idx].get("bonded_atom_indices", [])
                has_carbon_neighbor = any(
                    0 <= neighbor < len(atoms)
                    and str(atoms[neighbor].get("element", "")).upper() == "C"
                    for neighbor in oxygen_neighbors
                )
                score = 0 if has_carbon_neighbor else 10
            candidates.append((score, _distance_between_indices(atoms, (i, j)), i, j))
    if not candidates:
        return None
    _, _, left, right = min(candidates)
    return left, right


def _find_inter_substrate_pair(
    atoms: List[Dict],
    left_name: str,
    left_element: str,
    right_name: str,
    right_element: str,
) -> Optional[Tuple[int, int]]:
    left_indices = [
        i
        for i, atom in enumerate(atoms)
        if atom.get("substrate_name") == left_name
        and str(atom.get("element", "")).upper() == left_element.upper()
    ]
    right_indices = [
        i
        for i, atom in enumerate(atoms)
        if atom.get("substrate_name") == right_name
        and str(atom.get("element", "")).upper() == right_element.upper()
    ]
    pairs = [
        (_distance_between_indices(atoms, (left, right)), left, right)
        for left in left_indices
        for right in right_indices
        if left != right
    ]
    if not pairs:
        return None
    _, left, right = min(pairs)
    return left, right


def _electrophile_index(atoms: List[Dict], pair: Tuple[int, int]) -> int:
    for idx in pair:
        if str(atoms[idx].get("element", "")).upper() in {"P", "C"}:
            return idx
    return pair[0]


def _nearest_substrate_atom(
    atoms: List[Dict], substrate_name: str, element: str, reference_idx: int
) -> Optional[int]:
    candidates = [
        (_distance_between_indices(atoms, (reference_idx, i)), i)
        for i, atom in enumerate(atoms)
        if atom.get("substrate_name") == substrate_name
        and str(atom.get("element", "")).upper() == element.upper()
    ]
    return min(candidates)[1] if candidates else None


def _nearest_metal_substrate_pair(
    atoms: List[Dict], substrate_name: str, substrate_element: str
) -> Optional[Tuple[int, int]]:
    metals = [
        i
        for i, atom in enumerate(atoms)
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    substrate_atoms = [
        i
        for i, atom in enumerate(atoms)
        if atom.get("substrate_name") == substrate_name
        and str(atom.get("element", "")).upper() == substrate_element.upper()
    ]
    pairs = [
        (_distance_between_indices(atoms, (metal, substrate)), metal, substrate)
        for metal in metals
        for substrate in substrate_atoms
    ]
    if not pairs:
        return None
    _, metal, substrate = min(pairs)
    return metal, substrate


def _redox_coordinate_operations(atoms: List[Dict], task: ReactionTaskSpec) -> List[Dict]:
    operations: List[Dict] = []
    pair = _find_bonded_reactive_pair(
        atoms,
        task.transition_state.substrate_names,
        task.transition_state.reactive_bond_elements,
    )
    if pair is None and len(task.transition_state.substrate_names) >= 2:
        left_name, right_name = task.transition_state.substrate_names[:2]
        left_element, right_element = task.transition_state.reactive_bond_elements
        pair = _find_inter_substrate_pair(
            atoms, left_name, left_element, right_name, right_element
        )

    if task.nanozyme_type == "Superoxide Dismutase":
        metal_superoxide = _nearest_metal_substrate_pair(atoms, "superoxide", "O")
        if metal_superoxide:
            operations.append(
                _coordinate_operation(
                    atoms,
                    "bind",
                    metal_superoxide,
                    2.10,
                    move_component=True,
                )
            )
    elif task.nanozyme_type == "Glucose Oxidase" and pair:
        operations.append(_coordinate_operation(atoms, "stretch", pair, 1.60))
        hydrogen = pair[1]
        oxygen = _nearest_substrate_atom(atoms, "O2", "O", hydrogen)
        if oxygen is not None:
            operations.append(_coordinate_operation(atoms, "form", (oxygen, hydrogen), 1.10))
    elif task.nanozyme_type == "Glutathione Peroxidase" and pair:
        peroxide = _find_bonded_reactive_pair(atoms, ("H2O2",), ("O", "O"))
        if peroxide:
            operations.append(_coordinate_operation(atoms, "stretch", peroxide, 2.20))
        operations.append(_coordinate_operation(atoms, "bind", pair, 2.00, move_component=True))
    elif pair:
        target = (
            task.transition_state.final_bond_distance
            or _distance_between_indices(atoms, pair)
        )
        action = "bind" if task.transition_state.kind in {"proton_electron_transfer"} else "stretch"
        operations.append(
            _coordinate_operation(
                atoms,
                action,
                pair,
                target,
                move_component=action == "bind",
            )
        )
    return operations


def _apply_coordinate_operation(atoms: List[Dict], operation: Dict, fraction: float) -> None:
    left, right = operation["pair"]
    start = float(operation["initial_distance_a"])
    target = float(operation["target_distance_a"])
    desired = start + fraction * (target - start)
    anchor = np.array(atoms[left]["coords"], dtype=float)
    moving = np.array(atoms[right]["coords"], dtype=float)
    vector = moving - anchor
    norm = float(np.linalg.norm(vector))
    if norm < 1e-8:
        vector = np.array([1.0, 0.0, 0.0])
        norm = 1.0
    new_position = anchor + vector / norm * desired
    delta = new_position - moving
    moving_indices = _coordinate_moving_indices(atoms, left, right, operation)
    for idx in moving_indices:
        atoms[idx]["coords"] = (np.array(atoms[idx]["coords"], dtype=float) + delta).tolist()


def _coordinate_moving_indices(
    atoms: List[Dict], left: int, right: int, operation: Dict
) -> List[int]:
    if operation.get("move_fragment"):
        fragment = _bonded_component_after_cut(atoms, start=right, blocked_edge=(left, right))
        fragment.discard(left)
        if fragment:
            return sorted(fragment)

    molecule_id = atoms[right].get("molecule_id") if operation.get("move_component") else None
    if molecule_id:
        component = {
            idx
            for idx, atom in enumerate(atoms)
            if atom.get("molecule_id") == molecule_id
        }
        if left in component and len(component) > 1:
            component.discard(left)
        return sorted(component) if component else [right]

    return [right]


def _bonded_component_after_cut(
    atoms: List[Dict], start: int, blocked_edge: Tuple[int, int]
) -> set[int]:
    blocked = {tuple(blocked_edge), tuple(reversed(blocked_edge))}
    seen = set()
    stack = [start]
    while stack:
        idx = stack.pop()
        if idx in seen or idx < 0 or idx >= len(atoms):
            continue
        seen.add(idx)
        for neighbor in atoms[idx].get("bonded_atom_indices", []):
            if neighbor < 0 or neighbor >= len(atoms):
                continue
            if (idx, neighbor) in blocked:
                continue
            stack.append(neighbor)
    return seen


def _evaluate_electronic_states(
    atoms: List[Dict],
    config: PotentialEvaluationConfig,
    charge: int,
    multiplicities: Sequence[int],
) -> List[Dict]:
    states = []
    compatible_multiplicities = filter_spin_multiplicities_by_electron_parity(
        atoms,
        charge,
        multiplicities,
    )
    for multiplicity in compatible_multiplicities[:4]:
        try:
            result = evaluate_atoms_energy(
                atoms,
                replace(config, charge=charge, spin=multiplicity, relax=False),
            )
            states.append(
                {
                    "status": "success",
                    "multiplicity": multiplicity,
                    "energy_ev": result["energy_ev"],
                    "max_force_ev_per_a": result.get("max_force_ev_per_a"),
                }
            )
        except Exception as exc:
            states.append(
                {
                    "status": "failed",
                    "multiplicity": multiplicity,
                    "reason": str(exc),
                }
            )
    return states


coordinate_operation = _coordinate_operation
apply_coordinate_operation = _apply_coordinate_operation
evaluate_electronic_states = _evaluate_electronic_states


def _lowest_successful_state(states: Sequence[Dict]) -> Dict:
    successful = [state for state in states if state.get("status") == "success"]
    if not successful:
        reasons = "; ".join(str(state.get("reason", "failed")) for state in states)
        raise RuntimeError(f"all electronic-state calculations failed: {reasons}")
    return min(successful, key=lambda state: state["energy_ev"])


def _lowest_successful_state_or_none(states: Sequence[Dict]) -> Optional[Dict]:
    successful = [state for state in states if state.get("status") == "success"]
    if not successful:
        return None
    return min(successful, key=lambda state: state["energy_ev"])


def _combined_multiplicities(multiplicities: Sequence[int]) -> List[int]:
    totals = [1]
    for multiplicity in multiplicities:
        totals = _couple_multiplicity_sets(totals, [multiplicity])
    return totals


def _couple_multiplicity_sets(left: Sequence[int], right: Sequence[int]) -> List[int]:
    coupled = set()
    for left_mult in left:
        left_unpaired = max(int(left_mult) - 1, 0)
        for right_mult in right:
            right_unpaired = max(int(right_mult) - 1, 0)
            coupled.update(
                value + 1
                for value in range(
                    abs(left_unpaired - right_unpaired),
                    left_unpaired + right_unpaired + 1,
                    2,
                )
            )
    return sorted(coupled) or [1]


def _min_substrate_surface_distance(atoms: Sequence[Dict], surface_count: int) -> Optional[float]:
    if len(atoms) <= surface_count or surface_count == 0:
        return None
    surface = np.array([a["coords"] for a in atoms[:surface_count]], dtype=float)
    substrate = np.array([a["coords"] for a in atoms[surface_count:]], dtype=float)
    min_dist = float("inf")
    for coord in substrate:
        distances = np.linalg.norm(surface - coord, axis=1)
        min_dist = min(min_dist, float(np.min(distances)))
    return None if not np.isfinite(min_dist) else min_dist


def _distance_score(distance: Optional[float]) -> float:
    if distance is None:
        return 0.0
    # Penalize clashes below ~1.2 A and very remote poses above ~4.5 A.
    clash = 1.0 / (1.0 + math.exp(-6.0 * (distance - 1.2)))
    contact = math.exp(-((distance - 2.4) / 1.6) ** 2)
    return float(clash * contact)


def _find_reactive_pair(
    atoms: List[Dict],
    substrate_names: Tuple[str, ...],
    elements: Tuple[str, str],
) -> Optional[Tuple[int, int]]:
    if not elements[0] or not elements[1]:
        return None
    names = set(substrate_names)
    elem_a, elem_b = elements[0].upper(), elements[1].upper()
    best = None
    best_dist = float("inf")
    for i, atom_i in enumerate(atoms):
        if atom_i.get("substrate_name") not in names:
            continue
        if str(atom_i.get("element", "")).upper() != elem_a:
            continue
        for j, atom_j in enumerate(atoms):
            if i == j or atom_j.get("substrate_name") not in names:
                continue
            if str(atom_j.get("element", "")).upper() != elem_b:
                continue
            dist = float(
                np.linalg.norm(
                    np.array(atom_i["coords"], dtype=float)
                    - np.array(atom_j["coords"], dtype=float)
                )
            )
            if 0.1 < dist < best_dist:
                best = (i, j)
                best_dist = dist
    return best


def _atom_label(atom: Dict, idx: int) -> str:
    name = atom.get("substrate_name") or atom.get("residue_name") or "ATOM"
    return f"{idx}:{name}:{atom.get('atom_name', atom.get('element'))}"
