"""Scientific interpretation helpers for completed catalysis-screening results.

The functions in this module are deliberately dependency-light. They operate on
saved JSON-compatible payloads so that completed calculations can be audited
without rerunning MACE or tblite.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, Optional, Sequence

from .physchem_knowledge import get_screening_proxy_policy
from ..utils.constants import CATALYTIC_METAL_ELEMENTS


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def derive_scan_metrics(relative_energies_ev: Sequence[Any]) -> Dict[str, Optional[float]]:
    """Return distinct range, forward-peak, and endpoint metrics.

    ``scan_energy_range_ev`` describes the full sampled energy range and must
    not be interpreted as a forward barrier. ``forward_scan_peak_ev`` is the
    largest sampled energy relative to the first point. ``reaction_energy_ev``
    is the final sampled energy relative to the first point.
    """
    values = [_finite(value) for value in relative_energies_ev]
    finite = [value for value in values if value is not None]
    start = values[0] if values else None
    end = values[-1] if values else None
    if not finite or start is None:
        return {
            "scan_energy_range_ev": None,
            "forward_scan_peak_ev": None,
            "reaction_energy_ev": None,
        }
    return {
        "scan_energy_range_ev": max(finite) - min(finite),
        "forward_scan_peak_ev": max(value - start for value in finite),
        "reaction_energy_ev": end - start if end is not None else None,
    }


def profile_from_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Select the mechanism profile that represents the requested protocol.

    A profile may be ``insufficient_sampling`` or ``incomplete`` before it has
    relative energies. Returning only profiles with energies hides those audit
    states and can mislabel a redox attempt as ``not_applicable``.
    """
    redox = payload.get("redox_state_profile") or {}
    reaction = payload.get("reaction_profile") or {}
    redox_status = redox.get("status")
    reaction_status = reaction.get("status")
    if redox.get("relative_energies_ev"):
        return redox
    if reaction.get("relative_energies_ev"):
        return reaction
    if redox_status not in {None, "not_applicable", "not_run"}:
        return redox
    if reaction_status not in {None, "not_applicable", "not_run"}:
        return reaction
    return {}


def scan_metrics_from_payload(payload: Dict[str, Any]) -> Dict[str, Optional[float]]:
    profile = profile_from_payload(payload)
    return derive_scan_metrics(profile.get("relative_energies_ev") or [])


def _expanded_substrate_records(task: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    for substrate in task.get("substrates") or []:
        for copy_index in range(max(int(substrate.get("copies") or 1), 1)):
            yield {
                "name": substrate.get("name"),
                "role": substrate.get("role"),
                "charge": int(substrate.get("charge") or 0),
                "spin": int(substrate.get("spin") or 1),
                "copy_index": copy_index,
            }


def audit_charge_microstates(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Audit whether declared charge/protonation requirements are instantiated.

    This is an interpretation audit, not an attempt to assign chemically
    correct charges automatically. It reports what the saved calculation
    actually represented and flags assumptions requiring explicit validation.
    """
    task = payload.get("task") or {}
    calculation = task.get("calculation") or {}
    activity = str(task.get("nanozyme_type") or "")
    substrates = list(_expanded_substrate_records(task))
    substrate_charge = sum(item["charge"] for item in substrates)
    profile = profile_from_payload(payload)
    electronic_cluster = payload.get("electronic_cluster") or {}
    surface_charge = _finite(electronic_cluster.get("surface_charge"))
    profile_charge = _finite(profile.get("charge"))
    complex_charge = profile_charge
    if complex_charge is None:
        for pose in payload.get("adsorption_candidates") or []:
            electronic = pose.get("electronic_state_scan") or {}
            complex_charge = _finite(electronic.get("complex_charge"))
            if complex_charge is not None:
                break

    declared_microstates = [
        str(value) for value in (calculation.get("microstates") or [])
    ]
    selected_microstate = (
        payload.get("selected_microstate")
        or profile.get("selected_microstate")
        or electronic_cluster.get("selected_microstate")
    )
    charge_state_values = profile.get("charge_states") or profile.get("charges") or []
    charge_state_scan_performed = len(
        {value for value in (_finite(item) for item in charge_state_values) if value is not None}
    ) > 1
    frame_atom_counts = {
        len(frame.get("optimized_atoms") or [])
        for frame in (profile.get("frames") or [])
        if frame.get("optimized_atoms")
    }
    proton_count_changes = len(frame_atom_counts) > 1
    coordinate_text = " ".join(
        [
            str((task.get("transition_state") or {}).get("coordinate") or ""),
            str((task.get("transition_state") or {}).get("label") or ""),
            str(calculation.get("mechanism_family") or ""),
        ]
    ).lower()
    proton_coupled = "proton" in coordinate_text

    issues = []
    if calculation.get("requires_charge") and not charge_state_scan_performed:
        issues.append("fixed_total_charge_only")
    if declared_microstates and not selected_microstate:
        issues.append("declared_microstates_not_materialized")
    if proton_coupled and not proton_count_changes:
        issues.append("proton_transfer_without_explicit_proton_bookkeeping")
    if len(substrates) > 1:
        issues.append("joint_multi_substrate_reference")

    activity_specific = {
        "Catalase": (
            "high",
            [
                "two_peroxide_sequential_adsorption_not_resolved",
                "surface_aquo_hydroxo_state_not_instantiated",
            ],
        ),
        "Glutathione Peroxidase": (
            "critical",
            [
                "gsh_protonation_state_not_enumerated",
                "large_two_substrate_pose_space_under_sampled",
            ],
        ),
        "Glucose Oxidase": (
            "moderate",
            [
                "direct_metal_proxy_not_fad_mechanism",
                "oxygen_reduction_charge_transfer_not_explicit",
            ],
        ),
        "Oxidase": (
            "moderate",
            ["oxygen_reduction_charge_transfer_not_explicit"],
        ),
        "Peroxidase": (
            "moderate",
            ["peroxide_proton_transfer_microstates_not_enumerated"],
        ),
        "Phosphatase": (
            "high",
            ["pnpp_protonation_state_not_enumerated"],
        ),
        "DNase": (
            "moderate",
            ["nucleophile_protonation_state_not_enumerated"],
        ),
    }
    risk_level, specific_issues = activity_specific.get(activity, ("moderate", []))
    issues.extend(specific_issues)

    return {
        "activity": activity,
        "condition_id": calculation.get("condition_id"),
        "ph_range": list(calculation.get("ph_range") or []),
        "requires_charge": bool(calculation.get("requires_charge")),
        "requires_spin": bool(calculation.get("requires_spin")),
        "declared_microstates": declared_microstates,
        "selected_microstate": selected_microstate,
        "microstates_materialized": bool(selected_microstate),
        "substrates": substrates,
        "substrate_charge_total": substrate_charge,
        "surface_charge": int(surface_charge) if surface_charge is not None else None,
        "complex_charge": int(complex_charge) if complex_charge is not None else None,
        "charge_state_scan_performed": charge_state_scan_performed,
        "proton_count_changes_across_scan": proton_count_changes,
        "risk_level": risk_level,
        "issues": sorted(set(issues)),
    }


def classify_evidence_tier(
    *,
    activity: str,
    method: str,
    mode: str,
    adsorption_energy_ev: Any,
    forward_scan_peak_ev: Any,
) -> Dict[str, str]:
    """Assign a conservative interpretation tier to one completed profile."""
    adsorption = _finite(adsorption_energy_ev)
    forward_peak = _finite(forward_scan_peak_ev)
    forward_peak_max = get_screening_proxy_policy()["forward_scan_peak_max_ev"]
    independent = str(mode).startswith("independent_")
    representation_role = (
        "local_monometal_cluster"
        if independent
        else "direct_bimetallic_cluster"
        if mode == "bridged"
        else "unknown_representation"
    )

    if (
        activity in {"Oxidase", "Peroxidase"}
        and method == "GFN2-xTB"
        and adsorption is not None
        and adsorption < 0.0
        and forward_peak is not None
        and forward_peak <= forward_peak_max
    ):
        return {
            "evidence_tier": "A",
            "evidence_role": "core_screening_result",
            "representation_role": representation_role,
            "tier_reason": (
                "GFN2 profile with favorable adsorption and a forward scan peak "
                "at or below the screening threshold."
            ),
        }
    if (
        activity == "Glucose Oxidase"
        and forward_peak is not None
        and forward_peak <= forward_peak_max
    ):
        return {
            "evidence_tier": "B",
            "evidence_role": "topology_sensitivity_hypothesis",
            "representation_role": representation_role,
            "tier_reason": (
                "Low forward scan peak, but adsorption/reference-state and "
                "direct-metal proxy limitations prevent core quantitative use."
            ),
        }
    if (
        mode == "bridged"
        and activity in {"Oxidase", "Peroxidase"}
        and forward_peak is not None
        and forward_peak <= forward_peak_max
    ):
        return {
            "evidence_tier": "B",
            "evidence_role": "direct_bimetallic_representation",
            "representation_role": representation_role,
            "tier_reason": (
                "Both metals are present in the electronic cluster, but matched "
                "monometal controls are still required before claiming synergy."
            ),
        }
    if independent and method == "GFN1-xTB":
        return {
            "evidence_tier": "C",
            "evidence_role": "local_monometal_baseline",
            "representation_role": representation_role,
            "tier_reason": (
                "Independent task-local cluster evaluated with the GFN1 fallback; "
                "retain as a method-stratified local-site baseline."
            ),
        }
    return {
        "evidence_tier": "D",
        "evidence_role": "diagnostic_only",
        "representation_role": representation_role,
        "tier_reason": (
            "Retain for workflow diagnosis or limitations; do not use for "
            "cross-candidate quantitative ranking."
        ),
    }


def coordination_pairs_from_atoms(
    atoms: Sequence[Dict[str, Any]],
    bond_graph: Sequence[Sequence[Any]] = (),
) -> list[tuple[int, int]]:
    """Return unique metal-donor pairs represented by the saved assembly."""
    metal_indices = [
        index
        for index, atom in enumerate(atoms)
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    metal_by_site = {
        str(atoms[index].get("site_id")): index
        for index in metal_indices
        if atoms[index].get("site_id") is not None
    }
    pairs: set[tuple[int, int]] = set()
    for index, atom in enumerate(atoms):
        if not atom.get("is_coord_atom"):
            continue
        metal_index = metal_by_site.get(str(atom.get("site_id")))
        if metal_index is not None:
            pairs.add(tuple(sorted((metal_index, index))))
    for edge in bond_graph:
        if len(edge) < 3 or str(edge[2]) != "coordinate":
            continue
        left, right = int(edge[0]), int(edge[1])
        if left in metal_indices or right in metal_indices:
            pairs.add(tuple(sorted((left, right))))
    if metal_indices:
        import numpy as np

        positions = np.asarray([atom["coords"] for atom in atoms], dtype=float)
        for index, atom in enumerate(atoms):
            if not atom.get("is_bridge_atom") or index in metal_indices:
                continue
            distances = [
                (float(np.linalg.norm(positions[index] - positions[metal_index])), metal_index)
                for metal_index in metal_indices
            ]
            distance, metal_index = min(distances)
            if distance <= 2.70:
                pairs.add(tuple(sorted((metal_index, index))))
    return sorted(pairs)


def evaluate_unconstrained_topology(
    initial_atoms: Sequence[Dict[str, Any]],
    final_atoms: Sequence[Dict[str, Any]],
    *,
    bond_graph: Sequence[Sequence[Any]] = (),
    max_bond_delta_a: float = 0.30,
    max_aligned_metal_displacement_a: float = 0.75,
    relaxation_status: Optional[str] = None,
) -> Dict[str, Any]:
    """Evaluate whether key metal-donor topology survives free relaxation."""
    import numpy as np

    if len(initial_atoms) != len(final_atoms):
        return {
            "status": "failed",
            "topology_passed": False,
            "overall_passed": False,
            "reason": "atom_count_changed",
        }
    initial_positions = np.asarray(
        [atom["coords"] for atom in initial_atoms], dtype=float
    )
    final_positions = np.asarray(
        [atom["coords"] for atom in final_atoms], dtype=float
    )
    pairs = coordination_pairs_from_atoms(initial_atoms, bond_graph)
    records = []
    for left, right in pairs:
        left_is_metal = (
            str(initial_atoms[left].get("element", "")).upper()
            in CATALYTIC_METAL_ELEMENTS
        )
        metal_index = left if left_is_metal else right
        donor_index = right if left_is_metal else left
        initial_distance = float(
            np.linalg.norm(initial_positions[left] - initial_positions[right])
        )
        final_distance = float(
            np.linalg.norm(final_positions[left] - final_positions[right])
        )
        delta = final_distance - initial_distance
        records.append(
            {
                "pair": [left, right],
                "metal_element": str(
                    initial_atoms[metal_index].get("element")
                ).upper(),
                "donor_element": str(
                    initial_atoms[donor_index].get("element")
                ).upper(),
                "initial_distance_a": initial_distance,
                "final_distance_a": final_distance,
                "delta_a": delta,
                "passed": abs(delta) <= max_bond_delta_a,
            }
        )
    bond_deltas = [record["delta_a"] for record in records]
    bond_rmsd = (
        float(np.sqrt(np.mean(np.square(bond_deltas)))) if bond_deltas else None
    )
    max_abs_bond_delta = (
        float(max(abs(value) for value in bond_deltas)) if bond_deltas else None
    )

    support_indices = [
        index
        for index, atom in enumerate(initial_atoms)
        if str(atom.get("residue_name", "")).upper() == "GRA"
        and str(atom.get("element", "")).upper() != "H"
    ]
    aligned_final = final_positions.copy()
    if len(support_indices) >= 3:
        reference = initial_positions[support_indices]
        mobile = final_positions[support_indices]
        ref_center = reference.mean(axis=0)
        mobile_center = mobile.mean(axis=0)
        covariance = (mobile - mobile_center).T @ (reference - ref_center)
        left_svd, _, right_svd = np.linalg.svd(covariance)
        rotation = left_svd @ right_svd
        if np.linalg.det(rotation) < 0:
            left_svd[:, -1] *= -1
            rotation = left_svd @ right_svd
        aligned_final = (final_positions - mobile_center) @ rotation + ref_center
    metal_indices = [
        index
        for index, atom in enumerate(initial_atoms)
        if str(atom.get("element", "")).upper()
        in CATALYTIC_METAL_ELEMENTS
    ]
    metal_displacements = [
        float(np.linalg.norm(aligned_final[index] - initial_positions[index]))
        for index in metal_indices
    ]
    max_metal_displacement = max(metal_displacements) if metal_displacements else None
    topology_passed = bool(records) and all(record["passed"] for record in records)
    if max_metal_displacement is not None:
        topology_passed = (
            topology_passed
            and max_metal_displacement <= max_aligned_metal_displacement_a
        )
    optimization_converged = relaxation_status == "converged"
    overall_passed = topology_passed and optimization_converged
    status = (
        "passed"
        if overall_passed
        else "failed"
        if not topology_passed
        else "inconclusive_not_converged"
    )
    return {
        "status": status,
        "topology_passed": topology_passed,
        "optimization_converged": optimization_converged,
        "overall_passed": overall_passed,
        "coordination_pair_count": len(records),
        "coordination_bonds": records,
        "bond_length_rmsd_a": bond_rmsd,
        "max_abs_bond_delta_a": max_abs_bond_delta,
        "max_allowed_bond_delta_a": max_bond_delta_a,
        "aligned_metal_displacements_a": metal_displacements,
        "max_aligned_metal_displacement_a": max_metal_displacement,
        "max_allowed_aligned_metal_displacement_a": (
            max_aligned_metal_displacement_a
        ),
    }


def attach_structure_validation(
    row: Dict[str, Any],
    validation: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Attach candidate-level free-relaxation readiness to an activity row."""
    output = dict(row)
    if not validation:
        output.update(
            {
                "structure_validation_status": "not_tested",
                "structure_topology_passed": None,
                "publication_readiness": "pending_structure_validation",
            }
        )
        return output
    topology = validation.get("topology_validation") or {}
    status = str(topology.get("status") or validation.get("status") or "error")
    output.update(
        {
            "structure_validation_status": status,
            "structure_topology_passed": topology.get("topology_passed"),
            "structure_bond_length_rmsd_a": topology.get("bond_length_rmsd_a"),
            "structure_max_abs_bond_delta_a": topology.get(
                "max_abs_bond_delta_a"
            ),
            "structure_max_aligned_metal_displacement_a": topology.get(
                "max_aligned_metal_displacement_a"
            ),
            "publication_readiness": (
                "eligible_for_core_interpretation"
                if topology.get("overall_passed")
                else "blocked_by_structure_validation"
                if status == "failed"
                else "pending_structure_validation"
            ),
        }
    )
    return output
