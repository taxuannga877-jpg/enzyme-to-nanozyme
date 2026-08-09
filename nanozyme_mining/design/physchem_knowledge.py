"""Versioned physicochemical constraints used by generation and screening."""
from __future__ import annotations

from dataclasses import dataclass, field
import copy
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


_DATA_PATH = Path(__file__).with_name("data") / "physchem_constraints.v1.json"
_REQUIRED_SCREENING_PROXY_KEYS = (
    "sabatier_adsorption_optimum_ev",
    "sabatier_adsorption_width_ev",
    "forward_scan_peak_max_ev",
)


@dataclass(frozen=True)
class ConstructibilityDecision:
    status: str
    prototype_ids: Tuple[str, ...] = ()
    allowed_modes: Tuple[str, ...] = ()
    reason_codes: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    evidence_source_ids: Tuple[str, ...] = ()

    @property
    def constructible(self) -> bool:
        return self.status == "constructible"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "prototype_ids": list(self.prototype_ids),
            "allowed_modes": list(self.allowed_modes),
            "reason_codes": list(self.reason_codes),
            "warnings": list(self.warnings),
            "evidence_source_ids": list(self.evidence_source_ids),
        }


@lru_cache(maxsize=1)
def load_physchem_knowledge() -> Dict[str, Any]:
    with _DATA_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def knowledge_version() -> str:
    return str(load_physchem_knowledge()["schema_version"])


def get_metal_constraint(metal: str, oxidation_state: int) -> Optional[Dict[str, Any]]:
    key = f"{str(metal).upper()}:{int(oxidation_state)}"
    value = load_physchem_knowledge()["metal_constraints"].get(key)
    return dict(value) if value else None


def get_geometry_constraint(geometry: str) -> Optional[Dict[str, Any]]:
    key = normalize_geometry(geometry)
    value = load_physchem_knowledge()["geometry_families"].get(key)
    return dict(value) if value else None


def get_activity_prototype(activity: str) -> Optional[Dict[str, Any]]:
    value = load_physchem_knowledge()["activity_prototypes"].get(str(activity).strip())
    return dict(value) if value else None


def source_record(source_id: str) -> Optional[Dict[str, Any]]:
    value = load_physchem_knowledge()["sources"].get(source_id)
    return dict(value) if value else None


def get_screening_proxy_policy() -> Dict[str, Any]:
    policy = load_physchem_knowledge().get("policy", {})
    screening = dict(policy.get("screening_proxy", {}))
    missing = [key for key in _REQUIRED_SCREENING_PROXY_KEYS if key not in screening]
    if missing:
        raise KeyError(f"physchem screening_proxy missing required keys: {', '.join(missing)}")
    for key in _REQUIRED_SCREENING_PROXY_KEYS:
        screening[key] = float(screening[key])
    return screening


def donor_distance_range(
    metal: str,
    oxidation_state: int,
    donor: str,
    declared_range: Optional[Sequence[float]] = None,
) -> Optional[Tuple[float, float]]:
    if declared_range and len(declared_range) == 2:
        return float(declared_range[0]), float(declared_range[1])
    constraint = get_metal_constraint(metal, oxidation_state)
    if not constraint:
        return None
    value = constraint.get("donor_ranges", {}).get(str(donor).upper())
    if not value:
        return None
    return float(value[0]), float(value[1])


def allowed_spin_multiplicities(metal: str, oxidation_state: int) -> Tuple[int, ...]:
    constraint = get_metal_constraint(metal, oxidation_state) or {}
    return tuple(int(value) for value in constraint.get("allowed_spins", ()))


def evaluate_constructibility(design_spec) -> ConstructibilityDecision:
    """Reject unsupported activity/metal combinations before coordinates exist."""
    activities = _activities_for_spec(design_spec)
    if not activities or not getattr(design_spec, "metals", None):
        return ConstructibilityDecision(
            status="not_constructible",
            reason_codes=("missing_activity_or_metal",),
        )

    reasons: List[str] = []
    warnings: List[str] = []
    prototypes: List[str] = []
    evidence: List[str] = []
    mode_sets: List[set[str]] = []
    metal_symbols = [str(metal.metal_type).upper() for metal in design_spec.metals]

    for index, metal_spec in enumerate(design_spec.metals):
        activity = metal_spec.activity_type or activities[min(index, len(activities) - 1)]
        prototype = get_activity_prototype(activity)
        if prototype is None:
            reasons.append(f"missing_activity_prototype:{activity}")
            continue
        prototypes.append(str(prototype["prototype_id"]))
        evidence.extend(str(item) for item in prototype.get("sources", ()))
        allowed_metals = {str(item).upper() for item in prototype.get("allowed_metals", ())}
        if metal_symbols[index] not in allowed_metals:
            reasons.append(f"unsupported_metal_for_activity:{activity}:{metal_symbols[index]}")
        mode_sets.append(set(prototype.get("allowed_modes", ())))

        if prototype.get("zinc_requires_copper_partner") and metal_symbols[index] == "ZN":
            if "CU" not in metal_symbols:
                reasons.append("sod_zinc_requires_copper_partner")
        minimum_same = int(prototype.get("minimum_same_activity_centers", 0) or 0)
        if minimum_same:
            matching = sum(
                1
                for other_index, other in enumerate(design_spec.metals)
                if str(other.metal_type).upper() in allowed_metals
                and (
                    other.activity_type == activity
                    or (other.activity_type is None and activities[min(other_index, len(activities) - 1)] == activity)
                )
            )
            if matching < minimum_same:
                reasons.append(f"insufficient_same_activity_centers:{activity}:{minimum_same}")
        minimum_cooperative = int(prototype.get("minimum_cooperative_centers", 0) or 0)
        if minimum_cooperative and len(design_spec.metals) < minimum_cooperative:
            warnings.append(f"{activity} prefers at least {minimum_cooperative} cooperative centers")
        proxy_label = prototype.get("proxy_label")
        if proxy_label:
            warnings.append(f"{activity} is reported as {proxy_label}")

        if get_metal_constraint(metal_spec.metal_type, metal_spec.oxidation_state) is None:
            reasons.append(
                f"missing_metal_oxidation_constraint:{metal_symbols[index]}:{metal_spec.oxidation_state}"
            )

    allowed_modes = set.intersection(*mode_sets) if mode_sets else set()
    if len(design_spec.metals) > 1 and not allowed_modes:
        reasons.append("no_shared_bimetallic_mode")
    requested_mode = str(getattr(design_spec, "multi_metal_mode", "") or "")
    normalized_requested = normalize_mode(requested_mode)
    if len(design_spec.metals) > 1 and normalized_requested not in {"", "independent", "all"}:
        if normalized_requested not in allowed_modes:
            reasons.append(f"unsupported_requested_mode:{normalized_requested}")

    return ConstructibilityDecision(
        status="not_constructible" if reasons else "constructible",
        prototype_ids=tuple(dict.fromkeys(prototypes)),
        allowed_modes=tuple(sorted(allowed_modes)),
        reason_codes=tuple(dict.fromkeys(reasons)),
        warnings=tuple(dict.fromkeys(warnings)),
        evidence_source_ids=tuple(dict.fromkeys(evidence)),
    )


def enrich_design_spec(design_spec):
    """Return a copy populated with versioned evidence and range defaults."""
    enriched = copy.deepcopy(design_spec)
    enriched.knowledge_version = knowledge_version()
    activities = _activities_for_spec(enriched)
    for index, metal_spec in enumerate(enriched.metals):
        activity = metal_spec.activity_type or activities[min(index, len(activities) - 1)]
        prototype = get_activity_prototype(activity) or {}
        metal_constraint = get_metal_constraint(metal_spec.metal_type, metal_spec.oxidation_state) or {}
        metal_spec.prototype_id = metal_spec.prototype_id or prototype.get("prototype_id")
        metal_spec.geometry_family = metal_spec.geometry_family or normalize_geometry(
            metal_spec.coordination_geometry
        )
        supported_geometry_cn = {
            "linear": {2},
            "trigonal_planar": {3},
            "tetrahedral": {4},
            "square_planar": {4},
            "square_pyramidal": {5},
            "trigonal_bipyramidal": {5},
            "octahedral": {6},
        }
        if metal_spec.coordination_number not in supported_geometry_cn.get(metal_spec.geometry_family, set()):
            metal_spec.geometry_family = _preferred_geometry(
                activity,
                str(metal_spec.metal_type).upper(),
                int(metal_spec.coordination_number),
            )
            metal_spec.coordination_geometry = metal_spec.geometry_family
        if (
            metal_spec.geometry_family == "tetrahedral"
            and metal_spec.coordination_number == 4
            and not any(str(atom.residue_name).upper() in {"HOH", "OH"} for atom in metal_spec.coord_atoms)
        ):
            metal_spec.geometry_family = "square_planar"
            metal_spec.coordination_geometry = "square_planar"
        metal_spec.allowed_coordination_numbers = (
            list(metal_spec.allowed_coordination_numbers)
            or [int(value) for value in metal_constraint.get("allowed_cn", ())]
        )
        metal_spec.spin_candidates = (
            list(metal_spec.spin_candidates)
            or [int(value) for value in metal_constraint.get("allowed_spins", ())]
        )
        metal_spec.condition_id = metal_spec.condition_id or prototype.get("condition_id")
        microstates = list(prototype.get("microstates", ()))
        metal_spec.microstate_id = metal_spec.microstate_id or (microstates[0] if microstates else None)
        for donor_index, donor in enumerate(metal_spec.coord_atoms):
            donor.source_id = donor.source_id or metal_constraint.get("source_id")
            donor.bond_length_range = donor.bond_length_range or donor_distance_range(
                metal_spec.metal_type,
                metal_spec.oxidation_state,
                donor.donor_element,
            )
            inferred_role = _inferred_donor_role(
                metal_spec.geometry_family,
                donor_index,
                len(metal_spec.coord_atoms),
                donor.residue_name,
            )
            if donor.role == "equatorial_network" and inferred_role != "equatorial_network":
                donor.role = inferred_role
            if donor.role == "axial_labile":
                donor.labile = True
                donor.protonation_state = donor.protonation_state or (
                    "hydroxo" if str(donor.residue_name).upper() == "OH" else "aquo"
                )
    if not enriched.condition_id:
        condition_ids = [metal.condition_id for metal in enriched.metals if metal.condition_id]
        enriched.condition_id = condition_ids[0] if len(set(condition_ids)) == 1 else "activity_specific"
    return enriched


def normalize_geometry(value: str) -> str:
    key = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "squareplanar": "square_planar",
        "squarepyramidal": "square_pyramidal",
        "trigonalbipyramidal": "trigonal_bipyramidal",
        "trigonalplanar": "trigonal_planar",
    }
    return aliases.get(key, key)


def normalize_mode(value: str) -> str:
    key = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "adjacent": "independent_adjacent",
        "cooperative": "cooperative_adjacent",
        "separated": "independent_separated",
    }
    return aliases.get(key, key)


def _activities_for_spec(design_spec) -> List[str]:
    activities = [str(value) for value in getattr(design_spec, "activities", ()) if value]
    if not activities and getattr(design_spec, "nanozyme_type", None):
        activities = [
            part.strip()
            for part in str(design_spec.nanozyme_type).split("+")
            if part.strip()
        ]
    return activities


def _inferred_donor_role(
    geometry: str,
    index: int,
    count: int,
    residue_name: str,
) -> str:
    if str(residue_name).upper() in {"HOH", "OH"}:
        return "axial_labile"
    geometry = normalize_geometry(geometry)
    equatorial_count = {
        "octahedral": 4,
        "square_pyramidal": 4,
        "trigonal_bipyramidal": 3,
    }.get(geometry, count)
    return "equatorial_network" if index < equatorial_count else "axial_labile"


def _preferred_geometry(activity: str, metal: str, coordination_number: int) -> str:
    if coordination_number == 2:
        return "linear"
    if coordination_number == 3:
        return "trigonal_planar"
    if coordination_number == 4:
        if metal in {"CU", "NI", "PD", "PT"} or activity == "Peroxidase":
            return "square_planar"
        return "tetrahedral"
    if coordination_number == 5:
        return "square_pyramidal" if metal == "CU" else "trigonal_bipyramidal"
    if coordination_number == 6:
        return "octahedral"
    return "unsupported"
