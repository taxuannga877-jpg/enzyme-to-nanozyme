"""Shared physicochemical validation for generated and relaxed nanozymes."""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .physchem_knowledge import (
    donor_distance_range,
    get_geometry_constraint,
    get_metal_constraint,
    knowledge_version,
    normalize_geometry,
)


_IDEAL_ANGLE_MULTISETS = {
    "linear": [180.0],
    "trigonal_planar": [120.0, 120.0, 120.0],
    "tetrahedral": [109.47] * 6,
    "square_planar": [90.0] * 4 + [180.0] * 2,
    "square_pyramidal": [90.0] * 8 + [180.0] * 2,
    "trigonal_bipyramidal": [90.0] * 6 + [120.0] * 3 + [180.0],
    "octahedral": [90.0] * 12 + [180.0] * 3,
}

_COVALENT_RADII_A = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "FE": 1.24,
    "CO": 1.18,
    "NI": 1.17,
    "CU": 1.17,
    "ZN": 1.25,
    "MN": 1.39,
    "MO": 1.54,
    "W": 1.62,
    "V": 1.53,
}


def coordination_angle_rms_deg(
    geometry: str,
    metal_pos: np.ndarray,
    coord_atoms: Sequence[Dict[str, Any]],
) -> Optional[float]:
    """Return angular RMS against the shared ideal angle multiset, or None."""
    normalized = normalize_geometry(geometry)
    ideal = _IDEAL_ANGLE_MULTISETS.get(normalized)
    if ideal is None:
        return None

    vectors = []
    for donor in coord_atoms:
        vector = np.asarray(donor["coords"], dtype=float) - metal_pos
        norm = float(np.linalg.norm(vector))
        if norm > 1e-8:
            vectors.append(vector / norm)

    angles = sorted(
        math.degrees(math.acos(float(np.clip(np.dot(vectors[i], vectors[j]), -1.0, 1.0))))
        for i in range(len(vectors))
        for j in range(i + 1, len(vectors))
    )
    if len(angles) != len(ideal):
        return None

    ideal_sorted = sorted(ideal)
    return float(
        math.sqrt(sum((actual - target) ** 2 for actual, target in zip(angles, ideal_sorted)) / len(angles))
    )


@dataclass
class ValidationIssue:
    code: str
    message: str
    severity: str = "error"
    site_id: Optional[str] = None
    atom_indices: List[int] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "code": self.code,
            "message": self.message,
            "severity": self.severity,
            "site_id": self.site_id,
            "atom_indices": list(self.atom_indices),
            "details": dict(self.details),
        }


@dataclass
class CenterValidation:
    site_id: str
    metal: str
    oxidation_state: int
    expected_cn: int
    actual_cn: int
    geometry: str
    bond_diagnostics: List[Dict[str, Any]] = field(default_factory=list)
    angle_rms_deg: Optional[float] = None
    planarity_rms_a: Optional[float] = None
    issues: List[ValidationIssue] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "site_id": self.site_id,
            "metal": self.metal,
            "oxidation_state": self.oxidation_state,
            "expected_cn": self.expected_cn,
            "actual_cn": self.actual_cn,
            "geometry": self.geometry,
            "bond_diagnostics": list(self.bond_diagnostics),
            "angle_rms_deg": self.angle_rms_deg,
            "planarity_rms_a": self.planarity_rms_a,
            "passed": self.passed,
            "issues": [issue.to_dict() for issue in self.issues],
        }


@dataclass
class ValidationReport:
    stage: str
    centers: List[CenterValidation] = field(default_factory=list)
    issues: List[ValidationIssue] = field(default_factory=list)
    knowledge_version: str = field(default_factory=knowledge_version)

    @property
    def passed(self) -> bool:
        return not any(issue.severity == "error" for issue in self.all_issues)

    @property
    def all_issues(self) -> List[ValidationIssue]:
        return self.issues + [issue for center in self.centers for issue in center.issues]

    @property
    def reason_codes(self) -> List[str]:
        return list(dict.fromkeys(issue.code for issue in self.all_issues))

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "passed": self.passed,
            "knowledge_version": self.knowledge_version,
            "reason_codes": self.reason_codes,
            "centers": [center.to_dict() for center in self.centers],
            "issues": [issue.to_dict() for issue in self.issues],
        }


def validate_assembly(
    assembly: Dict[str, Any],
    design_spec,
    *,
    stage: str = "preflight",
    formal_charge: Optional[int] = None,
    spin_multiplicities: Optional[Sequence[int]] = None,
) -> ValidationReport:
    report = ValidationReport(stage=stage)
    atoms = list(assembly.get("atoms") or [])
    cores = list(assembly.get("cores") or [])
    if not cores:
        report.issues.append(ValidationIssue("missing_metal_core", "No metal cores were generated."))
        return report

    for index, core in enumerate(cores):
        metal_spec = design_spec.metals[index] if index < len(design_spec.metals) else None
        report.centers.append(_validate_center(core, atoms, metal_spec))

    _validate_bridge_accounting(report, cores)
    _validate_severe_clashes(report, atoms)
    if formal_charge is not None and spin_multiplicities:
        _validate_electron_parity(report, atoms, int(formal_charge), spin_multiplicities)
    return report


def _validate_center(core: Dict[str, Any], atoms: List[Dict[str, Any]], metal_spec) -> CenterValidation:
    metal = str(core.get("metal_type") or core.get("metal", {}).get("element") or "").upper()
    oxidation = int(core.get("oxidation_state", getattr(metal_spec, "oxidation_state", 0) or 0))
    site_id = str(core.get("site_id") or core.get("metal", {}).get("site_id") or "unknown")
    geometry = normalize_geometry(
        core.get("geometry") or getattr(metal_spec, "geometry_family", None)
        or getattr(metal_spec, "coordination_geometry", "")
    )
    expected_cn = int(getattr(metal_spec, "coordination_number", len(core.get("coord_atoms", ()))) or 0)
    coord_atoms = list(core.get("coord_atoms") or [])
    center = CenterValidation(
        site_id=site_id,
        metal=metal,
        oxidation_state=oxidation,
        expected_cn=expected_cn,
        actual_cn=len(coord_atoms),
        geometry=geometry,
    )

    if not coord_atoms:
        center.issues.append(
            ValidationIssue("coordination_number_zero", f"{site_id} has no declared donor atoms.", site_id=site_id)
        )
        return center
    if expected_cn and len(coord_atoms) != expected_cn:
        center.issues.append(
            ValidationIssue(
                "coordination_number_mismatch",
                f"{site_id} expected CN={expected_cn}, found CN={len(coord_atoms)}.",
                site_id=site_id,
                details={"expected": expected_cn, "actual": len(coord_atoms)},
            )
        )

    metal_constraint = get_metal_constraint(metal, oxidation)
    if metal_constraint is None:
        center.issues.append(
            ValidationIssue(
                "missing_metal_oxidation_constraint",
                f"No supported physicochemical constraint for {metal}({oxidation}+).",
                site_id=site_id,
            )
        )
    elif len(coord_atoms) not in set(int(value) for value in metal_constraint.get("allowed_cn", ())):
        center.issues.append(
            ValidationIssue(
                "coordination_number_outside_supported_set",
                f"{site_id} CN={len(coord_atoms)} is unsupported for {metal}({oxidation}+).",
                site_id=site_id,
                details={"allowed": metal_constraint.get("allowed_cn", [])},
            )
        )

    metal_pos = np.asarray(core["metal"]["coords"], dtype=float)
    seen = set()
    spec_atoms = list(getattr(metal_spec, "coord_atoms", ()) or ())
    for donor_index, donor in enumerate(coord_atoms):
        identity = _atom_identity(donor)
        if identity in seen:
            center.issues.append(
                ValidationIssue(
                    "duplicate_coordination_atom",
                    f"{site_id} counts donor {identity} more than once.",
                    site_id=site_id,
                )
            )
        seen.add(identity)
        donor_element = str(donor.get("donor_element") or donor.get("element") or "").upper()
        declared = spec_atoms[donor_index] if donor_index < len(spec_atoms) else None
        declared_range = getattr(declared, "bond_length_range", None)
        allowed = donor_distance_range(metal, oxidation, donor_element, declared_range)
        distance = float(np.linalg.norm(np.asarray(donor["coords"], dtype=float) - metal_pos))
        diagnostic = {
            "donor_index": donor_index,
            "donor_element": donor_element,
            "atom_name": donor.get("atom_name"),
            "role": donor.get("coordination_role") or getattr(declared, "role", None),
            "distance_a": distance,
            "allowed_range_a": list(allowed) if allowed else None,
            "source_id": donor.get("source_id") or getattr(declared, "source_id", None)
            or (metal_constraint or {}).get("source_id"),
        }
        center.bond_diagnostics.append(diagnostic)
        if allowed is None:
            center.issues.append(
                ValidationIssue(
                    "missing_bond_range",
                    f"No {metal}-{donor_element} distance range for {site_id}.",
                    site_id=site_id,
                )
            )
        elif distance < allowed[0]:
            center.issues.append(
                ValidationIssue(
                    "coordination_bond_too_short",
                    f"{site_id} {metal}-{donor_element} distance {distance:.3f} A is below {allowed[0]:.3f} A.",
                    site_id=site_id,
                    details=diagnostic,
                )
            )
        elif distance > allowed[1]:
            center.issues.append(
                ValidationIssue(
                    "coordination_bond_too_long",
                    f"{site_id} {metal}-{donor_element} distance {distance:.3f} A exceeds {allowed[1]:.3f} A.",
                    site_id=site_id,
                    details=diagnostic,
                )
            )

    _validate_geometry(center, metal_pos, coord_atoms)
    _validate_unexpected_donors(center, metal_pos, atoms, coord_atoms)
    return center


def _validate_geometry(
    center: CenterValidation,
    metal_pos: np.ndarray,
    coord_atoms: List[Dict[str, Any]],
) -> None:
    geometry_constraint = get_geometry_constraint(center.geometry)
    ideal = _IDEAL_ANGLE_MULTISETS.get(center.geometry)
    if geometry_constraint is None or ideal is None:
        center.issues.append(
            ValidationIssue(
                "unsupported_geometry_family",
                f"No shape constraint for geometry {center.geometry!r}.",
                site_id=center.site_id,
            )
        )
        return
    angle_rms = coordination_angle_rms_deg(center.geometry, metal_pos, coord_atoms)
    if angle_rms is None:
        vector_count = 0
        for donor in coord_atoms:
            vector = np.asarray(donor["coords"], dtype=float) - metal_pos
            if float(np.linalg.norm(vector)) > 1e-8:
                vector_count += 1
        center.issues.append(
            ValidationIssue(
                "geometry_coordination_mismatch",
                f"{center.site_id} geometry {center.geometry} is incompatible with CN={vector_count}.",
                site_id=center.site_id,
            )
        )
        return
    center.angle_rms_deg = angle_rms
    threshold = float(geometry_constraint.get("shape_rms_max_deg", 30.0))
    if center.angle_rms_deg > threshold:
        center.issues.append(
            ValidationIssue(
                "coordination_geometry_distorted",
                f"{center.site_id} angular RMS {center.angle_rms_deg:.1f} deg exceeds {threshold:.1f} deg.",
                site_id=center.site_id,
                details={"angle_rms_deg": center.angle_rms_deg, "threshold_deg": threshold},
            )
        )
    if "planarity_rms_max_a" in geometry_constraint and len(coord_atoms) >= 3:
        points = np.asarray([donor["coords"] for donor in coord_atoms], dtype=float)
        centroid = points.mean(axis=0)
        _, _, vh = np.linalg.svd(points - centroid, full_matrices=False)
        normal = vh[-1]
        donor_deviation = np.dot(points - centroid, normal)
        metal_deviation = float(np.dot(metal_pos - centroid, normal))
        deviations = np.append(donor_deviation, metal_deviation)
        center.planarity_rms_a = float(math.sqrt(float(np.mean(deviations ** 2))))
        threshold_a = float(geometry_constraint["planarity_rms_max_a"])
        if center.planarity_rms_a > threshold_a:
            center.issues.append(
                ValidationIssue(
                    "coordination_plane_distorted",
                    f"{center.site_id} plane RMS {center.planarity_rms_a:.3f} A exceeds {threshold_a:.3f} A.",
                    site_id=center.site_id,
                    details={"planarity_rms_a": center.planarity_rms_a, "threshold_a": threshold_a},
                )
            )


def _validate_unexpected_donors(
    center: CenterValidation,
    metal_pos: np.ndarray,
    atoms: List[Dict[str, Any]],
    declared: List[Dict[str, Any]],
) -> None:
    declared_ids = {_atom_identity(atom) for atom in declared}
    for index, atom in enumerate(atoms):
        element = str(atom.get("element") or "").upper()
        if element not in {"N", "O", "S"} or _atom_identity(atom) in declared_ids:
            continue
        if atom.get("substrate_name") or str(atom.get("residue_name", "")).upper() == "SUB":
            continue
        if atom.get("site_id") not in {None, center.site_id} and not atom.get("is_bridge_atom"):
            continue
        allowed = donor_distance_range(center.metal, center.oxidation_state, element)
        if allowed is None:
            continue
        distance = float(np.linalg.norm(np.asarray(atom["coords"], dtype=float) - metal_pos))
        if distance <= allowed[1] + 0.10:
            center.issues.append(
                ValidationIssue(
                    "unexpected_coordination_atom",
                    f"{center.site_id} has undeclared {element} donor at {distance:.3f} A.",
                    site_id=center.site_id,
                    atom_indices=[index],
                    details={"distance_a": distance},
                )
            )


def _validate_bridge_accounting(report: ValidationReport, cores: List[Dict[str, Any]]) -> None:
    bridge_counts: Dict[Tuple[Any, ...], List[str]] = {}
    for core in cores:
        site_id = str(core.get("site_id") or "unknown")
        for atom in core.get("coord_atoms", ()): 
            if atom.get("is_bridge_atom"):
                bridge_counts.setdefault(_atom_identity(atom), []).append(site_id)
    for identity, sites in bridge_counts.items():
        if len(sites) > 2 or len(sites) != len(set(sites)):
            report.issues.append(
                ValidationIssue(
                    "invalid_bridge_slot_accounting",
                    f"Bridge donor {identity} is counted by invalid centers {sites}.",
                    details={"sites": sites},
                )
            )


def _validate_severe_clashes(report: ValidationReport, atoms: List[Dict[str, Any]]) -> None:
    if len(atoms) < 2:
        return
    positions = np.asarray([atom["coords"] for atom in atoms], dtype=float)
    for left in range(len(atoms)):
        for right in range(left + 1, len(atoms)):
            distance = float(np.linalg.norm(positions[right] - positions[left]))
            left_element = str(atoms[left].get("element", "")).upper()
            right_element = str(atoms[right].get("element", "")).upper()
            radius_sum = (
                _COVALENT_RADII_A.get(left_element, 0.80)
                + _COVALENT_RADII_A.get(right_element, 0.80)
            )
            threshold = max(0.50, 0.70 * radius_sum)
            if distance < threshold:
                report.issues.append(
                    ValidationIssue(
                        "severe_atomic_clash",
                        (
                            f"Atoms {left} ({left_element}) and {right} ({right_element}) "
                            f"are only {distance:.3f} A apart; minimum clearance is "
                            f"{threshold:.3f} A."
                        ),
                        atom_indices=[left, right],
                        details={
                            "distance_a": distance,
                            "minimum_clearance_a": threshold,
                            "elements": [left_element, right_element],
                        },
                    )
                )


def _validate_electron_parity(
    report: ValidationReport,
    atoms: List[Dict[str, Any]],
    formal_charge: int,
    spin_multiplicities: Sequence[int],
) -> None:
    try:
        from ase.data import atomic_numbers
    except ImportError:
        return
    electron_count = sum(atomic_numbers.get(str(atom.get("element", "")).title(), 0) for atom in atoms)
    electron_count -= formal_charge
    compatible = [
        int(multiplicity)
        for multiplicity in spin_multiplicities
        if (electron_count - (int(multiplicity) - 1)) % 2 == 0
    ]
    if not compatible:
        report.issues.append(
            ValidationIssue(
                "electron_spin_parity_mismatch",
                f"Charge {formal_charge} and multiplicities {list(spin_multiplicities)} conflict with electron parity.",
                details={"electron_count": electron_count},
            )
        )


def _atom_identity(atom: Dict[str, Any]) -> Tuple[Any, ...]:
    coords = tuple(round(float(value), 5) for value in atom.get("coords", ()))
    return atom.get("site_id"), atom.get("atom_name"), atom.get("fragment_id"), coords
