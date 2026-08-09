"""
纳米酶结构评分引擎：基于真实化学/物理指标。
1. 配位距离评分（与实验统计值比较）
2. 配位角度评分（实际角度与理想几何比较）
3. 配位完整性评分（配位数满足程度）
"""
import numpy as np
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .validation import coordination_angle_rms_deg, validate_assembly

# 配位距离合理范围（Å）
# PR4-1 (M4 fix): added Mo / W / V / Cr / Ru / Pd / Pt / Au / Ag for completeness.
# Bond-length ranges from CSD median ± 2σ for the most populated coordination
# geometries of each metal in CCDC database (octahedral / square planar).
_COORD_DISTANCE_RANGES = {
    ("FE", "N"): (1.80, 2.30), ("FE", "O"): (1.70, 2.35), ("FE", "S"): (2.10, 2.50),
    ("CU", "N"): (1.85, 2.20), ("CU", "O"): (1.85, 2.20), ("CU", "S"): (2.00, 2.40),
    ("ZN", "N"): (1.90, 2.30), ("ZN", "O"): (1.85, 2.25), ("ZN", "S"): (2.10, 2.50),
    ("MN", "N"): (1.95, 2.35), ("MN", "O"): (1.80, 2.25),
    ("CO", "N"): (1.90, 2.20), ("CO", "O"): (1.85, 2.15),
    ("NI", "N"): (1.90, 2.20), ("NI", "O"): (1.85, 2.15),
    # PR4-1 (M4): added second-row + third-row + d8 noble metals
    ("MO", "N"): (1.95, 2.30), ("MO", "O"): (1.85, 2.25), ("MO", "S"): (2.30, 2.55),
    ("W",  "N"): (1.95, 2.30), ("W",  "O"): (1.85, 2.25), ("W",  "S"): (2.30, 2.55),
    ("V",  "N"): (1.95, 2.30), ("V",  "O"): (1.80, 2.20), ("V",  "S"): (2.25, 2.50),
    ("CR", "N"): (1.95, 2.20), ("CR", "O"): (1.80, 2.15),
    ("RU", "N"): (1.95, 2.20), ("RU", "O"): (1.95, 2.15), ("RU", "S"): (2.25, 2.50),
    ("PD", "N"): (1.95, 2.15), ("PD", "O"): (1.95, 2.15), ("PD", "S"): (2.20, 2.45),
    ("PT", "N"): (1.95, 2.15), ("PT", "O"): (1.95, 2.15), ("PT", "S"): (2.20, 2.45),
    ("AU", "N"): (1.95, 2.20), ("AU", "S"): (2.20, 2.45),
    ("AG", "N"): (2.05, 2.40), ("AG", "S"): (2.30, 2.55),
}
COORDINATION_DISTANCE_RANGES = _COORD_DISTANCE_RANGES

_VALID_COORDINATION = {
    ("FE", 2): (4, 6), ("FE", 3): (4, 6),
    ("CU", 1): (2, 4), ("CU", 2): (4, 6),
    ("ZN", 2): (4, 6),
    ("MN", 2): (4, 6), ("MN", 3): (4, 6), ("MN", 4): (4, 6),
    ("CO", 2): (4, 6), ("CO", 3): (6, 6),
    ("NI", 2): (4, 6),
    # PR4-1 (M5): added Mo / W / V / Cr / Ru / Pd / Pt / Au / Ag
    ("MO", 4): (4, 6), ("MO", 5): (5, 7), ("MO", 6): (4, 8),
    ("W",  4): (4, 6), ("W",  5): (5, 7), ("W",  6): (4, 8),
    ("V",  3): (4, 6), ("V",  4): (4, 6), ("V",  5): (4, 6),
    ("CR", 2): (4, 6), ("CR", 3): (6, 6), ("CR", 6): (4, 6),
    ("RU", 2): (6, 6), ("RU", 3): (6, 6), ("RU", 4): (6, 6),
    ("PD", 2): (4, 4),                       # d8 → square planar (CN=4)
    ("PT", 2): (4, 4), ("PT", 4): (6, 6),    # Pt(II) SP, Pt(IV) Oh
    ("AU", 1): (2, 2), ("AU", 3): (4, 4),    # Au(I) linear, Au(III) SP
    ("AG", 1): (2, 4),                        # Ag(I) linear/Td
}
VALID_COORDINATION_RANGES = _VALID_COORDINATION

SCORE_WEIGHTS = {"distance": 0.40, "angle": 0.40, "coordination": 0.20}


@dataclass
class ConstraintScore:
    geometry_score: float = 0.0   # 角度评分（前端显示"几何"）
    distance_score: float = 0.0
    energy_score: float = 0.0     # Backward-compatible alias for distance_score.
    coordination_score: float = 0.0
    steric_score: float = 0.0
    stability_score: float = 0.0
    total_score: float = 0.0
    passed_hard_constraints: bool = True
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    method: str = "coordination_rules"
    backend: str = "rules"
    raw_energy_ev: Optional[float] = None
    energy_per_atom_ev: Optional[float] = None
    relaxed_energy_ev: Optional[float] = None
    max_force_ev_per_a: Optional[float] = None
    details: Dict[str, Any] = field(default_factory=dict)


def score_assembly(assembly: dict, design_spec, geometry_constraints: List[Dict] = None) -> ConstraintScore:
    score = ConstraintScore()

    validation = validate_assembly(assembly, design_spec, stage="preflight")
    score.details["physchem_validation"] = validation.to_dict()
    score.passed_hard_constraints = validation.passed
    score.errors.extend(
        issue.message for issue in validation.all_issues if issue.severity == "error"
    )
    score.warnings.extend(
        issue.message for issue in validation.all_issues if issue.severity != "error"
    )

    if not score.passed_hard_constraints:
        return score

    dist_scores, angle_scores, coord_scores = [], [], []

    for i, core in enumerate(assembly.get("cores", [])):
        dist_scores.append(_distance_score(core))
        angle_scores.append(_angle_score(core))

        if hasattr(design_spec, "metals") and i < len(design_spec.metals):
            expected_cn = design_spec.metals[i].coordination_number
            actual_cn = len(core["coord_atoms"])
            # PR2-2 (M22 fix): expected_cn=0 used to return a full 1.0 score —
            # semantically nonsense ("no coordination required → perfect"). It
            # now contributes 0.7 (a neutral score with a small downward bias)
            # so that a missing/zero CN spec doesn't artificially boost the
            # overall ranking. When expected_cn > 0 the original ratio is used.
            coord_scores.append(_coordination_number_score(actual_cn, expected_cn))

    score.distance_score = float(np.mean(dist_scores)) if dist_scores else 0.7
    score.energy_score = score.distance_score
    score.geometry_score = float(np.mean(angle_scores)) if angle_scores else 0.7
    score.coordination_score = float(np.mean(coord_scores)) if coord_scores else 0.7
    score.stability_score = score.distance_score
    score.total_score = (
        SCORE_WEIGHTS["distance"] * score.distance_score
        + SCORE_WEIGHTS["angle"] * score.geometry_score
        + SCORE_WEIGHTS["coordination"] * score.coordination_score
    )
    score.details["score_semantics"] = {
        "distance_score": score.distance_score,
        "energy_score": "deprecated alias of distance_score unless an ML backend overwrites it",
    }
    return score


def _distance_score(core: dict) -> float:
    """配位距离与实验统计范围的符合程度。"""
    metal = core["metal_type"].upper()
    metal_pos = np.array(core["metal"]["coords"])
    scores = []
    for a in core["coord_atoms"]:
        elem = a["element"].upper()
        dist = float(np.linalg.norm(np.array(a["coords"]) - metal_pos))
        rng = _COORD_DISTANCE_RANGES.get((metal, elem))
        if rng:
            ideal = (rng[0] + rng[1]) / 2
            half_width = (rng[1] - rng[0]) / 2
            dev = max(0.0, abs(dist - ideal) - half_width * 0.3)
            scores.append(float(np.exp(-(dev ** 2) / (2 * 0.1 ** 2))))
    return float(np.mean(scores)) if scores else 0.7


def _coordination_number_score(actual_cn: int, expected_cn: int) -> float:
    if expected_cn <= 0:
        return 0.7
    if actual_cn <= 0:
        return 0.0
    ratio = actual_cn / expected_cn
    return float(min(ratio, 1.0 / ratio))


def _angle_score(core: dict) -> float:
    """实际配位角度与理想几何的符合程度。"""
    metal_pos = np.array(core["metal"]["coords"])
    angle_rms = coordination_angle_rms_deg(
        core.get("geometry", "square_planar"),
        metal_pos,
        core.get("coord_atoms", []),
    )
    if angle_rms is None:
        return 0.7
    return float(np.exp(-(angle_rms ** 2) / (2 * 15.0 ** 2)))


coordination_number_score = _coordination_number_score
angle_score = _angle_score


def _coord_cutoff(metal: str, donor: str) -> float:
    """
    PR2-3 (NEW-7 fix): per-metal-donor coordination cutoff.

    The hardcoded 2.6 Å was too tight for Mn-S (~2.4 + margin) and too loose
    for Fe-O / Cu-O (true contacts ≤ 2.2 Å). Use _COORD_DISTANCE_RANGES upper
    bound + 0.2 Å as the cutoff so distant atoms are not mis-counted as
    coordinated. Falls back to 2.6 Å for unknown (metal, donor) pairs.
    """
    key = (metal.upper(), donor.upper())
    if key in _COORD_DISTANCE_RANGES:
        _, hi = _COORD_DISTANCE_RANGES[key]
        return hi + 0.2
    return 2.6


def coordination_cutoff(metal: str, donor: str) -> float:
    """Public per-metal donor coordination cutoff used by design templates."""
    return _coord_cutoff(metal, donor)


def _donor_element_of(atom: dict) -> str:
    """Best-effort donor element from atom dict (handles both explicit donor_element
    and element fields)."""
    el = (atom.get("donor_element") or atom.get("element") or "").upper()
    # Strip multi-letter elements down to first char to match N/O/S keys
    if el in {"N", "O", "S"}:
        return el
    if el and el[0] in {"N", "O", "S"}:
        return el[0]
    return el or "N"


def _check_chemical_validity(core: dict) -> tuple:
    """Legacy single-core helper retained for downstream imports.

    New assembly scoring uses :func:`validate_assembly`, which checks both the
    lower and upper bond limits and rejects CN=0.
    """
    metal = core["metal_type"].upper()
    ox = core.get("oxidation_state", 2)
    metal_pos = np.array(core["metal"]["coords"])
    # PR2-3 (NEW-7): per-pair cutoff instead of one hardcoded 2.6 Å for all
    # metal-donor combinations. Without this, Mn-S contacts at 2.45 Å (within
    # the Mn-S valid range 2.10–2.50) are mis-counted as non-coordinated, and
    # over-stretched Fe-O contacts at 2.55 Å (well outside Fe-O 1.70–2.35) are
    # mis-counted as coordinated.
    actual_coord = []
    for a in core["coord_atoms"]:
        donor = _donor_element_of(a)
        cutoff = _coord_cutoff(metal, donor)
        if np.linalg.norm(np.array(a["coords"]) - metal_pos) < cutoff:
            actual_coord.append(a)
    cn = len(actual_coord)
    key = (metal, ox)
    if cn == 0:
        return False, f"{metal}({ox}+) CN=0 has no valid coordination contacts"
    if key in _VALID_COORDINATION:
        lo, hi = _VALID_COORDINATION[key]
        if not (lo <= cn <= hi):
            return False, f"{metal}({ox}+) CN={cn} out of valid range [{lo},{hi}]"
    for atom in actual_coord:
        donor = _donor_element_of(atom)
        distance = float(np.linalg.norm(np.array(atom["coords"]) - metal_pos))
        valid_range = _COORD_DISTANCE_RANGES.get((metal, donor))
        if valid_range and not (valid_range[0] <= distance <= valid_range[1]):
            return False, (
                f"{metal}({ox}+) {donor} distance {distance:.3f} A out of "
                f"range [{valid_range[0]:.3f},{valid_range[1]:.3f}]"
            )
    return True, ""
