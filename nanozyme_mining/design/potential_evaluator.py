"""
ML-ready candidate evaluation for nanozyme designs.

The default path is dependency-light and uses coordination quality plus a
steric clash proxy. If E2N_MLP_BACKEND is set to "fairchem" or "mace", this
module lazily loads the corresponding ASE calculator and records energy/force
diagnostics without making those heavy packages required for the web app.
"""
from __future__ import annotations

import math
import os
import copy
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .constraint_scorer import (
    COORDINATION_DISTANCE_RANGES,
    ConstraintScore,
    score_assembly,
)
from ..utils.config import env_float, env_int, parse_bool
from ..utils.constants import ALL_METAL_ELEMENTS, CATALYTIC_METAL_ELEMENTS  # PR4-1 (M12/M13)


TOTAL_SCORE_WEIGHTS = {
    "geometry": 0.30,
    "coordination": 0.25,
    "energy": 0.30,
    "steric": 0.15,
}

_COVALENT_RADII = {
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


@dataclass(frozen=True)
class PotentialEvaluationConfig:
    backend: str = "geometry_proxy"
    model: Optional[str] = None
    task: str = "oc20"
    device: str = "cpu"
    relax: bool = True
    steps: int = 120
    pre_relax_steps: int = 30
    extended_steps: int = 200
    fmax: float = 0.05
    usable_fmax: float = 0.15
    charge: int = 0
    spin: int = 1
    default_dtype: str = "float32"
    mace_head: Optional[str] = None
    constrained_relax: bool = True
    fix_coordination_bonds: bool = True
    fix_metal_metal_distance: bool = True
    fix_support_backbone: bool = True
    freeze_support: bool = True
    support_freeze_radius: float = 6.0
    coordination_restraint_k: float = 100.0
    metal_metal_restraint_k: float = 10.0
    support_backbone_restraint_k: float = 35.0
    production_restraint_scale: float = 1.0
    coordination_restraint_buffer: float = 0.05
    release_harmonic_restraints: bool = False
    solvent: Optional[str] = "water"
    tblite_solvation_model: str = "alpb"
    tblite_local_steps: int = 80
    tblite_local_fmax: float = 0.10

    @classmethod
    def from_env(cls) -> "PotentialEvaluationConfig":
        backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
        fmax = env_float("E2N_MLP_FMAX", 0.05, min_value=0.001, max_value=10.0)
        usable_fmax = env_float(
            "E2N_MLP_USABLE_FMAX",
            0.15,
            min_value=fmax,
            max_value=20.0,
        )
        return cls(
            backend=backend or "geometry_proxy",
            model=os.environ.get("E2N_MLP_MODEL") or None,
            task=os.environ.get("E2N_MLP_TASK", "oc20"),
            device=os.environ.get("E2N_MLP_DEVICE", "cpu"),
            relax=parse_bool(os.environ.get("E2N_MLP_RELAX"), default=True),
            steps=env_int("E2N_MLP_STEPS", 120, min_value=1, max_value=10_000),
            pre_relax_steps=env_int(
                "E2N_MLP_PRE_STEPS",
                30,
                min_value=0,
                max_value=5_000,
            ),
            extended_steps=env_int(
                "E2N_MLP_EXTENDED_STEPS",
                200,
                min_value=1,
                max_value=20_000,
            ),
            fmax=fmax,
            usable_fmax=usable_fmax,
            charge=env_int("E2N_MLP_CHARGE", 0, min_value=-20, max_value=20),
            spin=env_int("E2N_MLP_SPIN", 1, min_value=1, max_value=21),
            default_dtype=os.environ.get("E2N_MLP_DTYPE", "float32"),
            mace_head=os.environ.get("E2N_MACE_HEAD") or None,
            constrained_relax=parse_bool(
                os.environ.get("E2N_MLP_CONSTRAINED_RELAX"),
                default=True,
            ),
            fix_coordination_bonds=parse_bool(
                os.environ.get("E2N_MLP_FIX_COORDINATION"),
                default=True,
            ),
            fix_metal_metal_distance=parse_bool(
                os.environ.get("E2N_MLP_FIX_METAL_METAL"),
                default=True,
            ),
            fix_support_backbone=parse_bool(
                os.environ.get("E2N_MLP_FIX_SUPPORT_BACKBONE"),
                default=True,
            ),
            freeze_support=parse_bool(
                os.environ.get("E2N_MLP_FREEZE_SUPPORT"),
                default=True,
            ),
            support_freeze_radius=env_float(
                "E2N_MLP_SUPPORT_FREEZE_RADIUS",
                6.0,
                min_value=0.0,
                max_value=100.0,
            ),
            coordination_restraint_k=env_float(
                "E2N_MLP_COORDINATION_RESTRAINT_K",
                100.0,
                min_value=0.0,
                max_value=100_000.0,
            ),
            metal_metal_restraint_k=env_float(
                "E2N_MLP_METAL_METAL_RESTRAINT_K",
                10.0,
                min_value=0.0,
                max_value=100_000.0,
            ),
            support_backbone_restraint_k=env_float(
                "E2N_MLP_SUPPORT_BACKBONE_RESTRAINT_K",
                35.0,
                min_value=0.0,
                max_value=100_000.0,
            ),
            production_restraint_scale=env_float(
                "E2N_MLP_PRODUCTION_RESTRAINT_SCALE",
                1.0,
                min_value=0.0,
                max_value=100.0,
            ),
            coordination_restraint_buffer=env_float(
                "E2N_MLP_COORDINATION_RESTRAINT_BUFFER",
                0.05,
                min_value=0.0,
                max_value=5.0,
            ),
            release_harmonic_restraints=parse_bool(
                os.environ.get("E2N_MLP_RELEASE_RESTRAINTS"),
                default=False,
            ),
            solvent=os.environ.get("E2N_TBLITE_SOLVENT", "water") or None,
            tblite_solvation_model=os.environ.get("E2N_TBLITE_SOLVATION", "alpb"),
            tblite_local_steps=env_int(
                "E2N_TBLITE_LOCAL_STEPS",
                80,
                min_value=1,
                max_value=10_000,
            ),
            tblite_local_fmax=env_float(
                "E2N_TBLITE_LOCAL_FMAX",
                0.10,
                min_value=0.001,
                max_value=10.0,
            ),
        )


def evaluate_assembly(
    assembly: dict,
    design_spec,
    geometry_constraints: Optional[List[Dict]] = None,
    config: Optional[PotentialEvaluationConfig] = None,
) -> ConstraintScore:
    """Return an interpretable score, optionally enriched by an ASE ML potential."""
    config = config or PotentialEvaluationConfig.from_env()
    score = score_assembly(assembly, design_spec, geometry_constraints)
    atoms = assembly.get("atoms", [])

    steric_score, steric_details = _steric_score(atoms)
    score.steric_score = steric_score
    score.details["steric"] = steric_details

    if not score.passed_hard_constraints:
        score.total_score = 0.0
        score.backend = "rules"
        score.method = "failed_hard_constraints"
        return score

    proxy_energy = 0.65 * score.energy_score + 0.35 * steric_score
    score.energy_score = proxy_energy
    score.stability_score = proxy_energy
    score.backend = "geometry_proxy"
    score.method = "coordination_steric_proxy"

    if config.backend in {"fairchem", "mace"}:
        try:
            ml_result = _evaluate_with_ml_backend(atoms, config)
        except Exception as exc:
            score.warnings.append(
                f"{config.backend} evaluation unavailable; using geometry proxy ({exc})"
            )
            score.details["ml_error"] = str(exc)
        else:
            score.backend = config.backend
            score.method = f"{config.backend}_{'relaxed' if config.relax else 'single_point'}"
            score.raw_energy_ev = ml_result["energy_ev"]
            score.energy_per_atom_ev = ml_result["energy_per_atom_ev"]
            score.relaxed_energy_ev = ml_result.get("relaxed_energy_ev")
            score.max_force_ev_per_a = ml_result.get("max_force_ev_per_a")
            score.details["ml"] = ml_result
            score.warnings.append(
                "Raw ML total energies are recorded for diagnostics; total_score uses force/stability "
                "signals because different dopants/compositions need consistent references."
            )
            force_score = _force_score(score.max_force_ev_per_a)
            score.energy_score = 0.75 * force_score + 0.25 * steric_score
            score.stability_score = score.energy_score
    elif config.backend not in {"", "none", "geometry_proxy", "rules"}:
        score.warnings.append(f"Unknown E2N_MLP_BACKEND={config.backend}; using geometry proxy")

    score.total_score = _combine_scores(score)
    return score


def _combine_scores(score: ConstraintScore) -> float:
    return float(
        TOTAL_SCORE_WEIGHTS["geometry"] * score.geometry_score
        + TOTAL_SCORE_WEIGHTS["coordination"] * score.coordination_score
        + TOTAL_SCORE_WEIGHTS["energy"] * score.energy_score
        + TOTAL_SCORE_WEIGHTS["steric"] * score.steric_score
    )


def _force_score(max_force_ev_per_a: Optional[float]) -> float:
    """
    PR2-1 (M19 fix): replace the single-Gaussian decay with a piecewise multi-tier
    score that maps to the conventions ML practitioners actually use.

    Tier interpretation (matches ASE optimization conventions):
      - <0.05 eV/A: fully converged   → 1.00
      - <0.10 eV/A: tightly relaxed   → 0.90
      - <0.30 eV/A: loosely relaxed   → 0.75
      - <0.50 eV/A: starting structure→ 0.55
      - <1.00 eV/A: noticeable strain → 0.40
      - <2.00 eV/A: poor geometry     → 0.25
      - >=2.00:     unphysical        → 0.10 (does not vanish to 0; downstream
                                              steric/coordination scores still
                                              carry signal)
    The old `exp(-(F/1.5)^2)` form mapped 1.5 eV/A → 0.37, which gave a
    near-passing grade to clearly strained initial guesses; this form is
    explicit about what "good" looks like.
    """
    if max_force_ev_per_a is None or not np.isfinite(max_force_ev_per_a):
        return 0.55  # no info → middle of the road
    f = float(max_force_ev_per_a)
    if f < 0.05:
        return 1.00
    if f < 0.10:
        return 0.90
    if f < 0.30:
        return 0.75
    if f < 0.50:
        return 0.55
    if f < 1.00:
        return 0.40
    if f < 2.00:
        return 0.25
    return 0.10


def _steric_score(atoms: List[Dict]) -> tuple[float, Dict]:
    if len(atoms) < 2:
        return 1.0, {"min_distance": None, "clash_count": 0}

    positions = np.array([a["coords"] for a in atoms], dtype=float)
    elements = [str(a["element"]).upper() for a in atoms]
    penalties = []
    min_distance = float("inf")

    for i in range(len(atoms)):
        ri = _COVALENT_RADII.get(elements[i], 0.8)
        for j in range(i + 1, len(atoms)):
            d = float(np.linalg.norm(positions[i] - positions[j]))
            if d < 1e-8:
                continue
            min_distance = min(min_distance, d)
            rj = _COVALENT_RADII.get(elements[j], 0.8)
            clash_threshold = 0.55 * (ri + rj)
            if d < clash_threshold:
                penalties.append((clash_threshold - d) / clash_threshold)

    if penalties:
        mean_penalty = float(np.mean(penalties))
        max_penalty = float(np.max(penalties))
        score = math.exp(-(8.0 * mean_penalty + 2.0 * max_penalty))
    else:
        score = 1.0

    return float(score), {
        "min_distance": None if not np.isfinite(min_distance) else min_distance,
        "clash_count": len(penalties),
        "max_clash_fraction": float(np.max(penalties)) if penalties else 0.0,
    }


def _evaluate_with_ml_backend(atoms: List[Dict], config: PotentialEvaluationConfig) -> Dict:
    ase_atoms = _to_ase_atoms(atoms, config)
    calc = _get_calculator(config)
    ase_atoms.calc = calc

    energy = float(ase_atoms.get_potential_energy())
    forces = np.array(ase_atoms.get_forces(), dtype=float)
    result = {
        "energy_ev": energy,
        "energy_per_atom_ev": energy / max(len(ase_atoms), 1),
        "max_force_ev_per_a": _max_force(forces),
        "initial_energy_ev": energy,
        "initial_max_force_ev_per_a": _max_force(forces),
        "atom_count": len(ase_atoms),
    }

    if config.relax:
        from ase.optimize import FIRE

        relaxed = ase_atoms.copy()
        constraint_plan = _apply_relaxation_constraints(relaxed, atoms, config)
        relaxed.calc = _calculator_with_harmonic_restraints(calc, constraint_plan, config)
        opt = FIRE(relaxed, logfile=None)
        converged = bool(opt.run(fmax=config.fmax, steps=config.steps))
        relaxed_energy = float(relaxed.get_potential_energy())
        relaxed_forces = np.array(relaxed.get_forces(), dtype=float)
        restraint_diagnostics = _restraint_diagnostics(relaxed, constraint_plan)
        result.update(
            {
                "relaxed_energy_ev": relaxed_energy,
                "relaxed_energy_per_atom_ev": relaxed_energy / max(len(relaxed), 1),
                "max_force_ev_per_a": _max_force(relaxed_forces),
                "relaxed_max_force_ev_per_a": _max_force(relaxed_forces),
                "relax_converged": converged,
                "relaxed_positions": np.array(relaxed.get_positions(), dtype=float).tolist(),
                "relaxation_constraints": constraint_plan,
                "relaxation_diagnostics": restraint_diagnostics,
            }
        )
    return result


def relax_atoms_with_ml_backend(
    atoms: List[Dict],
    config: Optional[PotentialEvaluationConfig] = None,
) -> tuple[List[Dict], Dict]:
    """Run ASE/FIRE relaxation and return atom dicts with relaxed coordinates."""
    config = config or PotentialEvaluationConfig.from_env()
    if config.backend not in {"fairchem", "mace"}:
        raise RuntimeError("structure relaxation requires E2N_MLP_BACKEND=fairchem or mace")
    if not config.relax:
        raise RuntimeError("structure relaxation requires E2N_MLP_RELAX=1")

    ase_atoms = _to_ase_atoms(atoms, config)
    calc = _get_calculator(config)
    ase_atoms.calc = calc
    initial_energy = float(ase_atoms.get_potential_energy())
    initial_forces = np.array(ase_atoms.get_forces(), dtype=float)
    constraint_plan = _apply_relaxation_constraints(ase_atoms, atoms, config)

    from ase.optimize import FIRE

    ase_atoms.calc = _calculator_with_harmonic_restraints(
        calc, constraint_plan, config, restraint_scale=1.0
    )
    pre_opt = FIRE(ase_atoms, logfile=None)
    pre_opt.run(fmax=max(config.fmax * 2.0, 0.10), steps=config.pre_relax_steps)

    ase_atoms.calc = _calculator_with_harmonic_restraints(
        calc,
        constraint_plan,
        config,
        restraint_scale=config.production_restraint_scale,
    )
    opt = FIRE(ase_atoms, logfile=None)
    opt.run(fmax=config.fmax, steps=config.steps)
    restrained_energy = float(ase_atoms.get_potential_energy())
    restrained_forces = np.array(ase_atoms.get_forces(), dtype=float)

    # Topology-preserving MACE relaxation is the production default for
    # transition-metal coordination motifs. An optional release stage probes
    # whether the chosen MACE model also preserves that topology unrestrained.
    release_steps = max(0, config.extended_steps - config.pre_relax_steps - config.steps)
    release_opt = None
    release_converged = False
    if config.release_harmonic_restraints and release_steps:
        ase_atoms.calc = calc
        release_opt = FIRE(ase_atoms, logfile=None)
        release_converged = bool(
            release_opt.run(fmax=config.fmax, steps=release_steps)
        )
        relaxed_energy = float(ase_atoms.get_potential_energy())
        relaxed_forces = np.array(ase_atoms.get_forces(), dtype=float)
        criterion_max_force = _max_force(relaxed_forces)
        raw_max_force = criterion_max_force
        converged = criterion_max_force is not None and criterion_max_force <= config.fmax
        relaxation_status = _classify_relaxation_status(
            criterion_max_force,
            fmax=config.fmax,
            usable_fmax=config.usable_fmax,
            constrained=False,
        )
        relaxation_mode = "released_unrestrained"
    else:
        ase_atoms.calc = calc
        relaxed_energy = float(ase_atoms.get_potential_energy())
        raw_forces = np.array(ase_atoms.get_forces(), dtype=float)
        raw_max_force = _max_force(raw_forces)
        criterion_max_force = _max_force(restrained_forces)
        converged = criterion_max_force is not None and criterion_max_force <= config.fmax
        constrained = bool(constraint_plan.get("enabled"))
        relaxation_status = _classify_relaxation_status(
            criterion_max_force,
            fmax=config.fmax,
            usable_fmax=config.usable_fmax,
            constrained=constrained,
        )
        relaxation_mode = (
            "topology_preserving_constrained"
            if constrained
            else "unrestrained"
        )
    restraint_diagnostics = _restraint_diagnostics(ase_atoms, constraint_plan)
    relaxed = copy.deepcopy(atoms)
    for atom, pos in zip(relaxed, np.array(ase_atoms.get_positions(), dtype=float)):
        atom["coords"] = [float(pos[0]), float(pos[1]), float(pos[2])]

    return relaxed, {
        "status": "success",
        "backend": config.backend,
        "model": config.model,
        "task": config.task,
        "device": config.device,
        "steps_requested": config.steps,
        "pre_relax_steps_run": int(getattr(pre_opt, "nsteps", config.pre_relax_steps)),
        "steps_run": (
            int(getattr(pre_opt, "nsteps", config.pre_relax_steps))
            + int(getattr(opt, "nsteps", config.steps))
            + (
                int(getattr(release_opt, "nsteps", release_steps))
                if release_opt is not None
                else 0
            )
        ),
        "fmax_target_ev_per_a": config.fmax,
        "usable_fmax_ev_per_a": config.usable_fmax,
        "relax_converged": converged,
        "relaxation_status": relaxation_status,
        "relaxation_mode": relaxation_mode,
        "extended_relaxation": release_opt is not None,
        "release_relaxation": {
            "steps_requested": release_steps,
            "steps_run": (
                int(getattr(release_opt, "nsteps", release_steps))
                if release_opt is not None
                else 0
            ),
            "converged": release_converged,
            "harmonic_restraints": False,
        },
        "initial_energy_ev": initial_energy,
        "initial_max_force_ev_per_a": _max_force(initial_forces),
        "raw_mace_energy_ev": relaxed_energy,
        "raw_mace_max_force_ev_per_a": raw_max_force,
        "restrained_objective_energy_ev": restrained_energy,
        "restrained_objective_max_force_ev_per_a": _max_force(restrained_forces),
        "relaxed_energy_ev": relaxed_energy,
        "relaxed_max_force_ev_per_a": criterion_max_force,
        "relaxation_constraints": constraint_plan,
        "relaxation_diagnostics": restraint_diagnostics,
    }


_METAL_ELEMENTS = CATALYTIC_METAL_ELEMENTS


def _classify_relaxation_status(
    max_force_ev_per_a: Optional[float],
    *,
    fmax: float,
    usable_fmax: float,
    constrained: bool,
) -> str:
    if max_force_ev_per_a is None or not np.isfinite(max_force_ev_per_a):
        return "rejected"
    if max_force_ev_per_a <= fmax:
        return "converged_constrained" if constrained else "converged"
    if max_force_ev_per_a <= usable_fmax:
        return "usable_constrained" if constrained else "usable_not_converged"
    return "rejected"


def classify_relaxation_status(
    max_force_ev_per_a: Optional[float],
    *,
    fmax: float,
    usable_fmax: float,
    constrained: bool,
) -> str:
    """Public relaxation-status classifier shared by assembly and screening code."""
    return _classify_relaxation_status(
        max_force_ev_per_a,
        fmax=fmax,
        usable_fmax=usable_fmax,
        constrained=constrained,
    )


def _relaxation_plan_for_atoms(atoms: List[Dict], config: PotentialEvaluationConfig) -> Dict:
    """Build a dependency-light ASE constraint plan for chemically safe relaxation."""
    plan = {
        "enabled": bool(config.constrained_relax),
        "coordination_bonds": [],
        "metal_metal_bonds": [],
        "support_backbone_bonds": [],
        "frozen_atom_indices": [],
        "frozen_atom_count": 0,
        "warnings": [],
    }
    if not config.constrained_relax:
        return plan

    metal_indices = [
        idx for idx, atom in enumerate(atoms)
        if _element(atom) in _METAL_ELEMENTS
    ]
    metal_by_site = {
        str(atoms[idx].get("site_id")): idx
        for idx in metal_indices
        if atoms[idx].get("site_id") is not None
    }

    if config.fix_coordination_bonds:
        for idx, atom in enumerate(atoms):
            if not atom.get("is_coord_atom"):
                continue
            site_id = str(atom.get("site_id"))
            metal_idx = metal_by_site.get(site_id)
            if metal_idx is None:
                plan["warnings"].append(f"coordination atom {idx} has no matching metal site {site_id}")
                continue
            declared_range = atom.get("bond_length_range")
            if declared_range and len(declared_range) == 2:
                validation_lower_bound = float(declared_range[0])
                validation_upper_bound = float(declared_range[1])
                buffer = max(0.0, float(config.coordination_restraint_buffer))
                lower_bound = validation_lower_bound + buffer
                upper_bound = validation_upper_bound - buffer
                if lower_bound >= upper_bound:
                    lower_bound = validation_lower_bound
                    upper_bound = validation_upper_bound
                target = 0.5 * (lower_bound + upper_bound)
            else:
                target = float(atom.get("bond_length") or _distance(atoms[metal_idx], atom))
                lower_bound = target - 0.10
                upper_bound = target + 0.10
                validation_lower_bound = lower_bound
                validation_upper_bound = upper_bound
            plan["coordination_bonds"].append(
                {
                    "pair": [metal_idx, idx],
                    "site_id": site_id,
                    "target_distance": target,
                    "lower_bound": lower_bound,
                    "upper_bound": upper_bound,
                    "validation_lower_bound": validation_lower_bound,
                    "validation_upper_bound": validation_upper_bound,
                    "kind": "metal_donor",
                }
            )

        # Bridged topologies add donor atoms to both metal cores dynamically in
        # ``_core_from_atoms``.  Those atoms are intentionally not tagged with
        # a single site_id during scaffold construction, so the old relaxation
        # plan omitted them and MACE could collapse a bridge M-N contact below
        # the strict validation range.  Restrain each bridge donor to its
        # nearest metal using the same metal/donor distance table as validation.
        existing_pairs = {
            tuple(item["pair"]) for item in plan["coordination_bonds"]
        }
        for idx, atom in enumerate(atoms):
            if not atom.get("is_bridge_atom") or idx in metal_indices:
                continue
            donor_element = _element(atom)
            candidates = sorted(
                (
                    (_distance(atoms[metal_idx], atom), metal_idx)
                    for metal_idx in metal_indices
                ),
                key=lambda item: item[0],
            )
            if not candidates or candidates[0][0] > 2.70:
                continue
            current_distance, metal_idx = candidates[0]
            pair = (metal_idx, idx)
            if pair in existing_pairs:
                continue
            metal_element = _element(atoms[metal_idx])
            declared_range = COORDINATION_DISTANCE_RANGES.get(
                (metal_element, donor_element)
            )
            if declared_range is None:
                validation_lower_bound = max(0.5, current_distance - 0.20)
                validation_upper_bound = current_distance + 0.20
            else:
                validation_lower_bound = float(declared_range[0])
                validation_upper_bound = float(declared_range[1])
            buffer = max(0.0, float(config.coordination_restraint_buffer))
            lower_bound = validation_lower_bound + buffer
            upper_bound = validation_upper_bound - buffer
            if lower_bound >= upper_bound:
                lower_bound = validation_lower_bound
                upper_bound = validation_upper_bound
            target = 0.5 * (lower_bound + upper_bound)
            plan["coordination_bonds"].append(
                {
                    "pair": [metal_idx, idx],
                    "site_id": atoms[metal_idx].get("site_id"),
                    "target_distance": target,
                    "lower_bound": lower_bound,
                    "upper_bound": upper_bound,
                    "validation_lower_bound": validation_lower_bound,
                    "validation_upper_bound": validation_upper_bound,
                    "kind": "bridge_metal_donor",
                }
            )
            existing_pairs.add(pair)

    if config.fix_metal_metal_distance and len(metal_indices) >= 2:
        ordered = sorted(
            metal_indices,
            key=lambda idx: str(atoms[idx].get("site_id", f"M{idx}")),
        )
        for left, right in zip(ordered, ordered[1:]):
            target = float(
                atoms[left].get("target_metal_distance")
                or atoms[right].get("target_metal_distance")
                or _distance(atoms[left], atoms[right])
            )
            declared_range = atoms[left].get("metal_metal_range") or atoms[right].get("metal_metal_range")
            lower_bound, upper_bound = (
                (float(declared_range[0]), float(declared_range[1]))
                if declared_range and len(declared_range) == 2
                else (target - 0.50, target + 0.50)
            )
            if target <= 16.0:
                plan["metal_metal_bonds"].append(
                    {
                        "pair": [left, right],
                        "target_distance": target,
                        "lower_bound": lower_bound,
                        "upper_bound": upper_bound,
                        "kind": "bimetal_topology",
                    }
                )

    if config.fix_support_backbone:
        support_indices = [
            idx for idx, atom in enumerate(atoms)
            if _is_support_atom(atom) and _element(atom) != "H"
        ]
        support_positions = {
            idx: np.array(atoms[idx]["coords"], dtype=float)
            for idx in support_indices
        }
        for left_i, left in enumerate(support_indices):
            for right in support_indices[left_i + 1:]:
                distance = float(np.linalg.norm(support_positions[right] - support_positions[left]))
                if 1.15 <= distance <= 1.82:
                    plan["support_backbone_bonds"].append(
                        {
                            "pair": [left, right],
                            "target_distance": distance,
                            "kind": "carbon_support_backbone",
                        }
                    )

    if config.freeze_support and metal_indices:
        metal_positions = np.array([atoms[idx]["coords"] for idx in metal_indices], dtype=float)
        frozen = []
        for idx, atom in enumerate(atoms):
            if not _is_support_atom(atom):
                continue
            pos = np.array(atom["coords"], dtype=float)
            min_dist = float(np.min(np.linalg.norm(metal_positions - pos, axis=1)))
            if min_dist >= config.support_freeze_radius:
                frozen.append(idx)
        plan["frozen_atom_indices"] = frozen
        plan["frozen_atom_count"] = len(frozen)

    return plan


def relaxation_plan_for_atoms(atoms: List[Dict], config: PotentialEvaluationConfig) -> Dict:
    """Public wrapper for the chemically safe ASE relaxation-constraint plan."""
    return _relaxation_plan_for_atoms(atoms, config)


def _apply_relaxation_constraints(ase_atoms, atoms: List[Dict], config: PotentialEvaluationConfig) -> Dict:
    plan = _relaxation_plan_for_atoms(atoms, config)
    if not plan["enabled"]:
        return plan

    constraints = []
    if plan["frozen_atom_indices"]:
        from ase.constraints import FixAtoms

        constraints.append(FixAtoms(indices=plan["frozen_atom_indices"]))

    if constraints:
        ase_atoms.set_constraint(constraints)
    plan["ase_constraint_count"] = len(constraints)
    plan["restrained_pair_count"] = (
        len(plan["coordination_bonds"])
        + len(plan["metal_metal_bonds"])
        + len(plan["support_backbone_bonds"])
    )
    plan["fixed_bond_count"] = 0
    return plan


def _calculator_with_harmonic_restraints(
    calc,
    plan: Dict,
    config: PotentialEvaluationConfig,
    restraint_scale: float = 1.0,
):
    pairs = (
        (plan or {}).get("coordination_bonds", [])
        + (plan or {}).get("metal_metal_bonds", [])
        + (plan or {}).get("support_backbone_bonds", [])
    )
    if not (plan or {}).get("enabled") or not pairs:
        return calc

    from ase.calculators.calculator import Calculator, all_changes

    class HarmonicRestraintCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self, base_calc, restraint_plan, potential_config):
            super().__init__()
            self.base_calc = base_calc
            self.restraint_plan = restraint_plan
            self.potential_config = potential_config

        def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.base_calc.calculate(atoms, properties, system_changes)
            energy = float(self.base_calc.results["energy"])
            forces = np.array(self.base_calc.results["forces"], dtype=float).copy()
            positions = np.array(atoms.get_positions(), dtype=float)

            restraint_energy = 0.0
            max_delta = 0.0
            for item in self.restraint_plan.get("coordination_bonds", []):
                e, max_delta = _add_pair_restraint(
                    positions,
                    forces,
                    item,
                    self.potential_config.coordination_restraint_k * restraint_scale,
                    max_delta,
                )
                restraint_energy += e
            for item in self.restraint_plan.get("metal_metal_bonds", []):
                e, max_delta = _add_pair_restraint(
                    positions,
                    forces,
                    item,
                    self.potential_config.metal_metal_restraint_k * restraint_scale,
                    max_delta,
                )
                restraint_energy += e
            for item in self.restraint_plan.get("support_backbone_bonds", []):
                e, max_delta = _add_pair_restraint(
                    positions,
                    forces,
                    item,
                    self.potential_config.support_backbone_restraint_k * restraint_scale,
                    max_delta,
                )
                restraint_energy += e

            self.results["energy"] = energy + restraint_energy
            self.results["forces"] = forces
            self.results["raw_energy_ev"] = energy
            self.results["restraint_energy_ev"] = restraint_energy
            self.results["max_restraint_delta_a"] = max_delta

    return HarmonicRestraintCalculator(calc, plan, config)


def calculator_with_harmonic_restraints(
    calc,
    plan: Dict,
    config: PotentialEvaluationConfig,
    restraint_scale: float = 1.0,
):
    """Public wrapper that decorates an ASE calculator with E2N restraint forces."""
    return _calculator_with_harmonic_restraints(
        calc,
        plan,
        config,
        restraint_scale=restraint_scale,
    )


def _restraint_diagnostics(ase_atoms, plan: Dict) -> Dict:
    positions = np.array(ase_atoms.get_positions(), dtype=float)
    return _restraint_diagnostics_from_positions(positions, plan)


def _restraint_diagnostics_from_positions(positions: np.ndarray, plan: Dict) -> Dict:
    positions = np.array(positions, dtype=float)
    coordination = _pair_distance_stats(positions, (plan or {}).get("coordination_bonds", []))
    metal_metal = _pair_distance_stats(positions, (plan or {}).get("metal_metal_bonds", []))
    support_backbone = _pair_distance_stats(positions, (plan or {}).get("support_backbone_bonds", []))
    deltas = [
        abs(float(item.get("delta_a", 0.0)))
        for item in coordination["pairs"] + metal_metal["pairs"] + support_backbone["pairs"]
    ]
    return {
        "restrained_pair_count": len(deltas),
        "max_abs_delta_a": max(deltas) if deltas else 0.0,
        "rms_delta_a": float(math.sqrt(sum(d * d for d in deltas) / len(deltas))) if deltas else 0.0,
        "coordination": coordination,
        "metal_metal": metal_metal,
        "support_backbone": support_backbone,
    }


def restraint_diagnostics_from_positions(positions: np.ndarray, plan: Dict) -> Dict:
    """Public diagnostics helper for validated restraint-plan tests and reports."""
    return _restraint_diagnostics_from_positions(positions, plan)


def _pair_distance_stats(positions: np.ndarray, pairs: List[Dict]) -> Dict:
    rows = []
    abs_deltas = []
    for item in pairs:
        i, j = item["pair"]
        dist = float(np.linalg.norm(positions[j] - positions[i]))
        target = float(item["target_distance"])
        delta = dist - target
        rows.append(
            {
                "pair": list(item["pair"]),
                "kind": item.get("kind"),
                "target_distance_a": target,
                "allowed_range_a": [item.get("lower_bound"), item.get("upper_bound")]
                if item.get("lower_bound") is not None else None,
                "distance_a": dist,
                "delta_a": delta,
            }
        )
        abs_deltas.append(abs(delta))
    return {
        "count": len(rows),
        "max_abs_delta_a": max(abs_deltas) if abs_deltas else 0.0,
        "mean_abs_delta_a": float(np.mean(abs_deltas)) if abs_deltas else 0.0,
        "pairs": rows,
    }


def _add_pair_restraint(
    positions: np.ndarray,
    forces: np.ndarray,
    item: Dict,
    spring_k: float,
    max_delta: float,
) -> tuple[float, float]:
    i, j = item["pair"]
    target = float(item["target_distance"])
    vec = positions[j] - positions[i]
    dist = float(np.linalg.norm(vec))
    if dist < 1e-8:
        return 0.0, max_delta
    lower = float(item.get("lower_bound", target))
    upper = float(item.get("upper_bound", target))
    if lower <= dist <= upper:
        return 0.0, max_delta
    boundary = lower if dist < lower else upper
    delta = dist - boundary
    unit = vec / dist
    force = spring_k * delta * unit
    forces[i] += force
    forces[j] -= force
    return 0.5 * spring_k * delta * delta, max(max_delta, abs(delta))


def _element(atom: Dict) -> str:
    return str(atom.get("element", "")).upper()


def _distance(left: Dict, right: Dict) -> float:
    return float(np.linalg.norm(np.array(left["coords"], dtype=float) - np.array(right["coords"], dtype=float)))


def _is_support_atom(atom: Dict) -> bool:
    residue = str(atom.get("residue_name", "")).upper()
    return residue == "GRA" or residue.endswith("DP")


def evaluate_atoms_energy(
    atoms: List[Dict],
    config: Optional[PotentialEvaluationConfig] = None,
) -> Dict:
    """Evaluate an arbitrary atom list with the configured ML potential.

    This is intentionally stricter than ``evaluate_assembly``: geometry-only
    proxy scores are useful for ranking generated motifs, but adsorption
    energies and NEB barriers need a real, consistent potential energy surface.
    """
    config = config or PotentialEvaluationConfig.from_env()
    if config.backend not in {"fairchem", "mace", "tblite"}:
        raise RuntimeError(
            "adsorption/transition-state energies require E2N_MLP_BACKEND=fairchem, mace, or tblite"
        )
    return _evaluate_with_ml_backend(atoms, config)


def atoms_to_ase(atoms: List[Dict], config: Optional[PotentialEvaluationConfig] = None):
    """Public wrapper for scripts that need ASE images from E2N atom dicts."""
    return _to_ase_atoms(atoms, config or PotentialEvaluationConfig.from_env())


def get_ase_calculator(config: Optional[PotentialEvaluationConfig] = None):
    """Return the cached ASE calculator for the selected ML backend."""
    return _get_calculator(config or PotentialEvaluationConfig.from_env())


def calculator_capabilities(config: PotentialEvaluationConfig) -> Dict[str, bool]:
    """Return conservative electronic-state capabilities for a backend.

    The bundled MACE-MH-1 checkpoint is a ScaleShiftMACE without charge/spin
    embeddings. Its ASE adapter accepts these metadata fields, but the model
    energy is unchanged by them.
    """
    backend = config.backend.lower()
    if backend == "tblite":
        return {"energy": True, "forces": True, "charge": True, "spin": True}
    if backend == "mace":
        return {"energy": True, "forces": True, "charge": False, "spin": False}
    if backend == "fairchem":
        return {"energy": True, "forces": True, "charge": False, "spin": False}
    return {"energy": False, "forces": False, "charge": False, "spin": False}


def _to_ase_atoms(atoms: List[Dict], config: PotentialEvaluationConfig):
    try:
        from ase import Atoms
    except ImportError as exc:
        raise RuntimeError("ASE is required for ML potential evaluation") from exc

    symbols = [str(a["element"]).capitalize() for a in atoms]
    positions = np.array([a["coords"] for a in atoms], dtype=float)
    if len(positions) == 0:
        raise ValueError("Cannot evaluate an empty structure")

    span = np.ptp(positions, axis=0)
    cell_lengths = np.maximum(span + 16.0, 20.0)
    ase_atoms = Atoms(symbols=symbols, positions=positions, cell=cell_lengths, pbc=False)
    # PR2-1 (M18 nuance): MACE's MACECalculator DOES read atoms.info — the default
    # info_keys mapping (mace/calculators/mace.py:163-168) sends:
    #     "spin"   → total_spin
    #     "charge" → total_charge
    # so the four setdefault() calls below are the *primary* path that wires
    # config.spin / config.charge into the calculator. v4 audit M18's premise
    # ("ASE calculator doesn't read info") was wrong about MACE specifically;
    # we keep the info path and ALSO populate set_initial_magnetic_moments as
    # a belt-and-suspenders signal for non-MACE backends (LJ/EMT diagnostics,
    # FairChem variants) that consume only the per-atom array.
    metal_set = CATALYTIC_METAL_ELEMENTS
    metal_indices = [i for i, a in enumerate(atoms)
                      if str(a.get("element", "")).upper() in metal_set]
    moments = np.zeros(len(atoms), dtype=float)
    if metal_indices and config.spin and config.spin > 1:
        # Convention: config.spin is 2S+1 (multiplicity); convert to S then to
        # per-metal moment by splitting unpaired electrons across all metals.
        unpaired = float(config.spin - 1)  # S = (2S+1 - 1)/2 → unpaired = 2S
        per_metal = unpaired / len(metal_indices)
        for idx in metal_indices:
            moments[idx] = per_metal
        ase_atoms.set_initial_magnetic_moments(moments)
    if config.charge:
        # ASE doesn't carry total charge directly; calculators that need it read
        # info["total_charge"] or atoms.calc.parameters. We set both so downstream
        # backends find what they need.
        try:
            ase_atoms.set_initial_charges([0.0] * len(atoms))  # per-atom zero baseline
        except Exception:
            pass
    ase_atoms.info.setdefault("charge", config.charge)
    ase_atoms.info.setdefault("spin", config.spin)
    ase_atoms.info.setdefault("total_charge", config.charge)
    ase_atoms.info.setdefault("total_spin", config.spin)
    return ase_atoms


def _max_force(forces: np.ndarray) -> Optional[float]:
    if forces.size == 0:
        return None
    return float(np.max(np.linalg.norm(forces, axis=1)))


@lru_cache(maxsize=8)
def _get_calculator(config: PotentialEvaluationConfig):
    """
    PR2-1 (N-M8 note): lru_cache key is the frozen `config` dataclass itself,
    so any change to backend/device/model/dtype/head produces a different key.
    `device` is part of PotentialEvaluationConfig → switching cpu↔mps↔cuda
    automatically gets a fresh calculator. This was the original concern in
    N-M8 but the frozen-dataclass key already addresses it; documenting here
    so future readers don't re-raise the alarm.
    """
    import logging as _log
    log = _log.getLogger("e2n.potential_evaluator")

    if config.backend == "tblite":
        from tblite.ase import TBLite
        kwargs = dict(
            method=config.model or "GFN2-xTB",
            charge=config.charge,
            multiplicity=config.spin,
            verbosity=0,
        )
        if config.solvent:
            model = str(config.tblite_solvation_model or "alpb").lower()
            if model not in {"alpb", "gbsa", "cpcm", "gbe", "gb"}:
                raise ValueError(f"Unsupported tblite solvation model: {model}")
            kwargs[f"{model}_solvation"] = config.solvent
        return TBLite(**kwargs)

    if config.backend == "fairchem":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        model = config.model or "uma-s-1p2"
        predictor = pretrained_mlip.get_predict_unit(model, device=config.device)
        return FAIRChemCalculator(predictor, task_name=config.task)

    if config.backend == "mace":
        try:
            from mace.calculators import mace_mp
        except ImportError:
            mace_mp = None

        if config.model and Path(config.model).exists():
            from mace.calculators import MACECalculator

            kwargs = dict(
                model_paths=[config.model],
                device=config.device,
                default_dtype=config.default_dtype,
            )
            if config.mace_head:
                kwargs["head"] = config.mace_head
            try:
                return MACECalculator(**kwargs)
            except TypeError as e:
                # PR2-1 (M20 fix): the old code silently dropped `head` on TypeError,
                # so a misspelled MACE head or a version mismatch would produce
                # results computed against the *wrong* head with no warning.
                # Now we log a warning before retrying without head.
                if config.mace_head:
                    log.warning(
                        "MACECalculator rejected head=%r (%s); retrying without head — "
                        "energies/forces may differ from the intended head!",
                        config.mace_head, e,
                    )
                kwargs.pop("head", None)
                return MACECalculator(**kwargs)
        if mace_mp is None:
            raise RuntimeError("mace is not installed")
        kwargs = dict(
            model=config.model or "medium",
            device=config.device,
            default_dtype=config.default_dtype,
        )
        if config.mace_head:
            kwargs["head"] = config.mace_head
        return mace_mp(**kwargs)

    raise RuntimeError(f"Unsupported ML backend: {config.backend}")
