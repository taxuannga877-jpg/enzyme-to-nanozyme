"""
主组装入口：生成石墨烯片段纳米酶，支持批量变体生成。
"""
import uuid
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

from .design_spec import DesignSpec
from .constraint_scorer import ConstraintScore
from .carbon_scaffold import (
    build_graphene_fragment, embed_metal_in_graphene,
    embed_metal_in_carbon_network, embed_bridge_in_carbon_network,
    build_second_shell_on_graphene, passivate_graphene_edges, rdkit_optimize
)
from .bimetallic_topology import (
    BimetallicTopology,
    propose_bimetallic_topologies,
    score_metal_distance,
)
from .potential_evaluator import (
    PotentialEvaluationConfig,
    TOTAL_SCORE_WEIGHTS,
    evaluate_assembly,
    relax_atoms_with_ml_backend,
)
from .chemical_system import annotate_chemical_system
from .physchem_knowledge import enrich_design_spec, evaluate_constructibility
from .validation import validate_assembly

# PR4-1 (M14 fix): named cutoffs so the assemble_bimetallic call sites are
# self-documenting (previously magic numbers 1.10 and 0.95 looked like a bug).
SUPPORT_CLASH_CUTOFF_DEFAULT = 1.35  # metal/first-shell clearance before relaxation
SUPPORT_CLASH_CUTOFF_LINKER  = 0.95  # diimine / mu-bridge atoms (smaller)


@dataclass
class AssemblyResult:
    job_id: str
    atoms: List[Dict]
    cores: List[Dict]
    second_shell_atoms: List[Dict]
    score: ConstraintScore
    design_spec: DesignSpec
    bond_graph: List[Tuple[int, int, str]] = field(default_factory=list)
    formal_charge: int = 0
    spin_multiplicities: List[int] = field(default_factory=lambda: [1])
    chemistry_warnings: List[str] = field(default_factory=list)
    smiles: str = ""
    xyz: str = ""
    label: str = ""
    error: Optional[str] = None


class NanozymeAssembler:

    def assemble(self, spec: DesignSpec, geometry_constraints: List[Dict] = None) -> AssemblyResult:
        """单次组装（兼容旧接口）。"""
        results = self.assemble_batch(spec, geometry_constraints)
        return results[0] if results else _empty_result(spec)

    def assemble_batch(self, spec: DesignSpec, geometry_constraints: List[Dict] = None) -> List[AssemblyResult]:
        """
        批量生成所有变体：
        - 单金属：4种掺杂（none/N/S/NS）
        - 双金属：双活性拓扑图（邻位/协同桥联/级联双中心）
        """
        spec = enrich_design_spec(spec)
        decision = evaluate_constructibility(spec)
        if not decision.constructible:
            return [_not_constructible_result(spec, decision.to_dict())]

        results = []
        dopings = ["none", "N", "S", "NS"]

        if not spec.metals:
            return [_empty_result(spec)]

        if len(spec.metals) == 1:
            for doping in dopings:
                try:
                    r = self._assemble_single(spec, doping, geometry_constraints)
                except (ValueError, RuntimeError) as exc:
                    r = _generation_failure_result(spec, str(exc), label=f"single/{doping}")
                results.append(r)
        else:
            topologies = propose_bimetallic_topologies(spec)
            if not topologies:
                return [_not_constructible_result(spec, {
                    **decision.to_dict(),
                    "reason_codes": list(decision.reason_codes) + ["no_supported_topology"],
                })]
            for topology in topologies:
                try:
                    r = self._assemble_bimetallic(spec, topology, geometry_constraints)
                except (ValueError, RuntimeError) as exc:
                    r = _generation_failure_result(spec, str(exc), label=topology.edge.relation)
                results.append(r)
        return sorted(results, key=lambda r: r.score.total_score, reverse=True)

    def _assemble_single(self, spec: DesignSpec, doping: str, gc) -> AssemblyResult:
        metal_spec = spec.metals[0]
        coord_residues = [ca.residue_name for ca in metal_spec.coord_atoms]
        if not coord_residues:
            coord_residues = ["HIS"] * metal_spec.coordination_number

        graphene = build_graphene_fragment(radius=6.5)
        site_id = "M0"
        atoms = embed_metal_in_graphene(
            graphene, metal_spec.metal_type, coord_residues,
            metal_offset=np.zeros(3), doping=doping,
            coordination_number=metal_spec.coordination_number,
            coordination_geometry=metal_spec.coordination_geometry,
            coord_atoms=metal_spec.coord_atoms,
            site_id=site_id,
        )

        # 找金属原子坐标
        metal_center = np.array(next(a["coords"] for a in atoms if a["element"] == metal_spec.metal_type))

        # 第二配位层
        ss_specs = [{"residue_name": s.residue_name, "atom_name": s.atom_name,
                     "role": s.role, "distance_to_metal": s.distance_to_metal}
                    for s in spec.second_shell]
        ss_atoms = build_second_shell_on_graphene(atoms, metal_center, ss_specs)
        atoms = atoms + ss_atoms

        atoms = rdkit_optimize(atoms, metal_spec.metal_type)
        atoms = passivate_graphene_edges(atoms)
        atoms = _remove_isolated_support_artifacts(atoms)
        atoms = passivate_graphene_edges(atoms)
        pre_core = _core_from_atoms(atoms, metal_spec, site_id)
        pre_assembly = {"atoms": atoms, "cores": [pre_core], "linker_atoms": [], "mode": "single",
                    "second_shell": ss_atoms}
        pre_score = evaluate_assembly(pre_assembly, spec, gc, config=PotentialEvaluationConfig())
        if not pre_score.passed_hard_constraints:
            relax_details = {
                "status": "skipped",
                "backend": PotentialEvaluationConfig.from_env().backend,
                "reason": "preflight physicochemical validation failed",
            }
            _attach_relaxation_details(pre_score, relax_details)
            return _result_from_scored_assembly(
                spec, atoms, [pre_core], ss_atoms, pre_score,
                f"{metal_spec.metal_type}-invalid | {_doping_label(doping)}",
            )

        atoms, relax_details = _maybe_relax_atoms(atoms)
        core = _core_from_atoms(atoms, metal_spec, site_id)
        assembly = {"atoms": atoms, "cores": [core], "linker_atoms": [], "mode": "single",
                    "second_shell": ss_atoms}
        score = evaluate_assembly(assembly, spec, gc, config=_score_config_after_relax(relax_details))
        score.details["preflight_validation"] = pre_score.details.get("physchem_validation")
        _attach_relaxation_details(score, relax_details)
        label = f"{metal_spec.metal_type}-N{len(core['coord_atoms'])} | {_doping_label(doping)}"

        return _result_from_scored_assembly(spec, atoms, [core], ss_atoms, score, label)

    def _assemble_bimetallic(self, spec: DesignSpec, topology: BimetallicTopology, gc) -> AssemblyResult:
        radius = max(7.0, topology.edge.ideal_distance / 2.0 + 5.5)
        all_atoms = build_graphene_fragment(radius=radius)
        for node in topology.nodes:
            _apply_local_doping(all_atoms, np.array(node.position[:2], dtype=float), topology.edge.doping)

        center_records = []
        for node in topology.nodes:
            metal_spec = node.metal
            coord_residues = [ca.residue_name for ca in metal_spec.coord_atoms]
            if not coord_residues:
                coord_residues = ["HIS"] * metal_spec.coordination_number

            before_names = {atom.get("atom_name") for atom in all_atoms}
            partner = topology.nodes[1 - node.index]
            partner_direction = np.asarray(partner.position[:2], dtype=float) - np.asarray(node.position[:2], dtype=float)
            all_atoms = embed_metal_in_carbon_network(
                all_atoms, metal_spec.metal_type, coord_residues,
                metal_offset=np.array(node.position, dtype=float),
                coordination_number=metal_spec.coordination_number,
                coordination_geometry=metal_spec.coordination_geometry,
                coord_atoms=metal_spec.coord_atoms,
                site_id=node.site_id,
                reserve_coordination_slots=1 if topology.edge.relation == "bridged" else 0,
                reserved_direction_xy=partner_direction,
            )
            metal_center = np.array(
                next(
                    a["coords"] for a in all_atoms
                    if a.get("site_id") == node.site_id
                    and str(a["element"]).upper() == metal_spec.metal_type.upper()
                ),
                dtype=float,
            )
            _assert_no_detached_site_fragment(all_atoms, node.site_id, before_names)
            center_records.append((node, metal_center))

        metal_centers = [center for _, center in center_records]
        for node, _center in center_records:
            for atom in all_atoms:
                if atom.get("site_id") == node.site_id and str(atom.get("element", "")).upper() == node.metal.metal_type.upper():
                    atom["target_metal_distance"] = float(topology.edge.ideal_distance)
                    atom["metal_metal_range"] = list(topology.edge.distance_range)
                    atom["bimetallic_relation"] = topology.edge.relation
        metal_atoms = [
            atom for atom in all_atoms
            if str(atom.get("element", "")).upper() == str(atom.get("atom_name", "")).upper()
            and atom.get("is_embedded_metal")
        ]
        all_atoms = _drop_support_clashes(all_atoms, metal_atoms, cutoff=1.15)
        linker_atoms = embed_bridge_in_carbon_network(all_atoms, topology, metal_centers)

        all_atoms = passivate_graphene_edges(all_atoms)

        # 第二配位层
        ss_atoms = []
        for node, mc in center_records:
            ss_specs = [
                {
                    "residue_name": s.residue_name,
                    "atom_name": s.atom_name,
                    "role": s.role,
                    "distance_to_metal": s.distance_to_metal,
                }
                for s in spec.second_shell
                if s.target_metal_idx == node.index
            ]
            ss_atoms.extend(build_second_shell_on_graphene(all_atoms, mc, ss_specs, n_sites=2))
        all_atoms = all_atoms + ss_atoms

        all_atoms = rdkit_optimize(all_atoms, spec.metals[0].metal_type)
        all_atoms = passivate_graphene_edges(all_atoms)
        all_atoms = _remove_isolated_support_artifacts(all_atoms)
        all_atoms = passivate_graphene_edges(all_atoms)
        pre_cores = [_core_from_atoms(all_atoms, node.metal, node.site_id) for node in topology.nodes]
        pre_assembly = {"atoms": all_atoms, "cores": pre_cores, "linker_atoms": linker_atoms,
                        "mode": topology.edge.relation, "second_shell": ss_atoms}
        pre_score = evaluate_assembly(pre_assembly, spec, gc, config=PotentialEvaluationConfig())
        if not pre_score.passed_hard_constraints:
            relax_details = {
                "status": "skipped",
                "backend": PotentialEvaluationConfig.from_env().backend,
                "reason": "preflight physicochemical validation failed",
            }
            _attach_bimetallic_details(
                pre_score,
                topology,
                [np.asarray(core["metal"]["coords"], dtype=float) for core in pre_cores],
            )
            _attach_relaxation_details(pre_score, relax_details)
            return _result_from_scored_assembly(
                spec, all_atoms, pre_cores, ss_atoms, pre_score,
                f"not_constructible | {topology.edge.label}",
            )

        all_atoms, relax_details = _maybe_relax_atoms(all_atoms)
        cores = [_core_from_atoms(all_atoms, node.metal, node.site_id) for node in topology.nodes]
        relaxed_centers = [np.array(core["metal"]["coords"]) for core in cores]

        assembly = {"atoms": all_atoms, "cores": cores, "linker_atoms": linker_atoms,
                    "mode": topology.edge.relation, "second_shell": ss_atoms}
        score = evaluate_assembly(assembly, spec, gc, config=_score_config_after_relax(relax_details))
        score.details["preflight_validation"] = pre_score.details.get("physchem_validation")
        _attach_bimetallic_details(score, topology, relaxed_centers)
        _attach_relaxation_details(score, relax_details)
        bimetallic = score.details.get("bimetallic", {})
        center_label = " + ".join(
            f"{node.activity or 'Activity'}/{node.metal.metal_type}"
            for node in topology.nodes
        )
        distance = bimetallic.get("metal_distance") or topology.edge.ideal_distance
        label = (
            f"{center_label} | {topology.edge.label} | "
            f"{distance:.1f}Å | {_doping_label(topology.edge.doping)}"
        )

        return _result_from_scored_assembly(spec, all_atoms, cores, ss_atoms, score, label)


def _doping_label(d: str) -> str:
    return {
        "none": "Undoped",
        "N": "N-doped",
        "S": "S-doped",
        "NS": "N/S-doped",
    }.get(d, d)

def _mode_label(m: str) -> str:
    return {
        "bridged": "Bridged",
        "independent": "Independent",
        "independent_adjacent": "Independent adjacent",
        "independent_separated": "Independent non-adjacent",
    }.get(m, m)


def _orient_site_fragment(
    site_atoms: List[Dict], existing_atoms: List[Dict], metal_center: np.ndarray
) -> List[Dict]:
    """Rigidly rotate a new metal site to avoid the already placed site."""
    occupied = [
        atom
        for atom in existing_atoms
        if str(atom.get("residue_name", "")).upper() not in {"GRA", "NDP", "SDP"}
    ]
    if not occupied:
        return [dict(atom) for atom in site_atoms]

    occupied_coords = np.array([atom["coords"] for atom in occupied], dtype=float)
    coords = np.array([atom["coords"] for atom in site_atoms], dtype=float)
    centered = coords - metal_center
    best_coords = coords
    best_clearance = -float("inf")
    for angle in np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False):
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=float,
        )
        trial = centered @ rotation.T + metal_center
        clearance = min(
            float(np.min(np.linalg.norm(occupied_coords - coord, axis=1)))
            for coord in trial
        )
        if clearance > best_clearance:
            best_clearance = clearance
            best_coords = trial

    output = [dict(atom) for atom in site_atoms]
    for atom, coord in zip(output, best_coords):
        atom["coords"] = coord.tolist()
    return output


def _reserve_bridge_coordination_slot(
    site_atoms: List[Dict], metal_center: np.ndarray, partner_center: np.ndarray
) -> List[Dict]:
    """Remove the ligand fragment occupying the future inter-metal bridge axis."""
    direction = partner_center - metal_center
    direction /= np.linalg.norm(direction) or 1.0
    candidates = []
    for atom in site_atoms:
        if not atom.get("is_coord_atom") or not atom.get("fragment_id"):
            continue
        vector = np.array(atom["coords"], dtype=float) - metal_center
        vector /= np.linalg.norm(vector) or 1.0
        candidates.append((float(np.dot(vector, direction)), atom["fragment_id"]))
    if not candidates:
        return [dict(atom) for atom in site_atoms]
    _, fragment_id = max(candidates)
    return [dict(atom) for atom in site_atoms if atom.get("fragment_id") != fragment_id]


def _assert_no_detached_site_fragment(
    atoms: List[Dict], site_id: str, _before_names: set
) -> None:
    detached = [
        atom
        for atom in atoms
        if atom.get("site_id") == site_id
        and atom.get("is_coord_atom")
        and atom.get("fragment_id")
    ]
    if detached:
        raise RuntimeError(
            f"bimetallic site {site_id} generated detached ligand fragments "
            "instead of carbon-network donor atoms"
        )

def _empty_result(spec: DesignSpec) -> AssemblyResult:
    return AssemblyResult(job_id="error", atoms=[], cores=[], second_shell_atoms=[],
                          score=ConstraintScore(passed_hard_constraints=False, errors=["No metals"]),
                          design_spec=spec, error="No metals specified")


def _not_constructible_result(spec: DesignSpec, decision: Dict) -> AssemblyResult:
    reasons = list(decision.get("reason_codes") or ["not_constructible"])
    score = ConstraintScore(
        passed_hard_constraints=False,
        errors=reasons,
        method="physchem_constructibility_gate",
        backend="rules",
        details={"constructibility": decision},
    )
    return AssemblyResult(
        job_id=str(uuid.uuid4())[:8],
        atoms=[],
        cores=[],
        second_shell_atoms=[],
        score=score,
        design_spec=spec,
        label="not_constructible",
        error="; ".join(reasons),
    )


def _generation_failure_result(spec: DesignSpec, error: str, label: str) -> AssemblyResult:
    score = ConstraintScore(
        passed_hard_constraints=False,
        errors=[error],
        method="coordination_pocket_generation_failed",
        backend="rules",
        details={"constructibility": {"status": "not_constructible", "reason_codes": [error]}},
    )
    return AssemblyResult(
        job_id=str(uuid.uuid4())[:8], atoms=[], cores=[], second_shell_atoms=[],
        score=score, design_spec=spec, label=label, error=error,
    )


def _result_from_scored_assembly(
    spec: DesignSpec,
    atoms: List[Dict],
    cores: List[Dict],
    second_shell_atoms: List[Dict],
    score: ConstraintScore,
    label: str,
) -> AssemblyResult:
    chemistry = annotate_chemical_system(atoms, spec)
    validation = validate_assembly(
        {"atoms": atoms, "cores": cores},
        spec,
        stage="post_mace" if score.details.get("structure_relaxation", {}).get("status") == "success" else "preflight",
        formal_charge=chemistry.formal_charge,
        spin_multiplicities=chemistry.spin_multiplicities,
    )
    score.details["physchem_validation"] = validation.to_dict()
    if not validation.passed:
        score.passed_hard_constraints = False
        score.total_score = 0.0
        score.errors.extend(
            issue.message for issue in validation.all_issues
            if issue.severity == "error" and issue.message not in score.errors
        )
    return AssemblyResult(
        job_id=str(uuid.uuid4())[:8], atoms=atoms, cores=cores,
        second_shell_atoms=second_shell_atoms, score=score, design_spec=spec,
        bond_graph=chemistry.bond_graph,
        formal_charge=chemistry.formal_charge,
        spin_multiplicities=chemistry.spin_multiplicities,
        chemistry_warnings=chemistry.warnings,
        xyz=_to_xyz(atoms), label=label,
        error=None if score.passed_hard_constraints else "; ".join(score.errors[:3]),
    )

def _core_from_atoms(atoms: List[Dict], metal_spec, site_id: str) -> Dict:
    metal = next(
        a for a in atoms
        if a.get("site_id") == site_id and str(a["element"]).upper() == metal_spec.metal_type.upper()
    )
    coord_atoms = [
        a for a in atoms
        if a.get("site_id") == site_id and a.get("is_coord_atom")
    ]
    coord_atoms.extend(
        atom
        for atom in atoms
        if atom.get("is_bridge_atom")
        and float(
            np.linalg.norm(
                np.array(atom["coords"], dtype=float)
                - np.array(metal["coords"], dtype=float)
            )
        ) <= 2.55
    )
    return {
        "metal": metal,
        "coord_atoms": coord_atoms,
        "geometry": metal_spec.coordination_geometry,
        "metal_type": metal_spec.metal_type,
        "oxidation_state": metal_spec.oxidation_state,
        "site_id": site_id,
        "activity_type": metal_spec.activity_type,
    }


def _apply_local_doping(atoms: List[Dict], center_2d: np.ndarray, doping: str) -> None:
    counts = {"N": 2, "S": 1}
    for elem in ("N", "S"):
        if elem not in doping:
            continue
        candidates = sorted(
            (
                (float(np.linalg.norm(np.array(a["coords"][:2], dtype=float) - center_2d)), i)
                for i, a in enumerate(atoms)
                if str(a.get("element", "")).upper() == "C"
                and str(a.get("residue_name", "")).upper() == "GRA"
                and float(np.linalg.norm(np.array(a["coords"][:2], dtype=float) - center_2d)) >= 3.2
            ),
            key=lambda item: item[0],
        )
        start = min(6, max(0, len(candidates) - counts[elem]))
        for _, idx in candidates[start:start + counts[elem]]:
            atoms[idx]["element"] = elem
            atoms[idx]["residue_name"] = f"{elem}DP"
            atoms[idx]["atom_name"] = f"{elem}{idx}"


def _drop_support_clashes(
    existing_atoms: List[Dict],
    new_atoms: List[Dict],
    cutoff: float = SUPPORT_CLASH_CUTOFF_DEFAULT,
) -> List[Dict]:
    if not new_atoms:
        return existing_atoms
    new_positions = np.array([a["coords"] for a in new_atoms], dtype=float)
    kept = []
    for atom in existing_atoms:
        if not _is_support_atom(atom):
            kept.append(atom)
            continue
        if atom.get("is_coord_atom") or atom.get("is_bridge_atom"):
            kept.append(atom)
            continue
        pos = np.array(atom["coords"], dtype=float)
        distances = np.linalg.norm(new_positions - pos, axis=1)
        if len(distances) == 0 or float(np.min(distances)) >= cutoff:
            kept.append(atom)
    return kept


def _is_support_atom(atom: Dict) -> bool:
    residue = str(atom.get("residue_name", "")).upper()
    return residue == "GRA" or residue.endswith("DP")


def _remove_isolated_support_artifacts(atoms: List[Dict]) -> List[Dict]:
    """Drop support atoms left as isolated islands after clash pruning.

    Graphene clash pruning can occasionally leave a single support C/N/S atom
    that is then capped with hydrogens. It has bonds, but only to H, so the
    viewer shows it as a chemically misleading floating dot. Keep ordinary edge
    carbons with at least one heavy support neighbor and drop only true islands.
    """
    support_indices = [
        idx
        for idx, atom in enumerate(atoms)
        if _is_support_atom(atom) and str(atom.get("element", "")).upper() != "H"
    ]
    if not support_indices:
        return [dict(atom) for atom in atoms]

    positions = {idx: np.array(atoms[idx]["coords"], dtype=float) for idx in support_indices}
    remove = set()
    for idx in support_indices:
        atom = atoms[idx]
        if atom.get("is_coord_atom") or atom.get("is_bridge_atom"):
            continue
        support_heavy_neighbors = sum(
            1
            for other_idx in support_indices
            if other_idx != idx
            and 1.15 <= float(np.linalg.norm(positions[idx] - positions[other_idx])) <= 1.80
        )
        if support_heavy_neighbors == 0:
            remove.add(idx)
    if not remove:
        return [dict(atom) for atom in atoms]
    return [dict(atom) for idx, atom in enumerate(atoms) if idx not in remove]


def _attach_bimetallic_details(
    score: ConstraintScore,
    topology: BimetallicTopology,
    metal_centers: List[np.ndarray],
) -> None:
    details = score_metal_distance(topology, metal_centers)
    score.details["bimetallic"] = details
    distance_score = float(details.get("distance_score") or 0.0)
    score.geometry_score = 0.75 * score.geometry_score + 0.25 * distance_score
    if score.passed_hard_constraints:
        score.total_score = _weighted_total(score)


def _weighted_total(score: ConstraintScore) -> float:
    return float(
        TOTAL_SCORE_WEIGHTS["geometry"] * score.geometry_score
        + TOTAL_SCORE_WEIGHTS["coordination"] * score.coordination_score
        + TOTAL_SCORE_WEIGHTS["energy"] * score.energy_score
        + TOTAL_SCORE_WEIGHTS["steric"] * score.steric_score
    )


def _maybe_relax_atoms(atoms: List[Dict]) -> tuple[List[Dict], Dict]:
    config = PotentialEvaluationConfig.from_env()
    if config.backend not in {"mace", "fairchem"}:
        return atoms, {
            "status": "skipped",
            "backend": config.backend,
            "reason": "set E2N_MLP_BACKEND=mace or fairchem for ASE relaxation",
        }
    if not config.relax:
        return atoms, {
            "status": "skipped",
            "backend": config.backend,
            "reason": "set E2N_MLP_RELAX=1 to write relaxed coordinates to output",
        }
    try:
        return relax_atoms_with_ml_backend(atoms, config)
    except Exception as exc:
        return atoms, {
            "status": "failed",
            "backend": config.backend,
            "error": str(exc),
        }


def _score_config_after_relax(relax_details: Dict) -> PotentialEvaluationConfig:
    config = PotentialEvaluationConfig.from_env()
    if relax_details.get("status") == "success":
        from dataclasses import replace

        return replace(config, relax=False)
    return config


def _attach_relaxation_details(score: ConstraintScore, relax_details: Dict) -> None:
    score.details["structure_relaxation"] = relax_details
    status = relax_details.get("status")
    if status == "success":
        score.backend = relax_details.get("backend") or score.backend
        score.method = f"{score.backend}_relaxed_export"
        score.raw_energy_ev = relax_details.get("initial_energy_ev")
        score.relaxed_energy_ev = relax_details.get("relaxed_energy_ev")
        score.max_force_ev_per_a = relax_details.get("relaxed_max_force_ev_per_a")
        if "relaxed_energy_ev" in relax_details and len(score.details.get("ml", {})) == 0:
            score.details["ml"] = {
                "energy_ev": relax_details.get("relaxed_energy_ev"),
                "max_force_ev_per_a": relax_details.get("relaxed_max_force_ev_per_a"),
            }
        score.warnings.append(
            f"Exported coordinates were relaxed with {relax_details.get('backend')} "
            f"({relax_details.get('steps_run')} steps, converged={relax_details.get('relax_converged')})."
        )
        if relax_details.get("relaxation_status") == "rejected":
            score.passed_hard_constraints = False
            score.total_score = 0.0
            score.errors.append(
                "MACE relaxation rejected: raw maximum force exceeds the usable threshold."
            )
    elif status in {"failed", "skipped"}:
        reason = relax_details.get("error") or relax_details.get("reason")
        score.warnings.append(f"Exported coordinates use local clash cleanup only; ML relaxation {status}: {reason}")


def _to_xyz(atoms: List[Dict]) -> str:
    lines = [str(len(atoms)), "nanozyme structure"]
    for a in atoms:
        x, y, z = a["coords"]
        lines.append(f"{a['element']:<4} {x:10.4f} {y:10.4f} {z:10.4f}")
    return "\n".join(lines)
