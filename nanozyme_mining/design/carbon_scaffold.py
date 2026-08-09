"""
Build graphene-supported metal coordination motifs.

When buildamol is installed we use its MetalComplexer for the first shell.
Otherwise the fallback must still honor the requested coordination template and
produce a non-overlapping initial geometry for downstream relaxation.
"""
import itertools
import math
import numpy as np
from typing import List, Dict, Optional, Sequence
from ..utils.constants import ALL_METAL_ELEMENTS, CATALYTIC_METAL_ELEMENTS  # PR4-1 (M12/M13)

_CC = 1.42

# PR1-6 (H9): one-time warning flag for missing buildamol so we don't spam logs.
_BUILDAMOL_MISSING_WARNED = False

# Fallback defaults used only when the UI/database template does not provide CN.
#
# PR1-2 (v4 audit) — corrected to biologically/chemically reasonable defaults:
# - FE: was square_planar+cn=4 (rare for d5/d6 Fe). Octahedral CN=6 is the dominant
#   coordination for heme-Fe and most non-heme Fe enzymes (H8 fix).
# - MN: was tetrahedral+cn=4. Mn(III/IV) catalases are octahedral with bridging
#   oxos (M1/M2 fix). For octahedral-default we use cn=6, bl≈2.10 Å.
# - CO/NI: tightened CN/geom for typical d7/d8 cases.
# - PD/PT/AU: added for completeness (F-M2 fix) — d8 noble metals are square planar.
# - W/AG: added because dopant_modifier._METALS already references them.
_METAL_DEFAULTS = {
    "FE": {"cn": 6, "bl": 2.10, "geom": "octahedral"},      # H8 fix: heme-like default
    "CU": {"cn": 4, "bl": 2.00, "geom": "square_planar"},   # d9 Jahn-Teller → SP common
    "ZN": {"cn": 4, "bl": 2.10, "geom": "tetrahedral"},     # d10 closed-shell → Td common
    "MN": {"cn": 6, "bl": 2.10, "geom": "octahedral"},      # M1 fix: Mn-catalase Oh
    "CO": {"cn": 6, "bl": 2.00, "geom": "octahedral"},      # Co(II/III) Oh common
    "NI": {"cn": 4, "bl": 2.05, "geom": "square_planar"},   # Ni(II) d8 → SP common
    "MO": {"cn": 6, "bl": 2.10, "geom": "Octahedral",  "eq": 6},
    "V":  {"cn": 5, "bl": 2.10, "geom": "TrigonalBipyramidal", "eq": 5},
    "CR": {"cn": 6, "bl": 2.05, "geom": "Octahedral",  "eq": 6},
    "RU": {"cn": 6, "bl": 2.05, "geom": "Octahedral",  "eq": 6},
    # F-M2: d8 noble metals; also referenced by dopant_modifier._METALS
    "PD": {"cn": 4, "bl": 2.00, "geom": "square_planar"},
    "PT": {"cn": 4, "bl": 2.00, "geom": "square_planar"},
    "AU": {"cn": 4, "bl": 2.05, "geom": "square_planar"},
    "AG": {"cn": 2, "bl": 2.15, "geom": "linear"},          # Ag(I) d10 → linear typical
    "W":  {"cn": 6, "bl": 2.15, "geom": "Octahedral",  "eq": 6},
}

_BUILDAMOL_GEOMETRY = {
    "octahedral": "Octahedral",
    "tetrahedral": "Tetrahedral",
    "trigonal_bipyramidal": "TrigonalBipyramidal",
}

_DONOR_BY_RESIDUE = {
    "HIS": ("N", "NE2"),
    "CYS": ("S", "SG"),
    "MET": ("S", "SD"),
    "ASP": ("O", "OD1"),
    "GLU": ("O", "OE1"),
    "TYR": ("O", "OH"),
    "SER": ("O", "OG"),
    "THR": ("O", "OG1"),
    "LYS": ("N", "NZ"),
    "ARG": ("N", "NH1"),
    "HOH": ("O", "O"),
    "OH": ("O", "O"),
}

_TYPICAL_BOND_LENGTHS = {
    ("FE", "N"): 2.05, ("FE", "O"): 1.95, ("FE", "S"): 2.30,
    ("CU", "N"): 2.00, ("CU", "O"): 1.95, ("CU", "S"): 2.20,
    ("ZN", "N"): 2.10, ("ZN", "O"): 2.00, ("ZN", "S"): 2.30,
    ("MN", "N"): 2.15, ("MN", "O"): 1.90, ("MN", "S"): 2.40,
    ("CO", "N"): 2.00, ("CO", "O"): 1.95, ("CO", "S"): 2.25,
    ("NI", "N"): 2.05, ("NI", "O"): 2.00, ("NI", "S"): 2.25,
}

_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07,
    "FE": 1.24, "CU": 1.17, "ZN": 1.25, "MN": 1.39, "CO": 1.18,
    "NI": 1.17, "MO": 1.54, "V": 1.53, "CR": 1.39, "RU": 1.46,
}

# 残基 → (SMILES, 配位原子ID, 删除的H的ID)
_LIGAND_SMILES = {
    "HIS": ("c1cnc[nH]1", "N2", "H4"),   # 咪唑，N2 配位
    "CYS": ("[SH]CC",     "S1", "H1"),   # 硫醇，S 配位
    "ASP": ("OC(=O)C",   "O1", "H1"),   # 羧基，O 配位
    "GLU": ("OC(=O)CC",  "O1", "H1"),
    "TYR": ("Oc1ccccc1", "O1", "H1"),   # 酚羟基
    "SER": ("OCC",       "O1", "H1"),
    "THR": ("OC(C)C",    "O1", "H1"),
    "MET": ("CSC",       "S2", "H1"),   # 硫醚
    "LYS": ("NCCCC",     "N1", "H1"),
    "ARG": ("NC(=N)N",   "N1", "H1"),
    "HOH": ("O",         "O1", "H1"),
    "OH":  ("[OH-]",     "O1", "H1"),
}
_DEFAULT_LIGAND = ("c1cnc[nH]1", "N2", "H4")  # 默认咪唑

_FALLBACK_LIGAND_SMILES = {
    "HIS": "c1ncc[nH]1",
    "CYS": "[S-]CC",
    "ASP": "CC(=O)[O-]",
    "GLU": "CCC(=O)[O-]",
    "TYR": "[O-]c1ccccc1",
    "SER": "OCC",
    "THR": "OC(C)C",
    "MET": "CSC",
    "LYS": "NCCCC",
    "ARG": "NC(=[NH2+])N",
    "HOH": "O",
    "OH": "[OH-]",
}

# Octahedral 赤道位 H 标签（H1-H4）和轴向位（H5-H6）
_OCTAHEDRAL_EQ  = ["H1", "H2", "H3", "H4"]
_OCTAHEDRAL_AX  = ["H5", "H6"]
_TETRAHEDRAL_ALL = ["H1", "H2", "H3", "H4"]
_TRIGBIP_ALL    = ["H1", "H2", "H3", "H4", "H5"]


def build_graphene_fragment(radius: float = 6.5) -> List[Dict]:
    """生成以原点为中心的石墨烯蜂窝网格（z=0平面）。"""
    d = _CC
    a1 = np.array([np.sqrt(3) * d, 0.0])
    a2 = np.array([np.sqrt(3) * d / 2, 3 * d / 2])
    atoms, seen = [], set()
    R = int(radius / d) + 3
    for n in range(-R, R + 1):
        for m in range(-R, R + 1):
            for sub in (np.zeros(2), np.array([0.0, d])):
                pos = n * a1 + m * a2 + sub
                if np.linalg.norm(pos) > radius:
                    continue
                key = (round(pos[0], 3), round(pos[1], 3))
                if key in seen:
                    continue
                seen.add(key)
                atoms.append({"element": "C", "residue_name": "GRA",
                               "atom_name": f"C{len(atoms)+1}",
                               "coords": [pos[0], pos[1], 0.0]})
    return atoms


def embed_metal_in_graphene(
    graphene: List[Dict],
    metal_type: str,
    coord_residues: List[str],
    metal_offset: np.ndarray = None,
    doping: str = "none",
    n_dope: int = 2,
    s_dope: int = 1,
    coordination_number: Optional[int] = None,
    coordination_geometry: Optional[str] = None,
    coord_atoms: Optional[List] = None,
    site_id: Optional[str] = None,
) -> List[Dict]:
    """
    用 buildamol MetalComplexer 生成真实 3D 金属配合物，
    叠加到石墨烯片段上。
    """
    if metal_offset is None:
        metal_offset = np.zeros(3)

    site_id = site_id or "M0"

    # 用 buildamol 生成金属配合物核心；不可用时使用模板感知 fallback。
    complex_atoms = _buildamol_metal_complex(
        metal_type,
        coord_residues,
        metal_offset,
        coordination_number=coordination_number,
        coordination_geometry=coordination_geometry,
        coord_atoms=coord_atoms,
        site_id=site_id,
    )

    # 石墨烯背景：按实际距离移除与配合物冲突的 C，而不是按 XY 四舍五入。
    graphene = [dict(a) for a in graphene]
    graphene = _remove_graphene_clashes(graphene, complex_atoms)

    # N/S 掺杂
    center_2d = metal_offset[:2]
    if doping in ("N", "NS"):
        _apply_doping(graphene, center_2d, "N", n_dope, shell=2)
    if doping in ("S", "NS"):
        _apply_doping(graphene, center_2d, "S", s_dope, shell=3)

    return passivate_graphene_edges(graphene + complex_atoms)


def embed_metal_in_carbon_network(
    graphene: List[Dict],
    metal_type: str,
    coord_residues: List[str],
    metal_offset: np.ndarray = None,
    coordination_number: Optional[int] = None,
    coordination_geometry: Optional[str] = None,
    coord_atoms: Optional[List] = None,
    site_id: Optional[str] = None,
    reserve_coordination_slots: int = 0,
    reserved_direction_xy: Optional[np.ndarray] = None,
) -> List[Dict]:
    """Embed a metal into an existing graphene pocket.

    Unlike ``embed_metal_in_graphene`` this does not build a standalone metal
    complex. The first-shell N/O/S donors are converted from atoms already in
    the carbon support, so the initial graph handed to ML relaxation is a
    continuous heteroatom-doped carbon network with an embedded metal site.
    """
    if metal_offset is None:
        metal_offset = np.zeros(3)
    site_id = site_id or "M0"

    template = _coordination_template(
        metal_type,
        coord_residues,
        coordination_number,
        coordination_geometry,
        coord_atoms,
    )
    atoms = [dict(atom) for atom in graphene]
    donors = template["donors"][:template["cn"]]
    metal_xy = np.array(metal_offset[:2], dtype=float)

    layout = _network_geometry_layout(
        template["geometry"],
        donors,
        reserve_coordination_slots=reserve_coordination_slots,
        reserved_direction_xy=reserved_direction_xy,
    )
    donor_indices, donor_targets, fit = _fit_network_donor_pocket(
        atoms,
        metal_xy,
        layout["network_donors"],
        layout["network_vectors_xy"],
    )
    if fit["max_reposition_a"] > 1.25:
        raise ValueError(
            f"not_constructible:coordination_pocket_fit:{site_id}:"
            f"max_reposition={fit['max_reposition_a']:.3f}A"
        )
    _tag_network_donors(
        atoms,
        donor_indices,
        layout["network_donors"],
        site_id,
        target_positions=donor_targets,
    )
    atoms = _open_graphene_pocket(atoms, metal_xy, donor_indices)

    metal_coords = np.array(
        [float(metal_xy[0]), float(metal_xy[1]), float(layout["metal_z_a"])],
        dtype=float,
    )
    metal_atom = {
        "element": metal_type,
        "residue_name": "MET",
        "atom_name": metal_type,
        "coords": metal_coords.tolist(),
        "site_id": site_id,
        "formal_charge": 0,
        "is_embedded_metal": True,
        "coordination_pocket_fit": fit,
    }
    explicit_atoms = _build_explicit_coordination_donors(
        metal_coords,
        layout["explicit_donors"],
        layout["explicit_vectors"],
        site_id,
    )
    return atoms + [metal_atom] + explicit_atoms


def embed_bridge_in_carbon_network(
    atoms: List[Dict],
    topology,
    metal_centers: Sequence[np.ndarray],
) -> List[Dict]:
    """Convert existing support atoms into bimetallic bridge atoms."""
    if len(metal_centers) < 2:
        return []

    m0 = np.array(metal_centers[0], dtype=float)
    m1 = np.array(metal_centers[1], dtype=float)
    axis = m1 - m0
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-8:
        return []
    axis = axis / axis_norm
    normal = _in_plane_perpendicular(axis)
    midpoint = (m0 + m1) / 2.0
    relation = getattr(topology.edge, "relation", "")

    if relation == "bridged":
        targets = [
            ("N", "N1", m0 + axis * 2.0 + normal * 0.65, True),
            # Keep the terminal N-C links above the strict covalent-clearance
            # threshold.  A 3.0 A axial offset made both links exactly 1.00 A
            # for the default 5.8 A bridge, so preflight rejected the topology
            # before MACE had a chance to relax it.
            ("C", "C1", m0 + axis * 3.1 + normal * 0.65, False),
            ("C", "C2", m1 - axis * 3.1 - normal * 0.65, False),
            ("N", "N2", m1 - axis * 2.0 - normal * 0.65, True),
        ]
        tagged_atoms = []
        used = set()
        for element, atom_name, target, is_donor in targets:
            tagged = _tag_nearest_support_atom(
                atoms,
                target[:2],
                element=element,
                atom_name=atom_name,
                bridge_role="diimine_bridge" if is_donor else "bridge_backbone",
                bridge_relation=relation,
                exclude_coord=True,
                used_atom_names=used,
                mark_as_bridge=is_donor,
                reposition=True,
            )
            if tagged:
                used.add(tagged.get("atom_name"))
                tagged_atoms.append(tagged)
        return tagged_atoms

    return []


def _buildamol_metal_complex(
    metal_type: str,
    coord_residues: List[str],
    offset: np.ndarray,
    coordination_number: Optional[int] = None,
    coordination_geometry: Optional[str] = None,
    coord_atoms: Optional[List] = None,
    site_id: Optional[str] = None,
) -> List[Dict]:
    """用 buildamol 生成真实 3D 金属配合物，返回原子列表。

    PR1-6 (H9 fix): the catch-all `except Exception` block at the bottom used
    to silently fall back to the simple geometry for ANY error — including the
    case where buildamol is not installed. That made install issues invisible
    and made debugging legitimate API failures painful. Now ImportError
    surfaces a single one-time warning (and uses the fallback); other runtime
    exceptions are logged with their type so a regression in buildamol's API
    is visible in logs.
    """
    try:
        import buildamol as bam
        from buildamol.extensions import complexes
        import buildamol.structural as structural

        template = _coordination_template(
            metal_type,
            coord_residues,
            coordination_number,
            coordination_geometry,
            coord_atoms,
        )
        cn = template["cn"]
        bl = template["default_bond_length"]
        geom_key = _normalize_geometry(template["geometry"])
        geom_name = _BUILDAMOL_GEOMETRY.get(geom_key)
        if geom_name is None:
            raise RuntimeError(f"buildamol geometry not available for {geom_key}")

        metal_atom = bam.Atom.new(metal_type, pqr_charge=2)
        geom = getattr(structural.geometry, geom_name)()
        complexer = complexes.MetalComplexer(metal_atom, geom)
        complexer.make_core(bl)

        # 确定配体和 acceptor H 位置
        if geom_name == "Tetrahedral":
            acceptor_slots = _TETRAHEDRAL_ALL[:cn]
        elif geom_name == "TrigonalBipyramidal":
            acceptor_slots = _TRIGBIP_ALL[:cn]
        else:  # Octahedral
            acceptor_slots = (_OCTAHEDRAL_EQ + _OCTAHEDRAL_AX)[:cn]

        # 填充配体到 cn 个位置
        donors = template["donors"][:cn]

        ligands, binders, acceptors, deletes = [], [], [], []
        for i, donor in enumerate(donors):
            smi, bind_id, del_id = _LIGAND_SMILES.get(donor["residue_name"].upper(), _DEFAULT_LIGAND)
            lig = bam.Molecule.from_smiles(smi, id=f"{donor['residue_name']}{i}")
            ligands.append(lig)
            binders.append([bind_id])
            acceptors.append([acceptor_slots[i]])
            deletes.append([del_id])

        result = complexer.add_ligands(
            ligands=ligands,
            binders=binders,
            acceptors=acceptors,
            delete=deletes,
            optimize=False,
        )

        # buildamol MMFF 优化
        import buildamol.optimizers as bam_opt
        try:
            result = bam_opt.mmff_optimize(result, steps=500)
        except Exception:
            pass

        # 转为原子列表，平移到 offset。显式 H 对后续 tblite/xTB 是必需的。
        atoms = []
        for a in result.get_atoms():
            coord = np.array(a.coord) + offset
            res_name = "MET" if a.element == metal_type else a.parent.resname if hasattr(a, 'parent') else "LIG"
            atoms.append({
                "element": a.element,
                "residue_name": res_name,
                "atom_name": a.id,
                "coords": coord.tolist(),
                "site_id": site_id,
            })

        # 标记配位 N/O/S 原子（NVA）
        metal_pos = offset.copy()
        for a in atoms:
            if a["element"] in ("N", "O", "S"):
                dist = np.linalg.norm(np.array(a["coords"]) - metal_pos)
                if dist < bl * 1.15:
                    a["residue_name"] = "NVA"
                    a["is_coord_atom"] = True

        return atoms

    except ImportError as e:
        # PR1-6 (H9 fix): one-time WARNING if buildamol isn't installed. Use
        # module-global flag so we don't spam the log per-call.
        global _BUILDAMOL_MISSING_WARNED
        if not _BUILDAMOL_MISSING_WARNED:
            import logging
            logging.getLogger("e2n.carbon_scaffold").warning(
                "buildamol not installed (%s); using simple-geometry fallback. "
                "Install with: uv pip install --python .venv-mace/bin/python buildamol",
                e,
            )
            _BUILDAMOL_MISSING_WARNED = True
        return _fallback_metal_complex(
            metal_type,
            coord_residues,
            offset,
            coordination_number=coordination_number,
            coordination_geometry=coordination_geometry,
            coord_atoms=coord_atoms,
            site_id=site_id,
        )
    except Exception as e:
        # PR1-6 (H9 fix): log non-ImportError exceptions with type so a buildamol
        # API regression is visible. Continue with fallback so the design pipeline
        # still produces output.
        import logging
        logging.getLogger("e2n.carbon_scaffold").warning(
            "buildamol failed (%s: %s); using simple-geometry fallback",
            type(e).__name__, e,
        )
        return _fallback_metal_complex(
            metal_type,
            coord_residues,
            offset,
            coordination_number=coordination_number,
            coordination_geometry=coordination_geometry,
            coord_atoms=coord_atoms,
            site_id=site_id,
        )


def _fallback_metal_complex(
    metal_type: str,
    coord_residues: List[str],
    offset: np.ndarray,
    coordination_number: Optional[int] = None,
    coordination_geometry: Optional[str] = None,
    coord_atoms: Optional[List] = None,
    site_id: Optional[str] = None,
) -> List[Dict]:
    """Template-aware fallback when buildamol is unavailable."""
    template = _coordination_template(
        metal_type,
        coord_residues,
        coordination_number,
        coordination_geometry,
        coord_atoms,
    )
    cn = template["cn"]
    vectors = _coordination_vectors(template["geometry"], cn)
    atoms = [{"element": metal_type, "residue_name": "MET",
               "atom_name": metal_type, "coords": offset.tolist(),
               "site_id": site_id, "formal_charge": 0}]
    for i, donor in enumerate(template["donors"][:cn]):
        u = vectors[i]
        bl = donor["bond_length"] or template["default_bond_length"]
        donor_pos = offset + u * bl
        try:
            atoms.extend(_fallback_ligand_fragment(donor, donor_pos, u, i, site_id))
        except Exception:
            atoms.append({
                "element": donor["element"],
                "residue_name": "NVA",
                "source_residue_name": donor["residue_name"],
                "atom_name": donor["atom_name"],
                "coords": donor_pos.tolist(),
                "site_id": site_id,
                "is_coord_atom": True,
                "formal_charge": 0,
            })
            atoms.extend(_fallback_ligand_stub(donor, donor_pos, u, i, site_id))
    return atoms


def _apply_doping(graphene: List[Dict], center: np.ndarray,
                  elem: str, count: int, shell: int):
    candidates = sorted(
        [(np.linalg.norm(np.array(a["coords"][:2]) - center), i)
         for i, a in enumerate(graphene) if a["element"] == "C"]
    )
    start = min(shell * 6, max(0, len(candidates) - count))
    for _, i in candidates[start: start + count]:
        graphene[i]["element"] = elem
        graphene[i]["residue_name"] = f"{elem}DP"
        graphene[i]["atom_name"] = f"{elem}{i}"


def _select_network_donor_indices(
    atoms: List[Dict],
    metal_xy: np.ndarray,
    donors: List[Dict],
) -> List[int]:
    support = [
        (idx, atom, np.array(atom["coords"][:2], dtype=float))
        for idx, atom in enumerate(atoms)
        if _is_network_support(atom)
        and str(atom.get("element", "")).upper() != "H"
        and not atom.get("is_coord_atom")
        and not atom.get("is_bridge_atom")
    ]
    if not support:
        return []

    directions = _network_donor_directions(len(donors))
    selected: List[int] = []
    selected_set = set()
    for direction, donor in zip(directions, donors):
        donor_element = str(donor.get("element", "N")).upper()
        target_radius = _network_donor_projected_radius(donor)
        target_xy = metal_xy + direction * target_radius
        best = None
        best_score = float("inf")
        for idx, atom, pos in support:
            if idx in selected_set:
                continue
            radial = float(np.linalg.norm(pos - metal_xy))
            if radial < 0.95 or radial > 2.65:
                continue
            existing_element = str(atom.get("element", "")).upper()
            element_bonus = -0.30 if existing_element == donor_element else 0.0
            radial_penalty = 0.42 * abs(radial - target_radius)
            coord_penalty = float(np.linalg.norm(pos - target_xy))
            assigned_penalty = 1.25 if atom.get("site_id") else 0.0
            score = coord_penalty + radial_penalty + assigned_penalty + element_bonus
            if score < best_score:
                best_score = score
                best = idx
        if best is None:
            remaining = [
                (
                    abs(float(np.linalg.norm(pos - metal_xy)) - target_radius)
                    + 0.5 * float(np.linalg.norm(pos - target_xy)),
                    idx,
                )
                for idx, _atom, pos in support
                if idx not in selected_set
            ]
            if not remaining:
                continue
            best = min(remaining, key=lambda item: item[0])[1]
        selected.append(best)
        selected_set.add(best)
    return selected


def _fit_network_donor_pocket(
    atoms: List[Dict],
    metal_xy: np.ndarray,
    donors: List[Dict],
    vectors_xy: np.ndarray,
) -> tuple[List[int], List[np.ndarray], Dict]:
    """Jointly assign a complete donor pocket to ideal geometry targets."""
    if not donors:
        return [], [], {"rms_reposition_a": 0.0, "max_reposition_a": 0.0}
    support = [
        (idx, atom, np.asarray(atom["coords"], dtype=float))
        for idx, atom in enumerate(atoms)
        if _is_network_support(atom)
        and str(atom.get("element", "")).upper() != "H"
        and not atom.get("is_coord_atom")
        and not atom.get("is_bridge_atom")
    ]
    if len(support) < len(donors):
        raise ValueError("not_constructible:insufficient_support_atoms_for_coordination_pocket")

    targets = [
        np.array(
            [
                metal_xy[0] + vectors_xy[index][0] * float(donor["projected_radius_a"]),
                metal_xy[1] + vectors_xy[index][1] * float(donor["projected_radius_a"]),
                0.0,
            ],
            dtype=float,
        )
        for index, donor in enumerate(donors)
    ]
    candidate_ids = set()
    for target in targets:
        nearest = sorted(
            (
                (float(np.linalg.norm(position[:2] - target[:2])), idx)
                for idx, _atom, position in support
            ),
            key=lambda item: item[0],
        )[:10]
        candidate_ids.update(idx for _distance, idx in nearest)
    candidate_ids = set(
        idx
        for _score, idx in sorted(
            (
                (
                    min(float(np.linalg.norm(position[:2] - target[:2])) for target in targets),
                    idx,
                )
                for idx, _atom, position in support
                if idx in candidate_ids
            ),
            key=lambda item: item[0],
        )[:14]
    )
    positions = {idx: position for idx, _atom, position in support if idx in candidate_ids}
    elements = {idx: str(atom.get("element", "")).upper() for idx, atom, _position in support}

    best = None
    best_score = float("inf")
    for assignment in itertools.permutations(sorted(candidate_ids), len(donors)):
        displacements = [
            float(np.linalg.norm(positions[idx] - target))
            for idx, target in zip(assignment, targets)
        ]
        if max(displacements) > 1.35:
            continue
        element_penalty = sum(
            0.0 if elements[idx] == str(donor["element"]).upper() else 0.08
            for idx, donor in zip(assignment, donors)
        )
        score = sum(value * value for value in displacements) + element_penalty
        if score < best_score:
            best_score = score
            best = assignment, displacements
    if best is None:
        raise ValueError("not_constructible:no_joint_coordination_pocket_assignment")
    assignment, displacements = best
    return list(assignment), targets, {
        "method": "joint_assignment_with_local_pocket_reposition",
        "rms_reposition_a": float(math.sqrt(sum(d * d for d in displacements) / len(displacements))),
        "max_reposition_a": float(max(displacements)),
    }


def _tag_network_donors(
    atoms: List[Dict],
    donor_indices: List[int],
    donors: List[Dict],
    site_id: str,
    target_positions: Optional[List[np.ndarray]] = None,
) -> None:
    for order, (idx, donor) in enumerate(zip(donor_indices, donors), 1):
        element = str(donor.get("element", "N")).upper()
        atom = atoms[idx]
        atom["element"] = element
        atom["residue_name"] = f"{element}DP"
        atom["source_residue_name"] = donor.get("residue_name")
        atom["source_atom_name"] = donor.get("atom_name")
        atom["atom_name"] = f"{element}{site_id}_{order}"
        atom["site_id"] = site_id
        atom["is_coord_atom"] = True
        atom["is_network_donor"] = True
        atom["bond_length"] = float(donor.get("bond_length") or 0.0)
        atom["bond_length_range"] = list(donor.get("bond_length_range") or []) or None
        atom["coordination_role"] = donor.get("role", "equatorial_network")
        atom["labile"] = bool(donor.get("labile", False))
        atom["protonation_state"] = donor.get("protonation_state")
        atom["source_id"] = donor.get("source_id")
        if target_positions is not None and order - 1 < len(target_positions):
            atom["coords"] = [float(value) for value in target_positions[order - 1]]
        atom.setdefault("formal_charge", 0)


def _open_graphene_pocket(
    atoms: List[Dict],
    metal_xy: np.ndarray,
    donor_indices: List[int],
) -> List[Dict]:
    donor_names = {
        atoms[idx].get("atom_name")
        for idx in donor_indices
        if 0 <= idx < len(atoms)
    }
    mandatory_remove = []
    central_candidates = []
    donor_positions = [
        np.asarray(atoms[idx]["coords"], dtype=float)
        for idx in donor_indices
        if 0 <= idx < len(atoms)
    ]
    for idx, atom in enumerate(atoms):
        if atom.get("atom_name") in donor_names:
            continue
        if not _is_network_support(atom):
            continue
        if str(atom.get("element", "")).upper() == "H":
            continue
        dist = float(np.linalg.norm(np.array(atom["coords"][:2], dtype=float) - metal_xy))
        donor_clearance = min(
            (float(np.linalg.norm(np.asarray(atom["coords"], dtype=float) - donor_pos)) for donor_pos in donor_positions),
            default=float("inf"),
        )
        if donor_clearance < 1.10:
            mandatory_remove.append(idx)
        elif dist <= 0.92:
            central_candidates.append((dist, idx))

    if not mandatory_remove and not central_candidates:
        candidates = [
            (
                float(np.linalg.norm(np.array(atom["coords"][:2], dtype=float) - metal_xy)),
                idx,
            )
            for idx, atom in enumerate(atoms)
            if _is_network_support(atom)
            and str(atom.get("element", "")).upper() == "C"
            and atom.get("atom_name") not in donor_names
        ]
        candidates = [item for item in candidates if item[0] <= 1.22]
        central_candidates = sorted(candidates)[:1]
    remove = set(mandatory_remove)
    remove.update(idx for _dist, idx in sorted(central_candidates)[:2])
    if not remove:
        return [dict(atom) for atom in atoms]
    return [dict(atom) for idx, atom in enumerate(atoms) if idx not in remove]


def _network_metal_position(
    metal_xy: np.ndarray,
    donor_atoms: List[Dict],
    donor_specs: List[Dict],
) -> np.ndarray:
    if not donor_atoms:
        return np.array([float(metal_xy[0]), float(metal_xy[1]), 1.35], dtype=float)

    xy = np.array([atom["coords"][:2] for atom in donor_atoms], dtype=float)
    targets = np.array(
        [
            float(spec.get("bond_length") or 2.05)
            for spec in donor_specs[:len(donor_atoms)]
        ],
        dtype=float,
    )
    # Keep the metal over the intended topology node, but nudge it toward the
    # donor ring centroid if the discrete graphene lattice is off-center.
    centroid = np.mean(xy, axis=0)
    trial_centers = [
        metal_xy,
        0.85 * metal_xy + 0.15 * centroid,
        0.70 * metal_xy + 0.30 * centroid,
        0.50 * metal_xy + 0.50 * centroid,
        centroid,
    ]
    for radius in (0.18, 0.36):
        for direction in _network_donor_directions(6):
            trial_centers.append(0.70 * metal_xy + 0.30 * centroid + direction * radius)

    best_xy = trial_centers[0]
    best_z = 1.35
    best_score = float("inf")
    for trial_xy in trial_centers:
        radial = np.linalg.norm(xy - trial_xy, axis=1)
        heights = np.sqrt(np.maximum(targets * targets - radial * radial, 0.72 * 0.72))
        z = float(np.clip(np.mean(heights), 0.90, 1.85))
        distances = np.sqrt(radial * radial + z * z)
        deviations = np.abs(distances - targets)
        score = (
            float(np.mean(deviations))
            + 1.7 * float(np.max(deviations))
            + 0.10 * float(np.linalg.norm(trial_xy - metal_xy))
        )
        if score < best_score:
            best_score = score
            best_xy = trial_xy
            best_z = z
    return np.array([float(best_xy[0]), float(best_xy[1]), float(best_z)], dtype=float)


def _network_donor_directions(count: int) -> np.ndarray:
    count = max(1, int(count))
    if count == 1:
        angles = [0.0]
    elif count == 2:
        angles = [0.0, math.pi]
    elif count == 3:
        angles = [0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0]
    elif count == 4:
        angles = [0.0, math.pi / 2.0, math.pi, 3.0 * math.pi / 2.0]
    else:
        angles = [2.0 * math.pi * i / count for i in range(count)]
    return np.array([[math.cos(a), math.sin(a)] for a in angles], dtype=float)


def _network_donor_projected_radius(donor: Dict) -> float:
    bond_length = float(donor.get("bond_length") or 2.05)
    # A pyridinic graphene pore places donors about one C-C bond from the pore
    # center. The metal is above the sheet, so the projected radius is shorter
    # than the full metal-donor distance.
    return float(np.clip(0.68 * bond_length, 1.30, 1.58))


def _network_geometry_layout(
    geometry: str,
    donors: List[Dict],
    *,
    reserve_coordination_slots: int = 0,
    reserved_direction_xy: Optional[np.ndarray] = None,
) -> Dict:
    geometry = _normalize_geometry(geometry)
    all_donors = [dict(donor) for donor in donors]
    if geometry == "octahedral":
        network_count = min(4, len(all_donors))
        network_vectors = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], dtype=float)
        explicit_vectors = np.array([[0, 0, 1], [0, 0, -1]], dtype=float)
        metal_z = 0.0
    elif geometry == "square_pyramidal":
        network_count = min(4, len(all_donors))
        network_vectors = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], dtype=float)
        explicit_vectors = np.array([[0, 0, 1]], dtype=float)
        metal_z = 0.0
    elif geometry == "trigonal_bipyramidal":
        network_count = min(3, len(all_donors))
        network_vectors = _network_donor_directions(3)
        explicit_vectors = np.array([[0, 0, 1], [0, 0, -1]], dtype=float)
        metal_z = 0.0
    elif (
        geometry == "tetrahedral"
        and len(all_donors) >= 4
        and str(all_donors[-1].get("role", "")).lower() == "axial_labile"
    ):
        network_count = 3
        network_vectors = _network_donor_directions(3)
        explicit_vectors = np.array([[0, 0, 1]], dtype=float)
        metal_z = float(all_donors[0].get("bond_length") or 2.05) / 3.0
    else:
        network_count = len(all_donors)
        network_vectors = _network_donor_directions(max(network_count, 1))
        explicit_vectors = np.empty((0, 3), dtype=float)
        metal_z = 0.0

    network_donors = all_donors[:network_count]
    explicit_donors = all_donors[network_count:]
    for _ in range(min(int(reserve_coordination_slots), len(network_donors))):
        reserved = np.asarray(reserved_direction_xy if reserved_direction_xy is not None else [1.0, 0.0], dtype=float)
        reserved /= np.linalg.norm(reserved) or 1.0
        remove_index = int(np.argmax(network_vectors @ reserved))
        network_vectors = np.delete(network_vectors, remove_index, axis=0)
        del network_donors[remove_index]

    network_vectors = np.asarray(network_vectors[:len(network_donors)], dtype=float)
    if geometry == "tetrahedral" and len(network_donors) == 3:
        for donor in network_donors:
            donor["projected_radius_a"] = float(donor.get("bond_length") or 2.05) * math.sqrt(8.0 / 9.0)
    else:
        for donor in network_donors:
            donor["projected_radius_a"] = float(donor.get("bond_length") or 2.05)
    if len(explicit_vectors) < len(explicit_donors):
        explicit_vectors = _coordination_vectors(geometry, len(all_donors))[network_count:]
    return {
        "network_donors": network_donors,
        "network_vectors_xy": network_vectors,
        "explicit_donors": explicit_donors,
        "explicit_vectors": np.asarray(explicit_vectors[:len(explicit_donors)], dtype=float),
        "metal_z_a": metal_z,
    }


def _build_explicit_coordination_donors(
    metal_coords: np.ndarray,
    donors: List[Dict],
    vectors: np.ndarray,
    site_id: str,
) -> List[Dict]:
    atoms: List[Dict] = []
    for index, (donor, vector) in enumerate(zip(donors, vectors), 1):
        unit = np.asarray(vector, dtype=float)
        unit /= np.linalg.norm(unit) or 1.0
        distance = float(donor.get("bond_length") or 2.05)
        position = metal_coords + unit * distance
        residue = str(donor.get("residue_name") or "HOH").upper()
        atom = {
            "element": str(donor.get("element") or "O").upper(),
            "residue_name": residue,
            "atom_name": f"{donor.get('atom_name') or 'O'}{site_id}_AX{index}",
            "coords": [float(value) for value in position],
            "site_id": site_id,
            "is_coord_atom": True,
            "is_explicit_axial_ligand": True,
            "coordination_role": donor.get("role") or "axial_labile",
            "bond_length": distance,
            "bond_length_range": list(donor.get("bond_length_range") or []) or None,
            "labile": bool(donor.get("labile", True)),
            "protonation_state": donor.get("protonation_state") or ("hydroxo" if residue == "OH" else "aquo"),
            "source_id": donor.get("source_id"),
            "formal_charge": -1 if residue == "OH" else 0,
        }
        atoms.append(atom)
        if residue in {"HOH", "OH"}:
            hydrogen_count = 1 if residue == "OH" else 2
            tangent = np.array([1.0, 0.0, 0.0])
            if abs(float(np.dot(tangent, unit))) > 0.85:
                tangent = np.array([0.0, 1.0, 0.0])
            tangent -= unit * float(np.dot(tangent, unit))
            tangent /= np.linalg.norm(tangent) or 1.0
            for hydrogen_index in range(hydrogen_count):
                sign = -1.0 if hydrogen_index else 1.0
                h_pos = position + 0.82 * tangent * sign + 0.45 * unit
                atoms.append(
                    {
                        "element": "H",
                        "residue_name": residue,
                        "atom_name": f"H{site_id}_AX{index}_{hydrogen_index + 1}",
                        "coords": [float(value) for value in h_pos],
                        "site_id": site_id,
                        "support_parent": atom["atom_name"],
                        "formal_charge": 0,
                    }
                )
    return atoms


def _tag_nearest_support_atom(
    atoms: List[Dict],
    target_xy: np.ndarray,
    *,
    element: str,
    atom_name: str,
    bridge_role: str,
    bridge_relation: str,
    exclude_coord: bool = False,
    used_atom_names: Optional[set] = None,
    mark_as_bridge: bool = True,
    reposition: bool = False,
) -> Optional[Dict]:
    used_atom_names = used_atom_names or set()
    candidates = []
    for idx, atom in enumerate(atoms):
        if not _is_network_support(atom):
            continue
        if str(atom.get("element", "")).upper() == "H":
            continue
        if exclude_coord and atom.get("is_coord_atom"):
            continue
        if atom.get("atom_name") in used_atom_names:
            continue
        pos = np.array(atom["coords"][:2], dtype=float)
        candidates.append((float(np.linalg.norm(pos - target_xy)), idx))
    if not candidates:
        return None
    _dist, idx = min(candidates, key=lambda item: item[0])
    atom = atoms[idx]
    element = element.upper()
    atom["element"] = element
    atom["residue_name"] = f"{element}DP" if element != "C" else "GRA"
    atom["atom_name"] = atom_name
    atom["is_network_bridge"] = True
    atom["bridge_role"] = bridge_role
    atom["bridge_relation"] = bridge_relation
    if mark_as_bridge:
        atom["is_bridge_atom"] = True
        atom["is_network_donor"] = True
    if reposition:
        atom["coords"] = [float(target_xy[0]), float(target_xy[1]), 0.0]
    atom.setdefault("formal_charge", 0)
    return atom


def _in_plane_perpendicular(axis: np.ndarray) -> np.ndarray:
    axis = np.array([float(axis[0]), float(axis[1]), 0.0], dtype=float)
    if np.linalg.norm(axis) < 1e-8:
        return np.array([0.0, 1.0, 0.0], dtype=float)
    axis = axis / np.linalg.norm(axis)
    return np.array([-axis[1], axis[0], 0.0], dtype=float)


def _is_network_support(atom: Dict) -> bool:
    residue = str(atom.get("residue_name", "")).upper()
    return residue == "GRA" or residue.endswith("DP")


def _coordination_template(
    metal_type: str,
    coord_residues: List[str],
    coordination_number: Optional[int],
    coordination_geometry: Optional[str],
    coord_atoms: Optional[List],
) -> Dict:
    metal = metal_type.upper()
    default = _METAL_DEFAULTS.get(metal, {"cn": 4, "bl": 2.05, "geom": "tetrahedral"})
    cn = int(coordination_number or len(coord_atoms or []) or len(coord_residues or []) or default["cn"])
    cn = max(1, min(cn, 9))
    geometry = _normalize_geometry(coordination_geometry or default["geom"])
    donors = _donor_records(metal, coord_residues, coord_atoms, cn)
    return {
        "cn": cn,
        "geometry": geometry,
        "donors": donors,
        "default_bond_length": float(default["bl"]),
    }


def _donor_records(
    metal_type: str,
    coord_residues: List[str],
    coord_atoms: Optional[List],
    cn: int,
) -> List[Dict]:
    records = []
    if coord_atoms:
        for atom in coord_atoms:
            donor_element = str(getattr(atom, "donor_element", "") or _dict_get(atom, "donor_element", "N")).upper()
            residue_name = str(getattr(atom, "residue_name", "") or _dict_get(atom, "residue_name", "HIS")).upper()
            atom_name = str(getattr(atom, "atom_name", "") or _dict_get(atom, "atom_name", "") or _DONOR_BY_RESIDUE.get(residue_name, ("N", "NE2"))[1])
            bond_length = float(getattr(atom, "bond_length", 0.0) or _dict_get(atom, "bond_length", 0.0) or 0.0)
            records.append({
                "element": donor_element,
                "residue_name": residue_name,
                "atom_name": atom_name,
                "bond_length": bond_length or _typical_bond_length(metal_type, donor_element),
                "role": str(getattr(atom, "role", "") or _dict_get(atom, "role", "") or "equatorial_network"),
                "bond_length_range": getattr(atom, "bond_length_range", None) or _dict_get(atom, "bond_length_range"),
                "labile": bool(getattr(atom, "labile", False) or _dict_get(atom, "labile", False)),
                "protonation_state": getattr(atom, "protonation_state", None) or _dict_get(atom, "protonation_state"),
                "source_id": getattr(atom, "source_id", None) or _dict_get(atom, "source_id"),
            })
    else:
        for residue_name in coord_residues or []:
            residue_name = str(residue_name or "HIS").upper()
            donor_element, atom_name = _DONOR_BY_RESIDUE.get(residue_name, ("N", "NE2"))
            records.append({
                "element": donor_element,
                "residue_name": residue_name,
                "atom_name": atom_name,
                "bond_length": _typical_bond_length(metal_type, donor_element),
                "role": "equatorial_network",
                "bond_length_range": None,
                "labile": False,
                "protonation_state": None,
                "source_id": None,
            })

    while len(records) < cn:
        donor_element, atom_name = _DONOR_BY_RESIDUE["HIS"]
        records.append({
            "element": donor_element,
            "residue_name": "HIS",
            "atom_name": atom_name,
            "bond_length": _typical_bond_length(metal_type, donor_element),
            "role": "equatorial_network",
            "bond_length_range": None,
            "labile": False,
            "protonation_state": None,
            "source_id": None,
        })
    for index, record in enumerate(records[:cn]):
        residue = str(record.get("residue_name") or "").upper()
        if residue in {"HOH", "OH"}:
            record["role"] = "axial_labile"
            record["labile"] = True
    return records[:cn]


def _dict_get(obj, key: str, default=None):
    return obj.get(key, default) if isinstance(obj, dict) else default


def _typical_bond_length(metal_type: str, donor_element: str) -> float:
    return _TYPICAL_BOND_LENGTHS.get((metal_type.upper(), donor_element.upper()), 2.05)


def _normalize_geometry(geometry: Optional[str]) -> str:
    key = str(geometry or "").strip().lower().replace(" ", "_").replace("-", "_")
    aliases = {
        "squareplanar": "square_planar",
        "trigonalbipyramidal": "trigonal_bipyramidal",
        "trigonal_pyramidal": "tetrahedral",
        "unknown": "",
    }
    return aliases.get(key, key) or "tetrahedral"


def _coordination_vectors(geometry: str, cn: int) -> np.ndarray:
    geometry = _normalize_geometry(geometry)
    tables = {
        "linear": [[1, 0, 0], [-1, 0, 0]],
        "trigonal_planar": [[1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]],
        "square_planar": [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
        "tetrahedral": [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
        "square_pyramidal": [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        "trigonal_bipyramidal": [[1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0], [0, 0, 1], [0, 0, -1]],
        "octahedral": [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
    }
    raw = np.array(tables.get(geometry, []), dtype=float)
    if len(raw) < cn:
        raw = _fibonacci_vectors(cn)
    else:
        raw = raw[:cn]
    norms = np.linalg.norm(raw, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    return raw / norms


def _fibonacci_vectors(n: int) -> np.ndarray:
    vectors = []
    for i in range(max(n, 1)):
        theta = np.arccos(1 - 2 * (i + 0.5) / max(n, 1))
        phi = np.pi * (1 + 5 ** 0.5) * i
        vectors.append([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    return np.array(vectors, dtype=float)


def _fallback_ligand_stub(donor: Dict, donor_pos: np.ndarray, direction: np.ndarray, idx: int, site_id: Optional[str]) -> List[Dict]:
    u = direction / (np.linalg.norm(direction) or 1.0)
    p, q = _orthonormal_frame(u)
    res = donor["residue_name"]
    elem = donor["element"].upper()

    if elem == "N" and res == "HIS":
        offsets = [
            u * 1.32 + p * 0.72,
            u * 2.15 + p * 0.20 + q * 0.55,
            u * 2.10 - p * 0.75,
            u * 1.25 - p * 0.85 - q * 0.35,
        ]
        elems = ["C", "N", "C", "C"]
        res_name = "IMD"
    elif elem == "O":
        offsets = [u * 1.30, u * 1.30 + p * 1.22]
        elems = ["C", "O"]
        res_name = res[:3]
    elif elem == "S":
        offsets = [u * 1.80, u * 2.85 + p * 0.55]
        elems = ["C", "C"]
        res_name = res[:3]
    else:
        offsets = [u * 1.45, u * 2.50 + p * 0.45]
        elems = ["C", "C"]
        res_name = res[:3]

    atoms = []
    for j, (stub_elem, offset) in enumerate(zip(elems, offsets), 1):
        atoms.append({
            "element": stub_elem,
            "residue_name": res_name,
            "source_residue_name": res,
            "atom_name": f"{stub_elem}{idx + 1}{j}",
            "coords": (donor_pos + offset).tolist(),
            "site_id": site_id,
        })
    return atoms


def _fallback_ligand_fragment(
    donor: Dict,
    donor_pos: np.ndarray,
    direction: np.ndarray,
    idx: int,
    site_id: Optional[str],
) -> List[Dict]:
    """Generate an explicit-H, charged ligand fragment with RDKit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    residue = str(donor["residue_name"]).upper()
    smiles = _FALLBACK_LIGAND_SMILES.get(residue, _FALLBACK_LIGAND_SMILES["HIS"])
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 9107 + idx
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit embedding failed for {residue}")
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    donor_idx = _select_donor_atom(mol, donor["element"], residue)
    conf = mol.GetConformer()
    coords = np.array(
        [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ],
        dtype=float,
    )
    centered = coords - coords[donor_idx]
    heavy = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetAtomicNum() > 1 and i != donor_idx]
    substituents = [i for i in range(mol.GetNumAtoms()) if i != donor_idx]
    source_indices = heavy or substituents
    source = (
        centered[source_indices].mean(axis=0)
        if source_indices
        else np.array([1.0, 0.0, 0.0])
    )
    target = direction / (np.linalg.norm(direction) or 1.0)
    rotation = _rotation_between(source, target)
    transformed = centered @ rotation.T
    spin = _rotation_about_axis(target, 0.65 * idx)
    transformed = transformed @ spin.T + donor_pos

    adjacency = {i: [] for i in range(mol.GetNumAtoms())}
    orders = {i: {} for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        order = float(bond.GetBondTypeAsDouble())
        adjacency[left].append(right)
        adjacency[right].append(left)
        orders[left][right] = order
        orders[right][left] = order

    fragment_id = f"{site_id}:{residue}:{idx}"
    atoms = []
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        element = atom.GetSymbol().upper()
        atoms.append(
            {
                "element": element,
                "residue_name": residue,
                "source_residue_name": residue,
                "atom_name": donor["atom_name"] if atom_idx == donor_idx else f"{element}{idx + 1}_{atom_idx + 1}",
                "coords": transformed[atom_idx].tolist(),
                "site_id": site_id,
                "is_coord_atom": atom_idx == donor_idx,
                "formal_charge": int(atom.GetFormalCharge()),
                "fragment_id": fragment_id,
                "fragment_atom_index": atom_idx,
                "bonded_fragment_indices": sorted(adjacency[atom_idx]),
                "bond_orders": {str(k): v for k, v in orders[atom_idx].items()},
            }
        )
    return atoms


def _select_donor_atom(mol, element: str, residue: str) -> int:
    candidates = [
        atom for atom in mol.GetAtoms() if atom.GetSymbol().upper() == str(element).upper()
    ]
    if not candidates:
        raise ValueError(f"No {element} donor found for {residue}")
    if residue == "HIS":
        neutral_n = [a for a in candidates if a.GetTotalNumHs() == 0 and a.GetFormalCharge() == 0]
        if neutral_n:
            return neutral_n[0].GetIdx()
    charged = [a for a in candidates if a.GetFormalCharge() < 0]
    return (charged or candidates)[0].GetIdx()


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source = source / (np.linalg.norm(source) or 1.0)
    target = target / (np.linalg.norm(target) or 1.0)
    cross = np.cross(source, target)
    dot = float(np.clip(np.dot(source, target), -1.0, 1.0))
    norm = float(np.linalg.norm(cross))
    if norm < 1e-10:
        if dot > 0:
            return np.eye(3)
        axis = np.cross(source, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(axis) < 1e-8:
            axis = np.cross(source, np.array([0.0, 1.0, 0.0]))
        return _rotation_about_axis(axis, math.pi)
    axis = cross / norm
    return _rotation_about_axis(axis, math.acos(dot))


def _rotation_about_axis(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / (np.linalg.norm(axis) or 1.0)
    x, y, z = axis
    c, s = math.cos(angle), math.sin(angle)
    one = 1.0 - c
    return np.array(
        [
            [c + x * x * one, x * y * one - z * s, x * z * one + y * s],
            [y * x * one + z * s, c + y * y * one, y * z * one - x * s],
            [z * x * one - y * s, z * y * one + x * s, c + z * z * one],
        ],
        dtype=float,
    )


def passivate_graphene_edges(atoms: List[Dict], _depth: int = 0) -> List[Dict]:
    """Cap finite graphene carbon edges and newly cut holes with hydrogen."""
    atoms = [dict(atom) for atom in atoms if not atom.get("support_parent")]
    support = [
        (idx, atom)
        for idx, atom in enumerate(atoms)
        if _is_network_support(atom)
        and str(atom.get("element", "")).upper() != "H"
    ]
    if not support:
        return atoms
    positions = {idx: np.array(atom["coords"], dtype=float) for idx, atom in support}
    output = [dict(atom) for atom in atoms]
    h_count = 0
    for idx, atom in support:
        if str(atom.get("element", "")).upper() != "C":
            continue
        neighbors = [
            other_idx
            for other_idx, _ in support
            if other_idx != idx and 1.15 <= np.linalg.norm(positions[idx] - positions[other_idx]) <= 1.70
        ]
        missing = max(0, 3 - len(neighbors))
        if missing == 0:
            continue
        directions = _edge_hydrogen_directions(positions[idx], [positions[n] for n in neighbors], missing)
        for direction in directions:
            direction = _best_cap_direction(
                positions[idx],
                direction,
                output,
                atom.get("atom_name"),
            )
            h_count += 1
            output.append(
                {
                    "element": "H",
                    "residue_name": "GRA",
                    "atom_name": f"HG{h_count}",
                    "coords": (positions[idx] + 1.09 * direction).tolist(),
                    "formal_charge": 0,
                    "support_parent": atom.get("atom_name"),
                }
            )
    if _depth < 4:
        unsafe_parents = _unsafe_cap_parents(output, min_clearance=1.05)
        if unsafe_parents:
            pruned = [
                atom
                for atom in atoms
                if not (
                    atom.get("residue_name") == "GRA"
                    and atom.get("atom_name") in unsafe_parents
                )
            ]
            return passivate_graphene_edges(pruned, _depth=_depth + 1)
    return output


def _unsafe_cap_parents(atoms: List[Dict], min_clearance: float) -> set:
    unsafe = set()
    for i, atom in enumerate(atoms):
        parent = atom.get("support_parent")
        if not parent:
            continue
        pos = np.array(atom["coords"], dtype=float)
        for j, other in enumerate(atoms):
            if i == j:
                continue
            if other.get("atom_name") == parent and other.get("residue_name") == "GRA":
                continue
            if float(np.linalg.norm(pos - np.array(other["coords"], dtype=float))) < min_clearance:
                unsafe.add(parent)
                break
    return unsafe


def _best_cap_direction(
    parent: np.ndarray,
    preferred: np.ndarray,
    existing_atoms: List[Dict],
    parent_name: Optional[str],
) -> np.ndarray:
    """Rotate a cap direction in-plane to avoid H-H and H-ligand overlaps."""
    preferred = preferred / (np.linalg.norm(preferred) or 1.0)
    candidates = [
        _rotation_about_axis(np.array([0.0, 0.0, 1.0]), angle) @ preferred
        for angle in (0.0, math.pi / 6, -math.pi / 6, math.pi / 3, -math.pi / 3, math.pi / 2, -math.pi / 2, math.pi)
    ]
    best = preferred
    best_score = -float("inf")
    for candidate in candidates:
        trial = parent + 1.09 * candidate
        distances = []
        for atom in existing_atoms:
            if atom.get("atom_name") == parent_name and atom.get("residue_name") == "GRA":
                continue
            distances.append(float(np.linalg.norm(trial - np.array(atom["coords"], dtype=float))))
        clearance = min(distances) if distances else 10.0
        direction_fidelity = float(np.dot(candidate, preferred))
        score = clearance + 0.12 * direction_fidelity
        if score > best_score:
            best_score = score
            best = candidate
    return best / (np.linalg.norm(best) or 1.0)


def _edge_hydrogen_directions(
    center: np.ndarray, neighbor_positions: List[np.ndarray], missing: int
) -> List[np.ndarray]:
    if not neighbor_positions:
        return [np.array([math.cos(2 * math.pi * i / 3), math.sin(2 * math.pi * i / 3), 0.0]) for i in range(missing)]
    units = [(pos - center) / (np.linalg.norm(pos - center) or 1.0) for pos in neighbor_positions]
    base = -np.sum(units, axis=0)
    base[2] = 0.0
    if np.linalg.norm(base) < 1e-8:
        first = units[0]
        base = np.array([-first[1], first[0], 0.0])
    base = base / (np.linalg.norm(base) or 1.0)
    if missing == 1:
        return [base]
    angles = np.linspace(-math.pi / 3, math.pi / 3, missing)
    return [_rotation_about_axis(np.array([0.0, 0.0, 1.0]), float(angle)) @ base for angle in angles]


def _orthonormal_frame(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    u = direction / (np.linalg.norm(direction) or 1.0)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(u, ref))) > 0.85:
        ref = np.array([0.0, 1.0, 0.0])
    p = np.cross(u, ref)
    p = p / (np.linalg.norm(p) or 1.0)
    q = np.cross(u, p)
    q = q / (np.linalg.norm(q) or 1.0)
    return p, q


def _remove_graphene_clashes(graphene: List[Dict], complex_atoms: List[Dict]) -> List[Dict]:
    kept = []
    complex_positions = np.array([a["coords"] for a in complex_atoms], dtype=float)
    complex_elements = [str(a["element"]).upper() for a in complex_atoms]
    for atom in graphene:
        pos = np.array(atom["coords"], dtype=float)
        keep = True
        for cpos, elem in zip(complex_positions, complex_elements):
            dist = float(np.linalg.norm(pos - cpos))
            cutoff = 0.78 * (_COVALENT_RADII.get("C", 0.76) + _COVALENT_RADII.get(elem, 0.8))
            if elem in CATALYTIC_METAL_ELEMENTS:
                cutoff = max(cutoff, 1.35)
            if dist < cutoff:
                keep = False
                break
        if keep:
            kept.append(atom)
    return kept


def build_second_shell_on_graphene(
    graphene_atoms: List[Dict],
    metal_center: np.ndarray,
    second_shell_specs: List[Dict],
    n_sites: int = 4,
) -> List[Dict]:
    """Place explicit-H second-shell fragments without intersecting the site.

    Second-shell groups are optional environmental residues, so a rigid,
    collision-aware placement is safer than atom-wise clash relief that can
    destroy their covalent geometry.
    """
    specs = [s for s in second_shell_specs if not s.get("atom_name", "").endswith("CA")]
    added = []
    occupied = [dict(atom) for atom in graphene_atoms]
    for i, spec in enumerate(specs[:n_sites]):
        residue = str(spec.get("residue_name", "HIS")).upper()
        distance = float(spec.get("distance_to_metal", 5.0))
        fragment = _place_second_shell_fragment(
            residue,
            np.array(metal_center, dtype=float),
            distance,
            i,
            occupied + added,
        )
        added.extend(fragment)
    return added


def _fallback_second_shell_fragment(residue: str, anchor: np.ndarray, idx: int) -> List[Dict]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    smiles = _FALLBACK_LIGAND_SMILES.get(residue, _FALLBACK_LIGAND_SMILES["HIS"])
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 17011 + idx
    if AllChem.EmbedMolecule(mol, params) != 0:
        return []
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=150)
    except Exception as _exc:
        # PR4-1 (M6): UFF can fail on unusual topologies; log and continue with
        # the un-optimized ETKDG conformer (still chemically reasonable, just less
        # geometrically clean). Silent pass hid regressions in RDKit input handling.
        import logging as _log
        _log.getLogger("e2n.carbon_scaffold").debug(
            "UFFOptimizeMolecule skipped for %s residue idx=%d (%s)",
            residue, idx, type(_exc).__name__,
        )
    conf = mol.GetConformer()
    coords = np.array(
        [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ],
        dtype=float,
    )
    coords += anchor - coords.mean(axis=0)
    fragment_id = f"SS:{residue}:{idx}"
    adjacency = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    orders = {atom.GetIdx(): {} for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        order = float(bond.GetBondTypeAsDouble())
        adjacency[left].append(right)
        adjacency[right].append(left)
        orders[left][right] = order
        orders[right][left] = order
    return [
        {
            "element": atom.GetSymbol().upper(),
            "residue_name": f"SS{residue[:3]}",
            "source_residue_name": residue,
            "atom_name": f"{atom.GetSymbol().upper()}SS{idx + 1}_{atom.GetIdx() + 1}",
            "coords": coords[atom.GetIdx()].tolist(),
            "formal_charge": int(atom.GetFormalCharge()),
            "fragment_id": fragment_id,
            "fragment_atom_index": atom.GetIdx(),
            "bonded_fragment_indices": sorted(adjacency[atom.GetIdx()]),
            "bond_orders": {
                str(neighbor): order
                for neighbor, order in orders[atom.GetIdx()].items()
            },
        }
        for atom in mol.GetAtoms()
    ]


def _place_second_shell_fragment(
    residue: str,
    metal_center: np.ndarray,
    target_distance: float,
    idx: int,
    occupied: List[Dict],
) -> List[Dict]:
    base = _fallback_second_shell_fragment(residue, np.zeros(3), idx)
    if not base:
        return []
    base_coords = np.array([atom["coords"] for atom in base], dtype=float)
    base_coords -= base_coords.mean(axis=0)
    occupied_coords = np.array([atom["coords"] for atom in occupied], dtype=float)

    best = None
    best_clearance = -float("inf")
    best_anchor = None
    best_anchor_score = -float("inf")
    for radial_offset in (0.0, 0.6, 1.2, 1.8):
        radius = max(target_distance + radial_offset, 2.5)
        for azimuth_idx in range(16):
            azimuth = 2.0 * math.pi * (azimuth_idx + 0.37 * idx) / 16.0
            for z_fraction in (0.45, -0.45, 0.65, -0.65):
                z = radius * z_fraction
                xy = max(radius * radius - z * z, 0.0) ** 0.5
                center = metal_center + np.array(
                    [xy * math.cos(azimuth), xy * math.sin(azimuth), z]
                )
                for spin_idx in range(8):
                    rotation = _rotation_about_axis(
                        center - metal_center,
                        2.0 * math.pi * spin_idx / 8.0,
                    )
                    coords = base_coords @ rotation.T + center
                    clearance = _minimum_cross_distance(coords, occupied_coords)
                    hetero_distances = [
                        float(np.linalg.norm(coords[atom_index] - metal_center))
                        for atom_index, atom in enumerate(base)
                        if str(atom.get("element", "")).upper() in {"N", "O", "S"}
                    ]
                    if hetero_distances and min(hetero_distances) < 2.75:
                        continue
                    anchor = _second_shell_anchor_candidate(base, coords, occupied, metal_center)
                    if anchor:
                        anchor_score = clearance - 0.55 * abs(anchor["distance"] - anchor["ideal_distance"])
                        if anchor_score > best_anchor_score:
                            best_anchor_score = anchor_score
                            best_anchor = anchor
                            best = coords
                            best_clearance = clearance
                        if clearance >= 1.35 and abs(anchor["distance"] - anchor["ideal_distance"]) <= 0.35:
                            return _with_coordinates_and_anchor(base, coords, anchor)
                    if clearance > best_clearance:
                        best_clearance = clearance
                        best = coords
    if best is None or best_clearance < 0.90 or best_anchor is None:
        return []
    return _with_coordinates_and_anchor(base, best, best_anchor)


def _second_shell_anchor_candidate(
    fragment_atoms: List[Dict],
    fragment_coords: np.ndarray,
    occupied: List[Dict],
    metal_center: np.ndarray,
) -> Optional[Dict]:
    support = [
        atom
        for atom in occupied
        if _is_network_support(atom)
        and str(atom.get("element", "")).upper() != "H"
        and atom.get("atom_name")
    ]
    if not support:
        return None

    support_coords = np.array([atom["coords"] for atom in support], dtype=float)
    best = None
    best_score = float("inf")
    for frag_idx, atom in enumerate(fragment_atoms):
        element = str(atom.get("element", "")).upper()
        if element == "H":
            continue
        # Prefer carbon attachment handles; heteroatoms are often the catalytic
        # acid/base/nucleophile face that should remain oriented toward metal.
        element_penalty = 0.0 if element == "C" else 0.35
        distances = np.linalg.norm(support_coords - fragment_coords[frag_idx], axis=1)
        if not len(distances):
            continue
        support_idx = int(np.argmin(distances))
        distance = float(distances[support_idx])
        support_atom = support[support_idx]
        support_element = str(support_atom.get("element", "")).upper()
        ideal = _COVALENT_RADII.get(element, 0.77) + _COVALENT_RADII.get(support_element, 0.77)
        if not (max(1.20, ideal - 0.25) <= distance <= ideal + 0.35):
            continue
        metal_distance = float(np.linalg.norm(fragment_coords[frag_idx] - metal_center))
        score = abs(distance - ideal) + element_penalty - 0.015 * metal_distance
        if score < best_score:
            best_score = score
            best = {
                "fragment_atom_index": frag_idx,
                "support_atom_name": support_atom.get("atom_name"),
                "support_residue_name": support_atom.get("residue_name"),
                "support_element": support_element,
                "distance": distance,
                "ideal_distance": ideal,
            }
    return best


def _minimum_cross_distance(left: np.ndarray, right: np.ndarray) -> float:
    if not len(left) or not len(right):
        return float("inf")
    minimum = float("inf")
    for coord in left:
        minimum = min(minimum, float(np.min(np.linalg.norm(right - coord, axis=1))))
    return minimum


def _with_coordinates(atoms: List[Dict], coords: np.ndarray) -> List[Dict]:
    output = [dict(atom) for atom in atoms]
    for atom, coord in zip(output, coords):
        atom["coords"] = coord.tolist()
    return output


def _with_coordinates_and_anchor(atoms: List[Dict], coords: np.ndarray, anchor: Dict) -> List[Dict]:
    output = _with_coordinates(atoms, coords)
    try:
        atom = output[int(anchor["fragment_atom_index"])]
    except (KeyError, TypeError, ValueError, IndexError):
        return output
    atom["support_anchor_atom_name"] = anchor.get("support_atom_name")
    atom["support_anchor_residue_name"] = anchor.get("support_residue_name")
    atom["support_anchor_element"] = anchor.get("support_element")
    atom["support_anchor_distance"] = float(anchor.get("distance", 0.0) or 0.0)
    return output


def rdkit_optimize(atoms: List[Dict], metal_type: str = None) -> List[Dict]:
    """Lightweight geometry cleanup before optional ASE/ML relaxation.

    RDKit UFF/MMFF cannot reliably parameterize these metal/graphene adducts
    without a bond graph. Placement routines therefore keep every molecular
    fragment rigid and handle clashes before this point. Returning a copy here
    preserves valid ligand bond lengths until a real electronic relaxation is
    requested.
    """
    return [dict(atom) for atom in atoms]
