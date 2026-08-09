"""
基于 MetalloGen 几何方向向量 + BuildAMol MetalComplexer 构建金属配位核心。
输出：包含金属原子和第一配位层原子的 3D 坐标字典。
"""
import numpy as np

from .design_spec import MetalSpec
from .metallogen_bridge import MetalloGenConfig, load_known_geometry_vectors

# 直接从 MetalloGen globalvars 导入方向向量
_GEOM_VECTORS = load_known_geometry_vectors(MetalloGenConfig().metallogen_root)

# catalytic_metal_index.db 中的 coordination_geometry 名称 → MetalloGen 键名映射
_GEOM_MAP = {
    "octahedral":           "6_octahedral",
    "tetrahedral":          "4_tetrahedral",
    "square_planar":        "4_square_planar",
    "square_pyramidal":     "5_square_pyramidal",
    "trigonal_bipyramidal": "5_trigonal_bipyramidal",
    "linear":               "2_linear",
    "bent":                 "2_bent_135",
    "trigonal_planar":      "3_trigonal_planar",
    "seesaw":               "4_seesaw",
}

# 典型配位键长（Å），来自 coordination_distances_json 统计
_TYPICAL_BOND_LENGTHS = {
    ("FE", "N"): 2.05, ("FE", "O"): 1.95, ("FE", "S"): 2.30,
    ("CU", "N"): 2.00, ("CU", "O"): 1.95, ("CU", "S"): 2.20,
    ("ZN", "N"): 2.10, ("ZN", "O"): 2.00, ("ZN", "S"): 2.30,
    ("MN", "N"): 2.15, ("MN", "O"): 1.90, ("MN", "S"): 2.40,
    ("CO", "N"): 2.00, ("CO", "O"): 1.95,
    ("NI", "N"): 2.05, ("NI", "O"): 2.00,
}


def build_metal_core(metal_spec: MetalSpec, origin: np.ndarray = None) -> dict:
    """
    构建单个金属配位核心的 3D 坐标。

    返回:
        {
          "metal": {"element": "FE", "coords": [x,y,z]},
          "coord_atoms": [{"element": "N", "residue_name": "HIS", "atom_name": "NE2",
                           "coords": [x,y,z], "bond_length": 2.05}, ...],
          "geometry": "octahedral",
          "metal_type": "FE",
        }
    """
    if origin is None:
        origin = np.zeros(3)

    geom_key = _GEOM_MAP.get(metal_spec.coordination_geometry, "6_octahedral")
    if geom_key in _GEOM_VECTORS:
        raw_vectors = np.array(_GEOM_VECTORS[geom_key], dtype=float)
    else:
        raw_vectors = _fallback_vectors(metal_spec.coordination_number)

    # 归一化方向向量
    norms = np.linalg.norm(raw_vectors, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit_vectors = raw_vectors / norms

    # 自动补全配位原子到 coordination_number
    _DEFAULT_DONORS = {
        "FE": ("N", "HIS", "NE2"), "CU": ("N", "HIS", "NE2"),
        "ZN": ("N", "HIS", "NE2"), "MN": ("O", "ASP", "OD1"),
        "CO": ("N", "HIS", "NE2"), "NI": ("N", "HIS", "NE2"),
    }
    donor_elem, donor_res, donor_atom = _DEFAULT_DONORS.get(
        metal_spec.metal_type.upper(), ("N", "HIS", "NE2")
    )
    from .design_spec import CoordAtomSpec
    padded = list(metal_spec.coord_atoms)
    while len(padded) < metal_spec.coordination_number:
        padded.append(CoordAtomSpec(
            donor_element=donor_elem, residue_name=donor_res,
            atom_name=donor_atom, bond_length=0.0
        ))

    coord_atoms = []
    n_coord = min(len(padded), len(unit_vectors))
    for i in range(n_coord):
        atom_spec = padded[i]
        # 优先使用用户指定键长，否则查典型值
        bl = atom_spec.bond_length
        if bl <= 0:
            bl = _TYPICAL_BOND_LENGTHS.get(
                (metal_spec.metal_type.upper(), atom_spec.donor_element.upper()), 2.0
            )
        pos = origin + unit_vectors[i] * bl
        coord_atoms.append({
            "element": atom_spec.donor_element,
            "residue_name": atom_spec.residue_name,
            "atom_name": atom_spec.atom_name,
            "coords": pos.tolist(),
            "bond_length": bl,
        })

    return {
        "metal": {"element": metal_spec.metal_type, "coords": origin.tolist()},
        "coord_atoms": coord_atoms,
        "geometry": metal_spec.coordination_geometry,
        "metal_type": metal_spec.metal_type,
        "oxidation_state": metal_spec.oxidation_state,
        "coordination_number": metal_spec.coordination_number,
    }


def _fallback_vectors(cn: int) -> np.ndarray:
    """当 MetalloGen 不可用时的简单回退：均匀分布在球面上"""
    vectors = []
    for i in range(cn):
        theta = np.arccos(1 - 2 * (i + 0.5) / cn)
        phi = np.pi * (1 + 5 ** 0.5) * i
        vectors.append([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    return np.array(vectors)
