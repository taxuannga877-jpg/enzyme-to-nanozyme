"""
将第二配位层（催化残基/功能基团）附加到金属核心坐标上。
第二配位层不直接与金属成键，而是放置在配位原子外侧 3-5Å 处。
"""
import numpy as np
from typing import List, Dict
from .design_spec import SecondShellSpec

# 第二配位层原子到金属的典型距离范围（Å）
_SECOND_SHELL_DISTANCE = {
    "acid":          (3.5, 5.0),
    "base":          (3.0, 4.5),
    "nucleophile":   (3.0, 4.0),
    "electrostatic": (4.0, 6.0),
    "hydrogen_bond": (3.0, 4.5),
}


def attach_second_shell(metal_core: dict, second_shell_specs: List[SecondShellSpec]) -> List[Dict]:
    """
    在金属核心周围放置第二配位层原子。

    策略：在金属核心的配位原子延长线方向上，距金属 target_distance 处放置。
    每个第二配位层原子选择一个"空闲"方向（与已有配位原子方向夹角最大的方向）。

    返回：第二配位层原子列表，每项含 coords, residue_name, atom_name, role
    """
    metal_pos = np.array(metal_core["metal"]["coords"])
    occupied_dirs = [
        _unit(np.array(a["coords"]) - metal_pos)
        for a in metal_core["coord_atoms"]
    ]

    result = []
    used_dirs = list(occupied_dirs)

    for spec in second_shell_specs:
        dist_range = _SECOND_SHELL_DISTANCE.get(spec.role, (3.5, 5.0))
        target_dist = spec.distance_to_metal if spec.distance_to_metal > 0 else (dist_range[0] + dist_range[1]) / 2

        # 找与所有已用方向夹角最大的方向
        direction = _find_free_direction(used_dirs)
        pos = metal_pos + direction * target_dist
        used_dirs.append(direction)

        result.append({
            "element": _donor_element(spec.atom_name),
            "residue_name": spec.residue_name,
            "atom_name": spec.atom_name,
            "role": spec.role,
            "coords": pos.tolist(),
            "target_metal_idx": spec.target_metal_idx,
        })

    return result


def _find_free_direction(used_dirs: list) -> np.ndarray:
    """在球面上找与所有已用方向夹角最大的方向（贪心采样）"""
    if not used_dirs:
        return np.array([0.0, 0.0, 1.0])

    best_dir = None
    best_min_angle = -1.0
    rng = np.random.default_rng(42)

    for _ in range(200):
        candidate = _unit(rng.standard_normal(3))
        min_cos = max(abs(np.dot(candidate, d)) for d in used_dirs)
        min_angle = 1.0 - min_cos  # 越大越好
        if min_angle > best_min_angle:
            best_min_angle = min_angle
            best_dir = candidate

    return best_dir if best_dir is not None else np.array([0.0, 1.0, 0.0])


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-8 else np.array([0.0, 0.0, 1.0])


def _donor_element(atom_name: str) -> str:
    if atom_name.startswith("N"):
        return "N"
    if atom_name.startswith("O"):
        return "O"
    if atom_name.startswith("S"):
        return "S"
    return "C"
