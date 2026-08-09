"""杂原子掺杂：将碳骨架中的 C 原子替换为 N 或 S。"""
import copy

import numpy as np
from .structure_exporter import to_pdb_string
from ..utils.constants import ALL_METAL_ELEMENTS, CATALYTIC_METAL_ELEMENTS  # PR4-1 (M12/M13)

_METALS = ALL_METAL_ELEMENTS  # PR4-1 (M12/M13): shared constant


def _metal_center(atoms):
    """PR2-3 (NEW-5): find the active metal center for distance-aware doping."""
    for a in atoms:
        if str(a.get("element", "")).upper() in _METALS:
            return np.array(a["coords"], dtype=float)
    return None


def apply_doping(result, n_count: int = 0, s_count: int = 0):
    """返回深拷贝的 AssemblyResult，部分 C 原子替换为 N/S（跳过金属和已有配位杂原子）。

    PR2-3 (NEW-5 fix): replaceable C atoms are now ordered by ascending distance
    to the active metal center, so doping preferentially affects atoms in the
    first/second coordination shell where the electronic effect on catalysis is
    largest. Previously enumerate() order doped distant C atoms first, often
    placing N/S far from the metal where they had little effect on activity.

    Falls back to the original enumerate order when there is no metal center
    (e.g. supports without a metal — currently unused but defensive).
    """
    new = copy.deepcopy(result)
    metal_pos = _metal_center(new.atoms)
    replaceable = [
        i for i, a in enumerate(new.atoms)
        if a["element"].upper() == "C"
    ]
    if metal_pos is not None:
        # Order by ascending distance to metal → near-metal C first.
        replaceable.sort(
            key=lambda i: float(np.linalg.norm(
                np.array(new.atoms[i]["coords"], dtype=float) - metal_pos
            ))
        )
    idx = 0
    for _ in range(n_count):
        if idx >= len(replaceable): break
        new.atoms[replaceable[idx]]["element"] = "N"
        idx += 1
    for _ in range(s_count):
        if idx >= len(replaceable): break
        new.atoms[replaceable[idx]]["element"] = "S"
        idx += 1
    new.xyz = _to_xyz(new.atoms)
    return new


def result_to_pdb_string(result) -> str:
    """内存生成 PDB 字符串，供 3Dmol.js 直接消费。"""
    return to_pdb_string(result)


def _to_xyz(atoms) -> str:
    lines = [str(len(atoms)), "nanozyme structure"]
    for a in atoms:
        x, y, z = a["coords"]
        lines.append(f"{a['element']:<4} {x:10.4f} {y:10.4f} {z:10.4f}")
    return "\n".join(lines)
