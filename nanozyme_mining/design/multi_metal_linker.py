"""
多金属连接：三种模式
- independent: 独立并列，有机骨架桥连接
- bridged: 桥联共金属（共享配位原子）
- cooperative: 异核协同（共轭桥联配体）
"""
import numpy as np
from typing import List, Dict

# 有机连接体模板：(近似长度Å, SMILES, 连接原子名)
_LINKER_TEMPLATES = [
    (3.0,  "c1ccncc1",          "N1"),   # 吡啶
    (5.0,  "c1cnccn1",          "N1"),   # 吡嗪
    (7.0,  "c1ccnc(c1)-c1ccccn1", "N1"), # 联吡啶
    (9.0,  "c1ccc(cc1)-c1ccccc1", "C1"), # 联苯
    (12.0, "c1ccc(cc1)CCc1ccccc1", "C1"),# 联苄
    (15.0, "c1ccc(cc1)CCCc1ccccc1","C1"),# 苯丙基苯
]

# 桥联残基 SMILES（His 的 ND1-NE2 可同时配位两个金属）
_BRIDGE_SMILES = {
    "HIS": "c1cnc[nH]1",   # 咪唑，ND1 和 NE2 分别配位
    "CYS": "SCC",          # 硫醇，S 桥联
    "ASP": "OCC(=O)O",     # 天冬氨酸，两个 O 桥联
}

# 共轭桥联配体（cooperative 模式）
_COOPERATIVE_BRIDGES = {
    "pyrazine":     "c1cnccn1",
    "bipyridine":   "c1ccnc(c1)-c1ccccn1",
    "benzimidazole":"c1ccc2[nH]cnc2c1",
}


def link_metal_cores(cores: List[dict], mode: str, **kwargs) -> dict:
    """
    连接多个金属核心，返回包含所有原子的组合结构字典。

    cores: [build_metal_core() 的返回值, ...]
    mode: "independent" | "bridged" | "cooperative"
    kwargs:
        target_distance (float): independent/cooperative 模式下两金属目标间距
        bridge_residue (str): bridged 模式下桥联残基名
        bridge_type (str): cooperative 模式下桥联配体类型
    """
    if len(cores) == 1:
        return _single_core_to_assembly(cores[0])

    if mode == "bridged":
        return _link_bridged(cores, kwargs.get("bridge_residue", "HIS"))
    elif mode == "cooperative":
        return _link_cooperative(cores, kwargs.get("bridge_type", "pyrazine"),
                                 kwargs.get("target_distance", 12.0))
    else:  # independent（默认）
        return _link_independent(cores, kwargs.get("target_distance", 12.0))


def _single_core_to_assembly(core: dict) -> dict:
    atoms = [core["metal"]] + core["coord_atoms"]
    return {"atoms": atoms, "cores": [core], "linker_atoms": [], "mode": "single"}


def _link_independent(cores: List[dict], target_distance: float) -> dict:
    """
    独立并列：将各金属核心平移到合适位置，用有机连接体桥接。
    简化实现：沿 X 轴排列，间距 = target_distance。
    """
    all_atoms = []
    offset = np.zeros(3)

    for i, core in enumerate(cores):
        shift = np.array([i * target_distance, 0.0, 0.0])
        shifted_metal = {**core["metal"], "coords": (np.array(core["metal"]["coords"]) + shift).tolist()}
        shifted_coords = [{**a, "coords": (np.array(a["coords"]) + shift).tolist()} for a in core["coord_atoms"]]
        all_atoms.append(shifted_metal)
        all_atoms.extend(shifted_coords)

    # 在两金属之间插入连接体原子（简化：沿 X 轴均匀分布）
    linker_atoms = _place_linker_atoms(cores, target_distance)
    all_atoms.extend(linker_atoms)

    return {"atoms": all_atoms, "cores": cores, "linker_atoms": linker_atoms, "mode": "independent"}


def _link_bridged(cores: List[dict], bridge_residue: str) -> dict:
    """
    桥联共金属：两金属共享一个配位原子。
    Metal1 放在原点，Metal2 放在 (3.5, 0, 0)，桥联原子在中间。
    """
    bridge_dist = 3.5  # 桥联模式下两金属间距
    all_atoms = []

    # Metal1 在原点
    core1 = cores[0]
    all_atoms.append(core1["metal"])
    all_atoms.extend(core1["coord_atoms"])

    # Metal2 沿 X 轴偏移
    shift = np.array([bridge_dist, 0.0, 0.0])
    core2 = cores[1]
    shifted_metal2 = {**core2["metal"], "coords": (np.array(core2["metal"]["coords"]) + shift).tolist()}
    shifted_coords2 = [{**a, "coords": (np.array(a["coords"]) + shift).tolist()} for a in core2["coord_atoms"]]
    all_atoms.append(shifted_metal2)
    all_atoms.extend(shifted_coords2)

    # 桥联原子放在两金属中间
    bridge_pos = shift / 2
    bridge_atom = {
        "element": "N" if bridge_residue == "HIS" else "S",
        "residue_name": bridge_residue,
        "atom_name": "ND1" if bridge_residue == "HIS" else "SG",
        "coords": bridge_pos.tolist(),
        "is_bridge": True,
    }
    all_atoms.append(bridge_atom)

    return {"atoms": all_atoms, "cores": cores, "linker_atoms": [bridge_atom], "mode": "bridged",
            "bridge_residue": bridge_residue, "metal_distance": bridge_dist}


def _link_cooperative(cores: List[dict], bridge_type: str, target_distance: float) -> dict:
    """
    异核协同：共轭桥联配体连接两金属，保证电子传导。
    """
    all_atoms = []
    shift = np.array([target_distance, 0.0, 0.0])

    core1 = cores[0]
    all_atoms.append(core1["metal"])
    all_atoms.extend(core1["coord_atoms"])

    core2 = cores[1]
    shifted_metal2 = {**core2["metal"], "coords": (np.array(core2["metal"]["coords"]) + shift).tolist()}
    shifted_coords2 = [{**a, "coords": (np.array(a["coords"]) + shift).tolist()} for a in core2["coord_atoms"]]
    all_atoms.append(shifted_metal2)
    all_atoms.extend(shifted_coords2)

    # 共轭桥联配体原子（沿 X 轴均匀分布）
    linker_atoms = _place_linker_atoms(cores, target_distance, conjugated=True)
    all_atoms.extend(linker_atoms)

    return {"atoms": all_atoms, "cores": cores, "linker_atoms": linker_atoms,
            "mode": "cooperative", "bridge_type": bridge_type, "metal_distance": target_distance}


def _place_linker_atoms(cores: List[dict], target_distance: float, conjugated: bool = False) -> List[dict]:
    """在两金属之间沿 X 轴放置连接体原子（简化几何）"""
    # 选择合适的连接体模板
    template = min(_LINKER_TEMPLATES, key=lambda t: abs(t[0] - target_distance))
    smiles = template[1]

    # 简化：用 RDKit 生成连接体原子坐标
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        # 平移到两金属中间
        center = np.array([target_distance / 2, 0.0, 0.0])
        linker_atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() == 1:
                continue  # 跳过氢
            pos = np.array(conf.GetAtomPosition(i))
            pos = pos - pos.mean(axis=0) + center if pos.ndim == 1 else pos
            linker_atoms.append({
                "element": atom.GetSymbol(),
                "residue_name": "LNK",
                "atom_name": f"{atom.GetSymbol()}{i}",
                "coords": (pos + center).tolist(),
                "is_linker": True,
            })
        return linker_atoms
    except Exception:
        return []


def validate_metal_distances(assembly: dict) -> dict:
    """检查多金属间距是否在合理范围内"""
    mode = assembly.get("mode", "independent")
    cores = assembly.get("cores", [])
    if len(cores) < 2:
        return {"valid": True, "errors": []}

    errors = []
    ranges = {
        "independent":  (8.0, 20.0),
        "bridged":      (3.0, 5.0),
        "cooperative":  (10.0, 15.0),
    }
    lo, hi = ranges.get(mode, (0, 999))

    for i in range(len(cores)):
        for j in range(i + 1, len(cores)):
            p1 = np.array(cores[i]["metal"]["coords"])
            p2 = np.array(cores[j]["metal"]["coords"])
            dist = np.linalg.norm(p1 - p2)
            if not (lo <= dist <= hi):
                errors.append(f"Metal {i+1}-{j+1} distance {dist:.1f}Å out of range [{lo},{hi}] for {mode} mode")

    return {"valid": len(errors) == 0, "errors": errors}
