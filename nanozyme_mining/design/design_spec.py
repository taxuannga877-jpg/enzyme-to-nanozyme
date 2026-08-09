from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass
class CoordAtomSpec:
    """第一配位层原子规格"""
    donor_element: str   # "N", "O", "S"
    residue_name: str    # "HIS", "ASP", "CYS"
    atom_name: str       # "NE2", "OD1", "SG"
    bond_length: float   # Å
    role: str = "equatorial_network"
    bond_length_range: Optional[Tuple[float, float]] = None
    labile: bool = False
    protonation_state: Optional[str] = None
    source_id: Optional[str] = None

    def __post_init__(self) -> None:
        if self.bond_length_range is not None:
            self.bond_length_range = tuple(float(value) for value in self.bond_length_range)


@dataclass
class MetalSpec:
    """单个金属位点规格"""
    metal_type: str                          # "FE", "CU", "ZN"
    oxidation_state: int                     # +2, +3
    coordination_geometry: str               # "octahedral", "tetrahedral", ...
    coordination_number: int                 # 4, 5, 6
    coord_atoms: List[CoordAtomSpec]
    functional_role: str = "catalytic"       # "catalytic" | "structural"
    source_metal_site_id: Optional[str] = None
    activity_type: Optional[str] = None       # activity assigned to this center
    prototype_id: Optional[str] = None
    geometry_family: Optional[str] = None
    allowed_coordination_numbers: List[int] = field(default_factory=list)
    spin_candidates: List[int] = field(default_factory=list)
    condition_id: Optional[str] = None
    microstate_id: Optional[str] = None
    evidence_policy: str = "require_supported"


@dataclass
class SecondShellSpec:
    """第二配位层规格（催化残基/功能基团）"""
    residue_name: str        # "HIS", "ARG", "TYR"
    atom_name: str           # "NE2", "NH1", "OH"
    role: str                # "acid" | "base" | "nucleophile" | "electrostatic" | "hydrogen_bond"
    target_metal_idx: int = 0  # 对应第几个金属（多金属时）
    distance_to_metal: float = 4.0  # 到金属的目标距离 Å


@dataclass
class DesignSpec:
    """完整设计规格 — 汇总用户三步选择"""
    nanozyme_type: str
    ec_numbers: List[str]
    metals: List[MetalSpec]
    activities: List[str] = field(default_factory=list)
    second_shell: List[SecondShellSpec] = field(default_factory=list)
    multi_metal_mode: str = "independent"  # "independent" | "bridged" | "cooperative"
    bridge_residue: str = "HIS"            # bridged 模式下的桥联残基
    bridge_metal_indices: List[int] = field(default_factory=lambda: [0, 1])
    target_metal_distance: float = 12.0   # independent/cooperative 模式下两金属目标间距 Å
    condition_id: Optional[str] = None
    microstate_id: Optional[str] = None
    evidence_policy: str = "require_supported"
    knowledge_version: Optional[str] = None

    def to_dict(self) -> dict:
        import dataclasses
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "DesignSpec":
        activities = list(d.get("activities") or [])
        nanozyme_type = d.get("nanozyme_type") or " + ".join(activities) or "Unknown"
        if not activities and nanozyme_type:
            activities = [part.strip() for part in str(nanozyme_type).split("+") if part.strip()]
        metals = [MetalSpec(
            metal_type=m["metal_type"],
            oxidation_state=m["oxidation_state"],
            coordination_geometry=m["coordination_geometry"],
            coordination_number=m["coordination_number"],
            coord_atoms=[CoordAtomSpec(**a) for a in m["coord_atoms"]],
            functional_role=m.get("functional_role", "catalytic"),
            source_metal_site_id=m.get("source_metal_site_id"),
            activity_type=m.get("activity_type") or m.get("nanozyme_type") or m.get("source_activity"),
            prototype_id=m.get("prototype_id"),
            geometry_family=m.get("geometry_family"),
            allowed_coordination_numbers=list(m.get("allowed_coordination_numbers") or []),
            spin_candidates=list(m.get("spin_candidates") or []),
            condition_id=m.get("condition_id"),
            microstate_id=m.get("microstate_id"),
            evidence_policy=m.get("evidence_policy", "require_supported"),
        ) for m in d["metals"]]
        second_shell = [SecondShellSpec(**s) for s in d.get("second_shell", [])]
        return cls(
            nanozyme_type=nanozyme_type,
            ec_numbers=d.get("ec_numbers", []),
            metals=metals,
            activities=activities,
            second_shell=second_shell,
            multi_metal_mode=d.get("multi_metal_mode", "independent"),
            bridge_residue=d.get("bridge_residue", "HIS"),
            bridge_metal_indices=d.get("bridge_metal_indices", [0, 1]),
            target_metal_distance=d.get("target_metal_distance", 12.0),
            condition_id=d.get("condition_id"),
            microstate_id=d.get("microstate_id"),
            evidence_policy=d.get("evidence_policy", "require_supported"),
            knowledge_version=d.get("knowledge_version"),
        )
