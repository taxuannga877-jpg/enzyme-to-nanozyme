"""
Lightweight bimetallic topology layer for dual-activity nanozyme design.

The shape mirrors the useful architecture of metal-complex builders: metal
centers are vertices, center-center relationships are edges, and geometry
placement is kept separate from relaxation/scoring.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Sequence

import numpy as np

from .design_spec import DesignSpec, MetalSpec
from .physchem_knowledge import evaluate_constructibility


@dataclass(frozen=True)
class MetalCenterNode:
    index: int
    site_id: str
    activity: str
    metal: MetalSpec
    position: tuple[float, float, float]


@dataclass(frozen=True)
class MetalBridgeEdge:
    relation: str
    ideal_distance: float
    distance_range: tuple[float, float]
    bridge_kind: str
    doping: str
    label: str


@dataclass(frozen=True)
class BimetallicTopology:
    nodes: tuple[MetalCenterNode, MetalCenterNode]
    edge: MetalBridgeEdge


def propose_bimetallic_topologies(spec: DesignSpec) -> List[BimetallicTopology]:
    """Return chemically motivated dual-center candidates for the first two metals."""
    if len(spec.metals) < 2:
        return []
    decision = evaluate_constructibility(spec)
    if not decision.constructible:
        return []

    activities = _activities_for_spec(spec)
    metals = spec.metals[:2]
    centers = [
        _node_for_metal(i, metal, activities)
        for i, metal in enumerate(metals)
    ]

    edges = (
        MetalBridgeEdge(
            relation="bridged",
            ideal_distance=5.8,
            distance_range=(4.8, 7.2),
            bridge_kind="diimine_bridge",
            doping="NS",
            label="Bridged dual-metal",
        ),
        MetalBridgeEdge(
            relation="independent_adjacent",
            ideal_distance=7.2,
            distance_range=(6.0, 8.8),
            bridge_kind="none",
            doping="N",
            label="Independent adjacent",
        ),
        MetalBridgeEdge(
            relation="independent_separated",
            ideal_distance=10.5,
            distance_range=(8.0, 14.0),
            bridge_kind="none",
            doping="N",
            label="Independent non-adjacent",
        ),
    )

    return [
        BimetallicTopology(nodes=_position_nodes(centers, edge.ideal_distance), edge=edge)
        for edge in edges
        if edge.relation in set(decision.allowed_modes)
    ]


def bridge_atoms_for_topology(
    topology: BimetallicTopology,
    metal_centers: Sequence[np.ndarray],
) -> List[Dict]:
    """Place bridge atoms only for explicitly bridged bimetallic topologies."""
    if len(metal_centers) < 2:
        return []

    m0 = np.array(metal_centers[0], dtype=float)
    m1 = np.array(metal_centers[1], dtype=float)
    axis = _unit(m1 - m0)
    normal = _perpendicular(axis)
    relation = topology.edge.relation

    if relation == "bridged":
        # The four atoms below form a diimine-like bridge embedded into the
        # carbon support by carbon_scaffold.embed_bridge_in_carbon_network():
        #
        #     N1───C1═══C2───N2
        #     │                 │
        #     M0                M1
        #
        # N1/N2 are donor atoms on the pore edge; C1/C2 are existing support
        # atoms tagged as a conjugated bridge backbone.
        return [
            _bridge_atom("N", "BRG", "N1", m0 + axis * 2.0 + normal * 0.65, topology),
            _bridge_atom("C", "BRG", "C1", m0 + axis * 3.0 + normal * 0.65, topology),
            _bridge_atom("C", "BRG", "C2", m1 - axis * 3.0 - normal * 0.65, topology),
            _bridge_atom("N", "BRG", "N2", m1 - axis * 2.0 - normal * 0.65, topology),
        ]

    return []


def score_metal_distance(
    topology: BimetallicTopology,
    metal_centers: Sequence[np.ndarray],
) -> Dict:
    """Score the relaxed metal-metal distance against the topology target."""
    if len(metal_centers) < 2:
        distance = None
        distance_score = 0.0
    else:
        distance = float(np.linalg.norm(np.array(metal_centers[1]) - np.array(metal_centers[0])))
        lo, hi = topology.edge.distance_range
        _SIGMA_BY_RELATION = {
            "bridged": 0.6,
            "independent_adjacent": 0.8,
            "independent_separated": 1.2,
        }
        sigma = _SIGMA_BY_RELATION.get(topology.edge.relation,
                                         max((hi - lo) / 2.5, 0.5))
        distance_score = float(math.exp(-((distance - topology.edge.ideal_distance) / sigma) ** 2))

    return {
        "relation": topology.edge.relation,
        "label": topology.edge.label,
        "bridge_kind": topology.edge.bridge_kind,
        "doping": topology.edge.doping,
        "metal_distance": distance,
        "ideal_distance": topology.edge.ideal_distance,
        "distance_range": list(topology.edge.distance_range),
        "distance_score": distance_score,
        "activities": [node.activity for node in topology.nodes],
        "metals": [node.metal.metal_type for node in topology.nodes],
    }


def _activities_for_spec(spec: DesignSpec) -> List[str]:
    activities = [a for a in getattr(spec, "activities", []) if a]
    if not activities and spec.nanozyme_type:
        activities = [
            part.strip()
            for part in str(spec.nanozyme_type).split("+")
            if part.strip()
        ]
    while len(activities) < len(spec.metals):
        metal = spec.metals[len(activities)]
        activities.append(metal.activity_type or spec.nanozyme_type or f"activity_{len(activities) + 1}")
    return activities


def _node_for_metal(index: int, metal: MetalSpec, activities: Sequence[str]) -> MetalCenterNode:
    activity = metal.activity_type or (activities[index] if index < len(activities) else "")
    return MetalCenterNode(
        index=index,
        site_id=f"M{index}",
        activity=activity,
        metal=metal,
        position=(0.0, 0.0, 0.0),
    )


def _position_nodes(
    nodes: Sequence[MetalCenterNode],
    distance: float,
) -> tuple[MetalCenterNode, MetalCenterNode]:
    left = MetalCenterNode(
        index=nodes[0].index,
        site_id=nodes[0].site_id,
        activity=nodes[0].activity,
        metal=nodes[0].metal,
        position=(-distance / 2.0, 0.0, 0.0),
    )
    right = MetalCenterNode(
        index=nodes[1].index,
        site_id=nodes[1].site_id,
        activity=nodes[1].activity,
        metal=nodes[1].metal,
        position=(distance / 2.0, 0.0, 0.0),
    )
    return (left, right)


def _bridge_atom(
    element: str,
    residue_name: str,
    atom_name: str,
    coords: np.ndarray,
    topology: BimetallicTopology,
    role: str = "bridge",
) -> Dict:
    return {
        "element": element,
        "residue_name": residue_name,
        "atom_name": atom_name,
        "coords": [float(coords[0]), float(coords[1]), float(coords[2])],
        "is_bridge_atom": True,
        "bridge_role": role,
        "bridge_relation": topology.edge.relation,
    }


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-8:
        return np.array([1.0, 0.0, 0.0])
    return vector / norm


def _perpendicular(axis: np.ndarray) -> np.ndarray:
    trial = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(axis, trial))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    perp = np.cross(axis, trial)
    return _unit(perp)
