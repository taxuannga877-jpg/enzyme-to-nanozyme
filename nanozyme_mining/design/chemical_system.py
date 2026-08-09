"""Chemical annotations shared by structure generation and catalysis screening."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ..utils.constants import CATALYTIC_METAL_ELEMENTS


_COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "P": 1.07,
    "S": 1.05,
    "FE": 1.24,
    "CU": 1.17,
    "ZN": 1.25,
    "MN": 1.39,
    "CO": 1.18,
    "NI": 1.17,
    "MO": 1.54,
    "V": 1.53,
    "CR": 1.39,
    "RU": 1.46,
}


@dataclass(frozen=True)
class ChemicalAnnotation:
    bond_graph: List[Tuple[int, int, str]]
    formal_charge: int
    spin_multiplicities: List[int]
    warnings: List[str]


def annotate_chemical_system(atoms: List[Dict], design_spec) -> ChemicalAnnotation:
    """Build a conservative bond graph and whole-cluster electronic metadata."""
    for atom in atoms:
        atom.setdefault("formal_charge", 0)

    metal_by_site = {
        f"M{idx}": metal for idx, metal in enumerate(getattr(design_spec, "metals", []))
    }
    for atom in atoms:
        element = str(atom.get("element", "")).upper()
        site = atom.get("site_id")
        if element in CATALYTIC_METAL_ELEMENTS and site in metal_by_site:
            atom["formal_charge"] = int(metal_by_site[site].oxidation_state)

    charge = int(sum(int(atom.get("formal_charge", 0)) for atom in atoms))
    bond_graph = infer_bond_graph(atoms)
    multiplicities = infer_spin_multiplicities(getattr(design_spec, "metals", []))
    multiplicities = filter_spin_multiplicities_by_electron_parity(
        atoms, charge, multiplicities
    )
    warnings = _chemistry_warnings(atoms, bond_graph, multiplicities)
    return ChemicalAnnotation(
        bond_graph=bond_graph,
        formal_charge=charge,
        spin_multiplicities=multiplicities,
        warnings=warnings,
    )


def infer_bond_graph(atoms: Sequence[Dict]) -> List[Tuple[int, int, str]]:
    bonds: List[Tuple[int, int, str]] = []
    seen_bonds: set[tuple[int, int]] = set()

    def add_bond(left_idx: int, right_idx: int, kind: str) -> None:
        if left_idx == right_idx:
            return
        lo, hi = sorted((int(left_idx), int(right_idx)))
        if lo < 0 or hi >= len(atoms) or (lo, hi) in seen_bonds:
            return
        seen_bonds.add((lo, hi))
        bonds.append((lo, hi, kind))

    fragment_atoms: Dict[tuple[str, int], int] = {}
    for idx, atom in enumerate(atoms):
        fragment_id = atom.get("fragment_id")
        fragment_idx = atom.get("fragment_atom_index")
        if fragment_id is None or fragment_idx is None:
            continue
        try:
            fragment_atoms[(str(fragment_id), int(fragment_idx))] = idx
        except (TypeError, ValueError):
            continue

    for idx, atom in enumerate(atoms):
        fragment_id = atom.get("fragment_id")
        if fragment_id is None or atom.get("fragment_atom_index") is None:
            continue
        orders = atom.get("bond_orders") or {}
        for neighbor, order in orders.items():
            try:
                neighbor_key = (str(fragment_id), int(neighbor))
                order_value = float(order)
            except (TypeError, ValueError):
                continue
            if neighbor_key not in fragment_atoms:
                continue
            add_bond(idx, fragment_atoms[neighbor_key], _bond_kind_from_order(order_value))

    support_by_name = {
        str(atom.get("atom_name")): idx
        for idx, atom in enumerate(atoms)
        if _is_support_atom(atom) and atom.get("atom_name")
    }
    for idx, atom in enumerate(atoms):
        parent = atom.get("support_anchor_atom_name")
        if not parent:
            continue
        parent_idx = support_by_name.get(str(parent))
        if parent_idx is not None:
            add_bond(parent_idx, idx, "single")

    positions = np.array([atom["coords"] for atom in atoms], dtype=float)
    for i, left in enumerate(atoms):
        left_element = str(left.get("element", "")).upper()
        for j in range(i + 1, len(atoms)):
            if (i, j) in seen_bonds:
                continue
            right = atoms[j]
            right_element = str(right.get("element", "")).upper()
            distance = float(np.linalg.norm(positions[i] - positions[j]))
            if distance < 0.35:
                continue

            left_metal = left_element in CATALYTIC_METAL_ELEMENTS
            right_metal = right_element in CATALYTIC_METAL_ELEMENTS
            if left_metal or right_metal:
                metal, donor = (left, right) if left_metal else (right, left)
                same_site_donor = (
                    donor.get("is_coord_atom")
                    and donor.get("site_id") == metal.get("site_id")
                )
                bridge_donor = donor.get("is_bridge_atom") and distance <= 2.55
                if (same_site_donor and distance <= 3.0) or bridge_donor:
                    bonds.append((i, j, "coordinate"))
                continue

            left_support = _is_support_atom(left)
            right_support = _is_support_atom(right)
            if left_support or right_support:
                if not (left_support and right_support):
                    continue
                if left_element == "H" or right_element == "H":
                    hydrogen, carbon = (left, right) if left_element == "H" else (right, left)
                    if hydrogen.get("support_parent") != carbon.get("atom_name"):
                        continue

            radius = _COVALENT_RADII.get(left_element, 0.77) + _COVALENT_RADII.get(right_element, 0.77)
            if distance <= 1.20 * radius:
                add_bond(i, j, "single")
    return bonds


def _bond_kind_from_order(order: float) -> str:
    if order >= 2.5:
        return "triple"
    if order >= 1.75:
        return "double"
    if 1.35 <= order < 1.75:
        return "aromatic"
    return "single"


def _is_support_atom(atom: Dict) -> bool:
    residue = str(atom.get("residue_name", "")).upper()
    return residue == "GRA" or residue.endswith("DP")


def infer_spin_multiplicities(metals: Sequence) -> List[int]:
    if not metals:
        return [1]
    site_unpaired = [
        _site_unpaired_candidates(str(m.metal_type).upper(), int(m.oxidation_state))
        for m in metals
    ]
    totals = set(site_unpaired[0])
    for candidates in site_unpaired[1:]:
        coupled = set()
        for left in totals:
            for right in candidates:
                coupled.update(range(abs(left - right), left + right + 1, 2))
        totals = coupled
    return sorted({unpaired + 1 for unpaired in totals}) or [1]


def _site_unpaired_candidates(element: str, oxidation: int) -> List[int]:
    table = {
        ("FE", 2): [0, 2, 4],
        ("FE", 3): [1, 3, 5],
        ("CU", 1): [0],
        ("CU", 2): [1],
        ("ZN", 2): [0],
        ("MN", 2): [1, 3, 5],
        ("MN", 3): [0, 2, 4],
        ("MN", 4): [1, 3],
        ("CO", 2): [1, 3],
        ("CO", 3): [0, 2, 4],
        ("NI", 2): [0, 2],
    }
    return table.get((element, oxidation), [0])


def filter_spin_multiplicities_by_electron_parity(
    atoms: Sequence[Dict], charge: int, multiplicities: Sequence[int]
) -> List[int]:
    """Return spin multiplicities compatible with the system electron count.

    Local cluster extraction can add or remove an odd number of capping atoms,
    so metal-only spin guesses are not always compatible with the final capped
    atom list. When every supplied guess has the wrong parity, scan the nearest
    lower and higher physically allowed multiplicities instead of sending an
    impossible electronic state to the calculator.
    """
    try:
        from ase.data import atomic_numbers

        electrons = sum(atomic_numbers[str(atom["element"]).capitalize()] for atom in atoms) - charge
    except Exception:
        return list(multiplicities) or [1]
    candidates = sorted({max(int(m), 1) for m in multiplicities}) or [1]
    valid = [m for m in candidates if (electrons - (m - 1)) % 2 == 0]
    if valid:
        return valid

    adjusted = set()
    for multiplicity in candidates:
        for neighbor in (multiplicity - 1, multiplicity + 1):
            if neighbor >= 1 and (electrons - (neighbor - 1)) % 2 == 0:
                adjusted.add(neighbor)
    if adjusted:
        return sorted(adjusted)
    return [1 if electrons % 2 == 0 else 2]


def _chemistry_warnings(
    atoms: Sequence[Dict],
    bond_graph: Sequence[Tuple[int, int, str]],
    multiplicities: Sequence[int],
) -> List[str]:
    warnings: List[str] = []
    if not any(str(atom.get("element", "")).upper() == "H" for atom in atoms):
        warnings.append("structure_has_no_explicit_hydrogens")
    if len(multiplicities) > 1:
        warnings.append("multiple_spin_states_require_energy_comparison")
    if _minimum_interatomic_distance(atoms) < 0.90:
        warnings.append("hard_atomic_clash")

    degrees = {i: 0 for i in range(len(atoms))}
    for left, right, _ in bond_graph:
        degrees[left] += 1
        degrees[right] += 1
    if any(
        atom.get("residue_name") == "GRA"
        and str(atom.get("element", "")).upper() == "C"
        and degrees[i] != 3
        for i, atom in enumerate(atoms)
    ):
        warnings.append("under_or_over_coordinated_graphene_carbon")
    return warnings


def _minimum_interatomic_distance(atoms: Sequence[Dict]) -> float:
    if len(atoms) < 2:
        return float("inf")
    positions = np.array([atom["coords"] for atom in atoms], dtype=float)
    minimum = float("inf")
    for left in range(len(positions)):
        distances = np.linalg.norm(positions[left + 1 :] - positions[left], axis=1)
        if len(distances):
            minimum = min(minimum, float(np.min(distances)))
    return minimum
