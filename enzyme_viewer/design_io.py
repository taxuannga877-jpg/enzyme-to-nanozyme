"""Pure file/atom helpers for persisted nanozyme design results."""

import math
from pathlib import Path

from nanozyme_mining.design.design_spec import DesignSpec
from nanozyme_mining.utils.constants import CATALYTIC_METAL_ELEMENTS


def _parse_pdb_atoms(path: Path) -> list:
    atoms = []
    if not path.exists():
        return atoms
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            atom_name = line[12:16].strip() or "X"
            residue_name = line[17:20].strip() or "LIG"
            element = line[76:78].strip()
            if not element:
                element = "".join(ch for ch in atom_name if ch.isalpha())[:2] or "X"
            atoms.append(
                {
                    "element": element.upper(),
                    "atom_name": atom_name,
                    "residue_name": residue_name,
                    "coords": [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    "formal_charge": 0,
                }
            )
        except Exception:
            continue
    return atoms


def _parse_xyz_atoms(path: Path) -> tuple[list, str]:
    if not path.exists():
        return [], ""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    try:
        count = int(lines[0].strip())
    except Exception:
        return [], text
    atoms = []
    for idx, line in enumerate(lines[2 : 2 + count], 1):
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            element = parts[0].upper()
            atoms.append(
                {
                    "element": element,
                    "atom_name": f"{element}{idx}",
                    "residue_name": "LIG",
                    "coords": [float(parts[1]), float(parts[2]), float(parts[3])],
                    "formal_charge": 0,
                }
            )
        except Exception:
            continue
    return atoms, text


def _atoms_to_xyz(atoms: list) -> str:
    lines = [str(len(atoms)), "nanozyme structure"]
    for atom in atoms:
        x, y, z = atom.get("coords", [0.0, 0.0, 0.0])
        lines.append(
            f"{atom.get('element', 'X'):<2} {float(x):12.4f} {float(y):12.4f} {float(z):12.4f}"
        )
    return "\n".join(lines) + "\n"


def _coord_distance(left: dict, right: dict) -> float:
    lx, ly, lz = left.get("coords", [0.0, 0.0, 0.0])
    rx, ry, rz = right.get("coords", [0.0, 0.0, 0.0])
    return math.sqrt(
        (float(lx) - float(rx)) ** 2
        + (float(ly) - float(ry)) ** 2
        + (float(lz) - float(rz)) ** 2
    )


def _reconstruct_cores_for_loaded_design(atoms: list, spec: DesignSpec) -> list:
    cores = []
    used_metals = set()
    used_donors = set()
    metals = [
        idx
        for idx, atom in enumerate(atoms)
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]

    for metal_idx, metal_spec in enumerate(spec.metals):
        site_id = f"M{metal_idx}"
        desired = str(metal_spec.metal_type).upper()
        matching = [
            idx
            for idx in metals
            if idx not in used_metals and str(atoms[idx].get("element", "")).upper() == desired
        ]
        if not matching:
            matching = [idx for idx in metals if idx not in used_metals]
        if not matching:
            continue
        atom_index = matching[0]
        used_metals.add(atom_index)
        metal_atom = atoms[atom_index]
        metal_atom.update(
            {
                "element": desired,
                "atom_name": desired,
                "residue_name": "MET",
                "site_id": site_id,
                "formal_charge": int(metal_spec.oxidation_state),
            }
        )

        coord_atoms = []
        expected_donors = [
            str(coord.donor_element or "").upper() for coord in metal_spec.coord_atoms
        ] or ["N"] * int(metal_spec.coordination_number)
        for donor in expected_donors[: int(metal_spec.coordination_number)]:
            candidates = []
            for idx, atom in enumerate(atoms):
                if idx == atom_index or idx in used_metals or idx in used_donors:
                    continue
                element = str(atom.get("element", "")).upper()
                if element not in {"N", "O", "S"}:
                    continue
                if donor and element != donor:
                    continue
                candidates.append((idx, _coord_distance(metal_atom, atom)))
            if not candidates:
                for idx, atom in enumerate(atoms):
                    if idx == atom_index or idx in used_metals or idx in used_donors:
                        continue
                    element = str(atom.get("element", "")).upper()
                    if element in {"N", "O", "S"}:
                        candidates.append((idx, _coord_distance(metal_atom, atom)))
            if not candidates:
                continue
            donor_idx, _distance = min(candidates, key=lambda item: item[1])
            used_donors.add(donor_idx)
            donor_atom = atoms[donor_idx]
            donor_atom["site_id"] = site_id
            donor_atom["is_coord_atom"] = True
            coord_atoms.append(donor_atom)

        cores.append(
            {
                "site_id": site_id,
                "metal": metal_atom,
                "metal_type": desired,
                "oxidation_state": int(metal_spec.oxidation_state),
                "geometry": metal_spec.coordination_geometry,
                "activity_type": (
                    metal_spec.activity_type
                    or (spec.activities[metal_idx] if metal_idx < len(spec.activities) else None)
                    or spec.nanozyme_type
                ),
                "coord_atoms": coord_atoms,
            }
        )
    return cores
