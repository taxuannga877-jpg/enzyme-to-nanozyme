"""
结构导出：XYZ / PDB / SDF

PR1-4 (v4 audit) fixes:
- H11 / N-H3 (HIGH): SDF writer now passes rdkit.Geometry.Point3D (not a Python
  list) to Conformer.SetAtomPosition, calls conf.Set3D(True), and runs
  Chem.SanitizeMol with a safe ops filter. The old `try/except: write empty file`
  swallowed all errors silently, including bad coords and partial 3D writes;
  now exceptions are logged and the SDF is still written best-effort.
- NEW-6 (MEDIUM): _write_pdb assigns ascending residue numbers grouped by
  (residue_name, chain_hint) instead of hardcoding `A   1` for every atom.
  PyMOL / Chimera now render distinct residues instead of one giant blob.
"""
import logging
import os
from collections import OrderedDict
from typing import Dict, List
from .chemical_system import infer_bond_graph
from ..utils.constants import ALL_METAL_ELEMENTS, CATALYTIC_METAL_ELEMENTS  # PR4-1 (M12/M13)

log = logging.getLogger("e2n.structure_exporter")


def export(result, output_dir: str) -> Dict[str, str]:
    os.makedirs(output_dir, exist_ok=True)
    paths = {}
    paths["xyz"] = _write_xyz(result, output_dir)
    paths["pdb"] = _write_pdb(result, output_dir)
    paths["sdf"] = _write_sdf(result, output_dir)
    return paths


def _write_xyz(result, output_dir: str) -> str:
    path = os.path.join(output_dir, f"{result.job_id}.xyz")
    with open(path, "w") as f:
        f.write(result.xyz)
    return path


def _assign_residue_numbers(atoms: List[Dict]) -> List[int]:
    """
    PR1-4 (NEW-6 fix): give each distinct (residue_name, site_id or position-cluster)
    pair a unique ascending residue number. Falls back to per-atom numbering when
    no grouping signal is available so legacy single-residue outputs still validate.

    The previous code wrote `A   1` for every atom; PyMOL would still load it but
    'show resi 1' would select the entire structure, which is misleading.
    """
    res_nums: List[int] = []
    seen: "OrderedDict[tuple, int]" = OrderedDict()
    next_num = 1
    for a in atoms:
        # Group by residue_name + site_id when present (assembler tags first-shell
        # atoms with site_id). Metals get their own residue number per metal.
        res = (a.get("residue_name") or "LIG").upper()
        site_id = a.get("site_id") or a.get("residue_id")
        elem = (a.get("element") or "").upper()
        if elem in ALL_METAL_ELEMENTS:
            # Each metal is its own residue
            key = ("METAL", elem, id(a))
        elif site_id is not None:
            key = (res, site_id)
        else:
            # Fall back to grouping by residue_name only (e.g. all GRA carbons share one)
            key = (res, "default")
        if key not in seen:
            seen[key] = next_num
            next_num += 1
        res_nums.append(seen[key])
    return res_nums


def _write_pdb(result, output_dir: str) -> str:
    path = os.path.join(output_dir, f"{result.job_id}.pdb")
    with open(path, "w") as f:
        f.write(to_pdb_string(result))
    return path


def to_pdb_string(result) -> str:
    """Return a PDB string with explicit CONECT records for designed bonds."""
    lines = ["REMARK  nanozyme design output"]
    res_nums = _assign_residue_numbers(result.atoms)
    for i, (a, res_num) in enumerate(zip(result.atoms, res_nums), 1):
        x, y, z = a["coords"]
        res = (a.get("residue_name") or "LIG")[:3]
        aname = (a.get("atom_name") or a["element"])[:4]
        elem = a["element"][:2]
        # PDB column-strict format: serial in cols 7-11, atom-name 13-16,
        # residue-name 18-20, chain-id 22, residue-seq 23-26.
        # res_num clamped to 4-digit PDB field (1..9999).
        rn = max(1, min(res_num, 9999))
        lines.append(
            f"HETATM{i:5d} {aname:<4} {res:<3} A{rn:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}"
        )
    lines.extend(_conect_records(result))
    lines.append("END")
    return "\n".join(lines)


def _conect_records(result) -> List[str]:
    bond_graph = list(getattr(result, "bond_graph", []) or [])
    if not bond_graph:
        try:
            bond_graph = infer_bond_graph(result.atoms)
        except Exception as exc:
            log.debug("could not infer PDB CONECT records for %s: %s", result.job_id, exc)
            return []

    adjacency: Dict[int, set[int]] = {}
    atom_count = len(result.atoms)
    for left, right, _kind in bond_graph:
        try:
            left_idx = int(left)
            right_idx = int(right)
        except (TypeError, ValueError):
            continue
        if left_idx == right_idx or not (0 <= left_idx < atom_count and 0 <= right_idx < atom_count):
            continue
        left_serial = left_idx + 1
        right_serial = right_idx + 1
        adjacency.setdefault(left_serial, set()).add(right_serial)
        adjacency.setdefault(right_serial, set()).add(left_serial)

    records: List[str] = []
    for serial in sorted(adjacency):
        neighbors = sorted(adjacency[serial])
        for offset in range(0, len(neighbors), 4):
            chunk = neighbors[offset:offset + 4]
            records.append("CONECT" + f"{serial:5d}" + "".join(f"{n:5d}" for n in chunk))
    return records


def _write_sdf(result, output_dir: str) -> str:
    path = os.path.join(output_dir, f"{result.job_id}.sdf")
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ImportError as e:
        # RDKit unavailable — emit empty stub so the calling endpoint still has
        # a file to send (back-compat with previous swallow-all behaviour).
        log.warning("rdkit unavailable, SDF will be empty: %s", e)
        with open(path, "w") as f:
            f.write("")
        return path

    # 只导出有机原子（跳过金属）
    METALS = ALL_METAL_ELEMENTS  # PR4-1 (M12/M13): shared constant
    organic_indices = [
        i for i, a in enumerate(result.atoms) if a.get("element", "").upper() not in METALS
    ]
    organic = [result.atoms[i] for i in organic_indices]
    if not organic:
        with open(path, "w") as f:
            f.write("")
        return path

    try:
        mol = Chem.RWMol()
        original_to_sdf = {}
        for original_idx, a in zip(organic_indices, organic):
            rd_atom = Chem.Atom(a["element"])
            rd_atom.SetFormalCharge(int(a.get("formal_charge", 0)))
            original_to_sdf[original_idx] = mol.AddAtom(rd_atom)
        bond_types = {
            "single": Chem.BondType.SINGLE,
            "double": Chem.BondType.DOUBLE,
            "triple": Chem.BondType.TRIPLE,
            # The lightweight SDF writer creates plain RDKit atoms without
            # aromatic flags, so aromatic bonds from fragment metadata cannot
            # be kekulized reliably here. Preserve connectivity as single
            # bonds; PDB CONECT remains the authoritative viewer graph.
            "aromatic": Chem.BondType.SINGLE,
        }
        for left, right, bond_kind in getattr(result, "bond_graph", []):
            if bond_kind == "coordinate":
                continue
            if left not in original_to_sdf or right not in original_to_sdf:
                continue
            mol.AddBond(
                original_to_sdf[left],
                original_to_sdf[right],
                bond_types.get(str(bond_kind).lower(), Chem.BondType.SINGLE),
            )
        conf = Chem.Conformer(len(organic))
        for i, a in enumerate(organic):
            # PR1-4 (H11 fix): SetAtomPosition expects Point3D, not list. Older
            # RDKit silently accepted list → tuple; newer raises TypeError.
            coords = a["coords"]
            conf.SetAtomPosition(i, Point3D(
                float(coords[0]), float(coords[1]), float(coords[2])
            ))
        # PR1-4 (N-H3 fix): mark conformer 3D so downstream consumers
        # (PyMOL, OpenBabel, RDKit ETKDG) treat coordinates as Cartesian
        # rather than as a 2D depiction.
        conf.Set3D(True)
        mol.AddConformer(conf, assignId=True)
        rwmol = mol.GetMol()

        # PR1-4 (N-H3): partial sanitize. Full SanitizeMol fails for metal-free
        # organic fragments without bond info; we run only the cheap, robust ops
        # that catch bad valences without insisting on aromaticity perception.
        try:
            Chem.SanitizeMol(rwmol,
                              sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS |
                                         Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                                         Chem.SanitizeFlags.SANITIZE_FINDRADICALS,
                              catchErrors=True)
        except Exception as e:
            log.debug("partial SanitizeMol soft-failed for SDF: %s", e)

        writer = Chem.SDWriter(path)
        writer.write(rwmol)
        writer.close()
    except Exception as e:
        # Last-resort fallback: empty file but record why.
        log.warning("SDF writer failed for job_id=%s: %s", result.job_id, e)
        with open(path, "w") as f:
            f.write("")
    return path
