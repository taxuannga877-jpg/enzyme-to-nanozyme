"""Motif structure detail route."""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from flask import jsonify, request

from enzyme_viewer.security import error_response, is_valid_motif_id


_log = logging.getLogger("e2n.routes.motif_structure")


@dataclass(frozen=True)
class MotifStructureRouteServices:
    motif_output_dir: Callable[[], Path]
    motif_library_dir: Callable[[], Path]
    get_motif_db: Callable[[], Any]
    get_catalytic_metal_db: Callable[[], Any]
    resolve_motif_json_file: Callable[[str], Any]
    path_under_any_root: Callable[[Path, list], bool]


def generate_pdb_from_residues(residue_structures, anchor_atoms):
    """Generate a compact PDB string from motif residue structures."""
    pdb_lines = []
    atom_serial = 1

    standard_aa = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "SEC",
        "PYL",
    }

    for key_str, residue in residue_structures.items():
        parts = key_str.split("_")
        if len(parts) < 2:
            continue
        res_name = parts[0].upper()
        if res_name not in standard_aa:
            continue
        try:
            res_num = int(parts[1])
        except ValueError:
            continue

        chain_id = residue.get("chain_id", "A")
        for atom_info in residue.get("atoms", []):
            atom_name = atom_info.get("atom_name", "")
            element = atom_info.get("element", "")
            coords = atom_info.get("coordinates", [0, 0, 0])
            occupancy = atom_info.get("occupancy", 1.0)
            bfactor = atom_info.get("bfactor", 0.0)

            if len(coords) < 3:
                continue

            x, y, z = coords[0], coords[1], coords[2]
            atom_name_padded = atom_name[:4].ljust(4)
            res_name_padded = res_name[:3].ljust(3)
            element_padded = element[:2].rjust(2) if element else "  "
            pdb_lines.append(
                f"ATOM  {atom_serial:5d} {atom_name_padded:4s} "
                f"{res_name_padded:3s} {chain_id:1s}"
                f"{res_num:4d}   {x:8.3f}{y:8.3f}{z:8.3f}"
                f"{occupancy:6.2f}{bfactor:6.2f}          {element_padded:2s}  \n"
            )
            atom_serial += 1

    processed_ligands = set()
    for key_str, residue in residue_structures.items():
        parts = key_str.split("_")
        if len(parts) < 2:
            continue
        res_name = parts[0].upper()
        if res_name in standard_aa:
            continue
        if res_name in processed_ligands:
            continue
        processed_ligands.add(res_name)

        try:
            res_num = int(parts[1])
        except ValueError:
            continue

        chain_id = residue.get("chain_id", "A")
        for atom_info in residue.get("atoms", []):
            atom_name = atom_info.get("atom_name", "")
            element = atom_info.get("element", "")
            coords = atom_info.get("coordinates", [0, 0, 0])
            occupancy = atom_info.get("occupancy", 1.0)
            bfactor = atom_info.get("bfactor", 0.0)

            if len(coords) < 3:
                continue

            x, y, z = coords[0], coords[1], coords[2]
            atom_name_padded = atom_name[:4].ljust(4)
            res_name_padded = res_name[:3].ljust(3)
            element_padded = element[:2].rjust(2) if element else "  "
            pdb_lines.append(
                f"HETATM{atom_serial:5d} {atom_name_padded:4s} "
                f"{res_name_padded:3s} {chain_id:1s}"
                f"{res_num:4d}   {x:8.3f}{y:8.3f}{z:8.3f}"
                f"{occupancy:6.2f}{bfactor:6.2f}          {element_padded:2s}  \n"
            )
            atom_serial += 1

    pdb_lines.append("END\n")
    return "".join(pdb_lines)


def _is_ec_directory_name(name: str) -> bool:
    parts = name.split("_")
    if len(parts) != 4:
        return False
    try:
        [int(part) for part in parts]
    except ValueError:
        return False
    return True


def _fallback_find_motif_file(motif_id: str, nanozyme_type: str, motif_library_dir: Path):
    if nanozyme_type:
        nanozyme_dir = motif_library_dir / nanozyme_type
        if nanozyme_dir.exists() and nanozyme_dir.is_dir():
            potential_file = nanozyme_dir / f"{motif_id}.json"
            if potential_file.exists():
                return potential_file
            for motif_file_path in nanozyme_dir.rglob(f"{motif_id}.json"):
                return motif_file_path

    for nanozyme_dir in motif_library_dir.iterdir():
        if not nanozyme_dir.is_dir():
            continue
        if nanozyme_type and _is_ec_directory_name(nanozyme_dir.name):
            continue

        potential_file = nanozyme_dir / f"{motif_id}.json"
        if potential_file.exists():
            return potential_file
        for motif_file_path in nanozyme_dir.rglob(f"{motif_id}.json"):
            return motif_file_path
    return None


def _metal_site_payload(motif_id: str, nanozyme_type: str, services: MotifStructureRouteServices):
    metal_type = motif_id.split("_metal_")[0]
    metal_db = services.get_catalytic_metal_db()
    if not metal_db:
        return jsonify(
            {
                "status": "error",
                "error": (
                    "Catalytic metal-site database not found. "
                    "Build catalytic_metal_index.db first."
                ),
            }
        ), 500

    site_row = None
    if nanozyme_type:
        sites = metal_db.get_sites_for_metal(
            nanozyme_type=nanozyme_type,
            metal_type=metal_type,
            limit=1,
        )
        site_row = sites[0] if sites else None
    if site_row is None:
        site_row = metal_db.get_any_site_for_metal(metal_type=metal_type)

    if not site_row:
        return jsonify(
            {
                "status": "error",
                "error": (
                    "No catalytic metal-site record found in catalytic_metal_index "
                    f"for metal_type={metal_type}"
                ),
            }
        ), 404

    metal_name = site_row.get("metal_name") or metal_type
    chain = site_row.get("metal_chain") or ""
    residue_number = site_row.get("metal_residue_id") or 1
    coords = site_row.get("metal_coords") or [0.0, 0.0, 0.0]
    try:
        x, y, z = float(coords[0]), float(coords[1]), float(coords[2])
    except Exception as exc:
        logging.getLogger("e2n.app").debug(
            "coord parse failed (%s); using origin",
            type(exc).__name__,
        )
        x, y, z = 0.0, 0.0, 0.0

    coordinating_residues = site_row.get("coordinating_residues") or []
    coordinating_residue_stats = {}
    for residue in coordinating_residues:
        residue_name = (residue.get("residue_name") or "UNK").upper()
        coordinating_residue_stats[residue_name] = (
            coordinating_residue_stats.get(residue_name, 0) + 1
        )

    element = metal_type[:2].upper()
    atom_name_padded = element.rjust(2).ljust(4)
    res_name_padded = metal_type[:3].ljust(3)
    element_padded = element.rjust(2)
    pdb_string = (
        f"HETATM{1:5d} {atom_name_padded:4s} {res_name_padded:3s} "
        f"{chain[:1] if chain else 'A':1s}"
        f"{int(residue_number):4d}   {x:8.3f}{y:8.3f}{z:8.3f}  "
        f"1.00  0.00          {element_padded:2s}  \n"
        "END\n"
    )

    return jsonify(
        {
            "status": "success",
            "motif": {
                "motif_id": motif_id,
                "metal_type": metal_type,
                "metal_name": metal_name,
                "ligand_name": metal_name,
                "category": "metal_sites",
                "anchor_atoms": [],
                "anchor_atoms_count": 0,
                "geometry_constraints": [],
                "chemistry_tag": f"Metal site: {metal_type}",
                "extraction_method": "catalytic_metal_index",
                "is_metal_site": True,
                "uniprot_id": "",
                "ec_number": site_row.get("ec_number", ""),
                "nanozyme_type": site_row.get("nanozyme_type", nanozyme_type or ""),
                "coordination_number": site_row.get("coordination_number", 0),
                "coordination_geometry": site_row.get(
                    "coordination_geometry",
                    "unknown",
                ),
                "oxidation_state": site_row.get("oxidation_state", None),
                "coordinating_residues": coordinating_residues,
                "coordination_distances": site_row.get("coordination_distances", []),
                "coordination_angles": site_row.get("coordination_angles", []),
                "coordinating_residue_stats": coordinating_residue_stats,
                "evidence_score": site_row.get("evidence_score", 0),
                "evidence_level": site_row.get("evidence_level", "none"),
                "included_reason": site_row.get("included_reason", ""),
                "excluded_reason": site_row.get("excluded_reason", ""),
                "pdb_id": site_row.get("pdb_id", ""),
            },
            "pdb_string": pdb_string,
            "source": "catalytic_metal_index",
        }
    )


def register_motif_structure_routes(app, services: MotifStructureRouteServices) -> None:
    """Register motif structure APIs while preserving original endpoint names."""

    def get_motif_structure():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        motif_id = data.get("motif_id", "")
        nanozyme_type = data.get("nanozyme_type", "")

        if not motif_id:
            return jsonify({"error": "Missing motif_id"}), 400
        if not is_valid_motif_id(motif_id):
            return jsonify({"error": "invalid motif_id"}), 400

        try:
            motif_file = None
            db = services.get_motif_db()
            if db:
                db_motif = db.get_by_id(motif_id)
                if db_motif:
                    db_path = Path(db_motif.get("file_path") or "")
                    roots = [
                        services.motif_output_dir(),
                        services.motif_library_dir(),
                    ]
                    if db_path.exists() and services.path_under_any_root(db_path, roots):
                        motif_file = db_path
                        _log.info("found motif file from database: %s", motif_file)
                    elif db_path:
                        logging.getLogger("e2n.security").warning(
                            "ignored motif DB file_path outside allowed roots for motif_id=%s",
                            motif_id,
                        )

            if not motif_file or not motif_file.exists():
                motif_file = services.resolve_motif_json_file(motif_id)

            if not motif_file or not motif_file.exists():
                motif_file = _fallback_find_motif_file(
                    motif_id,
                    nanozyme_type,
                    services.motif_library_dir(),
                )

            if not motif_file or not motif_file.exists():
                if "_metal_" in motif_id:
                    return _metal_site_payload(motif_id, nanozyme_type, services)
                return jsonify(
                    {
                        "status": "error",
                        "error": f"Motif file not found: {motif_id}",
                    }
                ), 404

            with open(motif_file, "r", encoding="utf-8") as f:
                motif_data = json.load(f)

            formatted_motif = {
                "motif_id": motif_data.get("motif_id", ""),
                "uniprot_id": motif_data.get("source_uniprot_id", ""),
                "source_pdb_id": motif_data.get("source_pdb_id", ""),
                "ec_number": motif_data.get("source_ec_number", ""),
                "nanozyme_type": motif_data.get("nanozyme_type", ""),
                "anchor_atoms": motif_data.get("anchor_atoms", []),
                "geometry_constraints": [
                    {
                        "constraint_type": constraint.get("constraint_type", ""),
                        "atom_indices": constraint.get("atom_indices", []),
                        "value": f"{constraint.get('value', 0):.2f}",
                        "unit": constraint.get("unit", ""),
                    }
                    for constraint in motif_data.get("geometry_constraints", [])
                ],
                "chemistry_tag": motif_data.get("chemistry_tag", ""),
                "reaction_smiles": motif_data.get("reaction_smiles", ""),
                "extraction_method": motif_data.get("extraction_method", ""),
                "confidence_score": motif_data.get("confidence_score", 0.0),
                "residue_structures": motif_data.get("residue_structures", {}),
                "structure_2d_svg": motif_data.get("structure_2d_svg", ""),
            }

            pdb_string = ""
            if formatted_motif.get("residue_structures"):
                try:
                    pdb_string = generate_pdb_from_residues(
                        formatted_motif["residue_structures"],
                        formatted_motif["anchor_atoms"],
                    )
                except Exception as e:
                    _log.warning("failed to generate PDB string: %s", e)

            return jsonify(
                {
                    "status": "success",
                    "motif": formatted_motif,
                    "pdb_string": pdb_string,
                    "source": (
                        "database"
                        if db and db.get_by_id(motif_id)
                        else "filesystem"
                    ),
                }
            )
        except Exception as e:
            return error_response("failed to load motif structure", exc=e)

    app.add_url_rule(
        "/api/get_motif_structure",
        endpoint="get_motif_structure",
        view_func=get_motif_structure,
        methods=["POST"],
    )
