"""Ligand structure route."""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from flask import jsonify, request

from enzyme_viewer.security import (
    error_response,
    is_valid_ec_number,
    is_valid_pdb_id,
)
from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser


@dataclass(frozen=True)
class LigandRouteServices:
    resolve_pdb_library_file: Callable[[str, str], Any]


def register_ligand_routes(app, services: LigandRouteServices) -> None:
    """Register ligand APIs while preserving original endpoint names."""

    def get_ligand_structure():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        ligand_id = data.get("ligand_id", "")
        if data.get("pdb_path"):
            logging.getLogger("e2n.security").warning(
                "client provided pdb_path field to /api/get_ligand_structure ignored"
            )
        ec_number = str(data.get("ec_number", "") or "").strip()
        pdb_id = str(data.get("pdb_id", "") or "").upper().strip()
        ligand_name = data.get("ligand_name", "")
        chain = data.get("chain", "")
        residue_number = data.get("residue_number", "")

        if not is_valid_pdb_id(pdb_id):
            return jsonify({"error": "invalid pdb_id (expected 4-char alphanumeric)"}), 400
        if not is_valid_ec_number(ec_number):
            return jsonify({"error": "invalid ec_number"}), 400

        try:
            pdb_path = services.resolve_pdb_library_file(pdb_id, ec_number)
        except ValueError:
            return jsonify({"error": "invalid path"}), 400
        except FileNotFoundError:
            return jsonify({"error": f"PDB file not found: {pdb_id}"}), 404

        if not ligand_name:
            return jsonify({"error": "Ligand name is required"}), 400

        try:
            ligand_name_norm = str(ligand_name).strip().upper()
            chain_norm = "" if chain in [None, "null", "None"] else str(chain).strip()
            residue_number_norm = None
            if residue_number not in [None, "", "null", "None"]:
                try:
                    residue_number_norm = int(residue_number)
                except (TypeError, ValueError):
                    residue_number_norm = None

            pdb_parser = ComprehensivePDBParser()
            parsed_data = pdb_parser.parse_pdb_file(Path(pdb_path))
            ligands = parsed_data.get("ligands", [])

            matching_atoms = []
            for ligand in ligands:
                res_name = ligand.get("residue_name", "").upper()
                if res_name == "HOH":
                    continue

                ligand_chain = str(ligand.get("chain", "") or "").strip()
                ligand_resnum = ligand.get("residue_number")
                try:
                    ligand_resnum_int = (
                        int(ligand_resnum) if ligand_resnum is not None else None
                    )
                except (TypeError, ValueError):
                    ligand_resnum_int = None

                exact_ok = (
                    res_name == ligand_name_norm
                    and (not chain_norm or ligand_chain == chain_norm)
                    and (
                        residue_number_norm is None
                        or ligand_resnum_int == residue_number_norm
                    )
                )
                fallback_ok = res_name == ligand_name_norm

                if exact_ok or (not matching_atoms and fallback_ok):
                    matching_atoms.append(ligand)

            if not matching_atoms:
                return jsonify(
                    {
                        "status": "error",
                        "error": f"Ligand {ligand_name} not found in PDB file",
                    }
                ), 404

            pdb_lines = []
            atom_serial = 1

            for atom in matching_atoms:
                atom_name = atom.get("atom_name", "")
                element = atom.get("element", "")
                coords = atom.get("coordinates", [0, 0, 0])

                if len(coords) < 3:
                    continue

                x, y, z = coords[0], coords[1], coords[2]

                atom_name_padded = atom_name[:4].ljust(4)
                res_name_padded = ligand_name_norm[:3].ljust(3)
                element_padded = element[:2].rjust(2) if element else "  "

                resnum_for_line = residue_number_norm
                if resnum_for_line is None:
                    try:
                        resnum_for_line = (
                            int(atom.get("residue_number"))
                            if atom.get("residue_number") is not None
                            else 1
                        )
                    except (TypeError, ValueError):
                        resnum_for_line = 1
                chain_for_line = chain_norm or (atom.get("chain") or "A")

                pdb_line = (
                    f"HETATM{atom_serial:5d} {atom_name_padded:4s} "
                    f"{res_name_padded:3s} {str(chain_for_line)[:1]:1s}"
                    f"{int(resnum_for_line):4d}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                    f"{element_padded:2s}  \n"
                )
                pdb_lines.append(pdb_line)
                atom_serial += 1

            pdb_lines.append("END\n")
            pdb_string = "".join(pdb_lines)

            return jsonify(
                {
                    "status": "success",
                    "ligand_id": ligand_id,
                    "ligand_name": ligand_name,
                    "pdb_string": pdb_string,
                    "atom_count": len(matching_atoms),
                }
            )
        except Exception as e:
            return error_response("failed to load ligand structure", exc=e)

    app.add_url_rule(
        "/api/get_ligand_structure",
        endpoint="get_ligand_structure",
        view_func=get_ligand_structure,
        methods=["POST"],
    )
