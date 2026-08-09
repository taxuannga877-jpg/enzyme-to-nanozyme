"""Basic motif lookup and extraction routes."""

import json
import logging
from dataclasses import dataclass
from typing import Any, Callable

from flask import jsonify, request

from enzyme_viewer.security import (
    error_response,
    is_valid_ec_number,
    is_valid_motif_id,
    is_valid_pdb_id,
    require_json_csrf,
)


_log = logging.getLogger("e2n.routes.motif_basic")


@dataclass(frozen=True)
class MotifBasicRouteServices:
    resolve_motif_json_file: Callable[[str], Any]
    resolve_pdb_library_file: Callable[[str, str], Any]
    get_json_file_path: Callable[[str], Any]
    motif_extractor: Any
    motif_output_dir: Callable[[], Any]


def register_motif_basic_routes(app, services: MotifBasicRouteServices) -> None:
    """Register basic motif APIs while preserving original endpoint names."""

    def get_motif():
        motif_id = request.args.get("motif_id", "")

        if not is_valid_motif_id(motif_id):
            return jsonify({"error": "invalid motif_id"}), 400

        try:
            motif_file = services.resolve_motif_json_file(motif_id)

            if not motif_file or not motif_file.exists():
                return jsonify(
                    {
                        "status": "error",
                        "error": f"Motif file not found: {motif_id}",
                    }
                ), 404

            with open(motif_file, "r", encoding="utf-8") as f:
                motif_data = json.load(f)

            return jsonify({"status": "success", "motif": motif_data})
        except Exception as e:
            return error_response("failed to load motif", exc=e)

    @require_json_csrf
    def extract_motif():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        ec_number = data.get("ec_number", "")
        uniprot_id = data.get("uniprot_id", "")
        pdb_id = str(data.get("pdb_id", "") or "").upper().strip()
        nanozyme_type = data.get("nanozyme_type", "POD")

        if data.get("pdb_path"):
            logging.getLogger("e2n.security").warning(
                "client provided pdb_path field ignored (use pdb_id + ec_number)"
            )

        if not is_valid_pdb_id(pdb_id):
            return jsonify({"error": "invalid pdb_id (expected 4-char alphanumeric)"}), 400
        if not is_valid_ec_number(ec_number):
            return jsonify({"error": "invalid ec_number"}), 400

        try:
            pdb_path = services.resolve_pdb_library_file(pdb_id, ec_number)
        except ValueError as exc:
            return error_response("invalid pdb_id or ec_number", status=400, exc=exc)
        except FileNotFoundError:
            return jsonify({"error": f"PDB file not found: {pdb_id}"}), 404

        try:
            active_site_indices = []
            json_file = services.get_json_file_path(ec_number)

            if json_file.exists():
                with open(json_file, "r", encoding="utf-8") as f:
                    enzyme_data = json.load(f)

                for entry in enzyme_data:
                    entry_pdb_id = str(entry.get("pdb_id", "") or "").upper().strip()
                    if entry_pdb_id == pdb_id or (
                        uniprot_id and entry.get("uniprot_id") == uniprot_id
                    ):
                        active_sites = entry.get("active_sites", [])
                        for site in active_sites:
                            start = site.get("start", 0)
                            end = site.get("end", start)
                            active_site_indices.extend(range(start, end + 1))
                        break

            _log.info("extracting motif for uniprot_id=%s", uniprot_id)
            _log.debug("active site indices: %s", active_site_indices)

            motif = services.motif_extractor.extract_motif(
                pdb_path=str(pdb_path),
                uniprot_id=uniprot_id,
                ec_number=ec_number,
                nanozyme_type=nanozyme_type,
                active_site_indices=active_site_indices if active_site_indices else None,
            )

            if motif is None:
                return jsonify(
                    {
                        "status": "error",
                        "error": (
                            "Failed to extract catalytic motif, "
                            "no catalytic residues found"
                        ),
                    }
                ), 404

            motif_dict = motif.to_dict()
            motif_file = services.motif_output_dir() / f"{motif.motif_id}.json"
            with open(motif_file, "w", encoding="utf-8") as f:
                json.dump(motif_dict, f, indent=2)

            motif_info = {
                "motif_id": motif.motif_id,
                "uniprot_id": uniprot_id,
                "ec_number": ec_number,
                "nanozyme_type": nanozyme_type,
                "anchor_atoms": [
                    {
                        "atom_name": atom.atom_name,
                        "residue_name": atom.residue_name,
                        "residue_number": atom.residue_number,
                        "chain_id": atom.chain_id,
                        "coordinates": atom.coordinates,
                    }
                    for atom in motif.anchor_atoms
                ],
                "geometry_constraints": [
                    {
                        "type": constraint.constraint_type,
                        "atoms": constraint.atom_indices,
                        "value": f"{constraint.value:.2f}",
                        "unit": constraint.unit,
                    }
                    for constraint in motif.geometry_constraints
                ],
            }

            return jsonify(
                {
                    "status": "success",
                    "motif": motif_info,
                    "motif_file": str(motif_file),
                    "message": (
                        f"Successfully extracted {len(motif.anchor_atoms)} "
                        "catalytic residues"
                    ),
                }
            )
        except Exception as e:
            return error_response("failed to extract motif", exc=e)

    app.add_url_rule(
        "/api/get_motif",
        endpoint="get_motif",
        view_func=get_motif,
        methods=["GET"],
    )
    app.add_url_rule(
        "/api/extract_motif",
        endpoint="extract_motif",
        view_func=extract_motif,
        methods=["POST"],
    )
