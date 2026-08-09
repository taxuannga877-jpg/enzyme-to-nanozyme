"""PDB structure information and rendering routes."""

import glob
import json
import logging
import os
from dataclasses import dataclass
from typing import Any, Callable

import pandas as pd
from flask import jsonify, request

from enzyme_viewer.security import (
    error_response,
    is_valid_ec_number,
    is_valid_pdb_id,
    safe_join,
)
from enzyme_viewer.structure_info import _build_pdb_info_response
from nanozyme_mining.structure.pdb_metal_extractor import PDBMetalExtractor
from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser


@dataclass(frozen=True)
class StructureRouteServices:
    pdb_library_dir: Callable[[], Any]
    get_json_file_path: Callable[[str], Any]
    render_structure_cached: Callable[..., tuple]


def register_structure_routes(app, services: StructureRouteServices) -> None:
    """Register PDB structure APIs while preserving original endpoint names."""

    def get_pdb_full_info():
        if not request.is_json:
            return jsonify({"error": "Missing JSON"}), 400

        data = request.get_json()
        pdb_id = data.get("pdb_id", "")
        ec_number = data.get("ec_number", "")

        if not pdb_id or not ec_number:
            return jsonify({"error": "Missing pdb_id or ec_number"}), 400

        try:
            if not is_valid_pdb_id(pdb_id):
                return jsonify({"error": "invalid pdb_id (expected 4-char alphanumeric)"}), 400
            if not is_valid_ec_number(ec_number):
                return jsonify({"error": "invalid ec_number"}), 400

            ec_dir = ec_number.replace(".", "_")
            try:
                pdb_dir = safe_join(services.pdb_library_dir(), ec_dir)
            except ValueError:
                return jsonify({"error": "invalid path"}), 400

            pdb_file = None
            for candidate in pdb_dir.glob(f"*{glob.escape(pdb_id)}*.pdb"):
                pdb_file = candidate
                break

            if not pdb_file or not pdb_file.exists():
                return jsonify({"error": f"PDB file not found: {pdb_id}"}), 404

            with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
                pdb_lines = f.readlines()

            parser = ComprehensivePDBParser()
            parsed = parser.parse_pdb_lines(pdb_lines, source_path=pdb_file)

            metal_extractor = PDBMetalExtractor()
            metal_sites_extracted = metal_extractor.parse_pdb_lines(pdb_lines)

            result = _build_pdb_info_response(
                parsed,
                pdb_id,
                ec_number,
                pdb_file,
                metal_sites_extracted,
            )

            pdb_id_for_lookup = str(pdb_id or "").upper()
            json_file = services.get_json_file_path(ec_number)
            if json_file.exists():
                with open(json_file, "r", encoding="utf-8") as f:
                    enzyme_data = json.load(f)
                for entry in enzyme_data:
                    uid = entry.get("uniprot_id", "")
                    pid = str(entry.get("pdb_id", "") or "").upper()
                    if pdb_id_for_lookup == pid or pdb_id_for_lookup == uid:
                        raw_sites = entry.get("active_sites", [])
                        if raw_sites:
                            result["active_sites"] = [
                                {
                                    "site_id": site.get("type", "Site"),
                                    "residues": str(site.get("start", "")),
                                    "description": site.get("description", "") or "\u2014",
                                }
                                for site in raw_sites
                            ]
                        break

            return jsonify({"status": "success", "data": result})
        except Exception as e:
            return error_response("failed to load PDB details", exc=e)

    def get_structure():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        if data.get("pdb_path"):
            logging.getLogger("e2n.security").warning(
                "client provided pdb_path field ignored (use pdb_id + ec_number)"
            )
        ec_number = data.get("ec_number", "")
        uniprot_id = data.get("uniprot_id", "")
        pdb_id = str(data.get("pdb_id", "") or "").upper().strip()

        if not is_valid_pdb_id(pdb_id):
            return jsonify({"error": "invalid pdb_id (expected 4-char alphanumeric)"}), 400
        if not is_valid_ec_number(ec_number):
            return jsonify({"error": "invalid ec_number"}), 400

        try:
            ec_dir = ec_number.replace(".", "_")
            pdb_dir = safe_join(services.pdb_library_dir(), ec_dir)
        except ValueError:
            return jsonify({"error": "invalid path"}), 400

        pdb_path = None
        for candidate in pdb_dir.glob(f"*{glob.escape(pdb_id)}*.pdb"):
            pdb_path = str(candidate)
            break
        if not pdb_path or not os.path.exists(pdb_path):
            return jsonify({"error": f"PDB file not found: {pdb_id}"}), 404

        try:
            site_labels = None
            active_sites_data = []

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
                        if active_sites:
                            site_labels = {}
                            for site in active_sites:
                                site_type_str = site.get("type", "").lower()
                                start = site.get("start", 0)
                                end = site.get("end", start)

                                if "active site" in site_type_str:
                                    site_type = 2
                                elif "binding site" in site_type_str:
                                    site_type = 1
                                else:
                                    site_type = 3

                                for res_idx in range(start, end + 1):
                                    site_labels[res_idx] = site_type

                            for site in active_sites:
                                site_type = site.get("type", "Unknown")
                                start = site.get("start", 0)
                                end = site.get("end", start)
                                description = site.get("description", "")

                                if "active site" in site_type.lower():
                                    color = "#00B050"
                                elif "binding site" in site_type.lower():
                                    color = "#FF0000"
                                else:
                                    color = "#FFFF00"

                                if start == end:
                                    active_sites_data.append(
                                        {
                                            "Residue Index": start,
                                            "Residue Name": "",
                                            "Color": color,
                                            "Active Type": site_type,
                                            "Description": description,
                                        }
                                    )
                                else:
                                    active_sites_data.append(
                                        {
                                            "Residue Index": f"{start}-{end}",
                                            "Residue Name": "",
                                            "Color": color,
                                            "Active Type": site_type,
                                            "Description": description,
                                        }
                                    )
                        break

            view_size = (900, 900)
            show_active = site_labels is not None and len(site_labels) > 0
            structure_html, active_data = services.render_structure_cached(
                enzyme_structure_path=pdb_path,
                site_labels=site_labels,
                view_size=view_size,
                show_active=show_active,
            )

            active_data_html = ""
            if active_sites_data:
                active_data_df = pd.DataFrame(active_sites_data)
                active_data_html = active_data_df.to_html(
                    index=False,
                    escape=True,
                    classes="table table-striped",
                )
            elif active_data and len(active_data) > 0:
                active_data_df = pd.DataFrame(
                    active_data,
                    columns=["Residue Index", "Residue Name", "Color", "Active Type"],
                )
                active_data_html = active_data_df.to_html(
                    index=False,
                    escape=True,
                    classes="table table-striped",
                )
            else:
                active_data_html = "<p>No active site data found</p>"

            return jsonify(
                {
                    "status": "success",
                    "structure_html": structure_html,
                    "active_data_html": active_data_html,
                    "ec_number": ec_number,
                    "uniprot_id": uniprot_id,
                    "has_active_sites": len(active_sites_data) > 0,
                }
            )
        except Exception as e:
            return error_response("failed to render structure", exc=e)

    app.add_url_rule(
        "/api/get_pdb_full_info",
        endpoint="get_pdb_full_info",
        view_func=get_pdb_full_info,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/get_structure",
        endpoint="get_structure",
        view_func=get_structure,
        methods=["POST"],
    )
