"""Catalog and EC lookup routes."""

import json
import logging
from dataclasses import dataclass
from typing import Any, Callable

from flask import jsonify, request

from enzyme_viewer.security import error_response


_log = logging.getLogger("e2n.routes.catalog")


@dataclass(frozen=True)
class CatalogRouteServices:
    pdb_library_dir: Callable[[], Any]
    json_cache_dir: Callable[[], Any]
    motif_library_dir: Callable[[], Any]
    get_json_file_path: Callable[[str], Any]
    get_motif_db: Callable[[], Any]
    get_ec_activity_label: Callable[[str], str]


def register_catalog_routes(app, services: CatalogRouteServices) -> None:
    """Register catalog APIs while preserving original endpoint names."""

    def list_ec_with_labels():
        try:
            ec_list = []
            for json_file in services.pdb_library_dir().glob("*/*_sites.json"):
                ec_number = json_file.stem.replace("_sites", "")
                label = services.get_ec_activity_label(ec_number)
                ec_list.append(
                    {
                        "ec_number": ec_number,
                        "label": label,
                        "display": f"{ec_number} ({label})" if label else ec_number,
                    }
                )
            ec_list.sort(key=lambda x: x["ec_number"])
            return jsonify({"status": "success", "ec_list": ec_list})
        except Exception as e:
            return error_response("failed to list labeled EC entries", exc=e)

    def list_ec():
        try:
            ec_list = []

            for json_file in services.pdb_library_dir().glob("*/*_sites.json"):
                ec_number = json_file.stem.replace("_sites", "")
                ec_list.append(ec_number)

            json_cache_dir = services.json_cache_dir()
            if not ec_list and json_cache_dir.exists():
                for json_file in json_cache_dir.glob("*_sites.json"):
                    ec_number = json_file.stem.replace("_sites", "")
                    if ec_number not in ec_list:
                        ec_list.append(ec_number)

            ec_list.sort()
            return jsonify(
                {
                    "status": "success",
                    "ec_list": ec_list,
                    "total": len(ec_list),
                }
            )
        except Exception as e:
            return error_response("failed to list EC numbers", exc=e)

    def list_nanozyme_types():
        try:
            nanozyme_types_set = set()

            db = services.get_motif_db()
            if db:
                try:
                    if hasattr(db, "get_all_nanozyme_types"):
                        nanozyme_types_set = set(db.get_all_nanozyme_types())
                    else:
                        for motif in db.get_all():
                            nanozyme_type = motif.get("nanozyme_type", "")
                            if nanozyme_type:
                                nanozyme_types_set.add(nanozyme_type)
                    _log.info(
                        "database returned %d nanozyme types",
                        len(nanozyme_types_set),
                    )
                except Exception as db_error:
                    _log.warning(
                        "database query failed, falling back to files: %s",
                        db_error,
                    )
                    db = None

            if not db or len(nanozyme_types_set) == 0:
                motif_library_dir = services.motif_library_dir()

                if not motif_library_dir.exists():
                    return jsonify(
                        {
                            "status": "error",
                            "error": "motif library is not available",
                        }
                    ), 404

                _log.warning("scanning motif filesystem (%s)", motif_library_dir)

                for sub_dir in motif_library_dir.iterdir():
                    if not sub_dir.is_dir():
                        continue

                    dir_name = sub_dir.name
                    parts = dir_name.split("_")
                    is_ec_format = False
                    if len(parts) == 4:
                        try:
                            [int(p) for p in parts]
                            is_ec_format = True
                        except ValueError:
                            pass

                    if is_ec_format:
                        for category_dir in sub_dir.iterdir():
                            if not category_dir.is_dir():
                                continue
                            for motif_file in category_dir.glob("*.json"):
                                try:
                                    with open(motif_file, "r", encoding="utf-8") as f:
                                        motif_data = json.load(f)
                                    nanozyme_type = motif_data.get("nanozyme_type", "")
                                    if nanozyme_type:
                                        nanozyme_types_set.add(nanozyme_type)
                                except Exception as e:
                                    _log.warning(
                                        "failed to read motif file %s: %s",
                                        motif_file,
                                        e,
                                    )
                                    continue
                    else:
                        nanozyme_types_set.add(dir_name)

            try:
                from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE

                for ec_dir in services.pdb_library_dir().iterdir():
                    if not ec_dir.is_dir():
                        continue
                    parts = ec_dir.name.split("_")
                    if len(parts) == 4:
                        try:
                            ec_key = ".".join(parts)
                            nanozyme_type = EC_TO_NANOZYME_TYPE.get(ec_key)
                            if nanozyme_type:
                                nanozyme_types_set.add(
                                    nanozyme_type.value
                                    if hasattr(nanozyme_type, "value")
                                    else str(nanozyme_type)
                                )
                        except Exception as exc:
                            logging.getLogger("e2n.app").debug(
                                "EC parse skipped (%s)", type(exc).__name__
                            )
            except Exception as exc:
                logging.getLogger("e2n.app").debug(
                    "EC scan skipped (%s)", type(exc).__name__
                )

            nanozyme_types = sorted(list(nanozyme_types_set))
            return jsonify(
                {
                    "status": "success",
                    "nanozyme_types": nanozyme_types,
                    "total": len(nanozyme_types),
                }
            )
        except Exception as e:
            return error_response("failed to list nanozyme types", exc=e)

    def query_ec():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        ec_number = data.get("ec_number", "").strip()

        if not ec_number:
            return jsonify({"error": "EC number is required"}), 400

        try:
            json_file = services.get_json_file_path(ec_number)

            if not json_file.exists():
                return jsonify(
                    {
                        "status": "success",
                        "ec_number": ec_number,
                        "pdb_list": [],
                        "message": f"No local cache data found for EC {ec_number}",
                    }
                )

            with open(json_file, "r", encoding="utf-8") as f:
                enzyme_data = json.load(f)

            if not enzyme_data:
                return jsonify(
                    {
                        "status": "success",
                        "ec_number": ec_number,
                        "pdb_list": [],
                        "message": f"EC {ec_number} data is empty",
                    }
                )

            pdb_list = []
            ec_dir_name = ec_number.replace(".", "_")
            pdb_library_ec_dir = services.pdb_library_dir() / ec_dir_name

            for idx, entry in enumerate(enzyme_data):
                uniprot_id = entry.get("uniprot_id", "")
                pdb_id = str(entry.get("pdb_id", "")).upper().strip()
                sequence = entry.get("sequence", "")
                has_pdb = False
                if pdb_id and pdb_library_ec_dir.exists():
                    exp_pdb_path = pdb_library_ec_dir / f"{pdb_id}.pdb"
                    has_pdb = exp_pdb_path.is_file()

                pdb_list.append(
                    {
                        "id": idx + 1,
                        "pdb_id": pdb_id,
                        "uniprot_id": uniprot_id,
                        "ec_number": ec_number,
                        "has_pdb": has_pdb,
                        "sequence_length": len(sequence) if sequence else 0,
                        "has_active_sites": len(entry.get("active_sites", [])) > 0,
                    }
                )

            return jsonify(
                {
                    "status": "success",
                    "ec_number": ec_number,
                    "pdb_list": pdb_list,
                    "total": len(pdb_list),
                    "message": f"Found {len(pdb_list)} structures",
                }
            )
        except Exception as e:
            return error_response("failed to query EC structures", exc=e)

    app.add_url_rule(
        "/api/list_ec_with_labels",
        endpoint="list_ec_with_labels",
        view_func=list_ec_with_labels,
        methods=["GET"],
    )
    app.add_url_rule("/api/list_ec", endpoint="list_ec", view_func=list_ec, methods=["GET"])
    app.add_url_rule(
        "/api/list_nanozyme_types",
        endpoint="list_nanozyme_types",
        view_func=list_nanozyme_types,
        methods=["GET"],
    )
    app.add_url_rule("/api/query_ec", endpoint="query_ec", view_func=query_ec, methods=["POST"])
