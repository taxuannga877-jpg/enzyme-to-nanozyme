"""Motif library listing route."""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable

from flask import jsonify, request

from enzyme_viewer.security import error_response


MOTIF_CATEGORIES = [
    "metal_sites",
    "metal_motifs",
    "catalytic_sites",
    "binding_sites",
    "ligands_cofactors",
    "other",
]


_log = logging.getLogger("e2n.routes.motif_listing")
_PRIVATE_PATH_FIELDS = {
    "file_path",
    "file_paths",
    "pdb_path",
    "pdb_paths",
    "rep_pdb_path",
}


@dataclass(frozen=True)
class MotifListingRouteServices:
    motif_output_dir: Callable[[], Path]
    motif_library_dir: Callable[[], Path]
    get_motif_db: Callable[[], Any]
    get_catalytic_metal_db: Callable[[], Any]
    get_ligand_db: Callable[[], Any]
    classify_motif: Callable[[dict], str]
    path_under_any_root: Callable[[Path, list], bool]


def _initial_categories() -> dict[str, list]:
    return {category: [] for category in MOTIF_CATEGORIES}


def _response_category(category: str) -> str:
    return "metal_motifs" if category == "metal_sites" else category


def _safe_category(category: str) -> str:
    return category if category in MOTIF_CATEGORIES else "other"


def _public_record(record: dict) -> dict:
    return {
        key: value
        for key, value in record.items()
        if key not in _PRIVATE_PATH_FIELDS
    }


def _is_ec_directory_name(name: str) -> bool:
    parts = name.split("_")
    if len(parts) != 4:
        return False
    try:
        [int(part) for part in parts]
    except ValueError:
        return False
    return True


def _motif_info_from_payload(motif_data: dict, motif_file: Path, nanozyme_type: str) -> dict:
    category = _response_category(_safe_category(motif_data.get("_category", "other")))
    raw_category = "metal_sites" if category == "metal_motifs" else category
    return {
        "motif_id": motif_data.get("motif_id", ""),
        "uniprot_id": motif_data.get("source_uniprot_id", ""),
        "source_pdb_id": motif_data.get("source_pdb_id", ""),
        "ec_number": motif_data.get("source_ec_number", ""),
        "nanozyme_type": motif_data.get("nanozyme_type", nanozyme_type),
        "anchor_atoms_count": len(motif_data.get("anchor_atoms", [])),
        "category": category,
        "motif_category": raw_category,
        "source": "motif_library",
    }


def _iter_filesystem_motifs(
    motif_library_dir: Path,
    nanozyme_type: str,
    classify_motif: Callable[[dict], str],
) -> Iterable[tuple[str, dict]]:
    for sub_dir in motif_library_dir.iterdir():
        if not sub_dir.is_dir():
            continue

        dir_name = sub_dir.name
        require_type_match = False
        if _is_ec_directory_name(dir_name):
            require_type_match = True
            motif_files = []
            for category_dir in sub_dir.iterdir():
                if category_dir.is_dir():
                    motif_files.extend(category_dir.glob("*.json"))
        elif dir_name.upper() == nanozyme_type.upper():
            motif_files = list(sub_dir.rglob("*.json"))
        else:
            continue

        for motif_file in motif_files:
            try:
                with open(motif_file, "r", encoding="utf-8") as f:
                    motif_data = json.load(f)
            except Exception as e:
                _log.warning("error loading motif %s: %s", motif_file, e)
                continue

            file_nanozyme_type = motif_data.get("nanozyme_type", "").upper()
            if require_type_match and file_nanozyme_type != nanozyme_type.upper():
                continue
            if len(motif_data.get("anchor_atoms", [])) > 50:
                continue

            category = _safe_category(classify_motif(motif_data))
            motif_data["_category"] = category
            yield _response_category(category), _motif_info_from_payload(
                motif_data,
                motif_file,
                nanozyme_type,
            )


def register_motif_listing_routes(app, services: MotifListingRouteServices) -> None:
    """Register motif listing APIs while preserving original endpoint names."""

    def list_motifs():
        if not request.is_json:
            return jsonify({"error": "Missing JSON in request"}), 400

        data = request.get_json()
        nanozyme_type = data.get("nanozyme_type", "")

        if not nanozyme_type:
            return jsonify({"error": "Missing nanozyme type"}), 400

        try:
            motifs_by_category = _initial_categories()
            total_count = 0

            db = services.get_motif_db()
            if db:
                try:
                    if hasattr(db, "get_by_nanozyme_type"):
                        try:
                            db_motifs = db.get_by_nanozyme_type(
                                nanozyme_type,
                                max_anchor_atoms=50,
                            )
                        except TypeError:
                            db_motifs = [
                                motif
                                for motif in db.get_by_nanozyme_type(nanozyme_type)
                                if motif.get("anchor_atoms_count", 0) <= 50
                            ]
                    else:
                        db_motifs = [
                            motif
                            for motif in db.get_all()
                            if motif.get("nanozyme_type", "").upper()
                            == nanozyme_type.upper()
                            and motif.get("anchor_atoms_count", 0) <= 50
                        ]

                    for db_motif in db_motifs:
                        if db_motif.get("anchor_atoms_count", 0) > 50:
                            continue
                        category = db_motif.get("category", "other")
                        if category not in motifs_by_category:
                            file_path = db_motif.get("file_path", "")
                            candidate_path = Path(file_path) if file_path else None
                            roots = [
                                services.motif_output_dir(),
                                services.motif_library_dir(),
                            ]
                            if (
                                candidate_path
                                and candidate_path.exists()
                                and services.path_under_any_root(candidate_path, roots)
                            ):
                                try:
                                    with open(candidate_path, "r", encoding="utf-8") as f:
                                        category = services.classify_motif(json.load(f))
                                except Exception as e:
                                    _log.warning(
                                        "error reclassifying motif from %s: %s",
                                        file_path,
                                        e,
                                    )
                                    category = "other"
                            else:
                                if file_path:
                                    logging.getLogger("e2n.security").warning(
                                        "ignored motif DB file_path outside allowed roots "
                                        "for motif_id=%s",
                                        db_motif.get("motif_id", ""),
                                    )
                                category = "other"

                        category = _safe_category(category)
                        response_category = _response_category(category)
                        motifs_by_category[response_category].append(
                            {
                                "motif_id": db_motif.get("motif_id", ""),
                                "uniprot_id": db_motif.get("uniprot_id", ""),
                                "source_pdb_id": db_motif.get("source_pdb_id", ""),
                                "ec_number": db_motif.get("ec_number", ""),
                                "nanozyme_type": db_motif.get(
                                    "nanozyme_type",
                                    nanozyme_type,
                                ),
                                "anchor_atoms_count": db_motif.get(
                                    "anchor_atoms_count",
                                    0,
                                ),
                                "category": response_category,
                                "motif_category": category,
                                "source": "motif_library",
                            }
                        )
                        total_count += 1
                except Exception as db_error:
                    _log.warning(
                        "motif database query failed, falling back to filesystem: %s",
                        db_error,
                    )
                    db = None

            if not db:
                motif_library_dir = services.motif_library_dir()
                for category, motif_info in _iter_filesystem_motifs(
                    motif_library_dir,
                    nanozyme_type,
                    services.classify_motif,
                ):
                    motifs_by_category[category].append(motif_info)
                    total_count += 1

            try:
                metal_db = services.get_catalytic_metal_db()
                if metal_db:
                    metal_sites_summary = metal_db.get_metal_sites_summary_by_nanozyme_type(
                        nanozyme_type
                    )
                    motifs_by_category["metal_sites"].extend(metal_sites_summary)
                    total_count += len(metal_sites_summary)
                else:
                    motifs_by_category["metal_sites"] = []
            except Exception as e:
                _log.warning("catalytic metal query failed: %s", e)
                motifs_by_category["metal_sites"] = []

            motifs_by_category["metal_sites"].sort(
                key=lambda item: (
                    item.get("occurrence_count", 0),
                    item.get("evidence_score", 0),
                ),
                reverse=True,
            )

            try:
                ligand_db = services.get_ligand_db()
                if ligand_db:
                    ligands_summary = ligand_db.get_ligands_summary_by_nanozyme_type(
                        nanozyme_type
                    )
                    motifs_by_category["ligands_cofactors"] = ligands_summary
                    total_count += len(ligands_summary)
                else:
                    motifs_by_category["ligands_cofactors"] = []
            except Exception as e:
                _log.warning("ligand/cofactor query failed: %s", e)
                motifs_by_category["ligands_cofactors"] = []

            motifs_by_category["ligands_cofactors"].sort(
                key=lambda item: item.get("occurrence_count", 0),
                reverse=True,
            )

            other_count = len(motifs_by_category.get("other", []))
            total_count -= other_count
            motifs_by_category["other"] = []

            final_motifs_by_category = {}
            for category in MOTIF_CATEGORIES:
                entries = motifs_by_category.get(category, [])
                final_motifs_by_category[category] = (
                    [_public_record(entry) for entry in entries]
                    if isinstance(entries, list)
                    else []
                )

            return jsonify(
                {
                    "status": "success",
                    "nanozyme_type": nanozyme_type,
                    "motifs": final_motifs_by_category,
                    "total_count": total_count,
                    "source": "database" if db else "filesystem",
                }
            )
        except Exception as e:
            return error_response("failed to list motifs", exc=e)

    app.add_url_rule(
        "/api/list_motifs",
        endpoint="list_motifs",
        view_func=list_motifs,
        methods=["POST"],
    )
