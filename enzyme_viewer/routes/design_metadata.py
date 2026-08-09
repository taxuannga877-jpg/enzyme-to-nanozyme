"""Read-only nanozyme design metadata routes."""

from flask import jsonify, request

from enzyme_viewer.security import error_response
from nanozyme_mining.design.motif_selector import (
    get_activity_metals,
    get_coord_templates,
    get_second_shell,
)
from nanozyme_mining.design.substrate_catalog import list_reaction_tasks


def api_get_activity_metals():
    nanozyme_type = request.args.get("nanozyme_type", "")
    if not nanozyme_type:
        return jsonify({"status": "error", "error": "nanozyme_type required"}), 400
    try:
        metals = get_activity_metals(nanozyme_type)
        return jsonify({"status": "success", "metals": metals})
    except Exception as e:
        return error_response("failed to load activity metals", exc=e)


def api_get_coord_templates():
    metal_type = request.args.get("metal_type", "")
    nanozyme_type = request.args.get("nanozyme_type", "")
    oxidation_state_raw = request.args.get("oxidation_state", None)
    if not metal_type:
        return jsonify({"status": "error", "error": "metal_type required"}), 400
    try:
        oxidation_state = (
            int(oxidation_state_raw) if oxidation_state_raw not in (None, "") else None
        )
        templates = get_coord_templates(
            metal_type,
            nanozyme_type,
            oxidation_state=oxidation_state,
        )
        return jsonify({"status": "success", "templates": templates})
    except ValueError:
        return jsonify({"status": "error", "error": "oxidation_state must be an integer"}), 400
    except Exception as e:
        return error_response("failed to load coordination templates", exc=e)


def api_get_second_shell():
    nanozyme_type = request.args.get("nanozyme_type", "")
    metal_type = request.args.get("metal_type", "")
    try:
        residues = get_second_shell(nanozyme_type, metal_type)
        return jsonify({"status": "success", "residues": residues})
    except Exception as e:
        return error_response("failed to load second-shell recommendations", exc=e)


def api_design_substrate_tasks():
    nanozyme_type = request.args.get("nanozyme_type", "")
    tasks = list_reaction_tasks(nanozyme_type)
    return jsonify(
        {
            "status": "success",
            "nanozyme_type": nanozyme_type,
            "tasks": [task.to_dict() for task in tasks],
            "total": len(tasks),
        }
    )


def register_design_metadata_routes(app) -> None:
    """Register read-only design metadata APIs while preserving endpoint names."""
    app.add_url_rule(
        "/api/design/get_activity_metals",
        endpoint="api_get_activity_metals",
        view_func=api_get_activity_metals,
    )
    app.add_url_rule(
        "/api/design/get_coord_templates",
        endpoint="api_get_coord_templates",
        view_func=api_get_coord_templates,
    )
    app.add_url_rule(
        "/api/design/get_second_shell",
        endpoint="api_get_second_shell",
        view_func=api_get_second_shell,
    )
    app.add_url_rule(
        "/api/design/substrate_tasks",
        endpoint="api_design_substrate_tasks",
        view_func=api_design_substrate_tasks,
    )
