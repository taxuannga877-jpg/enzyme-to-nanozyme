"""Activity validation API routes.

The worker, cache, and persisted-result state stay owned by app.py; this module
only holds the HTTP route shell and receives the shared services by reference.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

from flask import jsonify, request, send_file

from enzyme_viewer.security import (
    error_response,
    is_valid_motif_id,
    require_json_csrf,
    safe_join,
)


@dataclass(frozen=True)
class ActivityValidationRouteServices:
    reference_figure_specs: Callable[[], Mapping[str, Mapping[str, Any]]]
    get_design_result: Callable[[str], Any]
    get_reaction_task: Callable[[str], Any]
    score_payload: Callable[[Any], dict]
    runtime_context: Callable[[Any], dict]
    structure_diagnostics: Callable[[Any], dict]
    reference_figures: Callable[[], list]
    resolve_activities: Callable[..., list]
    jobs: Any
    executor: Any
    run_worker: Callable[..., None]
    snapshot_for_client: Callable[[dict], dict]
    snapshot_from_disk: Callable[[str], Optional[dict]]


def register_activity_validation_routes(app, services: ActivityValidationRouteServices) -> None:
    """Register activity validation APIs while preserving original endpoint names."""

    def api_activity_validation_reference(key):
        spec = (services.reference_figure_specs() or {}).get(key)
        path = Path(spec.get("path") or "") if spec else None
        if not path or not path.exists():
            return jsonify({"status": "error", "error": "reference figure not found"}), 404
        return send_file(path)

    def api_activity_validation_context(job_id):
        if not is_valid_motif_id(job_id):
            return jsonify({"status": "error", "error": "invalid job_id"}), 400
        result = services.get_design_result(job_id)
        if not result:
            return jsonify({"status": "error", "error": "job not found"}), 404

        from nanozyme_mining.design.dopant_modifier import result_to_pdb_string

        activities = services.resolve_activities(result)
        tasks = []
        for activity in activities:
            task = services.get_reaction_task(activity)
            if task:
                tasks.append(task.to_dict())
        return jsonify(
            {
                "status": "success",
                "job_id": result.job_id,
                "label": result.label,
                "atom_count": len(result.atoms),
                "formal_charge": result.formal_charge,
                "spin_multiplicities": result.spin_multiplicities,
                "chemistry_warnings": result.chemistry_warnings,
                "score": services.score_payload(result.score),
                "runtime": services.runtime_context(result),
                "structure_diagnostics": services.structure_diagnostics(result),
                "reference_figures": services.reference_figures(),
                "activities": activities,
                "tasks": tasks,
                "pdb": result_to_pdb_string(result),
            }
        )

    @require_json_csrf
    def api_activity_validation_start(job_id):
        if not is_valid_motif_id(job_id):
            return jsonify({"status": "error", "error": "invalid job_id"}), 400
        result = services.get_design_result(job_id)
        if not result:
            return jsonify({"status": "error", "error": "job not found"}), 404

        try:
            data = request.get_json(silent=True)
            if data is None:
                data = {}
            if not isinstance(data, dict):
                return jsonify(
                    {"status": "error", "error": "request body must be a JSON object"}
                ), 400
            activities = services.resolve_activities(result, data)
            if not activities:
                return jsonify(
                    {"status": "error", "error": "no supported activity task found"}
                ), 400
            task_id = services.jobs.create(job_id=job_id, activities=activities)
            services.executor.submit(
                services.run_worker,
                task_id,
                result,
                data,
            )
            snapshot = services.jobs.get(task_id)
            return jsonify(
                {"status": "success", "task": services.snapshot_for_client(snapshot)}
            )
        except Exception as e:
            return error_response("failed to load activity validation context", exc=e)

    def api_activity_validation_status(task_id):
        if not is_valid_motif_id(task_id):
            return jsonify({"status": "error", "error": "invalid task_id"}), 400
        snapshot = services.jobs.get(task_id) or services.snapshot_from_disk(task_id)
        if not snapshot:
            return jsonify({"status": "error", "error": "task not found"}), 404
        return jsonify(
            {"status": "success", "task": services.snapshot_for_client(snapshot)}
        )

    def api_activity_validation_artifact(task_id, filename):
        if not is_valid_motif_id(task_id):
            return jsonify({"status": "error", "error": "invalid task_id"}), 400
        snapshot = services.jobs.get(task_id) or services.snapshot_from_disk(task_id)
        if not snapshot:
            return jsonify({"status": "error", "error": "task not found"}), 404
        if "/" in filename or "\\" in filename or not filename:
            return jsonify({"status": "error", "error": "invalid artifact filename"}), 400
        allowed = {".png", ".svg", ".json", ".pdb", ".xyz", ".sdf"}
        if Path(filename).suffix.lower() not in allowed:
            return jsonify({"status": "error", "error": "unsupported artifact format"}), 400
        try:
            task_dir = Path(snapshot.get("output_dir") or "")
            path = safe_join(task_dir, filename)
        except Exception as exc:
            return error_response("invalid artifact path", status=400, exc=exc)
        if not path.exists():
            return jsonify({"status": "error", "error": "artifact not found"}), 404
        return send_file(path)

    app.add_url_rule(
        "/api/design/activity_validation/reference/<key>",
        endpoint="api_activity_validation_reference",
        view_func=api_activity_validation_reference,
    )
    app.add_url_rule(
        "/api/design/activity_validation/context/<job_id>",
        endpoint="api_activity_validation_context",
        view_func=api_activity_validation_context,
    )
    app.add_url_rule(
        "/api/design/activity_validation/start/<job_id>",
        endpoint="api_activity_validation_start",
        view_func=api_activity_validation_start,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/design/activity_validation/status/<task_id>",
        endpoint="api_activity_validation_status",
        view_func=api_activity_validation_status,
    )
    app.add_url_rule(
        "/api/design/activity_validation/artifact/<task_id>/<path:filename>",
        endpoint="api_activity_validation_artifact",
        view_func=api_activity_validation_artifact,
    )
