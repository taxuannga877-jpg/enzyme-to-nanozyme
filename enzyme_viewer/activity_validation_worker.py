"""Background activity-validation orchestration."""

import json
import logging
import os
from pathlib import Path
import uuid

from enzyme_viewer.design_serialization import _json_safe
from enzyme_viewer.design_store import _assembly_result_payload
from enzyme_viewer.security import clamp_int, parse_bool
from nanozyme_mining.design import structure_exporter
from nanozyme_mining.design.activity_validation_report import (
    render_activity_validation_figures,
    summarize_activity_validation,
    write_activity_validation_report,
)
from nanozyme_mining.design.catalysis_screening import screen_catalysis
from nanozyme_mining.design.nanozyme_assembler import AssemblyResult
from nanozyme_mining.design.substrate_catalog import get_reaction_task

_CONFIG = None
_JOBS = None
_RUN_CATALYSIS_SCREEN_SUBPROCESS = None
log = logging.getLogger("e2n.activity_validation")


def configure_activity_validation_worker(
    *,
    config,
    jobs,
    run_catalysis_screen_subprocess,
) -> None:
    """Bind worker helpers to app-owned runtime state."""
    global _CONFIG, _JOBS, _RUN_CATALYSIS_SCREEN_SUBPROCESS
    _CONFIG = config
    _JOBS = jobs
    _RUN_CATALYSIS_SCREEN_SUBPROCESS = run_catalysis_screen_subprocess


def _config():
    if _CONFIG is None:
        raise RuntimeError("activity validation worker has not been configured")
    return _CONFIG


def _jobs():
    if _JOBS is None:
        raise RuntimeError("activity validation worker has not been configured")
    return _JOBS


def _run_catalysis_subprocess(result, payload_data):
    if _RUN_CATALYSIS_SCREEN_SUBPROCESS is None:
        raise RuntimeError("activity validation worker has not been configured")
    return _RUN_CATALYSIS_SCREEN_SUBPROCESS(result, payload_data)


def _resolve_validation_activities(result, data=None):
    data = data or {}
    requested = data.get("activities")
    if not requested:
        requested = getattr(result.design_spec, "activities", []) or []
    if not requested:
        requested = [
            getattr(metal, "activity_type", None)
            for metal in getattr(result.design_spec, "metals", [])
            if getattr(metal, "activity_type", None)
        ]
    if not requested:
        requested = [getattr(result.design_spec, "nanozyme_type", "")]

    resolved = []
    for item in requested:
        if not item:
            continue
        for part in str(item).split("+"):
            activity = part.strip()
            if activity and activity not in resolved and get_reaction_task(activity):
                resolved.append(activity)
    return resolved


def _validation_progress_callback(task_id: str, activity: str, start: float, span: float):
    def _callback(event: dict) -> None:
        local_progress = float(event.get("progress", 0.0))
        progress = start + span * (local_progress / 100.0)
        detail = {k: v for k, v in event.items() if k not in {"stage", "message", "progress"}}
        _jobs().event(
            task_id,
            event.get("stage", "running"),
            f"{activity}: {event.get('message', '')}",
            progress=progress,
            **detail,
        )

    return _callback


def _run_activity_validation_worker(task_id: str, result: AssemblyResult, data: dict) -> None:
    output_dir = _config()["ACTIVITY_VALIDATION_OUTPUT_DIR"] / task_id
    activity_payloads = []
    try:
        activities = _resolve_validation_activities(result, data)
        if not activities:
            raise ValueError("no supported activity task found for this design")

        _jobs().update(
            task_id,
            status="running",
            activities=activities,
            output_dir=str(output_dir),
        )
        _jobs().event(
            task_id,
            "preparing",
            f"prepared validation plan for {len(activities)} activity task(s)",
            progress=4,
        )

        max_poses = clamp_int(
            data.get("max_adsorption_poses", 3),
            default=3,
            min_value=1,
            max_value=8,
        )
        run_neb = parse_bool(data.get("run_neb"), default=False)
        run_reaction_scan = parse_bool(
            data.get("run_reaction_scan"),
            default=True,
        )
        reaction_scan_points = clamp_int(
            data.get("reaction_scan_points", 5),
            default=5,
            min_value=5,
            max_value=9,
        )
        backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
        use_subprocess = os.environ.get("E2N_CATALYSIS_SUBPROCESS", "1").lower() not in {
            "0",
            "false",
            "no",
            "off",
        }

        per_activity_span = 82.0 / max(len(activities), 1)
        for idx, activity in enumerate(activities):
            task = get_reaction_task(activity)
            activity_start = 8.0 + idx * per_activity_span
            _jobs().event(
                task_id,
                "activity_start",
                f"{activity}: starting substrate and reaction-coordinate validation",
                progress=activity_start,
            )
            payload_data = {
                **data,
                "nanozyme_type": activity,
                "max_adsorption_poses": max_poses,
                "run_neb": run_neb,
                "run_reaction_scan": run_reaction_scan,
                "reaction_scan_points": reaction_scan_points,
            }
            if use_subprocess and backend in {"mace", "fairchem"}:
                _jobs().event(
                    task_id,
                    "subprocess",
                    f"{activity}: running isolated ML process",
                    progress=activity_start + per_activity_span * 0.15,
                )
                payload = _run_catalysis_subprocess(result, payload_data)
                _jobs().event(
                    task_id,
                    "subprocess",
                    f"{activity}: isolated process finished with {payload.get('status')}",
                    progress=activity_start + per_activity_span * 0.95,
                )
            else:
                payload = screen_catalysis(
                    result,
                    task=task,
                    max_adsorption_poses=max_poses,
                    run_neb=run_neb,
                    run_reaction_scan=run_reaction_scan,
                    reaction_scan_points=reaction_scan_points,
                    progress_callback=_validation_progress_callback(
                        task_id,
                        activity,
                        activity_start,
                        per_activity_span,
                    ),
                )
            activity_payloads.append({"activity": activity, "payload": payload})
            partial = summarize_activity_validation(result, activity_payloads)
            _jobs().update(
                task_id,
                partial_results=partial.get("activity_results", []),
            )

        _jobs().event(
            task_id,
            "rendering",
            "rendering validation figures and report",
            progress=92,
        )
        summary = summarize_activity_validation(result, activity_payloads)
        output_dir.mkdir(parents=True, exist_ok=True)
        artifacts = render_activity_validation_figures(summary, output_dir)
        structure_paths = structure_exporter.export(result, str(output_dir))
        for fmt in ("pdb", "xyz", "sdf"):
            path = structure_paths.get(fmt)
            if path:
                artifacts.append(
                    {
                        "kind": "structure",
                        "label": f"Validated design structure ({fmt.upper()})",
                        fmt: Path(path).name,
                    }
                )
        design_spec_path = output_dir / "design_spec.json"
        design_spec_path.write_text(
            json.dumps(
                result.design_spec.to_dict(),
                indent=2,
                ensure_ascii=False,
                default=_json_safe,
            )
            + "\n",
            encoding="utf-8",
        )
        artifacts.append(
            {
                "kind": "design_spec",
                "label": "Design specification",
                "json": design_spec_path.name,
            }
        )
        assembly_result_path = output_dir / "assembly_result.json"
        assembly_result_path.write_text(
            json.dumps(
                _assembly_result_payload(result),
                indent=2,
                ensure_ascii=False,
                default=_json_safe,
            )
            + "\n",
            encoding="utf-8",
        )
        artifacts.append(
            {
                "kind": "assembly_result",
                "label": "Assembly result metadata",
                "json": assembly_result_path.name,
            }
        )
        report_path = write_activity_validation_report(summary, output_dir)
        artifacts.append(
            {
                "kind": "json",
                "label": "Machine-readable validation report",
                "json": Path(report_path).name,
            }
        )
        _jobs().update(
            task_id,
            status="complete",
            stage="complete",
            progress=100.0,
            result=summary,
            artifacts=artifacts,
            output_dir=str(output_dir),
        )
        _jobs().event(
            task_id,
            "complete",
            "activity validation complete",
            progress=100,
        )
    except Exception as exc:
        request_id = uuid.uuid4().hex[:12]
        log.exception(
            "activity validation failed task_id=%s request_id=%s",
            task_id,
            request_id,
        )
        _jobs().update(
            task_id,
            status="failed",
            error="activity validation failed",
            request_id=request_id,
            output_dir=str(output_dir),
        )
        _jobs().event(
            task_id,
            "failed",
            "activity validation failed",
            progress=100,
            request_id=request_id,
        )


def _validation_snapshot_for_client(snapshot: dict) -> dict:
    if not snapshot:
        return snapshot
    task_id = snapshot.get("task_id")
    enriched = dict(snapshot)
    enriched.pop("output_dir", None)
    artifacts = []
    for artifact in snapshot.get("artifacts") or []:
        item = dict(artifact)
        for key in ("png", "svg", "json", "pdb", "xyz", "sdf"):
            filename = item.get(key)
            if filename:
                item[f"{key}_url"] = (
                    f"/api/design/activity_validation/artifact/{task_id}/{filename}"
                )
        artifacts.append(item)
    enriched["artifacts"] = artifacts
    return enriched
