"""Nanozyme design job API routes."""

import os
import shutil
import tempfile
from dataclasses import dataclass
from typing import Any, Callable

from flask import after_this_request, jsonify, request, send_file

from enzyme_viewer.security import (
    clamp_float,
    clamp_int,
    error_response,
    is_valid_download_format,
    is_valid_motif_id,
    parse_bool,
    require_json_csrf,
    safe_join,
)
from nanozyme_mining.design import structure_exporter
from nanozyme_mining.design.catalysis_screening import screen_catalysis
from nanozyme_mining.design.design_spec import DesignSpec
from nanozyme_mining.design.nanozyme_assembler import NanozymeAssembler


@dataclass(frozen=True)
class DesignJobRouteServices:
    design_results: Any
    persist_design_result: Callable[[Any], Any]
    get_design_result: Callable[[str], Any]
    design_result_dir: Callable[[str], Any]
    score_payload: Callable[[Any], dict]
    get_reaction_task: Callable[[str], Any]
    run_assemble_subprocess: Callable[..., dict]
    run_catalysis_screen_subprocess: Callable[[Any, dict], dict]
    assembly_from_dict: Callable[[dict], Any]


def register_design_job_routes(app, services: DesignJobRouteServices) -> None:
    """Register design job APIs while preserving original endpoint names."""

    @require_json_csrf
    def api_design_assemble():
        try:
            data = request.get_json(silent=True)
            if not isinstance(data, dict):
                return jsonify(
                    {"status": "error", "error": "request body must be a JSON object"}
                ), 400
            spec = DesignSpec.from_dict(data)
            if not spec.metals:
                return jsonify(
                    {"status": "error", "error": "at least one metal site is required"}
                ), 400
            if len(spec.metals) > 2:
                return jsonify(
                    {"status": "error", "error": "at most two metal sites are supported"}
                ), 400
            backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
            use_subprocess = os.environ.get("E2N_ASSEMBLY_SUBPROCESS", "1").lower() not in {
                "0",
                "false",
                "no",
                "off",
            }
            if use_subprocess and backend in {"mace", "fairchem"}:
                assembled = services.run_assemble_subprocess(spec, mode="single")
                if assembled.get("status") != "success":
                    return jsonify(assembled), 500
                result = services.assembly_from_dict(assembled["results"][0])
            else:
                assembler = NanozymeAssembler()
                result = assembler.assemble(spec)
            services.design_results[result.job_id] = result
            services.persist_design_result(result)

            from nanozyme_mining.design.dopant_modifier import result_to_pdb_string

            return jsonify(
                {
                    "status": "success",
                    "job_id": result.job_id,
                    "xyz": result.xyz,
                    "pdb": result_to_pdb_string(result),
                    "smiles": result.smiles,
                    "score": services.score_payload(result.score),
                    "passed": result.score.passed_hard_constraints,
                    "errors": result.score.errors,
                    "multi_metal_mode": spec.multi_metal_mode,
                    "atom_count": len(result.atoms),
                }
            )
        except Exception as e:
            return error_response("failed to assemble nanozyme design", exc=e)

    @require_json_csrf
    def api_design_assemble_variants():
        try:
            from nanozyme_mining.design.dopant_modifier import result_to_pdb_string

            data = request.get_json(silent=True)
            if not isinstance(data, dict):
                return jsonify(
                    {"status": "error", "error": "request body must be a JSON object"}
                ), 400
            spec = DesignSpec.from_dict(data)
            if not spec.metals:
                return jsonify(
                    {"status": "error", "error": "at least one metal site is required"}
                ), 400
            if len(spec.metals) > 2:
                return jsonify(
                    {"status": "error", "error": "at most two metal sites are supported"}
                ), 400

            backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
            use_subprocess = os.environ.get("E2N_ASSEMBLY_SUBPROCESS", "1").lower() not in {
                "0",
                "false",
                "no",
                "off",
            }
            if use_subprocess and backend in {"mace", "fairchem"}:
                assembled = services.run_assemble_subprocess(spec, mode="variants")
                if assembled.get("status") != "success":
                    return jsonify(assembled), 500
                batch = [
                    services.assembly_from_dict(item)
                    for item in assembled["results"]
                ]
            else:
                assembler = NanozymeAssembler()
                batch = assembler.assemble_batch(spec)

            variants_out = {}
            variant_order = []
            for result in batch:
                services.design_results[result.job_id] = result
                services.persist_design_result(result)
                variant_order.append(result.job_id)
                variants_out[result.job_id] = {
                    "label": result.label,
                    "pdb": result_to_pdb_string(result),
                    "xyz": result.xyz,
                    "atom_count": len(result.atoms),
                    "score": services.score_payload(result.score),
                    "passed": result.score.passed_hard_constraints,
                    "errors": result.score.errors,
                }

            first = batch[0]
            return jsonify(
                {
                    "status": "success",
                    "job_id": first.job_id,
                    "total": len(batch),
                    "variant_order": variant_order,
                    "variants": variants_out,
                }
            )
        except Exception as e:
            return error_response("failed to assemble nanozyme variants", exc=e)

    @require_json_csrf
    def api_design_catalysis_screen(job_id):
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
            nanozyme_type = data.get("nanozyme_type") or result.design_spec.nanozyme_type
            task = services.get_reaction_task(nanozyme_type)
            max_poses = clamp_int(
                data.get("max_adsorption_poses", 6),
                default=6,
                min_value=1,
                max_value=20,
            )
            run_neb = parse_bool(data.get("run_neb"), default=False)
            backend = os.environ.get("E2N_MLP_BACKEND", "geometry_proxy").strip().lower()
            use_subprocess = os.environ.get("E2N_CATALYSIS_SUBPROCESS", "1").lower() not in {
                "0",
                "false",
                "no",
                "off",
            }
            if use_subprocess and backend in {"mace", "fairchem"}:
                payload = services.run_catalysis_screen_subprocess(result, data)
            else:
                payload = screen_catalysis(
                    result,
                    task=task,
                    max_adsorption_poses=max_poses,
                    run_neb=run_neb,
                    neb_steps=clamp_int(
                        data.get("neb_steps", 80),
                        default=80,
                        min_value=5,
                        max_value=500,
                    ),
                    neb_fmax=clamp_float(
                        data.get("neb_fmax", 0.08),
                        default=0.08,
                        min_value=0.01,
                        max_value=1.0,
                    ),
                    run_reaction_scan=parse_bool(
                        data.get("run_reaction_scan"),
                        default=False,
                    ),
                    reaction_scan_points=clamp_int(
                        data.get("reaction_scan_points", 5),
                        default=5,
                        min_value=3,
                        max_value=25,
                    ),
                )
            status_code = 200 if payload.get("status") == "success" else 400
            return jsonify(payload), status_code
        except Exception as e:
            return error_response("failed to screen catalysis", exc=e)

    def api_design_download(job_id, fmt):
        if not is_valid_motif_id(job_id):
            return jsonify({"status": "error", "error": "invalid job_id"}), 400
        if not is_valid_download_format(fmt):
            return jsonify({"status": "error", "error": f"format {fmt} not supported"}), 400
        result = services.design_results.get(job_id)
        if not result:
            try:
                persisted_path = safe_join(
                    services.design_result_dir(job_id),
                    f"{job_id}.{fmt}",
                )
            except ValueError:
                return jsonify({"status": "error", "error": "invalid job_id"}), 400
            if persisted_path.exists():
                return send_file(
                    persisted_path,
                    as_attachment=True,
                    download_name=f"nanozyme_{job_id}.{fmt}",
                )
            return jsonify({"status": "error", "error": "job not found"}), 404

        out_dir = tempfile.mkdtemp(prefix="e2n_dl_")

        @after_this_request
        def _cleanup_tempdir(response):
            try:
                shutil.rmtree(out_dir, ignore_errors=True)
            except Exception:
                pass
            return response

        paths = structure_exporter.export(result, out_dir)
        path = paths.get(fmt)
        if not path:
            return jsonify({"status": "error", "error": f"format {fmt} not supported"}), 400
        return send_file(
            path,
            as_attachment=True,
            download_name=f"nanozyme_{job_id}.{fmt}",
        )

    app.add_url_rule(
        "/api/design/assemble",
        endpoint="api_design_assemble",
        view_func=api_design_assemble,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/design/assemble_variants",
        endpoint="api_design_assemble_variants",
        view_func=api_design_assemble_variants,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/design/catalysis_screen/<job_id>",
        endpoint="api_design_catalysis_screen",
        view_func=api_design_catalysis_screen,
        methods=["POST"],
    )
    app.add_url_rule(
        "/api/design/download/<job_id>/<fmt>",
        endpoint="api_design_download",
        view_func=api_design_download,
    )
