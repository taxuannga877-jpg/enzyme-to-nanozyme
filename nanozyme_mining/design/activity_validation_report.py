"""Reporting utilities for designed-nanozyme activity validation."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from .physchem_knowledge import get_screening_proxy_policy
from .scientific_audit import profile_from_payload
from ..utils.constants import CATALYTIC_METAL_ELEMENTS


def summarize_activity_validation(assembly, activity_payloads: Sequence[Dict]) -> Dict:
    """Convert raw catalysis-screen payloads into a compact UI/report schema."""
    score = getattr(assembly, "score", None)
    score_details = getattr(score, "details", {}) if score is not None else {}
    structure_relaxation = dict(score_details.get("structure_relaxation") or {})
    activity_results = []
    for entry in activity_payloads:
        activity = entry.get("activity") or entry.get("nanozyme_type") or "Unknown"
        payload = entry.get("payload") or {}
        task = payload.get("task") or {}
        calculation = task.get("calculation") or {}
        poses = payload.get("adsorption_candidates") or []
        best = poses[0] if poses else {}
        reaction_profile = payload.get("reaction_profile") or {}
        redox_profile = payload.get("redox_state_profile") or {}
        mechanism_visualization = payload.get("mechanism_visualization") or {}
        adsorption_structure = payload.get("best_adsorption_structure") or {}
        transition_state = payload.get("transition_state") or {}
        profile_values = reaction_profile.get("relative_energies_ev") or []
        redox_values = redox_profile.get("relative_energies_ev") or []
        active_profile = profile_from_payload(payload)
        barrier = reaction_profile.get("proxy_barrier_ev")
        redox_activation = redox_profile.get("redox_activation_energy_ev")
        redox_energy_span = redox_profile.get("redox_energy_span_ev")
        activation_metric = barrier if barrier is not None else redox_activation
        activation_label = (
            "hydrolysis coordinate proxy barrier"
            if barrier is not None
            else "redox forward scan peak"
            if redox_activation is not None
            else None
        )
        result = {
            "activity": activity,
                "status": payload.get("status", "error"),
                "calculation_status": payload.get("calculation_status"),
                "error": payload.get("error"),
                "task_id": task.get("task_id"),
                "assay": task.get("assay"),
                "mechanism_family": calculation.get("mechanism_family"),
                "barrier_method": calculation.get("barrier_method"),
                "validation_level": calculation.get("validation_level", "screening_proxy"),
                "requires_charge": bool(calculation.get("requires_charge", False)),
                "requires_spin": bool(calculation.get("requires_spin", False)),
                "ml_backend": payload.get("ml_backend"),
                "ml_task": payload.get("ml_task"),
                "method_decision": payload.get("method_decision") or {},
                "active_center": payload.get("active_center") or {},
                "structure_electronic_state": payload.get("structure_electronic_state") or {},
                "pose_count": len(poses),
                "best_candidate_id": best.get("candidate_id"),
                "best_adsorption_energy_ev": _finite_or_none(best.get("adsorption_energy_ev")),
                "best_distance_score": _finite_or_none(best.get("distance_score")),
                "best_min_surface_distance_a": _finite_or_none(best.get("min_substrate_surface_distance")),
                "max_force_ev_per_a": _finite_or_none(best.get("max_force_ev_per_a")),
                "adsorption_local_optimization": best.get("local_optimization") or {},
                "transition_state_status": transition_state.get("status"),
                "transition_state_label": transition_state.get("label"),
                "transition_state_reason": transition_state.get("reason"),
                "reaction_profile_status": reaction_profile.get("status"),
                "reaction_profile_relative_energies_ev": [
                    _finite_or_none(v) for v in profile_values
                ],
                "proxy_barrier_ev": _finite_or_none(barrier),
                "redox_state_profile_status": redox_profile.get("status"),
                "redox_state_profile_relative_energies_ev": [
                    _finite_or_none(v) for v in redox_values
                ],
                "redox_activation_energy_ev": _finite_or_none(redox_activation),
                "redox_energy_span_ev": _finite_or_none(redox_energy_span),
                "activation_metric_ev": _finite_or_none(activation_metric),
                "activation_metric_label": activation_label,
                "mechanism_visualization": mechanism_visualization,
                "best_adsorption_structure": {
                    "candidate_id": adsorption_structure.get("candidate_id"),
                    "components": adsorption_structure.get("components") or [],
                    "atoms": _compact_structure_atoms(
                        adsorption_structure.get("atoms") or []
                    ),
                },
                "reaction_scan_gate": payload.get("reaction_scan_gate") or {},
                "interpretation": _interpret_activity(payload, best, reaction_profile),
            "scan_quality": active_profile.get("scan_quality") or {},
        }
        eligible, blockers = _figure_eligibility(
            structure_relaxation,
            result,
        )
        result["figure_eligible"] = eligible
        result["figure_blockers"] = blockers
        activity_results.append(result)

    successful = [r for r in activity_results if r["status"] == "success"]
    finite_adsorption = [
        r["best_adsorption_energy_ev"]
        for r in successful
        if r["best_adsorption_energy_ev"] is not None
    ]
    finite_barriers = [
        r["proxy_barrier_ev"]
        for r in successful
        if r["proxy_barrier_ev"] is not None
    ]
    finite_redox = [
        r["redox_activation_energy_ev"]
        for r in successful
        if r["redox_activation_energy_ev"] is not None
    ]
    finite_activation = [
        r["activation_metric_ev"]
        for r in successful
        if r["activation_metric_ev"] is not None
    ]
    summary = {
        "job_id": getattr(assembly, "job_id", ""),
        "label": getattr(assembly, "label", ""),
        "atom_count": len(getattr(assembly, "atoms", []) or []),
        "formal_charge": getattr(assembly, "formal_charge", 0),
        "spin_multiplicities": list(getattr(assembly, "spin_multiplicities", []) or []),
        "chemistry_warnings": list(getattr(assembly, "chemistry_warnings", []) or []),
        "structure_relaxation": structure_relaxation,
        "structure_atoms": _compact_structure_atoms(
            getattr(assembly, "atoms", []) or []
        ),
        "activity_count": len(activity_results),
        "completed_activity_count": len(successful),
        "mean_best_adsorption_energy_ev": _mean(finite_adsorption),
        "max_proxy_barrier_ev": max(finite_barriers) if finite_barriers else None,
        "max_redox_activation_energy_ev": max(finite_redox) if finite_redox else None,
        "max_activation_metric_ev": max(finite_activation) if finite_activation else None,
        "activity_results": activity_results,
        "caveat": (
            "These values are MACE + tblite screening proxies for prioritization. "
            "They should be compared only within the same task, charge/spin, and reference protocol."
        ),
    }
    return summary


def write_activity_validation_report(summary: Dict, output_dir: str | Path) -> str:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    path = output / "activity_validation_report.json"
    path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return str(path)


def render_activity_validation_figures(summary: Dict, output_dir: str | Path) -> List[Dict]:
    """Render publication-style PNG/SVG artifacts for the validation result."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:  # pragma: no cover - exercised only on missing matplotlib
        return [
            {
                "kind": "error",
                "label": "figure rendering unavailable",
                "error": str(exc),
            }
        ]

    _set_publication_style(plt)
    overview_png, overview_svg = _render_overview(summary, output, plt, np)
    artifacts = [
        {
            "kind": "overview",
            "label": "Activity validation overview",
            "png": Path(overview_png).name,
            "svg": Path(overview_svg).name,
        }
    ]
    if summary.get("structure_atoms"):
        structure_png, structure_svg = _render_relaxed_structure(
            summary,
            output,
            plt,
            np,
        )
        artifacts.append(
            {
                "kind": "relaxed_structure",
                "label": "MACE-relaxed dual-site structure",
                "png": Path(structure_png).name,
                "svg": Path(structure_svg).name,
            }
        )
    activity_rows = summary.get("activity_results", [])
    has_eligible_rows = any(row.get("figure_eligible", True) for row in activity_rows)
    has_mechanism_frames = any(
        (row.get("mechanism_visualization") or {}).get("frames")
        for row in activity_rows
    )
    has_adsorption_structures = any(
        (row.get("best_adsorption_structure") or {}).get("atoms")
        for row in activity_rows
    )
    if not has_eligible_rows and not has_mechanism_frames and not has_adsorption_structures:
        return artifacts
    if has_eligible_rows:
        profile_png, profile_svg = _render_reaction_profiles(summary, output, plt, np)
        artifacts.append(
            {
                "kind": "reaction_profile",
                "label": "Reaction-coordinate screening profiles",
                "png": Path(profile_png).name,
                "svg": Path(profile_svg).name,
            }
        )
    if has_mechanism_frames:
        artifacts.extend(_render_mechanism_panels(summary, output, plt, np))
        artifacts.extend(_render_mechanism_structure_frames(summary, output, plt, np))
    artifacts.extend(_render_adsorption_structure_frames(summary, output, plt, np))
    if any(
        row.get("figure_eligible", True)
        and (row.get("mechanism_visualization") or {}).get("kind")
        == "redox_electronic_state_path"
        for row in summary.get("activity_results", [])
    ):
        redox_png, redox_svg = _render_redox_state_map(summary, output, plt, np)
        artifacts.append({
            "kind": "redox_state_map",
            "label": "Redox charge/spin state map",
            "png": Path(redox_png).name,
            "svg": Path(redox_svg).name,
        })
    volcano_rows = _volcano_eligible_rows(summary)
    if volcano_rows:
        volcano_png, volcano_svg = _render_adsorption_volcano(
            {**summary, "activity_results": volcano_rows},
            output,
            plt,
            np,
        )
        artifacts.append({
            "kind": "adsorption_volcano",
            "label": "Adsorption-activation volcano",
            "png": Path(volcano_png).name,
            "svg": Path(volcano_svg).name,
        })
    return artifacts


def _render_overview(summary: Dict, output: Path, plt, np) -> tuple[str, str]:
    rows = summary.get("activity_results", [])
    labels = [r.get("activity") or "-" for r in rows] or ["No activity"]
    adsorption = [r.get("best_adsorption_energy_ev") for r in rows] or [None]
    barriers = [r.get("activation_metric_ev") for r in rows] or [None]
    y = np.arange(len(labels))

    fig, axes = plt.subplots(1, 2, figsize=(7.2, max(2.8, 0.55 * len(labels) + 1.4)))
    fig.patch.set_facecolor("#fbfaf6")
    policy = get_screening_proxy_policy()
    adsorption_optimum = policy["sabatier_adsorption_optimum_ev"]

    ax = axes[0]
    ads_values = [0.0 if v is None else v for v in adsorption]
    ads_colors = ["#b8c1cc" if v is None else "#2a9d8f" for v in adsorption]
    ax.axvspan(-2.3, 0.2, color="#2a9d8f", alpha=0.08, lw=0)
    ax.axvline(adsorption_optimum, color="#264653", lw=1.0, ls="--", alpha=0.65)
    ax.barh(y, ads_values, color=ads_colors, edgecolor="white", height=0.56)
    for yi, value in zip(y, adsorption):
        label = "not evaluated" if value is None else f"{value:.2f} eV"
        x = 0.05 if value is None else value
        ha = "left" if value is None or value <= 0 else "left"
        ax.text(x, yi, f"  {label}", va="center", ha=ha, fontsize=7, color="#3a4256")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("best adsorption energy (eV)")
    ax.set_title("Substrate Binding", loc="left")
    _symmetric_xlim(ax, ads_values + [adsorption_optimum], margin=0.35)

    ax = axes[1]
    barrier_values = [0.0 if v is None else v for v in barriers]
    barrier_colors = ["#b8c1cc" if v is None else "#e76f51" for v in barriers]
    ax.barh(y, barrier_values, color=barrier_colors, edgecolor="white", height=0.56)
    for yi, value in zip(y, barriers):
        label = "scan unavailable" if value is None else f"{value:.2f} eV"
        ax.text((0.02 if value is None else value) + 0.02, yi, label,
                va="center", ha="left", fontsize=7, color="#3a4256")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel("activation metric (eV)")
    ax.set_title("Mechanism Scan", loc="left")
    finite_barriers = [v for v in barrier_values if v and math.isfinite(v)]
    ax.set_xlim(0, max(finite_barriers + [1.0]) * 1.25)

    fig.suptitle(
        f"{summary.get('label') or summary.get('job_id')} — activity validation",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
        color="#264653",
    )
    fig.text(
        0.02,
        0.015,
        summary.get("caveat", ""),
        fontsize=7,
        color="#7c7a72",
    )
    fig.tight_layout(rect=[0.0, 0.05, 1.0, 0.92])
    png = output / "activity_validation_overview.png"
    svg = output / "activity_validation_overview.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_relaxed_structure(summary: Dict, output: Path, plt, np) -> tuple[str, str]:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")
    _draw_structure_projection(
        ax,
        summary.get("structure_atoms") or [],
        np,
        focus=False,
        label_metals=True,
    )
    relaxation = summary.get("structure_relaxation") or {}
    ax.set_title(
        f"{summary.get('label') or 'Dual-site candidate'}\n"
        f"MACE {relaxation.get('relaxation_status', 'unknown')} | "
        f"{summary.get('atom_count', 0)} atoms",
        loc="left",
    )
    fig.tight_layout()
    png = output / "mace_relaxed_structure.png"
    svg = output / "mace_relaxed_structure.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_reaction_profiles(summary: Dict, output: Path, plt, np) -> tuple[str, str]:
    rows = [
        row
        for row in summary.get("activity_results", [])
    ]
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    fig.patch.set_facecolor("#fbfaf6")
    palette = ["#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#0072b2"]
    plotted = 0
    for idx, row in enumerate(rows):
        values = row.get("reaction_profile_relative_energies_ev") or []
        profile_kind = "coordinate"
        if not any(v is not None for v in values):
            values = row.get("redox_state_profile_relative_energies_ev") or []
            profile_kind = "redox state"
        finite = [v for v in values if v is not None and math.isfinite(v)]
        if not finite or len(values) < 2:
            continue
        x = np.linspace(0.0, 1.0, len(values))
        y = np.array([np.nan if v is None else v for v in values], dtype=float)
        color = palette[idx % len(palette)]
        ax.plot(
            x,
            y,
            marker="o",
            lw=1.8,
            ms=4,
            ls="--" if profile_kind == "redox state" else "-",
            color=color,
            label=f"{row.get('activity')} ({profile_kind})",
        )
        ymax = np.nanmax(y)
        xmax = x[int(np.nanargmax(y))]
        ax.scatter([xmax], [ymax], s=42, color=color, edgecolor="white", zorder=4)
        plotted += 1

    if plotted:
        ax.set_xlabel("normalized reaction coordinate")
        ax.set_ylabel("relative energy (eV)")
        ax.set_title("Mechanism-Specific Screening Profiles", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    else:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.56,
            "Mechanism-specific scans were not available for this validation run.",
            ha="center",
            va="center",
            fontsize=11,
            color="#264653",
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.43,
            "Enable tblite/GFN2-xTB and mechanism scan to generate energy profiles.",
            ha="center",
            va="center",
            fontsize=6.6,
            color="#7c7a72",
        )

    fig.tight_layout()
    png = output / "reaction_coordinate_profiles.png"
    svg = output / "reaction_coordinate_profiles.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_mechanism_panels(summary: Dict, output: Path, plt, np) -> List[Dict]:
    rows = [
        row
        for row in summary.get("activity_results", [])
    ]
    panel_rows = [
        row for row in rows
        if (row.get("mechanism_visualization") or {}).get("frames")
    ]
    artifacts = []
    if not panel_rows:
        fig, ax = plt.subplots(figsize=(7.2, 3.4))
        fig.patch.set_facecolor("#fbfaf6")
        ax.set_axis_off()
        ax.text(
            0.5,
            0.56,
            "No mechanism frame summaries were produced.",
            ha="center",
            va="center",
            fontsize=11,
            color="#264653",
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.43,
            "Run a hydrolysis coordinate scan or redox state scan to populate this panel.",
            ha="center",
            va="center",
            fontsize=8,
            color="#7c7a72",
        )
        png = output / "mechanism_coordinate_panel_empty.png"
        svg = output / "mechanism_coordinate_panel_empty.svg"
        fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
        plt.close(fig)
        return [{
            "kind": "mechanism_panel",
            "label": "Mechanism-specific coordinate panel",
            "png": png.name,
            "svg": svg.name,
        }]

    for row in panel_rows:
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        fig.patch.set_facecolor("#fbfaf6")
        viz = row.get("mechanism_visualization") or {}
        frames = viz.get("frames") or []
        finite_frames = [
            frame for frame in frames
            if frame.get("relative_energy_ev") is not None
            and math.isfinite(float(frame.get("relative_energy_ev")))
        ]
        finite_frames = _representative_mechanism_frames(finite_frames)
        if not finite_frames:
            ax.set_axis_off()
            ax.text(0.5, 0.5, f"{row.get('activity')}: no finite mechanism path", ha="center", va="center")
        else:
            energies = [float(frame.get("relative_energy_ev")) for frame in finite_frames]
            x = np.arange(len(finite_frames), dtype=float)
            labels = _mechanism_node_labels(viz, row, len(finite_frames))
            span = _energy_span(energies)
            ymin = min(energies)
            ymax = max(energies)
            ax.set_xlim(-0.5, max(1.0, len(x) - 0.5))
            ax.set_ylim(ymin - 0.18 * span, ymax + 1.18 * span)

            for idx, (xi, yi) in enumerate(zip(x, energies)):
                ax.plot([xi - 0.18, xi + 0.18], [yi, yi], color="#111111", lw=1.5, solid_capstyle="round")
                if idx < len(x) - 1:
                    ax.plot(
                        [xi + 0.18, x[idx + 1] - 0.18],
                        [yi, energies[idx + 1]],
                        color="#111111",
                        lw=1.1,
                        ls=(0, (4, 3)),
                    )
                ax.text(
                    xi,
                    yi - 0.045 * span,
                    f"{yi:.3f}",
                    ha="center",
                    va="top",
                    fontsize=7,
                    color="#4b5563",
                )
                point_axes = ax.transAxes.inverted().transform(
                    ax.transData.transform((xi, yi + 0.11 * span))
                )
                inset = ax.inset_axes(
                    [point_axes[0] - 0.105, point_axes[1], 0.21, 0.29],
                    transform=ax.transAxes,
                    zorder=6,
                )
                _draw_ball_stick_thumbnail(inset, finite_frames[idx], np)

            activation = row.get("activation_metric_ev")
            if activation is not None and math.isfinite(float(activation)) and len(energies) >= 2:
                peak_idx = int(np.nanargmax(energies))
                peak = energies[peak_idx]
                if viz.get("kind") == "redox_electronic_state_path":
                    ax.text(
                        0.985,
                        0.965,
                        f"forward scan peak: {float(activation):.3f} eV",
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=7.5,
                        color="#e76f51",
                        fontweight="bold",
                    )
                else:
                    ax.annotate(
                        f"proxy barrier {float(activation):.3f} eV",
                        xy=(x[peak_idx], peak),
                        xytext=(x[peak_idx], peak + 0.34 * span),
                        ha="center",
                        color="#e76f51",
                        fontsize=7.5,
                        fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color="#e76f51", lw=0.8),
                    )
            ax.set_title(
                f"{row.get('activity')} — "
                f"{_display_label(viz.get('label') or viz.get('kind', 'mechanism'))}"
                f"{_diagnostic_title_suffix(row)}",
                loc="left",
            )
            ax.set_xlabel("constrained scan progress")
            ax.set_ylabel("relative energy (eV)")
            ax.set_xticks(x)
            ax.set_xticklabels(labels)
            ax.grid(axis="x", visible=False)

        fig.tight_layout()
        stem = f"{_filename_slug(row.get('activity'))}_mechanism_coordinate"
        png = output / f"{stem}.png"
        svg = output / f"{stem}.svg"
        fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
        plt.close(fig)
        artifacts.append({
            "kind": "mechanism_panel",
            "activity": row.get("activity"),
            "label": f"{row.get('activity')} mechanism coordinate",
            "png": png.name,
            "svg": svg.name,
        })
    return artifacts


def _render_mechanism_structure_frames(
    summary: Dict,
    output: Path,
    plt,
    np,
) -> List[Dict]:
    rows = [
        row
        for row in summary.get("activity_results", [])
        if (row.get("mechanism_visualization") or {}).get("frames")
    ]
    artifacts = []
    if not rows:
        return artifacts
    for row in rows:
        all_frames = (row.get("mechanism_visualization") or {}).get("frames") or []
        frames = _representative_mechanism_frames(all_frames)
        fig, axes = plt.subplots(
            1,
            max(len(frames), 1),
            figsize=(2.25 * max(len(frames), 1), 2.75),
            squeeze=False,
        )
        fig.patch.set_facecolor("#ffffff")
        view = _reference_structure_view(frames, np)
        bounds = _common_structure_projection_bounds(
            frames,
            view,
            np,
            focus=True,
        )
        for col_idx, frame in enumerate(frames):
            ax = axes[0][col_idx]
            _draw_structure_projection(
                ax,
                frame.get("atoms") or [],
                np,
                focus=True,
                label_metals=col_idx in {0, len(frames) - 1},
                view=view,
                bounds=bounds,
                show_regions=False,
                show_component_labels=False,
            )
            energy = frame.get("relative_energy_ev")
            coords = frame.get("coordinates") or []
            distance = coords[0].get("distance_a") if coords else None
            title = f"{row.get('activity')} | {col_idx + 1}/{len(frames)}"
            if energy is not None:
                title += f"\ndE {float(energy):.2f} eV"
            if distance is not None:
                title += f" | d {float(distance):.2f} A"
            ax.set_title(title, fontsize=8)
        if not row.get("figure_eligible", True):
            quality = row.get("scan_quality") or {}
            fig.suptitle(
                "Diagnostic preview — "
                f"{quality.get('converged_frame_count', 0)}/"
                f"{quality.get('frame_count', len(frames))} frames converged",
                x=0.01,
                ha="left",
                fontsize=8,
                color="#c45a43",
                fontweight="bold",
            )
            fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.91])
        else:
            fig.tight_layout()
        stem = f"{_filename_slug(row.get('activity'))}_mechanism_structures"
        png = output / f"{stem}.png"
        svg = output / f"{stem}.svg"
        fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
        plt.close(fig)
        artifacts.append({
            "kind": "mechanism_structure_frames",
            "activity": row.get("activity"),
            "label": f"{row.get('activity')} optimized mechanism structures",
            "png": png.name,
            "svg": svg.name,
        })
    return artifacts


def _representative_mechanism_frames(
    frames: Sequence[Dict],
    max_frames: int = 5,
) -> List[Dict]:
    """Keep calculation figures readable while preserving full frames in JSON.

    Some coupled hydrolysis coordinates contain a 5 x 5 grid (25 optimized
    structures).  The full grid remains in the calculation payload, but the
    publication/gallery panel should show five evenly spaced representatives,
    matching the compact Initial/25%/50%/75%/Final visual language.
    """
    frames = list(frames)
    if len(frames) <= max_frames:
        return frames
    if max_frames <= 1:
        return [frames[0]]
    indices = [
        round(index * (len(frames) - 1) / (max_frames - 1))
        for index in range(max_frames)
    ]
    return [frames[index] for index in indices]


representative_mechanism_frames = _representative_mechanism_frames


def _render_adsorption_structure_frames(
    summary: Dict,
    output: Path,
    plt,
    np,
) -> List[Dict]:
    rows = [
        row
        for row in summary.get("activity_results", [])
        if not (row.get("mechanism_visualization") or {}).get("frames")
        and (row.get("best_adsorption_structure") or {}).get("atoms")
    ]
    artifacts = []
    for row in rows:
        atoms = _adsorption_focus_atoms(
            (row.get("best_adsorption_structure") or {}).get("atoms") or [],
            np,
        )
        frame = {"atoms": atoms}
        view = _metal_axis_structure_view(atoms, np) or _reference_structure_view([frame], np)
        bounds = _common_structure_projection_bounds(
            [frame],
            view,
            np,
            focus=True,
        )
        fig, ax = plt.subplots(figsize=(7.2, 3.8))
        fig.patch.set_facecolor("#ffffff")
        _draw_structure_projection(
            ax,
            atoms,
            np,
            focus=True,
            label_metals=True,
            view=view,
            bounds=bounds,
            show_regions=False,
            show_component_labels=False,
        )
        adsorption = row.get("best_adsorption_energy_ev")
        energy_label = (
            f"Eads {float(adsorption):.2f} eV"
            if adsorption is not None
            else "Eads unavailable"
        )
        ax.set_title(
            f"{row.get('activity')} — adsorption pose | {energy_label}",
            loc="left",
            fontsize=10,
        )
        gate = row.get("reaction_scan_gate") or {}
        if gate.get("status") == "rejected":
            ax.text(
                0.99,
                0.98,
                "screened out",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7,
                color="#c45a43",
            )
            fig.tight_layout()
        else:
            fig.tight_layout()
        stem = f"{_filename_slug(row.get('activity'))}_adsorption_structure"
        png = output / f"{stem}.png"
        svg = output / f"{stem}.svg"
        fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
        plt.close(fig)
        artifacts.append({
            "kind": "adsorption_structure",
            "activity": row.get("activity"),
            "label": f"{row.get('activity')} adsorption structure",
            "png": png.name,
            "svg": svg.name,
        })
    return artifacts


def _adsorption_focus_atoms(atoms: Sequence[Dict], np) -> List[Dict]:
    metal_positions = [
        np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        for atom in atoms
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    if not metal_positions:
        return list(atoms)
    selected = []
    for atom in atoms:
        element = str(atom.get("element", "")).upper()
        if (
            element in CATALYTIC_METAL_ELEMENTS
            or atom.get("substrate_name")
            or atom.get("molecule_id")
        ):
            selected.append(atom)
            continue
        position = np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        nearest_metal = min(
            float(np.linalg.norm(position - metal))
            for metal in metal_positions
        )
        if nearest_metal <= 4.6:
            selected.append(atom)
    return selected


def _metal_axis_structure_view(atoms: Sequence[Dict], np) -> Dict | None:
    metals = [
        np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        for atom in atoms
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    if len(metals) < 2:
        return None
    center = np.mean(np.asarray(metals), axis=0)
    horizontal = np.asarray(metals[-1] - metals[0], dtype=float)
    horizontal_norm = float(np.linalg.norm(horizontal))
    if horizontal_norm <= 1e-12:
        return None
    horizontal = horizontal / horizontal_norm
    support = np.asarray(
        [
            atom.get("coords", (0.0, 0.0, 0.0))
            for atom in atoms
            if not atom.get("substrate_name")
            and not atom.get("molecule_id")
            and str(atom.get("element", "")).upper() not in CATALYTIC_METAL_ELEMENTS
        ],
        dtype=float,
    )
    try:
        _u, _s, vh = np.linalg.svd(support - np.mean(support, axis=0), full_matrices=False)
        surface_normal = np.asarray(vh[-1], dtype=float)
    except Exception:
        surface_normal = np.asarray([0.0, 0.0, 1.0], dtype=float)
    adsorbates = np.asarray(
        [
            atom.get("coords", (0.0, 0.0, 0.0))
            for atom in atoms
            if atom.get("substrate_name") or atom.get("molecule_id")
        ],
        dtype=float,
    )
    if len(adsorbates) and float(np.dot(np.mean(adsorbates, axis=0) - center, surface_normal)) < 0.0:
        surface_normal = -surface_normal
    horizontal = horizontal - float(np.dot(horizontal, surface_normal)) * surface_normal
    horizontal = horizontal / max(float(np.linalg.norm(horizontal)), 1e-12)
    in_plane_vertical = np.cross(surface_normal, horizontal)
    in_plane_vertical = in_plane_vertical / max(
        float(np.linalg.norm(in_plane_vertical)),
        1e-12,
    )
    tilt = math.radians(30.0)
    vertical = math.cos(tilt) * surface_normal + math.sin(tilt) * in_plane_vertical
    return {
        "center": center,
        "horizontal": horizontal,
        "vertical": vertical,
    }


def _render_redox_state_map(summary: Dict, output: Path, plt, np) -> tuple[str, str]:
    rows = [
        row
        for row in summary.get("activity_results", [])
        if row.get("figure_eligible", True)
    ]
    redox_rows = [
        row for row in rows
        if (row.get("mechanism_visualization") or {}).get("kind") == "redox_electronic_state_path"
    ]
    fig, axes = plt.subplots(
        max(len(redox_rows), 1),
        1,
        figsize=(7.2, max(2.8, 1.9 * max(len(redox_rows), 1))),
        squeeze=False,
    )
    fig.patch.set_facecolor("#fbfaf6")
    if not redox_rows:
        ax = axes[0][0]
        ax.set_axis_off()
        ax.text(
            0.5,
            0.55,
            "No redox state scans in this run.",
            ha="center",
            va="center",
            fontsize=11,
            color="#264653",
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.42,
            "Redox activities populate charge/spin heatmaps when tblite is available.",
            ha="center",
            va="center",
            fontsize=8,
            color="#7c7a72",
        )
    for ax, row in zip(axes[:, 0], redox_rows):
        viz = row.get("mechanism_visualization") or {}
        matrix = viz.get("relative_state_energy_matrix_ev") or []
        columns = viz.get("multiplicity_columns") or []
        if not matrix or not columns:
            ax.set_axis_off()
            ax.text(0.5, 0.5, f"{row.get('activity')}: no state matrix", ha="center", va="center")
            continue
        arr = np.array(
            [
                [np.nan if value is None else float(value) for value in row_values]
                for row_values in matrix
            ],
            dtype=float,
        )
        im = ax.imshow(arr, aspect="auto", cmap="viridis")
        ax.set_title(f"{row.get('activity')} redox state scan", loc="left")
        ax.set_xlabel("spin multiplicity")
        ax.set_ylabel("scan frame")
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels([str(col) for col in columns])
        ax.set_yticks(range(arr.shape[0]))
        ax.set_yticklabels([str(i) for i in range(arr.shape[0])])
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="relative energy (eV)")
    fig.tight_layout()
    png = output / "redox_state_heatmap.png"
    svg = output / "redox_state_heatmap.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _render_adsorption_volcano(summary: Dict, output: Path, plt, np) -> tuple[str, str]:
    rows = [
        row
        for row in summary.get("activity_results", [])
        if row.get("figure_eligible", True)
    ]
    points = [
        row for row in rows
        if row.get("best_adsorption_energy_ev") is not None
        and row.get("activation_metric_ev") is not None
    ]
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    fig.patch.set_facecolor("#fbfaf6")
    if points:
        colors = {
            "hydrolysis": "#2a9d8f",
            "peroxide_redox": "#e76f51",
            "peroxide_disproportionation": "#f4a261",
            "dioxygen_redox": "#315f8d",
            "redox_dismutation": "#8a5a44",
        }
        for row in points:
            family = row.get("mechanism_family") or "other"
            ax.scatter(
                [row["best_adsorption_energy_ev"]],
                [row["activation_metric_ev"]],
                s=70,
                color=colors.get(family, "#6b7280"),
                edgecolor="white",
                linewidth=0.8,
                label=family,
            )
            ax.text(
                row["best_adsorption_energy_ev"],
                row["activation_metric_ev"],
                f"  {row.get('activity')}",
                fontsize=7,
                va="center",
                color="#3a4256",
            )
        xs = [float(row["best_adsorption_energy_ev"]) for row in points]
        ys = [float(row["activation_metric_ev"]) for row in points]
        xmin, xmax = min(xs + [-1.6]), max(xs + [0.8])
        xgrid = np.linspace(xmin, xmax, 160)
        policy = get_screening_proxy_policy()
        optimum = policy["sabatier_adsorption_optimum_ev"]
        base = min(ys) if ys else 0.0
        scale = max((max(ys) - base) if ys else 1.0, 0.35)
        guide = base + scale * ((xgrid - optimum) / max(abs(xmax - xmin), 0.5) * 2.4) ** 2
        ax.plot(xgrid, guide, color="#cc79a7", lw=1.15, ls=":", alpha=0.9, label="Sabatier guide")
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="best")
        ax.axvline(optimum, color="#264653", lw=1.0, ls="--", alpha=0.55)
        ax.set_xlabel("best adsorption energy (eV)")
        ax.set_ylabel("activation metric (eV)")
        ax.set_title("Adsorption-Activation Volcano", loc="left")
    else:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.55,
            "Volcano plot needs both adsorption and activation metrics.",
            ha="center",
            va="center",
            fontsize=11,
            color="#264653",
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.42,
            "Run tblite adsorption plus hydrolysis/redox mechanism scans.",
            ha="center",
            va="center",
            fontsize=8,
            color="#7c7a72",
        )
    fig.tight_layout()
    png = output / "adsorption_activation_volcano.png"
    svg = output / "adsorption_activation_volcano.svg"
    fig.savefig(png, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(svg, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return str(png), str(svg)


def _set_publication_style(plt) -> None:
    font_family = "DejaVu Sans"
    try:
        from matplotlib import font_manager

        root = Path(__file__).resolve().parents[2]
        for font_path in [
            root / "fonts" / "NotoSansSC-VF.ttf",
            root / "fonts" / "NotoSansSC-Regular.otf",
        ]:
            if font_path.exists():
                font_manager.fontManager.addfont(str(font_path))
                font_family = font_manager.FontProperties(fname=str(font_path)).get_name()
                break
    except Exception:
        font_family = "DejaVu Sans"

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font_family, "DejaVu Sans", "Helvetica", "Arial"],
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.16,
            "grid.linestyle": "-",
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


_PLOT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "S": 1.05,
    "P": 1.07,
    "FE": 1.24,
    "CU": 1.17,
    "ZN": 1.25,
    "MN": 1.39,
    "CO": 1.18,
    "NI": 1.17,
}

_PLOT_COLORS = {
    "C": "#cbd1d6",
    "N": "#4a90d9",
    "O": "#e76f51",
    "S": "#d6a21f",
    "P": "#8a5a44",
    "CU": "#d55e00",
    "FE": "#0072b2",
    "ZN": "#7a6eb1",
    "MN": "#8c564b",
    "CO": "#6f4e7c",
    "NI": "#4d7c6f",
}


def _compact_structure_atoms(atoms: Sequence[Dict]) -> List[Dict]:
    compact = []
    for atom in atoms:
        row = {
            "element": str(atom.get("element", "")).upper(),
            "coords": [float(value) for value in atom.get("coords", (0.0, 0.0, 0.0))],
        }
        for key in (
            "atom_name",
            "residue_name",
            "site_id",
            "substrate_name",
            "molecule_id",
            "is_coord_atom",
            "is_embedded_metal",
        ):
            if atom.get(key) is not None:
                row[key] = atom.get(key)
        compact.append(row)
    return compact


def _draw_structure_projection(
    ax,
    atoms: Sequence[Dict],
    np,
    *,
    focus: bool,
    label_metals: bool,
    view: Dict | None = None,
    bounds: Sequence[float] | None = None,
    show_regions: bool = True,
    show_contact_guide: bool = True,
    show_component_labels: bool = True,
) -> None:
    selected, coords, projection = _project_structure_atoms(
        atoms,
        np,
        focus=focus,
        view=view,
    )
    if not selected:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No structure coordinates", ha="center", va="center")
        return

    for left in range(len(selected)):
        left_atom = selected[left][1]
        left_element = str(left_atom.get("element", "")).upper()
        left_component = _plot_component(left_atom)
        for right in range(left + 1, len(selected)):
            right_atom = selected[right][1]
            right_element = str(right_atom.get("element", "")).upper()
            right_component = _plot_component(right_atom)
            distance = float(np.linalg.norm(coords[right] - coords[left]))
            is_metal_bond = (
                left_element in CATALYTIC_METAL_ELEMENTS or right_element in CATALYTIC_METAL_ELEMENTS
            ) and (
                left_element in {"N", "O", "S"}
                or right_element in {"N", "O", "S"}
            )
            same_substrate_molecule = (
                left_component in {"tmb", "oxidant"}
                and left_component == right_component
                and left_atom.get("molecule_id") == right_atom.get("molecule_id")
            )
            same_support = left_component == right_component == "support"
            if not is_metal_bond and not same_substrate_molecule and not same_support:
                continue
            threshold = (
                2.55
                if is_metal_bond
                else 1.22
                * (
                    _PLOT_RADII.get(left_element, 0.8)
                    + _PLOT_RADII.get(right_element, 0.8)
                )
            )
            if 0.45 <= distance <= threshold:
                if is_metal_bond:
                    color, width = "#2a9d8f", 1.45
                elif left_component == "tmb":
                    color, width = "#59636e", 1.05
                elif left_component == "oxidant":
                    color, width = "#e4573d", 1.55
                else:
                    color, width = "#d8dde1", 0.55
                ax.plot(
                    [projection[left, 0], projection[right, 0]],
                    [projection[left, 1], projection[right, 1]],
                    color=color,
                    lw=width,
                    zorder=1,
                )

    if show_regions:
        _draw_component_regions(ax, selected, projection)

    for local, (_original, atom, _position) in enumerate(selected):
        element = str(atom.get("element", "")).upper()
        is_metal = element in CATALYTIC_METAL_ELEMENTS
        component = _plot_component(atom)
        if is_metal:
            color, size, edge_width = _PLOT_COLORS.get(element, "#9aa2aa"), 88, 0.55
        elif component == "tmb":
            organic_colors = {
                "C": "#65717d",
                "N": "#2f6fba",
                "O": "#cf5b46",
                "S": "#d6a21f",
                "P": "#8a5a44",
                "H": "#fff7f2",
            }
            color = organic_colors.get(element, "#65717d")
            size = 23 if element in {"N", "O", "S", "P"} else 15
            edge_width = 0.35
        elif component == "oxidant":
            if element == "H":
                color, size, edge_width = "#fff7f2", 18, 0.8
            else:
                color, size, edge_width = "#e4573d", 27, 0.45
        elif atom.get("is_coord_atom"):
            color, size, edge_width = "#4a90d9", 24, 0.25
        else:
            color, size, edge_width = _PLOT_COLORS.get(element, "#cfd5da"), 10, 0.2
        ax.scatter(
            [projection[local, 0]],
            [projection[local, 1]],
            s=size,
            color=color,
            edgecolor="white",
            linewidth=edge_width,
            zorder=3 if is_metal else 2,
        )
        if label_metals and is_metal:
            ax.text(
                projection[local, 0],
                projection[local, 1] + 0.28,
                element.title(),
                ha="center",
                va="bottom",
                fontsize=6.3,
                fontweight="bold",
                color="#111111",
                zorder=6,
            )
    if show_contact_guide:
        _draw_relative_position_guide(ax, selected, coords, projection, np)
    if show_component_labels:
        _draw_component_labels(ax, selected, projection)
    ax.set_aspect("equal")
    if bounds is None:
        xmin, xmax = float(np.min(projection[:, 0])), float(np.max(projection[:, 0]))
        ymin, ymax = float(np.min(projection[:, 1])), float(np.max(projection[:, 1]))
    else:
        xmin, xmax, ymin, ymax = [float(value) for value in bounds]
    xpad = max(0.06 * (xmax - xmin), 0.35)
    ypad = max(0.06 * (ymax - ymin), 0.35)
    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)
    ax.set_axis_off()


def _reference_structure_view(frames: Sequence[Dict], np) -> Dict | None:
    """Lock an oblique physical view: support tilted, adsorbates still above."""
    reference_atoms = next(
        (
            frame.get("atoms") or []
            for frame in frames
            if frame.get("atoms")
        ),
        [],
    )
    if not reference_atoms:
        return None

    support_positions = [
        np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        for atom in reference_atoms
        if str(atom.get("element", "")).upper() != "H"
        and not atom.get("substrate_name")
        and not atom.get("molecule_id")
        and str(atom.get("element", "")).upper() not in CATALYTIC_METAL_ELEMENTS
    ]
    if len(support_positions) < 3:
        return None
    support = np.asarray(support_positions, dtype=float)
    center = np.mean(support, axis=0)
    try:
        _u, _s, vh = np.linalg.svd(support - center, full_matrices=False)
        horizontal = np.asarray(vh[0], dtype=float)
        surface_normal = np.asarray(vh[-1], dtype=float)
    except Exception:
        horizontal = np.asarray([1.0, 0.0, 0.0], dtype=float)
        surface_normal = np.asarray([0.0, 0.0, 1.0], dtype=float)

    adsorbate_positions = [
        np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        for atom in reference_atoms
        if atom.get("substrate_name") or atom.get("molecule_id")
    ]
    if adsorbate_positions:
        adsorbate_center = np.mean(np.asarray(adsorbate_positions), axis=0)
        if float(np.dot(adsorbate_center - center, surface_normal)) < 0.0:
            surface_normal = -surface_normal
    dominant = int(np.argmax(np.abs(horizontal)))
    if float(horizontal[dominant]) < 0.0:
        horizontal = -horizontal
    horizontal = (
        horizontal
        - float(np.dot(horizontal, surface_normal)) * surface_normal
    )
    norm = float(np.linalg.norm(horizontal))
    if norm > 1e-12:
        horizontal = horizontal / norm
    in_plane_vertical = np.cross(surface_normal, horizontal)
    in_plane_norm = float(np.linalg.norm(in_plane_vertical))
    if in_plane_norm > 1e-12:
        in_plane_vertical = in_plane_vertical / in_plane_norm

    # A strict side view collapses the coordination pocket into one line.
    # Retain a 30-degree in-plane component so the metal-donor geometry remains
    # legible while the surface normal still points toward the adsorbates.
    tilt = math.radians(30.0)
    vertical = (
        math.cos(tilt) * surface_normal
        + math.sin(tilt) * in_plane_vertical
    )
    return {
        "center": center,
        "horizontal": horizontal,
        "vertical": vertical,
    }


def _substrate_plot_atoms(atoms: Sequence[Dict]) -> List[Dict]:
    """Keep only molecular reactants for compact energy-profile thumbnails."""
    return [
        atom
        for atom in atoms
        if atom.get("substrate_name") or atom.get("molecule_id")
    ]


def _plot_component(atom: Dict) -> str:
    substrate = str(atom.get("substrate_name") or "").upper()
    if substrate in {"O2", "H2O2", "SUPEROXIDE", "H2O"}:
        return "oxidant"
    if substrate or atom.get("molecule_id"):
        return "tmb"
    return "support"


def _display_formula(value: str) -> str:
    normalized = str(value or "").upper()
    if normalized == "O2":
        return r"$\mathrm{O_2}$"
    if normalized == "H2O2":
        return r"$\mathrm{H_2O_2}$"
    return str(value or "")


def _draw_component_regions(ax, selected, projection) -> None:
    try:
        from matplotlib.patches import Ellipse
    except Exception:
        return
    styles = {
        "tmb": ("#5b7188", 0.065),
        "oxidant": ("#e76f51", 0.10),
    }
    for component, (color, alpha) in styles.items():
        indices = [
            idx
            for idx, (_original, atom, _position) in enumerate(selected)
            if _plot_component(atom) == component
        ]
        if not indices:
            continue
        points = projection[indices]
        xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
        ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
        width = max(xmax - xmin + 0.55, 0.75)
        height = max(ymax - ymin + 0.45, 0.65)
        ax.add_patch(
            Ellipse(
                ((xmin + xmax) / 2, (ymin + ymax) / 2),
                width,
                height,
                facecolor=color,
                edgecolor=color,
                linewidth=0.8,
                linestyle=(0, (3, 2)),
                alpha=alpha,
                zorder=0,
            )
        )


def _draw_component_labels(ax, selected, projection) -> None:
    component_labels = {}
    for component in ("tmb", "oxidant"):
        indices = [
            idx
            for idx, (_original, atom, _position) in enumerate(selected)
            if _plot_component(atom) == component
        ]
        if not indices:
            continue
        points = projection[indices]
        names = [
            str(selected[idx][1].get("substrate_name") or "")
            for idx in indices
        ]
        label = "TMB" if component == "tmb" else _display_formula(next((n for n in names if n), "oxidant"))
        component_labels[component] = (
            label,
            float(points[:, 0].mean()),
            float(points[:, 1].max()),
        )
    for component, (label, x, y) in component_labels.items():
        color = "#425466" if component == "tmb" else "#c94f37"
        ax.text(
            x,
            y + 0.22,
            label,
            ha="center",
            va="bottom",
            fontsize=6.7,
            fontweight="bold",
            color=color,
            zorder=6,
        )


def _draw_relative_position_guide(ax, selected, coords, projection, np) -> None:
    tmb_indices = [
        idx
        for idx, (_original, atom, _position) in enumerate(selected)
        if _plot_component(atom) == "tmb"
    ]
    oxidant_indices = [
        idx
        for idx, (_original, atom, _position) in enumerate(selected)
        if _plot_component(atom) == "oxidant"
    ]
    if not tmb_indices or not oxidant_indices:
        return
    closest = min(
        (
            (float(np.linalg.norm(coords[tmb] - coords[oxidant])), tmb, oxidant)
            for tmb in tmb_indices
            for oxidant in oxidant_indices
        ),
        key=lambda item: item[0],
    )
    distance, tmb_idx, oxidant_idx = closest
    x_values = [projection[tmb_idx, 0], projection[oxidant_idx, 0]]
    y_values = [projection[tmb_idx, 1], projection[oxidant_idx, 1]]
    ax.plot(
        x_values,
        y_values,
        color="#c46a54",
        linewidth=0.75,
        linestyle=(0, (2, 2)),
        alpha=0.9,
        zorder=4,
    )


def _draw_ball_stick_thumbnail(ax, frame: Dict, np) -> None:
    """Render TMB and oxidant as separated ball-and-stick components."""
    atoms = frame.get("atoms") or []
    ax.set_axis_off()
    tmb_atoms = [atom for atom in atoms if _plot_component(atom) == "tmb"]
    oxidant_atoms = [atom for atom in atoms if _plot_component(atom) == "oxidant"]
    components = []
    if tmb_atoms and oxidant_atoms:
        components = [
            (ax.inset_axes([0.00, 0.02, 0.73, 0.96]), tmb_atoms),
            (ax.inset_axes([0.75, 0.22, 0.24, 0.58]), oxidant_atoms),
        ]
    elif tmb_atoms:
        components = [(ax.inset_axes([0.02, 0.04, 0.96, 0.92]), tmb_atoms)]
    elif oxidant_atoms:
        components = [(ax.inset_axes([0.08, 0.10, 0.84, 0.80]), oxidant_atoms)]
    for component_ax, component_atoms in components:
        _draw_structure_projection(
            component_ax,
            component_atoms,
            np,
            focus=False,
            label_metals=False,
            view=None,
            bounds=None,
            show_regions=False,
            show_contact_guide=False,
            show_component_labels=False,
        )


def _project_structure_atoms(
    atoms: Sequence[Dict],
    np,
    *,
    focus: bool,
    view: Dict | None,
):
    heavy = [
        (idx, atom)
        for idx, atom in enumerate(atoms)
        if str(atom.get("element", "")).upper() != "H"
        or (
            _plot_component(atom) == "oxidant"
            and str(atom.get("substrate_name") or "").upper() == "H2O2"
        )
    ]
    if not heavy:
        return [], np.empty((0, 3)), np.empty((0, 2))

    metal_positions = [
        np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        for _idx, atom in heavy
        if str(atom.get("element", "")).upper() in CATALYTIC_METAL_ELEMENTS
    ]
    focus_center = (
        np.mean(np.asarray(metal_positions), axis=0)
        if metal_positions
        else np.mean(
            np.asarray([atom.get("coords", (0.0, 0.0, 0.0)) for _idx, atom in heavy], dtype=float),
            axis=0,
        )
    )
    selected = []
    for idx, atom in heavy:
        element = str(atom.get("element", "")).upper()
        position = np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float)
        essential = (
            element in CATALYTIC_METAL_ELEMENTS
            or bool(atom.get("is_coord_atom"))
            or bool(atom.get("substrate_name"))
            or bool(atom.get("molecule_id"))
        )
        if not focus or essential or float(np.linalg.norm(position - focus_center)) <= 4.8:
            selected.append((idx, atom, position))
    if len(selected) < 2:
        selected = [
            (idx, atom, np.asarray(atom.get("coords", (0.0, 0.0, 0.0)), dtype=float))
            for idx, atom in heavy
        ]

    coords = np.asarray([position for _idx, _atom, position in selected], dtype=float)
    if view is not None:
        centered = coords - np.asarray(view["center"], dtype=float)
        projection = np.column_stack(
            [
                centered @ np.asarray(view["horizontal"], dtype=float),
                centered @ np.asarray(view["vertical"], dtype=float),
            ]
        )
    else:
        centered = coords - np.mean(coords, axis=0)
        try:
            _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
            projection = centered @ vh[:2].T
        except Exception:
            projection = centered[:, :2]
        if projection.ndim != 2 or projection.shape[1] < 2:
            projection = np.column_stack(
                [
                    projection[:, 0]
                    if projection.ndim == 2 and projection.shape[1]
                    else np.zeros(len(selected)),
                    np.zeros(len(selected)),
                ]
            )
    return selected, coords, projection


def _common_structure_projection_bounds(
    frames: Sequence[Dict],
    view: Dict | None,
    np,
    *,
    focus: bool,
) -> tuple[float, float, float, float] | None:
    projections = []
    for frame in frames:
        _selected, _coords, projection = _project_structure_atoms(
            frame.get("atoms") or [],
            np,
            focus=focus,
            view=view,
        )
        if len(projection):
            projections.append(projection)
    if not projections:
        return None
    merged = np.vstack(projections)
    return (
        float(np.min(merged[:, 0])),
        float(np.max(merged[:, 0])),
        float(np.min(merged[:, 1])),
        float(np.max(merged[:, 1])),
    )


def _display_label(value) -> str:
    return str(value or "").replace("_", " ").strip()


def _filename_slug(value) -> str:
    return "".join(
        char.lower() if char.isalnum() else "_"
        for char in str(value or "activity")
    ).strip("_")


def _diagnostic_title_suffix(row: Dict) -> str:
    if row.get("figure_eligible", True):
        return ""
    quality = row.get("scan_quality") or {}
    converged = quality.get("converged_frame_count")
    total = quality.get("frame_count")
    if converged is not None and total:
        return f" [diagnostic {converged}/{total} converged]"
    return " [diagnostic]"


def _mechanism_node_labels(viz: Dict, row: Dict, count: int) -> List[str]:
    if count == 1:
        return ["Initial"]
    if count == 5:
        return ["Initial", "25%", "50%", "75%", "Final"]
    base = ["Initial"] + [
        f"{100 * idx / (count - 1):.0f}%"
        for idx in range(1, count - 1)
    ] + ["Final"]
    if count <= len(base):
        return base[:count]
    return base + [f"S{i}" for i in range(len(base) + 1, count + 1)]


def _mechanism_path_colors(viz: Dict) -> Dict[str, str]:
    if viz.get("kind") == "redox_electronic_state_path":
        return {
            "metal": "#0072b2",
            "substrate": "#e76f51",
            "distance": "#cc79a7",
            "bond": "#2a9d8f",
        }
    return {
        "metal": "#d55e00",
        "substrate": "#4a90d9",
        "distance": "#2a9d8f",
        "bond": "#2a9d8f",
    }


def _energy_span(values: Sequence[float]) -> float:
    finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not finite:
        return 1.0
    return max(max(finite) - min(finite), 0.5)


def _draw_mechanism_thumbnail(ax, x: float, y: float, frame: Dict, colors: Dict[str, str], np) -> None:
    span = _energy_span([y, 0.0])
    y0 = y + 0.28 * span
    coords = frame.get("coordinates") or []
    primary = coords[0] if coords else {}
    action = str(primary.get("action") or "scan")
    distance = primary.get("distance_a")
    bond_alpha = 0.95
    if distance is not None:
        try:
            bond_alpha = max(0.28, min(0.95, 2.4 / max(float(distance), 0.2)))
        except Exception:
            bond_alpha = 0.75

    ax.scatter([x - 0.10], [y0], s=56, color=colors["metal"], edgecolor="white", linewidth=0.55, zorder=5, clip_on=False)
    ax.scatter([x + 0.08], [y0 + 0.06 * span], s=32, color=colors["substrate"], edgecolor="white", linewidth=0.45, zorder=5, clip_on=False)
    ax.scatter([x + 0.18], [y0 - 0.05 * span], s=26, color="#d0d6dc", edgecolor="white", linewidth=0.35, zorder=5, clip_on=False)
    ax.plot(
        [x - 0.06, x + 0.07],
        [y0 + 0.01 * span, y0 + 0.05 * span],
        color=colors["bond"],
        lw=1.2,
        alpha=bond_alpha,
        zorder=4,
        clip_on=False,
    )
    if action in {"break", "stretch"}:
        ax.plot([x + 0.10, x + 0.20], [y0 + 0.04 * span, y0 + 0.10 * span], color="#111111", lw=0.65, ls=":", zorder=4, clip_on=False)
    elif action == "form":
        ax.plot([x + 0.03, x + 0.19], [y0 - 0.01 * span, y0 - 0.03 * span], color="#111111", lw=0.65, ls="--", zorder=4, clip_on=False)


def _figure_eligibility(
    structure_relaxation: Dict,
    result: Dict,
) -> tuple[bool, List[str]]:
    blockers = []
    if structure_relaxation:
        if structure_relaxation.get("status") != "success":
            blockers.append("MACE relaxation was not completed")
        elif structure_relaxation.get("relaxation_status") not in {
            "converged",
            "converged_constrained",
        }:
            blockers.append(
                "MACE relaxation did not reach the configured force convergence threshold"
            )

    calculation_status = result.get("calculation_status")
    if calculation_status in {"reaction_scanned", "redox_state_scanned"}:
        values = (
            result.get("reaction_profile_relative_energies_ev")
            if calculation_status == "reaction_scanned"
            else result.get("redox_state_profile_relative_energies_ev")
        ) or []
        finite_count = sum(
            value is not None and math.isfinite(float(value))
            for value in values
        )
        if finite_count < 5:
            blockers.append("mechanism profile has fewer than 5 finite scan points")
        quality = result.get("scan_quality") or {}
        if quality and not quality.get("all_frames_converged", False):
            blockers.append("one or more constrained tblite scan frames did not converge")

    adsorption_optimization = result.get("adsorption_local_optimization") or {}
    if adsorption_optimization:
        if adsorption_optimization.get("status") != "success":
            blockers.append("tblite adsorption geometry optimization failed")
        elif adsorption_optimization.get("converged") is not True:
            blockers.append("tblite adsorption geometry optimization did not converge")
    return not blockers, blockers


def _volcano_eligible_rows(summary: Dict) -> List[Dict]:
    groups: Dict[str, List[Dict]] = {}
    for row in summary.get("activity_results", []):
        if not row.get("figure_eligible", True):
            continue
        if row.get("best_adsorption_energy_ev") is None or row.get("activation_metric_ev") is None:
            continue
        task_id = str(row.get("task_id") or "")
        if task_id:
            groups.setdefault(task_id, []).append(row)
    eligible = []
    for rows in groups.values():
        if len(rows) >= 5:
            eligible.extend(rows)
    return eligible


def _interpret_activity(payload: Dict, best_pose: Dict, reaction_profile: Dict) -> str:
    if payload.get("status") != "success":
        return "failed"
    if best_pose.get("adsorption_energy_ev") is None:
        return "geometry-only substrate pose; configure an energy backend for ranking"
    redox_profile = payload.get("redox_state_profile") or {}
    if redox_profile.get("status") == "success":
        return "adsorption and redox charge/spin activation profile calculated"
    if reaction_profile.get("status") == "success":
        return "adsorption and hydrolysis coordinate profile calculated"
    active_profile = profile_from_payload(payload)
    if active_profile.get("status") == "insufficient_sampling":
        return "adsorption calculated; mechanism profile was attempted but needs at least five scan points"
    if active_profile.get("status") == "incomplete":
        return "adsorption calculated; mechanism profile was attempted but one or more scan frames did not converge"
    return "adsorption calculated; mechanism profile requires tblite/GFN2-xTB and scan setting"


def _finite_or_none(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _mean(values: Iterable[float]) -> float | None:
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def _symmetric_xlim(ax, values: List[float], margin: float = 0.2) -> None:
    finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not finite:
        ax.set_xlim(-1.0, 1.0)
        return
    lo = min(finite + [0.0])
    hi = max(finite + [0.0])
    span = max(hi - lo, 1.0)
    ax.set_xlim(lo - margin * span, hi + margin * span)
