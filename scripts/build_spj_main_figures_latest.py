#!/usr/bin/env python3
"""Build five SPJ main figures from the canonical x1-x100 dataset."""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw

try:
    from _io_utils import pct_ratio as pct
    from build_x1_x100_dataset import ACTIVITY_ORDER, TOPOLOGY_ORDER
    from figure_style import (
        SPJ_WORKBENCH_PALETTE,
        VISUAL_ATLAS_ACTIVITY_COLORS,
        VISUAL_ATLAS_SOURCE_COLORS,
        apply_evidence_rcparams,
    )
except ModuleNotFoundError:
    from scripts._io_utils import pct_ratio as pct  # type: ignore
    from scripts.build_x1_x100_dataset import ACTIVITY_ORDER, TOPOLOGY_ORDER  # type: ignore
    from scripts.figure_style import (  # type: ignore
        SPJ_WORKBENCH_PALETTE,
        VISUAL_ATLAS_ACTIVITY_COLORS,
        VISUAL_ATLAS_SOURCE_COLORS,
        apply_evidence_rcparams,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "outputs" / "x1_x100_dataset"
DEFAULT_OUT_DIR = ROOT / "outputs" / "spj_main_figures_latest"


def _absolute(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def _env_path(name: str, default: Path) -> Path:
    return _absolute(os.environ.get(name, str(default)))


DATA_DIR = _env_path("E2N_SPJ_DATA_DIR", DEFAULT_DATA_DIR)
OUT_DIR = _env_path("E2N_SPJ_OUT_DIR", DEFAULT_OUT_DIR)
FIG_DIR = OUT_DIR / "figures"
SOURCE_DIR = OUT_DIR / "source_data"
QA_DIR = OUT_DIR / "qa"
MANUSCRIPT_DIR = OUT_DIR / "manuscript"


def _refresh_output_dirs() -> None:
    global FIG_DIR, SOURCE_DIR, QA_DIR, MANUSCRIPT_DIR
    FIG_DIR = OUT_DIR / "figures"
    SOURCE_DIR = OUT_DIR / "source_data"
    QA_DIR = OUT_DIR / "qa"
    MANUSCRIPT_DIR = OUT_DIR / "manuscript"


def configure_paths(
    *,
    data_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
) -> None:
    global DATA_DIR, OUT_DIR
    if data_dir is not None:
        DATA_DIR = _absolute(data_dir)
    if out_dir is not None:
        OUT_DIR = _absolute(out_dir)
        _refresh_output_dirs()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=str(DATA_DIR))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    return parser.parse_args(argv)


INK = SPJ_WORKBENCH_PALETTE["ink"]
MUTED = SPJ_WORKBENCH_PALETTE["muted"]
LINE = SPJ_WORKBENCH_PALETTE["line"]
SOFT = SPJ_WORKBENCH_PALETTE["soft"]
TEAL = SPJ_WORKBENCH_PALETTE["teal"]
BLUE = SPJ_WORKBENCH_PALETTE["blue"]
AMBER = SPJ_WORKBENCH_PALETTE["amber"]
GREEN = SPJ_WORKBENCH_PALETTE["green"]
GRAY = SPJ_WORKBENCH_PALETTE["gray"]

TOPOLOGY_COLORS = {
    "bridged": AMBER,
    "independent adjacent": BLUE,
    "independent separated": TEAL,
}
TOPOLOGY_MARKERS = {
    "bridged": "s",
    "independent adjacent": "^",
    "independent separated": "o",
}
ACTIVITY_SHORT = {
    "Catalase": "CAT",
    "Oxidase": "OXD",
    "Peroxidase": "POD",
    "Glutathione Peroxidase": "GPx",
    "Glucose Oxidase": "GOx",
    "DNase": "DNase",
}
METHOD_ORDER = [
    "First pass",
    "GFN1 extended",
    "GFN1 SCF fallback",
    "GFN2 deep",
    "GFN2 extended",
]
METHOD_KEY = {
    "First pass": "first_pass",
    "GFN1 extended": "gfn1_extended",
    "GFN1 SCF fallback": "gfn1_scf_fallback",
    "GFN2 deep": "gfn2_deep",
    "GFN2 extended": "gfn2_extended",
}
METHOD_COLORS = {
    label: VISUAL_ATLAS_SOURCE_COLORS[key]
    for label, key in METHOD_KEY.items()
}


def configure_style() -> None:
    apply_evidence_rcparams(
        mpl.rcParams,
        colors=SPJ_WORKBENCH_PALETTE,
        title_size=8,
        label_size=7,
        tick_size=6.5,
        legend_size=6.2,
        savefig_dpi=450,
        svg_fonttype_none=True,
    )
    mpl.rcParams.update(
        {
            "axes.linewidth": 0.7,
            "axes.titlepad": 6,
            "axes.unicode_minus": False,
            "figure.dpi": 150,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": False,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "lines.solid_capstyle": "round",
        }
    )


def reset_output() -> None:
    target = OUT_DIR.resolve()
    if target == ROOT.resolve() or target == Path(target.anchor):
        raise RuntimeError(f"Refusing to clear unsafe output path: {target}")
    if target.exists():
        shutil.rmtree(target)
    for path in [FIG_DIR, SOURCE_DIR, QA_DIR, MANUSCRIPT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def load_dataset() -> dict[str, pd.DataFrame]:
    required = {
        "candidates": DATA_DIR / "candidates.csv",
        "profiles": DATA_DIR / "profiles.csv",
        "designs": DATA_DIR / "designs.csv",
        "topology_tests": DATA_DIR / "topology_tests.csv",
        "representative": DATA_DIR / "representative.csv",
        "representative_scans": DATA_DIR / "representative_scans.csv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Canonical dataset is incomplete: " + "; ".join(missing))
    data = {name: pd.read_csv(path) for name, path in required.items()}
    candidates = data["candidates"]
    profiles = data["profiles"]
    tests = data["topology_tests"]
    if len(candidates) != 355 or candidates["candidate_id"].nunique() != 355:
        raise RuntimeError("Canonical candidate table failed the 355-row uniqueness check")
    if len(profiles) != 699 or profiles["profile_id"].nunique() != 699:
        raise RuntimeError("Canonical profile table failed the 699-row uniqueness check")
    if not profiles["frame_count"].eq(profiles["converged_frames"]).all():
        raise RuntimeError("Canonical profile table contains a non-converged frame")
    if len(tests) != 10 or int(tests["retained_at_fdr_0_05"].sum()) != 8:
        raise RuntimeError("Canonical topology-test table failed the 10-test / 8-retained check")
    return data


def clean_axis(ax: plt.Axes, *, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.set_axisbelow(True)


def rule_title(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=7.7, fontweight="bold", color=INK, pad=6)
    ax.plot([0, 1], [1.015, 1.015], transform=ax.transAxes, color=LINE, lw=0.7, clip_on=False)


def add_panel_labels(fig: plt.Figure, axes: list[plt.Axes], labels: str) -> None:
    for ax, label in zip(axes, labels):
        box = ax.get_position()
        fig.text(
            box.x0 - 0.026,
            box.y1 + 0.012,
            label,
            ha="left",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=INK,
        )


def save_figure(fig: plt.Figure, stem: str, *, height: float) -> dict[str, str]:
    fig.set_size_inches(7.0, height, forward=True)
    paths: dict[str, str] = {}
    for extension in ["png", "svg", "pdf"]:
        path = FIG_DIR / f"{stem}.{extension}"
        kwargs: dict[str, Any] = {"facecolor": "white", "edgecolor": "none"}
        if extension == "png":
            kwargs["dpi"] = 450
        fig.savefig(path, **kwargs)
        paths[extension] = str(path)
    plt.close(fig)
    with Image.open(paths["png"]) as image:
        if image.mode != "RGB":
            image.convert("RGB").save(paths["png"])
    return paths


def pair_short(value: str) -> str:
    return "+".join(ACTIVITY_SHORT.get(part.strip(), part.strip()) for part in value.split("+"))


def _jitter(values: np.ndarray, center: float, *, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(center, 0.055, size=len(values))


def draw_fig1(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    candidates = data["candidates"]
    profiles = data["profiles"]
    representative = data["representative"].copy()
    scans = data["representative_scans"].copy()

    design_topology = (
        candidates.groupby(["topology", "design_index"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    design_topology.to_csv(SOURCE_DIR / "fig1_design_topology.csv", index=False)

    representative.to_csv(SOURCE_DIR / "fig1_representative_candidate.csv", index=False)
    scans.to_csv(SOURCE_DIR / "fig1_representative_scans.csv", index=False)
    summary = pd.DataFrame(
        {
            "metric": ["retained_candidate_records", "complete_profiles", "converged_frames", "represented_design_indices"],
            "value": [
                len(candidates),
                len(profiles),
                int(profiles["converged_frames"].sum()),
                int(candidates["design_index"].nunique()),
            ],
        }
    )
    summary.to_csv(SOURCE_DIR / "fig1_canonical_summary.csv", index=False)

    contract = """# Fig. 1 Figure Contract

- **Claim:** The canonical x1-x100 release links a retained candidate record to a frozen structure snapshot, complete profile trajectories and corpus-wide design coverage.
- **Panel A:** A representative retained x57 record selected by image availability, two-profile linkage and proximity to the eligible-set median score.
- **Panel B:** The two complete activity-profile trajectories associated with that candidate.
- **Panel C:** Retained candidate-record counts at every design index and topology across the canonical release.
- **Decision rule:** Candidate rows satisfy source-specific retention rules; profile rows are complete with every recorded frame converged.
- **Boundary:** Empty design cells show absence from the retained public set, not a yield denominator.
- **Review risk:** The representative candidate illustrates traceability and is not presented as the best or experimentally validated catalyst.
"""
    (MANUSCRIPT_DIR / "fig1_figure_contract.md").write_text(contract, encoding="utf-8")

    fig = plt.figure(figsize=(7.0, 5.10))
    grid = fig.add_gridspec(
        2,
        12,
        left=0.08,
        right=0.93,
        bottom=0.095,
        top=0.93,
        height_ratios=[1.50, 0.68],
        hspace=0.42,
        wspace=0.48,
    )
    ax_a = fig.add_subplot(grid[0, :6])
    ax_b = fig.add_subplot(grid[0, 6:])
    ax_c = fig.add_subplot(grid[1, :])

    structure_paths = [
        DATA_DIR / "representative_structure.png",
        DATA_DIR / "assets" / "representative_structure.png",
    ]
    structure_path = next((path for path in structure_paths if path.exists()), None)
    if structure_path is None:
        checked = ", ".join(str(path) for path in structure_paths)
        raise FileNotFoundError(f"Representative structure image not found; checked: {checked}")
    shutil.copy2(structure_path, SOURCE_DIR / "fig1_representative_structure.png")
    with Image.open(structure_path) as structure_image:
        ax_a.imshow(structure_image.convert("RGB"))
    ax_a.set_axis_off()
    rep = representative.iloc[0]
    rule_title(ax_a, "Representative retained candidate")
    ax_a.text(
        0.02,
        0.03,
        (
            f"{rep.design_id} | {rep.topology}\n"
            f"d = {rep.distance_a:.1f} A | {rep.doping} | alpha = {int(rep.angle_deg)} degrees\n"
            f"{rep.activity_pair}"
        ),
        transform=ax_a.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.2,
        color=INK,
        linespacing=1.35,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 2.5},
    )

    scan_colors = {
        "Glucose Oxidase": VISUAL_ATLAS_ACTIVITY_COLORS["Glucose Oxidase"],
        "Peroxidase": VISUAL_ATLAS_ACTIVITY_COLORS["Peroxidase"],
    }
    for activity, rows in scans.groupby("activity", sort=False):
        rows = rows.sort_values("scan_step")
        color = scan_colors.get(activity, TEAL)
        ax_b.plot(
            rows["scan_step"],
            rows["relative_energy_ev"],
            color=color,
            lw=1.8,
            marker="o",
            markersize=4.8,
            markeredgecolor="white",
            markeredgewidth=0.6,
        )
        last = rows.iloc[-1]
        ax_b.text(
            float(last["scan_step"]) + 0.10,
            float(last["relative_energy_ev"])
            + {"Glucose Oxidase": 0.10, "Peroxidase": -0.10}.get(activity, 0.0),
            ACTIVITY_SHORT.get(activity, activity),
            va="center",
            ha="left",
            fontsize=6.4,
            color=color,
            fontweight="bold",
        )
    ax_b.axhline(0, color=LINE, lw=0.7, zorder=0)
    ax_b.set_xlim(0.8, 5.55)
    ax_b.set_xticks([1, 2, 3, 4, 5])
    ax_b.set_xlabel("Scan step")
    ax_b.set_ylabel("Relative energy (eV)")
    rule_title(ax_b, "Complete profile trajectories")
    clean_axis(ax_b)

    matrix = np.zeros((len(TOPOLOGY_ORDER), 100), dtype=float)
    for row in design_topology.itertuples(index=False):
        matrix[TOPOLOGY_ORDER.index(row.topology), int(row.design_index) - 1] = int(row.candidate_count)
    count_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "candidate_counts",
        ["#FFFFFF", "#DCEAE7", "#73B0A8", TEAL],
    )
    image = ax_c.imshow(matrix, aspect="auto", interpolation="nearest", cmap=count_cmap, vmin=0, vmax=max(1, matrix.max()))
    ax_c.set_yticks(range(3), ["bridged", "adjacent", "separated"])
    tick_values = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    ax_c.set_xticks([value - 1 for value in tick_values], tick_values)
    ax_c.set_xlabel("Design index")
    ax_c.axvline(int(rep.design_index) - 1, color=AMBER, lw=1.0)
    ax_c.text(
        int(rep.design_index) - 1,
        2.46,
        rep.design_id,
        ha="center",
        va="bottom",
        fontsize=6.2,
        color=AMBER,
        fontweight="bold",
    )
    colorbar = fig.colorbar(image, ax=ax_c, fraction=0.018, pad=0.012)
    colorbar.set_label("Candidate records", fontsize=6.2)
    rule_title(
        ax_c,
        (
            f"Canonical release | {len(candidates)} retained records | "
            f"{len(profiles)} complete profiles | {int(profiles['converged_frames'].sum())} frames"
        ),
    )

    add_panel_labels(fig, [ax_a, ax_b, ax_c], "ABC")
    return save_figure(fig, "fig1_canonical_evidence_trace", height=5.10)


def draw_fig2(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    candidates = data["candidates"]
    geometry = (
        candidates.groupby(["design_index", "topology", "distance_a", "doping", "angle_deg"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    geometry.to_csv(SOURCE_DIR / "fig2_geometry.csv", index=False)

    highlight_mask = (
        candidates["topology"].eq("independent separated")
        & candidates["distance_a"].eq(13.0)
        & candidates["doping"].eq("NS")
        & candidates["angle_deg"].isin([75, 105])
    )
    highlighted = (
        candidates[highlight_mask]
        .groupby(["design_index", "angle_deg"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    highlighted.to_csv(SOURCE_DIR / "fig2_highlighted_geometry.csv", index=False)
    candidates[["candidate_id", "topology", "distance_a", "doping", "angle_deg"]].to_csv(
        SOURCE_DIR / "fig2_distance_points.csv", index=False
    )

    contract = """# Fig. 2 Figure Contract

- **Claim:** Retained candidate records occupy discrete requested-geometry regimes, including a dense independent-separated d=13.0 A, NS, alpha=75/105 region.
- **Panel A:** Bubble area is retained candidate-record count at each distance-angle coordinate; shape is topology and color is doping.
- **Panel B:** The highlighted region is resolved into its contributing design indices.
- **Panel C:** Individual candidate distances are displayed by topology with median and IQR.
- **Decision rule:** The highlighted region requires independent-separated topology, d=13.0 A, NS and alpha=75 or 105 degrees.
- **Boundary:** Counts within the retained set are not yield estimates and do not establish catalytic or topology superiority.
- **Review risk:** Repeated coordinates can hide multiplicity; bubble area and direct counts make multiplicity explicit.
"""
    (MANUSCRIPT_DIR / "fig2_figure_contract.md").write_text(contract, encoding="utf-8")

    fig = plt.figure(figsize=(7.0, 4.75))
    grid = fig.add_gridspec(
        2,
        3,
        left=0.09,
        right=0.96,
        bottom=0.10,
        top=0.91,
        width_ratios=[1.12, 1.12, 0.92],
        height_ratios=[0.88, 1.12],
        wspace=0.60,
        hspace=0.58,
    )
    ax_a = fig.add_subplot(grid[:, :2])
    ax_b = fig.add_subplot(grid[0, 2])
    ax_c = fig.add_subplot(grid[1, 2])

    doping_colors = {"N": BLUE, "NS": TEAL, "S": AMBER}
    for topology in TOPOLOGY_ORDER:
        for doping, color in doping_colors.items():
            rows = geometry[geometry["topology"].eq(topology) & geometry["doping"].eq(doping)]
            if rows.empty:
                continue
            ax_a.scatter(
                rows["distance_a"],
                rows["angle_deg"],
                s=18 + rows["candidate_count"] * 3.0,
                marker=TOPOLOGY_MARKERS[topology],
                color=color,
                edgecolor=INK,
                linewidth=0.45,
                alpha=0.82,
            )
    ax_a.add_patch(
        Rectangle((12.72, 66), 0.56, 48, facecolor="#E7F3F1", edgecolor=TEAL, linewidth=1.1, linestyle="--")
    )
    selected_geometry = geometry[
        geometry["topology"].eq("independent separated")
        & geometry["distance_a"].eq(13.0)
        & geometry["doping"].eq("NS")
        & geometry["angle_deg"].isin([75, 105])
    ]
    ax_a.scatter(
        selected_geometry["distance_a"],
        selected_geometry["angle_deg"],
        s=28 + selected_geometry["candidate_count"] * 3.2,
        marker="o",
        facecolor="none",
        edgecolor=INK,
        linewidth=1.2,
        zorder=4,
    )
    for angle in [75, 105]:
        count = int(selected_geometry.loc[selected_geometry["angle_deg"].eq(angle), "candidate_count"].sum())
        ax_a.annotate(
            f"{count}",
            xy=(13.0, angle),
            xytext=(-36, -10 if angle == 75 else 10),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=6.3,
            color=INK,
            arrowprops=dict(arrowstyle="-", color=TEAL, lw=0.8),
        )
    ax_a.set_xlabel("Metal-pair distance (A)")
    ax_a.set_ylabel("Angle (degrees)")
    ax_a.set_xlim(5.1, 15.1)
    ax_a.set_ylim(-8, 160)
    rule_title(ax_a, "Retained-record requested geometry")
    clean_axis(ax_a, grid_axis="both")
    topology_legend = ax_a.legend(
        handles=[
            Line2D([0], [0], marker=TOPOLOGY_MARKERS[value], color="none", markerfacecolor="#D5E5E2", markeredgecolor=INK, label=value, markersize=6)
            for value in TOPOLOGY_ORDER
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.985),
        ncols=3,
        borderaxespad=0.2,
        columnspacing=0.8,
        handletextpad=0.3,
    )
    ax_a.add_artist(topology_legend)
    ax_a.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor=INK, label=doping, markersize=5.5)
            for doping, color in doping_colors.items()
        ],
        title="Doping",
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.985),
        ncols=3,
        borderaxespad=0.2,
        columnspacing=0.6,
        handletextpad=0.25,
        title_fontsize=6.0,
    )
    ax_a.text(0.01, 0.02, "Bubble area encodes candidate-record count", transform=ax_a.transAxes, fontsize=6.0, color=MUTED)

    rows = highlighted.sort_values("design_index")
    bars = ax_b.bar(np.arange(len(rows)), rows["candidate_count"], color=TEAL, width=0.68)
    ax_b.set_xticks(np.arange(len(rows)), [f"x{int(value)}" for value in rows["design_index"]])
    ax_b.set_ylabel("Candidate records")
    ax_b.set_ylim(0, rows["candidate_count"].max() * 1.30)
    for bar, value in zip(bars, rows["candidate_count"]):
        ax_b.text(bar.get_x() + bar.get_width() / 2, value + 1, str(int(value)), ha="center", fontsize=6.0)
    ax_b.text(
        0.02,
        0.98,
        "separated | d=13.0 A\nNS | alpha=75/105",
        transform=ax_b.transAxes,
        va="top",
        fontsize=5.8,
        color=MUTED,
        linespacing=1.25,
    )
    rule_title(ax_b, "Highlighted geometry")
    clean_axis(ax_b, grid_axis="y")

    y = np.arange(len(TOPOLOGY_ORDER))[::-1]
    for index, topology in enumerate(TOPOLOGY_ORDER):
        values = candidates.loc[candidates["topology"].eq(topology), "distance_a"].to_numpy(dtype=float)
        yi = y[index]
        ax_c.scatter(values, _jitter(values, yi, seed=100 + index), s=8, color=TOPOLOGY_COLORS[topology], alpha=0.25, edgecolor="none")
        q25, median, q75 = np.percentile(values, [25, 50, 75])
        ax_c.plot([q25, q75], [yi, yi], color=INK, lw=2.0)
        ax_c.scatter([median], [yi], s=28, color="white", edgecolor=INK, linewidth=0.8, zorder=4)
    ax_c.set_yticks(y, ["bridged", "adjacent", "separated"])
    ax_c.set_xlabel("Distance (A)")
    ax_c.set_ylim(-0.55, 2.55)
    rule_title(ax_c, "Distance distributions")
    clean_axis(ax_c, grid_axis="x")

    add_panel_labels(fig, [ax_a, ax_b, ax_c], "ABC")
    return save_figure(fig, "fig2_geometry_landscape", height=4.75)


def draw_fig3(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    candidates = data["candidates"]
    pair_topology = (
        candidates.groupby(["activity_pair", "topology"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    pair_topology.to_csv(SOURCE_DIR / "fig3_activity_pair_topology.csv", index=False)
    pair_counts = pair_topology.groupby("activity_pair", as_index=False)["candidate_count"].sum()
    pair_counts = pair_counts.sort_values("candidate_count", ascending=True).reset_index(drop=True)
    pair_counts.to_csv(SOURCE_DIR / "fig3_activity_pair_counts.csv", index=False)
    pair_topology_percent = pair_topology.merge(
        pair_counts.rename(columns={"candidate_count": "pair_candidate_count"}),
        on="activity_pair",
        how="left",
    )
    pair_topology_percent["topology_percent"] = (
        100.0 * pair_topology_percent["candidate_count"] / pair_topology_percent["pair_candidate_count"]
    )
    pair_topology_percent.to_csv(SOURCE_DIR / "fig3_activity_pair_topology_percent.csv", index=False)
    pair_design = (
        candidates.groupby(["activity_pair", "design_index"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    pair_design.to_csv(SOURCE_DIR / "fig3_activity_pair_design.csv", index=False)

    contract = """# Fig. 3 Figure Contract

- **Claim:** Thirteen activity pairs are represented across the same x1-x100 design grid, with pair-specific topology composition.
- **Panel A:** Retained candidate-record counts for every activity pair.
- **Panel B:** Within-pair percentages for bridged, independent adjacent and independent separated topology.
- **Panel C:** Candidate occupancy by activity pair and design index.
- **Decision rule:** Each candidate contributes once to its declared activity pair and topology.
- **Boundary:** Unequal counts reflect the displayed retained set and are not normalized yields.
- **Review risk:** Long activity names can obscure comparison; stable activity abbreviations are used throughout.
"""
    (MANUSCRIPT_DIR / "fig3_figure_contract.md").write_text(contract, encoding="utf-8")

    pair_order = pair_counts["activity_pair"].tolist()
    labels = [pair_short(value) for value in pair_order]
    fig = plt.figure(figsize=(7.0, 5.25))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.13,
        right=0.93,
        bottom=0.09,
        top=0.93,
        height_ratios=[1.12, 0.88],
        width_ratios=[0.78, 1.22],
        hspace=0.54,
        wspace=0.36,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    y = np.arange(len(pair_order))
    bars = ax_a.barh(y, pair_counts["candidate_count"], color=TEAL, height=0.62)
    ax_a.set_yticks(y, labels)
    ax_a.set_xlabel("Candidate records")
    ax_a.set_xlim(0, pair_counts["candidate_count"].max() * 1.18)
    for bar, value in zip(bars, pair_counts["candidate_count"]):
        ax_a.text(value + 1.2, bar.get_y() + bar.get_height() / 2, str(int(value)), va="center", fontsize=5.9)
    rule_title(ax_a, "Candidate coverage by activity pair")
    clean_axis(ax_a, grid_axis="x")

    pivot = pd.pivot_table(
        pair_topology,
        index="activity_pair",
        columns="topology",
        values="candidate_count",
        aggfunc="sum",
        fill_value=0,
    ).reindex(pair_order, fill_value=0)
    left = np.zeros(len(pair_order), dtype=float)
    pair_totals = pivot.sum(axis=1).to_numpy(dtype=float)
    topology_hatches = {"bridged": "//", "independent adjacent": "..", "independent separated": ""}
    for topology in TOPOLOGY_ORDER:
        counts = pivot.get(topology, pd.Series(0, index=pivot.index)).to_numpy(dtype=float)
        values = np.divide(100.0 * counts, pair_totals, out=np.zeros_like(counts), where=pair_totals > 0)
        ax_b.barh(
            y,
            values,
            left=left,
            height=0.62,
            color=TOPOLOGY_COLORS[topology],
            edgecolor="white",
            linewidth=0.35,
            hatch=topology_hatches[topology],
            label=topology,
        )
        left += values
    ax_b.set_yticks(y, [])
    ax_b.set_ylim(-0.6, len(pair_order) + 0.8)
    ax_b.set_xlim(0, 100)
    ax_b.set_xlabel("Candidate-record share within pair (%)")
    ax_b.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.99), ncols=3, borderaxespad=0, columnspacing=0.65)
    rule_title(ax_b, "Topology composition of each pair")
    clean_axis(ax_b, grid_axis="x")

    matrix = np.zeros((len(pair_order), 100), dtype=float)
    for row in pair_design.itertuples(index=False):
        matrix[pair_order.index(row.activity_pair), int(row.design_index) - 1] = int(row.candidate_count)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("pair_design", ["#FFFFFF", "#DDEAE8", "#75AEA7", TEAL])
    image = ax_c.imshow(matrix, aspect="auto", interpolation="nearest", origin="lower", cmap=cmap, vmin=0, vmax=max(1, matrix.max()))
    ax_c.set_yticks(range(len(pair_order)), labels)
    tick_values = [1, 20, 40, 60, 80, 100]
    ax_c.set_xticks([value - 1 for value in tick_values], tick_values)
    ax_c.set_xlabel("Design index")
    colorbar = fig.colorbar(image, ax=ax_c, fraction=0.018, pad=0.012)
    colorbar.set_label("Candidate records", fontsize=6.2)
    rule_title(ax_c, "Activity-pair occupancy across x1-x100")

    add_panel_labels(fig, [ax_a, ax_b, ax_c], "ABC")
    return save_figure(fig, "fig3_activity_pair_composition", height=5.25)


def _descriptor_panel(
    ax: plt.Axes,
    profiles: pd.DataFrame,
    column: str,
    xlabel: str,
    *,
    show_zero_fraction: bool,
) -> None:
    y = np.arange(len(ACTIVITY_ORDER))[::-1]
    x_max = float(profiles[column].max())
    x_min = float(profiles[column].min())
    for index, activity in enumerate(ACTIVITY_ORDER):
        values = profiles.loc[profiles["activity"].eq(activity), column].dropna().to_numpy(dtype=float)
        if len(values) == 0:
            continue
        yi = y[index]
        ax.scatter(
            values,
            _jitter(values, yi, seed=400 + index + (100 if show_zero_fraction else 0)),
            s=7,
            color=VISUAL_ATLAS_ACTIVITY_COLORS[activity],
            alpha=0.24,
            edgecolor="none",
        )
        if len(values) >= 10:
            q25, median, q75 = np.percentile(values, [25, 50, 75])
            ax.plot([q25, q75], [yi, yi], color=INK, lw=2.0, zorder=4)
            ax.scatter([median], [yi], s=26, color="white", edgecolor=INK, linewidth=0.8, zorder=5)
        else:
            ax.scatter(values, np.full(len(values), yi), s=34, marker="D", color="white", edgecolor=INK, linewidth=0.8, zorder=5)
        if show_zero_fraction:
            zero_count = int(np.sum(values == 0))
            ax.text(
                x_max * 1.04,
                yi,
                f"peakless {zero_count}/{len(values)}",
                ha="left",
                va="center",
                fontsize=5.7,
                color=MUTED,
            )
    ax.axvline(0, color=LINE, lw=0.7)
    ax.set_yticks(y, [ACTIVITY_SHORT[value] for value in ACTIVITY_ORDER])
    ax.set_ylim(-0.55, len(ACTIVITY_ORDER) - 0.45)
    if show_zero_fraction:
        ax.set_xlim(min(-0.15, x_min - 0.1), x_max * 1.28)
    else:
        ax.set_xlim(x_min - 0.25, x_max + 0.35)
    ax.set_xlabel(xlabel)
    clean_axis(ax, grid_axis="x")


def draw_fig4(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    profiles = data["profiles"]
    activity_counts = (
        profiles.groupby("activity", as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    activity_counts.to_csv(SOURCE_DIR / "fig4_activity_counts.csv", index=False)
    method_activity = (
        profiles.groupby(["activity", "selected_method"], as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    method_activity.to_csv(SOURCE_DIR / "fig4_method_activity.csv", index=False)
    descriptor_points = profiles[
        ["profile_id", "activity", "topology", "adsorption_energy_ev", "activation_metric_ev"]
    ].copy()
    descriptor_points.to_csv(SOURCE_DIR / "fig4_descriptor_points.csv", index=False)
    zero_activation = profiles.groupby("activity", as_index=False).agg(
        profile_count=("profile_id", "size"),
        zero_activation_count=("activation_metric_ev", lambda values: int(np.sum(np.asarray(values) == 0))),
    )
    zero_activation["zero_activation_fraction"] = zero_activation["zero_activation_count"] / zero_activation["profile_count"]
    zero_activation.to_csv(SOURCE_DIR / "fig4_zero_activation.csv", index=False)

    contract = """# Fig. 4 Figure Contract

- **Claim:** The 699 complete profiles have activity-specific method provenance and descriptor distributions.
- **Panel A:** Complete-profile count by activity.
- **Panel B:** Method composition within each activity.
- **Panel C:** Every adsorption descriptor, with median and IQR where n>=10.
- **Panel D:** Every activation descriptor, with zero-valued boundary counts and median/IQR where n>=10.
- **Decision rule:** Only complete profiles with all recorded frames converged are included.
- **Boundary:** Zero-valued activation metrics are peakless descriptor cases, not zero barriers; descriptors do not establish catalytic activity.
- **Review risk:** DNase has n=1 and is shown as one point without distributional inference.
"""
    (MANUSCRIPT_DIR / "fig4_figure_contract.md").write_text(contract, encoding="utf-8")

    fig = plt.figure(figsize=(7.0, 4.95))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.09,
        right=0.96,
        bottom=0.10,
        top=0.92,
        height_ratios=[0.86, 1.14],
        width_ratios=[0.92, 1.08],
        hspace=0.58,
        wspace=0.45,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    counts = activity_counts.set_index("activity").reindex(ACTIVITY_ORDER).reset_index().iloc[::-1]
    y = np.arange(len(counts))
    bars = ax_a.barh(y, counts["profile_count"], color=[VISUAL_ATLAS_ACTIVITY_COLORS[value] for value in counts["activity"]], height=0.58)
    ax_a.set_yticks(y, [ACTIVITY_SHORT[value] for value in counts["activity"]])
    ax_a.set_xlabel("Complete profiles")
    ax_a.set_xlim(0, counts["profile_count"].max() * 1.17)
    for bar, value in zip(bars, counts["profile_count"]):
        ax_a.text(value + 3, bar.get_y() + bar.get_height() / 2, str(int(value)), va="center", fontsize=6.1)
    rule_title(ax_a, "Profile coverage")
    clean_axis(ax_a, grid_axis="x")

    method_matrix = pd.pivot_table(
        method_activity,
        index="activity",
        columns="selected_method",
        values="profile_count",
        aggfunc="sum",
        fill_value=0,
    ).reindex(ACTIVITY_ORDER, fill_value=0)
    percentages = method_matrix.div(method_matrix.sum(axis=1), axis=0) * 100.0
    y = np.arange(len(ACTIVITY_ORDER))[::-1]
    left = np.zeros(len(ACTIVITY_ORDER), dtype=float)
    method_hatches = {
        "First pass": "",
        "GFN1 extended": "//",
        "GFN1 SCF fallback": "..",
        "GFN2 deep": "xx",
        "GFN2 extended": "\\\\",
    }
    for method in METHOD_ORDER:
        values = percentages.get(method, pd.Series(0, index=percentages.index)).to_numpy(dtype=float)
        ax_b.barh(
            y,
            values,
            left=left,
            height=0.58,
            color=METHOD_COLORS[method],
            edgecolor="white",
            linewidth=0.35,
            hatch=method_hatches[method],
            label=method,
        )
        left += values
    ax_b.set_yticks(y, [ACTIVITY_SHORT[value] for value in ACTIVITY_ORDER])
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(-0.5, len(ACTIVITY_ORDER) + 0.8)
    ax_b.set_xlabel("Profiles within activity (%)")
    ax_b.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.99), ncols=3, borderaxespad=0, columnspacing=0.55, fontsize=5.4)
    rule_title(ax_b, "Selected calculation method")
    clean_axis(ax_b, grid_axis="x")

    _descriptor_panel(ax_c, profiles, "adsorption_energy_ev", "Adsorption descriptor (eV)", show_zero_fraction=False)
    rule_title(ax_c, "Adsorption distributions")
    _descriptor_panel(ax_d, profiles, "activation_metric_ev", "Activation metric (eV)", show_zero_fraction=True)
    rule_title(ax_d, "Activation distributions | peakless boundary")

    add_panel_labels(fig, [ax_a, ax_b, ax_c, ax_d], "ABCD")
    return save_figure(fig, "fig4_profile_descriptors", height=4.95)


def draw_fig5(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    profiles = data["profiles"]
    tests = data["topology_tests"].copy()
    tests.to_csv(SOURCE_DIR / "fig5_topology_tests.csv", index=False)
    metric_column = {
        "adsorption": "adsorption_energy_ev",
        "activation": "activation_metric_ev",
    }
    median_rows: list[dict[str, Any]] = []
    retained_tests = tests[tests["retained_at_fdr_0_05"]].copy()
    for test in retained_tests.itertuples(index=False):
        for topology in TOPOLOGY_ORDER:
            values = profiles.loc[
                profiles["activity"].eq(test.activity) & profiles["topology"].eq(topology),
                metric_column[test.metric],
            ].dropna()
            if len(values) >= 2:
                median_rows.append(
                    {
                        "activity": test.activity,
                        "metric": test.metric,
                        "topology": topology,
                        "n": int(len(values)),
                        "median_ev": float(values.median()),
                    }
                )
    medians = pd.DataFrame(median_rows)
    medians.to_csv(SOURCE_DIR / "fig5_topology_medians.csv", index=False)

    contract = """# Fig. 5 Figure Contract

- **Claim:** Topology is associated with descriptor distributions for 8 of 10 valid activity-descriptor tests.
- **Panel A:** Exact Benjamini-Hochberg-adjusted q values.
- **Panel B:** Epsilon-squared effect sizes, independent of q value.
- **Panel C:** Topology medians for the eight associations retained at q<0.05.
- **Decision rule:** Kruskal-Wallis requires at least two topology groups with n>=2; BH correction is applied across all ten valid tests.
- **Boundary:** Statistical association does not establish causal topology effects, catalytic superiority or experimental activity.
- **Review risk:** Significance can obscure magnitude; the effect matrix and median directions are shown separately.
"""
    (MANUSCRIPT_DIR / "fig5_figure_contract.md").write_text(contract, encoding="utf-8")

    tested_activities = [value for value in ACTIVITY_ORDER if value != "DNase"]
    metric_order = ["adsorption", "activation"]
    q_matrix = np.full((len(tested_activities), 2), np.nan)
    effect_matrix = np.full_like(q_matrix, np.nan)
    retained_matrix = np.zeros_like(q_matrix, dtype=bool)
    for row in tests.itertuples(index=False):
        i = tested_activities.index(row.activity)
        j = metric_order.index(row.metric)
        q_matrix[i, j] = row.q_value
        effect_matrix[i, j] = row.epsilon_squared
        retained_matrix[i, j] = bool(row.retained_at_fdr_0_05)

    fig = plt.figure(figsize=(7.0, 4.85))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.10,
        right=0.965,
        bottom=0.10,
        top=0.92,
        height_ratios=[0.92, 1.08],
        width_ratios=[1, 1],
        hspace=0.58,
        wspace=0.52,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    q_display = -np.log10(np.clip(q_matrix, 1e-8, 1.0))
    q_cmap = mpl.colors.LinearSegmentedColormap.from_list("q_values", ["#F2F4F3", "#A9D2CC", TEAL])
    image_a = ax_a.imshow(q_display, cmap=q_cmap, aspect="auto", vmin=0, vmax=max(2.0, float(np.nanmax(q_display))))
    ax_a.set_xticks([0, 1], ["Adsorption", "Activation"])
    ax_a.set_yticks(range(len(tested_activities)), [ACTIVITY_SHORT[value] for value in tested_activities])
    for i in range(q_matrix.shape[0]):
        for j in range(q_matrix.shape[1]):
            value = q_matrix[i, j]
            label = "<0.001" if value < 0.001 else f"{value:.3f}"
            ax_a.text(j, i, label, ha="center", va="center", fontsize=6.2, color="white" if q_display[i, j] > 1.5 else INK, fontweight="bold" if retained_matrix[i, j] else "normal")
            if retained_matrix[i, j]:
                ax_a.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False, edgecolor=INK, linewidth=0.9))
    rule_title(ax_a, "BH q values | 8/10 retained")
    cbar_a = fig.colorbar(image_a, ax=ax_a, fraction=0.045, pad=0.03)
    cbar_a.set_label("-log10(q)", fontsize=6.2)

    effect_cmap = mpl.colors.LinearSegmentedColormap.from_list("effects", ["#F4F5F4", "#E7C68B", AMBER])
    image_b = ax_b.imshow(effect_matrix, cmap=effect_cmap, aspect="auto", vmin=0, vmax=max(0.35, float(np.nanmax(effect_matrix))))
    ax_b.set_xticks([0, 1], ["Adsorption", "Activation"])
    ax_b.set_yticks(range(len(tested_activities)), [ACTIVITY_SHORT[value] for value in tested_activities])
    for i in range(effect_matrix.shape[0]):
        for j in range(effect_matrix.shape[1]):
            value = effect_matrix[i, j]
            ax_b.text(j, i, f"{value:.3f}", ha="center", va="center", fontsize=6.2, color="white" if value > 0.25 else INK, fontweight="bold" if retained_matrix[i, j] else "normal")
    rule_title(ax_b, "Epsilon-squared effect size")
    cbar_b = fig.colorbar(image_b, ax=ax_b, fraction=0.045, pad=0.03)
    cbar_b.set_label("epsilon-squared", fontsize=6.2)

    retained = tests[tests["retained_at_fdr_0_05"]].copy()
    retained["activity_order"] = retained["activity"].map({value: index for index, value in enumerate(ACTIVITY_ORDER)})
    retained["metric_order"] = retained["metric"].map({"adsorption": 0, "activation": 1})
    retained = retained.sort_values(["activity_order", "metric_order"]).reset_index(drop=True)
    y = np.arange(len(retained))[::-1]
    labels: list[str] = []
    for yi, test in zip(y, retained.itertuples(index=False)):
        labels.append(f"{ACTIVITY_SHORT[test.activity]} | {test.metric[:3]}")
        values: list[float] = []
        for topology in TOPOLOGY_ORDER:
            hit = medians[
                medians["activity"].eq(test.activity)
                & medians["metric"].eq(test.metric)
                & medians["topology"].eq(topology)
            ]
            if hit.empty:
                continue
            value = float(hit.iloc[0]["median_ev"])
            values.append(value)
            ax_c.scatter([value], [yi], s=32, marker=TOPOLOGY_MARKERS[topology], color=TOPOLOGY_COLORS[topology], edgecolor="white", linewidth=0.6, zorder=3)
        if len(values) >= 2:
            ax_c.plot([min(values), max(values)], [yi, yi], color=LINE, lw=1.4, zorder=1)
    ax_c.axvline(0, color=LINE, lw=0.7)
    ax_c.set_yticks(y, labels)
    ax_c.set_xlabel("Topology median descriptor value (eV)")
    ax_c.set_ylim(-0.65, len(retained) + 0.15)
    ax_c.legend(
        handles=[
            Line2D([0], [0], marker=TOPOLOGY_MARKERS[value], color="none", markerfacecolor=TOPOLOGY_COLORS[value], markeredgecolor="white", label=value, markersize=6)
            for value in TOPOLOGY_ORDER
        ],
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.985),
        ncols=3,
        borderaxespad=0.2,
    )
    rule_title(ax_c, "Topology medians for retained activity-descriptor pairs")
    clean_axis(ax_c, grid_axis="x")

    add_panel_labels(fig, [ax_a, ax_b, ax_c], "ABC")
    return save_figure(fig, "fig5_topology_statistics", height=4.85)


def write_manuscript_support(data: dict[str, pd.DataFrame]) -> dict[str, str]:
    candidates = data["candidates"]
    profiles = data["profiles"]
    panel_rows = [
        ("Fig. 1", "A", "Representative retained candidate record and its canonical metadata.", "fig1_representative_candidate.csv; fig1_representative_structure.png"),
        ("Fig. 1", "B", "Complete activity-profile trajectories for the representative candidate.", "fig1_representative_scans.csv"),
        ("Fig. 1", "C", "Candidate-record occupancy and canonical release totals across x1-x100.", "fig1_design_topology.csv; fig1_canonical_summary.csv"),
        ("Fig. 2", "A", "Retained candidate-record counts at each requested geometry coordinate.", "fig2_geometry.csv"),
        ("Fig. 2", "B", "Design-index composition of the highlighted geometry.", "fig2_highlighted_geometry.csv"),
        ("Fig. 2", "C", "Individual candidate distances by topology.", "fig2_distance_points.csv"),
        ("Fig. 3", "A", "Retained candidate-record counts by declared activity pair.", "fig3_activity_pair_counts.csv"),
        ("Fig. 3", "B", "Within-pair topology composition.", "fig3_activity_pair_topology_percent.csv"),
        ("Fig. 3", "C", "Activity-pair occupancy by design index.", "fig3_activity_pair_design.csv"),
        ("Fig. 4", "A", "Complete-profile count by activity.", "fig4_activity_counts.csv"),
        ("Fig. 4", "B", "Selected-method composition by activity.", "fig4_method_activity.csv"),
        ("Fig. 4", "C", "Individual adsorption descriptor values.", "fig4_descriptor_points.csv"),
        ("Fig. 4", "D", "Individual activation values and zero-valued boundary counts.", "fig4_descriptor_points.csv; fig4_zero_activation.csv"),
        ("Fig. 5", "A", "BH-adjusted q values for topology tests.", "fig5_topology_tests.csv"),
        ("Fig. 5", "B", "Epsilon-squared effect sizes for topology tests.", "fig5_topology_tests.csv"),
        ("Fig. 5", "C", "Topology medians for retained associations.", "fig5_topology_medians.csv"),
    ]
    panel_index = pd.DataFrame(panel_rows, columns=["figure", "panel", "claim", "source_data"])
    panel_index["source_data_paths"] = panel_index["source_data"].map(
        lambda value: "; ".join(f"source_data/{item.strip()}" for item in value.split(";"))
    )
    panel_index_path = MANUSCRIPT_DIR / "panel_source_data_index.csv"
    panel_index.to_csv(panel_index_path, index=False)

    captions = f"""# Main Figure Captions

Fig. 1. Traceability within the retained canonical x1-x100 set. (A) A representative x57 independent-separated hypothesis is shown at a requested distance of 13.0 A, NS doping, and a requested angle of 105 degrees; the structure is a constraint-compatible construction rather than a ranked or experimentally validated catalyst. (B) The two complete activity-profile trajectories associated with this candidate are shown as relative energies across their recorded scan steps. (C) Candidate-record occupancy is mapped over one uninterrupted x1-x100 axis for bridged, independent adjacent and independent separated topologies. The canonical set contains {len(candidates)} retained canonical candidate records, {len(profiles)} complete activity-specific profiles and {int(profiles['converged_frames'].sum())} converged scan frames within complete profiles. Only retained candidate records and complete-profile trajectories are displayed; empty design cells are not yield denominators, and the trajectories are computational descriptors rather than experimental rates or density-functional barriers.

Fig. 2. Requested geometry-grid occupancy of retained candidate records. (A) Marker shape denotes topology, color denotes N, NS or S doping and bubble area denotes candidate-record count at each requested distance-angle coordinate. The outlined region contains independent-separated candidate records at a requested distance of 13.0 A, NS doping, and requested angles of 75 or 105 degrees. (B) x17, x41, x57, x81 and x97 contribute records to the highlighted region. (C) Every requested candidate distance is shown by topology; horizontal segments denote the interquartile range and open points denote the median. Requested-grid counts within the retained set do not establish yield, catalytic superiority or unconstrained topology retention.

Fig. 3. Declared activity-pair composition in the canonical candidate-record table. (A) Counts are shown for all {candidates['activity_pair'].nunique()} activity pairs. (B) Each pair is decomposed into the percentage contributed by bridged, independent adjacent and independent separated topologies; absolute pair size remains in panel A. (C) Activity-pair occupancy is mapped across the common x1-x100 design axis. Each retained candidate record contributes once to its declared pair and topology; unequal totals describe the retained set and are not normalized yields.

Fig. 4. Computational-route provenance and descriptor distributions for complete profiles. (A) Complete-profile counts are shown by activity. (B) Stacked bars report the selected calculation method as a percentage within each activity. (C) Every adsorption descriptor is shown; horizontal segments and open points denote the interquartile range and median when n>=10. (D) Every forward scan peak descriptor is shown with zero-valued counts. Zero-valued forward scan peak descriptors are peakless boundary cases, not zero barriers. DNase is represented by one complete profile and is shown without distributional inference.

Fig. 5. Observational topology-labeled comparisons in the complete-profile table. (A) Exact Benjamini-Hochberg-adjusted q values are shown for ten valid Kruskal-Wallis activity-descriptor tests; outlined cells are the eight associations retained at q<0.05. (B) Epsilon-squared values report effect magnitude independently of q value. (C) Topology medians show direction for retained associations. Tests require at least two topology groups with n>=2 and correction is applied across all ten valid tests. These success-conditioned, protocol-mixed, confounded omnibus associations do not establish causal topology effects, catalytic superiority or experimental activity.
"""
    captions_path = MANUSCRIPT_DIR / "figure_captions_latest.md"
    captions_path.write_text(captions, encoding="utf-8")
    return {
        "captions": "manuscript/figure_captions_latest.md",
        "panel_source_data_index": "manuscript/panel_source_data_index.csv",
    }


def qa_png(path: str | Path) -> dict[str, Any]:
    raw = Image.open(path)
    image = raw.convert("RGB")
    array = np.asarray(image)
    nonwhite = np.any(array < 245, axis=2)
    sample = array[
        :: max(1, array.shape[0] // 160),
        :: max(1, array.shape[1] // 160),
    ].reshape(-1, 3)
    return {
        "mode": raw.mode,
        "width": image.width,
        "height": image.height,
        "nonwhite_pixels": int(nonwhite.sum()),
        "nonwhite_fraction": float(nonwhite.mean()),
        "unique_color_sample": int(len(np.unique(sample, axis=0))),
        "nonblank": bool(nonwhite.sum() > 2000 and len(np.unique(sample, axis=0)) > 8),
    }


def build_contact_sheet(paths: list[Path], output_path: Path, *, columns: int = 2) -> str:
    width = 920
    label_height = 42
    padding = 18
    cells: list[Image.Image] = []
    for path in paths:
        image = Image.open(path).convert("RGB")
        scale = width / image.width
        image = image.resize((width, int(image.height * scale)), Image.Resampling.LANCZOS)
        cell = Image.new("RGB", (width, image.height + label_height), "white")
        ImageDraw.Draw(cell).text((10, 10), path.name, fill=(23, 33, 29))
        cell.paste(image, (0, label_height))
        cells.append(cell)
    rows = math.ceil(len(cells) / columns)
    cell_width = width + 2 * padding
    cell_height = max(cell.height for cell in cells) + 2 * padding
    sheet = Image.new("RGB", (columns * cell_width, rows * cell_height), (242, 244, 242))
    for index, cell in enumerate(cells):
        sheet.paste(cell, ((index % columns) * cell_width + padding, (index // columns) * cell_height + padding))
    sheet.save(output_path, quality=95)
    return str(output_path)


def render_pdf(pdf_path: Path, output_path: Path) -> tuple[bool, str]:
    try:
        import pypdfium2 as pdfium

        document = pdfium.PdfDocument(pdf_path)
        if len(document) < 1:
            return False, "no pages"
        document[0].render(scale=2.0).to_pil().convert("RGB").save(output_path)
        return True, "pypdfium2"
    except Exception as error:
        renderer = shutil.which("pdftoppm")
        if not renderer:
            return False, str(error)
        prefix = output_path.with_suffix("")
        process = subprocess.run(
            [renderer, "-f", "1", "-singlefile", "-png", "-r", "160", str(pdf_path), str(prefix)],
            capture_output=True,
            text=True,
            check=False,
        )
        rendered = prefix.with_suffix(".png")
        if process.returncode != 0 or not rendered.exists():
            return False, process.stderr.strip() or "PDF render failed"
        if rendered != output_path:
            rendered.replace(output_path)
        return True, "pdftoppm"


def write_qa(outputs: dict[str, dict[str, str]]) -> dict[str, str]:
    rows: list[dict[str, Any]] = []
    grayscale_paths: list[Path] = []
    pdf_paths: list[Path] = []
    for figure, paths in outputs.items():
        png_path = Path(paths["png"])
        png_qa = qa_png(png_path)
        svg_text = Path(paths["svg"]).read_text(encoding="utf-8", errors="ignore")
        grayscale_path = QA_DIR / f"{png_path.stem}_grayscale.png"
        Image.open(png_path).convert("L").convert("RGB").save(grayscale_path)
        grayscale_paths.append(grayscale_path)
        pdf_path = Path(paths["pdf"])
        pdf_render = QA_DIR / f"{pdf_path.stem}_pdf_readback.png"
        rendered, renderer = render_pdf(pdf_path, pdf_render)
        pdf_qa = qa_png(pdf_render) if rendered else {"nonblank": False}
        if rendered:
            pdf_paths.append(pdf_render)
        row = {
            "figure": figure,
            **png_qa,
            "estimated_dpi_at_7in": round(png_qa["width"] / 7.0, 1),
            "meets_300dpi_at_7in": bool(png_qa["width"] >= 2100),
            "rgb": bool(png_qa["mode"] == "RGB"),
            "svg_nonempty": bool(Path(paths["svg"]).stat().st_size > 1000),
            "svg_text_editable": bool("<text" in svg_text),
            "pdf_nonempty": bool(pdf_path.stat().st_size > 1000),
            "pdf_readback": bool(rendered and pdf_qa["nonblank"]),
            "pdf_renderer": renderer,
        }
        row["submission_ready"] = all(
            row[key]
            for key in [
                "nonblank",
                "meets_300dpi_at_7in",
                "rgb",
                "svg_nonempty",
                "svg_text_editable",
                "pdf_nonempty",
                "pdf_readback",
            ]
        )
        rows.append(row)
    qa = pd.DataFrame(rows)
    qa_path = QA_DIR / "submission_qa.csv"
    qa.to_csv(qa_path, index=False)
    if not qa["submission_ready"].all():
        raise RuntimeError("One or more figures failed submission QA")
    color_sheet = build_contact_sheet([Path(paths["png"]) for paths in outputs.values()], QA_DIR / "color_contact_sheet.png")
    grayscale_sheet = build_contact_sheet(grayscale_paths, QA_DIR / "grayscale_contact_sheet.png")
    pdf_sheet = build_contact_sheet(pdf_paths, QA_DIR / "pdf_readback_contact_sheet.png")
    return {
        "submission_qa": "qa/submission_qa.csv",
        "color_contact_sheet": "qa/color_contact_sheet.png",
        "grayscale_contact_sheet": "qa/grayscale_contact_sheet.png",
        "pdf_readback_contact_sheet": "qa/pdf_readback_contact_sheet.png",
    }


def write_manifest(
    data: dict[str, pd.DataFrame],
    outputs: dict[str, dict[str, str]],
    support: dict[str, str],
    qa: dict[str, str],
) -> None:
    profiles = data["profiles"]
    candidates = data["candidates"]
    relative_outputs = {
        figure: {extension: str(Path(path).relative_to(OUT_DIR)) for extension, path in paths.items()}
        for figure, paths in outputs.items()
    }
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "x1-x100",
        "data_dir": str(DATA_DIR.relative_to(ROOT)),
        "counts": {
            "candidates": int(len(candidates)),
            "profiles": int(len(profiles)),
            "frames": int(profiles["converged_frames"].sum()),
            "design_indices": int(candidates["design_index"].nunique()),
            "activity_pairs": int(candidates["activity_pair"].nunique()),
            "activities": int(profiles["activity"].nunique()),
            "topologies": int(candidates["topology"].nunique()),
        },
        "figures": relative_outputs,
        "manuscript": support,
        "qa": qa,
        "boundary": "Computational screening descriptors are not experimental activity measurements or DFT barriers.",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def write_readme(data: dict[str, pd.DataFrame], outputs: dict[str, dict[str, str]]) -> None:
    candidates = data["candidates"]
    profiles = data["profiles"]
    try:
        data_label = DATA_DIR.relative_to(ROOT).as_posix()
    except ValueError:
        data_label = str(DATA_DIR)
    lines = [
        "# SPJ Main Figures",
        "",
        f"All figures read only from `{data_label}`.",
        "",
        f"- Retained canonical candidate records: {len(candidates)}",
        f"- Complete profiles: {len(profiles)}",
        f"- Converged frames: {int(profiles['converged_frames'].sum())}",
        "- Figure width: 7 inches",
        "- PNG: 450 dpi RGB",
        "- SVG/PDF: hybrid vector containers with editable text and line work",
        "",
        "## Figures",
        "",
    ]
    for figure, paths in outputs.items():
        lines.append(f"- {figure}: {Path(paths['png']).name}")
    lines.extend(
        [
            "",
            "## Support",
            "",
            "- Captions: manuscript/figure_captions_latest.md",
            "- Panel source-data index: manuscript/panel_source_data_index.csv",
            "- Figure contracts: manuscript/fig1_figure_contract.md through manuscript/fig5_figure_contract.md",
        ]
    )
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    configure_paths(data_dir=args.data_dir, out_dir=args.out_dir)
    configure_style()
    reset_output()
    data = load_dataset()
    outputs = {
        "Fig. 1": draw_fig1(data),
        "Fig. 2": draw_fig2(data),
        "Fig. 3": draw_fig3(data),
        "Fig. 4": draw_fig4(data),
        "Fig. 5": draw_fig5(data),
    }
    support = write_manuscript_support(data)
    qa = write_qa(outputs)
    write_manifest(data, outputs, support, qa)
    write_readme(data, outputs)
    print(
        json.dumps(
            {
                "data_dir": str(DATA_DIR),
                "output_dir": str(OUT_DIR),
                "counts": {
                    "candidates": len(data["candidates"]),
                    "profiles": len(data["profiles"]),
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
