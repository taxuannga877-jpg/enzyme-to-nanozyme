from __future__ import annotations

from typing import Any, MutableMapping


NON_FE_CU_COLORS = {
    "ink": "#26312d",
    "muted": "#66716b",
    "grid": "#e4ebe6",
    "records": "#8aa399",
    "mace": "#4e7f9e",
    "calculable": "#d07035",
    "complete": "#0d7f73",
    "residual": "#b05a5a",
    "first": "#49656b",
    "gfn2": "#2f7f95",
    "gfn1": "#a45a3f",
    "extended": "#8f4d81",
    "boundary": "#8f6b2f",
    "planned": "#5b6f2a",
    "source": "#8f4d81",
}

NON_FE_CU_SOURCE_ORDER = [
    "first_pass",
    "gfn2_deep",
    "gfn1_scf_fallback",
    "gfn1_extended",
    "gfn2_extended",
]

NON_FE_CU_SOURCE_COLORS = {
    "first_pass": NON_FE_CU_COLORS["first"],
    "gfn2_deep": NON_FE_CU_COLORS["gfn2"],
    "gfn1_scf_fallback": NON_FE_CU_COLORS["gfn1"],
    "gfn1_extended": NON_FE_CU_COLORS["extended"],
    "gfn2_extended": NON_FE_CU_COLORS["boundary"],
}

X1_X30_COLORS = {
    "ink": "#1f2933",
    "muted": "#6b7280",
    "grid": "#d6d9df",
    "x1": "#2f6f73",
    "x26": "#b65c32",
    "total": "#334e68",
    "complete": "#2f6f73",
    "residual": "#b65c32",
    "first_pass": "#4c78a8",
    "gfn1_scf_fallback": "#59a14f",
    "gfn2_deep": "#f28e2b",
    "gfn1_extended": "#b07aa1",
    "gfn2_extended": "#e15759",
    "neutral": "#e8eaef",
}

MULTIVARIANT_RESULT_COLORS = {
    "complete": "#0072B2",
    "incomplete": "#D55E00",
    "success": "#009E73",
    "neutral": "#7A869A",
    "accent": "#CC79A7",
    "warning": "#E69F00",
    "line": "#D8DEE9",
    "text": "#1F2937",
}

PHYSCHEM_COMPARISON_COLORS = {
    "passed": "#2A9D8F",
    "failed": "#E76F51",
    "rejected": "#F4A261",
    "not_constructible": "#7B8794",
    "legacy": "#B0BEC5",
    "new": "#264653",
    "allowed": "#D4E6D4",
    "cu": "#D55E00",
    "fe": "#0072B2",
    "zn": "#009E73",
    "mn": "#CC79A7",
    "ni": "#E69F00",
    "support": "#C7CDD1",
    "donor": "#4A90D9",
}

VISUAL_ATLAS_SOURCE_COLORS = {
    "first_pass": "#275f8f",
    "gfn1_scf_fallback": "#0d7f73",
    "gfn2_deep": "#d07035",
    "gfn1_extended": "#b8860b",
    "gfn2_extended": "#8f4d81",
    "missing": "#9aa4a1",
}

VISUAL_ATLAS_ACTIVITY_COLORS = {
    "Catalase": "#0d7f73",
    "Oxidase": "#275f8f",
    "Peroxidase": "#d07035",
    "Glutathione Peroxidase": "#8f4d81",
    "Glucose Oxidase": "#b8860b",
    "DNase": "#9b4f42",
}

METAL_DIVERSE_METAL_COLORS = {
    "FE": "#275f8f",
    "CU": "#d07035",
    "MN": "#0d7f73",
    "CO": "#8f4d81",
    "NI": "#5b6f2a",
    "ZN": "#b8860b",
}

METAL_DIVERSE_ACTIVITY_COLORS = {
    **VISUAL_ATLAS_ACTIVITY_COLORS,
    "Phosphatase": "#5b6f2a",
    "Urease": "#6d6658",
}

INTEGRATED_ACTIVITY_COLORS = {
    "Catalase": "#D55E00",
    "DNase": "#CC79A7",
    "Glucose Oxidase": "#E69F00",
    "Glutathione Peroxidase": "#F0E442",
    "Oxidase": "#0072B2",
    "Peroxidase": "#009E73",
    "Phosphatase": "#56B4E9",
}

INTEGRATED_SOURCE_COLORS = {
    "first_pass_success": "#0072B2",
    "fresh_rerun_gfn2": "#009E73",
    "fresh_rerun_gfn1": "#D55E00",
}

PAPER_OKABE_COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "vermillion": "#D55E00",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "black": "#000000",
}

PAPER_TIER_COLORS = {
    "A": PAPER_OKABE_COLORS["green"],
    "B": PAPER_OKABE_COLORS["sky"],
    "C": PAPER_OKABE_COLORS["orange"],
    "D": PAPER_OKABE_COLORS["purple"],
}

PAPER_METHOD_COLORS = {
    "GFN2-xTB": PAPER_OKABE_COLORS["blue"],
    "GFN1-xTB": PAPER_OKABE_COLORS["red"],
}

PAPER_CATEGORY_COLORS = {
    "metal_sites": PAPER_OKABE_COLORS["blue"],
    "catalytic_sites": PAPER_OKABE_COLORS["orange"],
    "binding_sites": PAPER_OKABE_COLORS["green"],
}

SPJ_WORKBENCH_PALETTE = {
    "ink": "#17211d",
    "muted": "#5f6b65",
    "line": "#cfd8d2",
    "paper": "#ffffff",
    "soft": "#f5f7f4",
    "teal": "#087f73",
    "blue": "#305a8d",
    "amber": "#c56d2a",
    "green": "#27764d",
    "gray": "#8b9690",
}


def apply_evidence_rcparams(
    rcparams: MutableMapping[str, Any],
    *,
    colors: dict[str, str] = NON_FE_CU_COLORS,
    title_size: int = 13,
    label_size: int = 10,
    tick_size: int = 9,
    legend_size: int = 9,
    savefig_dpi: int | None = None,
    svg_fonttype_none: bool = True,
) -> None:
    params = {
        "font.family": "DejaVu Sans",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.edgecolor": colors["ink"],
        "axes.labelcolor": colors["ink"],
        "xtick.color": colors["ink"],
        "ytick.color": colors["ink"],
        "axes.titlecolor": colors["ink"],
        "axes.titleweight": "bold",
        "axes.titlesize": title_size,
        "axes.labelsize": label_size,
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        "legend.fontsize": legend_size,
    }
    if savefig_dpi is not None:
        params["savefig.dpi"] = savefig_dpi
    if svg_fonttype_none:
        params["svg.fonttype"] = "none"
    rcparams.update(params)


def style_evidence_axis(
    ax,
    *,
    colors: dict[str, str] = NON_FE_CU_COLORS,
    axis: str = "y",
    linewidth: float = 0.8,
) -> None:
    ax.grid(axis=axis, color=colors["grid"], linewidth=linewidth)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def apply_x1_x30_rcparams(rcparams: MutableMapping[str, Any]) -> None:
    rcparams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def apply_multivariant_result_rcparams(rcparams: MutableMapping[str, Any]) -> None:
    rcparams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": MULTIVARIANT_RESULT_COLORS["line"],
            "axes.labelcolor": MULTIVARIANT_RESULT_COLORS["text"],
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "svg.fonttype": "none",
        }
    )


def apply_physchem_comparison_rcparams(rcparams: MutableMapping[str, Any]) -> None:
    rcparams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "font.size": 10,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "legend.fontsize": 8.2,
            "legend.frameon": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.14,
            "grid.linestyle": "-",
        }
    )
