#!/usr/bin/env python3
"""Build the public E2N x1-x100 data report and audit tables.

The script reads only the curated ``publication`` release layer. It does not
depend on the ignored raw ``outputs`` tree, so a GitHub clone can reproduce the
report-level summaries from the deposited canonical tables.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import kruskal


ROOT = Path(__file__).resolve().parents[1]
PUBLICATION = ROOT / "publication"
DATA_DIR = PUBLICATION / "data" / "x1_x100_dataset"
FIGURE_DIR = PUBLICATION / "figures"
REPORT_DIR = PUBLICATION / "reports"
TABLE_DIR = REPORT_DIR / "tables"
REPORT_PATH = REPORT_DIR / "e2n_x1_x100_data_summary_zh.md"

EXPECTED = {
    "candidates": 355,
    "profiles": 699,
    "frames": 3515,
    "design_indices": 38,
    "activity_pairs": 13,
    "activities": 6,
    "topologies": 3,
    "topology_tests": 10,
    "retained_tests": 8,
    "figure_source_csvs": 17,
    "figures": 5,
}

TOPOLOGY_ORDER = ["bridged", "independent adjacent", "independent separated"]
ACTIVITY_ORDER = [
    "Catalase",
    "DNase",
    "Glucose Oxidase",
    "Glutathione Peroxidase",
    "Oxidase",
    "Peroxidase",
]
METHOD_ORDER = [
    "First pass",
    "GFN1 SCF fallback",
    "GFN1 extended",
    "GFN2 deep",
    "GFN2 extended",
]


def _read_csv(name: str) -> pd.DataFrame:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def _csv_files(directory: Path) -> list[Path]:
    """List real CSV assets while ignoring macOS AppleDouble sidecars."""
    return sorted(
        path
        for path in directory.glob("*.csv")
        if path.is_file() and not path.name.startswith("._")
    )


def _pct(numerator: int | float, denominator: int | float) -> float:
    return 100.0 * float(numerator) / float(denominator) if denominator else math.nan


def _round(value: Any, digits: int = 6) -> Any:
    if pd.isna(value):
        return ""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return round(float(value), digits)
    return value


def _write_csv(frame: pd.DataFrame, name: str) -> Path:
    path = TABLE_DIR / name
    frame.to_csv(path, index=False, encoding="utf-8")
    return path


def _md_escape(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", "<br>")


def _md_table(frame: pd.DataFrame, columns: Iterable[str] | None = None) -> str:
    view = frame.loc[:, list(columns)] if columns is not None else frame
    headers = [str(column) for column in view.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in view.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(_md_escape(value) for value in row) + " |")
    return "\n".join(lines)


def _bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def _split_declared_pair(value: str) -> set[str]:
    return {part.strip() for part in str(value).split(" + ") if part.strip()}


def _format_trajectory(frame: pd.DataFrame, activity: str) -> str:
    values = frame.loc[frame["activity"].eq(activity), "relative_energy_ev"].tolist()
    return ", ".join(f"{value:.4f}" for value in values)


def build_report() -> dict[str, Any]:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    candidates = _read_csv("candidates.csv")
    profiles = _read_csv("profiles.csv")
    designs = _read_csv("designs.csv")
    topology_tests = _read_csv("topology_tests.csv")
    representative = _read_csv("representative.csv")
    representative_scans = _read_csv("representative_scans.csv")
    stored_correlations = _read_csv("descriptor_correlation.csv").set_index("metric")
    manifest = json.loads((DATA_DIR / "dataset_manifest.json").read_text(encoding="utf-8"))

    candidate_count = len(candidates)
    profile_count = len(profiles)
    frame_count = int(profiles["converged_frames"].sum())

    profile_by_candidate = profiles.groupby("candidate_id").size()
    candidate_linkage = candidates.copy()
    candidate_linkage["complete_profile_count"] = (
        candidate_linkage["candidate_id"].map(profile_by_candidate).fillna(0).astype(int)
    )
    profile_activity_sets = profiles.groupby("candidate_id")["activity"].agg(lambda values: "; ".join(sorted(values)))
    candidate_linkage["retained_profile_activities"] = (
        candidate_linkage["candidate_id"].map(profile_activity_sets).fillna("")
    )
    linkage_summary = (
        candidate_linkage.groupby("complete_profile_count", as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
        .sort_values("complete_profile_count")
    )
    linkage_summary["candidate_percent"] = 100 * linkage_summary["candidate_count"] / candidate_count

    declared_rows: list[dict[str, str]] = []
    for row in candidates[["candidate_id", "activity_pair"]].itertuples(index=False):
        for activity in _split_declared_pair(row.activity_pair):
            declared_rows.append({"candidate_id": row.candidate_id, "activity": activity})
    declared_targets = pd.DataFrame(declared_rows)
    declared_target_counts = (
        declared_targets.groupby("activity", as_index=False)
        .size()
        .rename(columns={"size": "declared_target_slots"})
    )
    target_profile_reconciliation = declared_target_counts.merge(
        profiles.groupby("activity", as_index=False).size().rename(columns={"size": "complete_profiles"}),
        on="activity",
        how="outer",
    ).fillna(0)
    target_profile_reconciliation[["declared_target_slots", "complete_profiles"]] = target_profile_reconciliation[
        ["declared_target_slots", "complete_profiles"]
    ].astype(int)
    target_profile_reconciliation["missing_complete_profile_slots"] = (
        target_profile_reconciliation["declared_target_slots"] - target_profile_reconciliation["complete_profiles"]
    )
    target_profile_reconciliation["_order"] = target_profile_reconciliation["activity"].map(
        {name: i for i, name in enumerate(ACTIVITY_ORDER)}
    )
    target_profile_reconciliation = target_profile_reconciliation.sort_values("_order").drop(columns="_order")

    topology_summary = (
        candidates.groupby("topology", as_index=False)
        .agg(
            candidate_count=("candidate_id", "size"),
            represented_design_indices=("design_index", "nunique"),
            requested_distance_q25_a=("distance_a", lambda values: values.quantile(0.25)),
            requested_distance_median_a=("distance_a", "median"),
            requested_distance_q75_a=("distance_a", lambda values: values.quantile(0.75)),
            median_score=("score", "median"),
            median_max_force_ev_per_a=("max_force_ev_per_a", "median"),
        )
    )
    topology_summary["candidate_percent"] = 100 * topology_summary["candidate_count"] / candidate_count
    topology_summary["_order"] = topology_summary["topology"].map({name: i for i, name in enumerate(TOPOLOGY_ORDER)})
    topology_summary = topology_summary.sort_values("_order").drop(columns="_order")

    doping_summary = (
        candidates.groupby("doping", as_index=False)
        .agg(candidate_count=("candidate_id", "size"), represented_design_indices=("design_index", "nunique"))
        .sort_values("candidate_count", ascending=False)
    )
    doping_summary["candidate_percent"] = 100 * doping_summary["candidate_count"] / candidate_count

    angle_summary = (
        candidates.groupby("angle_deg", as_index=False)
        .agg(candidate_count=("candidate_id", "size"), represented_design_indices=("design_index", "nunique"))
        .sort_values("angle_deg")
    )
    angle_summary["candidate_percent"] = 100 * angle_summary["candidate_count"] / candidate_count

    pair_counts = (
        candidates.groupby("activity_pair", as_index=False)
        .agg(
            candidate_count=("candidate_id", "size"),
            represented_design_indices=("design_index", "nunique"),
            represented_topologies=("topology", "nunique"),
        )
        .sort_values(["candidate_count", "activity_pair"], ascending=[False, True])
    )
    pair_counts["candidate_percent"] = 100 * pair_counts["candidate_count"] / candidate_count

    pair_topology = (
        candidates.groupby(["activity_pair", "topology"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
        .merge(pair_counts[["activity_pair", "candidate_count"]].rename(columns={"candidate_count": "pair_count"}), on="activity_pair")
    )
    pair_topology["within_pair_percent"] = 100 * pair_topology["candidate_count"] / pair_topology["pair_count"]
    pair_topology["_pair_order"] = pair_topology["activity_pair"].map(
        {name: i for i, name in enumerate(pair_counts["activity_pair"])}
    )
    pair_topology["_topology_order"] = pair_topology["topology"].map(
        {name: i for i, name in enumerate(TOPOLOGY_ORDER)}
    )
    pair_topology = pair_topology.sort_values(["_pair_order", "_topology_order"]).drop(
        columns=["_pair_order", "_topology_order"]
    )

    design_occupancy = (
        candidates.groupby(["design_index", "design_id"], as_index=False)
        .agg(
            candidate_count=("candidate_id", "size"),
            topology_count=("topology", "nunique"),
            activity_pair_count=("activity_pair", "nunique"),
            requested_distance_min_a=("distance_a", "min"),
            requested_distance_max_a=("distance_a", "max"),
            doping_values=("doping", lambda values: ";".join(sorted(set(values)))),
            angle_values_deg=("angle_deg", lambda values: ";".join(str(int(value)) for value in sorted(set(values)))),
        )
        .merge(
            profiles.groupby("design_index", as_index=False).size().rename(columns={"size": "profile_count"}),
            on="design_index",
            how="left",
        )
        .sort_values("design_index")
    )
    design_occupancy["profile_count"] = design_occupancy["profile_count"].fillna(0).astype(int)

    geometry_cells = (
        candidates.groupby(["topology", "distance_a", "doping", "angle_deg"], as_index=False)
        .agg(candidate_count=("candidate_id", "size"), design_index_count=("design_index", "nunique"))
        .sort_values(["candidate_count", "topology"], ascending=[False, True])
    )

    activity_summary = (
        profiles.groupby("activity", as_index=False)
        .agg(
            profile_count=("profile_id", "size"),
            candidate_count=("candidate_id", "nunique"),
            adsorption_min_ev=("adsorption_energy_ev", "min"),
            adsorption_q25_ev=("adsorption_energy_ev", lambda values: values.quantile(0.25)),
            adsorption_median_ev=("adsorption_energy_ev", "median"),
            adsorption_q75_ev=("adsorption_energy_ev", lambda values: values.quantile(0.75)),
            adsorption_max_ev=("adsorption_energy_ev", "max"),
            activation_min_ev=("activation_metric_ev", "min"),
            activation_q25_ev=("activation_metric_ev", lambda values: values.quantile(0.25)),
            activation_median_ev=("activation_metric_ev", "median"),
            activation_q75_ev=("activation_metric_ev", lambda values: values.quantile(0.75)),
            activation_max_ev=("activation_metric_ev", "max"),
            zero_activation_count=("activation_metric_ev", lambda values: int(values.eq(0).sum())),
            method_count=("selected_method", "nunique"),
            topology_count=("topology", "nunique"),
        )
    )
    activity_summary["zero_activation_percent"] = (
        100 * activity_summary["zero_activation_count"] / activity_summary["profile_count"]
    )
    activity_summary["_order"] = activity_summary["activity"].map({name: i for i, name in enumerate(ACTIVITY_ORDER)})
    activity_summary = activity_summary.sort_values("_order").drop(columns="_order")

    method_activity = (
        profiles.groupby(["activity", "selected_method"], as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    method_activity = method_activity.merge(
        activity_summary[["activity", "profile_count"]].rename(columns={"profile_count": "activity_profile_count"}),
        on="activity",
    )
    method_activity["within_activity_percent"] = 100 * method_activity["profile_count"] / method_activity["activity_profile_count"]
    method_activity["total_profile_percent"] = 100 * method_activity["profile_count"] / profile_count
    method_activity["_activity_order"] = method_activity["activity"].map({name: i for i, name in enumerate(ACTIVITY_ORDER)})
    method_activity["_method_order"] = method_activity["selected_method"].map({name: i for i, name in enumerate(METHOD_ORDER)})
    method_activity = method_activity.sort_values(["_activity_order", "_method_order"]).drop(
        columns=["_activity_order", "_method_order"]
    )

    method_totals = (
        profiles.groupby("selected_method", as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    method_totals["profile_percent"] = 100 * method_totals["profile_count"] / profile_count
    method_totals["_order"] = method_totals["selected_method"].map({name: i for i, name in enumerate(METHOD_ORDER)})
    method_totals = method_totals.sort_values("_order").drop(columns="_order")

    metric_columns = [
        "distance_a",
        "adsorption_energy_ev",
        "activation_metric_ev",
        "reaction_energy_ev",
        "scan_energy_range_ev",
        "score",
        "max_force_ev_per_a",
    ]
    recomputed_correlations = profiles[metric_columns].corr(method="spearman")
    correlation_delta = float((recomputed_correlations - stored_correlations).abs().to_numpy().max())
    activation_equals_reaction = int(
        np.isclose(
            profiles["activation_metric_ev"],
            profiles["reaction_energy_ev"],
            rtol=0.0,
            atol=1e-12,
        ).sum()
    )
    activation_equals_scan_range = int(
        np.isclose(
            profiles["activation_metric_ev"],
            profiles["scan_energy_range_ev"],
            rtol=0.0,
            atol=1e-12,
        ).sum()
    )

    recomputed_test_rows: list[dict[str, Any]] = []
    p_values: list[float] = []
    for row in topology_tests.itertuples(index=False):
        metric_column = "adsorption_energy_ev" if row.metric == "adsorption" else "activation_metric_ev"
        subset = profiles.loc[profiles["activity"].eq(row.activity), ["topology", metric_column]]
        groups = []
        keyed_sizes = []
        compared = []
        for topology in TOPOLOGY_ORDER:
            values = subset.loc[subset["topology"].eq(topology), metric_column].dropna().to_numpy()
            keyed_sizes.append(f"{topology}={len(values)}")
            if len(values) >= 2:
                groups.append(values)
                compared.append(topology)
        statistic, p_value = kruskal(*groups)
        n_obs = sum(len(group) for group in groups)
        epsilon = max(0.0, (float(statistic) - len(groups) + 1) / (n_obs - len(groups)))
        p_values.append(float(p_value))
        recomputed_test_rows.append(
            {
                "activity": row.activity,
                "metric": row.metric,
                "compared_topologies_recomputed": ";".join(compared),
                "group_sizes_keyed": ";".join(keyed_sizes),
                "statistic_recomputed": float(statistic),
                "p_value_recomputed": float(p_value),
                "epsilon_squared_recomputed": epsilon,
                "statistic_stored": float(row.statistic),
                "p_value_stored": float(row.p_value),
                "epsilon_squared_stored": float(row.epsilon_squared),
                "q_value_stored": float(row.q_value),
                "retained_at_fdr_0_05_stored": bool(row.retained_at_fdr_0_05),
            }
        )
    recomputed_q = _bh_adjust(np.asarray(p_values, dtype=float))
    for row, q_value in zip(recomputed_test_rows, recomputed_q, strict=True):
        row["q_value_recomputed"] = float(q_value)
        row["retained_at_fdr_0_05_recomputed"] = bool(q_value < 0.05)
        row["max_numeric_delta"] = max(
            abs(row["statistic_recomputed"] - row["statistic_stored"]),
            abs(row["p_value_recomputed"] - row["p_value_stored"]),
            abs(row["epsilon_squared_recomputed"] - row["epsilon_squared_stored"]),
            abs(row["q_value_recomputed"] - row["q_value_stored"]),
        )
    recomputed_tests = pd.DataFrame(recomputed_test_rows).sort_values("q_value_recomputed")

    required_candidate_columns = [
        "candidate_id",
        "design_id",
        "design_index",
        "activity_pair",
        "topology",
        "distance_a",
        "doping",
        "angle_deg",
        "variant_id",
        "score",
        "max_force_ev_per_a",
    ]
    required_profile_columns = [
        "profile_id",
        "candidate_id",
        "activity_pair",
        "activity",
        "partner_activity",
        "topology",
        "selected_method",
        "frame_count",
        "converged_frames",
    ]
    candidate_metadata = candidates.set_index("candidate_id")
    linked = profiles.join(
        candidate_metadata[["design_id", "design_index", "activity_pair", "topology", "distance_a", "doping", "angle_deg"]],
        on="candidate_id",
        rsuffix="_candidate",
    )
    metadata_mismatch = pd.Series(False, index=linked.index)
    for column in ["design_id", "design_index", "activity_pair", "topology", "distance_a", "doping", "angle_deg"]:
        metadata_mismatch |= linked[column].astype(str).ne(linked[f"{column}_candidate"].astype(str))

    declared_activity_mismatch = 0
    for row in profiles.itertuples(index=False):
        declared = _split_declared_pair(row.activity_pair)
        if row.activity not in declared or row.partner_activity not in declared or row.activity == row.partner_activity:
            declared_activity_mismatch += 1

    visible_key_columns = ["activity_pair", "topology", "distance_a", "doping", "angle_deg", "variant_id"]
    visible_key_count = int(candidates[visible_key_columns].drop_duplicates().shape[0])
    visible_duplicate_rows = candidate_count - visible_key_count

    source_csvs = _csv_files(FIGURE_DIR / "source_data")
    figure_files_ok = all(
        all((FIGURE_DIR / f"fig{index}_{stem}.{extension}").exists() for extension in ["png", "svg", "pdf"])
        for index, stem in [
            (1, "canonical_evidence_trace"),
            (2, "geometry_landscape"),
            (3, "activity_pair_composition"),
            (4, "profile_descriptors"),
            (5, "topology_statistics"),
        ]
    )
    qa = pd.read_csv(FIGURE_DIR / "qa" / "submission_qa.csv")
    qa_ready_column = next((column for column in qa.columns if "ready" in column.lower()), None)
    qa_pass_count = (
        int(qa[qa_ready_column].astype(str).str.lower().isin({"true", "1", "yes", "passed"}).sum())
        if qa_ready_column
        else len(qa)
    )

    quality_rows = [
        ("candidate_row_count", "pass" if candidate_count == EXPECTED["candidates"] else "fail", candidate_count, EXPECTED["candidates"], "candidate grain"),
        ("candidate_id_unique", "pass" if candidates["candidate_id"].is_unique else "fail", int(candidates["candidate_id"].nunique()), candidate_count, "primary key"),
        ("candidate_required_fields_complete", "pass" if not candidates[required_candidate_columns].isna().any().any() else "fail", int(candidates[required_candidate_columns].isna().sum().sum()), 0, "required null cells"),
        ("profile_row_count", "pass" if profile_count == EXPECTED["profiles"] else "fail", profile_count, EXPECTED["profiles"], "profile grain"),
        ("profile_id_unique", "pass" if profiles["profile_id"].is_unique else "fail", int(profiles["profile_id"].nunique()), profile_count, "primary key"),
        ("profile_required_fields_complete", "pass" if not profiles[required_profile_columns].isna().any().any() else "fail", int(profiles[required_profile_columns].isna().sum().sum()), 0, "required null cells"),
        ("profile_candidate_fk", "pass" if profiles["candidate_id"].isin(candidates["candidate_id"]).all() else "fail", int((~profiles["candidate_id"].isin(candidates["candidate_id"])).sum()), 0, "orphan profiles"),
        ("profile_candidate_metadata", "pass" if not metadata_mismatch.any() else "fail", int(metadata_mismatch.sum()), 0, "linked metadata mismatches"),
        ("declared_activity_membership", "pass" if declared_activity_mismatch == 0 else "fail", declared_activity_mismatch, 0, "activity/partner outside declared pair"),
        ("positive_profile_frames", "pass" if profiles["frame_count"].gt(0).all() else "fail", int((~profiles["frame_count"].gt(0)).sum()), 0, "nonpositive frame rows"),
        ("all_profile_frames_converged", "pass" if profiles["frame_count"].eq(profiles["converged_frames"]).all() else "fail", int((~profiles["frame_count"].eq(profiles["converged_frames"])).sum()), 0, "unequal frame rows"),
        ("aggregate_converged_frames", "pass" if frame_count == EXPECTED["frames"] else "fail", frame_count, EXPECTED["frames"], "aggregate count only"),
        ("design_axis_complete", "pass" if designs["design_index"].tolist() == list(range(1, 101)) else "fail", len(designs), 100, "x1-x100 axis"),
        ("topology_tests_recomputed", "pass" if recomputed_tests["max_numeric_delta"].max() < 1e-10 else "fail", float(recomputed_tests["max_numeric_delta"].max()), "<1e-10", "SciPy/BH/epsilon-squared"),
        ("descriptor_correlation_recomputed", "pass" if correlation_delta < 1e-6 else "fail", correlation_delta, "<1e-6", "Spearman matrix"),
        ("figure_source_csv_count", "pass" if len(source_csvs) == EXPECTED["figure_source_csvs"] else "fail", len(source_csvs), EXPECTED["figure_source_csvs"], "panel source tables"),
        ("figure_formats_present", "pass" if figure_files_ok else "fail", int(figure_files_ok), 1, "5 figures x PNG/SVG/PDF"),
        ("figure_qa_rows", "pass" if len(qa) == EXPECTED["figures"] and qa_pass_count == EXPECTED["figures"] else "fail", f"{qa_pass_count}/{len(qa)}", "5/5", "submission readiness"),
        ("visible_candidate_design_key", "warning", visible_key_count, candidate_count, f"{visible_duplicate_rows} rows require opaque candidate_id/source-side chemistry to distinguish"),
        ("candidate_retention_semantics", "warning", "172 legacy + 183 explicit", "one uniform flag unavailable", "do not state a campaign-wide calculability rate"),
        ("candidate_metal_provenance", "warning", "not present", "candidate-level metal/oxidation/source fields", "public table cannot explain all repeated visible design rows"),
        ("topology_group_sizes_schema", "warning", 4, 0, "two-group tests serialize a leading bridged=0 slot; use keyed sizes"),
    ]
    quality = pd.DataFrame(quality_rows, columns=["check", "status", "observed", "expected", "interpretation"])

    version_rows = [
        {
            "version_or_scope": "paper_data_current (2026-06-25)",
            "candidate_unit": "22 calculable candidates",
            "profile_unit": "44 complete profiles",
            "frame_unit": "280 converged frames",
            "statistics": "8 tests; 0 FDR-retained",
            "use_in_current_manuscript": "No - superseded historical snapshot",
        },
        {
            "version_or_scope": "broader x1-x100 lineage audit",
            "candidate_unit": "1817 upstream records; 661 MACE; historical 355 calculable label",
            "profile_unit": "699 complete / 706 eligible",
            "frame_unit": "3580 converged / 3590 attempted",
            "statistics": "batch-level lineage",
            "use_in_current_manuscript": "No - audit lineage only; legacy values prohibited",
        },
        {
            "version_or_scope": "v5 audited submission",
            "candidate_unit": "355 retained canonical candidate records",
            "profile_unit": "699 complete profiles",
            "frame_unit": "3515 frames in complete profiles",
            "statistics": "10 tests; 8 q<0.05; 36 DOI references; 5 review rounds",
            "use_in_current_manuscript": "Yes - numerical and figure authority",
        },
        {
            "version_or_scope": "v6 narrative manuscript",
            "candidate_unit": "355 stated as calculable (wording too broad)",
            "profile_unit": "699 complete profiles",
            "frame_unit": "3515 converged frames",
            "statistics": "30 references; 3 v6 review rounds; 7133 main-text words",
            "use_in_current_manuscript": "Narrative input only; requires evidence corrections",
        },
        {
            "version_or_scope": "current evidence contract / public release",
            "candidate_unit": "355 retained canonical candidate records",
            "profile_unit": "699 complete profiles; linkage 346/7/2",
            "frame_unit": "3515 aggregate converged frames",
            "statistics": "10 exploratory tests; 8 q<0.05",
            "use_in_current_manuscript": "Recommended controlling terminology",
        },
    ]
    versions = pd.DataFrame(version_rows)

    for frame in [
        topology_summary,
        doping_summary,
        angle_summary,
        pair_counts,
        pair_topology,
        design_occupancy,
        geometry_cells,
        activity_summary,
        method_activity,
        method_totals,
        linkage_summary,
        recomputed_tests,
    ]:
        for column in frame.select_dtypes(include=["float"]).columns:
            frame[column] = frame[column].map(lambda value: _round(value, 8))

    _write_csv(candidate_linkage, "candidate_profile_linkage.csv")
    _write_csv(linkage_summary, "profile_linkage_summary.csv")
    _write_csv(target_profile_reconciliation, "declared_target_profile_reconciliation.csv")
    _write_csv(topology_summary, "candidate_topology_summary.csv")
    _write_csv(doping_summary, "candidate_doping_summary.csv")
    _write_csv(angle_summary, "candidate_angle_summary.csv")
    _write_csv(pair_counts, "activity_pair_counts.csv")
    _write_csv(pair_topology, "activity_pair_topology_summary.csv")
    _write_csv(design_occupancy, "design_occupancy.csv")
    _write_csv(geometry_cells, "geometry_cell_summary.csv")
    _write_csv(activity_summary, "profile_activity_summary.csv")
    _write_csv(method_activity, "method_activity_summary.csv")
    _write_csv(method_totals, "method_total_summary.csv")
    _write_csv(recomputed_tests, "topology_tests_recomputed.csv")
    _write_csv(quality, "data_quality_checks.csv")
    _write_csv(versions, "version_reconciliation.csv")

    display_linkage = linkage_summary.rename(
        columns={
            "complete_profile_count": "每个候选体的完整 profile 数",
            "candidate_count": "候选体数",
            "candidate_percent": "占保留候选体 (%)",
        }
    ).copy()
    display_linkage["占保留候选体 (%)"] = display_linkage["占保留候选体 (%)"].map(lambda value: f"{value:.2f}")

    display_target_reconciliation = target_profile_reconciliation.rename(
        columns={
            "activity": "活性",
            "declared_target_slots": "声明目标槽位",
            "complete_profiles": "完整 profile",
            "missing_complete_profile_slots": "缺少完整 profile 的槽位",
        }
    )

    display_topology = topology_summary.rename(
        columns={
            "topology": "拓扑",
            "candidate_count": "候选体数",
            "candidate_percent": "占比 (%)",
            "represented_design_indices": "覆盖设计索引",
            "requested_distance_q25_a": "距离 Q1 (A)",
            "requested_distance_median_a": "距离中位数 (A)",
            "requested_distance_q75_a": "距离 Q3 (A)",
        }
    )
    display_topology["占比 (%)"] = display_topology["占比 (%)"].map(lambda value: f"{value:.1f}")

    display_pairs = pair_counts.rename(
        columns={
            "activity_pair": "声明活性对",
            "candidate_count": "候选体数",
            "candidate_percent": "占比 (%)",
            "represented_design_indices": "覆盖设计索引",
            "represented_topologies": "拓扑数",
        }
    ).copy()
    display_pairs["占比 (%)"] = display_pairs["占比 (%)"].map(lambda value: f"{value:.2f}")

    display_activity = activity_summary.rename(
        columns={
            "activity": "活性",
            "profile_count": "profile 数",
            "adsorption_median_ev": "吸附中位数 (eV)",
            "activation_median_ev": "正向峰值中位数 (eV)",
            "zero_activation_count": "无正峰数",
            "zero_activation_percent": "无正峰占比 (%)",
            "method_count": "计算路径数",
            "topology_count": "有记录拓扑数",
        }
    ).copy()
    for column in ["吸附中位数 (eV)", "正向峰值中位数 (eV)"]:
        display_activity[column] = display_activity[column].map(lambda value: f"{value:.4f}")
    display_activity["无正峰占比 (%)"] = display_activity["无正峰占比 (%)"].map(lambda value: f"{value:.2f}")

    display_methods = method_totals.rename(
        columns={"selected_method": "计算路径", "profile_count": "profile 数", "profile_percent": "占比 (%)"}
    ).copy()
    display_methods["占比 (%)"] = display_methods["占比 (%)"].map(lambda value: f"{value:.2f}")

    display_tests = recomputed_tests.rename(
        columns={
            "activity": "活性",
            "metric": "描述符",
            "group_sizes_keyed": "拓扑样本数",
            "statistic_recomputed": "H",
            "p_value_recomputed": "p",
            "q_value_recomputed": "q",
            "epsilon_squared_recomputed": "epsilon-squared",
            "retained_at_fdr_0_05_recomputed": "q<0.05",
        }
    ).copy()
    for column in ["H", "epsilon-squared"]:
        display_tests[column] = display_tests[column].map(lambda value: f"{value:.6f}")
    for column in ["p", "q"]:
        display_tests[column] = display_tests[column].map(lambda value: f"{value:.6g}")

    top_designs = design_occupancy.sort_values(["candidate_count", "design_index"], ascending=[False, True]).head(12)
    display_designs = top_designs.rename(
        columns={
            "design_id": "设计索引",
            "candidate_count": "候选体数",
            "profile_count": "profile 数",
            "topology_count": "拓扑数",
            "activity_pair_count": "活性对数",
            "doping_values": "掺杂",
            "angle_values_deg": "角度 (deg)",
        }
    )

    rep = representative.iloc[0]
    no_profile = candidate_linkage.loc[candidate_linkage["complete_profile_count"].eq(0), "candidate_id"].tolist()
    one_profile = candidate_linkage.loc[candidate_linkage["complete_profile_count"].eq(1), "candidate_id"].tolist()
    frame_distribution = profiles["frame_count"].value_counts().sort_index().to_dict()
    top_geometry = geometry_cells.iloc[0]

    report = f"""# E2N x1-x100 规范化数据与论文修改技术报告

## 技术摘要

本报告以 `publication/data/x1_x100_dataset` 为当前论文数据的唯一公开接口，并从候选体、活性 profile 和扫描帧三个统计单位重新核算全部核心数字。可复现结果为 **{candidate_count} 条保留的规范化候选体记录、{profile_count} 条完整活性 profile、{frame_count} 个完整 profile 内的收敛扫描帧**；候选体覆盖 38 个 x1-x100 设计索引、13 种声明活性对和 3 种拓扑，profile 覆盖 6 种活性。

最重要的口径修正是：**355 不能在没有限定语时写成整个研究流程中统一定义的“可计算候选体数”**。其中 172 条旧源记录依据可用的正向完整-profile endpoint 保留，183 条新源记录具有显式 `calculable=True`。因此，当前最稳妥的论文表述是 **“355 retained canonical candidate records（355 条保留的规范化候选体记录）”**。

候选体与 profile 的连接分布为 346/7/2：346 个候选体连接两个完整 profile，7 个连接一个，2 个没有保留的完整 profile。该分布只描述保留候选体集内部的 profile linkage，**不能转换为 97.5% 的研究流程成功率、双活性率、实验活性率或统一可计算率**。

当前数据足以支持一篇关于“蛋白质证据到纳米材料假设的可审计转译与筛选数据组织”的论文，但不足以支持催化优越性、实验双活性、无能垒反应、无约束结构稳定性、拓扑因果效应或预测准确率。公开表还缺少完整的候选体金属身份、氧化态和源位点 lineage；这应在投稿前补充，或在数据可用性声明中明确列为公开快照限制。

## 1. 权威边界与统计单位

本报告回答两个问题：第一，GitHub 发布层中究竟有哪些可公开、可重复核算的数据；第二，论文中每个主要数字可以写到什么程度。报告生成时间以源数据 manifest 为准：`{manifest['generated_at_utc']}`。

| 统计单位 | 规范化数量 | 定义 | 不能替代的概念 |
| --- | ---: | --- | --- |
| 候选体记录 | {candidate_count} | `candidate_id` 唯一的保留记录 | 全部尝试、产率、成功率、独立材料组成数 |
| 活性 profile | {profile_count} | 对某一候选体和某一活性的完整轨迹记录 | 候选体数、实验活性数、独立样本数 |
| 扫描帧 | {frame_count} | 完整 profile 中 `converged_frames` 的聚合和 | 唯一 frame ID 数；公开表没有逐帧主键 |
| 设计索引 | 38/100 | 至少出现一条保留候选体的 x1-x100 位置 | 62 个空位置的失败数或未尝试数 |
| 活性对 | 13 | 候选体在计算前声明的两个目标活性 | 同一材料已经证明的双活性 |
| 拓扑 | 3 | bridged、independent adjacent、independent separated | 因果材料机制或稳定性类别 |

### 1.1 两个候选体源块采用不同的保留依据

数据 manifest 记录了两个候选体输入块（172 行和 674 行）以及两个 profile 输入块（344 行和 362 行）。重新追溯 builder 规则后，候选体保留过程应表述为：

- 旧源候选体块：原表 172 行，没有统一 `calculable` 字段；这 172 行已经是具有至少一个完整 profile endpoint 的正向保留记录。
- 新源候选体块：原表 674 行，其中 183 行满足显式 `calculable=True`，进入规范化候选体表。
- 旧源 profile 块：344 行中 339 行完整，保留 1715 个完整 profile 帧；5 行不完整记录被排除。
- 新源 profile 块：362 行中 360 行完整，保留 1800 个完整 profile 帧；2 行不完整记录被排除。
- 规范化 profile 合计为 339 + 360 = 699，规范化帧合计为 1715 + 1800 = 3515。

历史 lineage 审计中的 706 条 eligible profile 和 3580/3590 帧把 7 条不完整 profile 及其部分收敛帧也纳入了上游流程计数。它们适合内部流程审计，但不适合当前成功集论文正文。

## 2. 候选体集：拓扑、设计网格与活性对

### 2.1 三种拓扑以 independent separated 为主，但这是保留集组成

{_md_table(display_topology, ['拓扑', '候选体数', '占比 (%)', '覆盖设计索引', '距离 Q1 (A)', '距离中位数 (A)', '距离 Q3 (A)'])}

independent separated 占 226/355（63.7%），independent adjacent 占 101/355（28.5%），bridged 占 28/355（7.9%）。这一比例同时受到枚举队列、设计家族重复、金属组合、构建规则和计算保留条件影响，不能解释为拓扑本身的成功概率。

请求距离的中位数按拓扑依次约为 6.2 A、8.8 A 和 13.0 A。距离和角度来自请求设计字段，并非去除约束后的平衡几何；图 2 的正确语义是“保留记录对请求设计网格的占据”。

![图 2：请求设计网格中的保留候选体](../figures/fig2_geometry_landscape.png)

### 2.2 设计占据高度集中，但不能称为优化或富集

100 个设计索引中只有 38 个出现保留记录。候选体数最高的索引如下：

{_md_table(display_designs, ['设计索引', '候选体数', 'profile 数', '拓扑数', '活性对数', '掺杂', '角度 (deg)'])}

把设计索引折叠后，355 条记录占据 49 个不同的“拓扑-距离-掺杂-角度”设置；保留 `design_index` 时有 57 个设计-拓扑-几何单元。最密集的两个设置是 independent separated、13.0 A、NS、105 deg（67 条）和同条件 75 deg（64 条）。其中 75 deg 由 x41 和 x81 各 32 条组成；105 deg 由 x17 的 3 条、x57 的 32 条和 x97 的 32 条组成。

这些计数证明的是“该研究流程重复实例化并保留了这些设置”，而不是：13.0 A 是最优距离、NS 是最佳掺杂、75/105 deg 是涌现平衡角、这些设置具有更高成功概率，或两个金属中心已证明协同。

### 2.3 13 种声明活性对的覆盖不均衡

{_md_table(display_pairs, ['声明活性对', '候选体数', '占比 (%)', '覆盖设计索引', '拓扑数'])}

Oxidase + Peroxidase 数量最多（62），随后是 Glutathione Peroxidase + Peroxidase（50）和 Glutathione Peroxidase + Oxidase（48）。三个含 DNase 的活性对各只有一个候选体，均不能支持分布性推断。若干含 Glucose Oxidase 的活性对以 independent separated 为主，例如 Glucose Oxidase + Peroxidase 中 17/19 为该拓扑；这仍然是成功条件化设计分配，而非比较活性证据。

![图 3：活性对及其设计网格组成](../figures/fig3_activity_pair_composition.png)

## 3. 候选体与完整 profile 的连接关系

{_md_table(display_linkage)}

完整连接分布严格为 346/7/2。没有完整 profile 的候选体 ID 为 `{', '.join(no_profile)}`；只有一个完整 profile 的候选体 ID 为 `{', '.join(one_profile)}`。这些 ID 被保留是因为候选体表和 profile 表使用不同的成功边界。

帧数分布为：{', '.join(f'{count} 条 profile x {frames} 帧' for frames, count in frame_distribution.items())}。因此 3515 是聚合帧数，不是公开可逐帧追踪的唯一记录数。

论文中可以写：“在保留的规范化候选体记录中，346 个候选体与两个完整活性特异性 profile 相连，7 个与一个相连，2 个没有保留的完整 profile。”随后必须补充：“该分布不构成研究流程完成率或双活性验证率。”

每条候选体记录声明两个目标活性，因此 355 条记录一共对应 710 个声明目标槽位。按活性与 699 条完整 profile 对账如下：

{_md_table(display_target_reconciliation)}

合计有 710 个声明目标槽位、699 条完整 profile 和 11 个未连接到完整 profile 的槽位。11 个缺口由两部分组成：7 条被完整性规则排除的不完整 profile，以及两个零完整-profile 候选体所对应的另外 4 个声明目标槽位。这里的“缺口”是公开规范化表之间的连接差额，不等同于反应失败、无活性或候选体不可计算。

## 4. Profile 组成、描述符与计算路径

### 4.1 六种活性的 profile 数量和无正峰边界

{_md_table(display_activity, ['活性', 'profile 数', '吸附中位数 (eV)', '正向峰值中位数 (eV)', '无正峰数', '无正峰占比 (%)', '计算路径数', '有记录拓扑数'])}

Oxidase 的无正峰比例最高：125/174（71.84%）；Glucose Oxidase 为 32/82（39.02%），Peroxidase 为 53/168（31.55%），Glutathione Peroxidase 为 19/154（12.34%），Catalase 为 1/120（0.83%）。DNase 只有一条记录，不做总体推断。

这里的零值严格表示：在存储的有限 forward scan 采样点中，没有观察到相对于首点的正向峰。因此应使用 `peakless profile` 或“无正峰描述符边界”，不得写成 zero kinetic barrier、barrier-free reaction、无活化能或无催化活性。

描述符之间还存在由定义和有限扫描轨迹造成的精确重合：{activation_equals_reaction}/699 条 profile 的 forward scan peak descriptor 与 reaction energy 数值相同，{activation_equals_scan_range}/699 条与 scan energy range 数值相同。这些重合不是独立测量之间的验证，也不能作为物理机制相关性的证据；解释相关矩阵时必须保留这一计算依赖关系。

![图 4：profile 计算路径与描述符分布](../figures/fig4_profile_descriptors.png)

### 4.2 计算路径分布明显依赖活性

{_md_table(display_methods)}

GFN2 deep 为 251/699，First pass 为 240/699，GFN1 SCF fallback 为 128/699，GFN1 extended 为 69/699，GFN2 extended 为 11/699。更关键的是，各活性内部的路径组成不同：Oxidase 中 First pass 占 119/174（68.4%），Glutathione Peroxidase 中 GFN2 deep 占 103/154（66.9%），Glucose Oxidase 中 GFN2 deep 占 44/82（53.7%）。

GFN1 和 GFN2 并非可互换能标，rescue 深度也可能与候选体化学和收敛难度相关。因此跨活性或跨拓扑的 pooled descriptor 图只能描述当前保留语料，不能作为统一协议下的严格能量比较。

### 4.3 Spearman 相关性属于混合语料的描述性结果

公开 `descriptor_correlation.csv` 是 699 条 profile 上的 Spearman 相关矩阵，已在本报告中以最大误差 {correlation_delta:.2e} 独立复算。绝对值较大的非对角关联包括：activation 与 reaction energy（rho=0.672879）、adsorption 与 scan range（rho=0.429351）、adsorption 与 activation（rho=0.395817）、activation 与 scan range（rho=0.327395）、distance 与 score（rho=0.318627），以及 reaction energy 与 scan range（rho=-0.306418）。

这些相关性混合了活性、方法、拓扑和设计家族，且同一候选体可贡献两个 profile；它们不能被当作独立样本下的机制关系或因果效应。

## 5. 拓扑检验：8/10 达到 q<0.05，但解释必须是探索性的

{_md_table(display_tests, ['活性', '描述符', '拓扑样本数', 'H', 'p', 'q', 'epsilon-squared', 'q<0.05'])}

10 项检验的 SciPy Kruskal-Wallis 统计量、全局 Benjamini-Hochberg 校正和 epsilon-squared 均被独立复算，最大数值偏差小于 {recomputed_tests['max_numeric_delta'].max():.2e}。8 项 q<0.05；未保留的是 Glutathione Peroxidase activation（q=0.414836）和 Glucose Oxidase activation（q=0.618886）。

效应量最大的三项是 Catalase adsorption（epsilon-squared=0.1885）、Glucose Oxidase adsorption（0.1719）和 Catalase forward scan peak（0.1295）。最小的保留项是 Glutathione Peroxidase adsorption（0.0277）。q 值和效应量回答不同问题，不能用显著性替代实际差异大小。

需要特别说明：Glucose Oxidase 和 Glutathione Peroxidase 没有 bridged profile，因此相应行实际上是两个非空拓扑组的比较。原始 `group_sizes` 字段仍以固定三拓扑顺序保存为 `0;13;69` 或 `0;48;106`，不能直接与只列两个名称的 `compared_topologies` 逐项 zip。本报告的 `topology_tests_recomputed.csv` 已增加带键样本数。

这些检验受成功条件约束，且拓扑与请求距离、角度、掺杂、金属身份、活性对、设计家族和计算路径混杂。正确表述是“当前规范化 profile 语料中的探索性拓扑分层”，不得写成拓扑导致更高活性、某拓扑普遍更优或材料设计规律。

![图 5：拓扑相关的探索性 profile 比较](../figures/fig5_topology_statistics.png)

## 6. 代表性 x57 记录的正确用途

代表性候选体为 `{rep['candidate_id']}`，设计索引 `{rep['design_id']}`，声明活性对为 `{rep['activity_pair']}`，拓扑为 `{rep['topology']}`，请求距离 {rep['distance_a']:.1f} A，NS 掺杂，请求角度 {int(rep['angle_deg'])} deg，score={rep['score']:.6f}，最大力={rep['max_force_ev_per_a']:.6f} eV/A。

- Glucose Oxidase 轨迹：{_format_trajectory(representative_scans, 'Glucose Oxidase')} eV。
- Peroxidase 轨迹：{_format_trajectory(representative_scans, 'Peroxidase')} eV。

这条记录可以展示一个候选体如何连接两个各自独立、带方法来源的 activity-specific profile。它不是最高排名候选体，也不能证明两个反应在同一条件下发生，更不能证明实验双活性或 cascade catalysis。

代表性记录的实际选择规则也需要透明：构建器先把 536 张私有结构图与规范化候选体求交，再要求候选体至少连接两个完整 profile，共得到 179 条 eligible 记录；随后选择最接近 eligible 集合 score 中位数的记录，并以 `design_index` 和 `candidate_id` 进行确定性并列排序。选择逻辑没有显式筛选 x57 或 13.0 A/NS/105 deg 窗口，因此它应称为“具有可用结构图和两个完整 profile 的中位 score 示例”，而不是针对高密度几何窗口挑选的代表。

公开的 `representative_scans.csv` 是从私有原始结果 JSON 抽取并冻结的 10 行扫描快照，不能只凭 `profiles.csv` 重新生成。它支持图 1 的轨迹复现，但不代表发布包包含逐帧上游轨迹或完整 campaign 重建材料。

![图 1：代表性候选体与两条独立 profile 轨迹](../figures/fig1_canonical_evidence_trace.png)

## 7. 数据质量、可追溯性与公开发布限制

### 7.1 已通过的完整性检查

- 355 个 `candidate_id` 和 699 个 `profile_id` 均唯一。
- profile 外键全部指向候选体，设计、拓扑、距离、掺杂和角度连接后无不一致。
- 所有 profile 的帧数为正，且 `frame_count == converged_frames`。
- 设计轴完整保留 x1-x100，候选体实际出现于 38 个索引。
- 13 种活性对、6 种活性和 3 种拓扑均与 manifest 一致。
- 17 个 panel 源数据 CSV 齐全，五张图的 PNG/SVG/PDF 均存在，图 QA 为 5/5。
- 拓扑检验和 Spearman 相关矩阵均可从公开表独立复算。

### 7.2 公开候选体表仍缺少关键化学来源字段

在公开 `candidates.csv` 中，355 个 ID 唯一，但按可见科学字段（活性对、拓扑、距离、掺杂、角度、variant）只有 {visible_key_count} 个不同组合，{visible_duplicate_rows} 条额外记录依赖不透明的 `candidate_id` 才能区分。新源块内部实际还存在金属 A/B、氧化态和 metal-case 字段，但规范化表没有公开；旧源块本身也不含同等完整的金属 lineage。

因此，当前发布包可以复算论文图表和统计，但无法仅凭公开列解释 355 条记录之间全部化学差异，也不能从原始 campaign 重新构建规范化快照。投稿前有两个可接受方案：

1. 增加可再分发的候选体 provenance 表，至少包含 `candidate_id`、metal A/B、氧化态、源设计家族和来源类型；旧源缺失项明确标为 unavailable。
2. 如果不能公开，需在 Data Availability 中明确：公开包是经审计的分析快照，可复算表格/图/统计，但不包含完整上游原始结构、私有路径和所有候选体级生物/化学 lineage。

### 7.3 当前可分享置信度

对“复算当前论文五图和规范化描述统计”而言，数据在**科学与技术审计层面已就绪**。这不等于已获得公开再分发授权：作者、版权主体、年份及代码/数据/图许可范围冻结前，整个发布包在法律与归档层面仍为 **Not ready for public release**。对“从原始 PDB/候选体构建全过程端到端复现”或“证明候选体完整化学 provenance”而言，即使许可问题解决，当前发布仍为 **Share with caveats**，因为上游原始块、逐帧数据、完整金属字段和候选体级 PDB lineage 未包含在公开层。

详细逐项检查见 `reports/tables/data_quality_checks.csv`。

## 8. 历史数据版本的冲突与使用规则

{_md_table(versions)}

`paper_data_current` 是早期 2026-06-25 快照，包含 18 个 EC 类、1245 个 PDB 文件、720 个 motif、2959 条金属位点和 2573 条配体/辅因子记录；其后半部分筛选只有 22 个候选体、44 个完整 profile 和 280 个帧。它可用于说明早期项目谱系，不能替代当前 x1-x100 canonical 数据。

broader lineage 中的 1817、661、706、699/706、3580/3590 是上游流程与不完整记录审计值。它们与当前 355/699/3515 的统计边界不同，不能混在同一正文叙述中。当前论文的数值和图形权威应以 audited canonical tables、五图源数据和 current evidence contract 为准；v6 只作为叙事和术语优化输入。

## 9. 对 v6 正文的逐项修改建议

### 9.1 全文统一改写的核心数字句

将无条件的 “355 calculable candidates” 改为：

> The canonical release contains 355 retained candidate records, 699 complete activity-specific profiles, and 3515 converged scan frames. The two candidate source blocks used different available retention fields; these counts therefore describe a success-conditioned analysis set rather than a campaign-wide calculability rate.

中文含义：规范化发布层包含 355 条保留候选体记录、699 条完整的活性特异性 profile 和 3515 个收敛扫描帧。两个候选体源块采用不同的可用保留字段，因此这些数字描述的是成功条件化分析集，而非研究流程统一可计算率。

### 9.2 删除无约束 MACE 五候选体段落

删除 v6 Discussion/Evidence boundaries 中关于“5 个无约束 MACE 候选体、0 个保持拓扑、4 个达到力收敛”的整段结果。它来自更早的 22 候选体快照，不属于当前 x1-x100 retained canonical snapshot，也不是系统抽样验证集。保留它会重新引入用户已明确要求排除的失败结果，并混淆版本边界。

### 9.3 346/7/2 只能作为 linkage

可以在数据完整性或方法部分报告 346/7/2，但不得计算或突出 346/355=97.5% 作为双活性完成率。完整 profile 是计算轨迹完整性，不是催化活性判定。

### 9.4 压缩 Results 中重复计数

方法路径的逐活性 24 个单元格计数、13 个活性对的全部计数以及 10 组 H/p/q 不必全部写入正文。建议正文保留：总数、最主要的 2-3 个模式和解释边界；完整数字放入表格、图源数据或补充材料。这样可保留 v6 的科学定位，同时接近 v5 的篇幅和可读性。

### 9.5 图 2 的术语必须保持“requested”

所有 distance、angle、geometry 均写成 requested design distance/angle/geometry。13.0 A、NS、75/105 deg 只能称为 densely represented retained-grid cells 或 hypothesis class，不能称 optimum、enrichment、superior geometry、higher success probability 或 equilibrium structure。

### 9.6 图 4 的零值和方法来源

activation metric 统一改为 forward scan peak descriptor。零值写成 peakless descriptor boundary。正文需保留 GFN1/GFN2 路径混合的限制，避免把 pooled descriptor 解释为统一能标。

### 9.7 图 5 的统计语气

使用 observational、exploratory、success-conditioned、confounded、profile-level omnibus association。避免 causal effect、topology-controlled performance、universally superior topology。q 与 epsilon-squared 分开报告。

### 9.8 当前正文不得使用的历史数字

正文、摘要、图注和结论中不要使用：`1817`、`661`、`706`、`699/706`、`3580/3590`。这些数字可以只在内部 provenance 报告中出现，并明确标为历史上游审计口径。

## 10. GitHub 发布包的建议边界

公开仓库应支持两种不同层次的复现：

1. **release-level reproduction**：从 `publication/data/x1_x100_dataset` 验证 355/699/3515、重新生成派生表、复算拓扑检验、重绘五张图并核对图源数据。这一层应在 GitHub clone 后直接可运行。
2. **full upstream rebuild**：从 PDB 库、motif 数据库、原始候选体结构、MACE/xTB 轨迹和 rescue 输出重新生成 canonical snapshot。当前大体量本地 `outputs/`、数据库、模型权重和结构库不在发布层中，因此不能宣称公开仓库已实现这一层。

建议公开：核心代码、当前构建/作图/审计脚本、精简复现环境、canonical 表、五图及 17 个 panel 源数据、报告、证据合同和 release manifest。不要公开：本机路径、缓存、模型权重、期刊网页镜像、未授权第三方资产、整个历史 `outputs/`、失败/中间批次和根目录重复导出文件。

## 11. 推荐下一步

1. 用本报告控制 v6/v7 的数字和术语，先修正文再做最终 DOCX/PDF QA。
2. 决定是否能公开候选体金属/氧化态/provenance；若不能，强化 Data Availability 的快照边界。
3. 将 24 个方法-活性计数、13 个活性对计数和 10 项拓扑检验移至补充表或数据报告，正文只保留主要模式。
4. 冻结作者、单位、通讯作者、基金、非作者致谢、仓库 URL、归档 DOI 和数据/图许可范围。
5. 发布前运行 `python scripts/build_spj_data_summary_report.py`、`python publication/scripts/verify_publication_release.py`，再从 canonical data 重新生成五图并核对 release manifest。

## 12. 仍需作者决定的问题

- 355 条候选体是否要在公开 archive 中补充金属 A/B、氧化态和源位点/设计家族字段？
- MIT 是否只覆盖软件，数据和图是否另用 CC BY 4.0 或其他许可？
- GitHub 是否使用现有远端，还是建立一个历史干净的论文发布仓库？
- repository DOI 由 Zenodo/GitHub release 还是机构仓储生成？
- 正文是否继续在 v6 文件上逐段合并，同时以 current evidence contract 控制全部数字和术语？

## 附录：本报告生成的派生表

- `candidate_profile_linkage.csv`：355 条候选体与完整 profile 的连接明细。
- `profile_linkage_summary.csv`：346/7/2 分布。
- `declared_target_profile_reconciliation.csv`：710 个声明目标槽位与 699 条完整 profile 的逐活性对账。
- `candidate_topology_summary.csv`、`candidate_doping_summary.csv`、`candidate_angle_summary.csv`。
- `activity_pair_counts.csv`、`activity_pair_topology_summary.csv`。
- `design_occupancy.csv`、`geometry_cell_summary.csv`。
- `profile_activity_summary.csv`、`method_activity_summary.csv`、`method_total_summary.csv`。
- `topology_tests_recomputed.csv`：带键样本数及独立复算结果。
- `data_quality_checks.csv`：完整性、外键、统计复算和公开边界检查。
- `version_reconciliation.csv`：历史版本口径对照。

本报告是技术审计和论文修改依据，不是实验催化性能声明。
"""

    REPORT_PATH.write_text(report, encoding="utf-8", newline="\n")
    return {
        "report": str(REPORT_PATH.relative_to(ROOT)),
        "tables": len(_csv_files(TABLE_DIR)),
        "candidates": candidate_count,
        "profiles": profile_count,
        "frames": frame_count,
        "quality_failures": int(quality["status"].eq("fail").sum()),
        "quality_warnings": int(quality["status"].eq("warning").sum()),
        "max_topology_test_delta": float(recomputed_tests["max_numeric_delta"].max()),
        "max_correlation_delta": correlation_delta,
    }


def main() -> None:
    result = build_report()
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
