#!/usr/bin/env python3
"""Build the canonical x1-x100 success-only dataset.

The builder accepts generic candidate and profile tables, normalizes them into
one schema, and writes a self-contained dataset directory for downstream
analysis and figures. Source-specific labels and paths are not propagated into
the public tables.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageOps
from scipy.stats import kruskal


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = ROOT / "outputs" / "x1_x100_dataset"

CANDIDATE_PATTERNS = [
    "outputs/*x26_x30_with_figures/consolidated_candidate_profiles.csv",
    "outputs/*/official_visual_workbench/*candidate_table.csv",
]
PROFILE_PATTERNS = [
    "outputs/*x26_x30_with_figures/consolidated_activity_profiles.csv",
    "outputs/*/official_visual_workbench/*integrated_activity_table.csv",
]
STRUCTURE_PATTERN = "outputs/*/official_visual_workbench/assets/images/*_mace_relaxed_structure.png"

TOPOLOGY_ORDER = ["bridged", "independent adjacent", "independent separated"]
ACTIVITY_ORDER = [
    "Catalase",
    "Oxidase",
    "Peroxidase",
    "Glutathione Peroxidase",
    "Glucose Oxidase",
    "DNase",
]
METHOD_LABELS = {
    "first_pass": "First pass",
    "gfn1_extended": "GFN1 extended",
    "gfn1_scf_fallback": "GFN1 SCF fallback",
    "gfn2_deep": "GFN2 deep",
    "gfn2_extended": "GFN2 extended",
}

VARIANT_RE = re.compile(
    r"(?:^|_)(?P<kind>bridged|adjacent|separated)"
    r"(?:_x(?P<index>\d+))?_d(?P<distance>\d+)_"
    r"(?P<doping>ns|n|s)_a(?P<angle>\d+)$",
    re.IGNORECASE,
)
BASE_VARIANT_INDEX = {
    ("bridged", 5.6, "NS", 0): 1,
    ("bridged", 6.0, "NS", 30): 2,
    ("bridged", 6.4, "N", 60): 3,
    ("adjacent", 6.8, "N", 0): 1,
    ("adjacent", 7.6, "NS", 30): 2,
    ("adjacent", 8.4, "S", 60): 3,
    ("separated", 9.2, "N", 0): 1,
    ("separated", 10.8, "NS", 30): 2,
    ("separated", 12.4, "S", 60): 3,
}


def _absolute(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def _discover_one(pattern: str) -> Path:
    matches = sorted(ROOT.glob(pattern), key=lambda path: path.stat().st_mtime, reverse=True)
    if not matches:
        raise FileNotFoundError(f"No input matched {pattern}")
    return matches[0]


def discover_inputs() -> tuple[list[Path], list[Path]]:
    candidates = [_discover_one(pattern) for pattern in CANDIDATE_PATTERNS]
    profiles = [_discover_one(pattern) for pattern in PROFILE_PATTERNS]
    if len(set(candidates)) != len(candidates) or len(set(profiles)) != len(profiles):
        raise RuntimeError("Input discovery returned duplicate tables")
    return candidates, profiles


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    default_candidates, default_profiles = discover_inputs()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument(
        "--candidate-table",
        action="append",
        default=None,
        help="Candidate table; repeat for multiple sources",
    )
    parser.add_argument(
        "--profile-table",
        action="append",
        default=None,
        help="Profile table; repeat for multiple sources",
    )
    parser.set_defaults(
        default_candidate_tables=[str(path) for path in default_candidates],
        default_profile_tables=[str(path) for path in default_profiles],
    )
    return parser.parse_args(argv)


def _bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _successful_candidate_rows(frame: pd.DataFrame) -> pd.DataFrame:
    rows = frame.copy()
    if "calculable" in rows.columns:
        success = _bool_series(rows["calculable"])
    elif "candidate_status" in rows.columns:
        success = rows["candidate_status"].astype(str).str.strip().str.lower().eq("calculable")
    elif "complete_activity_profiles" in rows.columns:
        success = _numeric(rows["complete_activity_profiles"]).fillna(0).gt(0)
    else:
        raise ValueError(
            "Candidate table has no verifiable success field; expected calculable, "
            "candidate_status, or complete_activity_profiles"
        )
    return rows[success].copy()


def _successful_profile_rows(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"profile_status", "frame_count", "converged_frames"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Profile table is missing success fields: {', '.join(missing)}")

    frame_count = _numeric(frame["frame_count"])
    converged_frames = _numeric(frame["converged_frames"])
    success = (
        frame["profile_status"].astype(str).str.strip().str.lower().eq("success")
        & frame_count.gt(0)
        & frame_count.eq(converged_frames)
    )
    if "calculable" in frame.columns:
        success &= _bool_series(frame["calculable"])
    if "mace_status" in frame.columns:
        success &= frame["mace_status"].astype(str).str.strip().str.lower().eq("success")
    if "status" in frame.columns:
        success &= frame["status"].astype(str).str.strip().str.lower().eq("success")
    if "error" in frame.columns:
        errors = frame["error"]
        success &= errors.isna() | errors.astype(str).str.strip().eq("")
    return frame[success].copy()


def topology_label(value: Any) -> str:
    key = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    mapping = {
        "bridged": "bridged",
        "adjacent": "independent adjacent",
        "independent_adjacent": "independent adjacent",
        "separated": "independent separated",
        "independent_separated": "independent separated",
    }
    if key not in mapping:
        raise ValueError(f"Unknown topology value: {value}")
    return mapping[key]


def parse_variant(value: Any) -> dict[str, Any]:
    text = str(value).strip()
    match = VARIANT_RE.search(text)
    if match is None:
        raise ValueError(f"Cannot parse variant: {text}")
    kind = match.group("kind").lower()
    distance = int(match.group("distance")) / 10.0
    doping = match.group("doping").upper()
    angle = int(match.group("angle"))
    index_text = match.group("index")
    index = int(index_text) if index_text else BASE_VARIANT_INDEX.get((kind, distance, doping, angle))
    if index is None:
        raise ValueError(f"Cannot infer design index for variant: {text}")
    return {
        "design_id": f"x{int(index)}",
        "design_index": int(index),
        "topology": topology_label(kind),
        "distance_a": float(distance),
        "doping": doping,
        "angle_deg": int(angle),
        "variant_id": f"x{int(index)}_{kind}_d{distance:.1f}_{doping.lower()}_a{angle}",
    }


def _with_variant_fields(frame: pd.DataFrame) -> pd.DataFrame:
    parsed = pd.DataFrame([parse_variant(value) for value in frame["variant_id"]], index=frame.index)
    base = frame.drop(columns=[column for column in parsed.columns if column in frame.columns]).copy()
    return pd.concat([base, parsed], axis=1)


def _read_tables(paths: list[Path]) -> tuple[list[pd.DataFrame], list[dict[str, Any]]]:
    frames: list[pd.DataFrame] = []
    audit: list[dict[str, Any]] = []
    for index, path in enumerate(paths, start=1):
        frame = pd.read_csv(path)
        frames.append(frame)
        audit.append(
            {
                "input_id": f"input_{index}",
                "rows": int(len(frame)),
                "columns": int(len(frame.columns)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    return frames, audit


def _candidate_metrics(profile_rows: pd.DataFrame) -> pd.DataFrame:
    rows = profile_rows.copy()
    rows["score"] = _numeric(rows["score_total"])
    rows["max_force_ev_per_a"] = _numeric(rows["relaxed_max_force_ev_per_a"])
    return rows.groupby("candidate_id", as_index=False).agg(
        score=("score", "first"),
        max_force_ev_per_a=("max_force_ev_per_a", "first"),
    )


def normalize_tables(
    candidate_frames: list[pd.DataFrame],
    profile_frames: list[pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_profiles = pd.concat(profile_frames, ignore_index=True, sort=False)
    complete_private = _successful_profile_rows(raw_profiles)
    metrics = _candidate_metrics(complete_private)

    candidate_parts: list[pd.DataFrame] = []
    for frame in candidate_frames:
        rows = _successful_candidate_rows(frame)
        rows = _with_variant_fields(rows)
        rows["candidate_id"] = rows["candidate_id"].astype(str)
        rows["activity_pair"] = rows["combo"].astype(str)
        rows = rows.merge(metrics, on="candidate_id", how="left")
        if "score_total" in rows.columns:
            rows["score"] = _numeric(rows["score_total"]).fillna(rows["score"])
        if "relaxed_max_force_ev_per_a" in rows.columns:
            rows["max_force_ev_per_a"] = _numeric(rows["relaxed_max_force_ev_per_a"]).fillna(
                rows["max_force_ev_per_a"]
            )
        candidate_parts.append(rows)
    candidate_private = pd.concat(candidate_parts, ignore_index=True, sort=False)
    candidate_private = candidate_private.sort_values(
        ["design_index", "topology", "activity_pair", "candidate_id"]
    ).reset_index(drop=True)
    candidate_columns = [
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
    candidates = candidate_private[candidate_columns].copy()

    profile_private = _with_variant_fields(complete_private)
    profile_private["candidate_id"] = profile_private["candidate_id"].astype(str)
    profile_private["activity_pair"] = profile_private["combo"].astype(str)
    profile_private["selected_method"] = profile_private["selected_source"].map(METHOD_LABELS)
    profile_private["score"] = _numeric(profile_private["score_total"])
    profile_private["max_force_ev_per_a"] = _numeric(profile_private["relaxed_max_force_ev_per_a"])
    for column in [
        "adsorption_energy_ev",
        "activation_metric_ev",
        "reaction_energy_ev",
        "scan_energy_range_ev",
        "frame_count",
        "converged_frames",
    ]:
        profile_private[column] = _numeric(profile_private[column])
    profile_private["profile_id"] = (
        profile_private["candidate_id"]
        + "_"
        + profile_private["activity"].str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True).str.strip("_")
    )
    profile_private = profile_private.sort_values(
        ["design_index", "candidate_id", "activity"]
    ).reset_index(drop=True)
    profile_columns = [
        "profile_id",
        "candidate_id",
        "design_id",
        "design_index",
        "activity_pair",
        "activity",
        "partner_activity",
        "topology",
        "distance_a",
        "doping",
        "angle_deg",
        "variant_id",
        "selected_method",
        "adsorption_energy_ev",
        "activation_metric_ev",
        "reaction_energy_ev",
        "scan_energy_range_ev",
        "frame_count",
        "converged_frames",
        "score",
        "max_force_ev_per_a",
    ]
    profiles = profile_private[profile_columns].copy()
    return candidates, profiles, profile_private


def build_designs(candidates: pd.DataFrame, profiles: pd.DataFrame) -> pd.DataFrame:
    candidate_counts = candidates.groupby("design_index", as_index=False).agg(
        candidate_count=("candidate_id", "size"),
        activity_pair_count=("activity_pair", "nunique"),
        topology_count=("topology", "nunique"),
    )
    profile_counts = profiles.groupby("design_index", as_index=False).agg(
        profile_count=("profile_id", "size"),
        activity_count=("activity", "nunique"),
    )
    designs = pd.DataFrame({"design_index": np.arange(1, 101, dtype=int)})
    designs["design_id"] = designs["design_index"].map(lambda value: f"x{value}")
    designs = designs.merge(candidate_counts, on="design_index", how="left")
    designs = designs.merge(profile_counts, on="design_index", how="left")
    for column in [
        "candidate_count",
        "activity_pair_count",
        "topology_count",
        "profile_count",
        "activity_count",
    ]:
        designs[column] = designs[column].fillna(0).astype(int)
    return designs


def topology_tests(profiles: pd.DataFrame) -> pd.DataFrame:
    metric_columns = {
        "adsorption": "adsorption_energy_ev",
        "activation": "activation_metric_ev",
    }
    rows: list[dict[str, Any]] = []
    for activity in ACTIVITY_ORDER:
        activity_rows = profiles[profiles["activity"].eq(activity)]
        for metric, column in metric_columns.items():
            groups: list[np.ndarray] = []
            names: list[str] = []
            sizes: list[int] = []
            for topology in TOPOLOGY_ORDER:
                values = activity_rows.loc[activity_rows["topology"].eq(topology), column].dropna().to_numpy(dtype=float)
                sizes.append(len(values))
                if len(values) >= 2:
                    groups.append(values)
                    names.append(topology)
            if len(groups) < 2:
                continue
            statistic, p_value = kruskal(*groups)
            total_n = sum(len(group) for group in groups)
            group_count = len(groups)
            epsilon_squared = max(0.0, (float(statistic) - group_count + 1) / max(total_n - group_count, 1))
            rows.append(
                {
                    "activity": activity,
                    "metric": metric,
                    "compared_topologies": ";".join(names),
                    "group_sizes": ";".join(str(value) for value in sizes),
                    "statistic": float(statistic),
                    "p_value": float(p_value),
                    "epsilon_squared": float(epsilon_squared),
                }
            )
    rows.sort(key=lambda row: row["p_value"])
    count = len(rows)
    running = 1.0
    q_values = [1.0] * count
    for reverse_index in range(count - 1, -1, -1):
        rank = reverse_index + 1
        running = min(running, rows[reverse_index]["p_value"] * count / rank)
        q_values[reverse_index] = min(1.0, running)
    for row, q_value in zip(rows, q_values):
        row["q_value"] = float(q_value)
        row["retained_at_fdr_0_05"] = bool(q_value < 0.05)
    return pd.DataFrame(rows)


def build_summaries(
    candidates: pd.DataFrame,
    profiles: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    candidate_design_topology = (
        candidates.groupby(["design_index", "topology"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    activity_pair_topology = (
        candidates.groupby(["activity_pair", "topology"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    profile_activity_topology = (
        profiles.groupby(["activity", "topology"], as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    method_activity = (
        profiles.groupby(["activity", "selected_method"], as_index=False)
        .size()
        .rename(columns={"size": "profile_count"})
    )
    geometry = (
        candidates.groupby(["design_index", "topology", "distance_a", "doping", "angle_deg"], as_index=False)
        .size()
        .rename(columns={"size": "candidate_count"})
    )
    descriptor_summary = profiles.groupby("activity", as_index=False).agg(
        n=("profile_id", "size"),
        adsorption_q25=("adsorption_energy_ev", lambda values: float(np.percentile(values, 25))),
        adsorption_median=("adsorption_energy_ev", "median"),
        adsorption_q75=("adsorption_energy_ev", lambda values: float(np.percentile(values, 75))),
        activation_q25=("activation_metric_ev", lambda values: float(np.percentile(values, 25))),
        activation_median=("activation_metric_ev", "median"),
        activation_q75=("activation_metric_ev", lambda values: float(np.percentile(values, 75))),
        zero_activation_count=("activation_metric_ev", lambda values: int(np.sum(np.asarray(values) == 0))),
    )
    numeric = [
        "distance_a",
        "adsorption_energy_ev",
        "activation_metric_ev",
        "reaction_energy_ev",
        "scan_energy_range_ev",
        "score",
        "max_force_ev_per_a",
    ]
    correlation = profiles[numeric].corr(method="spearman").round(6)
    correlation.index.name = "metric"
    correlation = correlation.reset_index()
    return {
        "candidate_design_topology": candidate_design_topology,
        "activity_pair_topology": activity_pair_topology,
        "profile_activity_topology": profile_activity_topology,
        "method_activity": method_activity,
        "geometry": geometry,
        "descriptor_summary": descriptor_summary,
        "descriptor_correlation": correlation,
    }


def _crop_structure(path: Path) -> Image.Image:
    image = Image.open(path).convert("RGB")
    image = image.crop((0, int(image.height * 0.105), image.width, image.height))
    inverted = ImageOps.invert(image.convert("L"))
    bbox = inverted.point(lambda value: 255 if value > 10 else 0).getbbox()
    if bbox:
        left, top, right, bottom = bbox
        margin = 24
        image = image.crop(
            (
                max(0, left - margin),
                max(0, top - margin),
                min(image.width, right + margin),
                min(image.height, bottom + margin),
            )
        )
    return image


def _load_scan(row: pd.Series) -> pd.DataFrame:
    result_path = Path(str(row["selected_result_path"]))
    if not result_path.is_absolute():
        result_path = Path(str(row["run_dir"])) / result_path
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    values = payload.get("reaction_profile", {}).get("relative_energies_ev") or []
    if not values:
        values = payload.get("mechanism_visualization", {}).get("relative_energies_ev") or []
    energies = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if len(energies) == 0:
        raise RuntimeError(f"No scan energies for candidate {row['candidate_id']}")
    return pd.DataFrame(
        {
            "candidate_id": str(row["candidate_id"]),
            "activity": str(row["activity"]),
            "scan_step": np.arange(1, len(energies) + 1),
            "relative_energy_ev": energies,
        }
    )


def build_representative_assets(
    out_dir: Path,
    candidates: pd.DataFrame,
    profiles: pd.DataFrame,
    profile_private: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    images = sorted(ROOT.glob(STRUCTURE_PATTERN))
    image_by_candidate = {path.name.split("_", 1)[0]: path for path in images}
    profile_counts = profiles.groupby("candidate_id").size()
    eligible = candidates[
        candidates["candidate_id"].isin(image_by_candidate)
        & candidates["candidate_id"].map(profile_counts).fillna(0).ge(2)
    ].copy()
    if eligible.empty:
        raise RuntimeError("No successful candidate has both a structure image and two profiles")
    median_score = float(eligible["score"].median())
    eligible["selection_distance"] = (eligible["score"] - median_score).abs()
    representative = eligible.sort_values(["selection_distance", "design_index", "candidate_id"]).iloc[[0]].drop(
        columns="selection_distance"
    )
    candidate_id = str(representative.iloc[0]["candidate_id"])

    asset_dir = out_dir / "assets"
    asset_dir.mkdir(parents=True, exist_ok=True)
    _crop_structure(image_by_candidate[candidate_id]).save(asset_dir / "representative_structure.png")

    scan_rows = profile_private[profile_private["candidate_id"].astype(str).eq(candidate_id)]
    scans = pd.concat([_load_scan(row) for _, row in scan_rows.iterrows()], ignore_index=True)
    return representative, scans


def validate_dataset(
    candidates: pd.DataFrame,
    profiles: pd.DataFrame,
    designs: pd.DataFrame,
    tests: pd.DataFrame,
) -> dict[str, Any]:
    checks = {
        "candidate_rows": int(len(candidates)),
        "candidate_ids_unique": bool(candidates["candidate_id"].is_unique),
        "profile_rows": int(len(profiles)),
        "profile_ids_unique": bool(profiles["profile_id"].is_unique),
        "all_profile_frames_converged": bool(profiles["frame_count"].eq(profiles["converged_frames"]).all()),
        "profile_candidates_exist": bool(profiles["candidate_id"].isin(candidates["candidate_id"]).all()),
        "candidate_design_ids_match": bool(
            candidates["design_id"].eq(candidates["design_index"].map(lambda value: f"x{value}")).all()
        ),
        "profile_design_ids_match": bool(
            profiles["design_id"].eq(profiles["design_index"].map(lambda value: f"x{value}")).all()
        ),
        "design_ids_sequential": designs["design_id"].tolist() == [f"x{value}" for value in range(1, 101)],
        "design_rows": int(len(designs)),
        "design_index_min": int(designs["design_index"].min()),
        "design_index_max": int(designs["design_index"].max()),
        "topologies": sorted(candidates["topology"].unique().tolist()),
        "activities": sorted(profiles["activity"].unique().tolist()),
        "topology_tests": int(len(tests)),
        "retained_tests": int(tests["retained_at_fdr_0_05"].sum()),
    }
    expected = {
        "candidate_rows": 355,
        "profile_rows": 699,
        "design_rows": 100,
        "design_index_min": 1,
        "design_index_max": 100,
        "topology_tests": 10,
        "retained_tests": 8,
    }
    for key, value in expected.items():
        if checks[key] != value:
            raise RuntimeError(f"Dataset check failed for {key}: {checks[key]} != {value}")
    for key in [
        "candidate_ids_unique",
        "profile_ids_unique",
        "all_profile_frames_converged",
        "profile_candidates_exist",
        "candidate_design_ids_match",
        "profile_design_ids_match",
        "design_ids_sequential",
    ]:
        if not checks[key]:
            raise RuntimeError(f"Dataset check failed for {key}")
    if checks["topologies"] != TOPOLOGY_ORDER:
        raise RuntimeError(f"Unexpected topology vocabulary: {checks['topologies']}")
    return checks


def reset_output(out_dir: Path) -> None:
    target = out_dir.resolve()
    if target == ROOT.resolve() or target == Path(target.anchor):
        raise RuntimeError(f"Refusing to clear unsafe output path: {target}")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)


def build_dataset(
    candidate_paths: list[Path],
    profile_paths: list[Path],
    out_dir: Path,
) -> dict[str, Any]:
    reset_output(out_dir)
    candidate_frames, candidate_audit = _read_tables(candidate_paths)
    profile_frames, profile_audit = _read_tables(profile_paths)
    candidates, profiles, profile_private = normalize_tables(candidate_frames, profile_frames)
    designs = build_designs(candidates, profiles)
    tests = topology_tests(profiles)
    summaries = build_summaries(candidates, profiles)
    representative, scans = build_representative_assets(out_dir, candidates, profiles, profile_private)
    checks = validate_dataset(candidates, profiles, designs, tests)

    candidates.to_csv(out_dir / "candidates.csv", index=False)
    profiles.to_csv(out_dir / "profiles.csv", index=False)
    designs.to_csv(out_dir / "designs.csv", index=False)
    tests.to_csv(out_dir / "topology_tests.csv", index=False)
    representative.to_csv(out_dir / "representative.csv", index=False)
    scans.to_csv(out_dir / "representative_scans.csv", index=False)
    for name, frame in summaries.items():
        frame.to_csv(out_dir / f"{name}.csv", index=False)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "x1-x100",
        "selection": {
            "candidates": "calculable records only",
            "profiles": "complete records with all recorded frames converged",
        },
        "counts": {
            "candidates": int(len(candidates)),
            "profiles": int(len(profiles)),
            "frames": int(profiles["converged_frames"].sum()),
            "design_indices_with_candidates": int(candidates["design_index"].nunique()),
            "activity_pairs": int(candidates["activity_pair"].nunique()),
            "activities": int(profiles["activity"].nunique()),
            "topologies": int(candidates["topology"].nunique()),
        },
        "inputs": {
            "candidate_tables": candidate_audit,
            "profile_tables": profile_audit,
        },
        "checks": checks,
        "files": sorted(path.name for path in out_dir.iterdir() if path.is_file()),
    }
    (out_dir / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    readme = f"""# x1-x100 Dataset

This directory is the only data interface used by the main figure builder.

- Candidates: {len(candidates)}
- Complete profiles: {len(profiles)}
- Converged frames in complete profiles: {int(profiles['converged_frames'].sum())}
- Design indices with at least one candidate: {candidates['design_index'].nunique()}
- Activity pairs: {candidates['activity_pair'].nunique()}
- Activities: {profiles['activity'].nunique()}
- Topologies: {candidates['topology'].nunique()}

`candidates.csv`, `profiles.csv`, and `designs.csv` are the canonical record tables.
All remaining CSV files are deterministic summaries derived from those tables.
"""
    (out_dir / "README.md").write_text(readme, encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    candidate_values = args.candidate_table or args.default_candidate_tables
    profile_values = args.profile_table or args.default_profile_tables
    manifest = build_dataset(
        [_absolute(path) for path in candidate_values],
        [_absolute(path) for path in profile_values],
        _absolute(args.out_dir),
    )
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
