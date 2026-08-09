#!/usr/bin/env python3
"""Verify the self-contained E2N publication release.

The script uses only the Python standard library and resolves every release
path from its own location, so it can be invoked from any working directory.
It reports all detected failures before returning a non-zero exit status.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable
from urllib.parse import unquote, urlsplit


PUBLICATION_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PUBLICATION_ROOT / "data" / "x1_x100_dataset"
FIGURE_DIR = PUBLICATION_ROOT / "figures"
SOURCE_DATA_DIR = FIGURE_DIR / "source_data"
QA_DIR = FIGURE_DIR / "qa"
REPORT_DIR = PUBLICATION_ROOT / "reports"

EXPECTED_CANDIDATE_COUNT = 355
EXPECTED_PROFILE_COUNT = 699
EXPECTED_FRAME_COUNT = 3515
EXPECTED_DESIGN_INDEX_COUNT = 38
EXPECTED_ACTIVITY_PAIR_COUNT = 13
EXPECTED_REPRESENTATIVE_ID = "31967c28"

EXPECTED_ACTIVITIES = {
    "Catalase",
    "DNase",
    "Glucose Oxidase",
    "Glutathione Peroxidase",
    "Oxidase",
    "Peroxidase",
}
EXPECTED_TOPOLOGIES = {
    "bridged",
    "independent adjacent",
    "independent separated",
}
EXPECTED_LINKAGE = {2: 346, 1: 7, 0: 2}

EXPECTED_FIGURE_STEMS = {
    "Fig. 1": "fig1_canonical_evidence_trace",
    "Fig. 2": "fig2_geometry_landscape",
    "Fig. 3": "fig3_activity_pair_composition",
    "Fig. 4": "fig4_profile_descriptors",
    "Fig. 5": "fig5_topology_statistics",
}
EXPECTED_FIGURE_SOURCE_FILES = {
    "fig1_canonical_summary.csv",
    "fig1_design_topology.csv",
    "fig1_representative_candidate.csv",
    "fig1_representative_scans.csv",
    "fig2_distance_points.csv",
    "fig2_geometry.csv",
    "fig2_highlighted_geometry.csv",
    "fig3_activity_pair_counts.csv",
    "fig3_activity_pair_design.csv",
    "fig3_activity_pair_topology.csv",
    "fig3_activity_pair_topology_percent.csv",
    "fig4_activity_counts.csv",
    "fig4_descriptor_points.csv",
    "fig4_method_activity.csv",
    "fig4_zero_activation.csv",
    "fig5_topology_medians.csv",
    "fig5_topology_tests.csv",
}
EXPECTED_PANEL_ROWS = {
    ("Fig. 1", "A"): (
        "Representative retained candidate and its canonical metadata.",
        ("fig1_representative_candidate.csv",),
        "../data/x1_x100_dataset/representative_structure.png",
    ),
    ("Fig. 1", "B"): (
        "Complete activity-profile trajectories for the representative candidate.",
        ("fig1_representative_scans.csv",),
        "",
    ),
    ("Fig. 1", "C"): (
        "Candidate occupancy and canonical retained-set totals across x1-x100.",
        ("fig1_design_topology.csv", "fig1_canonical_summary.csv"),
        "",
    ),
    ("Fig. 2", "A"): (
        "Retained candidate-record counts at each requested geometry coordinate.",
        ("fig2_geometry.csv",),
        "",
    ),
    ("Fig. 2", "B"): (
        "Design-index composition of the highlighted requested geometry.",
        ("fig2_highlighted_geometry.csv",),
        "",
    ),
    ("Fig. 2", "C"): (
        "Individual requested candidate distances by topology.",
        ("fig2_distance_points.csv",),
        "",
    ),
    ("Fig. 3", "A"): (
        "Retained candidate-record counts by declared activity pair.",
        ("fig3_activity_pair_counts.csv",),
        "",
    ),
    ("Fig. 3", "B"): (
        "Absolute and within-pair topology composition of each activity pair.",
        (
            "fig3_activity_pair_topology.csv",
            "fig3_activity_pair_topology_percent.csv",
        ),
        "",
    ),
    ("Fig. 3", "C"): (
        "Declared activity-pair occupancy by design index.",
        ("fig3_activity_pair_design.csv",),
        "",
    ),
    ("Fig. 4", "A"): (
        "Complete-profile count by activity.",
        ("fig4_activity_counts.csv",),
        "",
    ),
    ("Fig. 4", "B"): (
        "Selected computational-route composition by activity.",
        ("fig4_method_activity.csv",),
        "",
    ),
    ("Fig. 4", "C"): (
        "Individual adsorption descriptor values.",
        ("fig4_descriptor_points.csv",),
        "",
    ),
    ("Fig. 4", "D"): (
        "Individual forward-scan peak descriptors and zero-valued boundary counts.",
        ("fig4_descriptor_points.csv", "fig4_zero_activation.csv"),
        "",
    ),
    ("Fig. 5", "A"): (
        "Benjamini-Hochberg-adjusted q values for topology tests.",
        ("fig5_topology_tests.csv",),
        "",
    ),
    ("Fig. 5", "B"): (
        "Epsilon-squared effect sizes for topology tests.",
        ("fig5_topology_tests.csv",),
        "",
    ),
    ("Fig. 5", "C"): (
        "Topology medians for retained associations.",
        ("fig5_topology_medians.csv",),
        "",
    ),
}
EXPECTED_REPORT_TABLES = {
    "activity_pair_counts.csv",
    "activity_pair_topology_summary.csv",
    "candidate_angle_summary.csv",
    "candidate_doping_summary.csv",
    "candidate_profile_linkage.csv",
    "candidate_topology_summary.csv",
    "data_quality_checks.csv",
    "declared_target_profile_reconciliation.csv",
    "design_occupancy.csv",
    "geometry_cell_summary.csv",
    "method_activity_summary.csv",
    "method_total_summary.csv",
    "profile_activity_summary.csv",
    "profile_linkage_summary.csv",
    "topology_tests_recomputed.csv",
    "version_reconciliation.csv",
}
EXPECTED_DATA_FILES = {
    "README.md",
    "activity_pair_topology.csv",
    "candidate_design_topology.csv",
    "candidates.csv",
    "dataset_manifest.json",
    "descriptor_correlation.csv",
    "descriptor_summary.csv",
    "designs.csv",
    "geometry.csv",
    "method_activity.csv",
    "profile_activity_topology.csv",
    "profiles.csv",
    "representative.csv",
    "representative_scans.csv",
    "representative_structure.png",
    "topology_tests.csv",
}
MANIFEST_RELATIVE_PATH = "RELEASE_MANIFEST.json"
EXPECTED_RELEASE_FILES = frozenset(
    {
        "CITATION.cff",
        "EVIDENCE_CONTRACT.md",
        "LICENSE_SCOPE.md",
        "README.md",
        MANIFEST_RELATIVE_PATH,
        "requirements-release.txt",
        "figures/README.md",
        "figures/figure_captions_v6.md",
        "figures/manifest_v6.json",
        "figures/panel_source_data_index.csv",
        "figures/qa/submission_qa.csv",
        "figures/qa/success_only_source_data_audit.json",
        "reports/e2n_x1_x100_data_summary_zh.md",
        "scripts/build_release_manifest.py",
        "scripts/verify_publication_release.py",
        *(f"data/x1_x100_dataset/{name}" for name in EXPECTED_DATA_FILES),
        *(
            f"figures/{stem}.{extension}"
            for stem in EXPECTED_FIGURE_STEMS.values()
            for extension in ("pdf", "png", "svg")
        ),
        *(f"figures/source_data/{name}" for name in EXPECTED_FIGURE_SOURCE_FILES),
        *(f"reports/tables/{name}" for name in EXPECTED_REPORT_TABLES),
    }
)
EXPECTED_MANIFEST_EVIDENCE_COUNTS = {
    "retained_canonical_candidate_records": EXPECTED_CANDIDATE_COUNT,
    "complete_activity_specific_profiles": EXPECTED_PROFILE_COUNT,
    "converged_scan_frames_within_complete_profiles": EXPECTED_FRAME_COUNT,
}
TOPOLOGY_ORDER = (
    "bridged",
    "independent adjacent",
    "independent separated",
)

if len(EXPECTED_RELEASE_FILES) != 79:
    raise RuntimeError(
        f"release allowlist has {len(EXPECTED_RELEASE_FILES)} paths; expected 79"
    )

CANDIDATE_PROFILE_METADATA_FIELDS = (
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
)
NUMERIC_METADATA_FIELDS = {
    "design_index",
    "distance_a",
    "angle_deg",
    "score",
    "max_force_ev_per_a",
}

MARKDOWN_LINK_PATTERN = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


@dataclass(frozen=True)
class CheckResult:
    passed: bool
    name: str
    detail: str


class Audit:
    def __init__(self) -> None:
        self.results: list[CheckResult] = []

    def check(self, name: str, condition: bool, detail: str) -> None:
        self.results.append(CheckResult(bool(condition), name, detail))

    def run(self, name: str, callback: Callable[[], None]) -> None:
        try:
            callback()
        except Exception as exc:  # Keep auditing independent release sections.
            self.check(name, False, f"unexpected {type(exc).__name__}: {exc}")

    @property
    def failure_count(self) -> int:
        return sum(not result.passed for result in self.results)

    def print_report(self) -> None:
        print(f"E2N publication release verification: {PUBLICATION_ROOT}")
        for result in self.results:
            label = "PASS" if result.passed else "FAIL"
            print(f"[{label}] {result.name}: {result.detail}")
        passed = len(self.results) - self.failure_count
        print(
            f"Summary: {passed} passed, {self.failure_count} failed, "
            f"{len(self.results)} total checks."
        )


def read_csv_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"CSV has duplicate headers: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    return read_csv_table(path)[1]


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def parse_int(value: str, field: str) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid integer in {field}: {value!r}") from exc
    if not math.isfinite(number) or not number.is_integer():
        raise ValueError(f"invalid integer in {field}: {value!r}")
    return int(number)


def parse_float(value: str, field: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid number in {field}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"non-finite number in {field}: {value!r}")
    return number


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().casefold()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid Boolean value: {value!r}")


def format_sample(values: Iterable[Any], limit: int = 5) -> str:
    items = [str(value) for value in values]
    if len(items) <= limit:
        return ", ".join(items)
    return f"{', '.join(items[:limit])}, ... ({len(items)} total)"


def duplicate_values(values: Iterable[str]) -> list[str]:
    return sorted(value for value, count in Counter(values).items() if count > 1)


def metadata_values_equal(field: str, left: str, right: str) -> bool:
    if field not in NUMERIC_METADATA_FIELDS:
        return left == right
    try:
        return math.isclose(
            parse_float(left, field),
            parse_float(right, field),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    except ValueError:
        return False


def require_dataset_files(audit: Audit) -> bool:
    required = {
        "candidates.csv",
        "profiles.csv",
        "representative.csv",
        "representative_scans.csv",
        "topology_tests.csv",
    }
    missing = sorted(
        name
        for name in required
        if not (DATA_DIR / name).is_file() or (DATA_DIR / name).stat().st_size == 0
    )
    audit.check(
        "canonical dataset files",
        not missing,
        "all five controlling tables are present"
        if not missing
        else f"missing or empty: {format_sample(missing)}",
    )
    return not missing


def audit_canonical_dataset(audit: Audit) -> None:
    if not require_dataset_files(audit):
        return

    candidates = read_csv_rows(DATA_DIR / "candidates.csv")
    profiles = read_csv_rows(DATA_DIR / "profiles.csv")
    tests = read_csv_rows(DATA_DIR / "topology_tests.csv")
    representative_rows = read_csv_rows(DATA_DIR / "representative.csv")
    scans = read_csv_rows(DATA_DIR / "representative_scans.csv")

    candidate_ids = [row.get("candidate_id", "").strip() for row in candidates]
    profile_ids = [row.get("profile_id", "").strip() for row in profiles]
    candidate_duplicates = duplicate_values(candidate_ids)
    profile_duplicates = duplicate_values(profile_ids)

    audit.check(
        "canonical row counts",
        len(candidates) == EXPECTED_CANDIDATE_COUNT
        and len(profiles) == EXPECTED_PROFILE_COUNT,
        f"candidates={len(candidates)} (expected 355); "
        f"profiles={len(profiles)} (expected 699)",
    )
    audit.check(
        "candidate primary keys",
        all(candidate_ids) and not candidate_duplicates,
        "355 non-empty unique candidate IDs"
        if all(candidate_ids) and not candidate_duplicates
        else f"empty={candidate_ids.count('')}; duplicates={format_sample(candidate_duplicates)}",
    )
    audit.check(
        "profile primary keys",
        all(profile_ids) and not profile_duplicates,
        "699 non-empty unique profile IDs"
        if all(profile_ids) and not profile_duplicates
        else f"empty={profile_ids.count('')}; duplicates={format_sample(profile_duplicates)}",
    )

    candidate_by_id = {
        row["candidate_id"].strip(): row
        for row in candidates
        if row.get("candidate_id", "").strip()
    }
    missing_foreign_keys = sorted(
        {
            row.get("candidate_id", "").strip()
            for row in profiles
            if row.get("candidate_id", "").strip() not in candidate_by_id
        }
    )
    audit.check(
        "profile candidate foreign keys",
        not missing_foreign_keys,
        "all 699 profiles link to canonical candidates"
        if not missing_foreign_keys
        else f"missing candidate IDs: {format_sample(missing_foreign_keys)}",
    )

    metadata_mismatches: list[str] = []
    activity_pair_mismatches: list[str] = []
    for profile in profiles:
        candidate = candidate_by_id.get(profile.get("candidate_id", "").strip())
        if candidate is None:
            continue
        for field in CANDIDATE_PROFILE_METADATA_FIELDS:
            if not metadata_values_equal(
                field,
                candidate.get(field, ""),
                profile.get(field, ""),
            ):
                metadata_mismatches.append(
                    f"{profile.get('profile_id', '?')}:{field}"
                )
        pair_activities = {
            value.strip()
            for value in candidate.get("activity_pair", "").split(" + ")
            if value.strip()
        }
        activity = profile.get("activity", "").strip()
        partner = profile.get("partner_activity", "").strip()
        if (
            len(pair_activities) != 2
            or {activity, partner} != pair_activities
            or activity == partner
        ):
            activity_pair_mismatches.append(profile.get("profile_id", "?"))
    audit.check(
        "linked candidate/profile metadata",
        not metadata_mismatches and not activity_pair_mismatches,
        "candidate metadata and declared activity pairs agree for every profile"
        if not metadata_mismatches and not activity_pair_mismatches
        else (
            f"field mismatches={format_sample(metadata_mismatches)}; "
            f"activity-pair mismatches={format_sample(activity_pair_mismatches)}"
        ),
    )

    invalid_frames: list[str] = []
    frame_total = 0
    converged_total = 0
    for profile in profiles:
        try:
            frame_count = parse_int(profile.get("frame_count", ""), "frame_count")
            converged = parse_int(
                profile.get("converged_frames", ""), "converged_frames"
            )
        except ValueError:
            invalid_frames.append(profile.get("profile_id", "?"))
            continue
        frame_total += frame_count
        converged_total += converged
        if frame_count <= 0 or converged <= 0 or frame_count != converged:
            invalid_frames.append(profile.get("profile_id", "?"))
    audit.check(
        "complete-profile frame accounting",
        not invalid_frames
        and frame_total == EXPECTED_FRAME_COUNT
        and converged_total == EXPECTED_FRAME_COUNT,
        f"frame_count sum={frame_total}; converged sum={converged_total}; "
        f"invalid profiles={format_sample(invalid_frames) if invalid_frames else 'none'}",
    )

    try:
        design_indices = {parse_int(row["design_index"], "design_index") for row in candidates}
    except (KeyError, ValueError) as exc:
        design_indices = set()
        design_error = str(exc)
    else:
        design_error = ""
    activity_pairs = {row.get("activity_pair", "").strip() for row in candidates}
    activities = {row.get("activity", "").strip() for row in profiles}
    candidate_topologies = {row.get("topology", "").strip() for row in candidates}
    profile_topologies = {row.get("topology", "").strip() for row in profiles}
    audit.check(
        "represented design indices",
        len(design_indices) == EXPECTED_DESIGN_INDEX_COUNT
        and all(1 <= value <= 100 for value in design_indices),
        f"represented={len(design_indices)} (expected 38)"
        + (f"; parse error={design_error}" if design_error else ""),
    )
    audit.check(
        "declared activity pairs",
        len(activity_pairs) == EXPECTED_ACTIVITY_PAIR_COUNT and "" not in activity_pairs,
        f"represented={len(activity_pairs)} (expected 13)",
    )
    audit.check(
        "profile activities",
        activities == EXPECTED_ACTIVITIES,
        f"observed={sorted(activities)!r}; expected six canonical activities",
    )
    audit.check(
        "candidate and profile topologies",
        candidate_topologies == EXPECTED_TOPOLOGIES
        and profile_topologies == EXPECTED_TOPOLOGIES,
        f"candidate={sorted(candidate_topologies)!r}; profile={sorted(profile_topologies)!r}",
    )

    linkage_counts = Counter({candidate_id: 0 for candidate_id in candidate_ids})
    for profile in profiles:
        candidate_id = profile.get("candidate_id", "").strip()
        if candidate_id in linkage_counts:
            linkage_counts[candidate_id] += 1
    linkage_distribution = Counter(linkage_counts.values())
    observed_linkage = {
        profile_count: linkage_distribution.get(profile_count, 0)
        for profile_count in EXPECTED_LINKAGE
    }
    unexpected_profile_counts = sorted(
        count for count in linkage_distribution if count not in EXPECTED_LINKAGE
    )
    audit.check(
        "candidate-to-profile linkage",
        observed_linkage == EXPECTED_LINKAGE and not unexpected_profile_counts,
        f"two={observed_linkage[2]}, one={observed_linkage[1]}, "
        f"zero={observed_linkage[0]}; expected 346/7/2"
        + (
            f"; unexpected profile counts={unexpected_profile_counts}"
            if unexpected_profile_counts
            else ""
        ),
    )

    representative_ok = len(representative_rows) == 1
    representative_id = (
        representative_rows[0].get("candidate_id", "").strip()
        if representative_ok
        else ""
    )
    representative_candidate = candidate_by_id.get(representative_id)
    representative_metadata_mismatches: list[str] = []
    if representative_ok and representative_candidate is not None:
        for field in ("candidate_id", *CANDIDATE_PROFILE_METADATA_FIELDS):
            if not metadata_values_equal(
                field,
                representative_candidate.get(field, ""),
                representative_rows[0].get(field, ""),
            ):
                representative_metadata_mismatches.append(field)
    representative_profiles = [
        row
        for row in profiles
        if row.get("candidate_id", "").strip() == representative_id
    ]
    audit.check(
        "representative candidate linkage",
        representative_ok
        and representative_id == EXPECTED_REPRESENTATIVE_ID
        and representative_candidate is not None
        and not representative_metadata_mismatches
        and len(representative_profiles) == 2,
        f"id={representative_id or 'missing'}; canonical_match="
        f"{representative_candidate is not None}; profiles={len(representative_profiles)}; "
        f"metadata_mismatches={format_sample(representative_metadata_mismatches) if representative_metadata_mismatches else 'none'}",
    )

    scan_ids = {row.get("candidate_id", "").strip() for row in scans}
    scan_activities = {row.get("activity", "").strip() for row in scans}
    representative_profile_activities = {
        row.get("activity", "").strip() for row in representative_profiles
    }
    scan_steps: dict[str, list[int]] = defaultdict(list)
    invalid_scan_rows: list[str] = []
    for row_number, row in enumerate(scans, start=2):
        activity = row.get("activity", "").strip()
        try:
            step = parse_int(row.get("scan_step", ""), "scan_step")
            parse_float(row.get("relative_energy_ev", ""), "relative_energy_ev")
        except ValueError:
            invalid_scan_rows.append(str(row_number))
            continue
        scan_steps[activity].append(step)
    step_sequences_ok = all(
        sorted(steps) == [1, 2, 3, 4, 5] for steps in scan_steps.values()
    )
    audit.check(
        "representative scan snapshot",
        len(scans) == 10
        and scan_ids == {representative_id}
        and len(scan_activities) == 2
        and scan_activities == representative_profile_activities
        and not invalid_scan_rows
        and step_sequences_ok,
        f"rows={len(scans)}; candidate_ids={sorted(scan_ids)!r}; "
        f"activities={sorted(scan_activities)!r}; invalid_rows="
        f"{format_sample(invalid_scan_rows) if invalid_scan_rows else 'none'}",
    )

    retained_flags: list[bool] = []
    invalid_test_rows: list[str] = []
    test_keys: list[str] = []
    for row_number, row in enumerate(tests, start=2):
        key = f"{row.get('activity', '')}::{row.get('metric', '')}"
        test_keys.append(key)
        try:
            q_value = parse_float(row.get("q_value", ""), "q_value")
            retained = parse_bool(row.get("retained_at_fdr_0_05", ""))
            parse_float(row.get("statistic", ""), "statistic")
            parse_float(row.get("p_value", ""), "p_value")
            parse_float(row.get("epsilon_squared", ""), "epsilon_squared")
        except ValueError:
            invalid_test_rows.append(str(row_number))
            continue
        retained_flags.append(retained)
        if retained != (q_value < 0.05):
            invalid_test_rows.append(str(row_number))
    audit.check(
        "topology tests",
        len(tests) == 10
        and sum(retained_flags) == 8
        and not duplicate_values(test_keys)
        and not invalid_test_rows,
        f"tests={len(tests)}; retained={sum(retained_flags)}; "
        f"invalid_rows={format_sample(invalid_test_rows) if invalid_test_rows else 'none'}",
    )


def csv_cell_paths(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def audit_panel_index(audit: Audit) -> None:
    index_path = FIGURE_DIR / "panel_source_data_index.csv"
    if not index_path.is_file():
        audit.check(
            "panel source-data index",
            False,
            "missing figures/panel_source_data_index.csv",
        )
        return
    headers, rows = read_csv_table(index_path)
    expected_headers = [
        "figure",
        "panel",
        "claim",
        "source_data",
        "source_data_paths",
        "asset_path",
    ]
    row_by_key: dict[tuple[str, str], dict[str, str]] = {}
    duplicate_keys: list[str] = []
    for row in rows:
        key = (row.get("figure", "").strip(), row.get("panel", "").strip())
        if key in row_by_key:
            duplicate_keys.append("::".join(key))
        row_by_key[key] = row

    mismatches: list[str] = []
    missing_paths: list[str] = []
    mapped_names: set[str] = set()
    for key, (expected_claim, expected_sources, expected_asset) in EXPECTED_PANEL_ROWS.items():
        row = row_by_key.get(key)
        label = "".join(key)
        if row is None:
            mismatches.append(f"{label}:missing row")
            continue
        source_names = tuple(csv_cell_paths(row.get("source_data", "")))
        source_paths = tuple(csv_cell_paths(row.get("source_data_paths", "")))
        expected_paths = tuple(f"source_data/{name}" for name in expected_sources)
        asset_path = row.get("asset_path", "").strip()
        mapped_names.update(source_names)
        if row.get("claim", "").strip() != expected_claim:
            mismatches.append(f"{label}:claim")
        if source_names != expected_sources:
            mismatches.append(f"{label}:source_data")
        if source_paths != expected_paths:
            mismatches.append(f"{label}:source_data_paths")
        if asset_path != expected_asset:
            mismatches.append(f"{label}:asset_path")
        for relative_path in expected_paths:
            if not (FIGURE_DIR / relative_path).is_file():
                missing_paths.append(relative_path)
        if expected_asset and not (FIGURE_DIR / expected_asset).resolve().is_file():
            missing_paths.append(expected_asset)

    expected_keys = set(EXPECTED_PANEL_ROWS)
    observed_keys = set(row_by_key)
    audit.check(
        "panel source-data index",
        headers == expected_headers
        and len(rows) == 16
        and observed_keys == expected_keys
        and not duplicate_keys
        and mapped_names == EXPECTED_FIGURE_SOURCE_FILES
        and not missing_paths
        and not mismatches,
        f"panels={len(rows)}; unique CSVs={len(mapped_names)}; "
        f"headers={'valid' if headers == expected_headers else repr(headers)}; "
        f"missing panels={format_sample(sorted('::'.join(key) for key in expected_keys - observed_keys)) if observed_keys != expected_keys else 'none'}; "
        f"unexpected panels={format_sample(sorted('::'.join(key) for key in observed_keys - expected_keys)) if observed_keys != expected_keys else 'none'}; "
        f"duplicates={format_sample(duplicate_keys) if duplicate_keys else 'none'}; "
        f"unmapped CSVs={format_sample(sorted(EXPECTED_FIGURE_SOURCE_FILES - mapped_names)) if mapped_names != EXPECTED_FIGURE_SOURCE_FILES else 'none'}; "
        f"missing paths={format_sample(missing_paths) if missing_paths else 'none'}; "
        f"row mismatches={format_sample(mismatches) if mismatches else 'none'}",
    )


@dataclass(frozen=True)
class CsvExpectation:
    headers: tuple[str, ...]
    rows: list[dict[str, Any]]
    keys: tuple[str, ...]
    numeric_fields: frozenset[str] = frozenset()
    bool_fields: frozenset[str] = frozenset()


def selected_fields(
    rows: Iterable[dict[str, str]], fields: tuple[str, ...]
) -> list[dict[str, Any]]:
    return [{field: row.get(field, "") for field in fields} for row in rows]


def grouped_counts(
    rows: Iterable[dict[str, str]],
    fields: tuple[str, ...],
    count_field: str,
) -> list[dict[str, Any]]:
    counts = Counter(tuple(row.get(field, "") for field in fields) for row in rows)
    return [
        {**dict(zip(fields, key)), count_field: count}
        for key, count in sorted(counts.items())
    ]


def expected_figure_source_tables() -> dict[str, CsvExpectation]:
    candidates = read_csv_rows(DATA_DIR / "candidates.csv")
    profiles = read_csv_rows(DATA_DIR / "profiles.csv")
    representative = read_csv_rows(DATA_DIR / "representative.csv")
    scans = read_csv_rows(DATA_DIR / "representative_scans.csv")
    tests = read_csv_rows(DATA_DIR / "topology_tests.csv")

    expectations: dict[str, CsvExpectation] = {}

    fig1_design_headers = ("topology", "design_index", "candidate_count")
    expectations["fig1_design_topology.csv"] = CsvExpectation(
        fig1_design_headers,
        grouped_counts(
            candidates, ("topology", "design_index"), "candidate_count"
        ),
        ("topology", "design_index"),
        frozenset({"design_index", "candidate_count"}),
    )

    representative_headers = (
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
    )
    expectations["fig1_representative_candidate.csv"] = CsvExpectation(
        representative_headers,
        selected_fields(representative, representative_headers),
        ("candidate_id",),
        frozenset(
            {
                "design_index",
                "distance_a",
                "angle_deg",
                "score",
                "max_force_ev_per_a",
            }
        ),
    )

    scan_headers = ("candidate_id", "activity", "scan_step", "relative_energy_ev")
    expectations["fig1_representative_scans.csv"] = CsvExpectation(
        scan_headers,
        selected_fields(scans, scan_headers),
        ("candidate_id", "activity", "scan_step"),
        frozenset({"scan_step", "relative_energy_ev"}),
    )

    summary_headers = ("metric", "value")
    expectations["fig1_canonical_summary.csv"] = CsvExpectation(
        summary_headers,
        [
            {"metric": "retained_candidate_records", "value": len(candidates)},
            {"metric": "complete_profiles", "value": len(profiles)},
            {
                "metric": "converged_frames",
                "value": sum(
                    parse_int(row.get("converged_frames", ""), "converged_frames")
                    for row in profiles
                ),
            },
            {
                "metric": "represented_design_indices",
                "value": len({row.get("design_index", "") for row in candidates}),
            },
        ],
        ("metric",),
        frozenset({"value"}),
    )

    geometry_headers = (
        "design_index",
        "topology",
        "distance_a",
        "doping",
        "angle_deg",
        "candidate_count",
    )
    expectations["fig2_geometry.csv"] = CsvExpectation(
        geometry_headers,
        grouped_counts(
            candidates,
            ("design_index", "topology", "distance_a", "doping", "angle_deg"),
            "candidate_count",
        ),
        ("design_index", "topology", "distance_a", "doping", "angle_deg"),
        frozenset({"design_index", "distance_a", "angle_deg", "candidate_count"}),
    )

    highlighted_candidates = [
        row
        for row in candidates
        if row.get("topology", "") == "independent separated"
        and math.isclose(
            parse_float(row.get("distance_a", ""), "distance_a"),
            13.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        and row.get("doping", "") == "NS"
        and parse_int(row.get("angle_deg", ""), "angle_deg") in {75, 105}
    ]
    highlighted_headers = ("design_index", "angle_deg", "candidate_count")
    expectations["fig2_highlighted_geometry.csv"] = CsvExpectation(
        highlighted_headers,
        grouped_counts(
            highlighted_candidates,
            ("design_index", "angle_deg"),
            "candidate_count",
        ),
        ("design_index", "angle_deg"),
        frozenset({"design_index", "angle_deg", "candidate_count"}),
    )

    distance_headers = (
        "candidate_id",
        "topology",
        "distance_a",
        "doping",
        "angle_deg",
    )
    expectations["fig2_distance_points.csv"] = CsvExpectation(
        distance_headers,
        selected_fields(candidates, distance_headers),
        ("candidate_id",),
        frozenset({"distance_a", "angle_deg"}),
    )

    pair_topology_rows = grouped_counts(
        candidates, ("activity_pair", "topology"), "candidate_count"
    )
    pair_topology_headers = ("activity_pair", "topology", "candidate_count")
    expectations["fig3_activity_pair_topology.csv"] = CsvExpectation(
        pair_topology_headers,
        pair_topology_rows,
        ("activity_pair", "topology"),
        frozenset({"candidate_count"}),
    )
    pair_totals = Counter(row.get("activity_pair", "") for row in candidates)
    pair_count_headers = ("activity_pair", "candidate_count")
    pair_count_rows = [
        {"activity_pair": activity_pair, "candidate_count": count}
        for activity_pair, count in pair_totals.items()
    ]
    expectations["fig3_activity_pair_counts.csv"] = CsvExpectation(
        pair_count_headers,
        pair_count_rows,
        ("activity_pair",),
        frozenset({"candidate_count"}),
    )
    pair_percent_headers = (
        "activity_pair",
        "topology",
        "candidate_count",
        "pair_candidate_count",
        "topology_percent",
    )
    pair_percent_rows: list[dict[str, Any]] = []
    for row in pair_topology_rows:
        pair_count = pair_totals[str(row["activity_pair"])]
        candidate_count = int(row["candidate_count"])
        pair_percent_rows.append(
            {
                **row,
                "pair_candidate_count": pair_count,
                "topology_percent": 100.0 * candidate_count / pair_count,
            }
        )
    expectations["fig3_activity_pair_topology_percent.csv"] = CsvExpectation(
        pair_percent_headers,
        pair_percent_rows,
        ("activity_pair", "topology"),
        frozenset(
            {"candidate_count", "pair_candidate_count", "topology_percent"}
        ),
    )
    pair_design_headers = ("activity_pair", "design_index", "candidate_count")
    expectations["fig3_activity_pair_design.csv"] = CsvExpectation(
        pair_design_headers,
        grouped_counts(
            candidates, ("activity_pair", "design_index"), "candidate_count"
        ),
        ("activity_pair", "design_index"),
        frozenset({"design_index", "candidate_count"}),
    )

    activity_count_headers = ("activity", "profile_count")
    activity_count_rows = grouped_counts(profiles, ("activity",), "profile_count")
    expectations["fig4_activity_counts.csv"] = CsvExpectation(
        activity_count_headers,
        activity_count_rows,
        ("activity",),
        frozenset({"profile_count"}),
    )
    method_headers = ("activity", "selected_method", "profile_count")
    expectations["fig4_method_activity.csv"] = CsvExpectation(
        method_headers,
        grouped_counts(
            profiles, ("activity", "selected_method"), "profile_count"
        ),
        ("activity", "selected_method"),
        frozenset({"profile_count"}),
    )
    descriptor_headers = (
        "profile_id",
        "activity",
        "topology",
        "adsorption_energy_ev",
        "activation_metric_ev",
    )
    expectations["fig4_descriptor_points.csv"] = CsvExpectation(
        descriptor_headers,
        selected_fields(profiles, descriptor_headers),
        ("profile_id",),
        frozenset({"adsorption_energy_ev", "activation_metric_ev"}),
    )
    activity_profile_counts = Counter(row.get("activity", "") for row in profiles)
    zero_counts = Counter(
        row.get("activity", "")
        for row in profiles
        if parse_float(row.get("activation_metric_ev", ""), "activation_metric_ev")
        == 0.0
    )
    zero_headers = (
        "activity",
        "profile_count",
        "zero_activation_count",
        "zero_activation_fraction",
    )
    zero_rows = [
        {
            "activity": activity,
            "profile_count": count,
            "zero_activation_count": zero_counts[activity],
            "zero_activation_fraction": zero_counts[activity] / count,
        }
        for activity, count in activity_profile_counts.items()
    ]
    expectations["fig4_zero_activation.csv"] = CsvExpectation(
        zero_headers,
        zero_rows,
        ("activity",),
        frozenset(
            {"profile_count", "zero_activation_count", "zero_activation_fraction"}
        ),
    )

    test_headers = (
        "activity",
        "metric",
        "compared_topologies",
        "group_sizes",
        "statistic",
        "p_value",
        "epsilon_squared",
        "q_value",
        "retained_at_fdr_0_05",
    )
    expectations["fig5_topology_tests.csv"] = CsvExpectation(
        test_headers,
        selected_fields(tests, test_headers),
        ("activity", "metric"),
        frozenset({"statistic", "p_value", "epsilon_squared", "q_value"}),
        frozenset({"retained_at_fdr_0_05"}),
    )
    metric_columns = {
        "adsorption": "adsorption_energy_ev",
        "activation": "activation_metric_ev",
    }
    median_rows: list[dict[str, Any]] = []
    for test in tests:
        if not parse_bool(test.get("retained_at_fdr_0_05", "")):
            continue
        activity = test.get("activity", "")
        metric = test.get("metric", "")
        metric_column = metric_columns.get(metric)
        if metric_column is None:
            raise ValueError(f"unsupported topology-test metric: {metric!r}")
        for topology in TOPOLOGY_ORDER:
            values = [
                parse_float(row.get(metric_column, ""), metric_column)
                for row in profiles
                if row.get("activity", "") == activity
                and row.get("topology", "") == topology
            ]
            if len(values) >= 2:
                median_rows.append(
                    {
                        "activity": activity,
                        "metric": metric,
                        "topology": topology,
                        "n": len(values),
                        "median_ev": statistics.median(values),
                    }
                )
    median_headers = ("activity", "metric", "topology", "n", "median_ev")
    expectations["fig5_topology_medians.csv"] = CsvExpectation(
        median_headers,
        median_rows,
        ("activity", "metric", "topology"),
        frozenset({"n", "median_ev"}),
    )

    if set(expectations) != EXPECTED_FIGURE_SOURCE_FILES:
        missing = sorted(EXPECTED_FIGURE_SOURCE_FILES - set(expectations))
        unexpected = sorted(set(expectations) - EXPECTED_FIGURE_SOURCE_FILES)
        raise RuntimeError(
            f"figure-source expectation mismatch: missing={missing!r}; "
            f"unexpected={unexpected!r}"
        )
    return expectations


def normalized_row_key(
    row: dict[str, Any], expectation: CsvExpectation
) -> tuple[Any, ...]:
    values: list[Any] = []
    for field in expectation.keys:
        value = row.get(field, "")
        if field in expectation.numeric_fields:
            values.append(parse_float(value, field))
        elif field in expectation.bool_fields:
            values.append(parse_bool(value))
        else:
            values.append(str(value))
    return tuple(values)


def compare_csv_expectation(
    path: Path, expectation: CsvExpectation
) -> tuple[bool, str]:
    try:
        headers, actual_rows = read_csv_table(path)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        return False, f"cannot read CSV: {exc}"
    expected_headers = list(expectation.headers)
    if headers != expected_headers:
        return False, f"headers={headers!r}; expected={expected_headers!r}"
    if len(actual_rows) != len(expectation.rows):
        return False, f"rows={len(actual_rows)}; expected={len(expectation.rows)}"
    try:
        actual_keys = [normalized_row_key(row, expectation) for row in actual_rows]
        expected_keys = [
            normalized_row_key(row, expectation) for row in expectation.rows
        ]
    except ValueError as exc:
        return False, str(exc)
    duplicate_actual = [
        repr(key) for key, count in Counter(actual_keys).items() if count > 1
    ]
    duplicate_expected = [
        repr(key) for key, count in Counter(expected_keys).items() if count > 1
    ]
    if duplicate_actual or duplicate_expected:
        return (
            False,
            f"duplicate actual keys={format_sample(duplicate_actual) if duplicate_actual else 'none'}; "
            f"duplicate expected keys={format_sample(duplicate_expected) if duplicate_expected else 'none'}",
        )
    actual_by_key = dict(zip(actual_keys, actual_rows))
    expected_by_key = dict(zip(expected_keys, expectation.rows))
    if set(actual_by_key) != set(expected_by_key):
        missing = sorted(set(expected_by_key) - set(actual_by_key), key=repr)
        unexpected = sorted(set(actual_by_key) - set(expected_by_key), key=repr)
        return (
            False,
            f"missing keys={format_sample(missing) if missing else 'none'}; "
            f"unexpected keys={format_sample(unexpected) if unexpected else 'none'}",
        )
    mismatches: list[str] = []
    for key in sorted(expected_by_key, key=repr):
        actual = actual_by_key[key]
        expected = expected_by_key[key]
        for field in expectation.headers:
            try:
                if field in expectation.numeric_fields:
                    equal = math.isclose(
                        parse_float(actual.get(field, ""), field),
                        parse_float(expected.get(field, ""), field),
                        rel_tol=1e-10,
                        abs_tol=1e-12,
                    )
                elif field in expectation.bool_fields:
                    equal = parse_bool(actual.get(field, "")) == parse_bool(
                        expected.get(field, "")
                    )
                else:
                    equal = str(actual.get(field, "")) == str(
                        expected.get(field, "")
                    )
            except ValueError:
                equal = False
            if not equal:
                mismatches.append(f"{key!r}:{field}")
    if mismatches:
        return False, f"value mismatches={format_sample(mismatches)}"
    return True, f"{len(actual_rows)} rows match canonical recomputation"


def audit_figure_source_data(audit: Audit) -> None:
    observed = (
        {path.name for path in SOURCE_DATA_DIR.glob("*.csv")}
        if SOURCE_DATA_DIR.is_dir()
        else set()
    )
    missing = sorted(EXPECTED_FIGURE_SOURCE_FILES - observed)
    unexpected = sorted(observed - EXPECTED_FIGURE_SOURCE_FILES)
    audit.check(
        "figure source-data inventory",
        len(observed) == 17 and not missing and not unexpected,
        f"CSV files={len(observed)}; missing={format_sample(missing) if missing else 'none'}; "
        f"unexpected={format_sample(unexpected) if unexpected else 'none'}",
    )
    expectations = expected_figure_source_tables()
    for name in sorted(EXPECTED_FIGURE_SOURCE_FILES):
        path = SOURCE_DATA_DIR / name
        if not path.is_file():
            audit.check(f"figure source data {name}", False, "file is missing")
            continue
        passed, detail = compare_csv_expectation(path, expectations[name])
        audit.check(f"figure source data {name}", passed, detail)


def valid_file_signature(path: Path, extension: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open("rb") as handle:
        prefix = handle.read(1024)
    if extension == "png":
        return prefix.startswith(b"\x89PNG\r\n\x1a\n")
    if extension == "pdf":
        return prefix.startswith(b"%PDF-")
    if extension == "svg":
        return b"<svg" in prefix.lower()
    return False


def audit_figure_files(audit: Audit) -> None:
    invalid: list[str] = []
    for stem in EXPECTED_FIGURE_STEMS.values():
        for extension in ("png", "svg", "pdf"):
            path = FIGURE_DIR / f"{stem}.{extension}"
            if not valid_file_signature(path, extension):
                invalid.append(path.name)
    audit.check(
        "main figure files",
        not invalid,
        "five figures are present as valid non-empty PNG/SVG/PDF files"
        if not invalid
        else f"missing or invalid: {format_sample(invalid)}",
    )


def audit_figure_qa(audit: Audit) -> None:
    qa_csv_path = QA_DIR / "submission_qa.csv"
    source_audit_path = QA_DIR / "success_only_source_data_audit.json"
    missing = [
        str(path.relative_to(PUBLICATION_ROOT)).replace("\\", "/")
        for path in (qa_csv_path, source_audit_path)
        if not path.is_file()
    ]
    if missing:
        audit.check("figure QA", False, f"missing: {format_sample(missing)}")
        return
    rows = read_csv_rows(qa_csv_path)
    qa_figures = {row.get("figure", "").strip() for row in rows}
    invalid_rows: list[str] = []
    for row in rows:
        figure = row.get("figure", "?").strip()
        try:
            ready = parse_bool(row.get("submission_ready", ""))
        except ValueError:
            ready = False
        if not ready:
            invalid_rows.append(figure)
    source_audit = read_json(source_audit_path)
    source_audit_ok = (
        isinstance(source_audit, dict)
        and str(source_audit.get("status", "")).casefold() == "passed"
        and source_audit.get("retained_candidate_records") == 355
        and source_audit.get("complete_profiles") == 699
        and source_audit.get("converged_frames") == 3515
        and str(source_audit.get("retained_topology_tests")) == "8/10"
        and source_audit.get("submission_ready_figures") == 5
    )
    audit.check(
        "figure QA",
        len(rows) == 5
        and qa_figures == set(EXPECTED_FIGURE_STEMS)
        and not invalid_rows
        and source_audit_ok,
        f"submission-ready rows={len(rows) - len(invalid_rows)}/5; "
        f"source-data audit={'passed' if source_audit_ok else 'invalid'}",
    )


def local_markdown_target(markdown_path: Path, raw_target: str) -> Path | None:
    target = raw_target.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1]
    elif " " in target:
        target = target.split(" ", 1)[0]
    parsed = urlsplit(target)
    if parsed.scheme or parsed.netloc or not parsed.path:
        return None
    decoded_path = unquote(parsed.path)
    return (markdown_path.parent / decoded_path).resolve()


def audit_report_assets(audit: Audit) -> None:
    report_path = REPORT_DIR / "e2n_x1_x100_data_summary_zh.md"
    if not report_path.is_file() or report_path.stat().st_size == 0:
        audit.check(
            "Chinese data summary report",
            False,
            "missing or empty reports/e2n_x1_x100_data_summary_zh.md",
        )
        return
    markdown = report_path.read_text(encoding="utf-8-sig")
    broken_links: list[str] = []
    for raw_target in MARKDOWN_LINK_PATTERN.findall(markdown):
        target_path = local_markdown_target(report_path, raw_target)
        if target_path is not None and not target_path.exists():
            broken_links.append(raw_target)
    tables_dir = REPORT_DIR / "tables"
    missing_tables = sorted(
        name
        for name in EXPECTED_REPORT_TABLES
        if not (tables_dir / name).is_file() or (tables_dir / name).stat().st_size == 0
    )
    audit.check(
        "Chinese data summary report",
        not broken_links and not missing_tables,
        f"local Markdown links={'valid' if not broken_links else format_sample(broken_links)}; "
        f"report tables={len(EXPECTED_REPORT_TABLES) - len(missing_tables)}/{len(EXPECTED_REPORT_TABLES)}"
        + (
            f"; missing={format_sample(missing_tables)}" if missing_tables else ""
        ),
    )


def normalize_manifest_path(raw_path: Any) -> str:
    raw_value = str(raw_path)
    value = raw_value.strip()
    if raw_value != value or "\\" in value or value.startswith("./"):
        raise ValueError(f"non-canonical manifest path: {raw_path!r}")
    pure_path = PurePosixPath(value)
    if (
        not value
        or pure_path.is_absolute()
        or ".." in pure_path.parts
        or ":" in pure_path.parts[0]
    ):
        raise ValueError(f"unsafe manifest path: {raw_path!r}")
    return pure_path.as_posix()


@dataclass(frozen=True)
class ManifestEntry:
    path: str
    byte_count: int
    sha256: str


def extract_manifest_entries(manifest: Any) -> list[ManifestEntry]:
    if not isinstance(manifest, dict):
        raise ValueError("manifest root must be a JSON object")
    container = manifest.get("files")
    if not isinstance(container, list):
        raise ValueError("manifest files must be a JSON list")

    entries: list[ManifestEntry] = []
    for item in container:
        if not isinstance(item, dict):
            raise ValueError("manifest list entries must be JSON objects")
        if set(item) != {"path", "bytes", "sha256"}:
            raise ValueError(
                "manifest entries require exactly path, bytes, and sha256: "
                f"{item!r}"
            )
        raw_path = item.get("path")
        byte_count = item.get("bytes")
        digest = item.get("sha256")
        if (
            not isinstance(raw_path, str)
            or isinstance(byte_count, bool)
            or not isinstance(byte_count, int)
            or byte_count < 0
            or not isinstance(digest, str)
        ):
            raise ValueError(f"invalid manifest entry: {item!r}")
        entries.append(
            ManifestEntry(
                normalize_manifest_path(raw_path),
                byte_count,
                digest.casefold(),
            )
        )
    return entries


def streamed_file_metadata(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            byte_count += len(chunk)
            digest.update(chunk)
    return byte_count, digest.hexdigest()


def audit_release_manifest(audit: Audit) -> None:
    manifest_path = PUBLICATION_ROOT / "RELEASE_MANIFEST.json"
    if not manifest_path.is_file():
        audit.check(
            "release manifest SHA-256",
            False,
            "missing publication/RELEASE_MANIFEST.json; generate it after all release files are final",
        )
        return
    try:
        manifest = read_json(manifest_path)
        entries = extract_manifest_entries(manifest)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        audit.check("release manifest SHA-256", False, f"invalid manifest: {exc}")
        return

    metadata_errors: list[str] = []
    if manifest.get("schema_version") != 1:
        metadata_errors.append("schema_version")
    if manifest.get("hash_algorithm") != "SHA-256":
        metadata_errors.append("hash_algorithm")
    if manifest.get("path_base") != "directory containing this manifest":
        metadata_errors.append("path_base")
    if manifest.get("evidence_counts") != EXPECTED_MANIFEST_EVIDENCE_COUNTS:
        metadata_errors.append("evidence_counts")

    paths = [entry.path for entry in entries]
    duplicates = duplicate_values(paths)
    malformed_digests = sorted(
        entry.path
        for entry in entries
        if re.fullmatch(r"[0-9a-f]{64}", entry.sha256) is None
    )
    missing_files: list[str] = []
    mismatched_sizes: list[str] = []
    mismatched_hashes: list[str] = []
    for entry in entries:
        path = PUBLICATION_ROOT / PurePosixPath(entry.path)
        if not path.is_file():
            missing_files.append(entry.path)
            continue
        if re.fullmatch(r"[0-9a-f]{64}", entry.sha256) is None:
            continue
        actual_size, actual_digest = streamed_file_metadata(path)
        if actual_size != entry.byte_count:
            mismatched_sizes.append(entry.path)
        if actual_digest != entry.sha256:
            mismatched_hashes.append(entry.path)

    expected_paths = EXPECTED_RELEASE_FILES - {MANIFEST_RELATIVE_PATH}
    manifest_paths = set(paths)
    unlisted = sorted(expected_paths - manifest_paths)
    unexpected = sorted(manifest_paths - expected_paths)
    audit.check(
        "release manifest SHA-256",
        len(entries) == len(expected_paths)
        and not metadata_errors
        and not duplicates
        and not malformed_digests
        and not missing_files
        and not mismatched_sizes
        and not mismatched_hashes
        and not unlisted
        and not unexpected,
        f"entries={len(entries)}/{len(expected_paths)}; "
        f"metadata={format_sample(metadata_errors) if metadata_errors else 'valid'}; "
        f"duplicate={format_sample(duplicates) if duplicates else 'none'}; "
        f"malformed={format_sample(malformed_digests) if malformed_digests else 'none'}; "
        f"missing={format_sample(missing_files) if missing_files else 'none'}; "
        f"size mismatches={format_sample(mismatched_sizes) if mismatched_sizes else 'none'}; "
        f"hash mismatches={format_sample(mismatched_hashes) if mismatched_hashes else 'none'}; "
        f"unlisted={format_sample(unlisted) if unlisted else 'none'}; "
        f"unexpected={format_sample(unexpected) if unexpected else 'none'}",
    )


def is_link_or_junction(path: Path) -> bool:
    is_junction = getattr(path, "is_junction", None)
    return path.is_symlink() or bool(is_junction and is_junction())


def audit_release_boundary(audit: Audit) -> None:
    observed: set[str] = set()
    links: list[str] = []
    for path in PUBLICATION_ROOT.rglob("*"):
        relative = path.relative_to(PUBLICATION_ROOT).as_posix()
        if is_link_or_junction(path):
            observed.add(relative)
            links.append(relative)
        elif path.is_file():
            observed.add(relative)
    missing = sorted(EXPECTED_RELEASE_FILES - observed)
    unexpected = sorted(observed - EXPECTED_RELEASE_FILES)
    audit.check(
        "public release boundary",
        len(observed) == 79 and not missing and not unexpected and not links,
        f"files={len(observed)}/79; "
        f"missing={format_sample(missing) if missing else 'none'}; "
        f"unexpected={format_sample(unexpected) if unexpected else 'none'}; "
        f"links/junctions={format_sample(sorted(links)) if links else 'none'}",
    )


def main() -> int:
    audit = Audit()
    audit.run("canonical dataset audit", lambda: audit_canonical_dataset(audit))
    audit.run("figure source-data audit", lambda: audit_figure_source_data(audit))
    audit.run("panel source-data audit", lambda: audit_panel_index(audit))
    audit.run("main figure audit", lambda: audit_figure_files(audit))
    audit.run("figure QA audit", lambda: audit_figure_qa(audit))
    audit.run("report asset audit", lambda: audit_report_assets(audit))
    audit.run("release manifest audit", lambda: audit_release_manifest(audit))
    audit.run("release boundary audit", lambda: audit_release_boundary(audit))
    audit.print_report()
    return 1 if audit.failure_count else 0


if __name__ == "__main__":
    sys.exit(main())
