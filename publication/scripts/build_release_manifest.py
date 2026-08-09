#!/usr/bin/env python3
"""Build the fixed-file SHA-256 inventory for the E2N publication release."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PUBLICATION_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = PUBLICATION_ROOT / "RELEASE_MANIFEST.json"
MANIFEST_RELATIVE_PATH = "RELEASE_MANIFEST.json"

DATA_FILES = {
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
FIGURE_STEMS = {
    "fig1_canonical_evidence_trace",
    "fig2_geometry_landscape",
    "fig3_activity_pair_composition",
    "fig4_profile_descriptors",
    "fig5_topology_statistics",
}
FIGURE_SOURCE_FILES = {
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
REPORT_TABLE_FILES = {
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
        *(f"data/x1_x100_dataset/{name}" for name in DATA_FILES),
        *(
            f"figures/{stem}.{extension}"
            for stem in FIGURE_STEMS
            for extension in ("pdf", "png", "svg")
        ),
        *(f"figures/source_data/{name}" for name in FIGURE_SOURCE_FILES),
        *(f"reports/tables/{name}" for name in REPORT_TABLE_FILES),
    }
)

if len(EXPECTED_RELEASE_FILES) != 79:
    raise RuntimeError(
        f"release allowlist has {len(EXPECTED_RELEASE_FILES)} paths; expected 79"
    )


def is_link_or_junction(path: Path) -> bool:
    is_junction = getattr(path, "is_junction", None)
    return path.is_symlink() or bool(is_junction and is_junction())


def observed_release_paths() -> tuple[set[str], list[str]]:
    observed: set[str] = set()
    links: list[str] = []
    for path in PUBLICATION_ROOT.rglob("*"):
        relative = path.relative_to(PUBLICATION_ROOT).as_posix()
        if is_link_or_junction(path):
            observed.add(relative)
            links.append(relative)
        elif path.is_file():
            observed.add(relative)
    return observed, sorted(links)


def release_files() -> list[Path]:
    observed, links = observed_release_paths()
    missing = sorted(EXPECTED_RELEASE_FILES - observed - {MANIFEST_RELATIVE_PATH})
    unexpected = sorted(observed - EXPECTED_RELEASE_FILES)
    if links or missing or unexpected:
        details = []
        if links:
            details.append(f"links/junctions={links!r}")
        if missing:
            details.append(f"missing={missing!r}")
        if unexpected:
            details.append(f"unexpected={unexpected!r}")
        raise RuntimeError("release boundary violation: " + "; ".join(details))
    hashed_paths = sorted(EXPECTED_RELEASE_FILES - {MANIFEST_RELATIVE_PATH})
    return [PUBLICATION_ROOT / Path(relative) for relative in hashed_paths]


def streamed_file_metadata(path: Path) -> tuple[int, str]:
    """Return byte count and SHA-256 from one file stream."""
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            byte_count += len(chunk)
            digest.update(chunk)
    return byte_count, digest.hexdigest()


def build_manifest() -> dict[str, object]:
    entries: list[dict[str, object]] = []
    for path in release_files():
        byte_count, digest = streamed_file_metadata(path)
        entries.append(
            {
                "path": path.relative_to(PUBLICATION_ROOT).as_posix(),
                "bytes": byte_count,
                "sha256": digest,
            }
        )
    return {
        "schema_version": 1,
        "release_status": "pre-release metadata blockers are documented in README.md and LICENSE_SCOPE.md",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "path_base": "directory containing this manifest",
        "hash_algorithm": "SHA-256",
        "evidence_counts": {
            "retained_canonical_candidate_records": 355,
            "complete_activity_specific_profiles": 699,
            "converged_scan_frames_within_complete_profiles": 3515,
        },
        "files": entries,
    }


def main() -> None:
    manifest = build_manifest()
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(
        json.dumps(
            {
                "manifest": str(MANIFEST_PATH),
                "files": len(manifest["files"]),
                "release_boundary_paths": len(EXPECTED_RELEASE_FILES),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
