# x1-x100 Public Analysis Dataset

This directory is the canonical public data interface for the main figures, report, and retained-set statistical summaries.

## Evidence scope

- 355 retained canonical candidate records
- 699 complete activity-specific profiles
- 3515 converged scan frames within complete profiles
- 38 of 100 design indices represented by at least one retained candidate record
- 13 declared activity pairs, 6 activities, and 3 topologies

The 62 unrepresented design indices are empty requested-grid positions, not failed candidates or a yield denominator.

## Canonical tables

- `candidates.csv`: one row per retained canonical candidate record, uniquely keyed by `candidate_id`
- `profiles.csv`: one row per complete activity-specific profile, linked to `candidates.csv` by `candidate_id`
- `designs.csv`: the complete requested x1-x100 design index

The 355 candidate IDs are unique, but the public candidate columns do not identify 355 unique material compositions. In particular, the release omits metal identities, oxidation states, source-site lineage, and source design family; only 240 distinct combinations remain after grouping the 355 rows by all visible scientific design fields. The opaque candidate IDs preserve record identity but do not replace the omitted chemistry.

## Derived tables and released snapshots

The remaining non-representative CSV files are deterministic summaries of the canonical tables. `representative.csv`, `representative_scans.csv`, and `representative_structure.png` are released snapshots associated with candidate `31967c28`. The builder intersected 536 private structure images with retained candidates having at least two complete profiles, yielding 179 eligible records, and then selected the record closest to the eligible-set median score with deterministic `design_index` and `candidate_id` tie-breaking. It did not explicitly filter for x57 or the highlighted 13.0 Å/NS/105° window. The scan trajectory was extracted from private raw result records and cannot be regenerated from `profiles.csv` alone, so the representative bundle should be treated as a fixed release artifact rather than a publicly rebuildable selection.

The five main figures use 17 panel-level source-data CSV files in `../../figures/source_data/`; their panel mapping is recorded in `../../figures/panel_source_data_index.csv`. The representative structure used in Fig. 1A is `representative_structure.png` in this directory.

## Reproducibility boundary

This package supports verification of the released counts, figure panels, descriptor summaries, and topology tests. It does not contain the full upstream candidate chemistry, raw scan results, or private structure-image pool and therefore does not support reconstruction of the entire upstream campaign from scratch. Requested distances and angles are design parameters rather than unconstrained equilibrium geometries, and computational descriptors are not experimental activity measurements or density-functional barriers.
