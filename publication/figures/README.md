# SPJ Main Figure Package

This directory contains the five main figures and the panel-level source data for the retained canonical x1-x100 release. Relative paths in this README and `manifest_v6.json` resolve from this directory.

## Evidence scope

- 355 retained canonical candidate records
- 699 complete activity-specific profiles
- 3515 converged scan frames within complete profiles
- 38 represented design indices, 13 declared activity pairs, 6 activities, and 3 topologies

These counts describe retained records. They do not establish 355 unique material compositions, an experimental success rate, or an unqualified set of calculable candidates.

## Figure files

| Figure | PNG | SVG | PDF |
|---|---|---|---|
| Fig. 1 | `fig1_canonical_evidence_trace.png` | `fig1_canonical_evidence_trace.svg` | `fig1_canonical_evidence_trace.pdf` |
| Fig. 2 | `fig2_geometry_landscape.png` | `fig2_geometry_landscape.svg` | `fig2_geometry_landscape.pdf` |
| Fig. 3 | `fig3_activity_pair_composition.png` | `fig3_activity_pair_composition.svg` | `fig3_activity_pair_composition.pdf` |
| Fig. 4 | `fig4_profile_descriptors.png` | `fig4_profile_descriptors.svg` | `fig4_profile_descriptors.pdf` |
| Fig. 5 | `fig5_topology_statistics.png` | `fig5_topology_statistics.svg` | `fig5_topology_statistics.pdf` |

The PNG files are 450 dpi RGB at a 7-inch target width. SVG and PDF files are
hybrid vector containers with editable text and line work; structure, heatmap,
or similar image layers may remain embedded raster content at the validated
final-size resolution.

## Source data and support

- Figure captions: `figure_captions_v6.md`
- Panel-to-source mapping: `panel_source_data_index.csv`
- Panel source data: 17 CSV files in `source_data/`, covering all 16 panels
- Representative structure used in Fig. 1A: `../data/x1_x100_dataset/representative_structure.png`
- Figure QA table: `qa/submission_qa.csv`
- Independent source-data audit: `qa/success_only_source_data_audit.json`
- Machine-readable package metadata: `manifest_v6.json`

The representative structure is a released snapshot. The representative scan trajectories are likewise released snapshots and cannot be reconstructed from `profiles.csv` alone.

## Interpretation boundary

Requested distances and angles are design parameters, not unconstrained equilibrium geometries. Computational screening descriptors are not experimental activity measurements or density-functional barriers. Empty design cells are not failure or yield denominators, and topology-labeled comparisons are observational, success-conditioned, protocol-mixed, and confounded.
