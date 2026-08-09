# E2N publication release

This directory is the public, paper-facing release layer for E2N. It packages
the canonical x1-x100 snapshot, panel-level figure source data, five main
figures, release checks, and a detailed Chinese data report.

The associated manuscript is being prepared as a Research Article for
`Research` (Science Partner Journal), with topical alignment to the special
issue `Nanozymes: Rational and Intelligent Design, Mechanistic Elucidation,
and Emerging Biological Applications`. The repository frames E2N as a
computational methods and resource contribution, not as an experimentally
validated catalyst report.

> [!IMPORTANT]
> This release is not ready for archival publication until the author list,
> affiliations, corresponding author, funding, repository DOI, copyright
> holder and year, artifact licenses, and non-author acknowledgments are
> completed. See `CITATION.cff` and `LICENSE_SCOPE.md`.

Files carrying a `v6` suffix identify the current figure-package generation;
they do not set the archival release version. The final tag, root
`CITATION.cff`, archive DOI, and release manifest metadata must use one version
after the author and license blockers are resolved.

## Canonical release scope

The controlling public evidence set contains:

- 355 retained canonical candidate records;
- 699 complete activity-specific profiles;
- 3515 converged frames within those complete profiles;
- 38 represented design indices from the x1-x100 axis;
- 13 declared activity pairs, 6 activities, and 3 topology modes;
- 10 profile-level topology tests, 8 retained at `q < 0.05` after
  Benjamini-Hochberg correction;
- 17 panel-level source-data CSV files; and
- 5 main figures, each supplied as PNG, SVG, and PDF.

These nouns are separate statistical units. In particular, the 355 records
must not be described without qualification as all attempted candidates,
355 unique material compositions, or a campaign success denominator. The
346/7/2 distribution describes only how retained candidate records link to
two, one, or zero complete profiles. It is not a success rate, dual-activity
rate, experimental activity rate, or unified completion rate.

The full interpretation rules are in [`EVIDENCE_CONTRACT.md`](EVIDENCE_CONTRACT.md).
If another release document uses less precise wording, the evidence contract
controls.

## Directory map

```text
publication/
|-- data/x1_x100_dataset/       Canonical and derived release tables
|-- figures/                    Five figures in PNG, SVG, and PDF
|   |-- source_data/            Seventeen panel-level CSV files
|   `-- qa/                     Figure and source-data QA records
|-- reports/                    Chinese technical report and audit tables
|-- scripts/                    Release verifier and manifest builder
|-- CITATION.cff                Citation metadata, currently blocked on authors/DOI
|-- EVIDENCE_CONTRACT.md        Controlling definitions and claim boundaries
|-- LICENSE_SCOPE.md            Artifact-by-artifact licensing status
|-- RELEASE_MANIFEST.json       Release inventory and SHA-256 checksums
`-- requirements-release.txt   Minimal Python dependencies for this layer
```

The canonical tables are under `data/x1_x100_dataset/`. The three primary
interfaces are:

- `candidates.csv`: one row per retained canonical candidate record;
- `profiles.csv`: one row per complete activity-specific profile; and
- `designs.csv`: the common x1-x100 design axis.

Other CSV files in that directory are deterministic summaries or released
representative assets. `representative_scans.csv` is a special provenance
case: it is a frozen snapshot extracted from private raw result JSON. The
per-step trajectory cannot be reconstructed from `profiles.csv`, which stores
only profile-level descriptors. The representative structure image is also a
released snapshot rather than a reconstruction target.

## Verify the release

Use Python 3.10 or newer from the repository root:

```bash
python -m pip install -r publication/requirements-release.txt
python publication/scripts/verify_publication_release.py
```

The verifier checks the canonical row counts and keys, candidate-profile
linkage, represented design dimensions, representative records, topology
tests, figure/source-data inventory, figure QA, manifest checksums, and the
absence of prohibited raw campaign artifacts in `publication/`.

After any intentional release-file change, rebuild the checksum inventory and
rerun the verifier:

```bash
python publication/scripts/build_release_manifest.py
python publication/scripts/verify_publication_release.py
```

Regenerate the Chinese technical report and its derived tables with:

```bash
python scripts/build_spj_data_summary_report.py
```

Rebuild the figure package into a disposable comparison directory with:

```bash
python scripts/build_spj_main_figures_latest.py \
  --data-dir publication/data/x1_x100_dataset \
  --out-dir .runtime/publication_figure_rebuild
```

Do not point `--out-dir` at `publication/figures`; the builder clears its
output directory before writing. Compare a verified temporary rebuild with
the released files before replacing any release artifact.

PDF readback QA uses `pypdfium2` from the release requirements, with the
external `pdftoppm` executable as a fallback when available.

## What this release reproduces

From a clean clone, the release is intended to support:

- validation of the 355/699/3515 canonical evidence counts;
- regeneration of the report-level derived tables;
- independent recalculation of the reported Spearman correlations and ten
  Kruskal-Wallis tests with global Benjamini-Hochberg correction;
- rebuilding of the five figures from the public canonical snapshot; and
- panel-level comparison against the 17 released source-data CSV files.

It does not support a full upstream campaign rebuild. The release omits large
raw output trees, private structure-image collections, raw result JSON,
source PDB/motif databases, model weights, calculation caches, and failed or
intermediate batches. The public candidate table also omits complete metal
identity, oxidation-state, source-site lineage, and source design-family
fields. Consequently, this repository supports release-level figure,
statistical, and report reproduction, not reconstruction of every upstream
enumeration and calculation step.

## Scientific interpretation

The release supports an auditable computational design and screening
workbench. It does not establish experimental catalytic activity, validated
multifunctionality, kinetic barriers, universal topology superiority,
equilibrium geometries, general material stability, clinical utility, or
predictive ranking accuracy.

Distances and angles in the public tables are requested design parameters.
Zero values of the forward scan peak descriptor mean that no positive peak
was observed among the stored scan points relative to the first point; they
do not mean a zero barrier or barrier-free catalysis. Topology statistics are
success-conditioned, method-mixed, profile-level, exploratory associations
and must not be interpreted causally.

## Reports and citation

The detailed Chinese audit is available at
[`reports/e2n_x1_x100_data_summary_zh.md`](reports/e2n_x1_x100_data_summary_zh.md).
It records data-quality checks, profile provenance, topology tests, limitations,
and manuscript revision recommendations.

Do not cite a placeholder DOI or placeholder author. Complete and validate
`CITATION.cff`, place the finalized copy at the repository root for GitHub
discovery, create the archival release, insert its DOI, and regenerate
`RELEASE_MANIFEST.json` before public deposition.
