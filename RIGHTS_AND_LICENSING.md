# License scope and release blockers

> [!CAUTION]
> Licensing is not finalized. This document records unresolved scope; it does
> not grant a new license. Confirm the real rights holder, copyright year, and
> licenses before publishing or archiving this release.

## Current repository state

The repository root contains an MIT license notice naming `Nanozyme Design
Team` and the year 2024. The project owner has not yet confirmed that name,
year, or the intended artifact scope for this paper release. They must not be
treated as verified author or rights-holder metadata.

The Research journal's published content is expected to be Open Access under
CC BY, but that journal policy does not by itself settle the pre-publication
license for repository code, canonical CSV data, figure source data, rendered
figures, reports, manuscript drafts, or third-party inputs.

## Artifact status

| Artifact class | Examples | Current status |
| --- | --- | --- |
| Source code | `scripts/*.py`, `publication/scripts/*.py` | License choice pending confirmation of rights holder/year and intended MIT scope |
| Canonical and derived data | `publication/data/**/*.csv`, `publication/figures/source_data/*.csv`, `publication/reports/tables/*.csv` | Data license not yet selected; no license is granted by this scope note |
| Original figures and captions | `publication/figures/*.png`, `*.svg`, `*.pdf`, caption files | Figure/content license not yet selected; no license is granted by this scope note |
| Reports and release documentation | `publication/reports/*.md`, this directory's Markdown files | Documentation/content license not yet selected |
| Manuscript and supplementary files | DOCX, PDF, Markdown manuscript drafts | Governed by author and journal decisions; not automatically covered by a software license |
| Third-party or source assets | Journal guides, externally sourced databases, PDB inputs, model weights, external images | Excluded unless a separate file explicitly documents permission and redistribution terms |

Absence of a selected license is not permission to redistribute or reuse an
artifact. Repository visibility and copyright permission are separate issues.

## Decisions required before public release

> [!IMPORTANT]
> Complete every item below before creating the archival tag or DOI.

1. Confirm the legal copyright holder name and correct copyright year or year
   range.
2. Decide whether the root MIT license applies only to source code or to a
   broader set of software/documentation files.
3. Select a data license for canonical CSV files and derived source data.
4. Select a content license for original figures, captions, and reports.
5. Verify that every released asset is original, licensed, public-domain, or
   covered by documented permission.
6. Add explicit per-directory license notices if code, data, and content use
   different licenses.
7. Update `CITATION.cff` with the confirmed SPDX identifier or identifiers.
8. Regenerate `RELEASE_MANIFEST.json` after all notices are final.

One common model, if the rights holder approves it, is MIT for original source
code and CC BY 4.0 for original data, figures, and documentation. This sentence
records an option only; it is not the current license selection.

## Exclusions

No license chosen for the public release should be interpreted as granting
rights to:

- journal templates, author guides, or publisher branding;
- third-party database content beyond its source license;
- model weights or software governed by separate upstream licenses;
- private raw calculations, caches, or local `outputs/` trees not included in
  the release; or
- trademarks, personal data, or author identity rights.

When an external asset is necessary for scientific provenance but cannot be
redistributed, cite its source and document how an authorized user can obtain
it instead of copying it into `publication/`.
