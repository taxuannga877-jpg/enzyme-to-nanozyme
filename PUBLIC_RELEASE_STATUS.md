# Public release status

**State: PRE-RELEASE — structurally organized, not approved for publication.**

This directory is a clean GitHub candidate generated from an allowlist. It is
safe to review and iterate on, but it must not receive a public archival tag,
Zenodo DOI, or “reproduces the paper” claim until every blocker below is
closed.

## Current strengths

- The public tree is separated from the large research workspace.
- Raw outputs, source structure libraries, caches, model weights, databases,
  runtime state, and machine-local paths are excluded.
- The canonical paper-facing data release has an internal checksum manifest
  and an executable verifier.
- Package metadata, CI, public-contract tests, reproducibility documentation,
  and a repository-level boundary audit are present.
- The Flask workbench source, templates, local browser assets, and package data
  are included while runtime databases and large research libraries remain excluded.
- Bundled browser assets now carry exact version notices, SHA-256 hashes, and
  the distributed license texts required by their upstream packages.
- Each generated export records source commit, dirty-state provenance, file
  sizes, and SHA-256 hashes in `BUILD_PROVENANCE.json`.

## Release blockers

| ID | Blocker | Exit condition |
| --- | --- | --- |
| R-01 | Export currently derives from a dirty working tree | Freeze a reviewed source commit and rebuild from a clean checkout |
| R-02 | Authors, affiliations, ORCIDs, release version/date, article citation, and DOI are placeholders | Validate `CITATION.cff` and GitHub/Zenodo metadata |
| R-03 | Code, data, figure, and documentation licenses are unresolved | Confirm rights holder/year and add approved license notices by artifact class |
| R-04 | Deposited figures are the verified `v6` generation, not the final V11 package | Freeze `publication/V11_ASSET_MANIFEST.json` with five figures, Table 1, captions, source data, and hashes |
| R-05 | The reviewed V11 Figure 5 contains unsupported DFT/NEB/rate/electrochemical-style claims | Replace it with a non-result design/validation roadmap or otherwise remove all unsupported outputs |
| R-06 | The broader code audit still has P1 correctness and collection blockers | Close the tracked P1 defects and run the supported Python matrix from a clean clone |
| R-07 | The public tests cover release contracts, not the full internal suite | Curate and pass a stable scientific core test suite without private fixtures |
| R-08 | Dependency ranges are declared but a final reproducible lock/environment file is not frozen | Add and verify the selected lockfile or conda environment on clean Linux/macOS installs |
| R-09 | Full upstream campaign reconstruction requires non-public raw inputs and model artifacts | Document acquisition/licensing or explicitly freeze the release as a release-level reproduction only |
| R-10 | Two rebuilt panel CSVs can differ at the final floating-point digit | Freeze V11 CSV `float_format` and require byte-stable source-data regeneration |
| R-12 | Flask workbench still has known UI/runtime audit findings | Close persistence, polling, Bootstrap-class, navigation, and claim-label issues or document them as unsupported |

## Required final gate

```bash
python tools/verify_public_repo.py --release-ready
python -m pytest -q
python -m build
```

The release is ready only when these commands pass from a clean clone, the
GitHub Actions matrix is green, and the generated V11 assets match the frozen
manifest byte-for-byte.
