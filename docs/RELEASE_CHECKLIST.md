# Archival release checklist

## Source and repository

- [ ] The public export was built from a reviewed, clean commit.
- [ ] The repository contains no raw outputs, caches, databases, model weights,
      local paths, secrets, or unlicensed third-party assets.
- [ ] Package installation, wheel build, source build, and clean-wheel import
      pass on the supported Python versions.
- [ ] The installed wheel contains Flask templates/static assets and the app
      smoke test writes only beneath the configured runtime directory.
- [ ] Public scientific-core tests and the GitHub Actions matrix are green.
- [ ] Every destructive output operation is bounded, validated, and covered by
      a regression test.

## Data and statistics

- [ ] Canonical schema, units, enums, nullable fields, and key relationships are
      frozen and documented.
- [ ] Candidate, profile, frame, and design counts are recomputed from the exact
      release tables.
- [ ] Failure and missingness semantics are preserved; no failed value reloads
      as passed or becomes numeric zero.
- [ ] Record-level, candidate-level, deduplicated, clustered, and method
      sensitivity analyses required by the manuscript are complete.
- [ ] All statistical families, including nonsignificant results, are released.
- [ ] The data and release manifests contain final SHA-256 hashes.

## V11 manuscript assets

- [ ] Exactly five main figures and one main table are frozen.
- [ ] Every panel/table maps to claim IDs, source data, code, caption, and hash.
- [ ] Figure 5 contains no unsupported DFT/NEB barriers, rates, onset
      potentials, pH optima, or simulated experimental readouts.
- [ ] Figure 2 chemistry and natural-site labels are verified.
- [ ] Figure 3 reports method, scale, frames, state labels, and selection logic.
- [ ] Figure 4 uses the final canonical statistical analysis and shows
      uncertainty/sensitivity honestly.
- [ ] Minimum font size, line width, color accessibility, panel labels, units,
      captions, and vector/raster exports pass visual QA.
- [ ] `publication/V11_ASSET_MANIFEST.json` matches all deposited bytes.

## Metadata and rights

- [ ] Contribution-ordered authors, affiliations, corresponding author, ORCIDs,
      funding, and acknowledgments are final.
- [ ] Code, data, figure, and documentation licenses are approved by the rights
      holder and expressed with appropriate notices.
- [ ] Third-party permissions and database/model licenses are documented.
- [ ] Bundled web libraries include exact versions, checksums, copyright
      notices, and complete license texts.
- [ ] `CITATION.cff` contains release version/date, software DOI, and preferred
      article citation with no placeholders.
- [ ] GitHub topics, description, README citation, and archive metadata agree.

## Final commands

```bash
python tools/verify_public_repo.py --release-ready
python publication/scripts/verify_publication_release.py
python -m pytest -q
python -m build
```

- [ ] The commands pass from a clean clone.
- [ ] A built wheel installs and passes tests in a second clean environment.
- [ ] The release tag points to the verified commit.
- [ ] The GitHub release, archive DOI, manuscript data-availability statement,
      and repository citation all reference the same version.
