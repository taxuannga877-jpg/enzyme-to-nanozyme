# Reproducibility guide

E2N exposes three different reproducibility levels. They must not be collapsed
into a single claim.

## Level A — release integrity and reported counts

**Supported now.** A clean clone can verify the curated paper-facing snapshot,
including canonical keys and counts, candidate-profile linkage, selected
statistics, panel source-data inventory, figure inventory, and SHA-256
checksums.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[release,test]"
python tools/verify_public_repo.py
python publication/scripts/verify_publication_release.py
python -m pytest -q
```

Expected canonical totals are 355 retained candidate records, 699 complete
activity-specific profiles, and 3515 converged frames within those profiles.
The verifier returns a nonzero status if the release inventory or checksums do
not match.

## Level B — report, statistics, and deposited figures

**Supported for the currently deposited release generation.** Rebuild the
derived audit tables and report:

```bash
python scripts/build_spj_data_summary_report.py
python publication/scripts/build_release_manifest.py
python publication/scripts/verify_publication_release.py
```

Rebuild figures into a new output directory:

```bash
python scripts/build_spj_main_figures_latest.py \
  --data-dir publication/data/x1_x100_dataset \
  --out-dir figure_rebuild
```

Compare regenerated values and panel source data with `publication/figures/`.
Rendering may vary slightly across operating systems because fonts and
graphics backends differ; numerical source tables and manifest hashes are the
primary comparison targets.

All deposited text assets use LF line endings so Git checkout does not alter
release hashes. The current figure generation can still produce last-digit
floating-point serialization differences (approximately `1e-16`) in two
panel CSVs while passing numerical recomputation. The final V11 builder must
freeze an explicit float format before claiming byte-identical regeneration.

The deposited figure files carry a `v6` generation suffix. They are not yet a
claim that the final V11 manuscript is reproducible. Final V11 release
readiness requires a new asset manifest connecting each of five figures and
Table 1 to its exact source data, caption, code entry point, and hash.

## Level C — full upstream campaign reconstruction

**Not supported by this public candidate.** The original workspace contains
large raw output trees, source PDB and motif libraries, model artifacts,
private/intermediate records, and representative raw trajectories that are
not redistributed here. Some public tables also omit metal identity,
oxidation state, source-site lineage, and design-family fields.

Consequently, this repository supports release-level verification and figure/
statistical reproduction, not exact regeneration of every upstream candidate,
geometry, atomistic relaxation, or scan frame.

## Optional atomistic backends

Install the lightweight optional group only when needed:

```bash
python -m pip install -e ".[atomistic,test]"
```

For every atomistic result, record at minimum:

- software and package versions;
- backend, model artifact, model head, and device;
- input structure identifier and SHA-256 hash;
- total charge and spin multiplicity;
- constraints/restraints and their force constants;
- optimizer, step limit, force threshold, and convergence result;
- every failed state and fallback route; and
- whether the reported geometry is requested, constrained-relaxed, or
  restraint-released.

MACE, FairChem, GPU environments, external model weights, and licensed source
databases are deliberately outside the minimal install. Add a frozen,
license-compatible environment specification before any archival claim that
depends on them.

## Clean-clone release procedure

1. Check out the exact reviewed commit in a new directory.
2. Create a new environment; do not reuse a workspace environment.
3. Install `.[release,test]`.
4. Run the default public verifier and tests.
5. Rebuild reports and figures into new directories.
6. Compare numerical source tables and V11 asset hashes.
7. Run `python tools/verify_public_repo.py --release-ready`.
8. Build both wheel and source distribution with `python -m build`.
9. Install the wheel in a second clean environment and rerun the public tests.
10. Only then create the GitHub release and archival DOI.
