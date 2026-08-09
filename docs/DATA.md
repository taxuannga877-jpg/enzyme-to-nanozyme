# Data documentation

## Canonical directory

The paper-facing snapshot is under
`publication/data/x1_x100_dataset/`. Its primary interfaces are:

| File | Unit | Key |
| --- | --- | --- |
| `candidates.csv` | Retained canonical candidate record | `candidate_id` |
| `profiles.csv` | Complete candidate-activity profile | `profile_id` |
| `designs.csv` | Requested x1-x100 design-axis position | `design_id` / `design_index` |
| `representative.csv` | Frozen representative candidate metadata | `candidate_id` |
| `representative_scans.csv` | Frozen representative scan points | candidate/activity/step fields |
| `topology_tests.csv` | Exploratory profile-level statistical test | activity/metric |

Derived tables summarize geometry occupancy, activity/method composition,
descriptor distributions, correlations, and candidate-profile linkage. The
authoritative field inventory and checksums are recorded in
`dataset_manifest.json` and `publication/RELEASE_MANIFEST.json`.

## Statistical units

- **Candidate record:** one retained `candidate_id` row in the canonical
  success-conditioned snapshot.
- **Profile:** one complete activity-specific scan summary associated with a
  candidate.
- **Frame:** one stored scan point counted within a complete profile.
- **Design:** one requested index on the x1-x100 design axis, whether or not a
  retained candidate is present.

The 346/7/2 linkage means 346 retained candidate records link to two complete
profiles, seven link to one, and two link to none. It is not an experimental
dual-activity rate or campaign-wide success fraction.

## Units and semantics

- Distances are reported in angstroms (`Å`) unless a field explicitly states
  another unit.
- Angles are reported in degrees.
- Energy-like computational descriptors are reported in electronvolts (`eV`).
- `activation_metric_ev` is a finite-coordinate screening descriptor derived
  from sampled points. It is not a transition-state free-energy barrier.
- Requested distance and angle fields are design parameters; they are not
  automatically equilibrium geometry measurements.
- A zero scan peak descriptor means no positive peak relative to the stored
  reference was present among available sampled values. It does not mean a
  barrier of zero.

## Selection and missingness

The canonical snapshot is success-conditioned. Current completeness is based
on stored profile status and converged-frame counts, not a fully public record
of every attempted profile. Missing upstream candidates and failed/intermediate
calculations are not represented as a campaign denominator.

Empty, null, and zero values must remain distinct. Never replace a failed or
missing energy with zero. Any downstream analysis should report the number of
records excluded and the exact exclusion rule.

## Public omissions

The release does not redistribute the full raw campaign, source PDB/motif
libraries, all trajectories, model weights, calculation caches, or all failed
states. Public candidate rows also omit selected metal identity, oxidation
state, source-site lineage, and design-family fields. These omissions protect
release boundaries but prevent complete upstream reconstruction.

## Analysis rules

1. Use a single declared analysis table per claim.
2. Do not treat multiple activity profiles from one candidate as independent
   materials without clustered or candidate-level sensitivity analysis.
3. Report all tested activity/metric families, including nonsignificant ones.
4. Treat topology tests as exploratory and confounded by selection and method.
5. Preserve method, backend, model/head, convergence, and failure provenance.
6. Regenerate the release manifest after any intentional data change.
