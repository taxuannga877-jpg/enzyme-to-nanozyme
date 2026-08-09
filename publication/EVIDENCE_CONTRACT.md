# E2N x1-x100 evidence contract

## Status and authority

This document defines the public evidence boundary for the E2N x1-x100 paper
release. It controls the terminology, denominators, provenance statements,
and scientific claims made from `publication/`.

The contract applies to the manuscript, figure captions, repository README,
data availability statement, supplementary material, release notes, and
archive metadata. If an older audit or manuscript version conflicts with this
document, the current canonical release tables and this contract take
precedence.

## 1. Canonical statistical units

| Unit | Contracted count | Definition | Must not be substituted with |
| --- | ---: | --- | --- |
| Retained canonical candidate record | 355 | A row with a unique `candidate_id` in `candidates.csv`, retained under the source-specific rules below | All attempted designs, campaign yield, successful experiments, or unique material compositions |
| Complete activity-specific profile | 699 | A row with a unique `profile_id`, a positive frame count, and equality of `frame_count` and `converged_frames` in `profiles.csv` | Candidate count, experimental assay, independent biological sample, or validated activity |
| Converged frame within a complete profile | 3515 | The sum of `converged_frames` across the 699 complete profiles | A table of 3515 independently keyed public frame records |
| Represented design index | 38 of 100 | An x1-x100 index occupied by at least one retained candidate record | Evidence that each empty index failed or was attempted |
| Declared activity pair | 13 | The two computational targets assigned to a candidate before profile interpretation | Experimentally demonstrated dual activity |
| Activity | 6 | An activity label represented by at least one complete profile | Six comparable assays under one uniform protocol |
| Topology | 3 | `bridged`, `independent adjacent`, or `independent separated` | A causal mechanism, stability class, or universal performance ordering |

Required headline wording:

> 355 retained canonical candidate records, 699 complete activity-specific
> profiles, and 3515 converged frames within complete profiles.

Do not shorten the first term to `355 calculable candidates`. The two source
blocks do not expose one uniform `calculable` field or one uniform selection
rule.

## 2. Candidate and profile retention

The 355 candidate records combine two source blocks with different available
metadata:

- 172 legacy-source records entered the canonical snapshot as an already
  filtered block with at least one available complete-profile endpoint. The
  common public candidate schema does not contain an explicit `calculable`
  field for these records.
- 183 newer-source records carry an explicit `calculable=True` state in their
  source table.

The two counts sum to 355, but this harmonization does not retroactively create
one campaign-wide calculability definition. The canonical release is a
retained analysis snapshot, not a denominator for every upstream attempt.

The profile snapshot contains 339 complete profiles from the legacy profile
block and 360 complete profiles from the newer profile block. Seven upstream
profile records that did not satisfy the complete-profile boundary were not
included in `profiles.csv`. Only frames belonging to retained complete
profiles contribute to the contracted total of 3515.

## 3. Candidate-profile linkage

The retained candidate records link to complete profiles as follows:

| Complete profiles linked to one candidate | Candidate records |
| ---: | ---: |
| 2 | 346 |
| 1 | 7 |
| 0 | 2 |

This 346/7/2 distribution is a foreign-key linkage audit inside the retained
snapshot. It is not:

- a campaign success or completion rate;
- a dual-functionality rate;
- an experimental activity rate;
- a profile eligibility denominator; or
- a unified calculability rate.

The two no-profile candidate IDs are `16e7f330` and `b5aad8dd`. The seven
single-profile candidate IDs are `644a5e9e`, `9034c955`, `77b47a61`,
`9b19eaa6`, `b876b6fa`, `e079d280`, and `74c8179e`. Their presence reflects
the difference between candidate retention and complete-profile retention;
it does not convert a missing profile into negative catalytic evidence.

## 4. Requested design grid

The public candidate records occupy 38 x1-x100 indices. The other 62 index
positions are blank on the common axis and have no attempted-design status in
the canonical public tables. They must not be counted as failures, rejected
candidates, or negative experiments.

Candidate rows are distributed across the three topology labels as follows:

| Topology | Retained candidate records |
| --- | ---: |
| `bridged` | 28 |
| `independent adjacent` | 101 |
| `independent separated` | 226 |

These are retained-set composition counts. They do not estimate topology
success probabilities because attempted counts, enumeration weights, source
chemistry, and selection processes are not controlled in the public snapshot.

All `distance_a` and `angle_deg` fields are requested design parameters. They
are not unconstrained equilibrium distances or angles. Dense occupancy at a
requested 13.0 A distance, NS co-doping, and requested angles of 75 or 105
degrees may be described as a represented hypothesis class or retained-grid
window. It must not be called an optimum, enrichment, yield advantage, or
equilibrium geometry.

## 5. Profiles and descriptor semantics

Each profile belongs to one candidate and one activity. A complete profile is
a computational trajectory record, not an assay result. The 699 profiles mix
five recorded computational routes:

| Computational route | Complete profiles |
| --- | ---: |
| GFN2 deep | 251 |
| First pass | 240 |
| GFN1 SCF fallback | 128 |
| GFN1 extended | 69 |
| GFN2 extended | 11 |

Method labels must remain attached to profile interpretation. Pooled energies
do not become uniform high-fidelity barriers merely because they share units.
The public descriptors are screening outputs and are not DFT-calibrated
kinetic barriers or experimental rates.

`activation_metric_ev` is a forward scan peak descriptor. A value of zero
means that no positive peak relative to the first stored point was observed
among the finite scan points. It does not establish a zero activation energy,
a barrier-free pathway, or catalytic inactivity.

## 6. Representative record and frozen snapshots

The released representative is candidate `31967c28` at x57. It has declared
targets Glucose Oxidase and Peroxidase, topology `independent separated`, a
requested distance of 13.0 A, NS co-doping, a requested angle of 105 degrees,
and exactly two complete canonical profiles.

Its selection was deterministic but not hypothesis-window filtered:

1. The dataset builder intersected 536 private structure images with canonical
   candidates linked to at least two complete profiles.
2. This produced 179 eligible records.
3. The builder selected the record closest to the eligible-set median score,
   breaking ties by design index and then candidate ID.
4. The selected record happens to be x57; the builder did not explicitly
   filter for x57 or the highlighted 13.0 A/NS/105-degree cell.

Therefore the representative illustrates traceability. It is not the best
candidate, a top-ranked catalyst, an experimentally validated material, or an
independent validation sample.

`representative_scans.csv` contains ten released rows: five scan points for
each of the two declared activities. These points were extracted from private
raw result JSON during dataset construction. `profiles.csv` contains only
profile-level descriptors, so the per-step trajectories cannot be regenerated
from `profiles.csv` alone. The scan table is a frozen, checksummed release
snapshot and must be preserved as such.

Likewise, `representative_structure.png` is a released representative asset
derived during dataset construction from a private structure-image
collection. The public package can verify its checksum and use it to rebuild
Fig. 1, but it cannot reconstruct that image from the canonical CSV tables.

## 7. Statistical tests

`topology_tests.csv` contains ten activity-by-descriptor Kruskal-Wallis omnibus
tests. P values were adjusted together by the Benjamini-Hochberg procedure;
eight tests have `q < 0.05` and two do not.

The tests are:

- profile-level rather than candidate-level;
- conditioned on the retained complete-profile snapshot;
- method-mixed;
- observational and exploratory;
- unbalanced across activity and topology; and
- potentially confounded by design family, requested geometry, support
  doping, activity, source block, and computational route.

They do not estimate a causal topology effect. A retained q value must be
reported separately from effect magnitude. No universal topology superiority,
mechanistic causality, or out-of-snapshot predictive claim is permitted.

Spearman descriptor correlations are also descriptive associations in the
same mixed profile set. They do not establish causal relationships or provide
independent-sample validation.

## 8. Public reproducibility boundary

The public release supports release-level reproduction:

- verify canonical counts, unique IDs, and foreign keys;
- regenerate deterministic summary tables;
- recalculate the released descriptive correlations and topology tests;
- rebuild the five figures from the canonical snapshot and frozen
  representative assets; and
- compare figure panels with their 17 released source-data CSV files.

It does not support a full upstream campaign rebuild. The public layer omits
raw result JSON, failed and intermediate batches, model weights, calculation
caches, the full structure-image collection, source PDB/motif databases, and
other large private output trees.

The candidate schema also omits metal A/B identities, oxidation states,
source-site lineage, and source design-family fields. Although all 355
`candidate_id` values are unique, only 240 records are unique on the visible
scientific design fields; the remaining 115 depend on opaque IDs and omitted
source-side chemistry for distinction. The release must not describe the
table as 355 unique material compositions.

## 9. Allowed and prohibited claims

Supported formulations include:

- auditable computational design and screening workbench;
- protein-to-nanomaterial translation as the paper's explicit conceptual
  framework;
- retained canonical candidate records and complete profiles;
- occupancy of the requested design grid;
- method-labeled computational descriptors;
- exploratory, success-conditioned topology associations; and
- prioritization of hypotheses for future synthesis and experimental testing.

The current release does not support claims of:

- experimentally validated catalytic activity or multifunctionality;
- direct catalytic rates or biological efficacy;
- DFT-calibrated kinetic barriers;
- barrier-free catalysis inferred from a zero forward peak descriptor;
- general material stability or unconstrained topology preservation;
- universal topology superiority or causal topology effects;
- campaign yield from the 38 occupied or 62 empty design indices;
- external benchmark outperformance;
- predictive ranking accuracy; or
- priority claims such as `first` or `unprecedented`.

Historical upstream lineage counts and partial-profile frame totals are
outside this contract. They must not be mixed with the 355/699/3515 public
snapshot in the manuscript, captions, abstract, conclusion, or repository
headline.

## 10. Change control

Any change to a canonical table, representative snapshot, figure source CSV,
figure, report, or release document requires:

1. rerunning the report generator and release verifier;
2. rebuilding affected figures into a temporary directory;
3. reviewing all changed scientific claims and denominators;
4. regenerating `RELEASE_MANIFEST.json` and every affected SHA-256 value; and
5. recording the new release version and archive DOI after deposition.

No count or definition in this contract should be changed merely to match an
older manuscript sentence. The source tables and auditable selection rules
must change first, with a documented version transition.

## 11. Citation and licensing blockers

Author names, contribution order, affiliations, corresponding author details,
funding, non-author acknowledgments, archival DOI, release version,
copyright holder/year, and final code/data/figure license choices remain to be
confirmed. Do not invent these fields. The current placeholders in
`CITATION.cff` and `LICENSE_SCOPE.md` are release blockers, not publication
metadata.
