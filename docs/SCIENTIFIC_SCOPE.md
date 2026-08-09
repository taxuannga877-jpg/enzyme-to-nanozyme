# Scientific scope and claim boundary

## What E2N is

E2N is a computational workflow and evidence resource for translating selected
metalloenzyme-derived constraints into explicit, topology-aware nanozyme
hypotheses. It integrates motif representation, candidate construction,
physicochemical constraints, optional atomistic evaluation, activity-specific
finite-coordinate screening, and provenance-aware release artifacts.

## Claims supported by the public release

Subject to the definitions in the evidence contract, the public snapshot can
support claims about:

- the structure and auditability of the computational workflow;
- the number and composition of retained canonical records;
- the distribution of requested design parameters in the retained snapshot;
- the availability and provenance of complete stored profiles;
- reproducible release-level summary statistics and visualizations; and
- explicit computational candidates for future validation.

## Claims not supported by the public release

The repository alone does not establish:

- experimental synthesis, catalytic activity, selectivity, stability, or
  biocompatibility;
- transition states, DFT/NEB free-energy barriers, kinetic rate constants, or
  reaction mechanisms;
- electrochemical onset potentials, pH optima, assay curves, spectroscopy, or
  microscopy observations;
- equilibrium geometries when coordinates were requested or relaxed under
  restraints;
- universal superiority of one topology, metal pair, or design family;
- validated multifunctionality for a material; or
- predictive accuracy against an independent experimental benchmark.

## Required vocabulary

Prefer:

- `enzyme-informed` or `enzyme-derived evidence` over lineage claims that have
  not been traced end-to-end;
- `candidate hypothesis` over `validated nanozyme`;
- `requested geometry` or `constrained-relaxed geometry` over `equilibrium
  structure`;
- `finite-coordinate scan descriptor` over `activation barrier`;
- `complete stored profile` over `successful catalyst`; and
- `exploratory association` over `topology effect`.

## Figure policy

Every result-bearing panel must have:

1. a unique claim identifier;
2. a canonical source-data file;
3. an executable generation entry point;
4. method and statistical unit in the caption;
5. a content hash in the V11 asset manifest; and
6. an explicit limitation when the visual could be mistaken for experiment,
   kinetics, DFT, or causal inference.

Prospective experiments may be shown only as clearly labelled design or
validation plans. They must not visually imitate collected microscopy,
spectroscopy, electrochemical, or assay results.

## Relationship to the manuscript

If manuscript wording and this document disagree, the numerical release must
be re-audited and the manuscript narrowed until each result-bearing statement
maps to a released source table. Caption changes cannot rescue an unsupported
quantitative panel.
