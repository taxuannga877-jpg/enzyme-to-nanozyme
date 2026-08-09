# Repository layout

```text
enzyme-to-nanozyme/
|-- .github/workflows/ci.yml       Lightweight clean-clone checks
|-- nanozyme_mining/               Installable scientific Python package
|-- enzyme_viewer/                 Optional Flask workbench, routes, templates, assets
|   |-- assembly/                  Motif-to-structure assembly primitives
|   |-- database/                  Database and reviewed annotation adapters
|   |-- design/                    Design schemas, scoring, validation, screening
|   |-- extraction/                Catalytic motif extraction
|   |-- structure/                 PDB and local-environment parsing
|   `-- utils/                     Shared constants, configuration, HTTP helpers
|-- scripts/                       Explicit paper/release entry points
|-- publication/
|   |-- data/x1_x100_dataset/      Canonical and derived release tables
|   |-- figures/                   Deposited figures, source data, and QA
|   |-- reports/                   Recomputed audit report and tables
|   `-- scripts/                   Manifest builder and release verifier
|-- examples/quickstart.py         Minimal package inspection example
|-- tests/test_public_contract.py  Public package/data contract tests
|-- docs/                          Reproducibility and scientific documentation
|-- tools/verify_public_repo.py    Repository-level boundary audit
|-- BUILD_PROVENANCE.json          Export source state and file hashes
|-- CITATION.cff                   Draft citation metadata
|-- RIGHTS_AND_LICENSING.md        Unresolved artifact licensing decisions
`-- PUBLIC_RELEASE_STATUS.md        Explicit blockers and final gate
```

## Intentionally excluded

- `outputs/`: raw and intermediate campaigns;
- `pdb_library/` and `motif_library/`: large source collections;
- `models/`: model cards, downloads, and weights;
- caches, local virtual environments, databases, and temporary renders;
- journal author guides, templates, and third-party reference projects;
- manuscript-building workspaces and local absolute paths; and
- runtime Flask databases, logs, caches, secrets, and generated job outputs.

The excluded workspace may remain useful for ongoing research, but it is not
the source of truth for the paper-associated GitHub release. Public membership
is controlled by the exporter allowlist.

## Adding a file to the public repository

1. Decide which scientific/reproducibility purpose the file serves.
2. Confirm it is original or redistributable and assign the correct artifact
   license after the licensing model is approved.
3. Remove machine paths, credentials, private identifiers, and raw restricted
   content.
4. Add it to the exporter manifest or template, not by manual copying into the
   generated directory.
5. Rebuild the candidate, run both verifiers, and review the generated diff.
6. Update citation, data, figure, or release documentation as appropriate.
