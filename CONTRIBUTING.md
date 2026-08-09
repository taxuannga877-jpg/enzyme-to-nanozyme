# Contributing to E2N

E2N combines scientific software with a paper-facing evidence release.
Contributions must therefore preserve both software behavior and the stated
evidence boundary.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[app,release,test]"
python tools/verify_public_repo.py
python -m pytest -q
```

Install `.[atomistic]` only when a change actually exercises ASE, RDKit, or
tblite. MACE, FairChem, GPU, and externally licensed model artifacts require a
separate documented environment and are not part of the lightweight CI path.

## Change requirements

- Keep each pull request focused and use a regression test for behavior fixes.
- Use explicit input and output paths in analysis scripts.
- Record backend, model/head, charge, spin, constraints, convergence criteria,
  and failure states for calculation changes.
- Distinguish requested geometry, relaxed geometry, screening descriptors,
  predicted properties, and experimental measurements.
- Update source data and hashes whenever a figure value changes.
- Do not commit raw campaigns, PDB libraries, model weights, caches, databases,
  credentials, journal templates, or private machine paths.

## Commit style

Use Conventional Commit prefixes such as `fix:`, `feat:`, `docs:`, `test:`,
`refactor:`, and `chore:`. Protect `main`; merge reviewed feature branches only
after CI passes.

## Pull request checklist

- [ ] The change has a narrow scientific/software purpose.
- [ ] New behavior has a regression test or a written reason why it cannot.
- [ ] Public verification and relevant tests pass.
- [ ] Claims remain within `docs/SCIENTIFIC_SCOPE.md`.
- [ ] Data and figure changes include provenance and updated checksums.
- [ ] No secret, private path, large artifact, or third-party restricted file
      was added.
- [ ] Any remaining limitation is stated in the pull request.
