#!/usr/bin/env python3
"""Audit an E2N public-repository candidate without modifying it."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


IGNORED_LOCAL_PARTS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "htmlcov",
}
FORBIDDEN_PARTS = {
    "cache",
    "models",
    "node_modules",
    "outputs",
    "pdb_library",
    "motif_library",
}
FORBIDDEN_SUFFIXES = {
    ".db",
    ".joblib",
    ".mdb",
    ".pkl",
    ".pyc",
    ".pyo",
    ".safetensors",
    ".sqlite",
    ".sqlite3",
}
TEXT_SUFFIXES = {
    "",
    ".cff",
    ".csv",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".toml",
    ".txt",
    ".yml",
    ".yaml",
}
REQUIRED_PATHS = {
    ".e2n-public-repo",
    ".env.example",
    ".github/workflows/ci.yml",
    "BUILD_PROVENANCE.json",
    "CHANGELOG.md",
    "CITATION.cff",
    "CONTRIBUTING.md",
    "RIGHTS_AND_LICENSING.md",
    "PUBLIC_RELEASE_STATUS.md",
    "README.md",
    "SECURITY.md",
    "THIRD_PARTY_NOTICES.md",
    "docs/DATA.md",
    "docs/REPOSITORY_LAYOUT.md",
    "docs/REPRODUCIBILITY.md",
    "docs/RELEASE_CHECKLIST.md",
    "docs/SCIENTIFIC_SCOPE.md",
    "examples/quickstart.py",
    "enzyme_viewer/__init__.py",
    "enzyme_viewer/app.py",
    "enzyme_viewer/static/js/3Dmol-min.js",
    "enzyme_viewer/static/js/3Dmol-min.js.LICENSE.txt",
    "enzyme_viewer/static/vendor/licenses/3Dmol-2.5.5-BSD-3-Clause.txt",
    "enzyme_viewer/static/vendor/licenses/Bootstrap-4.3.1-MIT.txt",
    "enzyme_viewer/static/vendor/licenses/Font-Awesome-Free-6.0.0-beta3.txt",
    "enzyme_viewer/static/vendor/licenses/Popper.js-1.14.7-MIT.txt",
    "enzyme_viewer/static/vendor/licenses/jQuery-3.3.1-MIT.txt",
    "enzyme_viewer/templates/index.html",
    "nanozyme_mining/__init__.py",
    "publication/EVIDENCE_CONTRACT.md",
    "publication/RELEASE_MANIFEST.json",
    "publication/data/x1_x100_dataset/candidates.csv",
    "publication/data/x1_x100_dataset/profiles.csv",
    "publication/scripts/verify_publication_release.py",
    "pyproject.toml",
    "scripts/build_spj_data_summary_report.py",
    "scripts/build_spj_main_figures_latest.py",
    "tests/test_public_contract.py",
    "tools/verify_public_repo.py",
}
ABSOLUTE_PATH_PATTERNS = {
    "macOS volume path": re.compile(r"/Volumes/[A-Za-z0-9_. -]+/"),
    "macOS user path": re.compile(r"/Users/[A-Za-z0-9_.-]+/"),
    "Windows user path": re.compile(r"[A-Za-z]:\\\\Users\\\\[^\\\\\s]+"),
}
SECRET_PATTERNS = {
    "AWS access key": re.compile(r"AKIA[0-9A-Z]{16}"),
    "GitHub classic token": re.compile(r"ghp_[A-Za-z0-9]{30,}"),
    "GitHub fine-grained token": re.compile(r"github_pat_[A-Za-z0-9_]{30,}"),
    "Hugging Face token": re.compile(r"hf_[A-Za-z0-9]{30,}"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
}
MAX_PUBLIC_FILE_BYTES = 20 * 1024 * 1024


@dataclass
class Finding:
    level: str
    check: str
    detail: str


class Audit:
    def __init__(self) -> None:
        self.findings: list[Finding] = []

    def pass_(self, check: str, detail: str) -> None:
        self.findings.append(Finding("PASS", check, detail))

    def warn(self, check: str, detail: str) -> None:
        self.findings.append(Finding("WARN", check, detail))

    def fail(self, check: str, detail: str) -> None:
        self.findings.append(Finding("FAIL", check, detail))

    @property
    def failed(self) -> bool:
        return any(item.level == "FAIL" for item in self.findings)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_ignored_local(path: Path, root: Path) -> bool:
    relative = path.relative_to(root)
    return any(
        part in IGNORED_LOCAL_PARTS or part.endswith(".egg-info")
        for part in relative.parts
    )


def _files(root: Path) -> list[Path]:
    return [
        path
        for path in sorted(root.rglob("*"))
        if path.is_file() and not _is_ignored_local(path, root)
    ]


def _read_text(path: Path) -> str | None:
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return None
    try:
        return path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None


def _check_required(root: Path, audit: Audit) -> None:
    missing = sorted(path for path in REQUIRED_PATHS if not (root / path).is_file())
    if missing:
        audit.fail("required paths", "missing: " + ", ".join(missing))
    else:
        audit.pass_("required paths", f"all {len(REQUIRED_PATHS)} required files are present")


def _check_boundary(root: Path, files: Iterable[Path], audit: Audit) -> None:
    problems: list[str] = []
    for path in files:
        relative = path.relative_to(root)
        if any(part in FORBIDDEN_PARTS for part in relative.parts):
            problems.append(f"forbidden directory: {relative.as_posix()}")
        if relative.name.startswith("._") or relative.name in {".DS_Store", "Thumbs.db"}:
            problems.append(f"OS metadata: {relative.as_posix()}")
        if relative.suffix.lower() in FORBIDDEN_SUFFIXES:
            problems.append(f"forbidden artifact: {relative.as_posix()}")
        if path.stat().st_size > MAX_PUBLIC_FILE_BYTES:
            problems.append(
                f"file exceeds 20 MiB: {relative.as_posix()} ({path.stat().st_size} bytes)"
            )
    if problems:
        audit.fail("public boundary", "; ".join(problems[:20]))
    else:
        audit.pass_("public boundary", "no caches, raw output trees, databases, weights, or oversized files")


def _check_text(root: Path, files: Iterable[Path], audit: Audit) -> None:
    local_paths: list[str] = []
    secrets: list[str] = []
    syntax_errors: list[str] = []
    for path in files:
        relative = path.relative_to(root).as_posix()
        text = _read_text(path)
        if text is None:
            continue
        for label, pattern in ABSOLUTE_PATH_PATTERNS.items():
            if pattern.search(text):
                local_paths.append(f"{relative} ({label})")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                secrets.append(f"{relative} ({label})")
        if path.suffix == ".py":
            try:
                ast.parse(text, filename=relative)
            except SyntaxError as exc:
                syntax_errors.append(f"{relative}:{exc.lineno}: {exc.msg}")

    if local_paths:
        audit.fail("portable paths", "local absolute paths found in: " + ", ".join(local_paths))
    else:
        audit.pass_("portable paths", "no macOS or Windows user-local paths detected")
    if secrets:
        audit.fail("secret scan", "possible credentials found in: " + ", ".join(secrets))
    else:
        audit.pass_("secret scan", "no supported credential signatures detected")
    if syntax_errors:
        audit.fail("Python syntax", "; ".join(syntax_errors))
    else:
        audit.pass_("Python syntax", "all exported Python files parse successfully")


def _check_provenance(root: Path, audit: Audit) -> None:
    path = root / "BUILD_PROVENANCE.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload["files"]
    except (OSError, KeyError, json.JSONDecodeError, TypeError) as exc:
        audit.fail("build provenance", f"invalid BUILD_PROVENANCE.json: {exc}")
        return

    mismatches: list[str] = []
    declared = {record["path"]: record for record in records}
    actual = {
        path.relative_to(root).as_posix(): path
        for path in _files(root)
        if path.name != "BUILD_PROVENANCE.json"
    }
    if set(declared) != set(actual):
        missing = sorted(set(declared) - set(actual))
        extra = sorted(set(actual) - set(declared))
        if missing:
            mismatches.append("missing: " + ", ".join(missing[:10]))
        if extra:
            mismatches.append("extra: " + ", ".join(extra[:10]))
    for relative in sorted(set(declared) & set(actual)):
        record = declared[relative]
        path = actual[relative]
        if record.get("bytes") != path.stat().st_size or record.get("sha256") != _sha256(path):
            mismatches.append(f"content mismatch: {relative}")
    if mismatches:
        audit.fail("build provenance", "; ".join(mismatches[:20]))
    else:
        audit.pass_("build provenance", f"{len(actual)} exported files match the build inventory")
    if payload.get("source_worktree_dirty"):
        audit.warn("source snapshot", "export was built from a dirty working tree; freeze a clean commit for archival release")
    else:
        audit.pass_("source snapshot", "export was built from a clean source commit")


def _check_publication_release(root: Path, audit: Audit) -> None:
    verifier = root / "publication" / "scripts" / "verify_publication_release.py"
    if not verifier.is_file():
        audit.fail("publication release", "publication verifier is missing")
        return
    result = subprocess.run(
        [sys.executable, str(verifier)],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        summary = next(
            (line.strip() for line in reversed(result.stdout.splitlines()) if line.strip()),
            "publication verifier passed",
        )
        audit.pass_("publication release", summary)
    else:
        detail = (result.stdout + "\n" + result.stderr).strip().replace("\n", " | ")
        audit.fail("publication release", detail[-2000:])


def _check_release_readiness(root: Path, audit: Audit, strict: bool) -> None:
    citation = (root / "CITATION.cff").read_text(encoding="utf-8", errors="replace")
    status = (root / "PUBLIC_RELEASE_STATUS.md").read_text(encoding="utf-8", errors="replace")
    blockers: list[str] = []
    if "REQUIRED BEFORE PUBLIC RELEASE" in citation:
        blockers.append("CITATION.cff still contains required placeholders")
    if not (root / "LICENSE").is_file():
        blockers.append("no finalized root LICENSE")
    if "PRE-RELEASE" in status:
        blockers.append("PUBLIC_RELEASE_STATUS is PRE-RELEASE")
    if not (root / "publication" / "V11_ASSET_MANIFEST.json").is_file():
        blockers.append("V11 figure/table asset manifest has not been frozen")
    if blockers:
        detail = "; ".join(blockers)
        if strict:
            audit.fail("archival readiness", detail)
        else:
            audit.warn("archival readiness", detail)
    else:
        audit.pass_("archival readiness", "license, citation metadata, status, and V11 asset manifest are finalized")


def verify(root: Path, *, release_ready: bool = False) -> Audit:
    root = root.resolve()
    audit = Audit()
    if not root.is_dir():
        audit.fail("repository root", f"not a directory: {root}")
        return audit
    files = _files(root)
    _check_required(root, audit)
    _check_boundary(root, files, audit)
    _check_text(root, files, audit)
    _check_provenance(root, audit)
    _check_publication_release(root, audit)
    _check_release_readiness(root, audit, release_ready)
    return audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--release-ready",
        action="store_true",
        help="Treat unresolved license, citation, V11 asset, and status items as failures.",
    )
    args = parser.parse_args(argv)
    audit = verify(args.root, release_ready=args.release_ready)
    for finding in audit.findings:
        print(f"[{finding.level}] {finding.check}: {finding.detail}")
    counts = {
        level: sum(item.level == level for item in audit.findings)
        for level in ("PASS", "WARN", "FAIL")
    }
    print(json.dumps(counts, ensure_ascii=False))
    return 1 if audit.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
