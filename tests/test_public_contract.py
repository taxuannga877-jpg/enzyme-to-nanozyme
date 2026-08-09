from __future__ import annotations

import csv
import importlib
import importlib.util
import json
from pathlib import Path

import nanozyme_mining
from nanozyme_mining.design.physchem_knowledge import knowledge_version
from nanozyme_mining.design.substrate_catalog import list_reaction_tasks


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "publication" / "data" / "x1_x100_dataset"


def _row_count(name: str) -> int:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def test_package_and_knowledge_schema_are_exposed() -> None:
    assert nanozyme_mining.__version__ == "0.2.0"
    assert knowledge_version()


def test_reaction_tasks_are_explicit_and_typed() -> None:
    tasks = list_reaction_tasks()
    assert tasks
    assert len({task.task_id for task in tasks}) == len(tasks)
    assert all(task.nanozyme_type and task.calculation.validation_level for task in tasks)


def test_canonical_release_counts_match_manifest() -> None:
    manifest = json.loads((DATA / "dataset_manifest.json").read_text(encoding="utf-8"))
    assert _row_count("candidates.csv") == 355
    assert _row_count("profiles.csv") == 699
    assert manifest["counts"]["candidates"] == 355
    assert manifest["counts"]["profiles"] == 699
    assert manifest["counts"]["frames"] == 3515


def test_public_candidate_has_no_final_license_claim() -> None:
    assert not (ROOT / "LICENSE").exists()
    scope = (ROOT / "RIGHTS_AND_LICENSING.md").read_text(encoding="utf-8")
    assert "Licensing is not finalized" in scope


def test_report_csv_discovery_ignores_appledouble(tmp_path: Path) -> None:
    script = ROOT / "scripts" / "build_spj_data_summary_report.py"
    spec = importlib.util.spec_from_file_location("e2n_public_report_builder", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    (tmp_path / "table.csv").write_text("key,value\na,1\n", encoding="utf-8")
    (tmp_path / "._table.csv").write_text("metadata", encoding="utf-8")
    assert [path.name for path in module._csv_files(tmp_path)] == ["table.csv"]


def test_flask_source_and_browser_assets_are_included() -> None:
    required = [
        ROOT / "enzyme_viewer" / "app.py",
        ROOT / "enzyme_viewer" / "templates" / "index.html",
        ROOT / "enzyme_viewer" / "static" / "css" / "e2n-theme.css",
        ROOT / "enzyme_viewer" / "static" / "js" / "3Dmol-min.js",
        ROOT / "enzyme_viewer" / "static" / "js" / "3Dmol-min.js.LICENSE.txt",
        ROOT / "enzyme_viewer" / "static" / "vendor" / "licenses" / "3Dmol-2.5.5-BSD-3-Clause.txt",
        ROOT / "enzyme_viewer" / "static" / "vendor" / "licenses" / "Bootstrap-4.3.1-MIT.txt",
        ROOT / "enzyme_viewer" / "static" / "vendor" / "licenses" / "Font-Awesome-Free-6.0.0-beta3.txt",
        ROOT / "enzyme_viewer" / "static" / "vendor" / "licenses" / "Popper.js-1.14.7-MIT.txt",
        ROOT / "enzyme_viewer" / "static" / "vendor" / "licenses" / "jQuery-3.3.1-MIT.txt",
        ROOT / "THIRD_PARTY_NOTICES.md",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in required)


def test_flask_app_uses_configured_runtime_root(tmp_path: Path, monkeypatch) -> None:
    data_root = tmp_path / "data"
    runtime_root = tmp_path / "runtime"
    monkeypatch.setenv("E2N_DATA_ROOT", str(data_root))
    monkeypatch.setenv("E2N_RUNTIME_DIR", str(runtime_root))
    monkeypatch.setenv("E2N_LOG_LEVEL", "WARNING")

    module = importlib.import_module("enzyme_viewer.app")
    module.app.config.update(TESTING=True)

    assert module.app.config["DATA_ROOT"] == data_root.resolve()
    assert module.app.config["RUNTIME_DIR"] == runtime_root.resolve()
    assert module.app.config["MOTIF_DB_PATH"].is_relative_to(runtime_root.resolve())
    assert module.app.config["DESIGN_OUTPUT_DIR"].is_relative_to(runtime_root.resolve())

    client = module.app.test_client()
    assert client.get("/").status_code == 200
    assert client.get("/static/css/e2n-theme.css").status_code == 200
    assert client.get("/static/js/3Dmol-min.js.LICENSE.txt").status_code == 200
