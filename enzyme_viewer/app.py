import os
import secrets
import threading
from pathlib import Path
from typing import Optional
from flask import Flask, g, jsonify, request, send_file
from flask_cors import CORS
import py3Dmol
from enzyme_viewer.images import get_structure_html_and_active_data
from enzyme_viewer.http_headers import apply_standard_response_headers
from enzyme_viewer.routes.activity_validation import (
    ActivityValidationRouteServices,
    register_activity_validation_routes,
)
from enzyme_viewer.routes.catalog import CatalogRouteServices, register_catalog_routes
from enzyme_viewer.routes.design_metadata import register_design_metadata_routes
from enzyme_viewer.routes.design_jobs import (
    DesignJobRouteServices,
    register_design_job_routes,
)
from enzyme_viewer.routes.motif_basic import (
    MotifBasicRouteServices,
    register_motif_basic_routes,
)
from enzyme_viewer.routes.motif_listing import (
    MotifListingRouteServices,
    register_motif_listing_routes,
)
from enzyme_viewer.routes.motif_structure import (
    MotifStructureRouteServices,
    generate_pdb_from_residues,
    register_motif_structure_routes,
)
from enzyme_viewer.routes.ligand import LigandRouteServices, register_ligand_routes
from enzyme_viewer.routes.pages import register_page_routes
from enzyme_viewer.routes.structure import (
    StructureRouteServices,
    register_structure_routes,
)
from enzyme_viewer.structure_info import (
    _build_pdb_info_response,
    _extract_active_sites,
    _extract_metal_sites,
    _extract_structure_info,
    _format_metal_sites_from_extractor,
    _infer_geometry,
)

from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher
from nanozyme_mining.extraction.extractor import MotifExtractor
from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser
from enzyme_viewer.motif_db import MotifDatabase, classify_motif
from enzyme_viewer.catalytic_metal_db import CatalyticMetalDatabase
from enzyme_viewer.ligand_db import LigandDatabase

app = Flask(__name__)

# PR0-2 (v4 audit): security configuration migrated to enzyme_viewer/security.py.
# - N-H4: CORS limited to local origins by default; override via FLASK_CORS_ORIGINS=*
# - C4: install unified error handler (no traceback to client; rid for log correlation)
from enzyme_viewer.security import (
    CORS_ORIGINS,
    DEBUG as _SEC_DEBUG,
    HOST  as _SEC_HOST,
    PORT  as _SEC_PORT,
    install_global_error_handler,
    error_response,
    env_int,                       # PR4-1 (M28): safe int env-var parser
    is_valid_pdb_id,
    is_valid_ec_number,
    is_valid_motif_id,
    install_basic_auth,
    load_secret_key,
    require_json_csrf,
    safe_join,
)

_GZIP_MIN_BYTES = env_int(
    "E2N_GZIP_MIN_BYTES",
    1024,
    min_value=0,
    max_value=16 * 1024 * 1024,
)


if CORS_ORIGINS == "*":
    CORS(app)
else:
    CORS(app, origins=CORS_ORIGINS)

install_global_error_handler(app)
install_basic_auth(app, host=_SEC_HOST)

@app.before_request
def set_csp_nonce():
    g.csp_nonce = secrets.token_urlsafe(16)


@app.context_processor
def inject_csp_nonce():
    return {"csp_nonce": getattr(g, "csp_nonce", "")}


@app.after_request
def add_security_headers(response):
    return apply_standard_response_headers(
        response,
        request=request,
        gzip_min_bytes=_GZIP_MIN_BYTES,
        script_nonce=getattr(g, "csp_nonce", None),
    )

# Configurable data/runtime roots keep an installed wheel from writing into
# ``site-packages``. Defaults preserve the source-checkout behavior when the
# app is launched from the repository root.
DATA_ROOT = Path(os.environ.get("E2N_DATA_ROOT", Path.cwd())).expanduser().resolve()
RUNTIME_DIR = Path(
    os.environ.get("E2N_RUNTIME_DIR", DATA_ROOT / ".runtime")
).expanduser().resolve()
DB_DIR = RUNTIME_DIR / "db"

app.config['DATA_ROOT'] = DATA_ROOT
app.config['RUNTIME_DIR'] = RUNTIME_DIR
load_secret_key(app, app.config['RUNTIME_DIR'] / 'secret_key')
app.config['PDB_LIBRARY_DIR'] = Path(
    os.environ.get("E2N_PDB_LIBRARY_DIR", DATA_ROOT / "pdb_library")
).expanduser().resolve()
app.config['MOTIF_LIBRARY_DIR'] = Path(
    os.environ.get("E2N_MOTIF_LIBRARY_DIR", DATA_ROOT / "motif_library")
).expanduser().resolve()
app.config['MOTIF_OUTPUT_DIR'] = Path(
    os.environ.get("E2N_MOTIF_OUTPUT_DIR", RUNTIME_DIR / "motifs")
).expanduser().resolve()
app.config['MOTIF_DB_PATH'] = Path(
    os.environ.get("E2N_MOTIF_DB_PATH", DB_DIR / "motif_index.db")
).expanduser().resolve()
app.config['CATALYTIC_METAL_DB_PATH'] = Path(
    os.environ.get("E2N_CATALYTIC_METAL_DB_PATH", DB_DIR / "catalytic_metal_index.db")
).expanduser().resolve()
app.config['LIGAND_DB_PATH'] = Path(
    os.environ.get("E2N_LIGAND_DB_PATH", DB_DIR / "ligand_index.db")
).expanduser().resolve()
app.config['DESIGN_OUTPUT_DIR'] = Path(
    os.environ.get("E2N_DESIGN_OUTPUT_DIR", RUNTIME_DIR / "outputs" / "design")
).expanduser().resolve()
app.config['ACTIVITY_VALIDATION_OUTPUT_DIR'] = Path(
    os.environ.get(
        "E2N_ACTIVITY_VALIDATION_OUTPUT_DIR",
        RUNTIME_DIR / "outputs" / "activity_validation",
    )
).expanduser().resolve()
app.config['ACTIVITY_VALIDATION_REFERENCE_DIR'] = Path(
    os.environ.get("E2N_ACTIVITY_VALIDATION_REFERENCE_DIR", DATA_ROOT / "参考图示")
).expanduser().resolve()
app.config['ACTIVITY_VALIDATION_REFERENCE_FIGURES'] = {
    "structure_comparison": {
        "label": "Strict structure validation reference",
        "path": Path(
            os.environ.get(
                "E2N_STRUCTURE_COMPARISON_REFERENCE",
                DATA_ROOT
                / "outputs"
                / "physchem_comparison_reports"
                / "report_20260616_135828"
                / "fig_structure_comparison.png",
            )
        ),
    },
    "adsorption_volcano_reference": {
        "label": "Adsorption volcano reference",
        "path": Path(
            os.environ.get(
                "E2N_ADSORPTION_VOLCANO_REFERENCE",
                app.config['ACTIVITY_VALIDATION_REFERENCE_DIR'] / "吸附火山图.jpg",
            )
        ),
    },
    "barrier_profile_reference": {
        "label": "Barrier profile reference",
        "path": Path(
            os.environ.get(
                "E2N_BARRIER_PROFILE_REFERENCE",
                app.config['ACTIVITY_VALIDATION_REFERENCE_DIR'] / "过渡态能垒.jpg",
            )
        ),
    },
}

# EC号到酶活性标签的映射 — 已迁移到 nanozyme_mining.utils.ec_mappings 作为单一真源
# (PR0-1: NEW-2 fix)
# 此处 re-export 保持原有 get_ec_activity_label / EC_ACTIVITY_LABELS 调用者向后兼容。
# 新增 1.11.1.21 KatG=CAT/POD 双功能精确条目，避免 prefix "1.11.1" 误判为 POD only。
from nanozyme_mining.utils.ec_mappings import (
    EC_ACTIVITY_LABELS,
    get_ec_activity_label,
)

# 向后兼容：保留旧路径（但不再使用）
app.config['CACHE_DIR'] = Path(
    os.environ.get("E2N_CACHE_DIR", RUNTIME_DIR / "cache")
).expanduser().resolve()
app.config['JSON_CACHE_DIR'] = app.config['CACHE_DIR'] / 'json'  # 仅用于向后兼容
app.config['PDB_CACHE_DIR'] = app.config['CACHE_DIR'] / 'pdb'  # 仅用于向后兼容

# 确保文件夹存在
DB_DIR.mkdir(parents=True, exist_ok=True)
app.config['PDB_LIBRARY_DIR'].mkdir(parents=True, exist_ok=True)
app.config['MOTIF_OUTPUT_DIR'].mkdir(parents=True, exist_ok=True)
app.config['DESIGN_OUTPUT_DIR'].mkdir(parents=True, exist_ok=True)
app.config['ACTIVITY_VALIDATION_OUTPUT_DIR'].mkdir(parents=True, exist_ok=True)

# PR4-1 (M27 fix): replace startup print() with structured logging. The init
# block runs once on import — keeping it as print made it impossible to
# silence under pytest / suppress in production. logger respects
# E2N_LOG_LEVEL env var (default INFO).
import logging as _app_logging
_app_logging.basicConfig(
    level=os.environ.get("E2N_LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
)
_log = _app_logging.getLogger("e2n.app")

_log.info("Using local data: PDB library = %s", app.config['PDB_LIBRARY_DIR'])

# 初始化功能模块
_log.info("Initializing modules ...")
uniprot_fetcher = UniProtFetcher(
    cache_dir=str(app.config['CACHE_DIR']),
    pdb_library_dir=str(app.config['PDB_LIBRARY_DIR'])
)
_log.info("  - UniProt Fetcher ready")

# Motif提取器
motif_extractor = MotifExtractor(output_dir=str(app.config['MOTIF_OUTPUT_DIR']))
_log.info("  - Motif Extractor ready")

# Motif数据库（延迟初始化）
_db_init_lock = threading.RLock()
motif_db = None


def _get_lazy_db(store_name, config_key, factory, loaded_message, missing_logger):
    """Initialize one module-level database store under the shared lock."""

    store = globals()[store_name]
    if store is not None:
        return store

    with _db_init_lock:
        store = globals()[store_name]
        if store is not None:
            return store

        db_path = app.config[config_key]
        if not db_path.exists():
            missing_logger(db_path)
            return None

        store = factory(db_path)
        globals()[store_name] = store
        _log.info(loaded_message, db_path)
        return store


def get_motif_db():
    """获取Motif数据库实例（延迟初始化）"""
    return _get_lazy_db(
        "motif_db",
        "MOTIF_DB_PATH",
        MotifDatabase,
        "Motif database loaded: %s",
        lambda db_path: _log.warning(
            "Motif database not found: %s; build it with "
            "python enzyme_viewer/motif_db.py",
            db_path,
        ),
    )

# 催化金属位点库（延迟初始化）
catalytic_metal_db = None
def get_catalytic_metal_db():
    """获取高精度催化金属位点库实例（延迟初始化）"""
    return _get_lazy_db(
        "catalytic_metal_db",
        "CATALYTIC_METAL_DB_PATH",
        CatalyticMetalDatabase,
        "Catalytic metal-site database loaded: %s",
        lambda db_path: _log.warning(
            "Catalytic metal-site database not found: %s; build it with "
            "python -m enzyme_viewer.catalytic_metal_db "
            "--pdb-library %s --motif-library %s --out %s",
            db_path,
            app.config['PDB_LIBRARY_DIR'],
            app.config['MOTIF_LIBRARY_DIR'],
            db_path,
        ),
    )

# 配体/辅因子库（延迟初始化）
ligand_db = None
def get_ligand_db():
    """获取配体/辅因子数据库实例（延迟初始化）"""
    return _get_lazy_db(
        "ligand_db",
        "LIGAND_DB_PATH",
        LigandDatabase,
        "Ligand/cofactor database loaded: %s",
        lambda db_path: _log.warning(
            "Ligand/cofactor database not found: %s; build it with "
            "python -m enzyme_viewer.ligand_db --pdb-library %s --out %s",
            db_path,
            app.config['PDB_LIBRARY_DIR'],
            db_path,
        ),
    )

@app.teardown_appcontext
def close_db_connections(_exc=None):
    """Close per-thread SQLite handles held by lazy global stores."""
    for store in (motif_db, catalytic_metal_db, ligand_db):
        if store is not None and hasattr(store, "close"):
            store.close()

register_page_routes(app)
register_design_metadata_routes(app)

def get_json_file_path(ec_number: str) -> Path:
    """获取JSON文件路径（优先从pdb_library，向后兼容旧路径）。

    PR0-3 (H1 fix): the legacy fallback used to call `get_json_file_path(ec_number)`
    recursively with the same argument, which is an infinite loop when both
    PDB_LIBRARY_DIR/<ec>/...sites.json AND JSON_CACHE_DIR exist but the file is
    not under PDB_LIBRARY_DIR. Now the fallback resolves directly to
    JSON_CACHE_DIR without recursion.
    """
    ec_dir_name = ec_number.replace(".", "_")
    json_file = app.config['PDB_LIBRARY_DIR'] / ec_dir_name / f"{ec_number}_sites.json"
    if json_file.exists():
        return json_file
    # Legacy fallback: try old cache location, but DO NOT recurse.
    legacy = app.config['JSON_CACHE_DIR'] / f"{ec_number}_sites.json"
    return legacy


def _resolve_pdb_library_file(pdb_id: str, ec_number: str) -> Path:
    """Resolve a PDB ID to a server-owned file under PDB_LIBRARY_DIR."""
    normalized_pdb_id = str(pdb_id or "").upper().strip()
    normalized_ec = str(ec_number or "").strip()
    if not is_valid_pdb_id(normalized_pdb_id):
        raise ValueError("invalid pdb_id (expected 4-char alphanumeric)")
    if not is_valid_ec_number(normalized_ec):
        raise ValueError("invalid ec_number")

    pdb_dir = safe_join(app.config['PDB_LIBRARY_DIR'], normalized_ec.replace(".", "_"))
    import glob as _glob
    for candidate in pdb_dir.glob(f"*{_glob.escape(normalized_pdb_id)}*.pdb"):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(normalized_pdb_id)

register_catalog_routes(
    app,
    CatalogRouteServices(
        pdb_library_dir=lambda: app.config['PDB_LIBRARY_DIR'],
        json_cache_dir=lambda: app.config['JSON_CACHE_DIR'],
        motif_library_dir=lambda: app.config['MOTIF_LIBRARY_DIR'],
        get_json_file_path=get_json_file_path,
        get_motif_db=get_motif_db,
        get_ec_activity_label=get_ec_activity_label,
    ),
)



def _path_under_any_root(path: Path, roots) -> bool:
    resolved = path.resolve()
    for root in roots:
        try:
            resolved.relative_to(Path(root).resolve())
            return True
        except ValueError:
            continue
    return False


def _resolve_motif_json_file(motif_id: str) -> Optional[Path]:
    roots = [app.config['MOTIF_OUTPUT_DIR'], app.config['MOTIF_LIBRARY_DIR']]

    output_file = safe_join(app.config['MOTIF_OUTPUT_DIR'], f"{motif_id}.json")
    if output_file.exists():
        return output_file

    db = get_motif_db()
    if db:
        db_motif = db.get_by_id(motif_id)
        if db_motif:
            db_path = Path(db_motif.get("file_path") or "")
            if db_path.exists() and _path_under_any_root(db_path, roots):
                return db_path

    library_root = app.config['MOTIF_LIBRARY_DIR']
    for candidate in library_root.rglob(f"{motif_id}.json"):
        if candidate.is_file() and _path_under_any_root(candidate, roots):
            return candidate
    return None

register_motif_basic_routes(
    app,
    MotifBasicRouteServices(
        resolve_motif_json_file=_resolve_motif_json_file,
        resolve_pdb_library_file=_resolve_pdb_library_file,
        get_json_file_path=get_json_file_path,
        motif_extractor=motif_extractor,
        motif_output_dir=lambda: app.config['MOTIF_OUTPUT_DIR'],
    ),
)

register_ligand_routes(
    app,
    LigandRouteServices(
        resolve_pdb_library_file=_resolve_pdb_library_file,
    ),
)

register_motif_listing_routes(
    app,
    MotifListingRouteServices(
        motif_output_dir=lambda: app.config['MOTIF_OUTPUT_DIR'],
        motif_library_dir=lambda: app.config['MOTIF_LIBRARY_DIR'],
        get_motif_db=lambda: get_motif_db(),
        get_catalytic_metal_db=lambda: get_catalytic_metal_db(),
        get_ligand_db=lambda: get_ligand_db(),
        classify_motif=lambda motif_data: classify_motif(motif_data),
        path_under_any_root=_path_under_any_root,
    ),
)

register_motif_structure_routes(
    app,
    MotifStructureRouteServices(
        motif_output_dir=lambda: app.config['MOTIF_OUTPUT_DIR'],
        motif_library_dir=lambda: app.config['MOTIF_LIBRARY_DIR'],
        get_motif_db=lambda: get_motif_db(),
        get_catalytic_metal_db=lambda: get_catalytic_metal_db(),
        resolve_motif_json_file=_resolve_motif_json_file,
        path_under_any_root=_path_under_any_root,
    ),
)



from enzyme_viewer.design_serialization import (
    _assembly_from_dict,
    _json_safe,
    _loads_subprocess_json,
    _score_from_dict,
)
from enzyme_viewer.design_io import (
    _atoms_to_xyz,
    _coord_distance,
    _parse_pdb_atoms,
    _parse_xyz_atoms,
    _reconstruct_cores_for_loaded_design,
)
from enzyme_viewer.design_store import (
    _activity_validation_reference_figures,
    _activity_validation_runtime_context,
    _activity_validation_structure_diagnostics,
    _activity_validation_task_dir,
    _assembly_result_payload,
    _design_result_dir,
    _get_design_result,
    _load_design_result_from_disk,
    _persist_design_result,
    _persisted_validation_artifacts,
    _score_payload,
    _validation_snapshot_from_disk,
    configure_design_store,
)
from enzyme_viewer.design_subprocess import (
    _run_assemble_subprocess,
    _run_catalysis_screen_subprocess,
)
from enzyme_viewer.activity_validation_worker import (
    _resolve_validation_activities,
    _run_activity_validation_worker,
    _validation_progress_callback,
    _validation_snapshot_for_client,
    configure_activity_validation_worker,
)
from enzyme_viewer.runtime_cache import _ActivityValidationCache, _DesignResultCache
from nanozyme_mining.design.substrate_catalog import get_reaction_task

from concurrent.futures import ThreadPoolExecutor

_DESIGN_CACHE_SIZE = env_int(
    "E2N_DESIGN_CACHE_SIZE",
    128,
    min_value=1,
    max_value=4096,
)
_DESIGN_CACHE_TTL = env_int(
    "E2N_DESIGN_CACHE_TTL",
    3600,
    min_value=1,
    max_value=30 * 24 * 3600,
)

_design_results = _DesignResultCache(_DESIGN_CACHE_SIZE, _DESIGN_CACHE_TTL)
configure_design_store(config=app.config, design_results=_design_results)

_STRUCTURE_RENDER_CACHE_SIZE = env_int(
    "E2N_STRUCTURE_RENDER_CACHE_SIZE",
    64,
    min_value=1,
    max_value=4096,
)
_STRUCTURE_RENDER_CACHE_TTL = env_int(
    "E2N_STRUCTURE_RENDER_CACHE_TTL",
    1800,
    min_value=1,
    max_value=30 * 24 * 3600,
)
_STRUCTURE_RENDER_CACHE_MAX_BYTES = env_int(
    "E2N_STRUCTURE_RENDER_CACHE_MAX_BYTES",
    50 * 1024 * 1024,
    min_value=1024,
    max_value=1024 * 1024 * 1024,
)
_structure_render_cache = _DesignResultCache(
    _STRUCTURE_RENDER_CACHE_SIZE,
    _STRUCTURE_RENDER_CACHE_TTL,
)


def _site_labels_fingerprint(site_labels):
    if not site_labels:
        return ()
    normalized = []
    for residue_index, site_type in site_labels.items():
        try:
            normalized.append((int(residue_index), int(site_type)))
        except (TypeError, ValueError):
            normalized.append((str(residue_index), str(site_type)))
    return tuple(sorted(normalized))


def _structure_render_cache_key(enzyme_structure_path, site_labels, view_size, show_active):
    path = Path(enzyme_structure_path).resolve()
    stat = path.stat()
    return (
        str(path),
        stat.st_mtime_ns,
        stat.st_size,
        tuple(view_size),
        bool(show_active),
        _site_labels_fingerprint(site_labels),
    )


def _get_structure_render_cached(*, enzyme_structure_path, site_labels, view_size, show_active):
    key = None
    try:
        key = _structure_render_cache_key(
            enzyme_structure_path,
            site_labels,
            view_size,
            show_active,
        )
    except OSError:
        key = None

    if key is not None and key[2] <= _STRUCTURE_RENDER_CACHE_MAX_BYTES:
        cached = _structure_render_cache.get(key)
        if cached is not None:
            return cached

    result = get_structure_html_and_active_data(
        enzyme_structure_path=enzyme_structure_path,
        site_labels=site_labels,
        view_size=view_size,
        show_active=show_active,
    )
    if key is not None and key[2] <= _STRUCTURE_RENDER_CACHE_MAX_BYTES:
        _structure_render_cache[key] = result
    return result


register_structure_routes(
    app,
    StructureRouteServices(
        pdb_library_dir=lambda: app.config['PDB_LIBRARY_DIR'],
        get_json_file_path=get_json_file_path,
        render_structure_cached=_get_structure_render_cached,
    ),
)

_VALIDATION_CACHE_SIZE = env_int(
    "E2N_VALIDATION_CACHE_SIZE",
    64,
    min_value=1,
    max_value=4096,
)
_VALIDATION_CACHE_TTL = env_int(
    "E2N_VALIDATION_CACHE_TTL",
    7200,
    min_value=1,
    max_value=30 * 24 * 3600,
)
_VALIDATION_MAX_WORKERS = env_int(
    "E2N_VALIDATION_MAX_WORKERS",
    2,
    min_value=1,
    max_value=32,
)


_activity_validation_jobs = _ActivityValidationCache(
    _VALIDATION_CACHE_SIZE,
    _VALIDATION_CACHE_TTL,
)
_activity_validation_executor = ThreadPoolExecutor(
    max_workers=max(1, _VALIDATION_MAX_WORKERS),
    thread_name_prefix="e2n-validation",
)
configure_activity_validation_worker(
    config=app.config,
    jobs=_activity_validation_jobs,
    run_catalysis_screen_subprocess=_run_catalysis_screen_subprocess,
)


register_activity_validation_routes(
    app,
    ActivityValidationRouteServices(
        reference_figure_specs=lambda: app.config.get(
            'ACTIVITY_VALIDATION_REFERENCE_FIGURES'
        ) or {},
        get_design_result=_get_design_result,
        get_reaction_task=get_reaction_task,
        score_payload=_score_payload,
        runtime_context=_activity_validation_runtime_context,
        structure_diagnostics=_activity_validation_structure_diagnostics,
        reference_figures=_activity_validation_reference_figures,
        resolve_activities=_resolve_validation_activities,
        jobs=_activity_validation_jobs,
        executor=_activity_validation_executor,
        run_worker=_run_activity_validation_worker,
        snapshot_for_client=_validation_snapshot_for_client,
        snapshot_from_disk=_validation_snapshot_from_disk,
    ),
)


register_design_job_routes(
    app,
    DesignJobRouteServices(
        design_results=_design_results,
        persist_design_result=_persist_design_result,
        get_design_result=_get_design_result,
        design_result_dir=_design_result_dir,
        score_payload=_score_payload,
        get_reaction_task=get_reaction_task,
        run_assemble_subprocess=_run_assemble_subprocess,
        run_catalysis_screen_subprocess=_run_catalysis_screen_subprocess,
        assembly_from_dict=_assembly_from_dict,
    ),
)


if __name__ == '__main__':
    # PR0-2 (C3 fix): debug + host driven by FLASK_DEBUG / FLASK_HOST env vars
    # Default is host=127.0.0.1 + debug=off (production-safe).
    # For development with live reload, set FLASK_DEBUG=1 FLASK_HOST=0.0.0.0
    app.run(debug=_SEC_DEBUG, host=_SEC_HOST, port=_SEC_PORT)
