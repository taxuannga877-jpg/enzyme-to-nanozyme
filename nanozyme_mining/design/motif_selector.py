"""
从 catalytic_metal_index.db / motif_index.db 查询辅助数据，
支持前端三步向导的三个 GET API。

ligand_index.db 保留为Web只读参考，不作为设计推荐依据。

PR1-5 (v4 audit):
- H13 (HIGH): get_second_shell 用 motif_index DB 查询替代 glob 全树扫描。
  motif_library 当前 313 MB / 数千 JSON，glob+json.load 是 O(N) IO；
  改为 SQL 查询是 O(log N) 索引命中。如果 motif_index.db 缺失则回退到 glob。
- H14 (HIGH): 模块级单例 + threading.local 连接复用，替代每次调用都
  sqlite3.connect + close 的模式。
- N-M2 (MEDIUM): _db_path 用 pathlib + parents 替代字符串拼接 '..'，
  更健壮（不依赖 cwd，跨平台正确）。
"""
import sqlite3
import threading
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

from .constraint_scorer import VALID_COORDINATION_RANGES, coordination_cutoff
from .physchem_knowledge import donor_distance_range, get_activity_prototype


# PR1-5 (N-M2 fix): pathlib + parents instead of os.path.join('..', '..').
# Path(__file__).parents[2] is the project root regardless of cwd.
_ROOT = Path(__file__).resolve().parents[2]
_ENZYME_VIEWER = _ROOT / "enzyme_viewer"
_MOTIF_LIBRARY = _ROOT / "motif_library"


def _db_path(name: str) -> str:
    return str(_ENZYME_VIEWER / name)


# PR1-5 (H14 fix): per-thread connection cache, opened once per (thread, db_name).
# Flask debug server is threaded=True; without per-thread caching, two simultaneous
# requests would race on a single shared sqlite3 connection.
_LOCAL = threading.local()


def _get_conn(db_name: str) -> sqlite3.Connection:
    cache = getattr(_LOCAL, "conns", None)
    if cache is None:
        cache = _LOCAL.conns = {}
    conn = cache.get(db_name)
    if conn is None:
        conn = sqlite3.connect(
            _db_path(db_name),
            check_same_thread=False,
        )
        # WAL + busy_timeout (mirror motif_db.py settings).
        try:
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA busy_timeout=5000")
        except sqlite3.OperationalError:
            pass
        cache[db_name] = conn
    return conn


def get_activity_metals(nanozyme_type: str) -> List[Dict]:
    """Step1 辅助：返回该活性下金属类型及出现频率"""
    conn = _get_conn("catalytic_metal_index.db")
    rows = conn.execute(
        "SELECT metal_type, COUNT(*) as cnt FROM catalytic_metal_site WHERE nanozyme_type=? GROUP BY metal_type ORDER BY cnt DESC",
        (nanozyme_type,)
    ).fetchall()
    total = sum(r[1] for r in rows)
    return [{"metal_type": r[0], "count": r[1], "percentage": round(r[1] / total * 100, 1) if total else 0} for r in rows]


def get_coord_templates(
    metal_type: str,
    nanozyme_type: str,
    oxidation_state: Optional[int] = None,
) -> List[Dict]:
    """Step2 辅助：返回该金属在该活性下的典型配位层模板（按配位几何聚类）"""
    import json
    metal_key = (metal_type or "").upper()
    ox = _normalize_oxidation_state(oxidation_state)
    conn = _get_conn("catalytic_metal_index.db")
    rows = conn.execute(
        """SELECT coordination_geometry, coordination_number, coordinating_residues_json,
                  coordination_distances_json, COUNT(*) as cnt
           FROM catalytic_metal_site
           WHERE metal_type=? AND nanozyme_type=?
           GROUP BY coordination_geometry, coordination_number
           ORDER BY cnt DESC""",
        (metal_type, nanozyme_type)
    ).fetchall()

    templates = []
    for row in rows:
        geom, cn, res_json, dist_json, cnt = row
        try:
            residues = json.loads(res_json) if res_json else []
            distances = json.loads(dist_json) if dist_json else []
        except Exception:
            continue
        # 聚合配位原子规格
        coord_atoms = []
        for i, res in enumerate(residues):
            coord_atoms.append({
                "donor_element": _donor_element(res.get("atom_name", ""), res.get("residue_name", "")),
                "residue_name": res.get("residue_name", ""),
                "atom_name": res.get("atom_name", ""),
                "bond_length": round(distances[i], 3) if i < len(distances) else 2.0,
            })
        label = _make_template_label(metal_type, geom, coord_atoms)
        templates.append({
            "name": label,
            "geometry": geom,
            "coordination_number": cn,
            "coord_atoms": coord_atoms,
            "source_pdb_count": cnt,
            "source": "catalytic_metal_index",
        })

    curated = _curated_coord_template(metal_key, nanozyme_type, ox)
    if curated and (
        not templates
        or not _is_reasonable_template(templates[0], metal_key, ox)
        or (ox is not None and not _matches_template_profile(templates[0], curated))
    ):
        templates = [curated] + templates
    return templates


# PR1-5 (H13 fix): cache the result of the (expensive) second-shell query.
# Inputs are bounded (~11 nanozyme types × ~10 metals = ~110 combos), so an
# LRU of 128 covers every realistic working set without unbounded growth.
@lru_cache(maxsize=128)
def _get_second_shell_cached(nanozyme_type: str, metal_type: str) -> tuple:
    """Internal: returns tuple of dict-like tuples so lru_cache stays hashable."""
    rows = _second_shell_from_db(nanozyme_type, metal_type)
    if rows is None:
        rows = _second_shell_from_glob(nanozyme_type)
    return tuple(
        (r["residue_name"], r["atom_name"], r["role"], r["frequency"])
        for r in rows
    )


def get_second_shell(nanozyme_type: str, metal_type: str) -> List[Dict]:
    """Step3 辅助：返回该活性+金属组合下的典型第二配位层残基"""
    return [
        {"residue_name": r, "atom_name": a, "role": role, "frequency": cnt}
        for r, a, role, cnt in _get_second_shell_cached(nanozyme_type, metal_type)
    ]


def _second_shell_from_db(nanozyme_type: str, metal_type: str) -> Optional[List[Dict]]:
    """
    PR1-5 (H13 fix): query motif_index.db instead of glob-scanning motif_library/.

    motif_index.db has (nanozyme_type, file_path) indexed; we read only the
    matching JSON files instead of opening every JSON in the library. Returns
    None if the DB is missing or has no matching rows, so the caller can fall
    back to the (slow) glob path.
    """
    import json
    db_path = _db_path("motif_index.db")
    if not Path(db_path).exists():
        return None
    try:
        conn = _get_conn("motif_index.db")
        rows = conn.execute(
            "SELECT file_path FROM motif_index WHERE UPPER(nanozyme_type) = UPPER(?)",
            (nanozyme_type,),
        ).fetchall()
    except sqlite3.Error:
        return None
    if not rows:
        return None
    counts: Dict[tuple, int] = defaultdict(int)
    for (file_path,) in rows:
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            continue
        for atom in data.get("anchor_atoms", []):
            atom_name = str(atom.get("atom_name", "") or "").upper()
            if atom_name == "CA":
                continue
            key = (atom.get("residue_name", ""),
                    atom.get("atom_name", ""),
                    atom.get("role", "unknown"))
            counts[key] += 1
    if not counts:
        return None
    result = [
        {"residue_name": r, "atom_name": a, "role": role, "frequency": cnt}
        for (r, a, role), cnt in sorted(counts.items(), key=lambda x: -x[1])
    ]
    return result[:20]


def _second_shell_from_glob(nanozyme_type: str) -> List[Dict]:
    """Legacy fallback path: glob-scan motif_library. Used only when
    motif_index.db is unavailable."""
    import glob
    import json
    import os
    counts: Dict[tuple, int] = defaultdict(int)
    for path in glob.glob(str(_MOTIF_LIBRARY / "**" / "*.json"), recursive=True):
        try:
            with open(path) as f:
                data = json.load(f)
            if data.get("nanozyme_type") != nanozyme_type:
                continue
            for atom in data.get("anchor_atoms", []):
                atom_name = str(atom.get("atom_name", "") or "").upper()
                if atom_name == "CA":
                    continue
                key = (atom.get("residue_name", ""),
                        atom.get("atom_name", ""),
                        atom.get("role", "unknown"))
                counts[key] += 1
        except Exception:
            continue
    result = [
        {"residue_name": r, "atom_name": a, "role": role, "frequency": cnt}
        for (r, a, role), cnt in sorted(counts.items(), key=lambda x: -x[1])
    ]
    return result[:20]


def _donor_element(atom_name: str, residue_name: str) -> str:
    if atom_name.startswith("N"):
        return "N"
    if atom_name.startswith("O"):
        return "O"
    if atom_name.startswith("S"):
        return "S"
    return "N"


def _make_template_label(metal: str, geom: str, coord_atoms: List[Dict]) -> str:
    from collections import Counter
    res_counts = Counter(a["residue_name"] for a in coord_atoms)
    res_str = "".join(f"{r}{n}" for r, n in sorted(res_counts.items()))
    geom_short = {"octahedral": "Oct", "tetrahedral": "Tet", "square_planar": "SqPl",
                  "square_pyramidal": "SqPy", "trigonal_bipyramidal": "TBP"}.get(geom, geom)
    return f"{metal}({geom_short}) {res_str}"


_COMMON_OXIDATION_STATE = {
    "FE": 3,
    "CU": 2,
    "ZN": 2,
    "MN": 2,
    "CO": 2,
    "NI": 2,
    "MO": 4,
    "W": 4,
    "V": 4,
    "CR": 3,
    "RU": 3,
    "PD": 2,
    "PT": 2,
    "AU": 1,
    "AG": 1,
}


def _normalize_oxidation_state(value) -> Optional[int]:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _is_reasonable_template(
    template: Dict,
    metal_type: str = "",
    oxidation_state: Optional[int] = None,
) -> bool:
    geom = str(template.get("geometry") or "").strip().lower()
    cn = int(template.get("coordination_number") or 0)
    atoms = template.get("coord_atoms") or []
    metal = str(metal_type or "").upper()
    if geom in {"", "unknown", "none"} or cn <= 0 or len(atoms) != cn:
        return False

    ox = oxidation_state if oxidation_state is not None else _COMMON_OXIDATION_STATE.get(metal)
    valid_range = VALID_COORDINATION_RANGES.get((metal, ox)) if ox is not None else None
    if valid_range:
        lo, hi = valid_range
        if not (lo <= cn <= hi):
            return False

    for atom in atoms:
        donor = _donor_element(atom.get("atom_name", ""), atom.get("residue_name", ""))
        try:
            bond_length = float(atom.get("bond_length", 0.0))
        except (TypeError, ValueError):
            return False
        if bond_length <= 0:
            return False
        if metal:
            supported_range = donor_distance_range(metal, ox or 0, donor)
            if supported_range and not (supported_range[0] <= bond_length <= supported_range[1]):
                return False
            if not supported_range and bond_length >= coordination_cutoff(metal, donor):
                return False
    return True


def _matches_template_profile(template: Dict, curated: Dict) -> bool:
    return (
        str(template.get("geometry") or "").lower()
        == str(curated.get("geometry") or "").lower()
        and int(template.get("coordination_number") or 0)
        == int(curated.get("coordination_number") or 0)
    )


def _curated_coord_template(
    metal_type: str,
    nanozyme_type: str,
    oxidation_state: Optional[int] = None,
) -> Dict:
    activity = str(nanozyme_type or "").lower()
    metal = str(metal_type or "").upper()
    ox = oxidation_state if oxidation_state is not None else _COMMON_OXIDATION_STATE.get(metal)

    if "superoxide" in activity and metal == "ZN":
        atoms = [
            _coord_atom("N", "HIS", "NE2", 2.05),
            _coord_atom("N", "HIS", "NE2", 2.05),
            _coord_atom("N", "HIS", "NE2", 2.05),
            _coord_atom("O", "ASP", "OD1", 2.00),
        ]
        geom = "tetrahedral"
    elif activity == "peroxidase" and metal == "FE":
        # The mined Fe-peroxidase records are often under-specified heme contacts.
        # A compact His4 square-planar proxy is the most stable first-shell default
        # for the current adsorption/reaction-coordinate screening path; users can
        # still add axial O/N ligands manually when they need a heme-like model.
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)]
        geom = "square_planar"
    elif metal == "CU" and ox == 1:
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(2)]
        geom = "linear"
    elif metal == "AU" and ox == 1:
        atoms = [
            _coord_atom("S", "CYS", "SG", _bond_length(metal, "S")),
            _coord_atom("S", "CYS", "SG", _bond_length(metal, "S")),
        ]
        geom = "linear"
    elif metal == "AG" and ox == 1:
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(2)]
        geom = "linear"
    elif "catalase" in activity and metal in {"MN", "FE"}:
        # PR1-2 (M1 fix): Mn(III/IV) catalase active sites are octahedral with
        # bridging oxos / waters, not tetrahedral. 4 oxygens + 2 water/hydroxo
        # is a reasonable octahedral approximation when no curated template exists.
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "ASP", "OD2", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE2", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O",   _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O",   _bond_length(metal, "O")),
        ]
        geom = "octahedral"
    elif metal == "CU":
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)]
        geom = "square_planar"
    elif metal == "ZN":
        atoms = [
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
        ]
        geom = "tetrahedral"
    elif metal == "FE":
        # PR1-2 (M2 fix): FE specifically — heme-like octahedral, not square_planar.
        atoms = [
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
        ]
        geom = "octahedral"
    elif metal == "MN":
        # PR1-2 (M2 fix): Mn — octahedral default.
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O",   _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O",   _bond_length(metal, "O")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
        ]
        geom = "octahedral"
    elif metal == "CO":
        # PR1-2 (M2): Co(II/III) octahedral default.
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)] + [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O",   _bond_length(metal, "O")),
        ]
        geom = "octahedral"
    elif metal == "NI":
        # PR1-2 (M2): Ni(II) d8 square-planar default.
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)]
        geom = "square_planar"
    elif metal == "MO":
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "ASP", "OD2", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
        ]
        geom = "octahedral"
    elif metal == "V" and ox in {3, 4}:
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "ASP", "OD2", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
        ]
        geom = "octahedral"
    elif metal == "V" and ox == 5:
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
        ]
        geom = "tetrahedral"
    elif metal == "CR" and ox in {2, 3}:
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "ASP", "OD2", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
            _coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")),
        ]
        geom = "octahedral"
    elif metal == "CR" and ox == 6:
        atoms = [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "GLU", "OE1", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
        ]
        geom = "tetrahedral"
    elif metal in {"PD", "PT"} and ox == 2:
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)]
        geom = "square_planar"
    elif metal == "PT" and ox == 4:
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)] + [
            _coord_atom("O", "ASP", "OD1", _bond_length(metal, "O")),
            _coord_atom("O", "HOH", "O", _bond_length(metal, "O")),
        ]
        geom = "octahedral"
    else:
        atoms = [_coord_atom("N", "HIS", "NE2", _bond_length(metal, "N")) for _ in range(4)]
        geom = "tetrahedral"

    equatorial_count = {
        "octahedral": 4,
        "square_pyramidal": 4,
        "trigonal_bipyramidal": 3,
        "tetrahedral": 3,
    }.get(geom, len(atoms))
    for index, atom in enumerate(atoms):
        supported_range = donor_distance_range(metal, ox or 0, atom["donor_element"])
        atom["role"] = "equatorial_network" if index < equatorial_count else "axial_labile"
        atom["bond_length_range"] = list(supported_range) if supported_range else None
        atom["labile"] = atom["role"] == "axial_labile"
        atom["protonation_state"] = (
            "hydroxo" if atom["residue_name"] == "OH" else "aquo"
            if atom["role"] == "axial_labile" else None
        )
        atom["source_id"] = "curated_coordination_ranges"
    prototype = get_activity_prototype(nanozyme_type) or {}
    return {
        "name": f"{metal}{f'(+{ox})' if ox is not None else ''} {geom} curated {''.join(a['donor_element'] for a in atoms)}",
        "geometry": geom,
        "coordination_number": len(atoms),
        "coord_atoms": atoms,
        "source_pdb_count": 0,
        "source": "curated_fallback",
        "reason": "database template was missing, unknown, or under-coordinated",
        "prototype_id": prototype.get("prototype_id"),
        "allowed_modes": list(prototype.get("allowed_modes", ())),
        "condition_id": prototype.get("condition_id"),
        "microstates": list(prototype.get("microstates", ())),
        "proxy_label": prototype.get("proxy_label"),
    }


def _coord_atom(element: str, residue: str, atom: str, bond_length: float) -> Dict:
    return {
        "donor_element": element,
        "residue_name": residue,
        "atom_name": atom,
        "bond_length": bond_length,
    }


def _bond_length(metal: str, donor: str) -> float:
    ranges = {
        ("FE", "N"): 2.05, ("FE", "O"): 2.02,
        ("MN", "N"): 2.15, ("MN", "O"): 2.10,
        ("CU", "N"): 2.00, ("CU", "O"): 1.98,
        ("ZN", "N"): 2.05, ("ZN", "O"): 2.00,
        ("CO", "N"): 1.98, ("CO", "O"): 1.95,
        ("NI", "N"): 1.95, ("NI", "O"): 1.95,
    }
    return ranges.get((metal.upper(), donor.upper()), 2.05)
