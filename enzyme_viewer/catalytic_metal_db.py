"""
Catalytic metal site library (SQLite)
====================================

This is the high-precision metal index used by the PDB-only mining pipeline.
It intentionally does not store every metal from a PDB file. A site is included
only when the metal has plausible catalytic chemistry plus direct evidence from
UniProt active-site annotations, PDB records, or extracted motif anchors.
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from nanozyme_mining.utils.constants import (
    CATALYTIC_METAL_RESIDUE_NAMES,
    CONDITIONAL_CATALYTIC_METAL_ELEMENTS,
    NON_CATALYTIC_METAL_ELEMENTS,
)

CORE_CATALYTIC_METALS = CATALYTIC_METAL_RESIDUE_NAMES
CONDITIONAL_METALS = CONDITIONAL_CATALYTIC_METAL_ELEMENTS
NON_CATALYTIC_METALS = NON_CATALYTIC_METAL_ELEMENTS
DONOR_RESIDUES = {"HIS", "ASP", "GLU", "CYS", "TYR", "SER", "THR", "MET", "ASN", "GLN"}


def _json_dumps(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def _json_loads(s: str) -> Any:
    return json.loads(s) if s else None


def _normalize_metal_type(metal_type: str) -> str:
    return (metal_type or "").upper().strip()


def _coords(value: Any) -> Optional[List[float]]:
    if not isinstance(value, (list, tuple)) or len(value) < 3:
        return None
    try:
        return [float(value[0]), float(value[1]), float(value[2])]
    except (TypeError, ValueError):
        return None


def _distance(a: List[float], b: List[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _residue_number(item: Dict[str, Any]) -> Optional[int]:
    raw = item.get("residue_number", item.get("residue_id", item.get("res_id", item.get("start"))))
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _chain_id(item: Dict[str, Any]) -> str:
    return str(item.get("chain_id", item.get("chain", "")) or "")


def _residue_keys(item: Dict[str, Any]) -> set:
    number = _residue_number(item)
    if number is None:
        return set()
    chain = _chain_id(item)
    return {(chain, number), ("", number)}


def _expanded_annotation_keys(active_sites: Optional[List[Dict[str, Any]]]) -> set:
    keys = set()
    for site in active_sites or []:
        start = site.get("start", site.get("residue_number", site.get("position")))
        end = site.get("end", start)
        chain = _chain_id(site)
        try:
            start_i = int(start)
            end_i = int(end)
        except (TypeError, ValueError):
            continue
        for number in range(start_i, end_i + 1):
            keys.add((chain, number))
            keys.add(("", number))
    return keys


def _annotation_is_strong(site: Dict[str, Any]) -> bool:
    text = " ".join(
        str(site.get(k, ""))
        for k in ("type", "description", "feature_type", "metal_role", "evidence")
    ).lower()
    return any(token in text for token in ("active", "catalytic", "metal", "binding", "cofactor"))


def _has_active_site_evidence(site: Dict[str, Any], active_sites: Optional[List[Dict[str, Any]]]) -> bool:
    if not active_sites:
        return False
    coord_keys = set()
    for residue in site.get("coordinating_residues") or []:
        coord_keys.update(_residue_keys(residue))
    metal_keys = _residue_keys(site)
    annotation_keys = _expanded_annotation_keys(active_sites)
    strong_annotation = any(_annotation_is_strong(annotation) for annotation in active_sites)
    return bool(strong_annotation and (coord_keys | metal_keys) & annotation_keys)


def _has_anchor_evidence(site: Dict[str, Any], motif_anchor_atoms: Optional[List[Dict[str, Any]]]) -> bool:
    anchors = motif_anchor_atoms or []
    if not anchors:
        return False

    anchor_keys = set()
    for atom in anchors:
        anchor_keys.update(_residue_keys(atom))
    for residue in site.get("coordinating_residues") or []:
        if anchor_keys & _residue_keys(residue):
            return True

    metal_coords = _coords(site.get("metal_coords") or site.get("coordinates"))
    if metal_coords:
        for atom in anchors:
            atom_coords = _coords(atom.get("coordinates"))
            if atom_coords and _distance(metal_coords, atom_coords) <= 6.0:
                return True
    return False


def _has_pdb_record_evidence(site: Dict[str, Any], pdb_record_residue_keys: Optional[set]) -> bool:
    if not pdb_record_residue_keys:
        return False
    site_keys = set()
    site_keys.update(_residue_keys(site))
    for residue in site.get("coordinating_residues") or []:
        site_keys.update(_residue_keys(residue))
    return bool(site_keys & pdb_record_residue_keys)


def _donor_profile(site: Dict[str, Any]) -> Tuple[bool, List[str]]:
    donors = []
    for residue in site.get("coordinating_residues") or []:
        residue_name = str(residue.get("residue_name", "")).upper()
        atom_name = str(residue.get("atom_name", "")).upper()
        if residue_name in DONOR_RESIDUES or atom_name[:1] in {"N", "O", "S"}:
            donors.append(residue_name or atom_name[:1])
    return len(donors) >= 1, donors


def evaluate_catalytic_metal_site(
    site: Dict[str, Any],
    *,
    ec_number: str = "",
    nanozyme_type: str = "",
    active_sites: Optional[List[Dict[str, Any]]] = None,
    pdb_record_residue_keys: Optional[set] = None,
    motif_anchor_atoms: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """
    Score one metal site using a high-precision catalytic policy.

    Returns a dict containing is_catalytic, evidence_score, evidence_level,
    included_reason, excluded_reason, and evidence_tags.
    """
    metal_type = _normalize_metal_type(site.get("metal_type"))
    coord_residues = site.get("coordinating_residues") or []
    res_names = [str(r.get("residue_name", "")).upper() for r in coord_residues]
    res_set = set(res_names)

    evidence_tags: List[str] = []
    excluded_reason = ""
    score = 0

    if metal_type in NON_CATALYTIC_METALS:
        excluded_reason = f"non-catalytic metal type: {metal_type}"
        return _evaluation(False, 0, evidence_tags, "", excluded_reason)

    is_core_metal = metal_type in CORE_CATALYTIC_METALS
    is_conditional = metal_type in CONDITIONAL_METALS
    if not is_core_metal and not is_conditional:
        excluded_reason = f"metal type lacks catalytic prior: {metal_type or 'UNK'}"
        return _evaluation(False, 0, evidence_tags, "", excluded_reason)

    active_evidence = _has_active_site_evidence(site, active_sites)
    pdb_record_evidence = _has_pdb_record_evidence(site, pdb_record_residue_keys)
    anchor_evidence = _has_anchor_evidence(site, motif_anchor_atoms)
    donor_ok, donors = _donor_profile(site)
    role_text = str(site.get("functional_role", "")).lower()
    role_evidence = "catalytic" in role_text

    has_his_asp_glu_active = active_evidence and bool(res_set & {"HIS", "ASP", "GLU"})
    if metal_type in {"ZN", "ZN2"} and res_names.count("CYS") >= 3 and not has_his_asp_glu_active:
        excluded_reason = "Cys-rich structural zinc without active-site evidence"
        return _evaluation(False, 0, evidence_tags, "", excluded_reason)

    if is_core_metal:
        score += 20
        evidence_tags.append("catalytic-metal-prior")
    if role_evidence:
        score += 25
        evidence_tags.append("extractor-catalytic-role")
    if active_evidence:
        score += 40
        evidence_tags.append("uniprot-active-site")
    if pdb_record_evidence:
        score += 35
        evidence_tags.append("pdb-site-or-link-record")
    if anchor_evidence:
        score += 35
        evidence_tags.append("motif-anchor-proximity")
    if donor_ok:
        score += 15
        evidence_tags.append("N/O/S-coordination")

    strong_evidence = active_evidence or pdb_record_evidence or anchor_evidence or (
        is_core_metal and role_evidence
    )
    if is_conditional and not (active_evidence and strong_evidence):
        excluded_reason = f"{metal_type} requires strong active-site evidence"
        return _evaluation(False, score, evidence_tags, "", excluded_reason)

    if not donor_ok:
        excluded_reason = "no protein N/O/S donor coordination"
        return _evaluation(False, score, evidence_tags, "", excluded_reason)

    if not strong_evidence:
        excluded_reason = "isolated metal without active-site, PDB-record, or motif-anchor evidence"
        return _evaluation(False, score, evidence_tags, "", excluded_reason)

    included = score >= 55
    included_reason = "; ".join(evidence_tags) if included else ""
    if not included:
        excluded_reason = "insufficient catalytic evidence score"
    return _evaluation(included, score, evidence_tags, included_reason, excluded_reason)


def _evaluation(
    is_catalytic: bool,
    score: int,
    evidence_tags: List[str],
    included_reason: str,
    excluded_reason: str,
) -> Dict[str, Any]:
    if score >= 70:
        level = "strong"
    elif score >= 55:
        level = "moderate"
    elif score > 0:
        level = "weak"
    else:
        level = "none"
    return {
        "is_catalytic": bool(is_catalytic),
        "evidence_score": int(score),
        "evidence_level": level,
        "evidence_tags": evidence_tags,
        "included_reason": included_reason,
        "excluded_reason": excluded_reason,
    }


def _calculate_angles_degrees(
    metal_coords: Tuple[float, float, float],
    coordinating_residues: List[Dict[str, Any]],
) -> List[float]:
    coords = []
    for residue in coordinating_residues or []:
        c = _coords(residue.get("coordinates"))
        if c:
            coords.append(c)
    if len(coords) < 2:
        return []

    angles: List[float] = []
    mx, my, mz = metal_coords
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            v1 = [coords[i][0] - mx, coords[i][1] - my, coords[i][2] - mz]
            v2 = [coords[j][0] - mx, coords[j][1] - my, coords[j][2] - mz]
            n1 = math.sqrt(sum(x * x for x in v1))
            n2 = math.sqrt(sum(x * x for x in v2))
            if n1 == 0 or n2 == 0:
                continue
            cosang = sum(a * b for a, b in zip(v1, v2)) / (n1 * n2)
            cosang = max(-1.0, min(1.0, cosang))
            angles.append(round(math.degrees(math.acos(cosang)), 2))
    return angles


class CatalyticMetalDatabase:
    """Thread-safe SQLite access for high-precision catalytic metal sites."""

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self._local = threading.local()
        self._init_database()

    def _get_connection(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            self._local.conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
            self._local.conn.row_factory = sqlite3.Row
            try:
                self._local.conn.execute("PRAGMA journal_mode=WAL")
                self._local.conn.execute("PRAGMA busy_timeout=5000")
                self._local.conn.execute("PRAGMA synchronous=NORMAL")
            except sqlite3.OperationalError:
                pass
        return self._local.conn

    def close(self) -> None:
        """Close the current thread's cached SQLite connection."""
        if hasattr(self._local, "conn") and self._local.conn:
            self._local.conn.close()
            self._local.conn = None

    def _init_database(self) -> None:
        conn = self._get_connection()
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS catalytic_metal_site (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                nanozyme_type TEXT NOT NULL,
                ec_number TEXT NOT NULL,
                pdb_id TEXT NOT NULL,
                pdb_path TEXT NOT NULL,
                metal_type TEXT NOT NULL,
                metal_name TEXT NOT NULL,
                metal_chain TEXT,
                metal_residue_id INTEGER,
                metal_coords_json TEXT,
                coordination_number INTEGER,
                coordination_geometry TEXT,
                oxidation_state TEXT,
                coordinating_residues_json TEXT,
                coordination_distances_json TEXT,
                coordination_angles_json TEXT,
                evidence_score INTEGER,
                evidence_level TEXT,
                evidence_tags_json TEXT,
                included_reason TEXT,
                excluded_reason TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        cur.execute("CREATE INDEX IF NOT EXISTS idx_cat_metal_nanozyme ON catalytic_metal_site(nanozyme_type)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_cat_metal_ec ON catalytic_metal_site(ec_number)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_cat_metal_type ON catalytic_metal_site(metal_type)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_cat_metal_pdb ON catalytic_metal_site(pdb_id)")
        cur.execute(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS idx_cat_metal_unique_site
            ON catalytic_metal_site(
                nanozyme_type, ec_number, pdb_id, metal_type,
                metal_residue_id
            )
            """
        )
        conn.commit()

    def clear(self) -> None:
        conn = self._get_connection()
        conn.execute("DELETE FROM catalytic_metal_site")
        conn.commit()

    def commit(self) -> None:
        self._get_connection().commit()

    def add_site(
        self,
        nanozyme_type: str,
        ec_number: str,
        pdb_id: str,
        pdb_path: str,
        site: Dict[str, Any],
        angles: List[float],
        evaluation: Dict[str, Any],
        *,
        commit: bool = True,
        existing_coords: Optional[List[List[float]]] = None,
    ) -> bool:
        conn = self._get_connection()
        # 对称链去重：同 pdb+metal_type 已有坐标距离 < 2.0 Å 的位点则跳过
        coords = site.get("metal_coords")
        if coords:
            import math
            if existing_coords is None:
                existing = conn.execute(
                    "SELECT metal_coords_json FROM catalytic_metal_site "
                    "WHERE nanozyme_type=? AND ec_number=? AND pdb_id=? AND metal_type=?",
                    (nanozyme_type, ec_number, pdb_id, site.get("metal_type", ""))
                ).fetchall()
                coords_to_check = [_json_loads(cj) for (cj,) in existing]
            else:
                coords_to_check = existing_coords
            for c2 in coords_to_check:
                if c2 and math.dist(coords[:3], c2[:3]) < 2.0:
                    return False
        cur = conn.execute(
            """
            INSERT OR IGNORE INTO catalytic_metal_site (
                nanozyme_type, ec_number, pdb_id, pdb_path,
                metal_type, metal_name, metal_chain, metal_residue_id, metal_coords_json,
                coordination_number, coordination_geometry, oxidation_state,
                coordinating_residues_json, coordination_distances_json, coordination_angles_json,
                evidence_score, evidence_level, evidence_tags_json, included_reason, excluded_reason
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                nanozyme_type,
                ec_number,
                pdb_id,
                pdb_path,
                site.get("metal_type", "") or "UNK",
                site.get("metal_name", "") or "Unknown",
                site.get("metal_chain", ""),
                site.get("metal_residue_id", None),
                _json_dumps(site.get("metal_coords")),
                int(site.get("coordination_number") or 0),
                site.get("coordination_geometry", "") or "unknown",
                site.get("oxidation_state", None),
                _json_dumps(site.get("coordinating_residues") or []),
                _json_dumps(site.get("coordination_distances") or []),
                _json_dumps(angles or []),
                int(evaluation.get("evidence_score") or 0),
                evaluation.get("evidence_level", "none"),
                _json_dumps(evaluation.get("evidence_tags") or []),
                evaluation.get("included_reason", ""),
                evaluation.get("excluded_reason", ""),
            ),
        )
        inserted = cur.rowcount > 0
        if inserted and coords and existing_coords is not None:
            existing_coords.append(list(coords[:3]))
        if commit:
            conn.commit()
        return inserted

    def get_metal_sites_summary_by_nanozyme_type(self, nanozyme_type: str) -> List[Dict[str, Any]]:
        if not nanozyme_type:
            return []
        conn = self._get_connection()
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                metal_type,
                MIN(metal_name) AS metal_name,
                COUNT(*) AS occurrence_count,
                GROUP_CONCAT(DISTINCT pdb_id) AS pdb_ids,
                GROUP_CONCAT(DISTINCT pdb_path) AS pdb_paths,
                MAX(evidence_score) AS evidence_score,
                MIN(evidence_level) AS evidence_level,
                MIN(included_reason) AS included_reason,
                MIN(coordination_number) AS coordination_number,
                MIN(coordination_geometry) AS coordination_geometry
            FROM catalytic_metal_site
            WHERE lower(nanozyme_type) = lower(?)
            GROUP BY metal_type
            ORDER BY occurrence_count DESC, evidence_score DESC
            """,
            (nanozyme_type,),
        )
        rows = cur.fetchall()

        result: List[Dict[str, Any]] = []
        for row in rows:
            metal_type = row["metal_type"]
            cur.execute(
                """
                SELECT coordinating_residues_json
                FROM catalytic_metal_site
                WHERE lower(nanozyme_type)=lower(?) AND metal_type=?
                ORDER BY evidence_score DESC
                LIMIT 1
                """,
                (nanozyme_type, metal_type),
            )
            rep = cur.fetchone()
            coords = _json_loads(rep[0]) if rep and rep[0] else []
            pdb_ids = sorted(set(p.upper() for p in (row["pdb_ids"] or "").split(",") if p))
            pdb_paths = sorted(set(p for p in (row["pdb_paths"] or "").split(",") if p))
            occurrence_count = int(row["occurrence_count"] or 0)
            metal_id = f"{metal_type}_metal_{occurrence_count}occ"
            result.append(
                {
                    "metal_type": metal_type,
                    "metal_name": row["metal_name"],
                    "ligand_name": row["metal_name"] or metal_type,
                    "ligand_id": f"{metal_type}_metal",
                    "metal_id": metal_id,
                    "motif_id": metal_id,
                    "anchor_atoms_count": 0,
                    "occurrence_count": occurrence_count,
                    "pdb_ids": pdb_ids,
                    "file_paths": pdb_paths,
                    "atom_count": 1,
                    "category": "metal_sites",
                    "source": "catalytic_metal_index",
                    "functional_role": "catalytic",
                    "coordinating_residues": coords,
                    "coordination_number": row["coordination_number"],
                    "coordination_geometry": row["coordination_geometry"],
                    "evidence_score": row["evidence_score"],
                    "evidence_level": row["evidence_level"],
                    "included_reason": row["included_reason"],
                }
            )
        return result

    def get_sites_for_metal(self, nanozyme_type: str, metal_type: str, limit: int = 50) -> List[Dict[str, Any]]:
        conn = self._get_connection()
        cur = conn.cursor()
        cur.execute(
            """
            SELECT * FROM catalytic_metal_site
            WHERE lower(nanozyme_type) = lower(?) AND upper(metal_type) = upper(?)
            ORDER BY evidence_score DESC, pdb_id
            LIMIT ?
            """,
            (nanozyme_type, metal_type, int(limit)),
        )
        return [self._row_to_dict(row) for row in cur.fetchall()]

    def get_any_site_for_metal(self, metal_type: str) -> Optional[Dict[str, Any]]:
        conn = self._get_connection()
        cur = conn.cursor()
        cur.execute(
            """
            SELECT * FROM catalytic_metal_site
            WHERE upper(metal_type) = upper(?)
            ORDER BY evidence_score DESC, pdb_id
            LIMIT 1
            """,
            (metal_type,),
        )
        row = cur.fetchone()
        return self._row_to_dict(row) if row else None

    @staticmethod
    def _row_to_dict(row: sqlite3.Row) -> Dict[str, Any]:
        data = dict(row)
        data["metal_coords"] = _json_loads(data.pop("metal_coords_json", "") or "")
        data["coordinating_residues"] = _json_loads(data.pop("coordinating_residues_json", "") or "")
        data["coordination_distances"] = _json_loads(data.pop("coordination_distances_json", "") or "")
        data["coordination_angles"] = _json_loads(data.pop("coordination_angles_json", "") or "")
        data["evidence_tags"] = _json_loads(data.pop("evidence_tags_json", "") or "")
        return data


def _ec_dir_to_ec_number(ec_dir_name: str) -> Optional[str]:
    parts = ec_dir_name.split("_")
    if len(parts) != 4:
        return None
    try:
        [int(part) for part in parts]
    except ValueError:
        return None
    return ec_dir_name.replace("_", ".")


def _load_library_index(pdb_library_dir: Path) -> Dict[str, str]:
    index_path = pdb_library_dir / "library_index.json"
    if not index_path.exists():
        return {}
    with open(index_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    mapping: Dict[str, str] = {}
    for ec_number, entry in (data.get("ec_entries") or {}).items():
        if isinstance(entry, dict) and entry.get("nanozyme_type"):
            mapping[str(ec_number)] = str(entry.get("nanozyme_type"))
    return mapping


def _load_sites_json(ec_dir: Path, ec_number: str) -> Dict[str, Dict[str, Any]]:
    path = ec_dir / f"{ec_number}_sites.json"
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        entries = json.load(handle)
    out: Dict[str, Dict[str, Any]] = {}
    for entry in entries or []:
        pdb_id = str(entry.get("pdb_id", "")).upper().strip()
        if pdb_id:
            out[pdb_id] = entry
    return out


def _load_motif_anchors(motif_library_dir: Optional[Path]) -> Dict[str, List[Dict[str, Any]]]:
    if not motif_library_dir or not motif_library_dir.exists():
        return {}
    anchors_by_pdb: Dict[str, List[Dict[str, Any]]] = {}
    for motif_file in motif_library_dir.glob("*/*/*.json"):
        try:
            with open(motif_file, "r", encoding="utf-8") as handle:
                data = json.load(handle)
        except Exception:
            continue
        pdb_id = str(data.get("source_pdb_id", "")).upper().strip()
        if not pdb_id:
            continue
        anchors_by_pdb.setdefault(pdb_id, []).extend(data.get("anchor_atoms") or [])
    return anchors_by_pdb


def _extract_pdb_record_context_from_lines(lines) -> Tuple[set, List[str]]:
    """Extract reusable PDB record context for catalytic-evidence checks."""
    site_keys = set()
    evidence_lines: List[str] = []
    for line in lines:
        if line.startswith("SITE"):
            for i in range(4):
                start = 18 + i * 11
                res_name = line[start:start + 3].strip()
                chain = line[start + 4:start + 5].strip()
                res_num = line[start + 5:start + 10].strip()
                if not res_name or not res_num:
                    continue
                try:
                    num = int(res_num)
                except ValueError:
                    continue
                site_keys.add((chain, num))
                site_keys.add(("", num))
        elif line.startswith(("LINK", "CONECT", "REMARK 620", "REMARK 800")):
            evidence_lines.append(line.upper())
    return site_keys, evidence_lines


def _extract_pdb_record_context(pdb_path: Path) -> Tuple[set, List[str]]:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            return _extract_pdb_record_context_from_lines(handle)
    except Exception:
        return set(), []


def _pdb_record_residue_keys_from_context(
    site: Dict[str, Any],
    context: Tuple[set, List[str]],
) -> set:
    """Project reusable PDB record context onto one metal site."""
    keys = set(context[0])
    evidence_lines = context[1]
    metal_type = _normalize_metal_type(site.get("metal_type"))
    metal_residue_id = site.get("metal_residue_id")
    try:
        metal_residue_id = int(metal_residue_id)
    except (TypeError, ValueError):
        metal_residue_id = None
    for line in evidence_lines:
        if metal_type and metal_type not in line:
            continue
        if metal_residue_id is not None:
            keys.add(("", metal_residue_id))
    return keys


def _extract_pdb_record_residue_keys(pdb_path: Path, site: Dict[str, Any]) -> set:
    """Extract a conservative residue-key set from SITE/LINK/CONECT/REMARK records."""
    return _pdb_record_residue_keys_from_context(site, _extract_pdb_record_context(pdb_path))


def build_catalytic_metal_db(
    pdb_library_dir: Path,
    out_db_path: Path,
    motif_library_dir: Optional[Path] = None,
    clear_existing: bool = True,
    limit_pdb_files: Optional[int] = None,
) -> None:
    from nanozyme_mining.structure.pdb_metal_extractor import PDBMetalExtractor, MetalSite

    pdb_library_dir = Path(pdb_library_dir)
    out_db_path = Path(out_db_path)
    out_db_path.parent.mkdir(parents=True, exist_ok=True)

    ec_to_type = _load_library_index(pdb_library_dir)
    anchors_by_pdb = _load_motif_anchors(Path(motif_library_dir) if motif_library_dir else None)
    extractor = PDBMetalExtractor()
    db = CatalyticMetalDatabase(out_db_path)
    if clear_existing:
        db.clear()

    seen_coords_by_site: Dict[Tuple[str, str, str, str], List[List[float]]] = {}
    pdb_files: List[Tuple[str, str, Path, Dict[str, Any]]] = []
    for ec_dir in sorted(pdb_library_dir.iterdir()):
        if not ec_dir.is_dir():
            continue
        ec_number = _ec_dir_to_ec_number(ec_dir.name)
        if not ec_number:
            continue
        nanozyme_type = ec_to_type.get(ec_number, "") or ""
        if not nanozyme_type:
            continue
        site_entries = _load_sites_json(ec_dir, ec_number)
        for pdb_path in sorted(ec_dir.glob("*.pdb")):
            pdb_id = pdb_path.stem.upper()
            pdb_files.append((nanozyme_type, ec_number, pdb_path, site_entries.get(pdb_id, {})))

    if limit_pdb_files is not None:
        pdb_files = pdb_files[: int(limit_pdb_files)]

    total = len(pdb_files)
    included_count = 0
    examined_count = 0
    for idx, (nanozyme_type, ec_number, pdb_path, site_entry) in enumerate(pdb_files, start=1):
        try:
            with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
                pdb_lines = handle.readlines()
            metal_sites: List[MetalSite] = extractor.parse_pdb_lines(pdb_lines)
        except Exception:
            continue
        if not metal_sites:
            continue

        pdb_id = pdb_path.stem.upper()
        active_sites = site_entry.get("active_sites", []) if isinstance(site_entry, dict) else []
        motif_anchor_atoms = anchors_by_pdb.get(pdb_id, [])
        pdb_record_context = _extract_pdb_record_context_from_lines(pdb_lines)

        for site_obj in metal_sites:
            examined_count += 1
            site = site_obj.to_dict()
            metal_coords = _coords(site.get("metal_coords")) or [0.0, 0.0, 0.0]
            site["metal_coords"] = metal_coords
            pdb_record_keys = _pdb_record_residue_keys_from_context(site, pdb_record_context)
            evaluation = evaluate_catalytic_metal_site(
                site,
                ec_number=ec_number,
                nanozyme_type=nanozyme_type,
                active_sites=active_sites,
                pdb_record_residue_keys=pdb_record_keys,
                motif_anchor_atoms=motif_anchor_atoms,
            )
            if not evaluation["is_catalytic"]:
                continue
            angles = _calculate_angles_degrees(tuple(metal_coords), site.get("coordinating_residues") or [])
            coords_key = (
                nanozyme_type,
                ec_number,
                pdb_id,
                str(site.get("metal_type", "")),
            )
            existing_coords = seen_coords_by_site.setdefault(coords_key, []) if clear_existing else None
            inserted = db.add_site(
                nanozyme_type=nanozyme_type,
                ec_number=ec_number,
                pdb_id=pdb_id,
                pdb_path=str(pdb_path),
                site=site,
                angles=angles,
                evaluation=evaluation,
                commit=False,
                existing_coords=existing_coords,
            )
            if inserted:
                included_count += 1

        if idx % 50 == 0:
            db.commit()
            print(f"[catalytic-metal] processed {idx}/{total} PDB files...")

    db.commit()

    print(
        f"[catalytic-metal] done. db={out_db_path}; "
        f"examined={examined_count}; included={included_count}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build high-precision catalytic metal SQLite index from pdb_library")
    parser.add_argument("--pdb-library", required=True, help="Path to pdb_library directory")
    parser.add_argument("--out", required=True, help="Output SQLite DB path")
    parser.add_argument("--motif-library", default=None, help="Optional motif_library path for anchor evidence")
    parser.add_argument("--no-clear", action="store_true", help="Do not clear existing rows before building")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of PDB files (debug)")
    args = parser.parse_args()

    build_catalytic_metal_db(
        pdb_library_dir=Path(args.pdb_library),
        out_db_path=Path(args.out),
        motif_library_dir=Path(args.motif_library) if args.motif_library else None,
        clear_existing=not args.no_clear,
        limit_pdb_files=args.limit,
    )


if __name__ == "__main__":
    main()
