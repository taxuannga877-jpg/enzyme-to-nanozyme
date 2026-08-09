"""
Ligand / cofactor library (SQLite, offline)
===========================================

Purpose
-------
- Extract ligands/cofactors from local PDB files under `pdb_library`.
- Store lightweight metadata for fast listing in `/api/list_motifs` without
  reparsing PDBs on each request.

Build
-----
python -m enzyme_viewer.ligand_db --pdb-library /abs/path/pdb_library --out /abs/path/enzyme_viewer/ligand_index.db
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import threading
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser
from nanozyme_mining.utils.constants import KNOWN_METAL_ELEMENTS, METAL_RESIDUE_NAMES


def _json_dumps(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def _json_loads(s: str) -> Any:
    return json.loads(s) if s else None


def _ec_dir_to_ec_number(ec_dir_name: str) -> Optional[str]:
    parts = ec_dir_name.split("_")
    if len(parts) != 4:
        return None
    try:
        [int(p) for p in parts]
    except ValueError:
        return None
    return ec_dir_name.replace("_", ".")


def _load_library_index(pdb_library_dir: Path) -> Dict[str, str]:
    """
    Returns mapping: ec_number -> nanozyme_type
    """
    index_path = pdb_library_dir / "library_index.json"
    if not index_path.exists():
        return {}
    with open(index_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    ec_entries = data.get("ec_entries", {}) or {}
    mapping: Dict[str, str] = {}
    for ec_number, entry in ec_entries.items():
        if not isinstance(entry, dict):
            continue
        nanozyme_type = entry.get("nanozyme_type", "") or ""
        if ec_number and nanozyme_type:
            mapping[str(ec_number)] = str(nanozyme_type)
    return mapping


def _is_water(res_name: str) -> bool:
    return res_name.strip().upper() == "HOH"


def _is_metal_ligand(res_name: str, element: str) -> bool:
    """Skip metals (handled by catalytic_metal_db / motif_db).

    PR2 (NEW-8 / v4 audit): docstring previously referenced 'mental_db' (a typo
    for metal_db that was the legacy precursor to catalytic_metal_index.db).
    Updated to the current module names; behaviour unchanged.
    """
    res_upper = res_name.strip().upper()
    elem_upper = element.strip().upper()
    if elem_upper in KNOWN_METAL_ELEMENTS:
        return True
    return res_upper in METAL_RESIDUE_NAMES


class LigandDatabase:
    """Thread-safe SQLite access for ligand/cofactor library."""

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
            CREATE TABLE IF NOT EXISTS ligand_entry (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                nanozyme_type TEXT NOT NULL,
                ec_number TEXT NOT NULL,
                pdb_id TEXT NOT NULL,
                pdb_path TEXT NOT NULL,
                ligand_name TEXT NOT NULL,
                chain TEXT,
                residue_number INTEGER,
                atom_count INTEGER,
                het_id TEXT,
                het_names_json TEXT,
                het_synonyms_json TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        cur.execute("CREATE INDEX IF NOT EXISTS idx_ligand_nanozyme ON ligand_entry(nanozyme_type)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_ligand_ec ON ligand_entry(ec_number)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_ligand_name ON ligand_entry(ligand_name)")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_ligand_pdb ON ligand_entry(pdb_id)")
        conn.commit()

    def clear(self) -> None:
        conn = self._get_connection()
        conn.execute("DELETE FROM ligand_entry")
        conn.commit()

    def add_ligand(
        self,
        nanozyme_type: str,
        ec_number: str,
        pdb_id: str,
        pdb_path: str,
        ligand_name: str,
        chain: str,
        residue_number: Optional[int],
        atom_count: int,
        het_id: Optional[str],
        het_names: Iterable[str],
        het_synonyms: Iterable[str],
    ) -> None:
        conn = self._get_connection()
        conn.execute(
            """
            INSERT INTO ligand_entry (
                nanozyme_type, ec_number, pdb_id, pdb_path,
                ligand_name, chain, residue_number, atom_count,
                het_id, het_names_json, het_synonyms_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                nanozyme_type,
                ec_number,
                pdb_id,
                pdb_path,
                ligand_name,
                chain,
                residue_number,
                atom_count,
                het_id,
                _json_dumps(list(het_names) if het_names else []),
                _json_dumps(list(het_synonyms) if het_synonyms else []),
            ),
        )
        conn.commit()

    def get_ligands_summary_by_nanozyme_type(self, nanozyme_type: str) -> List[Dict[str, Any]]:
        """
        Return compact ligand/cofactor summary for the UI.
        Grouped by ligand_name (residue name).
        """
        if not nanozyme_type:
            return []
        conn = self._get_connection()
        cur = conn.cursor()
        # IMPORTANT:
        # - The UI needs a *consistent* instance (pdb_path + chain + residue_number) that actually exists together.
        # - Using MIN(chain)/MIN(residue_number)/MIN(pdb_path) independently can create impossible combinations,
        #   causing "Ligand XXX not found in PDB file" when fetching structure.
        cur.execute(
            """
            WITH filtered AS (
                SELECT
                    id,
                    UPPER(ligand_name) AS ligand_key,
                    ligand_name,
                    pdb_id,
                    pdb_path,
                    chain,
                    residue_number,
                    atom_count
                FROM ligand_entry
                WHERE lower(nanozyme_type) = lower(?)
            ),
            agg AS (
                SELECT
                    ligand_key,
                    MIN(id) AS rep_id,
                    COUNT(*) AS occurrence_count,
                    GROUP_CONCAT(DISTINCT pdb_id) AS pdb_ids,
                    GROUP_CONCAT(DISTINCT pdb_path) AS pdb_paths
                FROM filtered
                GROUP BY ligand_key
            )
            SELECT
                a.ligand_key AS ligand_key,
                f.ligand_name AS ligand_name,
                a.occurrence_count AS occurrence_count,
                a.pdb_ids AS pdb_ids,
                a.pdb_paths AS pdb_paths,
                f.pdb_id AS rep_pdb_id,
                f.pdb_path AS rep_pdb_path,
                f.chain AS rep_chain,
                f.residue_number AS rep_residue_number,
                f.atom_count AS rep_atom_count
            FROM agg a
            JOIN filtered f ON f.id = a.rep_id
            ORDER BY a.occurrence_count DESC
            """,
            (nanozyme_type,),
        )
        rows = cur.fetchall()
        result: List[Dict[str, Any]] = []
        for r in rows:
            pdb_ids = [x for x in (r["pdb_ids"] or "").split(",") if x]
            pdb_paths = [x for x in (r["pdb_paths"] or "").split(",") if x]
            ligand_name = r["ligand_name"] or r["ligand_key"] or "UNK"
            ligand_id = f"{ligand_name}_ligand"
            result.append(
                {
                    "ligand_name": ligand_name,
                    "ligand_id": ligand_id,
                    # Use representative instance fields (same row), not aggregated mins.
                    "pdb_id": (r["rep_pdb_id"] or (pdb_ids[0] if pdb_ids else "")),
                    "pdb_ids": sorted(set([p.upper() for p in pdb_ids])),
                    "file_path": (r["rep_pdb_path"] or (pdb_paths[0] if pdb_paths else "")),
                    "file_paths": sorted(set(pdb_paths)),
                    "chain": r["rep_chain"] or "",
                    "residue_number": r["rep_residue_number"],
                    "atom_count": int(r["rep_atom_count"] or 0),
                    "occurrence_count": int(r["occurrence_count"] or 0),
                    "category": "ligands_cofactors",
                    "source": "@ligand",
                    "nanozyme_type": nanozyme_type,
                    "ec_number": None,
                    "uniprot_id": "",
                }
            )
        return result


def _extract_ligands_from_pdb(pdb_path: Path) -> List[Dict[str, Any]]:
    """
    Parse a single PDB and return grouped ligand entries (excluding water/metals).
    """
    parser = ComprehensivePDBParser()
    parsed = parser.parse_pdb_file(pdb_path)
    het_info = parsed.get("het_info") or {}
    grouped = parser.extract_ligands_and_cofactors(parsed, het_info)
    ligands: List[Dict[str, Any]] = []
    for group in grouped:
        res_name = (group.get("residue_name") or "").strip()
        if not res_name or _is_water(res_name):
            continue
        atoms = group.get("atoms") or []
        if not atoms:
            continue
        # Use first atom to check element for metal screening
        first_atom = atoms[0] if atoms else {}
        element = (first_atom.get("element") or "").strip()
        if _is_metal_ligand(res_name, element):
            continue
        ligands.append(
            {
                "ligand_name": res_name,
                "chain": group.get("chain") or "",
                "residue_number": group.get("residue_number"),
                "atom_count": len(atoms),
                "het_id": group.get("het_info", {}).get("het_id") if group.get("het_info") else res_name,
                "het_names": (group.get("het_info") or {}).get("names") or [],
                "het_synonyms": (group.get("het_info") or {}).get("synonyms") or [],
            }
        )
    return ligands


def build_ligand_db(
    pdb_library_dir: Path,
    out_db_path: Path,
    clear_existing: bool = True,
    limit_pdb_files: Optional[int] = None,
) -> None:
    pdb_library_dir = Path(pdb_library_dir)
    out_db_path = Path(out_db_path)
    out_db_path.parent.mkdir(parents=True, exist_ok=True)

    ec_to_type = _load_library_index(pdb_library_dir)
    db = LigandDatabase(out_db_path)
    if clear_existing:
        db.clear()

    pdb_files: List[Tuple[str, str, Path]] = []
    for ec_dir in sorted(pdb_library_dir.iterdir()):
        if not ec_dir.is_dir():
            continue
        ec_number = _ec_dir_to_ec_number(ec_dir.name)
        if not ec_number:
            continue
        nanozyme_type = ec_to_type.get(ec_number, "") or ""
        if not nanozyme_type:
            continue
        for pdb_path in sorted(ec_dir.glob("*.pdb")):
            pdb_id = pdb_path.stem.upper()
            pdb_files.append((nanozyme_type, ec_number, pdb_path))

    if limit_pdb_files is not None:
        pdb_files = pdb_files[: int(limit_pdb_files)]

    total = len(pdb_files)
    for idx, (nanozyme_type, ec_number, pdb_path) in enumerate(pdb_files, start=1):
        try:
            ligands = _extract_ligands_from_pdb(pdb_path)
        except Exception:
            continue
        if not ligands:
            continue
        for lig in ligands:
            db.add_ligand(
                nanozyme_type=nanozyme_type,
                ec_number=ec_number,
                pdb_id=pdb_path.stem.upper(),
                pdb_path=str(pdb_path),
                ligand_name=lig.get("ligand_name", "UNK"),
                chain=lig.get("chain", ""),
                residue_number=lig.get("residue_number"),
                atom_count=int(lig.get("atom_count") or 0),
                het_id=lig.get("het_id"),
                het_names=lig.get("het_names") or [],
                het_synonyms=lig.get("het_synonyms") or [],
            )

        if idx % 50 == 0:
            print(f"[ligand] processed {idx}/{total} PDB files...")

    print(f"[ligand] done. db={out_db_path} (rows stored in ligand_entry table)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build ligand/cofactor SQLite index from pdb_library")
    parser.add_argument("--pdb-library", required=True, help="Path to pdb_library directory")
    parser.add_argument("--out", required=True, help="Output SQLite DB path")
    parser.add_argument("--no-clear", action="store_true", help="Do not clear existing rows before building")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of PDB files (debug)")
    args = parser.parse_args()

    build_ligand_db(
        pdb_library_dir=Path(args.pdb_library),
        out_db_path=Path(args.out),
        clear_existing=not args.no_clear,
        limit_pdb_files=args.limit,
    )


if __name__ == "__main__":
    main()
