"""
Nanozyme Database - Stage 1 Core Module
========================================

Manages EC number to nanozyme function type mapping database.
Based on ChemEnzyRetroPlanner's UniProtParserEC architecture.
"""

import os
import json
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from datetime import datetime

from ..utils.constants import NanozymeType, EC_TO_NANOZYME_TYPE
from ..utils.ec_mappings import EC_PATTERNS


@dataclass
class EnzymeEntry:
    """Single enzyme entry in the database."""
    uniprot_id: str
    ec_number: str
    nanozyme_type: str
    sequence: str
    sequence_length: int
    pdb_id: Optional[str] = None
    pdb_path: Optional[str] = None
    active_sites: Optional[str] = None  # JSON string
    created_at: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EnzymeEntry":
        return cls(**data)


class NanozymeDatabase:
    """
    Stage 1: EC -> Nanozyme Function Type Mapping Database

    Manages a SQLite database that maps EC numbers to nanozyme function types
    and organizes structure entries for downstream processing.

    Based on ChemEnzyRetroPlanner's data organization patterns.
    """

    def __init__(self, db_path: str = "nanozyme_db.sqlite"):
        """
        Initialize the nanozyme database.

        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_database()

    def _init_database(self):
        """Initialize database schema."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Create enzymes table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS enzymes (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                uniprot_id TEXT UNIQUE NOT NULL,
                ec_number TEXT NOT NULL,
                nanozyme_type TEXT NOT NULL,
                sequence TEXT,
                sequence_length INTEGER,
                pdb_id TEXT,
                pdb_path TEXT,
                active_sites TEXT,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Create EC mapping table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS ec_mappings (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ec_number TEXT UNIQUE NOT NULL,
                nanozyme_type TEXT NOT NULL,
                description TEXT,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Create indexes
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_ec_number
            ON enzymes(ec_number)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_nanozyme_type
            ON enzymes(nanozyme_type)
        """)

        conn.commit()
        conn.close()

        # Initialize EC mappings
        self._init_ec_mappings()

    def _init_ec_mappings(self):
        """Initialize EC number to nanozyme type mappings."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Insert predefined EC mappings
        for nanozyme_type, ec_list in EC_PATTERNS.items():
            for ec_number in ec_list:
                try:
                    cursor.execute("""
                        INSERT OR IGNORE INTO ec_mappings
                        (ec_number, nanozyme_type, description)
                        VALUES (?, ?, ?)
                    """, (ec_number, nanozyme_type.value, nanozyme_type.name))
                except sqlite3.IntegrityError:
                    pass

        conn.commit()
        conn.close()

    def add_enzyme(self, entry: EnzymeEntry) -> bool:
        """
        Add an enzyme entry to the database.

        Args:
            entry: EnzymeEntry object

        Returns:
            True if successful, False otherwise
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO enzymes
                (uniprot_id, ec_number, nanozyme_type, sequence,
                 sequence_length, pdb_id, pdb_path, active_sites)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                entry.uniprot_id,
                entry.ec_number,
                entry.nanozyme_type,
                entry.sequence,
                entry.sequence_length,
                entry.pdb_id,
                entry.pdb_path,
                entry.active_sites
            ))
            conn.commit()
            return True
        except Exception as e:
            print(f"Error adding enzyme: {e}")
            return False
        finally:
            conn.close()

    def add_enzymes_batch(self, entries: List[EnzymeEntry]) -> int:
        """
        Add multiple enzyme entries in batch.

        Args:
            entries: List of EnzymeEntry objects

        Returns:
            Number of successfully added entries
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        count = 0

        for entry in entries:
            try:
                cursor.execute("""
                    INSERT OR REPLACE INTO enzymes
                    (uniprot_id, ec_number, nanozyme_type, sequence,
                     sequence_length, pdb_id, pdb_path, active_sites)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    entry.uniprot_id,
                    entry.ec_number,
                    entry.nanozyme_type,
                    entry.sequence,
                    entry.sequence_length,
                    entry.pdb_id,
                    entry.pdb_path,
                    entry.active_sites
                ))
                count += 1
            except Exception as e:
                print(f"Error adding enzyme {entry.uniprot_id}: {e}")

        conn.commit()
        conn.close()
        return count

    def get_nanozyme_type(self, ec_number: str) -> Optional[NanozymeType]:
        """
        Get nanozyme type for an EC number.

        Args:
            ec_number: EC number string (e.g., "1.11.1.7")

        Returns:
            NanozymeType enum or None if not found
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT nanozyme_type FROM ec_mappings
            WHERE ec_number = ?
        """, (ec_number,))

        result = cursor.fetchone()
        conn.close()

        if result:
            try:
                return NanozymeType(result[0])
            except ValueError:
                return NanozymeType.UNKNOWN
        return None

    def query_by_ec(self, ec_number: str) -> List[EnzymeEntry]:
        """
        Query enzymes by EC number.

        Args:
            ec_number: EC number string

        Returns:
            List of EnzymeEntry objects
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT uniprot_id, ec_number, nanozyme_type, sequence,
                   sequence_length, pdb_id, pdb_path, active_sites
            FROM enzymes WHERE ec_number = ?
        """, (ec_number,))

        results = cursor.fetchall()
        conn.close()

        return [
            EnzymeEntry(
                uniprot_id=r[0],
                ec_number=r[1],
                nanozyme_type=r[2],
                sequence=r[3],
                sequence_length=r[4],
                pdb_id=r[5],
                pdb_path=r[6],
                active_sites=r[7]
            ) for r in results
        ]

    def query_by_nanozyme_type(
        self, nanozyme_type: NanozymeType
    ) -> List[EnzymeEntry]:
        """
        Query enzymes by nanozyme type.

        Args:
            nanozyme_type: NanozymeType enum

        Returns:
            List of EnzymeEntry objects
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT uniprot_id, ec_number, nanozyme_type, sequence,
                   sequence_length, pdb_id, pdb_path, active_sites
            FROM enzymes WHERE nanozyme_type = ?
        """, (nanozyme_type.value,))

        results = cursor.fetchall()
        conn.close()

        return [
            EnzymeEntry(
                uniprot_id=r[0],
                ec_number=r[1],
                nanozyme_type=r[2],
                sequence=r[3],
                sequence_length=r[4],
                pdb_id=r[5],
                pdb_path=r[6],
                active_sites=r[7]
            ) for r in results
        ]

    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("SELECT COUNT(*) FROM enzymes")
        total = cursor.fetchone()[0]

        cursor.execute("""
            SELECT nanozyme_type, COUNT(*)
            FROM enzymes GROUP BY nanozyme_type
        """)
        by_type = dict(cursor.fetchall())

        conn.close()

        return {
            "total_enzymes": total,
            "by_nanozyme_type": by_type
        }
