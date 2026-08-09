"""
Motif数据库索引模块
使用SQLite存储motif的索引信息，加速查询
"""

import json
import math
import sqlite3
import threading
from pathlib import Path
from typing import Any, List, Dict, Optional
from dataclasses import dataclass

from nanozyme_mining.utils.constants import KNOWN_METAL_ELEMENTS


@dataclass
class MotifIndex:
    """Motif索引数据类"""
    motif_id: str
    uniprot_id: str
    source_pdb_id: str
    ec_number: str
    nanozyme_type: str
    category: str  # metal_sites, catalytic_sites, binding_sites, other
    anchor_atoms_count: int
    file_path: str
    confidence_score: float = 0.0
    chemistry_tag: str = ""
    reaction_smiles: str = ""


class MotifDatabase:
    """Motif数据库管理类（线程安全）"""

    def __init__(self, db_path: Path):
        """
        初始化数据库

        Args:
            db_path: SQLite数据库文件路径
        """
        self.db_path = Path(db_path)
        # 使用线程本地存储，为每个线程创建独立的连接
        self._local = threading.local()
        self._init_database()

    def _get_connection(self):
        """获取当前线程的数据库连接（线程安全）"""
        if not hasattr(self._local, 'conn') or self._local.conn is None:
            # 为当前线程创建新连接
            self._local.conn = sqlite3.connect(
                str(self.db_path),
                # PR3 (NEW-4 docstring fix): the previous comment said
                # "allow use in different threads", which made the intent
                # sound like cross-thread sharing — that would be unsafe.
                # The real purpose: this class uses threading.local() above
                # so each thread gets its OWN connection. sqlite3's default
                # check_same_thread=True would raise ProgrammingError when
                # Flask's threaded=True dev server hands a request off to a
                # different worker thread, even though threading.local
                # ensures that thread has its own conn. Setting False here
                # is the standard pattern for "per-thread cached conn".
                check_same_thread=False,
            )
            self._local.conn.row_factory = sqlite3.Row  # 返回字典格式
            # PR1-1 (M42 fix): enable WAL + reasonable busy_timeout for concurrent
            # reads/writes. WAL lets readers run while a writer holds the lock.
            try:
                self._local.conn.execute("PRAGMA journal_mode=WAL")
                self._local.conn.execute("PRAGMA busy_timeout=5000")
                self._local.conn.execute("PRAGMA synchronous=NORMAL")
            except sqlite3.OperationalError:
                # WAL may be unavailable on some filesystems (e.g. network shares);
                # gracefully degrade to default journal mode.
                pass
        return self._local.conn

    def _init_database(self):
        """初始化数据库表结构"""
        conn = self._get_connection()
        cursor = conn.cursor()

        # 创建motif索引表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS motif_index (
                motif_id TEXT PRIMARY KEY,
                uniprot_id TEXT NOT NULL,
                source_pdb_id TEXT DEFAULT '',
                ec_number TEXT NOT NULL,
                nanozyme_type TEXT NOT NULL,
                category TEXT NOT NULL,
                anchor_atoms_count INTEGER NOT NULL,
                file_path TEXT NOT NULL,
                confidence_score REAL DEFAULT 0.0,
                chemistry_tag TEXT DEFAULT '',
                reaction_smiles TEXT DEFAULT '',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # 创建索引以加速查询
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_ec_number ON motif_index(ec_number)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_uniprot_id ON motif_index(uniprot_id)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_category ON motif_index(category)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_nanozyme_type ON motif_index(nanozyme_type)
        """)
        # PR1-1 (H4 fix): expression index matches the case-insensitive query
        # below (`WHERE UPPER(nanozyme_type) = UPPER(?)`). Without this expression
        # index, SQLite linearly scans even though idx_nanozyme_type exists.
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_nanozyme_type_upper
            ON motif_index(UPPER(nanozyme_type))
        """)
        self._ensure_column(cursor, "motif_index", "source_pdb_id", "TEXT DEFAULT ''")
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_source_pdb_id ON motif_index(source_pdb_id)
        """)

        conn.commit()

    def _ensure_column(self, cursor, table_name: str, column_name: str, column_def: str):
        # PR1-1 (H2 fix): defensive validation — even though all current call sites
        # pass hardcoded literals, accept only [A-Za-z0-9_] in table/column names so
        # this method can't become an injection vector if reused with user input.
        import re
        if not re.match(r"^[A-Za-z0-9_]+$", table_name):
            raise ValueError(f"invalid table name: {table_name!r}")
        if not re.match(r"^[A-Za-z0-9_]+$", column_name):
            raise ValueError(f"invalid column name: {column_name!r}")
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = {row[1] for row in cursor.fetchall()}
        if column_name not in columns:
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {column_def}")

    def add_motif(self, motif: MotifIndex) -> bool:
        """
        添加或更新motif索引

        Args:
            motif: Motif索引对象

        Returns:
            是否成功
        """
        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            cursor.execute("""
                INSERT OR REPLACE INTO motif_index
                (motif_id, uniprot_id, source_pdb_id, ec_number, nanozyme_type, category,
                 anchor_atoms_count, file_path, confidence_score, chemistry_tag,
                 reaction_smiles, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """, (
                motif.motif_id,
                motif.uniprot_id,
                motif.source_pdb_id,
                motif.ec_number,
                motif.nanozyme_type,
                motif.category,
                motif.anchor_atoms_count,
                motif.file_path,
                motif.confidence_score,
                motif.chemistry_tag,
                motif.reaction_smiles
            ))
            conn.commit()
            return True
        except Exception as e:
            print(f"  ⚠️  添加motif失败 {motif.motif_id}: {e}")
            return False

    def add_motifs_batch(self, motifs: List[MotifIndex]) -> int:
        """
        PR1-1 (H5 fix): batch INSERT inside a single transaction.

        Single-row add_motif() commits on each call — that's correct for a
        random web write but ~10–100x too slow when ingesting thousands of
        motifs at indexing time. This API wraps the whole batch in one
        transaction; SQLite groups the disk syncs and the latency drops
        proportionally.

        Returns the count of successfully inserted rows.
        """
        if not motifs:
            return 0
        conn = self._get_connection()
        rows = [
            (m.motif_id, m.uniprot_id, m.source_pdb_id, m.ec_number,
             m.nanozyme_type, m.category, m.anchor_atoms_count, m.file_path,
             m.confidence_score, m.chemistry_tag, m.reaction_smiles)
            for m in motifs
        ]
        try:
            with conn:  # `with conn:` opens a transaction and commits on success
                conn.executemany("""
                    INSERT OR REPLACE INTO motif_index
                    (motif_id, uniprot_id, source_pdb_id, ec_number, nanozyme_type,
                     category, anchor_atoms_count, file_path, confidence_score,
                     chemistry_tag, reaction_smiles, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
                """, rows)
            return len(rows)
        except Exception as e:
            print(f"  ⚠️  批量添加motif失败: {e}")
            return 0

    def get_by_ec(self, ec_number: str) -> List[Dict]:
        """
        根据EC号查询motif

        Args:
            ec_number: EC号

        Returns:
            Motif列表（字典格式）
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT * FROM motif_index
            WHERE ec_number = ?
            ORDER BY confidence_score DESC, anchor_atoms_count DESC
        """, (ec_number,))

        rows = cursor.fetchall()
        return [dict(row) for row in rows]

    def get_by_category(self, ec_number: str, category: str) -> List[Dict]:
        """
        根据EC号和分类查询motif

        Args:
            ec_number: EC号
            category: 分类（metal_sites, catalytic_sites, binding_sites, other）

        Returns:
            Motif列表
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT * FROM motif_index
            WHERE ec_number = ? AND category = ?
            ORDER BY confidence_score DESC, anchor_atoms_count DESC
        """, (ec_number, category))

        rows = cursor.fetchall()
        return [dict(row) for row in rows]

    def get_by_id(self, motif_id: str) -> Optional[Dict]:
        """
        根据motif_id查询

        Args:
            motif_id: Motif ID

        Returns:
            Motif信息或None
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT * FROM motif_index
            WHERE motif_id = ?
        """, (motif_id,))

        row = cursor.fetchone()
        return dict(row) if row else None

    def get_all_ec_numbers(self) -> List[str]:
        """
        获取所有EC号列表

        Returns:
            EC号列表
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT DISTINCT ec_number FROM motif_index
            ORDER BY ec_number
        """)

        return [row[0] for row in cursor.fetchall()]

    def get_all_nanozyme_types(self) -> List[str]:
        """
        获取所有纳米酶类型列表

        Returns:
            纳米酶类型列表
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT DISTINCT nanozyme_type FROM motif_index
            WHERE nanozyme_type IS NOT NULL AND nanozyme_type != ''
            ORDER BY nanozyme_type
        """)

        return [row[0] for row in cursor.fetchall()]

    def get_by_nanozyme_type(self, nanozyme_type: str,
                              max_anchor_atoms: Optional[int] = None) -> List[Dict]:
        """
        根据纳米酶类型查询motif（使用SQL过滤，避免全表拉取导致卡顿）

        Args:
            nanozyme_type: 纳米酶类型（不区分大小写）
            max_anchor_atoms: PR1-1 (M30 fix) — 在 SQL 层过滤超大 motif。
                              之前调用方在 Python 层做 `if anchor_atoms_count > 50: continue`，
                              先把数千行拉到内存再过滤；改为 SQL WHERE 后直接减少 IO。
                              传 None 不过滤。

        Returns:
            Motif列表（字典格式）
        """
        if not nanozyme_type:
            return []
        conn = self._get_connection()
        cursor = conn.cursor()
        # PR1-1 (H4 fix): UPPER(nanozyme_type) hits idx_nanozyme_type_upper expression
        # index, so this is no longer a full table scan on case-insensitive match.
        sql = """
            SELECT * FROM motif_index
            WHERE UPPER(nanozyme_type) = UPPER(?)
        """
        params: list = [nanozyme_type]
        if max_anchor_atoms is not None:
            sql += " AND anchor_atoms_count <= ?"
            params.append(max_anchor_atoms)
        sql += " ORDER BY confidence_score DESC, anchor_atoms_count DESC"
        cursor.execute(sql, tuple(params))
        rows = cursor.fetchall()
        return [dict(row) for row in rows]

    def get_all(self) -> List[Dict]:
        """
        获取所有motif

        Returns:
            Motif列表（字典格式）
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT * FROM motif_index
            ORDER BY nanozyme_type, ec_number, confidence_score DESC
        """)

        rows = cursor.fetchall()
        return [dict(row) for row in rows]

    def count_by_ec(self, ec_number: str) -> Dict[str, int]:
        """
        统计指定EC号的motif数量（按分类）

        Args:
            ec_number: EC号

        Returns:
            各分类的数量统计
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("""
            SELECT category, COUNT(*) as count
            FROM motif_index
            WHERE ec_number = ?
            GROUP BY category
        """, (ec_number,))

        result = {}
        for row in cursor.fetchall():
            result[row[0]] = row[1]

        return result

    def clear(self):
        """清空数据库"""
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("DELETE FROM motif_index")
        conn.commit()

    def close(self):
        """关闭当前线程的数据库连接"""
        if hasattr(self._local, 'conn') and self._local.conn:
            self._local.conn.close()
            self._local.conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def _as_float_coords(value: Any) -> Optional[List[float]]:
    if not isinstance(value, (list, tuple)) or len(value) < 3:
        return None
    try:
        return [float(value[0]), float(value[1]), float(value[2])]
    except (TypeError, ValueError):
        return None


def _distance(a: List[float], b: List[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _residue_key(item: Dict) -> tuple:
    chain = item.get("chain_id", item.get("chain", ""))
    number = item.get("residue_number", item.get("residue_id", item.get("res_id")))
    try:
        number = int(number)
    except (TypeError, ValueError):
        number = None
    return (str(chain or ""), number)


def _residue_keys(item: Dict) -> set:
    chain, number = _residue_key(item)
    if number is None:
        return set()
    return {(chain, number), ("", number)}


def _active_site_keys(active_sites: Optional[List[Dict]]) -> set:
    keys = set()
    if not active_sites:
        return keys
    for site in active_sites:
        start = site.get("start", site.get("residue_number", site.get("position")))
        end = site.get("end", start)
        chain = str(site.get("chain_id", site.get("chain", "")) or "")
        try:
            start_i = int(start)
            end_i = int(end)
        except (TypeError, ValueError):
            continue
        for num in range(start_i, end_i + 1):
            keys.add((chain, num))
            keys.add(("", num))
    return keys


def _metal_site_is_linked_to_motif(site: Dict, anchor_atoms: List[Dict], active_sites: Optional[List[Dict]]) -> bool:
    """Return True only when the metal has residue/space/annotation evidence tied to the motif."""
    if not site.get("coordinating_residues"):
        return False

    anchor_keys = set()
    for atom in anchor_atoms:
        anchor_keys.update(_residue_keys(atom))
    coord_keys = set()
    for residue in site.get("coordinating_residues", []):
        coord_keys.update(_residue_keys(residue))
    if anchor_keys & coord_keys:
        return True

    active_keys = _active_site_keys(active_sites)
    if active_keys and coord_keys & active_keys:
        return True

    metal_coords = _as_float_coords(site.get("metal_coords") or site.get("coordinates"))
    if metal_coords:
        for atom in anchor_atoms:
            atom_coords = _as_float_coords(atom.get("coordinates"))
            if atom_coords and _distance(metal_coords, atom_coords) <= 6.0:
                return True

    for residue in site.get("coordinating_residues", []):
        residue_coords = _as_float_coords(residue.get("coordinates"))
        if not residue_coords:
            continue
        for atom in anchor_atoms:
            atom_coords = _as_float_coords(atom.get("coordinates"))
            if atom_coords and _distance(residue_coords, atom_coords) <= 2.0:
                return True

    return False


def classify_motif(motif_data: Dict, active_sites: Optional[List[Dict]] = None) -> str:
    """
    Classify motif content without promoting every PDB metal into metal_sites.

    A motif is metal-related only when at least one metal site is tied to the motif
    anchors by residue identity, active-site annotation, or close spatial contact.
    """
    anchor_atoms = motif_data.get('anchor_atoms', []) or []
    nanozyme_type = str(motif_data.get('nanozyme_type', '') or '').upper()
    for site in motif_data.get('metal_sites', []) or []:
        if _metal_site_is_linked_to_motif(site, anchor_atoms, active_sites):
            return 'metal_sites'

    has_catalytic = False
    has_binding = False

    # 常见的金属配位残基
    metal_coordinating_residues = ['HIS', 'CYS', 'ASP', 'GLU', 'TYR', 'SER', 'THR']

    has_metal = False
    for atom in anchor_atoms:
        role = atom.get('role', '').lower()
        element = atom.get('element', '').upper()
        residue_name = atom.get('residue_name', '').upper()

        # 检查金属元素
        if element in KNOWN_METAL_ELEMENTS and active_sites:
            has_metal = True

        # 检查角色
        if 'metal' in role or 'ligand' in role:
            has_metal = True
        if 'catalytic' in role or 'active' in role or 'nucleophile' in role or 'base' in role or 'acid' in role:
            has_catalytic = True
        if 'binding' in role:
            has_binding = True

        # 如果没有role信息，根据残基类型推断
        if not role:
            # 金属配位残基通常用于金属位点
            if residue_name in metal_coordinating_residues:
                # 检查是否是常见的催化残基
                catalytic_residues = ['HIS', 'CYS', 'SER', 'ASP', 'GLU', 'LYS', 'ARG', 'TYR']
                if residue_name in catalytic_residues:
                    has_catalytic = True
                else:
                    has_binding = True

    # 根据nanozyme类型推断（如果还没有明确分类）
    if not has_catalytic and not has_binding:
        # 某些nanozyme类型通常涉及金属
        metal_related_types = ['SOD', 'SUPEROXIDE DISMUTASE', 'CATALASE', 'CAT', 'LACCASE', 'LAC']
        if any(mt in nanozyme_type for mt in metal_related_types):
            has_catalytic = True  # 这些类型通常是催化位点

    # 优先级：金属位点 > 催化位点 > 结合位点 > 其他
    if has_metal:
        return 'metal_sites'
    elif has_catalytic:
        return 'catalytic_sites'
    elif has_binding:
        return 'binding_sites'
    else:
        # 默认：如果有锚点原子，至少是结合位点
        if len(anchor_atoms) > 0:
            return 'binding_sites'
        return 'other'


def build_motif_index(motif_library_dir: Path, db_path: Path, clear_existing: bool = False):
    """
    构建motif索引数据库

    Args:
        motif_library_dir: Motif库目录
        db_path: 数据库文件路径
        clear_existing: 是否清空现有数据
    """
    print(f"构建Motif索引数据库...")
    print(f"  Motif库目录: {motif_library_dir}")
    print(f"  数据库路径: {db_path}")

    with MotifDatabase(db_path) as db:
        if clear_existing:
            print("  清空现有数据...")
            db.clear()

        # 遍历目录结构：支持两种组织方式
        # 方式1: motif_library/{EC_number}/{category}/*.json (新格式)
        # 方式2: motif_library/{nanozyme_type}/*.json (旧格式)
        motif_count = 0

        for top_dir in motif_library_dir.iterdir():
            if not top_dir.is_dir():
                continue

            print(f"  扫描目录: {top_dir.name}")

            # 检查是否是EC号目录（包含分类子目录）
            category_dirs = ['metal_sites', 'catalytic_sites', 'binding_sites', 'other']
            has_category_dirs = any((top_dir / cat).is_dir() for cat in category_dirs)

            if has_category_dirs:
                # 方式1: EC号目录结构
                ec_number = top_dir.name.replace("_", ".")

                for category in category_dirs:
                    category_dir = top_dir / category
                    if not category_dir.is_dir():
                        continue

                    for motif_file in category_dir.glob("*.json"):
                        try:
                            with open(motif_file, 'r') as f:
                                motif_data = json.load(f)

                            # 提取motif信息
                            motif_id = motif_data.get('motif_id', motif_file.stem)
                            uniprot_id = motif_data.get('source_uniprot_id', '')
                            ec_number_from_data = motif_data.get('source_ec_number', ec_number)
                            nanozyme_type = motif_data.get('nanozyme_type', '')
                            anchor_atoms = motif_data.get('anchor_atoms', [])

                            # 使用目录中的分类（更准确）
                            # 也可以从数据中重新分类验证
                            category_from_data = classify_motif(motif_data)

                            # 创建索引对象
                            motif_index = MotifIndex(
                                motif_id=motif_id,
                                uniprot_id=uniprot_id,
                                source_pdb_id=motif_data.get('source_pdb_id', ''),
                                ec_number=ec_number_from_data,
                                nanozyme_type=nanozyme_type,
                                category=category,  # 使用目录分类
                                anchor_atoms_count=len(anchor_atoms),
                                file_path=str(motif_file),
                                confidence_score=motif_data.get('confidence_score', 0.0),
                                chemistry_tag=motif_data.get('chemistry_tag', ''),
                                reaction_smiles=motif_data.get('reaction_smiles', '')
                            )

                            # 添加到数据库
                            if db.add_motif(motif_index):
                                motif_count += 1

                        except Exception as e:
                            print(f"    ⚠️  处理 {motif_file.name} 失败: {e}")
                            continue
            else:
                # 方式2: 旧格式（nanozyme类型目录）
                for motif_file in top_dir.glob("*.json"):
                    try:
                        with open(motif_file, 'r') as f:
                            motif_data = json.load(f)

                        # 提取motif信息
                        motif_id = motif_data.get('motif_id', motif_file.stem)
                        uniprot_id = motif_data.get('source_uniprot_id', '')
                        ec_number = motif_data.get('source_ec_number', '')
                        nanozyme_type = motif_data.get('nanozyme_type', top_dir.name)
                        anchor_atoms = motif_data.get('anchor_atoms', [])

                        # 分类（尝试从文件路径获取active_sites信息）
                        # 注意：这里无法直接获取active_sites，使用motif数据推断
                        category = classify_motif(motif_data)

                        # 创建索引对象
                        motif_index = MotifIndex(
                            motif_id=motif_id,
                            uniprot_id=uniprot_id,
                            source_pdb_id=motif_data.get('source_pdb_id', ''),
                            ec_number=ec_number,
                            nanozyme_type=nanozyme_type,
                            category=category,
                            anchor_atoms_count=len(anchor_atoms),
                            file_path=str(motif_file),
                            confidence_score=motif_data.get('confidence_score', 0.0),
                            chemistry_tag=motif_data.get('chemistry_tag', ''),
                            reaction_smiles=motif_data.get('reaction_smiles', '')
                        )

                        # 添加到数据库
                        if db.add_motif(motif_index):
                            motif_count += 1

                    except Exception as e:
                        print(f"    ⚠️  处理 {motif_file.name} 失败: {e}")
                        continue

        print(f"\n✓ 索引构建完成，共索引 {motif_count} 个motif")

        # 显示统计信息
        all_ec = db.get_all_ec_numbers()
        print(f"  覆盖 {len(all_ec)} 个EC号")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='构建Motif索引数据库')
    parser.add_argument(
        '--motif-library',
        type=str,
        default='motif_library',
        help='Motif库目录路径（默认: motif_library）'
    )
    parser.add_argument(
        '--db-path',
        type=str,
        default='enzyme_viewer/motif_index.db',
        help='数据库文件路径（默认: enzyme_viewer/motif_index.db）'
    )
    parser.add_argument(
        '--clear',
        action='store_true',
        help='清空现有数据后重建'
    )

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    motif_library_dir = project_root / args.motif_library
    db_path = project_root / args.db_path

    if not motif_library_dir.exists():
        print(f"错误: Motif库目录不存在: {motif_library_dir}")
        exit(1)

    build_motif_index(motif_library_dir, db_path, clear_existing=args.clear)
