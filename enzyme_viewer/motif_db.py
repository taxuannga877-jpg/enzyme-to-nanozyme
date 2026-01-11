"""
Motif数据库索引模块
使用SQLite存储motif的索引信息，加速查询
"""

import json
import sqlite3
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class MotifIndex:
    """Motif索引数据类"""
    motif_id: str
    uniprot_id: str
    ec_number: str
    nanozyme_type: str
    category: str  # metal_sites, catalytic_sites, binding_sites, other
    anchor_atoms_count: int
    file_path: str
    confidence_score: float = 0.0
    chemistry_tag: str = ""
    reaction_smiles: str = ""


class MotifDatabase:
    """Motif数据库管理类"""
    
    def __init__(self, db_path: Path):
        """
        初始化数据库
        
        Args:
            db_path: SQLite数据库文件路径
        """
        self.db_path = Path(db_path)
        self.conn = None
        self._init_database()
    
    def _init_database(self):
        """初始化数据库表结构"""
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row  # 返回字典格式
        
        cursor = self.conn.cursor()
        
        # 创建motif索引表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS motif_index (
                motif_id TEXT PRIMARY KEY,
                uniprot_id TEXT NOT NULL,
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
        
        self.conn.commit()
    
    def add_motif(self, motif: MotifIndex) -> bool:
        """
        添加或更新motif索引
        
        Args:
            motif: Motif索引对象
            
        Returns:
            是否成功
        """
        try:
            cursor = self.conn.cursor()
            cursor.execute("""
                INSERT OR REPLACE INTO motif_index 
                (motif_id, uniprot_id, ec_number, nanozyme_type, category, 
                 anchor_atoms_count, file_path, confidence_score, chemistry_tag, 
                 reaction_smiles, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """, (
                motif.motif_id,
                motif.uniprot_id,
                motif.ec_number,
                motif.nanozyme_type,
                motif.category,
                motif.anchor_atoms_count,
                motif.file_path,
                motif.confidence_score,
                motif.chemistry_tag,
                motif.reaction_smiles
            ))
            self.conn.commit()
            return True
        except Exception as e:
            print(f"  ⚠️  添加motif失败 {motif.motif_id}: {e}")
            return False
    
    def get_by_ec(self, ec_number: str) -> List[Dict]:
        """
        根据EC号查询motif
        
        Args:
            ec_number: EC号
            
        Returns:
            Motif列表（字典格式）
        """
        cursor = self.conn.cursor()
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
        cursor = self.conn.cursor()
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
        cursor = self.conn.cursor()
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
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT DISTINCT ec_number FROM motif_index
            ORDER BY ec_number
        """)
        
        return [row[0] for row in cursor.fetchall()]
    
    def count_by_ec(self, ec_number: str) -> Dict[str, int]:
        """
        统计指定EC号的motif数量（按分类）
        
        Args:
            ec_number: EC号
            
        Returns:
            各分类的数量统计
        """
        cursor = self.conn.cursor()
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
        cursor = self.conn.cursor()
        cursor.execute("DELETE FROM motif_index")
        self.conn.commit()
    
    def close(self):
        """关闭数据库连接"""
        if self.conn:
            self.conn.close()
            self.conn = None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def classify_motif(motif_data: Dict, active_sites: Optional[List[Dict]] = None) -> str:
    """
    根据motif的特征分类
    
    Args:
        motif_data: Motif数据字典
        active_sites: 活性位点信息列表（可选，用于更准确的分类）
        
    Returns:
        分类名称（metal_sites, catalytic_sites, binding_sites, other）
    """
    anchor_atoms = motif_data.get('anchor_atoms', [])
    nanozyme_type = motif_data.get('nanozyme_type', '').upper()
    
    # 首先检查active_sites中的type信息（最准确）
    if active_sites:
        has_metal_site = False
        has_active_site = False
        has_binding_site = False
        
        # 获取motif涉及的残基编号
        motif_residues = {atom.get('residue_number') for atom in anchor_atoms}
        
        for site in active_sites:
            site_type = site.get('type', '').lower()
            start = site.get('start', 0)
            end = site.get('end', start)
            site_residues = set(range(start, end + 1))
            
            # 检查是否有重叠
            if motif_residues & site_residues:
                if 'metal' in site_type or 'metal binding' in site_type:
                    has_metal_site = True
                elif 'active' in site_type or 'catalytic' in site_type:
                    has_active_site = True
                elif 'binding' in site_type:
                    has_binding_site = True
        
        # 优先级：金属位点 > 活性位点 > 结合位点
        if has_metal_site:
            return 'metal_sites'
        elif has_active_site:
            return 'catalytic_sites'
        elif has_binding_site:
            return 'binding_sites'
    
    # 如果没有active_sites信息，使用原子特征推断
    has_metal = False
    has_catalytic = False
    has_binding = False
    
    # 常见的金属配位残基
    metal_coordinating_residues = ['HIS', 'CYS', 'ASP', 'GLU', 'TYR', 'SER', 'THR']
    
    for atom in anchor_atoms:
        role = atom.get('role', '').lower()
        element = atom.get('element', '').upper()
        residue_name = atom.get('residue_name', '').upper()
        
        # 检查金属元素
        metal_elements = ['FE', 'CU', 'ZN', 'MG', 'MN', 'CO', 'NI', 'CA', 'MO', 'W']
        if element in metal_elements:
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
    if not has_metal and not has_catalytic and not has_binding:
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

