#!/usr/bin/env python
"""
从JSON缓存批量提取Motif并覆盖motif_library
===========================================

流程：
1. 读取cache/json/*_sites.json文件
2. 确保PDB文件已下载
3. 提取motif（按nanozyme类型组织）
4. 覆盖motif_library目录
5. 建立数据库索引

使用方法:
    python scripts/extract_motifs_from_cache.py
    python scripts/extract_motifs_from_cache.py --ec 1.11.1.7
    python scripts/extract_motifs_from_cache.py --clear
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from nanozyme_mining.extraction.extractor import MotifExtractor
from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE, NanozymeType
from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher
from enzyme_viewer.motif_db import build_motif_index, classify_motif


# EC号到Nanozyme类型名称的映射（用于目录命名）
EC_TO_NANOZYME_NAME = {
    "1.11.1.7": "Peroxidase",
    "1.11.1.11": "Peroxidase",
    "1.11.1.21": "Peroxidase",
    "1.11.1.6": "Catalase",
    "1.15.1.1": "Superoxide Dismutase",
    "1.11.1.9": "Glutathione Peroxidase",
    "1.11.1.12": "Glutathione Peroxidase",
    "1.1.3.4": "Glucose Oxidase",
    "1.10.3.2": "Laccase",
    "1.4.3.4": "Oxidase",
    "1.3.3.4": "Oxidase",
    "3.1.3.1": "Phosphatase",
    "3.1.21.1": "DNase",
}


def get_nanozyme_name(ec_number: str) -> str:
    """根据EC号获取nanozyme类型名称"""
    return EC_TO_NANOZYME_NAME.get(ec_number, "Other")


# 不再使用黑名单 - 通过进程隔离和错误处理来防止崩溃


def ensure_pdb_downloaded(
    uniprot_id: str,
    alphafold_id: str,
    pdb_cache_dir: Path,
    fetcher: UniProtFetcher,
    ec_number: Optional[str] = None,
    pdb_library_dir: Optional[Path] = None,
    pdb_id: Optional[str] = None
) -> Optional[Path]:
    """
    确保PDB文件已下载 - 优先从pdb_library查找，然后从cache查找，最后下载
    优先查找实验PDB，如果没有再查找AlphaFold PDB
    
    Args:
        uniprot_id: UniProt ID
        alphafold_id: AlphaFold ID
        pdb_cache_dir: PDB缓存目录（旧目录，兼容性）
        fetcher: UniProtFetcher实例
        ec_number: EC号（用于在pdb_library中查找）
        pdb_library_dir: PDB库目录（新目录，按EC号组织）
        pdb_id: 实验PDB ID（可选，优先使用）
        
    Returns:
        PDB文件路径，如果下载失败返回None
    """
    # 尝试多种可能的文件名格式
    AFDB_VERSION = "6"
    
    # ========== 优先从pdb_library查找（按EC号组织）==========
    if pdb_library_dir and ec_number:
        ec_dir_name = ec_number.replace('.', '_')
        pdb_library_ec_dir = pdb_library_dir / ec_dir_name
        
        if pdb_library_ec_dir.exists():
            # 优先级1: 查找实验PDB (格式: {PDB_ID}.pdb)
            if pdb_id:
                pdb_id = pdb_id.upper().strip()
                exp_pdb_filename = f"{pdb_id}.pdb"
                exp_pdb_path = pdb_library_ec_dir / exp_pdb_filename
                if exp_pdb_path.exists():
                    return exp_pdb_path
            
            # 优先级2: AlphaFold PDB - 标准格式
            pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
            pdb_path = pdb_library_ec_dir / pdb_filename
            
            if pdb_path.exists():
                return pdb_path
            
            # 优先级3: AlphaFold PDB - 查找任何包含该ID的PDB文件
            matching_pdb = list(pdb_library_ec_dir.glob(f"AF-{alphafold_id}-*.pdb"))
            if matching_pdb:
                return matching_pdb[0]
            
            # 优先级4: 如果alphafold_id和uniprot_id不同，尝试uniprot_id
            if uniprot_id and uniprot_id != alphafold_id:
                pdb_filename = f"AF-{uniprot_id}-F1-model_v{AFDB_VERSION}.pdb"
                pdb_path = pdb_library_ec_dir / pdb_filename
                if pdb_path.exists():
                    return pdb_path
                
                matching_pdb = list(pdb_library_ec_dir.glob(f"AF-{uniprot_id}-*.pdb"))
                if matching_pdb:
                    return matching_pdb[0]
    
    # ========== 回退到旧缓存目录（兼容性）==========
    # 方法1: 标准格式
    pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
    pdb_path = pdb_cache_dir / pdb_filename
    
    if pdb_path.exists():
        return pdb_path
    
    # 方法2: 查找任何包含该ID的PDB文件
    matching_pdb = list(pdb_cache_dir.glob(f"AF-{alphafold_id}-*.pdb"))
    if matching_pdb:
        return matching_pdb[0]
    
    # 方法3: 如果alphafold_id和uniprot_id不同，尝试uniprot_id
    if uniprot_id and uniprot_id != alphafold_id:
        pdb_filename = f"AF-{uniprot_id}-F1-model_v{AFDB_VERSION}.pdb"
        pdb_path = pdb_cache_dir / pdb_filename
        if pdb_path.exists():
            return pdb_path
        
        matching_pdb = list(pdb_cache_dir.glob(f"AF-{uniprot_id}-*.pdb"))
        if matching_pdb:
            return matching_pdb[0]
    
    # 方法4: 尝试下载
    print(f"    ⚠️  PDB文件不存在，尝试下载 {alphafold_id}...")
    try:
        download_id = alphafold_id or uniprot_id
        pdb_path = fetcher.download_pdb(download_id)
        if pdb_path and pdb_path.exists():
            print(f"    ✓ 下载成功: {pdb_path.name}")
            return pdb_path
    except Exception as e:
        print(f"    ✗ 下载失败: {e}")
    
    return None


def extract_motif_from_entry(
    entry: Dict,
    ec_number: str,
    extractor: MotifExtractor,
    fetcher: UniProtFetcher,
    pdb_cache_dir: Path,
    pdb_library_dir: Optional[Path] = None
) -> Optional[Dict]:
    """
    从单个条目提取motif
    
    Args:
        entry: JSON缓存中的条目
        ec_number: EC号
        extractor: MotifExtractor实例
        fetcher: UniProtFetcher实例
        pdb_cache_dir: PDB缓存目录（旧目录，兼容性）
        pdb_library_dir: PDB库目录（新目录，按EC号组织）
        
    Returns:
        提取结果字典，失败返回None
    """
    uniprot_id = entry.get('uniprot_id', '')
    alphafold_id = entry.get('alphafold_id', uniprot_id)
    pdb_id = entry.get('pdb_id', '')
    
    active_sites = entry.get('active_sites', [])
    
    # 确保PDB文件存在 - 优先从pdb_library查找（优先使用实验PDB）
    pdb_path = ensure_pdb_downloaded(
        uniprot_id, alphafold_id, pdb_cache_dir, fetcher,
        ec_number=ec_number, pdb_library_dir=pdb_library_dir,
        pdb_id=pdb_id
    )
    if not pdb_path:
        return None
    
    # 提取活性位点残基索引
    active_site_indices = []
    for site in active_sites:
        start = site.get('start', 0)
        end = site.get('end', start)
        active_site_indices.extend(range(start, end + 1))
    
    # 获取nanozyme类型
    nanozyme_type_enum = EC_TO_NANOZYME_TYPE.get(ec_number)
    if not nanozyme_type_enum:
        print(f"    ⚠️  未找到EC号 {ec_number} 对应的nanozyme类型")
        return None
    
    nanozyme_name = get_nanozyme_name(ec_number)
    
    # 提取motif（仅从PDB文件提取，不使用模型预测）
    try:
        motif = extractor.extract_motif(
            pdb_path=str(pdb_path),
            uniprot_id=uniprot_id,
            ec_number=ec_number,
            nanozyme_type=nanozyme_name,
            active_site_indices=active_site_indices if active_site_indices else None
        )
        
        if motif:
            return {
                'motif': motif,
                'nanozyme_name': nanozyme_name,
                'pdb_path': str(pdb_path)
            }
        else:
            return None
            
    except Exception as e:
        print(f"    ✗ 提取失败: {e}")
        return None


def classify_motif_from_entry(motif_data: Dict, active_sites: List[Dict]) -> str:
    """
    根据motif和active_sites信息分类
    
    Args:
        motif_data: Motif数据字典
        active_sites: 活性位点信息列表
        
    Returns:
        分类名称（metal_sites, catalytic_sites, binding_sites, other）
    """
    return classify_motif(motif_data, active_sites=active_sites)


def process_ec_number(
    ec_number: str,
    json_cache_dir: Path,
    pdb_cache_dir: Path,
    motif_library_dir: Path,
    fetcher: UniProtFetcher,
    clear_existing: bool = False,
    pdb_library_dir: Optional[Path] = None
):
    """
    处理单个EC号的所有motif提取，按EC号和类型分类组织
    
    目录结构：
    motif_library/
      {EC_number}/
        metal_sites/
        catalytic_sites/
        binding_sites/
        other/
    
    Args:
        ec_number: EC号
        json_cache_dir: JSON缓存目录
        pdb_cache_dir: PDB缓存目录
        motif_library_dir: Motif库目录
        fetcher: UniProtFetcher实例
        clear_existing: 是否清空现有motif
    """
    print(f"\n{'='*80}")
    print(f"处理 EC {ec_number}")
    print(f"{'='*80}")
    
    # 读取JSON缓存（优先从pdb_library读取）
    if pdb_library_dir:
        ec_dir_name = ec_number.replace(".", "_")
        json_file = pdb_library_dir / ec_dir_name / f"{ec_number}_sites.json"
    else:
        json_file = json_cache_dir / f"{ec_number}_sites.json"
    
    if not json_file.exists():
        print(f"  ⚠️  JSON缓存文件不存在: {json_file}")
        return
    
    with open(json_file, 'r', encoding='utf-8') as f:
        enzyme_data = json.load(f)
    
    if not enzyme_data:
        print(f"  ⚠️  EC {ec_number} 的数据为空")
        return
    
    print(f"  ✓ 找到 {len(enzyme_data)} 个酶结构")
    
    # 按EC号组织目录结构
    ec_output_dir = motif_library_dir / ec_number.replace(".", "_")
    
    # 创建分类子目录
    category_dirs = {
        'metal_sites': ec_output_dir / 'metal_sites',
        'catalytic_sites': ec_output_dir / 'catalytic_sites',
        'binding_sites': ec_output_dir / 'binding_sites',
        'other': ec_output_dir / 'other'
    }
    
    # 清空现有motif（如果指定）
    if clear_existing and ec_output_dir.exists():
        print(f"  🗑️  清空现有motif目录: {ec_output_dir}")
        import shutil
        shutil.rmtree(ec_output_dir)
    
    # 创建所有分类目录
    for cat_dir in category_dirs.values():
        cat_dir.mkdir(parents=True, exist_ok=True)
    
    # 初始化提取器（使用临时目录，实际保存到分类目录）
    extractor = MotifExtractor(output_dir=str(ec_output_dir))
    
    # 处理每个酶
    success_count = 0
    fail_count = 0
    skip_count = 0
    category_counts = {cat: 0 for cat in category_dirs.keys()}
    
    for idx, entry in enumerate(enzyme_data, 1):
        uniprot_id = entry.get('uniprot_id', 'Unknown')
        
        # 检查是否已处理过（检查所有分类目录）
        nanozyme_name = get_nanozyme_name(ec_number)
        expected_motif_id = f"{uniprot_id}_{ec_number}_{nanozyme_name}"
        already_processed = False
        
        for cat_dir in category_dirs.values():
            if (cat_dir / f"{expected_motif_id}.json").exists():
                already_processed = True
                break
        
        if already_processed:
            skip_count += 1
            print(f"  [{idx}/{len(enzyme_data)}] 跳过 {uniprot_id} (已处理)")
            continue
        
        print(f"  [{idx}/{len(enzyme_data)}] 处理 {uniprot_id}...")
        
        result = extract_motif_from_entry(
            entry, ec_number, extractor, fetcher, pdb_cache_dir,
            pdb_library_dir=pdb_library_dir
        )
        
        if result:
            motif = result['motif']
            
            # 获取active_sites信息用于分类
            active_sites = entry.get('active_sites', [])
            
            # 分类motif
            motif_dict = motif.to_dict()
            category = classify_motif_from_entry(motif_dict, active_sites)
            
            # 保存到对应的分类目录
            category_dir = category_dirs[category]
            output_file = category_dir / f"{motif.motif_id}.json"
            motif.to_json(str(output_file))
            
            success_count += 1
            category_counts[category] += 1
            print(f"    ✓ 提取成功: {motif.motif_id} ({len(motif.anchor_atoms)} 个锚点原子) → {category}")
        else:
            fail_count += 1
            print(f"    ✗ 提取失败")
    
    print(f"\n  [完成] EC {ec_number}:")
    print(f"    ✓ 成功: {success_count}")
    print(f"    ⏭️  跳过: {skip_count}")
    print(f"    ✗ 失败: {fail_count}")
    print(f"    📊 分类统计:")
    for cat, count in category_counts.items():
        if count > 0:
            print(f"      - {cat}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="从JSON缓存批量提取Motif并覆盖motif_library"
    )
    parser.add_argument(
        "--ec",
        type=str,
        help="要处理的EC号（如 1.11.1.7），可以指定多个，用逗号分隔。不指定则处理所有EC号"
    )
    parser.add_argument(
        "--clear",
        action="store_true",
        help="清空现有motif_library后重建"
    )
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="./cache",
        help="缓存目录（默认: ./cache）"
    )
    parser.add_argument(
        "--motif-dir",
        type=str,
        default="./motif_library",
        help="Motif库目录（默认: ./motif_library）"
    )
    parser.add_argument(
        "--skip-index",
        action="store_true",
        help="跳过数据库索引构建"
    )
    
    args = parser.parse_args()
    
    # 设置路径
    project_root = Path(__file__).parent.parent
    pdb_library_dir = project_root / "pdb_library"  # PDB库目录（按EC号组织，包含JSON和PDB）
    motif_library_dir = project_root / args.motif_dir
    
    # 保留旧路径用于兼容性（但不再使用）
    cache_dir = project_root / args.cache_dir if args.cache_dir else None
    json_cache_dir = cache_dir / "json" if cache_dir else None
    pdb_cache_dir = cache_dir / "pdb" if cache_dir else None
    
    # 清空motif_library（如果指定）
    if args.clear and motif_library_dir.exists():
        print(f"🗑️  清空现有motif_library: {motif_library_dir}")
        import shutil
        shutil.rmtree(motif_library_dir)
    
    motif_library_dir.mkdir(parents=True, exist_ok=True)
    
    # 初始化UniProtFetcher（使用pdb_library）
    fetcher = UniProtFetcher(
        cache_dir=str(cache_dir) if cache_dir else "./cache",
        pdb_library_dir=str(pdb_library_dir)
    )
    
    # 确定要处理的EC号列表
    if args.ec:
        ec_numbers = [ec.strip() for ec in args.ec.split(",")]
    else:
        # 扫描所有_sites.json文件（从pdb_library）
        if pdb_library_dir.exists():
            json_files = list(pdb_library_dir.glob("*/*_sites.json"))
            ec_numbers = [f.stem.replace("_sites", "") for f in json_files]
            ec_numbers.sort()
        elif json_cache_dir and json_cache_dir.exists():
            # 向后兼容：如果pdb_library不存在，尝试旧路径
            json_files = list(json_cache_dir.glob("*_sites.json"))
            ec_numbers = [f.stem.replace("_sites", "") for f in json_files]
            ec_numbers.sort()
        else:
            ec_numbers = []
    
    if not ec_numbers:
        print("⚠️  未找到任何EC号数据")
        return
    
    print("\n" + "="*80)
    print("Motif提取流程")
    print("="*80)
    print(f"EC号列表: {ec_numbers}")
    print(f"PDB库: {pdb_library_dir} (包含JSON和PDB文件)")
    print(f"Motif库: {motif_library_dir}")
    print("="*80)
    
    # 处理每个EC号
    total_success = 0
    total_fail = 0
    
    for i, ec_number in enumerate(ec_numbers, 1):
        print(f"\n[{i}/{len(ec_numbers)}] 处理 EC {ec_number}...")
        try:
            process_ec_number(
                ec_number=ec_number,
                json_cache_dir=json_cache_dir,  # 保留用于兼容性，但不再使用
                pdb_cache_dir=pdb_cache_dir,  # 保留用于兼容性，但不再使用
                motif_library_dir=motif_library_dir,
                fetcher=fetcher,
                clear_existing=args.clear and i == 1,  # 只在第一次清空
                pdb_library_dir=pdb_library_dir
            )
        except KeyboardInterrupt:
            print("\n\n用户中断，退出...")
            break
        except Exception as e:
            print(f"\n  ✗ 处理 EC {ec_number} 时出错: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*80)
    print("Motif提取完成！")
    print("="*80)
    print(f"Motif库位置: {motif_library_dir}")
    
    # 构建数据库索引
    if not args.skip_index:
        print("\n" + "="*80)
        print("构建Motif数据库索引...")
        print("="*80)
        
        db_path = project_root / "enzyme_viewer" / "motif_index.db"
        build_motif_index(
            motif_library_dir=motif_library_dir,
            db_path=db_path,
            clear_existing=True
        )
        
        print("\n✓ 数据库索引构建完成！")


if __name__ == "__main__":
    main()

