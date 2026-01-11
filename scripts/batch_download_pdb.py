#!/usr/bin/env python3
"""
批量下载PDB文件脚本 - 纳米酶EC号专用
支持两种模式：
1. 从JSON缓存扫描（已有数据）
2. 直接从UniProt查询所有纳米酶EC号并下载（大量下载）
"""

import json
import sys
from pathlib import Path
from typing import Set, List, Dict
import time

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher, FETCH_ALL_RESULTS
from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE
from nanozyme_mining.utils.ec_mappings import EC_PATTERNS


def scan_json_cache(json_cache_dir: Path) -> Dict[str, List[Dict]]:
    """
    扫描JSON缓存目录，提取所有需要下载的PDB信息
    
    Returns:
        Dict mapping EC number to list of entries with alphafold_id
    """
    ec_entries = {}
    
    print(f"扫描JSON缓存目录: {json_cache_dir}")
    
    for json_file in json_cache_dir.glob("*_sites.json"):
        ec_number = json_file.stem.replace("_sites", "")
        
        try:
            with open(json_file, 'r') as f:
                entries = json.load(f)
            
            if ec_number not in ec_entries:
                ec_entries[ec_number] = []
            
            for entry in entries:
                alphafold_id = entry.get('alphafold_id', '').strip()
                uniprot_id = entry.get('uniprot_id', '').strip()
                
                if alphafold_id:
                    ec_entries[ec_number].append({
                        'alphafold_id': alphafold_id,
                        'uniprot_id': uniprot_id
                    })
        
        except Exception as e:
            print(f"  ⚠️  读取 {json_file.name} 失败: {e}")
            continue
    
    return ec_entries


def check_existing_pdb(pdb_cache_dir: Path, alphafold_id: str) -> bool:
    """检查PDB文件是否已存在"""
    # AlphaFold DB v6格式
    pdb_file = pdb_cache_dir / f"AF-{alphafold_id}-F1-model_v6.pdb"
    return pdb_file.exists() and pdb_file.stat().st_size > 0


def get_all_nanozyme_ec_numbers() -> List[str]:
    """获取所有纳米酶EC号列表"""
    all_ecs = []
    for nanozyme_type, ec_list in EC_PATTERNS.items():
        all_ecs.extend(ec_list)
    return sorted(set(all_ecs))


def query_ec_from_uniprot(
    fetcher: UniProtFetcher,
    ec_number: str,
    max_per_ec: int = FETCH_ALL_RESULTS
) -> List[Dict]:
    """
    从UniProt查询EC号的所有条目（包含AlphaFold ID）
    
    Args:
        fetcher: UniProtFetcher实例
        ec_number: EC号
        max_per_ec: 每个EC号最大查询数量（-1表示全部）
    
    Returns:
        包含alphafold_id和uniprot_id的条目列表
    """
    print(f"  查询 EC {ec_number}...", end=' ')
    
    # 使用query_with_active_sites获取详细信息（包含alphafold_id）
    entries = fetcher.query_with_active_sites(ec_number, max_results=max_per_ec)
    
    # 提取alphafold_id
    result = []
    for entry in entries:
        alphafold_id = entry.get('alphafold_id', '').strip()
        uniprot_id = entry.get('uniprot_id', '').strip()
        
        if alphafold_id:  # 只保留有AlphaFold结构的
            result.append({
                'alphafold_id': alphafold_id,
                'uniprot_id': uniprot_id
            })
    
    print(f"找到 {len(result)} 个有AlphaFold结构的条目")
    return result


def batch_download_pdb(
    json_cache_dir: Path = None,
    pdb_cache_dir: Path = None,
    cache_dir: str = "./cache",
    max_entries: int = None,
    max_per_ec: int = FETCH_ALL_RESULTS,
    delay: float = 0.5,
    use_cache: bool = True,
    ec_numbers: List[str] = None
):
    """
    批量下载PDB文件
    
    Args:
        json_cache_dir: JSON缓存目录（如果使用缓存模式）
        pdb_cache_dir: PDB缓存目录
        cache_dir: 缓存根目录
        max_entries: 全局最大下载条目数（None表示全部）
        max_per_ec: 每个EC号最大查询数量（-1表示全部，默认全部）
        delay: 每次下载之间的延迟（秒）
        use_cache: 是否使用JSON缓存（False则直接从UniProt查询）
        ec_numbers: 要处理的EC号列表（None表示所有纳米酶EC号）
    """
    # 初始化UniProtFetcher
    fetcher = UniProtFetcher(cache_dir=cache_dir)
    
    if pdb_cache_dir is None:
        pdb_cache_dir = Path(cache_dir) / "pdb"
    pdb_cache_dir = Path(pdb_cache_dir)
    pdb_cache_dir.mkdir(parents=True, exist_ok=True)
    
    # 确定要处理的EC号列表
    if ec_numbers is None:
        ec_numbers = get_all_nanozyme_ec_numbers()
    
    print("=" * 80)
    print("纳米酶PDB批量下载工具")
    print("=" * 80)
    print(f"处理 {len(ec_numbers)} 个纳米酶EC号")
    print(f"每个EC号最大查询: {'全部' if max_per_ec == FETCH_ALL_RESULTS else max_per_ec}")
    print(f"使用缓存模式: {use_cache}")
    print("=" * 80)
    
    # 步骤1: 获取所有需要下载的条目
    print("\n[步骤1] 获取所有需要下载的条目...")
    all_entries = []
    
    if use_cache and json_cache_dir and json_cache_dir.exists():
        # 模式1: 从JSON缓存扫描
        print("  模式: 从JSON缓存扫描")
        ec_entries = scan_json_cache(json_cache_dir)
        for ec_number, entries in ec_entries.items():
            for entry in entries:
                all_entries.append({
                    'ec_number': ec_number,
                    'alphafold_id': entry['alphafold_id'],
                    'uniprot_id': entry['uniprot_id']
                })
    else:
        # 模式2: 直接从UniProt查询
        print("  模式: 直接从UniProt查询所有纳米酶EC号")
        for i, ec_number in enumerate(ec_numbers, 1):
            print(f"  [{i}/{len(ec_numbers)}] ", end='')
            entries = query_ec_from_uniprot(fetcher, ec_number, max_per_ec)
            for entry in entries:
                all_entries.append({
                    'ec_number': ec_number,
                    'alphafold_id': entry['alphafold_id'],
                    'uniprot_id': entry['uniprot_id']
                })
            # 查询之间延迟
            if i < len(ec_numbers):
                time.sleep(1.0)
    
    total_entries = len(all_entries)
    print(f"\n  总计找到 {total_entries} 个条目")
    
    # 步骤2: 检查已存在的PDB文件
    print("\n[步骤2] 检查已存在的PDB文件...")
    to_download = []
    already_exists = 0
    
    for entry in all_entries:
        alphafold_id = entry['alphafold_id']
        if check_existing_pdb(pdb_cache_dir, alphafold_id):
            already_exists += 1
        else:
            to_download.append(entry)
    
    print(f"  已存在: {already_exists} 个")
    print(f"  需要下载: {len(to_download)} 个")
    
    if not to_download:
        print("\n✓ 所有PDB文件已存在，无需下载！")
        return
    
    # 不限制下载数量 - 下载所有可用文件
    # 移除所有限制，下载全部
    
    # 步骤3: 开始下载
    print("\n[步骤3] 开始批量下载PDB文件...")
    print("=" * 80)
    
    success_count = 0
    fail_count = 0
    failed_ids = []
    
    for i, entry in enumerate(to_download, 1):
        ec_number = entry['ec_number']
        alphafold_id = entry['alphafold_id']
        uniprot_id = entry['uniprot_id']
        
        print(f"[{i}/{len(to_download)}] EC {ec_number} - {alphafold_id} ({uniprot_id})...", end=' ', flush=True)
        
        try:
            pdb_path = fetcher.download_pdb(alphafold_id)
            
            if pdb_path and pdb_path.exists():
                print("✓ 成功")
                success_count += 1
            else:
                print("✗ 失败（文件不存在）")
                fail_count += 1
                failed_ids.append(alphafold_id)
        
        except Exception as e:
            print(f"✗ 失败: {e}")
            fail_count += 1
            failed_ids.append(alphafold_id)
        
        # 延迟，避免请求过快
        if i < len(to_download):
            time.sleep(delay)
    
    # 输出统计信息
    print("\n" + "=" * 80)
    print("下载完成统计:")
    print(f"  成功: {success_count} 个")
    print(f"  失败: {fail_count} 个")
    print(f"  总计: {len(to_download)} 个")
    print(f"  PDB缓存目录: {pdb_cache_dir}")
    
    if failed_ids:
        print(f"\n失败的AlphaFold ID (前10个):")
        for fid in failed_ids[:10]:
            print(f"  - {fid}")
        if len(failed_ids) > 10:
            print(f"  ... 还有 {len(failed_ids) - 10} 个")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='批量下载PDB文件')
    parser.add_argument(
        '--json-cache',
        type=str,
        default='cache/json',
        help='JSON缓存目录路径（默认: cache/json）'
    )
    parser.add_argument(
        '--pdb-cache',
        type=str,
        default='cache/pdb',
        help='PDB缓存目录路径（默认: cache/pdb）'
    )
    parser.add_argument(
        '--max-entries',
        type=int,
        default=None,
        help='最大下载条目数（默认: 无限制，下载全部）'
    )
    parser.add_argument(
        '--max-per-ec',
        type=int,
        default=-1,
        help='每个EC号最大查询数量（默认: -1表示全部，无限制）'
    )
    parser.add_argument(
        '--no-cache',
        action='store_true',
        help='不使用JSON缓存，直接从UniProt查询所有EC号'
    )
    parser.add_argument(
        '--cache-dir',
        type=str,
        default='./cache',
        help='缓存根目录（默认: ./cache）'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=0.5,
        help='每次下载之间的延迟（秒，默认: 0.5）'
    )
    
    args = parser.parse_args()
    
    # 转换为Path对象
    project_root = Path(__file__).parent.parent
    json_cache_dir = project_root / args.json_cache if args.json_cache else None
    pdb_cache_dir = project_root / args.pdb_cache if args.pdb_cache else None
    cache_dir = project_root / args.cache_dir if args.cache_dir else project_root / "./cache"
    
    # 如果使用缓存模式，检查目录
    if not args.no_cache and json_cache_dir and not json_cache_dir.exists():
        print(f"警告: JSON缓存目录不存在: {json_cache_dir}")
        print("  将切换到直接从UniProt查询模式")
        args.no_cache = True
    
    # 执行批量下载（无限制模式）
    batch_download_pdb(
        json_cache_dir=json_cache_dir,
        pdb_cache_dir=pdb_cache_dir,
        cache_dir=str(cache_dir),
        max_entries=args.max_entries,  # None表示无限制
        max_per_ec=args.max_per_ec if args.max_per_ec > 0 else FETCH_ALL_RESULTS,
        delay=args.delay,
        use_cache=not args.no_cache
    )


if __name__ == '__main__':
    main()


