#!/usr/bin/env python3
"""
升级PDB库脚本 - 将AlphaFold PDB替换为实验PDB
===========================================
遍历所有EC号的JSON文件，检查哪些条目有实验PDB ID，
然后下载实验PDB文件替换现有的AlphaFold PDB文件
"""

import json
import sys
import shutil
from pathlib import Path
from typing import Dict, List, Optional
import time
import requests

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE
from nanozyme_mining.utils.ec_mappings import EC_PATTERNS
from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher

# RCSB PDB下载URL
RCSB_PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"

def get_all_nanozyme_ec_numbers() -> List[str]:
    """获取所有纳米酶EC号列表"""
    all_ecs = []
    for nanozyme_type, ec_list in EC_PATTERNS.items():
        all_ecs.extend(ec_list)
    return sorted(set(all_ecs))


def download_experimental_pdb(
    pdb_id: str,
    target_dir: Path,
    timeout: int = 15
) -> Optional[Path]:
    """
    下载实验PDB文件（从RCSB PDB）到指定的库目录
    
    Args:
        pdb_id: PDB ID (e.g., "1ABC")
        target_dir: 目标目录
        timeout: 下载超时时间（秒）
    
    Returns:
        下载的PDB文件路径，失败返回None
    """
    pdb_id = pdb_id.upper().strip()
    if not pdb_id:
        return None
        
    pdb_filename = f"{pdb_id}.pdb"
    pdb_file = target_dir / pdb_filename
    
    # 如果已存在，直接返回
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        return pdb_file
    
    # 下载PDB文件
    url = RCSB_PDB_URL_TEMPLATE.format(pdb_id=pdb_id)
    
    try:
        # 使用较短的超时时间，避免卡住
        response = requests.get(url, timeout=timeout, stream=False)
        response.raise_for_status()
        
        # 检查响应是否有效（不是404页面）
        if len(response.text) < 1000:
            return None
        
        # 写入文件
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        
        return pdb_file
    
    except requests.Timeout:
        return None
    except requests.RequestException as e:
        return None
    except Exception as e:
        return None


def find_alphafold_pdb(
    uniprot_id: str,
    alphafold_id: str,
    ec_dir: Path
) -> Optional[Path]:
    """
    在EC目录中查找对应的AlphaFold PDB文件
    
    Args:
        uniprot_id: UniProt ID
        alphafold_id: AlphaFold ID
        ec_dir: EC号目录
    
    Returns:
        找到的PDB文件路径，未找到返回None
    """
    AFDB_VERSION = "6"
    
    # 方法1: 使用alphafold_id查找
    if alphafold_id:
        pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
        pdb_path = ec_dir / pdb_filename
        if pdb_path.exists():
            return pdb_path
        
        # 模糊匹配
        matching_pdb = list(ec_dir.glob(f"AF-{alphafold_id}-*.pdb"))
        if matching_pdb:
            return matching_pdb[0]
    
    # 方法2: 使用uniprot_id查找
    if uniprot_id and uniprot_id != alphafold_id:
        pdb_filename = f"AF-{uniprot_id}-F1-model_v{AFDB_VERSION}.pdb"
        pdb_path = ec_dir / pdb_filename
        if pdb_path.exists():
            return pdb_path
        
        # 模糊匹配
        matching_pdb = list(ec_dir.glob(f"AF-{uniprot_id}-*.pdb"))
        if matching_pdb:
            return matching_pdb[0]
    
    return None


def get_pdb_id_from_uniprot(uniprot_id: str, fetcher: UniProtFetcher) -> Optional[str]:
    """
    从UniProt查询单个UniProt ID的PDB ID
    
    Args:
        uniprot_id: UniProt ID
        fetcher: UniProtFetcher实例
    
    Returns:
        PDB ID，如果不存在则返回None
    """
    try:
        # 使用UniProt REST API查询单个条目
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url, timeout=10)  # 减少超时时间
        if response.status_code == 200:
            data = response.json()
            # 使用fetcher的_extract_pdb_id方法提取PDB ID
            pdb_id = fetcher._extract_pdb_id(data)
            return pdb_id
    except requests.Timeout:
        pass  # 超时，静默失败
    except Exception as e:
        pass  # 其他错误，静默失败
    return None


def upgrade_ec_directory(
    ec_number: str,
    library_dir: Path,
    fetcher: Optional[UniProtFetcher] = None,
    delay: float = 0.5,
    refetch_missing_pdb_id: bool = True
) -> Dict:
    """
    升级单个EC号目录的PDB文件
    
    Args:
        ec_number: EC号
        library_dir: PDB库目录
        delay: 每次下载之间的延迟（秒）
    
    Returns:
        升级统计信息
    """
    ec_dir_name = ec_number.replace('.', '_')
    ec_dir = library_dir / ec_dir_name
    json_file = ec_dir / f"{ec_number}_sites.json"
    
    stats = {
        'ec_number': ec_number,
        'total_entries': 0,
        'entries_with_pdb_id': 0,
        'experimental_pdb_downloaded': 0,
        'experimental_pdb_already_exists': 0,
        'experimental_pdb_failed': 0,
        'alphafold_pdb_replaced': 0,
        'alphafold_pdb_not_found': 0
    }
    
    if not json_file.exists():
        print(f"  ⚠️  JSON文件不存在: {json_file}")
        return stats
    
    # 读取JSON文件
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            entries = json.load(f)
    except Exception as e:
        print(f"  ✗ 读取JSON文件失败: {e}")
        return stats
    
    stats['total_entries'] = len(entries)
    
    print(f"\n  EC {ec_number} - 处理 {len(entries)} 个条目")
    print("  " + "-" * 76)
    
    for i, entry in enumerate(entries, 1):
        pdb_id = entry.get('pdb_id', '')
        if pdb_id:
            pdb_id = str(pdb_id).strip().upper()
        uniprot_id = entry.get('uniprot_id', '')
        alphafold_id = entry.get('alphafold_id', '')
        name = entry.get('name', '')[:50]  # 限制显示长度
        
        # 如果JSON中没有pdb_id，尝试从UniProt查询（不输出，避免过多日志）
        if not pdb_id and refetch_missing_pdb_id and fetcher and uniprot_id:
            try:
                pdb_id = get_pdb_id_from_uniprot(uniprot_id, fetcher)
                if pdb_id:
                    pdb_id = str(pdb_id).strip().upper()
                    # 更新JSON条目（可选，但不写入文件以避免修改原始数据）
                    entry['pdb_id'] = pdb_id
            except Exception:
                # 查询失败，静默继续
                pass
        
        if not pdb_id:
            continue  # 没有实验PDB ID，跳过
        
        stats['entries_with_pdb_id'] += 1
        
        print(f"  [{i}/{len(entries)}] {pdb_id} ({uniprot_id or alphafold_id}) {name}...", end=' ', flush=True)
        
        # 检查实验PDB是否已存在
        exp_pdb_file = ec_dir / f"{pdb_id}.pdb"
        if exp_pdb_file.exists() and exp_pdb_file.stat().st_size > 0:
            print("✓ 实验PDB已存在")
            stats['experimental_pdb_already_exists'] += 1
            continue
        
        # 下载实验PDB（添加进度提示）
        print("下载中...", end=' ', flush=True)
        downloaded_file = download_experimental_pdb(pdb_id, ec_dir, timeout=10)
        
        if downloaded_file and downloaded_file.exists():
            print("✓ 下载成功", end='')
            stats['experimental_pdb_downloaded'] += 1
            
            # 查找并删除对应的AlphaFold PDB（如果存在）
            af_pdb_file = find_alphafold_pdb(uniprot_id, alphafold_id, ec_dir)
            if af_pdb_file and af_pdb_file.exists():
                try:
                    # 备份AlphaFold PDB到backup目录（可选）
                    # backup_dir = ec_dir / "backup_alphafold"
                    # backup_dir.mkdir(exist_ok=True)
                    # shutil.copy2(af_pdb_file, backup_dir / af_pdb_file.name)
                    
                    # 删除AlphaFold PDB
                    af_pdb_file.unlink()
                    print(" (已删除AlphaFold PDB)")
                    stats['alphafold_pdb_replaced'] += 1
                except Exception as e:
                    print(f" (删除AlphaFold PDB失败: {e})")
            else:
                print()
                stats['alphafold_pdb_not_found'] += 1
        else:
            print("✗ 下载失败")
            stats['experimental_pdb_failed'] += 1
        
        # 延迟，避免请求过快
        if i < len(entries):
            time.sleep(delay)
    
    return stats


def main():
    """主函数"""
    library_dir = project_root / "pdb_library"
    
    if not library_dir.exists():
        print(f"✗ PDB库目录不存在: {library_dir}")
        return
    
    print("=" * 80)
    print("PDB库升级工具 - 将AlphaFold PDB替换为实验PDB")
    print("=" * 80)
    print(f"PDB库目录: {library_dir}")
    
    # 初始化UniProtFetcher（用于查询缺失的PDB ID）
    print("\n正在初始化UniProtFetcher...")
    cache_dir = project_root / "cache"
    try:
        fetcher = UniProtFetcher(
            cache_dir=str(cache_dir),
            pdb_library_dir=str(library_dir),
            use_mcsa=False  # 禁用M-CSA以加快初始化速度
        )
        print("✓ UniProtFetcher初始化完成")
    except Exception as e:
        print(f"⚠️  UniProtFetcher初始化失败: {e}")
        print("   继续运行，但不会查询缺失的PDB ID")
        fetcher = None
    
    # 获取所有EC号
    ec_numbers = get_all_nanozyme_ec_numbers()
    print(f"处理 {len(ec_numbers)} 个纳米酶EC号")
    print("注意: 如果JSON文件中缺少pdb_id字段，将尝试从UniProt查询获取")
    print("=" * 80)
    
    all_stats = []
    
    for i, ec_number in enumerate(ec_numbers, 1):
        print(f"\n[{i}/{len(ec_numbers)}] ", end='')
        stats = upgrade_ec_directory(ec_number, library_dir, fetcher=fetcher, delay=0.5, refetch_missing_pdb_id=True)
        all_stats.append(stats)
        
        # 打印统计
        if stats['entries_with_pdb_id'] > 0:
            print(f"\n    统计: {stats['entries_with_pdb_id']} 个有实验PDB ID的条目")
            print(f"      - 已下载: {stats['experimental_pdb_downloaded']}")
            print(f"      - 已存在: {stats['experimental_pdb_already_exists']}")
            print(f"      - 失败: {stats['experimental_pdb_failed']}")
            print(f"      - 替换AlphaFold: {stats['alphafold_pdb_replaced']}")
    
    # 汇总统计
    print("\n" + "=" * 80)
    print("升级完成 - 汇总统计")
    print("=" * 80)
    
    total_entries = sum(s['total_entries'] for s in all_stats)
    total_with_pdb_id = sum(s['entries_with_pdb_id'] for s in all_stats)
    total_downloaded = sum(s['experimental_pdb_downloaded'] for s in all_stats)
    total_already_exists = sum(s['experimental_pdb_already_exists'] for s in all_stats)
    total_failed = sum(s['experimental_pdb_failed'] for s in all_stats)
    total_replaced = sum(s['alphafold_pdb_replaced'] for s in all_stats)
    
    print(f"总条目数: {total_entries}")
    print(f"有实验PDB ID的条目: {total_with_pdb_id}")
    print(f"新下载的实验PDB: {total_downloaded}")
    print(f"已存在的实验PDB: {total_already_exists}")
    print(f"下载失败: {total_failed}")
    print(f"替换的AlphaFold PDB: {total_replaced}")
    print("=" * 80)


if __name__ == "__main__":
    main()

