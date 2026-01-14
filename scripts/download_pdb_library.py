#!/usr/bin/env python3
"""
纳米酶PDB库下载脚本
==================
专门用于下载和整理所有纳米酶EC号的PDB文件到独立的PDB库文件夹
不存到缓存目录，单独管理
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

from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher, FETCH_ALL_RESULTS
from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE
from nanozyme_mining.utils.ec_mappings import EC_PATTERNS

# AlphaFold DB配置
AFDB_VERSION = "6"
ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{alphafold_id}-F1-model_v{version}.pdb"
RCSB_PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"


def get_all_nanozyme_ec_numbers() -> List[str]:
    """获取所有纳米酶EC号列表"""
    all_ecs = []
    for nanozyme_type, ec_list in EC_PATTERNS.items():
        all_ecs.extend(ec_list)
    return sorted(set(all_ecs))


def download_experimental_pdb_to_library(
    pdb_id: str,
    target_dir: Path,
    timeout: int = 30
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
    pdb_id = pdb_id.upper()
    pdb_filename = f"{pdb_id}.pdb"
    pdb_file = target_dir / pdb_filename
    
    # 如果已存在，直接返回
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        return pdb_file
    
    # 下载PDB文件
    url = RCSB_PDB_URL_TEMPLATE.format(pdb_id=pdb_id)
    
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        
        # 检查响应是否有效（不是404页面）
        if len(response.text) < 1000:
            print(f"    ⚠️  文件可能无效（太小）")
        
        # 写入文件
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        
        return pdb_file
    
    except Exception as e:
        print(f"    ✗ 下载失败: {e}")
        return None


def download_alphafold_pdb_to_library(
    alphafold_id: str,
    target_dir: Path,
    timeout: int = 30
) -> Optional[Path]:
    """
    下载AlphaFold PDB文件到指定的库目录
    
    Args:
        alphafold_id: AlphaFold ID
        target_dir: 目标目录
        timeout: 下载超时时间（秒）
    
    Returns:
        下载的PDB文件路径，失败返回None
    """
    pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
    pdb_file = target_dir / pdb_filename
    
    # 如果已存在，直接返回
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        return pdb_file
    
    # 下载PDB文件
    url = ALPHAFOLD_URL_TEMPLATE.format(
        alphafold_id=alphafold_id,
        version=AFDB_VERSION
    )
    
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        
        # 写入文件
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        
        return pdb_file
    
    except Exception as e:
        print(f"    ✗ 下载失败: {e}")
        return None


def download_pdb_to_library(
    pdb_id: Optional[str] = None,
    alphafold_id: Optional[str] = None,
    target_dir: Path = None,
    timeout: int = 30,
    replace_alphafold: bool = True
) -> Optional[Path]:
    """
    下载PDB文件到指定的库目录（优先使用实验PDB，如果没有则使用AlphaFold）
    如果成功下载了实验PDB，可以选择删除对应的AlphaFold PDB文件
    
    Args:
        pdb_id: 实验PDB ID (可选，优先使用)
        alphafold_id: AlphaFold ID (可选，回退选项)
        target_dir: 目标目录
        timeout: 下载超时时间（秒）
        replace_alphafold: 如果下载了实验PDB，是否删除对应的AlphaFold PDB（默认True）
    
    Returns:
        下载的PDB文件路径，失败返回None
    """
    # 优先级1: 尝试下载实验PDB
    if pdb_id:
        pdb_path = download_experimental_pdb_to_library(pdb_id, target_dir, timeout)
        if pdb_path and pdb_path.exists():
            # 如果成功下载了实验PDB，删除对应的AlphaFold PDB（如果存在）
            if replace_alphafold and alphafold_id:
                af_pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
                af_pdb_path = target_dir / af_pdb_filename
                if af_pdb_path.exists():
                    try:
                        af_pdb_path.unlink()
                    except Exception as e:
                        pass  # 忽略删除错误，继续执行
            return pdb_path
    
    # 优先级2: 回退到AlphaFold
    if alphafold_id:
        return download_alphafold_pdb_to_library(alphafold_id, target_dir, timeout)
    
        return None


def query_ec_from_uniprot(
    fetcher: UniProtFetcher,
    ec_number: str,
    library_dir: Path,
    cache_dir: str = "./cache",
    max_per_ec: int = FETCH_ALL_RESULTS
) -> List[Dict]:
    """
    从缓存文件或UniProt查询EC号的所有条目（包含PDB ID和AlphaFold ID）
    优先从缓存文件读取，如果缓存不存在或为空，才查询UniProt
    
    Args:
        fetcher: UniProtFetcher实例
        ec_number: EC号
        library_dir: PDB库目录
        cache_dir: 缓存目录路径
        max_per_ec: 每个EC号最大查询数量（-1表示全部）
    
    Returns:
        包含pdb_id、alphafold_id和uniprot_id的条目列表
    """
    print(f"  查询 EC {ec_number}...", end=' ', flush=True)
    
    # 先尝试从缓存文件读取
    # 优先从pdb_library读取，向后兼容旧路径
    ec_dir_name = ec_number.replace(".", "_")
    cache_file = library_dir / ec_dir_name / f"{ec_number}_sites.json"
    if not cache_file.exists() and cache_dir:
        cache_file = Path(cache_dir) / "json" / f"{ec_number}_sites.json"
    if cache_file.exists():
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                cached_entries = json.load(f)
            
            # 从缓存数据提取pdb_id和alphafold_id
            result = []
            for entry in cached_entries:
                pdb_id = entry.get('pdb_id', '')
                if pdb_id:
                    pdb_id = str(pdb_id).strip().upper()
                alphafold_id = entry.get('alphafold_id', '')
                if alphafold_id:
                    alphafold_id = str(alphafold_id).strip()
                uniprot_id = entry.get('uniprot_id', '')
                if uniprot_id:
                    uniprot_id = str(uniprot_id).strip()
                
                # 保留有PDB ID或AlphaFold结构的条目
                if pdb_id or alphafold_id:
                    result.append({
                        'pdb_id': pdb_id,
                        'alphafold_id': alphafold_id,
                        'uniprot_id': uniprot_id,
                        'name': entry.get('name', ''),
                        'organism': entry.get('organism', '')
                    })
            
            if result:
                pdb_count = sum(1 for e in result if e['pdb_id'])
                af_count = sum(1 for e in result if e['alphafold_id'])
                print(f"从缓存读取: {len(result)} 个条目（{pdb_count} 个有实验PDB，{af_count} 个有AlphaFold）")
                return result
        except Exception as e:
            print(f"读取缓存失败: {e}, 将查询UniProt...", end=' ', flush=True)
    
    # 缓存不存在或为空，从UniProt查询
    entries = fetcher.query_with_active_sites(ec_number, max_results=max_per_ec)
    
    # 提取pdb_id和alphafold_id
    result = []
    for entry in entries:
        pdb_id = entry.get('pdb_id', '')
        if pdb_id:
            pdb_id = str(pdb_id).strip().upper()
        alphafold_id = entry.get('alphafold_id', '')
        if alphafold_id:
            alphafold_id = str(alphafold_id).strip()
        uniprot_id = entry.get('uniprot_id', '')
        if uniprot_id:
            uniprot_id = str(uniprot_id).strip()
        
        # 保留有PDB ID或AlphaFold结构的条目
        if pdb_id or alphafold_id:
            result.append({
                'pdb_id': pdb_id,
                'alphafold_id': alphafold_id,
                'uniprot_id': uniprot_id,
                'name': entry.get('name', ''),
                'organism': entry.get('organism', '')
            })
    
    pdb_count = sum(1 for e in result if e['pdb_id'])
    af_count = sum(1 for e in result if e['alphafold_id'])
    print(f"从UniProt查询: {len(result)} 个条目（{pdb_count} 个有实验PDB，{af_count} 个有AlphaFold）")
    return result


def download_pdb_library(
    library_dir: Path = None,
    cache_dir: str = "./cache",
    delay: float = 0.5
):
    """
    下载所有纳米酶EC号的PDB文件到专门的库文件夹
    
    Args:
        library_dir: PDB库目录（默认: ./pdb_library）
        cache_dir: 临时缓存目录（用于UniProt查询）
        delay: 每次下载之间的延迟（秒）
    """
    # 设置PDB库目录
    if library_dir is None:
        library_dir = project_root / "pdb_library"
    library_dir = Path(library_dir)
    
    # 创建库目录结构
    library_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建按EC号组织的子目录
    ec_numbers = get_all_nanozyme_ec_numbers()
    
    print("=" * 80)
    print("纳米酶PDB库下载工具")
    print("=" * 80)
    print(f"PDB库目录: {library_dir}")
    print(f"处理 {len(ec_numbers)} 个纳米酶EC号")
    print("=" * 80)
    
    # 初始化UniProtFetcher（仅用于查询，不用于下载）
    fetcher = UniProtFetcher(cache_dir=cache_dir)
    
    # 存储所有下载信息
    library_index = {
        'total_ec_numbers': len(ec_numbers),
        'total_pdb_files': 0,
        'ec_entries': {}
    }
    
    # 步骤1: 查询所有EC号的条目
    print("\n[步骤1] 查询所有EC号的条目...")
    all_entries = {}
    
    for i, ec_number in enumerate(ec_numbers, 1):
        print(f"  [{i}/{len(ec_numbers)}] ", end='')
        entries = query_ec_from_uniprot(fetcher, ec_number, library_dir, cache_dir, FETCH_ALL_RESULTS)
        all_entries[ec_number] = entries
        
        # 查询之间延迟
        if i < len(ec_numbers):
            time.sleep(1.0)
    
    total_entries = sum(len(entries) for entries in all_entries.values())
    print(f"\n  总计找到 {total_entries} 个有AlphaFold结构的条目")
    
    # 步骤2: 下载所有PDB文件
    print("\n[步骤2] 下载所有PDB文件...")
    print("=" * 80)
    
    success_count = 0
    fail_count = 0
    failed_ids = []
    
    for ec_number, entries in all_entries.items():
        # 创建EC号子目录
        ec_dir = library_dir / ec_number.replace('.', '_')
        ec_dir.mkdir(parents=True, exist_ok=True)
        
        nanozyme_type = EC_TO_NANOZYME_TYPE.get(ec_number, 'UNKNOWN')
        print(f"\nEC {ec_number} ({nanozyme_type.value}) - {len(entries)} 个条目")
        print("-" * 80)
        
        ec_info = {
            'ec_number': ec_number,
            'nanozyme_type': nanozyme_type.value if hasattr(nanozyme_type, 'value') else str(nanozyme_type),
            'total_entries': len(entries),
            'successful_downloads': 0,
            'failed_downloads': 0,
            'pdb_files': []
        }
        
        for i, entry in enumerate(entries, 1):
            pdb_id = entry.get('pdb_id', '')
            alphafold_id = entry.get('alphafold_id', '')
            uniprot_id = entry['uniprot_id']
            
            # 显示标识符
            if pdb_id:
                display_id = f"{pdb_id} (exp)"
            elif alphafold_id:
                display_id = f"{alphafold_id} (AF)"
            else:
                display_id = uniprot_id
            
            print(f"  [{i}/{len(entries)}] {display_id} ({uniprot_id})...", end=' ', flush=True)
            
            # 优先下载实验PDB，如果没有则下载AlphaFold
            pdb_path = download_pdb_to_library(
                pdb_id=pdb_id,
                alphafold_id=alphafold_id,
                target_dir=ec_dir,
                timeout=30
            )
            
            if pdb_path and pdb_path.exists():
                source_type = "experimental" if pdb_id else "alphafold"
                print(f"✓ ({source_type})")
                success_count += 1
                ec_info['successful_downloads'] += 1
                ec_info['pdb_files'].append({
                    'pdb_id': pdb_id,
                    'alphafold_id': alphafold_id,
                    'uniprot_id': uniprot_id,
                    'source_type': source_type,
                    'name': entry.get('name', ''),
                    'organism': entry.get('organism', ''),
                    'pdb_file': pdb_path.name,
                    'pdb_path': str(pdb_path.relative_to(library_dir))
                })
            else:
                print("✗")
                fail_count += 1
                failed_ids.append(display_id)
                ec_info['failed_downloads'] += 1
            
            # 下载之间延迟
            if i < len(entries):
                time.sleep(delay)
        
        library_index['ec_entries'][ec_number] = ec_info
        library_index['total_pdb_files'] += ec_info['successful_downloads']
    
    # 步骤3: 保存索引文件
    print("\n[步骤3] 保存索引文件...")
    index_file = library_dir / "library_index.json"
    with open(index_file, 'w', encoding='utf-8') as f:
        json.dump(library_index, f, indent=2, ensure_ascii=False)
    print(f"  ✓ 索引文件已保存: {index_file}")
    
    # 步骤4: 创建README文件
    readme_file = library_dir / "README.md"
    readme_content = f"""# 纳米酶PDB库

## 概述

本目录包含所有纳米酶EC号的AlphaFold结构文件（PDB格式）。

## 统计信息

- **EC号总数**: {library_index['total_ec_numbers']}
- **PDB文件总数**: {library_index['total_pdb_files']}
- **成功下载**: {success_count}
- **下载失败**: {fail_count}

## 目录结构

```
pdb_library/
├── README.md                    # 本文件
├── library_index.json           # 详细索引文件
├── 1_1_3_4/                     # EC 1.1.3.4 (Glucose Oxidase)
├── 1_10_3_2/                    # EC 1.10.3.2 (Laccase)
├── 1_11_1_6/                    # EC 1.11.1.6 (Catalase)
├── ...
```

## EC号列表

"""
    
    for ec_number, ec_info in sorted(library_index['ec_entries'].items()):
        nanozyme_type = ec_info['nanozyme_type']
        total = ec_info['total_entries']
        success = ec_info['successful_downloads']
        failed = ec_info['failed_downloads']
        readme_content += f"- **{ec_number}** ({nanozyme_type}): {success}/{total} 成功, {failed} 失败\n"
    
    readme_content += f"""
## 使用说明

1. 所有PDB文件按EC号组织在各自的子目录中
2. 文件名格式: `AF-{{alphafold_id}}-F1-model_v{AFDB_VERSION}.pdb`
3. 详细索引信息请查看 `library_index.json`

## 更新日期

{time.strftime('%Y-%m-%d %H:%M:%S')}
"""
    
    with open(readme_file, 'w', encoding='utf-8') as f:
        f.write(readme_content)
    print(f"  ✓ README文件已保存: {readme_file}")
    
    # 输出最终统计
    print("\n" + "=" * 80)
    print("下载完成统计:")
    print(f"  成功: {success_count} 个")
    print(f"  失败: {fail_count} 个")
    print(f"  总计: {total_entries} 个")
    print(f"  PDB库目录: {library_dir}")
    print("=" * 80)
    
    if failed_ids:
        print(f"\n失败的AlphaFold ID (前10个):")
        for fid in failed_ids[:10]:
            print(f"  - {fid}")
        if len(failed_ids) > 10:
            print(f"  ... 还有 {len(failed_ids) - 10} 个")
    
    print(f"\n✓ PDB库已创建完成！")
    print(f"  所有文件保存在: {library_dir}")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='下载纳米酶PDB库到专门的文件夹')
    parser.add_argument(
        '--library-dir',
        type=str,
        default='./pdb_library',
        help='PDB库目录路径（默认: ./pdb_library）'
    )
    parser.add_argument(
        '--cache-dir',
        type=str,
        default='./cache',
        help='临时缓存目录（用于UniProt查询，默认: ./cache）'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=0.5,
        help='每次下载之间的延迟（秒，默认: 0.5）'
    )
    
    args = parser.parse_args()
    
    # 转换为Path对象
    library_dir = project_root / args.library_dir if not Path(args.library_dir).is_absolute() else Path(args.library_dir)
    cache_dir = project_root / args.cache_dir if not Path(args.cache_dir).is_absolute() else Path(args.cache_dir)
    
    # 执行下载
    download_pdb_library(
        library_dir=library_dir,
        cache_dir=str(cache_dir),
        delay=args.delay
    )


if __name__ == '__main__':
    main()

