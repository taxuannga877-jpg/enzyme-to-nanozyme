#!/usr/bin/env python3
"""
构建实验PDB库脚本
=================
整合多个数据源（UniProt、RCSB PDB等）查询所有纳米酶EC号的实验PDB结构
自动去重后，只下载有详细结构信息的实验PDB文件到对应的EC目录
生成只包含实验PDB信息的JSON文件

注意：
- 不区分数据来源，只要是有详细信息的酶PDB就下载
- 只下载包含足够结构信息（ATOM/HETATM记录）的PDB文件
- 实验PDB信息量通常多于AlphaFold预测结构（包含配体、辅因子等）
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
import requests

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher, FETCH_ALL_RESULTS
from nanozyme_mining.utils.constants import EC_TO_NANOZYME_TYPE
from nanozyme_mining.utils.ec_mappings import EC_PATTERNS

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
    timeout: int = 60,
    max_retries: int = 3
) -> Optional[Path]:
    """
    下载实验PDB文件（从RCSB PDB）到指定的库目录
    只下载有详细结构信息的PDB文件（包含ATOM记录）
    
    Args:
        pdb_id: PDB ID (e.g., "1ABC")
        target_dir: 目标目录
        timeout: 下载超时时间（秒，默认60秒）
        max_retries: 最大重试次数（默认3次）
    
    Returns:
        下载的PDB文件路径，失败或无效返回None
    """
    pdb_id = pdb_id.upper().strip()
    if not pdb_id:
        return None
        
    pdb_filename = f"{pdb_id}.pdb"
    pdb_file = target_dir / pdb_filename
    
    # 如果已存在，验证文件是否有效
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        # 验证文件是否包含结构信息
        if _is_valid_pdb_file(pdb_file):
            return pdb_file
        else:
            # 文件无效，删除并重新下载
            pdb_file.unlink()
    
    # 下载PDB文件（带重试机制）
    url = RCSB_PDB_URL_TEMPLATE.format(pdb_id=pdb_id)
    
    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            # 使用stream模式下载，避免大文件超时
            response = requests.get(url, timeout=timeout, stream=True)
            response.raise_for_status()
            
            # 流式读取内容
            content_parts = []
            for chunk in response.iter_content(chunk_size=8192, decode_unicode=True):
                if chunk:
                    content_parts.append(chunk)
            
            content = ''.join(content_parts)
            
            # 检查响应是否有效（不是404页面）
            if len(content) < 1000:
                last_error = f"文件太小（{len(content)} 字节），可能无效"
                continue
            
            # 验证内容是否包含结构信息
            if not _has_structure_info(content):
                last_error = "文件不包含足够的ATOM/HETATM记录"
                continue
            
            # 写入文件
            with open(pdb_file, 'w', encoding='utf-8') as f:
                f.write(content)
            
            # 再次验证写入的文件
            if _is_valid_pdb_file(pdb_file):
                return pdb_file
            else:
                # 文件无效，删除
                if pdb_file.exists():
                    pdb_file.unlink()
                last_error = "文件验证失败"
                continue
        
        except requests.Timeout:
            last_error = f"请求超时（{timeout}秒）"
            if attempt < max_retries:
                time.sleep(2 * attempt)  # 指数退避
            continue
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                last_error = f"PDB ID不存在（404）"
                break  # 404不需要重试
            else:
                last_error = f"HTTP错误 {e.response.status_code}"
                if attempt < max_retries:
                    time.sleep(2 * attempt)
                continue
        except requests.ConnectionError as e:
            last_error = f"连接错误: {str(e)[:100]}"
            if attempt < max_retries:
                time.sleep(2 * attempt)
            continue
        except Exception as e:
            last_error = f"未知错误: {type(e).__name__}: {str(e)[:100]}"
            if attempt < max_retries:
                time.sleep(2 * attempt)
            continue
    
    # 所有重试都失败，记录错误（但不打印，由调用者决定）
    return None


def _has_structure_info(content: str) -> bool:
    """
    检查PDB内容是否包含结构信息（ATOM或HETATM记录）
    
    Args:
        content: PDB文件内容
    
    Returns:
        如果包含结构信息返回True
    """
    # 原则上只要包含 ATOM 或 HETATM 记录就认为是“有结构信息”
    # 之前要求至少 10 条，导致一些真实但较小/截断的结构被误判为失败
    return ('ATOM  ' in content) or ('HETATM' in content)


def _is_valid_pdb_file(pdb_file: Path) -> bool:
    """
    验证PDB文件是否有效（包含结构信息）
    
    Args:
        pdb_file: PDB文件路径
    
    Returns:
        如果文件有效返回True
    """
    try:
        # 只要文件存在且非空就进一步检查内容
        if not pdb_file.exists() or pdb_file.stat().st_size <= 0:
            return False
        
        # 读取文件前几KB检查
        with open(pdb_file, 'r') as f:
            content = f.read(50000)  # 读取前50KB
        
        return _has_structure_info(content)
    except Exception:
        return False


def query_experimental_pdb_entries(
    fetcher: UniProtFetcher,
    ec_number: str,
    max_per_ec: int = FETCH_ALL_RESULTS
) -> List[Dict]:
    """
    从多个数据源查询EC号的所有实验PDB结构（UniProt + RCSB PDB）
    合并所有来源的PDB ID，去重后返回
    
    Args:
        fetcher: UniProtFetcher实例
        ec_number: EC号
        max_per_ec: 每个EC号最大查询数量（-1表示全部）
    
    Returns:
        包含所有实验PDB ID的条目列表（已去重）
    """
    print(f"  查询 EC {ec_number}...", end=' ', flush=True)
    
    # 步骤1: 从UniProt查询所有条目（包含所有PDB ID）
    entries = fetcher.query_with_active_sites(ec_number, max_results=max_per_ec)
    
    # 收集所有PDB ID（从UniProt）
    uniprot_pdb_map = {}  # pdb_id -> entry with sequence and active sites
    
    for entry in entries:
        # 获取所有PDB ID（包括pdb_ids字段）
        all_pdb_ids = entry.get('pdb_ids', [])
        if not all_pdb_ids:
            # 回退到单个pdb_id（兼容旧格式）
            single_pdb_id = entry.get('pdb_id', '')
            if single_pdb_id:
                all_pdb_ids = [single_pdb_id.strip().upper()]
        
        for pdb_id in all_pdb_ids:
            pdb_id = str(pdb_id).strip().upper()
            if pdb_id:
                # 如果还没有映射，创建新的条目
                if pdb_id not in uniprot_pdb_map:
                    uniprot_pdb_map[pdb_id] = {
                        'pdb_id': pdb_id,
                        'uniprot_id': entry.get('uniprot_id', ''),
                        'sequence': entry.get('sequence', ''),
                        'active_sites': entry.get('active_sites', [])
                    }
    
    # 步骤2: 从RCSB PDB直接查询该EC号的所有实验结构
    rcsb_pdb_ids = fetcher.query_rcsb_pdb_by_ec(ec_number)
    
    # 步骤3: 合并结果（去重）
    all_pdb_ids = set(uniprot_pdb_map.keys()).union(set(rcsb_pdb_ids))
    
    # 构建最终条目列表
    experimental_entries = []
    for pdb_id in sorted(all_pdb_ids):
        if pdb_id in uniprot_pdb_map:
            # 有UniProt信息的条目
            experimental_entries.append(uniprot_pdb_map[pdb_id])
        else:
            # 只有RCSB PDB信息，没有UniProt序列（但仍需要下载）
            experimental_entries.append({
                'pdb_id': pdb_id,
                'uniprot_id': '',  # 未知
                'sequence': '',  # 未知（可以从PDB文件中提取）
                'active_sites': []
            })
    
    print(f"✓ 找到 {len(experimental_entries)} 个实验PDB结构")
    return experimental_entries


def build_experimental_pdb_library(
    library_dir: Path = None,
    cache_dir: str = "./cache",
    delay: float = 1.0
):
    """
    构建只包含实验PDB的数据库
    
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
    
    # 获取所有EC号
    ec_numbers = get_all_nanozyme_ec_numbers()
    
    print("=" * 80)
    print("构建实验PDB库")
    print("=" * 80)
    print(f"PDB库目录: {library_dir}")
    print(f"处理 {len(ec_numbers)} 个纳米酶EC号")
    print("注意: 只下载和保存有实验PDB ID的条目")
    print("=" * 80)
    
    # 初始化UniProtFetcher
    fetcher = UniProtFetcher(
        cache_dir=cache_dir,
        pdb_library_dir=str(library_dir),
        use_mcsa=False  # 禁用M-CSA以加快速度
    )
    
    # 存储所有统计信息
    library_index = {
        'total_ec_numbers': len(ec_numbers),
        'total_pdb_files': 0,
        'total_entries': 0,
        'ec_entries': {}
    }
    
    all_stats = {
        'total_queried': 0,
        'with_experimental_pdb': 0,
        'successful_downloads': 0,
        'failed_downloads': 0,
        'already_exists': 0
    }
    
    # 步骤1: 查询所有EC号的有实验PDB的条目
    print("\n[步骤1] 查询所有EC号的有实验PDB的条目...")
    all_experimental_entries = {}
    
    for i, ec_number in enumerate(ec_numbers, 1):
        print(f"  [{i}/{len(ec_numbers)}] ", end='')
        entries = query_experimental_pdb_entries(fetcher, ec_number, FETCH_ALL_RESULTS)
        all_experimental_entries[ec_number] = entries
        
        all_stats['total_queried'] += len(entries)
        all_stats['with_experimental_pdb'] += len(entries)
        
        # 查询之间延迟
        if i < len(ec_numbers):
            time.sleep(1.0)
    
    print(f"\n  总计找到 {all_stats['with_experimental_pdb']} 个有实验PDB的条目")
    
    # 步骤2: 下载所有实验PDB文件并生成JSON
    print("\n[步骤2] 下载实验PDB文件并生成JSON...")
    print("=" * 80)
    
    for ec_number, entries in all_experimental_entries.items():
        if not entries:
            print(f"\nEC {ec_number}: 没有实验PDB条目，跳过")
            continue
        
        # 创建EC号子目录
        ec_dir = library_dir / ec_number.replace('.', '_')
        ec_dir.mkdir(parents=True, exist_ok=True)
        
        nanozyme_type = EC_TO_NANOZYME_TYPE.get(ec_number, 'UNKNOWN')
        print(f"\nEC {ec_number} ({nanozyme_type.value if hasattr(nanozyme_type, 'value') else str(nanozyme_type)}) - {len(entries)} 个条目")
        print("-" * 80)
        
        ec_info = {
            'ec_number': ec_number,
            'nanozyme_type': nanozyme_type.value if hasattr(nanozyme_type, 'value') else str(nanozyme_type),
            'total_entries': len(entries),
            'successful_downloads': 0,
            'failed_downloads': 0,
            'already_exists': 0
        }
        
        # 存储每个EC号的JSON条目（只包含实验PDB）
        json_entries = []
        
        for i, entry in enumerate(entries, 1):
            pdb_id = entry.get('pdb_id', '').strip().upper()
            uniprot_id = entry.get('uniprot_id', '')
            
            if not pdb_id:
                continue  # 跳过没有PDB ID的条目（理论上不应该有）
            
            # 显示信息
            if uniprot_id:
                print(f"  [{i}/{len(entries)}] {pdb_id} (UniProt: {uniprot_id})...", end=' ', flush=True)
            else:
                print(f"  [{i}/{len(entries)}] {pdb_id}...", end=' ', flush=True)
            
            # 检查PDB是否已存在
            exp_pdb_file = ec_dir / f"{pdb_id}.pdb"
            if exp_pdb_file.exists() and exp_pdb_file.stat().st_size > 0:
                print("✓ 已存在")
                all_stats['already_exists'] += 1
                ec_info['already_exists'] += 1
            else:
                # 下载实验PDB（使用更长的超时时间和重试机制）
                print("下载中...", end=' ', flush=True)
                downloaded_file = download_experimental_pdb(pdb_id, ec_dir, timeout=60, max_retries=3)
                
                if downloaded_file and downloaded_file.exists():
                    print("✓ 成功")
                    all_stats['successful_downloads'] += 1
                    ec_info['successful_downloads'] += 1
                else:
                    print("✗ 失败（可能原因：超时/网络问题/PDB不存在）")
                    all_stats['failed_downloads'] += 1
                    ec_info['failed_downloads'] += 1
                    # 即使下载失败，也添加到JSON中（标记下载状态）
                    continue  # 跳过下载失败的条目，不加入JSON
            
            # 添加到JSON条目列表（只包含成功下载或已存在的）
            json_entries.append({
                'uniprot_id': uniprot_id,
                'pdb_id': pdb_id,
                'sequence': entry.get('sequence', ''),
                'active_sites': entry.get('active_sites', [])
            })
            
            # 下载之间延迟
            if i < len(entries):
                time.sleep(delay)
        
        # 保存该EC号的JSON文件（只包含有实验PDB的条目）
        json_file = ec_dir / f"{ec_number}_sites.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_entries, f, indent=2, ensure_ascii=False)
        
        print(f"\n  ✓ JSON文件已保存: {json_file} ({len(json_entries)} 个条目)")
        
        library_index['ec_entries'][ec_number] = ec_info
        library_index['total_pdb_files'] += ec_info['successful_downloads'] + ec_info['already_exists']
        library_index['total_entries'] += len(json_entries)
    
    # 步骤3: 保存索引文件
    print("\n[步骤3] 保存索引文件...")
    index_file = library_dir / "library_index.json"
    with open(index_file, 'w', encoding='utf-8') as f:
        json.dump(library_index, f, indent=2, ensure_ascii=False)
    print(f"  ✓ 索引文件已保存: {index_file}")
    
    # 步骤4: 创建README文件
    readme_file = library_dir / "README.md"
    readme_content = f"""# 纳米酶实验PDB库

## 概述

本目录包含所有纳米酶EC号的**实验PDB结构文件**（从RCSB PDB下载）。
**只包含有实验PDB ID的条目**，不包含AlphaFold预测结构。

## 统计信息

- **EC号总数**: {library_index['total_ec_numbers']}
- **有实验PDB的条目总数**: {library_index['total_entries']}
- **PDB文件总数**: {library_index['total_pdb_files']}
- **新下载**: {all_stats['successful_downloads']}
- **已存在**: {all_stats['already_exists']}
- **下载失败**: {all_stats['failed_downloads']}

**数据来源**: 整合了UniProt和RCSB PDB等多个数据源，自动去重后下载所有有详细信息的酶PDB结构。

## 目录结构

```
pdb_library/
├── README.md                    # 本文件
├── library_index.json           # 详细索引文件
├── 1_1_3_4/                     # EC 1.1.3.4
│   ├── 1.1.3.4_sites.json       # 只包含有实验PDB的条目
│   ├── 1ABC.pdb                 # 实验PDB文件
│   ├── 2DEF.pdb
│   └── ...
├── 1_10_3_2/                    # EC 1.10.3.2
│   └── ...
```

## EC号列表

"""
    
    for ec_number, ec_info in sorted(library_index['ec_entries'].items()):
        if ec_info['total_entries'] > 0:  # 只列出有数据的EC号
            nanozyme_type = ec_info['nanozyme_type']
            total = ec_info['total_entries']
            success = ec_info['successful_downloads'] + ec_info['already_exists']
            failed = ec_info['failed_downloads']
            readme_content += f"- **{ec_number}** ({nanozyme_type}): {success} 个PDB文件 ({total} 个条目)\n"
    
    readme_content += f"""
## 使用说明

1. 所有PDB文件按EC号组织在各自的子目录中
2. **只包含实验PDB文件**（格式: `{{PDB_ID}}.pdb`），不包含AlphaFold预测结构
3. 每个EC号目录下有对应的JSON文件（`{{EC_number}}_sites.json`），**只包含有实验PDB的条目**
4. 详细索引信息请查看 `library_index.json`

## JSON文件格式

每个EC号的JSON文件包含以下字段：
- `uniprot_id`: UniProt ID
- `pdb_id`: 实验PDB ID（必须字段）
- `sequence`: 蛋白质序列
- `active_sites`: 活性位点信息数组

## 更新日期

{time.strftime('%Y-%m-%d %H:%M:%S')}
"""
    
    with open(readme_file, 'w', encoding='utf-8') as f:
        f.write(readme_content)
    print(f"  ✓ README文件已保存: {readme_file}")
    
    # 输出最终统计
    print("\n" + "=" * 80)
    print("构建完成统计:")
    print(f"  查询到的条目: {all_stats['total_queried']}")
    print(f"  有实验PDB的条目: {all_stats['with_experimental_pdb']}")
    print(f"  新下载: {all_stats['successful_downloads']} 个")
    print(f"  已存在: {all_stats['already_exists']} 个")
    print(f"  下载失败: {all_stats['failed_downloads']} 个")
    print(f"  总计PDB文件: {library_index['total_pdb_files']} 个")
    print(f"  PDB库目录: {library_dir}")
    print("=" * 80)
    
    print(f"\n✓ 实验PDB库已构建完成！")
    print(f"  所有文件保存在: {library_dir}")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='构建只包含实验PDB的数据库')
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
        default=1.0,
        help='每次下载之间的延迟（秒，默认: 1.0，建议至少1秒以避免服务器限制）'
    )
    
    args = parser.parse_args()
    
    # 转换为Path对象
    library_dir = project_root / args.library_dir if not Path(args.library_dir).is_absolute() else Path(args.library_dir)
    cache_dir = project_root / args.cache_dir if not Path(args.cache_dir).is_absolute() else Path(args.cache_dir)
    
    # 执行构建
    build_experimental_pdb_library(
        library_dir=library_dir,
        cache_dir=str(cache_dir),
        delay=args.delay
    )


if __name__ == '__main__':
    main()

