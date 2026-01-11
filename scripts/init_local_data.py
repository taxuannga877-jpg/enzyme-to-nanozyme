#!/usr/bin/env python3
"""
本地数据初始化脚本
一键完成：
1. 批量下载所有缺失的PDB文件
2. 构建Motif数据库索引
"""

import sys
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.batch_download_pdb import batch_download_pdb
from enzyme_viewer.motif_db import build_motif_index


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='初始化本地数据（PDB下载 + Motif索引）')
    parser.add_argument(
        '--skip-pdb',
        action='store_true',
        help='跳过PDB下载步骤'
    )
    parser.add_argument(
        '--skip-motif',
        action='store_true',
        help='跳过Motif索引构建步骤'
    )
    parser.add_argument(
        '--max-pdb',
        type=int,
        default=None,
        help='最大PDB下载数量（默认: 全部）'
    )
    parser.add_argument(
        '--pdb-delay',
        type=float,
        default=0.5,
        help='PDB下载延迟（秒，默认: 0.5）'
    )
    parser.add_argument(
        '--clear-motif-db',
        action='store_true',
        help='清空现有Motif数据库后重建'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("本地数据初始化")
    print("=" * 60)
    
    # 步骤1: 批量下载PDB
    if not args.skip_pdb:
        print("\n" + "=" * 60)
        print("步骤1: 批量下载PDB文件")
        print("=" * 60)
        
        json_cache_dir = project_root / 'cache' / 'json'
        pdb_cache_dir = project_root / 'cache' / 'pdb'
        
        if not json_cache_dir.exists():
            print(f"错误: JSON缓存目录不存在: {json_cache_dir}")
            print("请先运行EC号查询以生成JSON缓存文件")
            sys.exit(1)
        
        batch_download_pdb(
            json_cache_dir=json_cache_dir,
            pdb_cache_dir=pdb_cache_dir,
            max_entries=args.max_pdb,
            delay=args.pdb_delay
        )
    else:
        print("\n跳过PDB下载步骤")
    
    # 步骤2: 构建Motif索引
    if not args.skip_motif:
        print("\n" + "=" * 60)
        print("步骤2: 构建Motif数据库索引")
        print("=" * 60)
        
        motif_library_dir = project_root / 'motif_library'
        db_path = project_root / 'enzyme_viewer' / 'motif_index.db'
        
        if not motif_library_dir.exists():
            print(f"警告: Motif库目录不存在: {motif_library_dir}")
            print("跳过Motif索引构建")
        else:
            build_motif_index(
                motif_library_dir=motif_library_dir,
                db_path=db_path,
                clear_existing=args.clear_motif_db
            )
    else:
        print("\n跳过Motif索引构建步骤")
    
    print("\n" + "=" * 60)
    print("✓ 本地数据初始化完成！")
    print("=" * 60)
    print("\n现在可以启动Flask应用，系统将优先使用本地数据：")
    print("  - PDB文件: cache/pdb/")
    print("  - Motif索引: enzyme_viewer/motif_index.db")
    print("\n启动应用:")
    print("  cd enzyme_viewer && python app.py")


if __name__ == '__main__':
    main()


