#!/usr/bin/env python
"""
实时监控Motif提取进度
====================

批量处理所有pdb_library中的EC号，确保motif_library对齐
实时显示每个EC号的处理进度和统计信息
"""

import os
import sys
import subprocess
import time
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

def get_ec_numbers_from_pdb_library(pdb_library_dir: Path) -> List[str]:
    """从pdb_library获取所有EC号"""
    ec_numbers = []
    if pdb_library_dir.exists():
        for ec_dir in sorted(pdb_library_dir.glob("*/")):
            if ec_dir.is_dir() and not ec_dir.name.startswith('.'):
                ec_number = ec_dir.name.replace("_", ".")
                ec_numbers.append(ec_number)
    return sorted(ec_numbers)


def get_ec_numbers_from_motif_library(motif_library_dir: Path) -> List[str]:
    """从motif_library获取所有EC号"""
    ec_numbers = []
    if motif_library_dir.exists():
        for ec_dir in sorted(motif_library_dir.glob("*/")):
            if ec_dir.is_dir() and not ec_dir.name.startswith('.'):
                ec_number = ec_dir.name.replace("_", ".")
                ec_numbers.append(ec_number)
    return sorted(ec_numbers)


def count_motifs_in_ec(ec_dir: Path) -> Dict[str, int]:
    """统计EC目录下的motif数量"""
    counts = {
        'metal_sites': 0,
        'catalytic_sites': 0,
        'binding_sites': 0,
        'other': 0,
        'total': 0
    }
    
    if not ec_dir.exists():
        return counts
    
    for category in ['metal_sites', 'catalytic_sites', 'binding_sites', 'other']:
        cat_dir = ec_dir / category
        if cat_dir.exists():
            count = len(list(cat_dir.glob("*.json")))
            counts[category] = count
            counts['total'] += count
    
    return counts


def format_time(seconds: float) -> str:
    """格式化时间"""
    if seconds < 60:
        return f"{seconds:.1f}秒"
    elif seconds < 3600:
        return f"{int(seconds // 60)}分{int(seconds % 60)}秒"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}小时{minutes}分"


def print_status_header():
    """打印状态表头"""
    print("\n" + "="*100)
    print(f"{'EC号':<15} {'状态':<12} {'进度':<20} {'Metal':<8} {'Catalytic':<10} {'Binding':<8} {'Other':<8} {'总计':<8}")
    print("="*100)


def print_status_line(ec_number: str, status: str, progress: str, counts: Dict[str, int]):
    """打印状态行"""
    ec_display = ec_number.replace(".", "_")
    print(f"{ec_display:<15} {status:<12} {progress:<20} "
          f"{counts['metal_sites']:<8} {counts['catalytic_sites']:<10} "
          f"{counts['binding_sites']:<8} {counts['other']:<8} {counts['total']:<8}")


def process_ec_with_monitoring(
    ec_number: str,
    project_root: Path,
    motif_library_dir: Path
) -> Tuple[bool, str, Dict[str, int]]:
    """
    处理单个EC号并返回结果
    
    Returns:
        (success, message, counts)
    """
    ec_dir_name = ec_number.replace(".", "_")
    ec_dir = motif_library_dir / ec_dir_name
    
    # 运行提取脚本
    script_path = project_root / "scripts" / "extract_motifs_from_cache.py"
    
    try:
        # 使用subprocess运行，捕获输出
        process = subprocess.Popen(
            [sys.executable, str(script_path), "--ec", ec_number, "--skip-index"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # 实时输出（可选，这里我们只等待完成）
        output_lines = []
        for line in process.stdout:
            output_lines.append(line)
            # 可以在这里实时打印，但为了简洁，我们只收集
        
        process.wait()
        
        # 统计结果
        counts = count_motifs_in_ec(ec_dir)
        
        if process.returncode == 0:
            return True, "完成", counts
        else:
            return False, f"错误(退出码{process.returncode})", counts
            
    except Exception as e:
        counts = count_motifs_in_ec(ec_dir)
        return False, f"异常: {str(e)[:30]}", counts


def main():
    project_root = Path(__file__).parent.parent
    pdb_library_dir = project_root / "pdb_library"
    motif_library_dir = project_root / "motif_library"
    
    # 获取所有EC号
    pdb_ec_numbers = get_ec_numbers_from_pdb_library(pdb_library_dir)
    motif_ec_numbers = get_ec_numbers_from_motif_library(motif_library_dir)
    
    print("\n" + "="*100)
    print("Motif提取监控系统")
    print("="*100)
    print(f"PDB库EC号数量: {len(pdb_ec_numbers)}")
    print(f"Motif库EC号数量: {len(motif_ec_numbers)}")
    print(f"缺失EC号数量: {len(pdb_ec_numbers) - len(motif_ec_numbers)}")
    print("="*100)
    
    # 找出缺失的EC号
    missing_ec = [ec for ec in pdb_ec_numbers if ec not in motif_ec_numbers]
    
    if not missing_ec:
        print("\n✓ 所有EC号都已处理完成！")
        print("\n当前状态:")
        print_status_header()
        for ec_number in pdb_ec_numbers:
            ec_dir_name = ec_number.replace(".", "_")
            ec_dir = motif_library_dir / ec_dir_name
            counts = count_motifs_in_ec(ec_dir)
            print_status_line(ec_number, "已完成", "100%", counts)
        return
    
    print(f"\n需要处理的EC号: {len(missing_ec)}")
    print(f"EC号列表: {[ec.replace('.', '_') for ec in missing_ec]}")
    
    # 显示初始状态
    print("\n初始状态:")
    print_status_header()
    
    all_ec = sorted(set(pdb_ec_numbers))
    for ec_number in all_ec:
        ec_dir_name = ec_number.replace(".", "_")
        ec_dir = motif_library_dir / ec_dir_name
        counts = count_motifs_in_ec(ec_dir)
        
        if ec_number in missing_ec:
            status = "待处理"
            progress = "0%"
        else:
            status = "已完成"
            progress = "100%"
        
        print_status_line(ec_number, status, progress, counts)
    
    # 开始处理
    print("\n" + "="*100)
    print("开始批量处理...")
    print("="*100)
    
    start_time = time.time()
    results = {}
    
    for i, ec_number in enumerate(missing_ec, 1):
        ec_dir_name = ec_number.replace(".", "_")
        ec_display = ec_number.replace(".", "_")
        
        print(f"\n[{i}/{len(missing_ec)}] 处理 EC {ec_number} ({ec_display})...")
        print("-" * 100)
        
        ec_start_time = time.time()
        success, message, counts = process_ec_with_monitoring(
            ec_number, project_root, motif_library_dir
        )
        ec_elapsed = time.time() - ec_start_time
        
        results[ec_number] = {
            'success': success,
            'message': message,
            'counts': counts,
            'time': ec_elapsed
        }
        
        status = "✓ 成功" if success else "✗ 失败"
        print(f"\n{status}: {message}")
        print(f"耗时: {format_time(ec_elapsed)}")
        print(f"提取结果: Metal={counts['metal_sites']}, Catalytic={counts['catalytic_sites']}, "
              f"Binding={counts['binding_sites']}, Other={counts['other']}, 总计={counts['total']}")
        
        # 更新状态显示
        print("\n当前进度:")
        print_status_header()
        for ec in all_ec:
            if ec in results:
                r = results[ec]
                status = "✓ 完成" if r['success'] else f"✗ {r['message']}"
                progress = "100%"
            elif ec in missing_ec:
                status = "处理中..."
                progress = "..."
            else:
                status = "已完成"
                progress = "100%"
            
            ec_dir = motif_library_dir / ec.replace(".", "_")
            counts = count_motifs_in_ec(ec_dir)
            print_status_line(ec, status, progress, counts)
    
    # 最终统计
    total_time = time.time() - start_time
    success_count = sum(1 for r in results.values() if r['success'])
    fail_count = len(results) - success_count
    total_motifs = sum(r['counts']['total'] for r in results.values())
    
    print("\n" + "="*100)
    print("处理完成！")
    print("="*100)
    print(f"总耗时: {format_time(total_time)}")
    print(f"成功: {success_count}/{len(missing_ec)}")
    print(f"失败: {fail_count}/{len(missing_ec)}")
    print(f"新增Motif总数: {total_motifs}")
    
    # 最终状态
    print("\n最终状态:")
    print_status_header()
    for ec_number in all_ec:
        ec_dir = motif_library_dir / ec_number.replace(".", "_")
        counts = count_motifs_in_ec(ec_dir)
        status = "已完成"
        progress = "100%"
        print_status_line(ec_number, status, progress, counts)
    
    # 构建索引
    print("\n" + "="*100)
    print("构建数据库索引...")
    print("="*100)
    
    from enzyme_viewer.motif_db import build_motif_index
    db_path = project_root / "enzyme_viewer" / "motif_index.db"
    build_motif_index(
        motif_library_dir=motif_library_dir,
        db_path=db_path,
        clear_existing=True
    )
    
    print("\n✓ 数据库索引构建完成！")


if __name__ == "__main__":
    main()

