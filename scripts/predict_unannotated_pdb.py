#!/usr/bin/env python
"""
使用 EasIFA 预测未标注的 PDB 文件
==================================

此脚本会：
1. 扫描 cache/unannotated/ 目录下的所有 JSON 文件
2. 使用 EasIFA 模型预测每个 PDB 文件的活性位点
3. 保存预测结果到 processed/predicted/ 目录

使用方法:
    python scripts/predict_unannotated_pdb.py
    python scripts/predict_unannotated_pdb.py --device cuda
    python scripts/predict_unannotated_pdb.py --cache_dir ./cache --output_dir ./processed
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from nanozyme_mining.core import DualTrackProcessor
from nanozyme_mining.prediction import EasIFAPredictor


def load_unannotated_entries(cache_dir: str = "./cache") -> List[Dict]:
    """
    加载所有未标注的条目
    
    Args:
        cache_dir: 缓存目录路径
        
    Returns:
        未标注条目列表
    """
    unannotated_dir = Path(cache_dir) / "unannotated"
    
    if not unannotated_dir.exists():
        print(f"[错误] 未标注数据目录不存在: {unannotated_dir}")
        return []
    
    entries = []
    json_files = list(unannotated_dir.glob("*.json"))
    
    print(f"[信息] 找到 {len(json_files)} 个未标注的 JSON 文件")
    
    for json_file in json_files:
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                entry = json.load(f)
                
            # 检查是否有 PDB 文件
            pdb_path = entry.get("pdb_path")
            if not pdb_path:
                print(f"  ⚠️  {entry.get('uniprot_id', 'Unknown')}: 没有 PDB 路径")
                continue
                
            # 检查 PDB 文件是否存在
            pdb_full_path = Path(cache_dir).parent / pdb_path if not Path(pdb_path).is_absolute() else Path(pdb_path)
            if not pdb_full_path.exists():
                print(f"  ⚠️  {entry.get('uniprot_id', 'Unknown')}: PDB 文件不存在: {pdb_path}")
                continue
            
            # 更新为绝对路径
            entry["pdb_path"] = str(pdb_full_path)
            entries.append(entry)
            
        except Exception as e:
            print(f"  ✗ 读取 {json_file} 失败: {e}")
            continue
    
    print(f"[信息] 有效条目: {len(entries)} 个（有 PDB 文件）")
    return entries


def predict_single_pdb(
    predictor: EasIFAPredictor,
    entry: Dict,
    reaction_smiles: str = "C>>C"
) -> Optional[Dict]:
    """
    预测单个 PDB 文件的活性位点
    
    Args:
        predictor: EasIFA 预测器
        entry: 条目数据
        reaction_smiles: 反应 SMILES
        
    Returns:
        预测结果字典，失败返回 None
    """
    uniprot_id = entry.get("uniprot_id", "Unknown")
    pdb_path = entry.get("pdb_path")
    
    if not pdb_path:
        print(f"  ✗ {uniprot_id}: 没有 PDB 路径")
        return None
    
    try:
        # 使用 EasIFA 预测
        result = predictor.predict_with_details(
            pdb_path=pdb_path,
            uniprot_id=uniprot_id,
            reaction_smiles=reaction_smiles
        )
        
        if result is None:
            print(f"  ✗ {uniprot_id}: 预测失败")
            return None
        
        # 转换为字典格式
        return {
            "uniprot_id": result.uniprot_id,
            "ec_number": entry.get("ec_number", ""),
            "nanozyme_type": entry.get("nanozyme_type", ""),
            "pdb_path": result.pdb_path,
            "sequence": entry.get("sequence", ""),
            "active_sites": [
                {
                    "residue_index": s.residue_index,
                    "residue_name": s.residue_name,
                    "site_type": s.site_type,
                    "coordinates": s.coordinates
                }
                for s in result.sites
            ],
            "source": "predicted",
            "num_sites": len(result.sites),
            "num_catalytic": sum(1 for s in result.sites if s.site_type == "Catalytic"),
            "num_binding": sum(1 for s in result.sites if s.site_type == "Binding")
        }
        
    except Exception as e:
        print(f"  ✗ {uniprot_id}: 预测异常 - {e}")
        return None


def save_prediction_result(result: Dict, output_dir: str):
    """
    保存预测结果到 JSON 文件
    
    Args:
        result: 预测结果字典
        output_dir: 输出目录
    """
    output_path = Path(output_dir) / "predicted" / f"{result['uniprot_id']}.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    
    print(f"  ✓ 已保存: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="使用 EasIFA 预测未标注的 PDB 文件",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 使用默认设置（CPU）
  python scripts/predict_unannotated_pdb.py
  
  # 使用 GPU
  python scripts/predict_unannotated_pdb.py --device cuda
  
  # 指定目录
  python scripts/predict_unannotated_pdb.py --cache_dir ./cache --output_dir ./processed
  
  # 指定反应 SMILES
  python scripts/predict_unannotated_pdb.py --reaction "C>>C"
        """
    )
    
    parser.add_argument(
        "--cache_dir",
        type=str,
        default="./cache",
        help="缓存目录路径（默认: ./cache）"
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./processed",
        help="输出目录路径（默认: ./processed）"
    )
    
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="运行设备（默认: cpu）"
    )
    
    parser.add_argument(
        "--reaction",
        type=str,
        default="C>>C",
        help="反应 SMILES（默认: C>>C）"
    )
    
    parser.add_argument(
        "--batch_size",
        type=int,
        default=10,
        help="批处理大小（默认: 10，设置为 0 表示不显示进度）"
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("EasIFA 未标注 PDB 文件批量预测")
    print("=" * 60)
    print(f"缓存目录: {args.cache_dir}")
    print(f"输出目录: {args.output_dir}")
    print(f"运行设备: {args.device}")
    print(f"反应 SMILES: {args.reaction}")
    print()
    
    # Step 1: 加载未标注条目
    print("[Step 1] 加载未标注条目...")
    entries = load_unannotated_entries(args.cache_dir)
    
    if len(entries) == 0:
        print("  ⚠️  没有找到未标注的条目，退出")
        return
    
    print(f"  ✓ 加载完成: {len(entries)} 个条目")
    print()
    
    # Step 2: 初始化 EasIFA 预测器
    print("[Step 2] 初始化 EasIFA 预测器...")
    try:
        predictor = EasIFAPredictor(device=args.device)
        print("  ✓ EasIFA 预测器初始化成功")
    except Exception as e:
        print(f"  ✗ EasIFA 预测器初始化失败: {e}")
        print("\n提示:")
        print("  1. 确保 EasIFA 模型文件已正确放置")
        print("  2. 确保 ChemEnzyRetroPlanner 已正确安装")
        print("  3. 检查环境变量 CHEMENZY_PATH 是否正确设置")
        return
    print()
    
    # Step 3: 批量预测
    print(f"[Step 3] 开始批量预测 {len(entries)} 个 PDB 文件...")
    print()
    
    success_count = 0
    fail_count = 0
    
    for i, entry in enumerate(entries, 1):
        uniprot_id = entry.get("uniprot_id", "Unknown")
        print(f"[{i}/{len(entries)}] 预测 {uniprot_id}...")
        
        result = predict_single_pdb(
            predictor=predictor,
            entry=entry,
            reaction_smiles=args.reaction
        )
        
        if result:
            save_prediction_result(result, args.output_dir)
            success_count += 1
            print(f"  ✓ 预测成功: {result['num_sites']} 个活性位点 "
                  f"({result['num_catalytic']} 催化, {result['num_binding']} 结合)")
        else:
            fail_count += 1
        
        print()
    
    # Step 4: 总结
    print("=" * 60)
    print("预测完成")
    print("=" * 60)
    print(f"总条目数: {len(entries)}")
    print(f"成功: {success_count}")
    print(f"失败: {fail_count}")
    print(f"成功率: {success_count/len(entries)*100:.1f}%")
    print()
    print(f"预测结果已保存到: {Path(args.output_dir) / 'predicted'}")


if __name__ == "__main__":
    main()




