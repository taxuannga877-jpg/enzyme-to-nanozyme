#!/usr/bin/env python
"""
纳米酶PDB数据库建立与Motif提取完整流程
========================================

流程：
1. 根据EC号从UniProt获取酶数据
2. 下载AlphaFold PDB结构
3. 分类为有标注/无标注数据
4. 对无标注数据使用EasIFA预测活性位点
5. 提取催化motif并保存到motif_library

使用方法:
    python scripts/run_pipeline.py --ec 1.11.1.7 --max_results 100
    python scripts/run_pipeline.py --all --max_results 50
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from nanozyme_mining.database import UniProtFetcher, NanozymeDatabase
from nanozyme_mining.core import DualTrackProcessor
from nanozyme_mining.extraction import MotifExtractor
from nanozyme_mining.utils.constants import NanozymeType, EC_TO_NANOZYME_TYPE
from nanozyme_mining.utils.ec_mappings import EC_PATTERNS


def get_all_ec_numbers() -> List[str]:
    """获取所有支持的EC号列表"""
    all_ecs = []
    for nanozyme_type, ec_list in EC_PATTERNS.items():
        all_ecs.extend(ec_list)
    return sorted(set(all_ecs))


def process_ec_number(
    ec_number: str,
    max_results: int = 100,
    cache_dir: str = "./cache",
    processed_dir: str = "./processed",
    motif_dir: str = "./motif_library",
    device: str = "cpu",
    skip_prediction: bool = False
):
    """
    处理单个EC号的完整流程
    
    Args:
        ec_number: EC号
        max_results: 最大获取数量
        cache_dir: 缓存目录
        processed_dir: 处理结果目录
        motif_dir: motif输出目录
        device: EasIFA运行设备
        skip_prediction: 是否跳过预测（仅处理有标注数据）
    """
    print("\n" + "="*80)
    print(f"处理 EC {ec_number}")
    print("="*80)
    
    # 获取纳米酶类型
    nanozyme_type = EC_TO_NANOZYME_TYPE.get(ec_number)
    if not nanozyme_type:
        print(f"  ⚠️  未找到EC号 {ec_number} 对应的纳米酶类型，跳过")
        return
    
    print(f"  纳米酶类型: {nanozyme_type.value}")
    
    # Step 1: 获取并分类数据
    print(f"\n[Step 1] 从UniProt获取数据 (max_results={max_results})...")
    fetcher = UniProtFetcher(cache_dir=cache_dir)
    
    try:
        annotated, unannotated = fetcher.fetch_and_classify(
            ec_number=ec_number,
            nanozyme_type=nanozyme_type,
            max_results=max_results
        )
        print(f"  ✓ 获取完成: {len(annotated)} 有标注, {len(unannotated)} 无标注")
    except Exception as e:
        print(f"  ✗ 获取数据失败: {e}")
        return
    
    if len(annotated) == 0 and len(unannotated) == 0:
        print(f"  ⚠️  没有获取到任何数据，跳过")
        return
    
    # Step 2: 处理有标注数据（直接使用）
    annotated_processed = []
    if len(annotated) > 0:
        print(f"\n[Step 2] 处理 {len(annotated)} 条有标注数据...")
        for entry in annotated:
            if entry.get("pdb_path") and os.path.exists(entry["pdb_path"]):
                annotated_processed.append(entry)
            else:
                print(f"  ⚠️  {entry.get('uniprot_id')}: PDB文件不存在")
        print(f"  ✓ 有效有标注数据: {len(annotated_processed)}")
    
    # Step 3: 处理无标注数据（使用EasIFA预测）
    unannotated_processed = []
    if len(unannotated) > 0 and not skip_prediction:
        print(f"\n[Step 3] 使用EasIFA预测 {len(unannotated)} 条无标注数据...")
        try:
            processor = DualTrackProcessor(
                output_dir=processed_dir,
                device=device
            )
            
            # 只处理有PDB文件的条目
            valid_unannotated = [
                e for e in unannotated 
                if e.get("pdb_path") and os.path.exists(e["pdb_path"])
            ]
            
            if valid_unannotated:
                predicted_results = processor.predict_unannotated_batch(
                    valid_unannotated,
                    reaction_smiles="C>>C"  # 默认反应，可以根据EC号定制
                )
                unannotated_processed = [
                    {
                        "uniprot_id": r.uniprot_id,
                        "ec_number": r.ec_number,
                        "nanozyme_type": r.nanozyme_type,
                        "pdb_path": r.pdb_path,
                        "sequence": r.sequence,
                        "active_sites": r.active_sites,
                        "source": r.source
                    }
                    for r in predicted_results
                ]
                print(f"  ✓ 预测完成: {len(unannotated_processed)} 条")
            else:
                print(f"  ⚠️  没有有效的PDB文件用于预测")
        except Exception as e:
            print(f"  ✗ 预测失败: {e}")
            print(f"    提示: 可能需要配置EasIFA模型路径")
    elif len(unannotated) > 0 and skip_prediction:
        print(f"\n[Step 3] 跳过预测，仅处理有标注数据")
    
    # Step 4: 提取Motif
    print(f"\n[Step 4] 提取催化Motif...")
    
    # 创建motif输出目录（按纳米酶类型组织）
    type_motif_dir = Path(motif_dir) / nanozyme_type.value
    type_motif_dir.mkdir(parents=True, exist_ok=True)
    
    extractor = MotifExtractor(output_dir=str(type_motif_dir))
    
    all_entries = annotated_processed + unannotated_processed
    success_count = 0
    fail_count = 0
    
    for entry in all_entries:
        uniprot_id = entry.get("uniprot_id", "")
        pdb_path = entry.get("pdb_path", "")
        
        if not pdb_path or not os.path.exists(pdb_path):
            print(f"  ⚠️  {uniprot_id}: PDB文件不存在")
            fail_count += 1
            continue
        
        # 提取活性位点残基索引
        active_sites = entry.get("active_sites", [])
        site_indices = []
        for site in active_sites:
            # 根据不同的数据源格式提取残基索引
            if "start" in site:
                site_indices.append(site["start"])
            elif "residue_index" in site:
                site_indices.append(site["residue_index"])
            elif "residue_number" in site:
                site_indices.append(site["residue_number"])
        
        try:
            motif = extractor.extract_motif(
                pdb_path=pdb_path,
                uniprot_id=uniprot_id,
                ec_number=ec_number,
                nanozyme_type=nanozyme_type.value,
                active_site_indices=site_indices if site_indices else None
            )
            
            if motif:
                output_file = type_motif_dir / f"{motif.motif_id}.json"
                motif.to_json(str(output_file))
                success_count += 1
                print(f"  ✓ {uniprot_id}: 提取成功")
            else:
                fail_count += 1
                print(f"  ✗ {uniprot_id}: 提取失败（未找到催化残基）")
                
        except Exception as e:
            fail_count += 1
            print(f"  ✗ {uniprot_id}: 提取失败 - {e}")
    
    print(f"\n[完成] EC {ec_number} 处理完成:")
    print(f"  - 有标注数据: {len(annotated_processed)}")
    print(f"  - 预测数据: {len(unannotated_processed)}")
    print(f"  - Motif提取成功: {success_count}")
    print(f"  - Motif提取失败: {fail_count}")


def main():
    parser = argparse.ArgumentParser(
        description="纳米酶PDB数据库建立与Motif提取完整流程"
    )
    parser.add_argument(
        "--ec",
        type=str,
        help="要处理的EC号（如 1.11.1.7），可以指定多个，用逗号分隔"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="处理所有支持的EC号"
    )
    parser.add_argument(
        "--max_results",
        type=int,
        default=100,
        help="每个EC号最大获取数量（默认: 100）"
    )
    parser.add_argument(
        "--cache_dir",
        type=str,
        default="./cache",
        help="缓存目录（默认: ./cache）"
    )
    parser.add_argument(
        "--processed_dir",
        type=str,
        default="./processed",
        help="处理结果目录（默认: ./processed）"
    )
    parser.add_argument(
        "--motif_dir",
        type=str,
        default="./motif_library",
        help="Motif输出目录（默认: ./motif_library）"
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="EasIFA运行设备（默认: cpu）"
    )
    parser.add_argument(
        "--skip_prediction",
        action="store_true",
        help="跳过预测步骤，仅处理有标注数据"
    )
    
    args = parser.parse_args()
    
    # 确定要处理的EC号列表
    if args.all:
        ec_numbers = get_all_ec_numbers()
        print(f"将处理所有 {len(ec_numbers)} 个EC号")
    elif args.ec:
        ec_numbers = [ec.strip() for ec in args.ec.split(",")]
    else:
        print("错误: 必须指定 --ec 或 --all")
        parser.print_help()
        return
    
    print("\n" + "="*80)
    print("纳米酶PDB数据库建立与Motif提取流程")
    print("="*80)
    print(f"EC号列表: {ec_numbers}")
    print(f"每个EC号最大数量: {args.max_results}")
    print(f"设备: {args.device}")
    print(f"跳过预测: {args.skip_prediction}")
    print("="*80)
    
    # 处理每个EC号
    total_success = 0
    total_fail = 0
    
    for i, ec_number in enumerate(ec_numbers, 1):
        print(f"\n[{i}/{len(ec_numbers)}] 处理 EC {ec_number}...")
        try:
            process_ec_number(
                ec_number=ec_number,
                max_results=args.max_results,
                cache_dir=args.cache_dir,
                processed_dir=args.processed_dir,
                motif_dir=args.motif_dir,
                device=args.device,
                skip_prediction=args.skip_prediction
            )
        except KeyboardInterrupt:
            print("\n\n用户中断，退出...")
            break
        except Exception as e:
            print(f"\n  ✗ 处理 EC {ec_number} 时出错: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*80)
    print("所有EC号处理完成！")
    print("="*80)
    print(f"Motif库位置: {args.motif_dir}")
    print(f"处理结果位置: {args.processed_dir}")


if __name__ == "__main__":
    main()

