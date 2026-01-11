#!/usr/bin/env python
"""
批量提取催化motif脚本
按照纳米酶类型组织motif库
"""

import json
import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from nanozyme_mining import MotifExtractor


def extract_motifs_by_type(nanozyme_type: str, annotated_dir: Path, output_dir: Path):
    """
    按照纳米酶类型批量提取motif

    Args:
        nanozyme_type: 纳米酶类型（如 "Peroxidase", "Catalase"）
        annotated_dir: 有标注数据目录
        output_dir: motif输出目录
    """
    # 创建类型专属目录
    type_output_dir = output_dir / nanozyme_type.replace(" ", "_")
    type_output_dir.mkdir(parents=True, exist_ok=True)

    # 初始化提取器
    extractor = MotifExtractor(output_dir=str(type_output_dir))

    # 查找该类型的所有数据文件
    json_files = list(annotated_dir.glob("*.json"))

    success_count = 0
    fail_count = 0

    for json_file in json_files:
        with open(json_file, 'r') as f:
            data = json.load(f)

        # 检查是否是目标类型
        if data.get('nanozyme_type') != nanozyme_type:
            continue

        # 提取活性位点索引
        site_indices = [
            site.get('start')
            for site in data.get('active_sites', [])
            if site.get('start')
        ]

        # 提取motif
        try:
            motif = extractor.extract_motif(
                pdb_path=data['pdb_path'],
                uniprot_id=data['uniprot_id'],
                ec_number=data['ec_number'],
                nanozyme_type=data['nanozyme_type'],
                active_site_indices=site_indices
            )

            if motif:
                # 保存motif
                output_file = type_output_dir / f"{motif.motif_id}.json"
                motif.to_json(str(output_file))
                success_count += 1
            else:
                fail_count += 1

        except Exception as e:
            print(f"  ✗ {data['uniprot_id']}: {e}")
            fail_count += 1

    return success_count, fail_count


if __name__ == "__main__":
    # 配置路径
    annotated_dir = Path("cache/annotated")
    output_dir = Path("motif_library")

    # 纳米酶类型列表
    nanozyme_types = [
        "Peroxidase",
        "Catalase",
        "Superoxide Dismutase",
        "Glutathione Peroxidase",
        "Laccase",
        "Oxidase"
    ]

    print("========================================")
    print("批量提取催化Motif - 按酶活性分类")
    print("========================================\n")

    total_success = 0
    total_fail = 0

    for i, ntype in enumerate(nanozyme_types, 1):
        print(f"[{i}/{len(nanozyme_types)}] {ntype}")
        success, fail = extract_motifs_by_type(ntype, annotated_dir, output_dir)
        print(f"  ✓ 成功: {success}  ✗ 失败: {fail}\n")

        total_success += success
        total_fail += fail

    print("========================================")
    print(f"总计: {total_success} 成功, {total_fail} 失败")
    print("========================================")
