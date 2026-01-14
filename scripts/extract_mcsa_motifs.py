#!/usr/bin/env python3
"""
M-CSA Motif提取脚本

从M-CSA数据库中提取motif，为每个EC号提取对应的催化残基信息并生成motif。
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# 添加项目根目录到Python路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from nanozyme_mining.database.mcsa_query import MCSAQuery
from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher
from nanozyme_mining.extraction.extractor import MotifExtractor
from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser


def get_ec_numbers_from_motif_library(motif_library_dir: Path) -> List[str]:
    """从motif_library目录获取所有EC号"""
    ec_numbers = set()
    
    for sub_dir in motif_library_dir.iterdir():
        if not sub_dir.is_dir():
            continue
        
        # 检查是否为EC号格式目录（如 1_11_1_6）
        dir_name = sub_dir.name
        parts = dir_name.split('_')
        if len(parts) == 4:
            try:
                [int(p) for p in parts]
                ec_number = dir_name.replace('_', '.')
                ec_numbers.add(ec_number)
            except ValueError:
                pass
    
    return sorted(list(ec_numbers))


def get_nanozyme_name(ec_number: str) -> str:
    """根据EC号获取nanozyme名称"""
    # 这里可以根据EC号映射到nanozyme名称
    # 暂时返回EC号作为名称
    return ec_number


def convert_mcsa_residue_to_active_site_indices(
    pdb_path: Path,
    catalytic_residues: List[Dict]
) -> List[int]:
    """
    将M-CSA的催化残基信息转换为active_site_indices
    
    Args:
        pdb_path: PDB文件路径
        catalytic_residues: M-CSA的催化残基列表，每个包含pdb_id, chain, residue_number等
    
    Returns:
        残基索引列表（UniProt序列位置）
    """
    active_site_indices = []
    
    try:
        # 解析PDB文件，获取残基映射信息
        pdb_parser = ComprehensivePDBParser()
        parsed_data = pdb_parser.parse_pdb_file(pdb_path)
        
        # 从PDB文件中提取ATOM记录，构建残基编号到序列位置的映射
        # 注意：这里需要根据实际情况处理PDB编号到UniProt编号的映射
        # 简化处理：直接使用PDB残基编号（如果PDB文件中的编号就是序列位置）
        
        for residue in catalytic_residues:
            pdb_residue_number = residue.get("residue_number")
            chain = residue.get("chain", "")
            
            if pdb_residue_number is not None:
                # 简化处理：假设PDB残基编号就是序列位置
                # 实际应用中可能需要更复杂的映射逻辑
                active_site_indices.append(pdb_residue_number)
    
    except Exception as e:
        print(f"  ⚠️  转换残基编号时出错: {e}")
    
    return active_site_indices


def extract_mcsa_motif_for_ec(
    ec_number: str,
    mcsa_query: MCSAQuery,
    uniprot_fetcher: UniProtFetcher,
    motif_extractor: MotifExtractor,
    pdb_library_dir: Path,
    motif_library_dir: Path
) -> int:
    """
    为单个EC号提取M-CSA motif
    
    Returns:
        成功提取的motif数量
    """
    print(f"\n{'='*80}")
    print(f"处理 EC {ec_number}")
    print(f"{'='*80}")
    
    # 查询M-CSA数据库
    mcsa_entries = mcsa_query.query_by_ec(ec_number)
    
    if not mcsa_entries:
        print(f"  ⚠️  未找到M-CSA条目")
        return 0
    
    print(f"  ✓ 找到 {len(mcsa_entries)} 个M-CSA条目")
    
    extracted_count = 0
    
    for entry in mcsa_entries:
        mcsa_id = entry.get("mcsa_id")
        uniprot_id = entry.get("uniprot_id", "")
        catalytic_residues = entry.get("catalytic_residues", [])
        
        if not catalytic_residues:
            print(f"  ⚠️  M-CSA条目 {mcsa_id} 没有催化残基信息")
            continue
        
        print(f"\n  处理 M-CSA条目 {mcsa_id} (UniProt: {uniprot_id})")
        print(f"    催化残基数: {len(catalytic_residues)}")
        
        # 从催化残基中提取PDB ID（取第一个残基的PDB ID）
        pdb_id = None
        for residue in catalytic_residues:
            pdb_id = residue.get("pdb_id", "").lower()
            if pdb_id:
                break
        
        if not pdb_id:
            print(f"    ⚠️  未找到PDB ID")
            continue
        
        # 查找PDB文件
        ec_dir_name = ec_number.replace(".", "_")
        pdb_library_ec_dir = pdb_library_dir / ec_dir_name
        
        pdb_path = None
        # 优先查找实验PDB
        exp_pdb_path = pdb_library_ec_dir / f"{pdb_id.upper()}.pdb"
        if exp_pdb_path.exists():
            pdb_path = exp_pdb_path
        else:
            # 如果不存在，尝试下载
            print(f"    📥 下载PDB文件 {pdb_id}...")
            try:
                downloaded_path = uniprot_fetcher.download_pdb(pdb_id)
                if downloaded_path and downloaded_path.exists():
                    # 移动到pdb_library目录
                    target_path = pdb_library_ec_dir / downloaded_path.name
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    import shutil
                    shutil.copy2(downloaded_path, target_path)
                    pdb_path = target_path
            except Exception as e:
                print(f"    ⚠️  下载PDB文件失败: {e}")
        
        if not pdb_path or not pdb_path.exists():
            print(f"    ⚠️  PDB文件不存在: {pdb_id}")
            continue
        
        print(f"    ✓ 找到PDB文件: {pdb_path}")
        
        # 转换M-CSA残基信息为active_site_indices
        active_site_indices = convert_mcsa_residue_to_active_site_indices(
            pdb_path, catalytic_residues
        )
        
        if not active_site_indices:
            print(f"    ⚠️  无法转换残基编号")
            continue
        
        print(f"    ✓ 提取到 {len(active_site_indices)} 个活性位点残基")
        
        # 构建functional_roles字典
        functional_roles = {}
        for residue in catalytic_residues:
            res_name = residue.get("residue_type", "")
            res_num = residue.get("residue_number")
            roles = residue.get("roles", [])
            roles_summary = residue.get("roles_summary", "")
            
            if res_num is not None:
                key = (res_name, res_num)
                functional_roles[key] = roles_summary if roles_summary else ", ".join(roles)
        
        # 获取nanozyme类型（从EC号推断）
        nanozyme_name = get_nanozyme_name(ec_number)
        nanozyme_type = "POD"  # 默认值，可以根据EC号映射
        
        # 提取motif
        try:
            print(f"    🔧 提取motif...")
            motif = motif_extractor.extract_motif(
                pdb_path=str(pdb_path),
                uniprot_id=uniprot_id,
                ec_number=ec_number,
                nanozyme_type=nanozyme_type,
                active_site_indices=active_site_indices if active_site_indices else None,
                functional_roles=functional_roles if functional_roles else None
            )
            
            if motif is None:
                print(f"    ⚠️  未能提取到motif")
                continue
            
            # 生成motif_id（格式：{uniprot_id}_{ec_number}_{nanozyme_name}_mcsa_{mcsa_id}）
            motif_id = f"{uniprot_id}_{ec_number.replace('.', '_')}_{nanozyme_name}_mcsa_{mcsa_id}"
            motif.motif_id = motif_id
            
            # 保存motif到motif_library目录
            ec_dir_name = ec_number.replace(".", "_")
            motif_library_ec_dir = motif_library_dir / ec_dir_name
            
            # 分类motif
            from enzyme_viewer.motif_db import classify_motif
            motif_dict = motif.to_dict()
            motif_dict['source'] = 'M-CSA'
            motif_dict['mcsa_id'] = mcsa_id
            motif_dict['extraction_method'] = 'mcsa'
            
            category = classify_motif(motif_dict)
            category_dir = motif_library_ec_dir / category
            category_dir.mkdir(parents=True, exist_ok=True)
            
            motif_file = category_dir / f"{motif_id}.json"
            with open(motif_file, 'w', encoding='utf-8') as f:
                json.dump(motif_dict, f, indent=2, ensure_ascii=False)
            
            print(f"    ✅ 成功提取并保存motif: {motif_file}")
            extracted_count += 1
        
        except Exception as e:
            print(f"    ❌ 提取motif失败: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    return extracted_count


def main():
    """主函数"""
    print("="*80)
    print("M-CSA Motif提取脚本")
    print("="*80)
    
    # 配置路径
    project_root = Path(__file__).parent.parent
    pdb_library_dir = project_root / 'pdb_library'
    motif_library_dir = project_root / 'motif_library'
    mcsa_file = project_root / 'data' / 'mcsa_database' / 'mcsa_processed.json'
    
    # 初始化组件
    print("\n初始化组件...")
    mcsa_query = MCSAQuery(mcsa_file=str(mcsa_file))
    mcsa_query.load()
    
    uniprot_fetcher = UniProtFetcher(
        cache_dir=str(project_root / 'cache'),
        pdb_library_dir=str(pdb_library_dir)
    )
    
    motif_extractor = MotifExtractor(output_dir=str(project_root / 'motifs'))
    
    print("✓ 组件初始化完成\n")
    
    # 获取所有EC号
    print("扫描motif_library目录获取EC号...")
    ec_numbers = get_ec_numbers_from_motif_library(motif_library_dir)
    print(f"✓ 找到 {len(ec_numbers)} 个EC号\n")
    
    if not ec_numbers:
        print("⚠️  未找到EC号，请先运行motif提取脚本")
        return
    
    # 为每个EC号提取M-CSA motif
    total_extracted = 0
    for ec_number in ec_numbers:
        try:
            count = extract_mcsa_motif_for_ec(
                ec_number,
                mcsa_query,
                uniprot_fetcher,
                motif_extractor,
                pdb_library_dir,
                motif_library_dir
            )
            total_extracted += count
        except Exception as e:
            print(f"\n❌ 处理EC {ec_number} 时出错: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "="*80)
    print(f"提取完成！共提取 {total_extracted} 个M-CSA motif")
    print("="*80)


if __name__ == "__main__":
    main()


