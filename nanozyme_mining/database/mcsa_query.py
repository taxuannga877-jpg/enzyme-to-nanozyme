"""
M-CSA 数据库查询工具

提供快速查询 M-CSA 数据库中的金属位点、催化残基等信息
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict


class MCSAQuery:
    """M-CSA 数据库查询器"""
    
    def __init__(self, mcsa_file: str = "data/mcsa_database/mcsa_processed.json"):
        """
        初始化查询器
        
        Args:
            mcsa_file: M-CSA 处理后的 JSON 文件路径
        """
        self.mcsa_file = Path(mcsa_file)
        self.ec_to_entries = {}
        self.mcsa_id_to_entry = {}
        self.loaded = False
        
    def load(self):
        """加载 M-CSA 数据库"""
        if self.loaded:
            return
        
        if not self.mcsa_file.exists():
            print(f"⚠️  M-CSA 文件不存在: {self.mcsa_file}")
            print("   请先运行: python nanozyme_mining/database/download_mcsa.py")
            return
        
        with open(self.mcsa_file, 'r') as f:
            entries = json.load(f)
        
        # 按 EC 号索引
        self.ec_to_entries = defaultdict(list)
        for entry in entries:
            ec_numbers = entry.get("ec_numbers", [])
            for ec in ec_numbers:
                self.ec_to_entries[ec].append(entry)
            
            # 按 M-CSA ID 索引
            mcsa_id = entry.get("mcsa_id")
            if mcsa_id:
                self.mcsa_id_to_entry[mcsa_id] = entry
        
        self.loaded = True
        print(f"✓ 加载 M-CSA 数据库: {len(entries)} 条目, {len(self.ec_to_entries)} 个 EC 号")
    
    def query_by_ec(self, ec_number: str) -> List[Dict]:
        """
        根据 EC 号查询
        
        Args:
            ec_number: EC 号 (e.g., "1.15.1.1")
            
        Returns:
            M-CSA 条目列表
        """
        if not self.loaded:
            self.load()
        
        return self.ec_to_entries.get(ec_number, [])
    
    def get_metal_sites(self, ec_number: str) -> Dict:
        """
        获取 EC 号对应的金属位点信息
        
        Args:
            ec_number: EC 号
            
        Returns:
            金属位点信息字典
        """
        entries = self.query_by_ec(ec_number)
        
        if not entries:
            return {
                "has_metal": False,
                "metal_types": [],
                "metal_coordination": [],
                "catalytic_residues": [],
                "mcsa_references": [],
            }
        
        # 整合所有条目的金属信息
        metal_types_set = set()
        metal_coordination = []
        catalytic_residues = []
        mcsa_references = []
        has_metal = False
        
        for entry in entries:
            metal_info = entry.get("metal_info", {})
            
            if metal_info.get("has_metal", False):
                has_metal = True
                
                # 金属类型
                entry_metal_types = metal_info.get("metal_types", [])
                for metal in entry_metal_types:
                    metal_types_set.add(metal)
                
                # 金属配体 - 添加该条目对应的金属类型信息
                for ligand in metal_info.get("metal_ligands", []):
                    metal_coordination.append({
                        "residue_type": ligand.get("residue_type"),
                        "residue_number": ligand.get("residue_number"),
                        "chain": ligand.get("chain"),
                        "pdb_id": ligand.get("pdb_id"),
                        "roles": ligand.get("roles", []),
                        "mcsa_id": entry.get("mcsa_id"),
                        "metal_types": entry_metal_types,  # 添加该残基配位的金属类型
                    })
            
            # 催化残基 - 如果是金属配体，也添加金属类型信息
            entry_metal_types = metal_info.get("metal_types", []) if metal_info.get("has_metal", False) else []
            for residue in entry.get("catalytic_residues", []):
                is_metal_ligand = residue.get("is_metal_ligand", False)
                catalytic_residues.append({
                    "residue_type": residue.get("residue_type"),
                    "residue_number": residue.get("residue_number"),
                    "chain": residue.get("chain"),
                    "pdb_id": residue.get("pdb_id"),
                    "roles": residue.get("roles", []),
                    "roles_summary": residue.get("roles_summary", ""),
                    "is_metal_ligand": is_metal_ligand,
                    "mcsa_id": entry.get("mcsa_id"),
                    "metal_types": entry_metal_types if is_metal_ligand else [],  # 如果是金属配体，添加金属类型
                })
            
            # 参考信息
            mcsa_references.append({
                "mcsa_id": entry.get("mcsa_id"),
                "enzyme_name": entry.get("enzyme_name"),
                "description": entry.get("description"),
                "url": entry.get("url"),
                "uniprot_id": entry.get("uniprot_id"),
            })
        
        return {
            "has_metal": has_metal,
            "metal_types": sorted(list(metal_types_set)),
            "metal_coordination": metal_coordination,
            "catalytic_residues": catalytic_residues,
            "mcsa_references": mcsa_references,
        }
    
    def get_catalytic_residues(self, ec_number: str) -> List[Dict]:
        """
        获取催化残基列表
        
        Args:
            ec_number: EC 号
            
        Returns:
            催化残基列表
        """
        metal_sites = self.get_metal_sites(ec_number)
        return metal_sites["catalytic_residues"]
    
    def has_metal(self, ec_number: str) -> bool:
        """
        检查是否含有金属
        
        Args:
            ec_number: EC 号
            
        Returns:
            是否含有金属
        """
        metal_sites = self.get_metal_sites(ec_number)
        return metal_sites["has_metal"]
    
    def get_metal_types(self, ec_number: str) -> List[str]:
        """
        获取金属类型列表
        
        Args:
            ec_number: EC 号
            
        Returns:
            金属类型列表
        """
        metal_sites = self.get_metal_sites(ec_number)
        return metal_sites["metal_types"]
    
    def summary(self) -> Dict:
        """
        获取数据库摘要统计
        
        Returns:
            统计信息字典
        """
        if not self.loaded:
            self.load()
        
        total_entries = len(self.mcsa_id_to_entry)
        entries_with_metal = sum(
            1 for entry in self.mcsa_id_to_entry.values()
            if entry.get("metal_info", {}).get("has_metal", False)
        )
        
        return {
            "total_entries": total_entries,
            "total_ec_numbers": len(self.ec_to_entries),
            "entries_with_metal": entries_with_metal,
            "metal_percentage": f"{entries_with_metal/total_entries*100:.1f}%" if total_entries > 0 else "0%",
        }


# 全局单例
_mcsa_query = None


def get_mcsa_query() -> MCSAQuery:
    """获取全局 M-CSA 查询器单例"""
    global _mcsa_query
    if _mcsa_query is None:
        _mcsa_query = MCSAQuery()
        _mcsa_query.load()
    return _mcsa_query


def query_metal_sites(ec_number: str) -> Dict:
    """
    快速查询金属位点（便捷函数）
    
    Args:
        ec_number: EC 号
        
    Returns:
        金属位点信息
    """
    return get_mcsa_query().get_metal_sites(ec_number)


if __name__ == "__main__":
    # 测试查询
    query = MCSAQuery()
    query.load()
    
    # 示例：查询超氧化物歧化酶 (SOD)
    test_ecs = ["1.15.1.1", "3.4.21.4", "1.11.1.6"]
    
    print("\n" + "=" * 80)
    print("M-CSA 查询测试")
    print("=" * 80)
    
    for ec in test_ecs:
        print(f"\n🔍 查询 EC {ec}:")
        metal_sites = query.get_metal_sites(ec)
        
        print(f"   含金属: {metal_sites['has_metal']}")
        print(f"   金属类型: {', '.join(metal_sites['metal_types']) if metal_sites['metal_types'] else '无'}")
        print(f"   催化残基数: {len(metal_sites['catalytic_residues'])}")
        print(f"   M-CSA 参考: {len(metal_sites['mcsa_references'])}")
        
        if metal_sites['mcsa_references']:
            for ref in metal_sites['mcsa_references'][:2]:
                print(f"      - {ref['enzyme_name']} (MCSA-{ref['mcsa_id']})")
    
    print("\n" + "=" * 80)
    print("数据库摘要:")
    summary = query.summary()
    for key, value in summary.items():
        print(f"   {key}: {value}")
    print("=" * 80)


