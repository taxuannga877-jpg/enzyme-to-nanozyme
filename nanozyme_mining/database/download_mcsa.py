#!/usr/bin/env python3
"""
完整下载 M-CSA (Mechanism and Catalytic Site Atlas) 数据库

该脚本会：
1. 遍历所有 M-CSA 条目（通过分页 API）
2. 下载每个条目的详细信息
3. 提取催化残基、金属配位、机制描述
4. 保存为结构化的 JSON 数据库

M-CSA 包含约 1000+ 个酶的催化机制数据
"""

import os
import json
import time
import requests
from pathlib import Path
from typing import Dict, List, Optional
from tqdm import tqdm


class MCSADownloader:
    """M-CSA 数据库完整下载器"""
    
    def __init__(self, output_dir: str = "data/mcsa_database"):
        self.base_url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"
        self.entries_url = f"{self.base_url}/entries/"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储路径
        self.raw_dir = self.output_dir / "raw_entries"
        self.raw_dir.mkdir(exist_ok=True)
        
        self.processed_file = self.output_dir / "mcsa_processed.json"
        self.summary_file = self.output_dir / "mcsa_summary.json"
        
    def download_all_entries(self) -> List[Dict]:
        """
        下载所有 M-CSA 条目（分页获取）
        
        Returns:
            所有条目的基本信息列表
        """
        print("=" * 80)
        print("开始下载完整 M-CSA 数据库")
        print("=" * 80)
        
        all_entries = []
        next_url = self.entries_url
        page = 1
        
        while next_url:
            print(f"\n📥 正在下载第 {page} 页...")
            
            try:
                response = requests.get(next_url, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                # 获取结果
                results = data.get("results", [])
                all_entries.extend(results)
                
                # 显示进度
                count = data.get("count", 0)
                print(f"   已获取 {len(all_entries)}/{count} 个条目")
                
                # 获取下一页 URL
                next_url = data.get("next")
                page += 1
                
                # 礼貌等待，避免过载服务器
                time.sleep(0.5)
                
            except Exception as e:
                print(f"❌ 下载第 {page} 页失败: {e}")
                break
        
        print(f"\n✅ 成功下载 {len(all_entries)} 个条目的基本信息")
        
        # 保存基本信息
        basic_file = self.output_dir / "mcsa_basic_info.json"
        with open(basic_file, 'w') as f:
            json.dump(all_entries, f, indent=2)
        print(f"   已保存到: {basic_file}")
        
        return all_entries
    
    def download_detailed_entry(self, mcsa_id: int) -> Optional[Dict]:
        """
        下载单个条目的详细信息
        
        Args:
            mcsa_id: M-CSA 条目 ID
            
        Returns:
            详细信息字典，失败返回 None
        """
        url = f"{self.entries_url}{mcsa_id}/"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            print(f"   ⚠️  下载 MCSA-{mcsa_id} 详细信息失败: {e}")
            return None
    
    def download_all_detailed_entries(self, basic_entries: List[Dict]) -> List[Dict]:
        """
        下载所有条目的详细信息
        
        Args:
            basic_entries: 基本信息列表
            
        Returns:
            详细信息列表
        """
        print("\n" + "=" * 80)
        print(f"开始下载 {len(basic_entries)} 个条目的详细信息")
        print("=" * 80)
        
        detailed_entries = []
        
        for entry in tqdm(basic_entries, desc="下载详细信息"):
            mcsa_id = entry.get("mcsa_id")
            if not mcsa_id:
                continue
            
            # 检查是否已下载
            raw_file = self.raw_dir / f"{mcsa_id}.json"
            if raw_file.exists():
                with open(raw_file, 'r') as f:
                    detailed = json.load(f)
            else:
                detailed = self.download_detailed_entry(mcsa_id)
                if detailed:
                    # 保存原始数据
                    with open(raw_file, 'w') as f:
                        json.dump(detailed, f, indent=2)
                    
                    # 礼貌等待
                    time.sleep(0.3)
            
            if detailed:
                detailed_entries.append(detailed)
        
        print(f"\n✅ 成功下载 {len(detailed_entries)} 个详细条目")
        return detailed_entries
    
    def parse_catalytic_residues(self, residues: List[Dict]) -> List[Dict]:
        """
        解析催化残基信息
        
        Args:
            residues: 残基列表
            
        Returns:
            解析后的催化残基信息
        """
        parsed_residues = []
        
        for residue in residues:
            # 提取残基链信息
            residue_chains = residue.get("residue_chains", [])
            if not residue_chains:
                continue
            
            # 获取参考链（通常是第一个或标记为 reference 的）
            ref_chain = None
            for chain in residue_chains:
                if chain.get("is_reference", False):
                    ref_chain = chain
                    break
            if not ref_chain:
                ref_chain = residue_chains[0]
            
            # 提取角色信息
            roles = residue.get("roles", [])
            role_names = [r.get("function", "") for r in roles]
            
            parsed_residues.append({
                "residue_type": ref_chain.get("code", ""),
                "residue_number": ref_chain.get("auth_resid", ref_chain.get("resid")),
                "chain": ref_chain.get("chain_name", ""),
                "pdb_id": ref_chain.get("pdb_id", ""),
                "roles": role_names,
                "roles_summary": residue.get("roles_summary", ""),
                "is_metal_ligand": any("metal" in r.lower() for r in role_names),
            })
        
        return parsed_residues
    
    def extract_metal_info(self, entry: Dict) -> Dict:
        """
        提取金属相关信息
        
        Args:
            entry: M-CSA 条目
            
        Returns:
            金属信息字典
        """
        metal_info = {
            "has_metal": False,
            "metal_ligands": [],
            "metal_types": [],
        }
        
        # 从描述中查找金属关键词
        description = entry.get("description", "").lower()
        metal_keywords = ["zinc", "iron", "copper", "manganese", "magnesium", 
                         "calcium", "nickel", "cobalt", "molybdenum", "metal"]
        
        found_metals = [kw for kw in metal_keywords if kw in description]
        if found_metals:
            metal_info["has_metal"] = True
            metal_info["metal_types"] = found_metals
        
        # 从催化残基中查找金属配体
        residues = entry.get("residues", [])
        parsed_residues = self.parse_catalytic_residues(residues)
        
        metal_ligands = [r for r in parsed_residues if r["is_metal_ligand"]]
        if metal_ligands:
            metal_info["has_metal"] = True
            metal_info["metal_ligands"] = metal_ligands
        
        # 从角色中查找金属相关功能
        for residue in residues:
            roles_summary = residue.get("roles_summary", "").lower()
            if "metal" in roles_summary:
                metal_info["has_metal"] = True
        
        return metal_info
    
    def process_entry(self, entry: Dict) -> Dict:
        """
        处理单个 M-CSA 条目，提取关键信息
        
        Args:
            entry: 原始条目数据
            
        Returns:
            处理后的条目信息
        """
        mcsa_id = entry.get("mcsa_id")
        
        # 基本信息
        processed = {
            "mcsa_id": mcsa_id,
            "enzyme_name": entry.get("enzyme_name", ""),
            "description": entry.get("description", ""),
            "uniprot_id": entry.get("reference_uniprot_id", ""),
            "ec_numbers": entry.get("all_ecs", []),
            "url": f"https://www.ebi.ac.uk/thornton-srv/m-csa/entry/{mcsa_id}/",
        }
        
        # 催化残基
        residues = entry.get("residues", [])
        processed["catalytic_residues"] = self.parse_catalytic_residues(residues)
        processed["num_catalytic_residues"] = len(processed["catalytic_residues"])
        
        # 金属信息
        processed["metal_info"] = self.extract_metal_info(entry)
        
        # 反应信息
        reaction = entry.get("reaction", {})
        if reaction:
            processed["reaction"] = {
                "name": reaction.get("name", ""),
                "equation": reaction.get("equation", ""),
            }
        
        return processed
    
    def process_all_entries(self, detailed_entries: List[Dict]) -> List[Dict]:
        """
        处理所有条目
        
        Args:
            detailed_entries: 详细条目列表
            
        Returns:
            处理后的条目列表
        """
        print("\n" + "=" * 80)
        print("处理和提取关键信息")
        print("=" * 80)
        
        processed_entries = []
        
        for entry in tqdm(detailed_entries, desc="处理条目"):
            try:
                processed = self.process_entry(entry)
                processed_entries.append(processed)
            except Exception as e:
                mcsa_id = entry.get("mcsa_id", "unknown")
                print(f"\n⚠️  处理 MCSA-{mcsa_id} 失败: {e}")
        
        print(f"\n✅ 成功处理 {len(processed_entries)} 个条目")
        
        # 保存处理后的数据
        with open(self.processed_file, 'w') as f:
            json.dump(processed_entries, f, indent=2)
        print(f"   已保存到: {self.processed_file}")
        
        return processed_entries
    
    def generate_summary(self, processed_entries: List[Dict]) -> Dict:
        """
        生成数据库摘要统计
        
        Args:
            processed_entries: 处理后的条目列表
            
        Returns:
            摘要统计字典
        """
        print("\n" + "=" * 80)
        print("生成数据库摘要")
        print("=" * 80)
        
        # 统计信息
        total_entries = len(processed_entries)
        entries_with_metal = sum(1 for e in processed_entries if e["metal_info"]["has_metal"])
        
        # EC 号统计
        ec_counter = {}
        for entry in processed_entries:
            for ec in entry["ec_numbers"]:
                ec_counter[ec] = ec_counter.get(ec, 0) + 1
        
        # 金属类型统计
        metal_types_counter = {}
        for entry in processed_entries:
            for metal in entry["metal_info"]["metal_types"]:
                metal_types_counter[metal] = metal_types_counter.get(metal, 0) + 1
        
        summary = {
            "total_entries": total_entries,
            "entries_with_metal": entries_with_metal,
            "metal_percentage": f"{entries_with_metal/total_entries*100:.1f}%",
            "total_ec_numbers": len(ec_counter),
            "top_ec_numbers": sorted(ec_counter.items(), key=lambda x: x[1], reverse=True)[:10],
            "metal_types_distribution": sorted(metal_types_counter.items(), key=lambda x: x[1], reverse=True),
            "total_catalytic_residues": sum(e["num_catalytic_residues"] for e in processed_entries),
        }
        
        # 打印摘要
        print(f"\n📊 M-CSA 数据库统计:")
        print(f"   总条目数: {summary['total_entries']}")
        print(f"   含金属的条目: {summary['entries_with_metal']} ({summary['metal_percentage']})")
        print(f"   覆盖的 EC 号: {summary['total_ec_numbers']}")
        print(f"   催化残基总数: {summary['total_catalytic_residues']}")
        
        print(f"\n   金属类型分布:")
        for metal, count in summary["metal_types_distribution"][:10]:
            print(f"     {metal}: {count}")
        
        print(f"\n   Top 10 EC 号:")
        for ec, count in summary["top_ec_numbers"]:
            print(f"     {ec}: {count}")
        
        # 保存摘要
        with open(self.summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"\n   已保存摘要到: {self.summary_file}")
        
        return summary
    
    def run(self):
        """执行完整的下载和处理流程"""
        start_time = time.time()
        
        # 1. 下载基本信息
        basic_entries = self.download_all_entries()
        
        # 2. 下载详细信息
        detailed_entries = self.download_all_detailed_entries(basic_entries)
        
        # 3. 处理和提取信息
        processed_entries = self.process_all_entries(detailed_entries)
        
        # 4. 生成摘要
        summary = self.generate_summary(processed_entries)
        
        elapsed = time.time() - start_time
        print("\n" + "=" * 80)
        print(f"✅ M-CSA 数据库下载完成！")
        print(f"   耗时: {elapsed/60:.1f} 分钟")
        print(f"   数据保存在: {self.output_dir}")
        print("=" * 80)
        
        return processed_entries, summary


def main():
    """主函数"""
    downloader = MCSADownloader()
    downloader.run()


if __name__ == "__main__":
    main()

