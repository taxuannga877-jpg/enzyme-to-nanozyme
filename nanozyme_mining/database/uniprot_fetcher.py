"""
UniProt Data Fetcher - Stage 1 Data Acquisition
================================================

Fetches enzyme data from UniProt and downloads AlphaFold structures.
Based on ChemEnzyRetroPlanner's UniProtParserEC architecture.
"""

import os
import json
import time
import subprocess
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .nanozyme_db import NanozymeDatabase, EnzymeEntry
from ..utils.constants import (
    NanozymeType,
    DEFAULT_UNIPROT_QUERY_SIZE,
    DEFAULT_MAX_ENTRIES,
    FETCH_ALL_RESULTS,
    UNIPROT_PAGE_SIZE
)
from ..utils.ec_mappings import EC_PATTERNS

# 导入 M-CSA 查询工具
try:
    from .mcsa_query import get_mcsa_query
    MCSA_AVAILABLE = True
except ImportError:
    MCSA_AVAILABLE = False
    print("⚠️  M-CSA 查询模块不可用")


AFDB_VERSION = "6"  # AlphaFold DB version (updated to v6)


class UniProtFetcher:
    """
    Fetches enzyme data from UniProt REST API and downloads structures.

    Based on ChemEnzyRetroPlanner's UniProtParserEC class.

    NOTE: This class does NOT use any ML models.
    It directly queries UniProt API to get:
    - Enzyme entries by EC number
    - Known active site annotations (ft_act_site, ft_binding)
    - AlphaFold structure downloads
    """

    def __init__(
        self,
        cache_dir: str = "./cache",
        pdb_library_dir: Optional[str] = None,
        download_timeout: int = 60,
        max_sequence_length: int = 600,
        use_mcsa: bool = True
    ):
        """
        Initialize UniProt fetcher.

        Args:
            cache_dir: Directory for caching downloaded data (legacy, for CSV and other caches)
            pdb_library_dir: PDB library directory organized by EC numbers (new structure)
            download_timeout: Timeout for downloads in seconds
            max_sequence_length: Maximum protein sequence length
            use_mcsa: Whether to query M-CSA for metal site information
        """
        self.cache_dir = Path(cache_dir)
        self.pdb_library_dir = Path(pdb_library_dir) if pdb_library_dir else None
        self.download_timeout = download_timeout
        self.max_sequence_length = max_sequence_length
        self.use_mcsa = use_mcsa and MCSA_AVAILABLE

        # Create cache directories (for CSV and other legacy caches)
        self.csv_cache = self.cache_dir / "csv"
        # JSON cache: use pdb_library if provided, otherwise use legacy cache/json
        self.json_cache = self.cache_dir / "json"  # Legacy path (for backward compatibility)
        self.pdb_cache = self.cache_dir / "pdb"  # Legacy PDB cache (for backward compatibility)

        # Separate folders for annotated vs unannotated
        self.annotated_dir = self.cache_dir / "annotated"
        self.unannotated_dir = self.cache_dir / "unannotated"
        
        # Initialize M-CSA query tool
        self.mcsa_query = None
        if self.use_mcsa:
            try:
                self.mcsa_query = get_mcsa_query()
                print("✓ M-CSA 数据库已启用")
            except Exception as e:
                print(f"⚠️  M-CSA 数据库加载失败: {e}")
                self.use_mcsa = False

        for d in [self.csv_cache, self.json_cache, self.pdb_cache,
                  self.annotated_dir, self.unannotated_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # API templates - TSV format for basic info
        self.uniprot_query_template = (
            "https://rest.uniprot.org/uniprotkb/stream?"
            "query=ec:{ec}+AND+reviewed:true+AND+database:(alphafolddb)"
            "&fields=accession,ec,sequence,xref_alphafolddb"
            "&format=tsv&size={size}"
        )

        # JSON format for detailed info including active sites and metal binding
        # Note: Removed database:(alphafolddb) filter as it may cause 400 errors
        # We'll filter for AlphaFold entries after fetching
        # Note: ft_metal is not a valid field, removed it
        self.uniprot_json_template = (
            "https://rest.uniprot.org/uniprotkb/search?"
            "query=ec:{ec}+AND+reviewed:true"
            "&fields=accession,ec,sequence,ft_act_site,ft_binding,ft_site,xref_alphafolddb,xref_pdb"
            "&format=json"
        )

        self.alphafold_url_template = (
            "https://alphafold.ebi.ac.uk/files/"
            "AF-{alphafold_id}-F1-model_v{version}.pdb"
        )
        
        # RCSB PDB URL template for experimental structures
        self.rcsb_pdb_url_template = "https://files.rcsb.org/download/{pdb_id}.pdb"

    def query_by_ec(
        self,
        ec_number: str,
        size: int = DEFAULT_UNIPROT_QUERY_SIZE
    ) -> List[Dict]:
        """
        Query UniProt for enzymes by EC number.

        Args:
            ec_number: EC number (e.g., "1.11.1.7")
            size: Maximum number of results (default: 500, use -1 for all results)

        Returns:
            List of enzyme data dictionaries
        """
        cache_file = self.csv_cache / f"{ec_number}.tsv"

        # Check cache first
        if cache_file.exists():
            return self._parse_tsv(cache_file)

        # Query UniProt API
        # Use large size for /stream endpoint (supports up to 500+ results)
        # If size is FETCH_ALL_RESULTS (-1), try to get all results by using a very large number
        query_size = 10000 if size == FETCH_ALL_RESULTS else size
        url = self.uniprot_query_template.format(
            ec=ec_number, size=query_size
        )

        try:
            response = requests.get(url, timeout=self.download_timeout)
            response.raise_for_status()

            # Save to cache
            with open(cache_file, 'w') as f:
                f.write(response.text)

            results = self._parse_tsv(cache_file)
            
            # If size was specified and we got fewer results, that's fine
            # If size was FETCH_ALL_RESULTS (-1), return all results
            if size > 0 and len(results) > size:
                results = results[:size]
            
            print(f"  ✓ Retrieved {len(results)} entries for EC {ec_number}")
            return results

        except Exception as e:
            print(f"Error querying UniProt for EC {ec_number}: {e}")
            return []

    def _parse_tsv(self, tsv_file: Path) -> List[Dict]:
        """Parse TSV file from UniProt."""
        results = []

        with open(tsv_file, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            return results

        headers = lines[0].strip().split('\t')

        for line in lines[1:]:
            values = line.strip().split('\t')
            if len(values) >= len(headers):
                entry = dict(zip(headers, values))
                results.append(entry)

        return results

    def _get_json_cache_path(self, ec_number: str) -> Path:
        """
        Get the JSON cache file path for an EC number.
        Uses pdb_library if available, otherwise falls back to legacy cache/json.
        
        Args:
            ec_number: EC number (e.g., "1.4.3.4")
            
        Returns:
            Path to JSON cache file
        """
        if self.pdb_library_dir:
            # New structure: pdb_library/{EC号}/{EC号}_sites.json
            ec_dir_name = ec_number.replace(".", "_")
            ec_dir = self.pdb_library_dir / ec_dir_name
            ec_dir.mkdir(parents=True, exist_ok=True)
            return ec_dir / f"{ec_number}_sites.json"
        else:
            # Legacy structure: cache/json/{EC号}_sites.json
            return self.json_cache / f"{ec_number}_sites.json"

    def query_with_active_sites(
        self, 
        ec_number: str,
        max_results: int = DEFAULT_UNIPROT_QUERY_SIZE
    ) -> List[Dict]:
        """
        Query UniProt with active site annotations (JSON format) with pagination.

        This method gets KNOWN active sites from UniProt annotations,
        NOT predicted by any ML model.

        Args:
            ec_number: EC number
            max_results: Maximum number of results to retrieve (default: 500, use FETCH_ALL_RESULTS=-1 for all)

        Returns:
            List of enzyme data with active site info
        """
        cache_file = self._get_json_cache_path(ec_number)

        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    cached_results = json.load(f)
                    # Check if cache is empty or invalid
                    if not cached_results or len(cached_results) == 0:
                        print(f"  ⚠️  Cache file is empty, re-fetching...")
                    # Check if cache has required fields (alphafold_id)
                    elif cached_results and len(cached_results) > 0 and "alphafold_id" not in cached_results[0]:
                        print(f"  ⚠️  Cache file missing alphafold_id, re-fetching...")
                    # If requesting all results and cache has fewer, re-fetch
                    elif max_results == FETCH_ALL_RESULTS and len(cached_results) < 100:
                        print(f"  ⚠️  Cache has only {len(cached_results)} entries, re-fetching...")
                    else:
                        return cached_results
            except (json.JSONDecodeError, KeyError, IndexError) as e:
                print(f"  ⚠️  Cache file corrupted, re-fetching... ({e})")

        results = []
        page_size = UNIPROT_PAGE_SIZE
        from_index = 0
        
        print(f"  Querying UniProt for EC {ec_number} (with pagination)...")

        while True:
            # Build URL with pagination
            url = (
                f"{self.uniprot_json_template.format(ec=ec_number)}"
                f"&size={page_size}&from={from_index}"
            )

            try:
                response = requests.get(url, timeout=self.download_timeout)
                response.raise_for_status()
                data = response.json()

                page_results = data.get("results", [])
                if not page_results:
                    break

                for entry in page_results:
                    active_sites = self._extract_active_sites(entry)
                    # Extract AlphaFold ID from xref_alphafolddb
                    alphafold_id = self._extract_alphafold_id(entry)
                    # Extract PDB ID from xref_pdb (experimental structures)
                    # 优先使用所有PDB ID，如果没有则使用第一个
                    all_pdb_ids = self._extract_all_pdb_ids(entry)
                    pdb_id = all_pdb_ids[0] if all_pdb_ids else None
                    results.append({
                        "uniprot_id": entry.get("primaryAccession", ""),
                        "alphafold_id": alphafold_id,
                        "pdb_id": pdb_id,  # Experimental PDB ID (first one for compatibility)
                        "pdb_ids": all_pdb_ids,  # All PDB IDs
                        "sequence": entry.get("sequence", {}).get("value", ""),
                        "active_sites": active_sites
                    })

                # Check if we've got all results
                total_results = data.get("totalResults", len(page_results))
                from_index += len(page_results)
                
                print(f"    Retrieved {from_index}/{total_results} entries...", end='\r')

                # Stop if we've got all results or reached max_results
                if from_index >= total_results:
                    break
                if max_results > 0 and len(results) >= max_results:
                    results = results[:max_results]
                    break

            except Exception as e:
                print(f"\n  ⚠️  Error querying page (from={from_index}): {e}")
                break

        print(f"\n  ✓ Retrieved {len(results)} entries with active site info")

        # Save to cache
        with open(cache_file, 'w') as f:
            json.dump(results, f, indent=2)

        return results

    def _extract_active_sites(self, entry: Dict) -> List[Dict]:
        """Extract active site annotations from UniProt entry, including metal binding sites."""
        sites = []
        features = entry.get("features", [])

        for feat in features:
            feat_type = feat.get("type", "")
            # Include metal binding sites for nanozyme design
            if feat_type in ["Active site", "Binding site", "Site", "Metal binding"]:
                location = feat.get("location", {})
                sites.append({
                    "type": feat_type,
                    "start": location.get("start", {}).get("value"),
                    "end": location.get("end", {}).get("value"),
                    "description": feat.get("description", "")
                })

        return sites

    def _extract_alphafold_id(self, entry: Dict) -> Optional[str]:
        """Extract AlphaFold database ID from UniProt entry."""
        cross_references = entry.get("uniProtKBCrossReferences", [])
        for xref in cross_references:
            if xref.get("database") == "AlphaFoldDB":
                # Extract ID from the URL or id field
                xref_id = xref.get("id", "")
                if xref_id:
                    return xref_id
                # Alternative: extract from properties
                properties = xref.get("properties", {})
                if "id" in properties:
                    return properties["id"]
        return None

    def _extract_pdb_id(self, entry: Dict) -> Optional[str]:
        """Extract PDB database ID from UniProt entry (experimental structures).
        
        Returns the first available PDB ID, or None if no experimental structure exists.
        """
        cross_references = entry.get("uniProtKBCrossReferences", [])
        for xref in cross_references:
            if xref.get("database") == "PDB":
                pdb_id = xref.get("id")
                if pdb_id:
                    return pdb_id.upper()  # PDB IDs are uppercase

        return None

    def _extract_all_pdb_ids(self, entry: Dict) -> List[str]:
        """Extract ALL PDB database IDs from UniProt entry (experimental structures).
        
        Returns a list of all available PDB IDs, or empty list if no experimental structure exists.
        """
        pdb_ids = []
        cross_references = entry.get("uniProtKBCrossReferences", [])
        for xref in cross_references:
            if xref.get("database") == "PDB":
                pdb_id = xref.get("id")
                if pdb_id:
                    pdb_ids.append(pdb_id.upper())  # PDB IDs are uppercase

        return pdb_ids

    def query_rcsb_pdb_by_ec(self, ec_number: str, timeout: int = 30) -> List[str]:
        """
        直接从RCSB PDB查询指定EC号的所有实验结构PDB ID
        
        Args:
            ec_number: EC号 (e.g., "1.1.1.1")
            timeout: 请求超时时间（秒）
        
        Returns:
            PDB ID列表
        """
        # RCSB PDB Search API v2 endpoint
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        # 构建查询JSON - 使用full_text服务查询EC号
        # 注意：text服务不支持rcsb_ec_lineage.id属性，需要使用full_text服务
        query_json = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": f"EC {ec_number}"
                }
            },
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }
        
        try:
            response = requests.post(
                search_url,
                json=query_json,
                headers={"Content-Type": "application/json"},
                timeout=timeout
            )
            response.raise_for_status()
            
            result = response.json()
            pdb_ids = []
            
            # 解析结果
            if "result_set" in result:
                for item in result["result_set"]:
                    if "identifier" in item:
                        pdb_ids.append(item["identifier"].upper())
            
            return sorted(set(pdb_ids))  # 去重并排序
            
        except Exception as e:
            print(f"  ⚠️  查询RCSB PDB失败 (EC {ec_number}): {e}")
            return []

    def download_experimental_pdb(self, pdb_id: str, target_dir: Optional[Path] = None) -> Optional[Path]:
        """
        Download experimental PDB structure from RCSB PDB.

        Args:
            pdb_id: PDB database ID (e.g., "1ABC")
            target_dir: Target directory for download (default: pdb_cache)

        Returns:
            Path to downloaded PDB file or None
        """
        if target_dir is None:
            target_dir = self.pdb_cache
        target_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_id = pdb_id.upper()
        pdb_file = target_dir / f"{pdb_id}.pdb"

        if pdb_file.exists() and pdb_file.stat().st_size > 0:
            return pdb_file

        url = self.rcsb_pdb_url_template.format(pdb_id=pdb_id)

        try:
            response = requests.get(url, timeout=self.download_timeout)
            response.raise_for_status()

            # Check if response is valid (not a 404 page)
            if len(response.text) < 1000:  # RCSB 404 pages are usually small
                print(f"Warning: Downloaded file for {pdb_id} seems too small, might be invalid")
            
            with open(pdb_file, 'w') as f:
                f.write(response.text)

            return pdb_file

        except Exception as e:
            print(f"Error downloading experimental PDB for {pdb_id}: {e}")
            return None

    def download_pdb(self, alphafold_id: str, target_dir: Optional[Path] = None) -> Optional[Path]:
        """
        Download AlphaFold PDB structure.

        Args:
            alphafold_id: AlphaFold database ID
            target_dir: Target directory for download (default: pdb_cache)

        Returns:
            Path to downloaded PDB file or None
        """
        if target_dir is None:
            target_dir = self.pdb_cache
        target_dir.mkdir(parents=True, exist_ok=True)

        pdb_file = target_dir / f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"

        if pdb_file.exists() and pdb_file.stat().st_size > 0:
            return pdb_file

        url = self.alphafold_url_template.format(
            alphafold_id=alphafold_id,
            version=AFDB_VERSION
        )

        try:
            response = requests.get(url, timeout=self.download_timeout)
            response.raise_for_status()

            with open(pdb_file, 'w') as f:
                f.write(response.text)

            return pdb_file

        except Exception as e:
            print(f"Error downloading PDB for {alphafold_id}: {e}")
            return None

    def download_pdb_prioritized(
        self, 
        pdb_id: Optional[str] = None,
        alphafold_id: Optional[str] = None,
        target_dir: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Download PDB structure with priority: experimental PDB > AlphaFold.
        
        Args:
            pdb_id: Experimental PDB ID (optional, highest priority)
            alphafold_id: AlphaFold ID (optional, fallback)
            target_dir: Target directory for download (default: pdb_cache)
            
        Returns:
            Path to downloaded PDB file or None
        """
        # Priority 1: Try experimental PDB
        if pdb_id:
            pdb_path = self.download_experimental_pdb(pdb_id, target_dir)
            if pdb_path:
                return pdb_path
        
        # Priority 2: Fallback to AlphaFold
        if alphafold_id:
            return self.download_pdb(alphafold_id, target_dir)
        
            return None

    def fetch_and_populate(
        self,
        db: NanozymeDatabase,
        ec_number: str,
        nanozyme_type: NanozymeType,
        max_entries: int = DEFAULT_MAX_ENTRIES
    ) -> int:
        """
        Fetch enzymes by EC and populate database.

        Args:
            db: NanozymeDatabase instance
            ec_number: EC number to query
            nanozyme_type: Nanozyme type for this EC
            max_entries: Maximum entries to fetch

        Returns:
            Number of entries added
        """
        entries_data = self.query_by_ec(ec_number, max_entries)
        count = 0

        for data in entries_data:
            try:
                # Extract AlphaFold ID
                alphafold_id = data.get('AlphaFoldDB', '').split(';')[0]
                if not alphafold_id:
                    continue

                sequence = data.get('Sequence', '')
                if len(sequence) > self.max_sequence_length:
                    continue

                # Download PDB
                pdb_path = self.download_pdb(alphafold_id)

                entry = EnzymeEntry(
                    uniprot_id=data.get('Entry', ''),
                    ec_number=ec_number,
                    nanozyme_type=nanozyme_type.value,
                    sequence=sequence,
                    sequence_length=len(sequence),
                    alphafold_id=alphafold_id,
                    pdb_path=str(pdb_path) if pdb_path else None
                )

                if db.add_enzyme(entry):
                    count += 1

            except Exception as e:
                print(f"Error processing entry: {e}")

        return count

    def _enrich_with_mcsa(self, entry_data: Dict, ec_number: str) -> Dict:
        """
        Enrich entry data with M-CSA metal site information.
        
        Args:
            entry_data: Entry data dictionary
            ec_number: EC number
            
        Returns:
            Enriched entry data
        """
        try:
            # Query M-CSA for metal sites
            mcsa_data = self.mcsa_query.get_metal_sites(ec_number)
            
            if mcsa_data["has_metal"]:
                # Add M-CSA metal information
                if "metal_info" not in entry_data:
                    entry_data["metal_info"] = {}
                
                entry_data["metal_info"].update({
                    "has_metal": True,
                    "metal_types": mcsa_data["metal_types"],
                    "metal_coordination": mcsa_data["metal_coordination"],
                    "mcsa_references": mcsa_data["mcsa_references"],
                })
                
                # Merge M-CSA catalytic residues into active_sites
                active_sites = entry_data.get("active_sites", [])
                
                for residue in mcsa_data["catalytic_residues"]:
                    # Convert M-CSA residue format to active_site format
                    site = {
                        "position": residue.get("residue_number"),
                        "type": residue.get("residue_type"),
                        "description": f"M-CSA catalytic residue ({', '.join(residue.get('roles', []))})",
                        "source": "M-CSA",
                        "is_metal_ligand": residue.get("is_metal_ligand", False),
                        "mcsa_id": residue.get("mcsa_id"),
                    }
                    active_sites.append(site)
                
                entry_data["active_sites"] = active_sites
                
        except Exception as e:
            print(f"  ⚠️  M-CSA enrichment failed for {ec_number}: {e}")
        
        return entry_data
    
    def fetch_and_classify(
        self,
        ec_number: str,
        nanozyme_type: NanozymeType,
        max_results: int = DEFAULT_UNIPROT_QUERY_SIZE
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Fetch enzymes and classify into annotated vs unannotated.

        Two-pronged strategy:
        - Annotated: Has active site info from UniProt/M-CSA
        - Unannotated: No active site info, needs model prediction

        Args:
            ec_number: EC number to query
            nanozyme_type: Nanozyme type
            max_results: Maximum number of results to fetch (default: 500, use -1 for all)

        Returns:
            Tuple of (annotated_list, unannotated_list)
        """
        # Get detailed info with active sites (with pagination)
        entries = self.query_with_active_sites(ec_number, max_results=max_results)

        annotated = []
        unannotated = []

        for entry in entries:
            uniprot_id = entry.get("uniprot_id", "")
            alphafold_id = entry.get("alphafold_id", "")
            active_sites = entry.get("active_sites", [])
            sequence = entry.get("sequence", "")

            if len(sequence) > self.max_sequence_length:
                continue

            # Download PDB using AlphaFold ID
            pdb_path = None
            if alphafold_id:
                pdb_path = self.download_pdb(alphafold_id)
            else:
                print(f"  ⚠️  No AlphaFold ID for {uniprot_id}, skipping PDB download")

            entry_data = {
                "uniprot_id": uniprot_id,
                "ec_number": ec_number,
                "nanozyme_type": nanozyme_type.value,
                "sequence": sequence,
                "pdb_path": str(pdb_path) if pdb_path else None,
                "active_sites": active_sites
            }
            
            # Enrich with M-CSA metal site information
            if self.use_mcsa and self.mcsa_query:
                entry_data = self._enrich_with_mcsa(entry_data, ec_number)

            # Classify based on annotation
            if entry_data["active_sites"] and len(entry_data["active_sites"]) > 0:
                annotated.append(entry_data)
                self._save_to_folder(entry_data, self.annotated_dir)
            else:
                unannotated.append(entry_data)
                self._save_to_folder(entry_data, self.unannotated_dir)

        print(f"EC {ec_number}: {len(annotated)} annotated, {len(unannotated)} unannotated")
        return annotated, unannotated

    def _save_to_folder(self, entry_data: Dict, folder: Path):
        """Save entry data to JSON file in specified folder."""
        uniprot_id = entry_data["uniprot_id"]
        output_file = folder / f"{uniprot_id}.json"
        with open(output_file, 'w') as f:
            json.dump(entry_data, f, indent=2)
