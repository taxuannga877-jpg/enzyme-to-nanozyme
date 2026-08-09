"""
UniProt Data Fetcher - Stage 1 Data Acquisition
================================================

Fetches enzyme data from UniProt and downloads experimental PDB structures.
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
from ..utils.http_utils import make_session, RateLimiter

# PR1-3 (H7 fix): single shared Session w/ User-Agent + automatic retry, plus
# polite rate limit (~3 req/s). UniProt's anonymous quotas are not formally
# documented but they explicitly request a UA + considerate cadence.
_UNIPROT_SESSION = make_session()
_UNIPROT_RATE = RateLimiter(min_interval=0.3)


class UniProtFetcher:
    """
    Fetches enzyme data from UniProt REST API and downloads structures.

    NOTE: This class does NOT use any ML models.
    It directly queries UniProt API to get:
    - Enzyme entries by EC number
    - Known active site annotations (ft_act_site, ft_binding)
    - Experimental PDB cross references
    """

    def __init__(
        self,
        cache_dir: str = "./cache",
        pdb_library_dir: Optional[str] = None,
        download_timeout: int = 60,
        max_sequence_length: int = 600,
    ):
        """
        Initialize UniProt fetcher.

        Args:
            cache_dir: Directory for caching downloaded data
            pdb_library_dir: PDB library directory organized by EC numbers
            download_timeout: Timeout for downloads in seconds
            max_sequence_length: Maximum protein sequence length
        """
        self.cache_dir = Path(cache_dir)
        self.pdb_library_dir = Path(pdb_library_dir) if pdb_library_dir else None
        self.download_timeout = download_timeout
        self.max_sequence_length = max_sequence_length

        self.csv_cache = self.cache_dir / "csv"
        self.json_cache = self.cache_dir / "json"
        self.pdb_cache = self.cache_dir / "pdb"
        self.annotated_dir = self.cache_dir / "annotated"
        self.unannotated_dir = self.cache_dir / "unannotated"

        for d in [self.csv_cache, self.json_cache, self.pdb_cache,
                  self.annotated_dir, self.unannotated_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # API templates - TSV format for basic info
        self.uniprot_query_template = (
            "https://rest.uniprot.org/uniprotkb/stream?"
            "query=ec:{ec}+AND+reviewed:true"
            "&fields=accession,ec,sequence,xref_pdb"
            "&format=tsv&size={size}"
        )

        # JSON format for detailed info including active sites and metal binding
        # Note: ft_metal is not a valid field, removed it
        self.uniprot_json_template = (
            "https://rest.uniprot.org/uniprotkb/search?"
            "query=ec:{ec}+AND+reviewed:true"
            "&fields=accession,ec,sequence,ft_act_site,ft_binding,ft_site,xref_pdb"
            "&format=json"
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
            # PR1-3 (H7): polite throttle + shared Session
            _UNIPROT_RATE.wait()
            response = _UNIPROT_SESSION.get(url, timeout=self.download_timeout)
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
                # PR1-3 (H7): polite throttle + shared Session
                _UNIPROT_RATE.wait()
                response = _UNIPROT_SESSION.get(url, timeout=self.download_timeout)
                response.raise_for_status()
                data = response.json()

                page_results = data.get("results", [])
                if not page_results:
                    break

                for entry in page_results:
                    active_sites = self._extract_active_sites(entry)
                    # Extract PDB ID from xref_pdb (experimental structures)
                    # 优先使用所有PDB ID，如果没有则使用第一个
                    all_pdb_ids = self._extract_all_pdb_ids(entry)
                    pdb_id = all_pdb_ids[0] if all_pdb_ids else None
                    results.append({
                        "uniprot_id": entry.get("primaryAccession", ""),
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
                    "description": feat.get("description", ""),
                    "metal_role": self._parse_metal_role(feat.get("description", "")) if feat_type == "Metal binding" else None,
                })

        return sites

    @staticmethod
    def _parse_metal_role(description: str) -> str:
        desc_lower = description.lower()
        if "catalytic" in desc_lower:
            return "catalytic"
        if "structural" in desc_lower:
            return "structural"
        if "substrate" in desc_lower:
            return "substrate_binding"
        return "unknown"

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
            # PR1-3 (H7): polite throttle + shared Session. RCSB search API also
            # benefits from connection reuse + automatic retry on 5xx.
            _UNIPROT_RATE.wait()
            response = _UNIPROT_SESSION.post(
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
            # PR1-3 (H7): polite throttle + shared Session
            _UNIPROT_RATE.wait()
            response = _UNIPROT_SESSION.get(url, timeout=self.download_timeout)
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