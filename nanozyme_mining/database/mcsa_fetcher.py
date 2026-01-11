"""
M-CSA (Mechanism and Catalytic Site Atlas) Data Fetcher
========================================================

Fetches catalytic mechanism and metal binding site information from M-CSA database.
M-CSA contains detailed annotations about enzyme catalytic sites including:
- Metal ions involved in catalysis
- Coordinating residues
- Catalytic mechanism descriptions

Reference: https://www.ebi.ac.uk/thornton-srv/m-csa/
"""

import requests
import json
import time
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict


@dataclass
class CatalyticSite:
    """Represents a catalytic site from M-CSA."""
    ec_number: str
    pdb_id: str
    chain_id: str
    uniprot_id: Optional[str]
    residues: List[Dict]  # List of catalytic residues
    cofactors: List[Dict]  # Metal ions and other cofactors
    mechanism: Optional[str]
    
    def to_dict(self):
        return asdict(self)


class MCSAFetcher:
    """Fetch catalytic site data from M-CSA database."""
    
    def __init__(self, cache_dir: str = "./cache/mcsa"):
        """
        Initialize M-CSA fetcher.
        
        Args:
            cache_dir: Directory for caching M-CSA data
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # M-CSA REST API endpoints
        self.base_url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"
        self.search_by_ec_url = f"{self.base_url}/entries"
        
        # Rate limiting
        self.request_delay = 0.5  # seconds between requests
    
    def query_by_ec(self, ec_number: str) -> List[CatalyticSite]:
        """
        Query M-CSA for catalytic sites by EC number.
        
        Args:
            ec_number: EC number (e.g., "1.11.1.6")
            
        Returns:
            List of CatalyticSite objects
        """
        cache_file = self.cache_dir / f"{ec_number}_mcsa.json"
        
        # Check cache first
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                cached_data = json.load(f)
                return [CatalyticSite(**site) for site in cached_data]
        
        print(f"  Querying M-CSA for EC {ec_number}...")
        
        try:
            # Query M-CSA API (with pagination)
            params = {"ec_number": ec_number}
            response = requests.get(
                self.search_by_ec_url, 
                params=params,
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            sites = []
            # M-CSA returns paginated results
            results = data.get("results", [])
            if not results and isinstance(data, list):
                # Fallback: sometimes API returns a list directly
                results = data
            
            for entry in results:
                site = self._parse_entry(entry, ec_number)
                if site:
                    sites.append(site)
            
            # Save to cache
            with open(cache_file, 'w') as f:
                json.dump([site.to_dict() for site in sites], f, indent=2)
            
            print(f"  ✓ Found {len(sites)} catalytic sites in M-CSA")
            time.sleep(self.request_delay)
            
            return sites
        
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                print(f"  ⚠️  No M-CSA entries found for EC {ec_number}")
                # Save empty cache to avoid repeated queries
                with open(cache_file, 'w') as f:
                    json.dump([], f)
                return []
            else:
                print(f"  ❌ Error querying M-CSA: {e}")
                return []
        
        except Exception as e:
            print(f"  ❌ Error querying M-CSA for EC {ec_number}: {e}")
            return []
    
    def _parse_entry(self, entry: Dict, ec_number: str) -> Optional[CatalyticSite]:
        """Parse a single M-CSA entry."""
        try:
            pdb_id = entry.get("pdb_code", "")
            chain_id = entry.get("chain_id", "")
            uniprot_id = entry.get("uniprot_accession")
            
            # Extract catalytic residues
            residues = []
            for res in entry.get("residues", []):
                residues.append({
                    "residue_name": res.get("residue_name"),
                    "residue_number": res.get("residue_number"),
                    "chain_id": res.get("chain_id"),
                    "role": res.get("chemical_function", "")
                })
            
            # Extract cofactors (including metals)
            cofactors = []
            for cof in entry.get("cofactors", []):
                cofactor_data = {
                    "name": cof.get("name"),
                    "type": cof.get("type"),
                    "formula": cof.get("formula", ""),
                }
                
                # Check if it's a metal ion
                metal_keywords = ["iron", "copper", "zinc", "manganese", "cobalt", 
                                 "nickel", "heme", "metal"]
                is_metal = any(kw in cof.get("name", "").lower() for kw in metal_keywords)
                cofactor_data["is_metal"] = is_metal
                
                cofactors.append(cofactor_data)
            
            mechanism = entry.get("mechanism_description")
            
            return CatalyticSite(
                ec_number=ec_number,
                pdb_id=pdb_id,
                chain_id=chain_id,
                uniprot_id=uniprot_id,
                residues=residues,
                cofactors=cofactors,
                mechanism=mechanism
            )
        
        except Exception as e:
            print(f"  ⚠️  Error parsing M-CSA entry: {e}")
            return None
    
    def fetch_for_nanozyme_types(
        self, 
        ec_to_nanozyme: Dict[str, str]
    ) -> Dict[str, List[CatalyticSite]]:
        """
        Fetch M-CSA data for multiple EC numbers.
        
        Args:
            ec_to_nanozyme: Dictionary mapping EC numbers to nanozyme types
            
        Returns:
            Dictionary mapping EC numbers to lists of catalytic sites
        """
        results = {}
        
        print(f"\n{'='*80}")
        print(f"Querying M-CSA for {len(ec_to_nanozyme)} EC numbers...")
        print(f"{'='*80}\n")
        
        for i, (ec_number, nanozyme_type) in enumerate(ec_to_nanozyme.items()):
            print(f"[{i+1}/{len(ec_to_nanozyme)}] {nanozyme_type} (EC {ec_number})")
            sites = self.query_by_ec(ec_number)
            results[ec_number] = sites
            
            # Print metal cofactors found
            for site in sites:
                metal_cofactors = [c for c in site.cofactors if c.get("is_metal")]
                if metal_cofactors:
                    print(f"    ✓ PDB {site.pdb_id}: {len(metal_cofactors)} metal cofactor(s)")
                    for metal in metal_cofactors:
                        print(f"      - {metal['name']}")
        
        return results
    
    def generate_metal_summary(
        self, 
        mcsa_results: Dict[str, List[CatalyticSite]]
    ) -> Dict:
        """Generate summary statistics of metal cofactors."""
        summary = {
            "total_sites": 0,
            "sites_with_metals": 0,
            "metal_types": {},
            "ec_metal_mapping": {}
        }
        
        for ec_number, sites in mcsa_results.items():
            summary["total_sites"] += len(sites)
            metal_types_for_ec = set()
            
            for site in sites:
                metal_cofactors = [c for c in site.cofactors if c.get("is_metal")]
                
                if metal_cofactors:
                    summary["sites_with_metals"] += 1
                    
                    for metal in metal_cofactors:
                        metal_name = metal["name"]
                        summary["metal_types"][metal_name] = \
                            summary["metal_types"].get(metal_name, 0) + 1
                        metal_types_for_ec.add(metal_name)
            
            if metal_types_for_ec:
                summary["ec_metal_mapping"][ec_number] = sorted(metal_types_for_ec)
        
        return summary


def main():
    """Example usage."""
    from ..utils.ec_mappings import EC_PATTERNS
    
    fetcher = MCSAFetcher()
    
    # Fetch for core nanozyme types
    core_ec_numbers = {
        "1.15.1.1": "SOD",
        "1.11.1.6": "CAT",
        "1.11.1.7": "POD",
        "1.10.3.2": "LAC"
    }
    
    results = fetcher.fetch_for_nanozyme_types(core_ec_numbers)
    
    # Generate summary
    summary = fetcher.generate_metal_summary(results)
    
    print("\n" + "="*80)
    print("M-CSA METAL COFACTOR SUMMARY")
    print("="*80)
    print(f"Total catalytic sites: {summary['total_sites']}")
    print(f"Sites with metal cofactors: {summary['sites_with_metals']}")
    print(f"\nMetal types found:")
    for metal, count in sorted(summary["metal_types"].items(), key=lambda x: -x[1]):
        print(f"  {metal:30} : {count}")
    print(f"\nEC to Metal mapping:")
    for ec, metals in summary["ec_metal_mapping"].items():
        print(f"  EC {ec}: {', '.join(metals)}")
    
    # Save summary
    output_file = Path("metal_cofactor_summary.json")
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n✓ Saved summary to {output_file}")


if __name__ == "__main__":
    main()
