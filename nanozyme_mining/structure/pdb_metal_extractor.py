"""
PDB Metal Extractor - Extract Metal Binding Sites from PDB Files
==================================================================

Parses PDB structure files to identify:
- Metal ions (Fe, Cu, Zn, Mn, Ca, Mg, Co, Ni, etc.)
- Metal-coordinating residues (within 3.5 Å)
- Metal coordination geometry
"""

import os
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict


# Common metal ions in proteins
METAL_IONS = {
    # Transition metals (critical for catalysis)
    "FE": "Iron",
    "FE2": "Iron(II)",  # Ferrous
    "FE3": "Iron(III)",  # Ferric
    "CU": "Copper",
    "CU1": "Copper(I)",  # Cuprous
    "CU2": "Copper(II)",  # Cupric
    "ZN": "Zinc",
    "ZN2": "Zinc(II)",
    "MN": "Manganese",
    "MN2": "Manganese(II)",
    "CO": "Cobalt",
    "CO2": "Cobalt(II)",
    "NI": "Nickel",
    "NI2": "Nickel(II)",
    
    # Alkali and alkaline earth metals
    "CA": "Calcium",
    "MG": "Magnesium",
    "K": "Potassium",
    "NA": "Sodium",
    
    # Other metals
    "MO": "Molybdenum",
    "W": "Tungsten",
    "V": "Vanadium",
    "SE": "Selenium",
    
    # Heme groups (contain Fe)
    "HEM": "Heme",
    "HEC": "Heme C",
}

# Metal-coordinating amino acids
COORDINATING_RESIDUES = {
    "HIS": ["NE2", "ND1"],  # Histidine - nitrogen
    "CYS": ["SG"],          # Cysteine - sulfur
    "MET": ["SD"],          # Methionine - sulfur
    "ASP": ["OD1", "OD2"],  # Aspartate - oxygen
    "GLU": ["OE1", "OE2"],  # Glutamate - oxygen
    "SER": ["OG"],          # Serine - oxygen
    "THR": ["OG1"],         # Threonine - oxygen
    "TYR": ["OH"],          # Tyrosine - oxygen
}

# Distance cutoff for metal coordination (Angstroms)
COORDINATION_DISTANCE = 3.5


@dataclass
class MetalSite:
    """Represents a metal binding site in a protein."""
    metal_type: str
    metal_name: str
    metal_residue_id: int
    metal_chain: str
    metal_coords: Tuple[float, float, float]
    coordinating_residues: List[Dict]
    coordination_number: int
    
    def to_dict(self):
        return asdict(self)


class PDBMetalExtractor:
    """Extract metal binding sites from PDB files."""
    
    def __init__(self):
        """Initialize the PDB metal extractor."""
        self.metal_ions = METAL_IONS
        self.coordinating_residues = COORDINATING_RESIDUES
        self.coordination_distance = COORDINATION_DISTANCE
    
    def parse_pdb_file(self, pdb_file: Path) -> List[MetalSite]:
        """
        Parse a PDB file and extract metal binding sites.
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            List of MetalSite objects
        """
        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        # Parse PDB file
        atoms = self._parse_atoms(pdb_file)
        
        # Separate metal ions and protein atoms
        metal_atoms = [a for a in atoms if a["res_name"] in self.metal_ions]
        protein_atoms = [a for a in atoms if a["res_name"] not in self.metal_ions]
        
        # Find coordinating residues for each metal
        metal_sites = []
        for metal in metal_atoms:
            coordinating = self._find_coordinating_residues(metal, protein_atoms)
            
            if coordinating:  # Only include if there are coordinating residues
                site = MetalSite(
                    metal_type=metal["res_name"],
                    metal_name=self.metal_ions.get(metal["res_name"], "Unknown"),
                    metal_residue_id=metal["res_id"],
                    metal_chain=metal["chain"],
                    metal_coords=(metal["x"], metal["y"], metal["z"]),
                    coordinating_residues=coordinating,
                    coordination_number=len(coordinating)
                )
                metal_sites.append(site)
        
        return metal_sites
    
    def _parse_atoms(self, pdb_file: Path) -> List[Dict]:
        """Parse ATOM and HETATM records from PDB file."""
        atoms = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        atom = {
                            "record": line[0:6].strip(),
                            "atom_id": int(line[6:11].strip()),
                            "atom_name": line[12:16].strip(),
                            "res_name": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "res_id": int(line[22:26].strip()),
                            "x": float(line[30:38].strip()),
                            "y": float(line[38:46].strip()),
                            "z": float(line[46:54].strip()),
                            "element": line[76:78].strip() if len(line) > 76 else ""
                        }
                        atoms.append(atom)
                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue
        
        return atoms
    
    def _find_coordinating_residues(
        self, 
        metal_atom: Dict, 
        protein_atoms: List[Dict]
    ) -> List[Dict]:
        """Find residues coordinating a metal ion."""
        coordinating = []
        metal_coords = np.array([metal_atom["x"], metal_atom["y"], metal_atom["z"]])
        
        # Group atoms by residue
        residues = {}
        for atom in protein_atoms:
            res_key = (atom["chain"], atom["res_id"], atom["res_name"])
            if res_key not in residues:
                residues[res_key] = []
            residues[res_key].append(atom)
        
        # Check each residue for coordination
        for (chain, res_id, res_name), atoms in residues.items():
            if res_name not in self.coordinating_residues:
                continue
            
            # Check specific coordinating atoms
            coord_atoms = self.coordinating_residues[res_name]
            for atom in atoms:
                if atom["atom_name"] in coord_atoms:
                    atom_coords = np.array([atom["x"], atom["y"], atom["z"]])
                    distance = np.linalg.norm(metal_coords - atom_coords)
                    
                    if distance <= self.coordination_distance:
                        coordinating.append({
                            "residue_name": res_name,
                            "residue_id": res_id,
                            "chain": chain,
                            "atom_name": atom["atom_name"],
                            "distance": round(distance, 2)
                        })
                        break  # Only one coordinating atom per residue
        
        return coordinating
    
    def extract_from_directory(
        self, 
        pdb_dir: Path, 
        output_dir: Optional[Path] = None
    ) -> Dict[str, List[MetalSite]]:
        """
        Extract metal sites from all PDB files in a directory.
        
        Args:
            pdb_dir: Directory containing PDB files
            output_dir: Optional directory to save JSON output
            
        Returns:
            Dictionary mapping PDB filename to list of metal sites
        """
        if not pdb_dir.exists():
            raise FileNotFoundError(f"Directory not found: {pdb_dir}")
        
        pdb_files = sorted(pdb_dir.glob("*.pdb"))
        print(f"Found {len(pdb_files)} PDB files")
        
        results = {}
        metals_found = 0
        files_with_metals = 0
        
        for i, pdb_file in enumerate(pdb_files):
            try:
                metal_sites = self.parse_pdb_file(pdb_file)
                
                if metal_sites:
                    results[pdb_file.name] = [site.to_dict() for site in metal_sites]
                    metals_found += len(metal_sites)
                    files_with_metals += 1
                
                if (i + 1) % 100 == 0:
                    print(f"  Processed {i+1}/{len(pdb_files)} files...", end='\r')
            
            except Exception as e:
                print(f"\n  ⚠️  Error processing {pdb_file.name}: {e}")
        
        print(f"\n✓ Found {metals_found} metal sites in {files_with_metals}/{len(pdb_files)} files")
        
        # Save results
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / "metal_sites_summary.json"
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"  Saved results to {output_file}")
        
        return results
    
    def generate_statistics(self, results: Dict[str, List[MetalSite]]) -> Dict:
        """Generate statistics about metal sites."""
        stats = {
            "total_files": len(results),
            "total_sites": sum(len(sites) for sites in results.values()),
            "metal_types": {},
            "coordination_numbers": {},
            "coordinating_residues": {}
        }
        
        for pdb_file, sites in results.items():
            for site in sites:
                # Count metal types
                metal = site["metal_type"]
                stats["metal_types"][metal] = stats["metal_types"].get(metal, 0) + 1
                
                # Count coordination numbers
                coord_num = site["coordination_number"]
                stats["coordination_numbers"][coord_num] = \
                    stats["coordination_numbers"].get(coord_num, 0) + 1
                
                # Count coordinating residues
                for res in site["coordinating_residues"]:
                    res_name = res["residue_name"]
                    stats["coordinating_residues"][res_name] = \
                        stats["coordinating_residues"].get(res_name, 0) + 1
        
        return stats


def main():
    """Example usage."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python pdb_metal_extractor.py <pdb_directory>")
        sys.exit(1)
    
    pdb_dir = Path(sys.argv[1])
    output_dir = Path("metal_sites_output")
    
    extractor = PDBMetalExtractor()
    results = extractor.extract_from_directory(pdb_dir, output_dir)
    
    # Generate and print statistics
    stats = extractor.generate_statistics(results)
    print("\n" + "="*80)
    print("STATISTICS:")
    print("="*80)
    print(f"Total PDB files with metals: {stats['total_files']}")
    print(f"Total metal sites: {stats['total_sites']}")
    print(f"\nMetal types:")
    for metal, count in sorted(stats['metal_types'].items(), key=lambda x: -x[1]):
        metal_name = METAL_IONS.get(metal, "Unknown")
        print(f"  {metal:6} ({metal_name:20}): {count}")
    print(f"\nCoordination numbers:")
    for coord_num, count in sorted(stats['coordination_numbers'].items()):
        print(f"  {coord_num}: {count}")
    print(f"\nCoordinating residues:")
    for res, count in sorted(stats['coordinating_residues'].items(), key=lambda x: -x[1]):
        print(f"  {res}: {count}")


if __name__ == "__main__":
    main()

