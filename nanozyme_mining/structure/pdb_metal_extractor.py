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
from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass, asdict, field


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
    coordination_geometry: str = "unknown"  # e.g., "tetrahedral", "octahedral", "square_planar"
    oxidation_state: Optional[str] = None  # e.g., "Fe(II)", "Cu(II)"
    coordination_distances: List[float] = field(default_factory=list)  # Distances to coordinating atoms
    
    def to_dict(self):
        """Convert to dictionary with numpy type conversion."""
        result = asdict(self)
        # Convert numpy types to Python native types
        return self._convert_numpy_types(result)
    
    def _convert_numpy_types(self, obj: Any) -> Any:
        """Recursively convert numpy types to Python native types."""
        if isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._convert_numpy_types(item) for item in obj]
        else:
            return obj


class PDBMetalExtractor:
    """Extract metal binding sites from PDB files."""
    
    def __init__(self):
        """Initialize the PDB metal extractor."""
        self.metal_ions = METAL_IONS
        self.coordinating_residues = COORDINATING_RESIDUES
        self.coordination_distance = COORDINATION_DISTANCE
    
    def _is_metal_element(self, element: str) -> bool:
        """
        Check if an element symbol represents a metal.
        
        Args:
            element: Element symbol (e.g., "FE", "CU", "ZN")
            
        Returns:
            True if the element is a metal
        """
        if not element:
            return False
        
        element_upper = element.upper().strip()
        
        # Common metal elements
        metal_elements = {
            "FE", "CU", "ZN", "MN", "CA", "MG", "CO", "NI",
            "MO", "W", "V", "CR", "TI", "AL", "PB", "HG",
            "CD", "AG", "AU", "PT", "PD", "RU", "RH", "IR",
            "OS", "RE", "TC", "NB", "TA", "HF", "ZR", "Y",
            "SC", "LA", "CE", "PR", "ND", "PM", "SM", "EU",
            "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
            "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM",
            "BK", "CF", "ES", "FM", "MD", "NO", "LR"
        }
        
        return element_upper in metal_elements
    
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
        # Check both res_name (for HETATM records) and element (for ATOM records)
        metal_atoms = []
        protein_atoms = []
        
        for atom in atoms:
            res_name = atom.get("res_name", "")
            element = atom.get("element", "")
            
            # Check if it's a metal: either res_name is in metal_ions dict, or element is a metal
            is_metal = (res_name in self.metal_ions) or self._is_metal_element(element)
            
            if is_metal:
                metal_atoms.append(atom)
            else:
                protein_atoms.append(atom)
        
        # Find coordinating residues for each metal
        metal_sites = []
        for metal in metal_atoms:
            coordinating, distances = self._find_coordinating_residues(metal, protein_atoms)
            
            if coordinating:  # Only include if there are coordinating residues
                # Determine metal type and name
                # Prefer res_name if it's in metal_ions dict, otherwise use element
                res_name = metal.get("res_name", "")
                element = metal.get("element", "").upper()
                
                if res_name in self.metal_ions:
                    metal_type = res_name
                    metal_name = self.metal_ions[res_name]
                elif element:
                    # Use element as metal_type, and try to get a readable name
                    metal_type = element
                    # Try to find a matching metal name
                    metal_name = self.metal_ions.get(element, f"{element} (metal)")
                else:
                    metal_type = res_name if res_name else "UNK"
                    metal_name = "Unknown Metal"
                
                # Analyze coordination geometry
                geometry = self._analyze_coordination_geometry(
                    metal, coordinating, distances
                )
                
                # Determine oxidation state (use metal_type which might be res_name or element)
                oxidation_state = self._determine_oxidation_state(
                    metal_type, len(coordinating), geometry
                )
                
                site = MetalSite(
                    metal_type=metal_type,
                    metal_name=metal_name,
                    metal_residue_id=metal["res_id"],
                    metal_chain=metal["chain"],
                    metal_coords=(metal["x"], metal["y"], metal["z"]),
                    coordinating_residues=coordinating,
                    coordination_number=len(coordinating),
                    coordination_geometry=geometry,
                    oxidation_state=oxidation_state,
                    coordination_distances=distances
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
    ) -> Tuple[List[Dict], List[float]]:
        """
        Find residues coordinating a metal ion.
        
        Returns:
            Tuple of (coordinating_residues, distances)
        """
        coordinating = []
        distances = []
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
                        # Convert numpy distance to Python float
                        distance_float = float(distance)
                        coordinating.append({
                            "residue_name": res_name,
                            "residue_id": res_id,
                            "chain": chain,
                            "atom_name": atom["atom_name"],
                            "distance": round(distance_float, 2),
                            "coordinates": [atom["x"], atom["y"], atom["z"]]
                        })
                        distances.append(distance_float)
                        break  # Only one coordinating atom per residue
        
        return coordinating, distances
    
    def _analyze_coordination_geometry(
        self,
        metal: Dict,
        coordinating: List[Dict],
        distances: List[float]
    ) -> str:
        """
        Analyze coordination geometry based on coordination number and angles.
        
        Args:
            metal: Metal atom dictionary
            coordinating: List of coordinating residue dictionaries
            distances: List of coordination distances
            
        Returns:
            Geometry type string
        """
        coord_num = len(coordinating)
        
        if coord_num == 2:
            return "linear"
        elif coord_num == 3:
            # Could be trigonal planar or T-shaped
            return "trigonal_planar"
        elif coord_num == 4:
            # Check angles to distinguish tetrahedral from square planar
            if len(coordinating) >= 4:
                angles = self._calculate_coordination_angles(metal, coordinating)
                # Square planar typically has angles around 90° and 180°
                # Tetrahedral has angles around 109.5°
                if angles:
                    avg_angle = float(np.mean(angles))
                    if 85 <= avg_angle <= 95:
                        return "square_planar"
                    elif 105 <= avg_angle <= 115:
                        return "tetrahedral"
            return "tetrahedral"  # Default
        elif coord_num == 5:
            return "trigonal_bipyramidal"
        elif coord_num == 6:
            # Check angles to distinguish octahedral from trigonal prismatic
            if len(coordinating) >= 6:
                angles = self._calculate_coordination_angles(metal, coordinating)
                if angles:
                    avg_angle = float(np.mean(angles))
                    if 85 <= avg_angle <= 95:
                        return "octahedral"
            return "octahedral"  # Default
        elif coord_num == 7:
            return "pentagonal_bipyramidal"
        elif coord_num == 8:
            return "square_antiprismatic"
        else:
            return "unknown"
    
    def _calculate_coordination_angles(
        self,
        metal: Dict,
        coordinating: List[Dict]
    ) -> List[float]:
        """
        Calculate angles between coordinating atoms.
        
        Args:
            metal: Metal atom dictionary
            coordinating: List of coordinating residue dictionaries
            
        Returns:
            List of angles in degrees
        """
        if len(coordinating) < 3:
            return []
        
        metal_coords = np.array([metal["x"], metal["y"], metal["z"]])
        angles = []
        
        # Calculate angles between pairs of coordinating atoms
        for i in range(len(coordinating)):
            for j in range(i + 1, len(coordinating)):
                coord1 = np.array(coordinating[i].get("coordinates", [0, 0, 0]))
                coord2 = np.array(coordinating[j].get("coordinates", [0, 0, 0]))
                
                if np.all(coord1 == 0) or np.all(coord2 == 0):
                    continue
                
                v1 = coord1 - metal_coords
                v2 = coord2 - metal_coords
                
                cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                cos_angle = np.clip(cos_angle, -1, 1)
                angle = float(np.degrees(np.arccos(cos_angle)))
                angles.append(angle)
        
        return angles
    
    def _determine_oxidation_state(
        self,
        metal_type: str,
        coordination_number: int,
        geometry: str
    ) -> Optional[str]:
        """
        Determine oxidation state based on metal type and coordination.
        
        Args:
            metal_type: Metal type code (e.g., "FE", "CU", "ZN")
            coordination_number: Number of coordinating atoms
            geometry: Coordination geometry
            
        Returns:
            Oxidation state string or None
        """
        # Common oxidation states for different metals
        oxidation_states = {
            "FE": {
                4: "Fe(II)",  # Common for heme
                5: "Fe(II)",
                6: "Fe(II)" if geometry == "octahedral" else "Fe(III)"
            },
            "FE2": "Fe(II)",
            "FE3": "Fe(III)",
            "CU": {
                2: "Cu(I)",
                4: "Cu(II)",
                5: "Cu(II)",
                6: "Cu(II)"
            },
            "CU1": "Cu(I)",
            "CU2": "Cu(II)",
            "ZN": "Zn(II)",  # Always +2
            "ZN2": "Zn(II)",
            "MN": "Mn(II)",
            "MN2": "Mn(II)",
            "CO": "Co(II)",
            "CO2": "Co(II)",
            "NI": "Ni(II)",
            "NI2": "Ni(II)",
            "CA": "Ca(II)",
            "MG": "Mg(II)",
            "K": "K(I)",
            "NA": "Na(I)"
        }
        
        if metal_type in oxidation_states:
            state_info = oxidation_states[metal_type]
            if isinstance(state_info, dict):
                return state_info.get(coordination_number, None)
            else:
                return state_info
        
        return None
    
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

