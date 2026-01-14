"""
Motif Extractor - Stage 2 Core Module
======================================

Extracts catalytic motifs from enzyme structures.
Based on ChemEnzyRetroPlanner's active site extraction patterns.
"""

import os
import json
import math
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .motif import CatalyticMotif, AnchorAtom, GeometryConstraint
from ..structure.pdb_parser import PDBParser as ComprehensivePDBParser
from ..structure.pdb_metal_extractor import PDBMetalExtractor
from ..structure.environment_analyzer import EnvironmentAnalyzer

# Try to import RDKit for 2D structure generation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("⚠️  RDKit not available. 2D structure generation will be disabled.")

# Try to import BioPython for residue structure extraction
try:
    from Bio.PDB import PDBParser, Select
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("⚠️  BioPython not available. Full residue structure extraction will be limited.")


# Note: Catalytic residues are now extracted from PDB SITE records,
# not hardcoded. This ensures we extract what's actually in the PDB file.

# Key donor atoms for each residue type
DONOR_ATOMS = {
    "HIS": ["NE2", "ND1"],
    "SER": ["OG"],
    "CYS": ["SG"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "TYR": ["OH"],
    "ARG": ["NH1", "NH2"],
    "LYS": ["NZ"],
    "ASN": ["OD1", "ND2"],
    "GLN": ["OE1", "NE2"],
    "SEC": ["SE"],
}


class MotifExtractor:
    """
    Extracts catalytic motifs from PDB structures.

    Based on ChemEnzyRetroPlanner's MyProtein and active site
    extraction patterns.
    """

    def __init__(self, output_dir: str = "./motifs"):
        """
        Initialize extractor.

        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize PDB parser if BioPython is available
        if BIOPYTHON_AVAILABLE:
            self.pdb_parser = PDBParser(QUIET=True)
        
        # Initialize comprehensive PDB parser and analyzers
        self.comprehensive_parser = ComprehensivePDBParser()
        self.metal_extractor = PDBMetalExtractor()
        self.env_analyzer = EnvironmentAnalyzer()

    def parse_pdb(self, pdb_path: str) -> List[Dict]:
        """
        Parse PDB file and extract atom information.

        Args:
            pdb_path: Path to PDB file

        Returns:
            List of atom dictionaries
        """
        atoms = []

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom = {
                        "atom_name": line[12:16].strip(),
                        "residue_name": line[17:20].strip(),
                        "chain_id": line[21].strip(),
                        "residue_number": int(line[22:26].strip()),
                        "x": float(line[30:38].strip()),
                        "y": float(line[38:46].strip()),
                        "z": float(line[46:54].strip()),
                        "element": line[76:78].strip() if len(line) > 76 else ""
                    }
                    atoms.append(atom)

        return atoms

    def calculate_distance(self, coord1: List[float], coord2: List[float]) -> float:
        """Calculate Euclidean distance between two points."""
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

    def calculate_angle(self, c1: List[float], c2: List[float], c3: List[float]) -> float:
        """Calculate angle (degrees) between three points."""
        v1 = np.array(c1) - np.array(c2)
        v2 = np.array(c3) - np.array(c2)
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return math.degrees(math.acos(np.clip(cos_angle, -1, 1)))

    def find_catalytic_residues(
        self,
        atoms: List[Dict],
        pdb_path: str,
        active_site_indices: Optional[List[int]] = None
    ) -> Tuple[List[Dict], str]:
        """
        Find catalytic residues from PDB file.

        Priority:
        1. PDB SITE records (highest priority)
        2. active_site_indices parameter (from JSON file)

        Args:
            atoms: List of atom dictionaries
            pdb_path: Path to PDB file (to extract SITE records)
            active_site_indices: Known active site residue indices (optional)

        Returns:
            Tuple of (catalytic_atoms: List[Dict], source: str)
            - catalytic_atoms: List of atom dictionaries in unified format
            - source: "pdb_site" or "none"
        """
        # 1. Try PDB SITE records first (highest priority)
        site_residues = self._extract_site_residues(pdb_path)
        
        # If no SITE records found, fall back to active_site_indices if provided
        has_explicit_active_sites = active_site_indices is not None and len(active_site_indices) > 0
        if not site_residues and has_explicit_active_sites:
            site_residues = {(None, idx) for idx in active_site_indices}
        
        if site_residues:
            # Find atoms matching SITE residues
            catalytic_atoms = []
            for atom in atoms:
                res_name = atom["residue_name"]
                res_num = atom["residue_number"]
                chain_id = atom["chain_id"]
                
                # Check if this residue is in SITE records
                is_catalytic = False
                for site_chain, site_num in site_residues:
                    if site_chain is None:
                        # Match by residue number only
                        if res_num == site_num:
                            is_catalytic = True
                            break
                    else:
                        # Match by chain and residue number
                        if chain_id == site_chain and res_num == site_num:
                            is_catalytic = True
                            break
                
                if not is_catalytic:
                    continue
                
                # Check if atom is a donor atom (catalytic atoms are usually donor/acceptor atoms)
                if res_name in DONOR_ATOMS:
                    if atom["atom_name"] in DONOR_ATOMS[res_name]:
                        catalytic_atoms.append(atom)
                else:
                    # If residue type not in DONOR_ATOMS, still include key atoms
                    catalytic_atoms.append(atom)
            
            # If we have explicit active_site information, try to extract atoms even if no match found
            if has_explicit_active_sites:
                if catalytic_atoms:
                    return catalytic_atoms, "pdb_site"
                else:
                    # Try to extract atoms from active_site_indices
                    # This handles cases where residue numbering doesn't match
                    fallback_atoms = []
                    for idx in active_site_indices:
                        for atom in atoms:
                            if atom["residue_number"] == idx:
                                # Prefer key atoms (CA, or donor atoms)
                                if atom["atom_name"] == "CA":
                                    fallback_atoms.append(atom)
                                    break
                                elif atom["residue_name"] in DONOR_ATOMS:
                                    if atom["atom_name"] in DONOR_ATOMS[atom["residue_name"]]:
                                        fallback_atoms.append(atom)
                                        break
                    
                    # If still no atoms found, at least include CA atoms
                    if not fallback_atoms:
                        for idx in active_site_indices:
                            for atom in atoms:
                                if atom["residue_number"] == idx and atom["atom_name"] == "CA":
                                    fallback_atoms.append(atom)
                                    break
                    
                    return fallback_atoms, "pdb_site" if fallback_atoms else "none"
            
            # If we found atoms from PDB SITE records, return them
            if catalytic_atoms:
                return catalytic_atoms, "pdb_site"
        
        # No catalytic residues found
        return [], "none"
    
    def _extract_site_residues(self, pdb_path: str) -> set:
        """
        Extract catalytic residues from PDB SITE records.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            Set of (chain_id, residue_number) tuples from SITE records
        """
        site_residues = set()
        
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith("SITE"):
                        # Parse SITE record
                        # Format: SITE    1 AC1 4 HIS A 146  HIS A  57  HIS A  87  HIS A 119
                        try:
                            # SITE can have up to 4 residues per line
                            for i in range(4):
                                start = 18 + i * 11
                                if len(line) > start + 10:
                                    res_name = line[start:start+3].strip()
                                    chain = line[start+4:start+5].strip() if len(line) > start+4 else ""
                                    res_num_str = line[start+5:start+10].strip()
                                    if res_name and res_num_str:
                                        try:
                                            res_num = int(res_num_str)
                                            chain_id = chain if chain else None
                                            site_residues.add((chain_id, res_num))
                                        except ValueError:
                                            pass
                        except (ValueError, IndexError):
                            continue
        except Exception as e:
            print(f"⚠️  Warning: Error reading SITE records from {pdb_path}: {e}")
        
        return site_residues


    def extract_motif(
        self,
        pdb_path: str,
        uniprot_id: str,
        ec_number: str,
        nanozyme_type: str,
        active_site_indices: Optional[List[int]] = None,
        functional_roles: Optional[Dict[Tuple[str, int], str]] = None
    ) -> Optional[CatalyticMotif]:
        """
        Extract catalytic motif from PDB structure with full residue information.

        Args:
            pdb_path: Path to PDB file
            uniprot_id: UniProt ID
            ec_number: EC number
            nanozyme_type: Nanozyme type
            active_site_indices: Known active site residue indices
            functional_roles: Dict mapping (residue_name, residue_number) to functional role (optional)

        Returns:
            CatalyticMotif object or None
        """
        atoms = self.parse_pdb(pdb_path)
        if not atoms:
            return None

        # Extract catalytic residues from PDB file
        catalytic_atoms, extraction_source = self.find_catalytic_residues(
            atoms, pdb_path, active_site_indices
        )

        if not catalytic_atoms:
            return None

        # Extract full residue structures
        residue_structures = self._extract_residue_structures(
            pdb_path, catalytic_atoms
        )

        # Create anchor atoms with functional roles
        anchor_atoms = []
        for atom in catalytic_atoms:
            # Get functional role
            role = ""
            if functional_roles:
                key = (atom["residue_name"], atom["residue_number"])
                role = functional_roles.get(key, "")
            
            # If no role from functional_roles, try to infer from residue type
            if not role:
                role = self._infer_functional_role(atom["residue_name"], atom["atom_name"])

            anchor = AnchorAtom(
                atom_name=atom["atom_name"],
                residue_name=atom["residue_name"],
                residue_number=atom["residue_number"],
                chain_id=atom["chain_id"],
                element=atom["element"],
                coordinates=[atom["x"], atom["y"], atom["z"]],
                is_donor=True,
                role=role
            )
            anchor_atoms.append(anchor)

        # Calculate geometry constraints
        geometry = self._calculate_geometry(anchor_atoms)

        motif_id = f"{uniprot_id}_{ec_number}_{nanozyme_type}"

        motif = CatalyticMotif(
            motif_id=motif_id,
            source_uniprot_id=uniprot_id,
            source_ec_number=ec_number,
            nanozyme_type=nanozyme_type,
            anchor_atoms=anchor_atoms,
            geometry_constraints=geometry,
            extraction_method=extraction_source,
            reaction_smiles=""  # No longer used, but kept for compatibility
        )

        # Add residue structure information
        motif.residue_structures = residue_structures

        # Generate 2D structure if RDKit is available
        if RDKIT_AVAILABLE:
            motif.structure_2d_svg = self._generate_2d_structure(motif, residue_structures)
        
        # Extract extended PDB information
        try:
            self._extract_extended_pdb_info(motif, pdb_path, catalytic_atoms)
        except Exception as e:
            print(f"⚠️  Warning: Error extracting extended PDB information: {e}")

        return motif
    
    def _extract_extended_pdb_info(
        self,
        motif: CatalyticMotif,
        pdb_path: str,
        catalytic_atoms: List[Dict]
    ) -> None:
        """
        Extract extended PDB information and add to motif.
        
        Args:
            motif: CatalyticMotif object to update
            pdb_path: Path to PDB file
            catalytic_atoms: List of catalytic atom dictionaries
        """
        pdb_file = Path(pdb_path)
        
        # 1. Parse comprehensive PDB information
        parsed_data = self.comprehensive_parser.parse_pdb_file(pdb_file)
        
        # Check if any data was calculated (fallback) rather than parsed
        calculated_flags = []
        if parsed_data.get("_calculated_ssbonds"):
            calculated_flags.append("disulfide bonds")
        if parsed_data.get("_calculated_ligands"):
            calculated_flags.append("ligands/cofactors")
        
        # 2. Extract metal sites first (needed for calculating metal coordination links)
        metal_sites = []
        try:
            metal_sites = self.metal_extractor.parse_pdb_file(pdb_file)
            motif.metal_sites = [site.to_dict() for site in metal_sites]
        except Exception as e:
            print(f"⚠️  Warning: Error extracting metal sites: {e}")
            motif.metal_sites = []
        
        # 3. Extract chemical bonds (with fallback calculation support)
        motif.chemical_bonds = self.comprehensive_parser.extract_chemical_bonds(
            parsed_data, pdb_file=pdb_file, metal_sites=[site.to_dict() for site in metal_sites]
        )
        
        # Check if links were calculated
        if motif.chemical_bonds.get("links") and any(
            link.get("calculated_from_coords", False) 
            for link in motif.chemical_bonds.get("links", [])
        ):
            calculated_flags.append("metal coordination links")
        
        # 4. Extract ligands and cofactors
        motif.ligands_and_cofactors = self.comprehensive_parser.extract_ligands_and_cofactors(
            parsed_data, parsed_data.get("het_info")
        )
        
        # 5. Extract site annotations
        motif.site_annotations = parsed_data.get("sites", [])
        
        # 6. Extract secondary structure
        motif.secondary_structure = self.comprehensive_parser.extract_secondary_structure(parsed_data)
        
        # 7. Print warnings if data was calculated rather than parsed
        if calculated_flags:
            print(f"ℹ️  Note: The following information was calculated from coordinates "
                  f"(not parsed from PDB records): {', '.join(calculated_flags)}. "
                  f"This is common for AlphaFold-predicted structures.")
        
        # 8. Analyze residue environment for catalytic residues
        try:
            target_residues = [
                (atom["residue_name"], atom["residue_number"])
                for atom in catalytic_atoms
            ]
            # Remove duplicates
            target_residues = list(set(target_residues))
            
            residue_env = self.env_analyzer.analyze_residue_environment(
                pdb_file, target_residues
            )
            motif.residue_environment = residue_env
        except Exception as e:
            print(f"⚠️  Warning: Error analyzing residue environment: {e}")
            motif.residue_environment = {}

    def _calculate_geometry(self, anchor_atoms: List[AnchorAtom]) -> List[GeometryConstraint]:
        """Calculate geometry constraints between anchor atoms."""
        constraints = []

        # Calculate pairwise distances
        for i in range(len(anchor_atoms)):
            for j in range(i + 1, len(anchor_atoms)):
                dist = self.calculate_distance(
                    anchor_atoms[i].coordinates,
                    anchor_atoms[j].coordinates
                )
                constraints.append(GeometryConstraint(
                    constraint_type="distance",
                    atom_indices=[i, j],
                    value=dist,
                    unit="angstrom"
                ))

        # Calculate angles for triplets
        if len(anchor_atoms) >= 3:
            for i in range(len(anchor_atoms) - 2):
                angle = self.calculate_angle(
                    anchor_atoms[i].coordinates,
                    anchor_atoms[i + 1].coordinates,
                    anchor_atoms[i + 2].coordinates
                )
                constraints.append(GeometryConstraint(
                    constraint_type="angle",
                    atom_indices=[i, i + 1, i + 2],
                    value=angle,
                    unit="degree"
                ))

        return constraints

    def _extract_residue_structures(
        self,
        pdb_path: str,
        catalytic_atoms: List[Dict]
    ) -> Dict[Tuple[str, int], Dict]:
        """
        Extract full residue structures from PDB file.

        Args:
            pdb_path: Path to PDB file
            catalytic_atoms: List of catalytic atom dictionaries

        Returns:
            Dictionary mapping (residue_name, residue_number) to residue structure
        """
        residue_structures = {}

        if BIOPYTHON_AVAILABLE:
            try:
                structure = self.pdb_parser.get_structure('motif', pdb_path)
                
                # Get unique residues from catalytic atoms
                unique_residues = {}
                for atom in catalytic_atoms:
                    key = (atom["residue_name"], atom["residue_number"], atom["chain_id"])
                    if key not in unique_residues:
                        unique_residues[key] = atom

                # Extract full residue structures
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            res_name = residue.get_resname()
                            res_num = residue.get_id()[1]
                            chain_id = chain.get_id()
                            
                            key = (res_name, res_num, chain_id)
                            if key in unique_residues:
                                # Extract all atoms in this residue
                                atoms_in_residue = []
                                for atom in residue:
                                    atoms_in_residue.append({
                                        "atom_name": atom.get_name(),
                                        "element": atom.element,
                                        "coordinates": list(atom.get_coord()),
                                        "occupancy": atom.get_occupancy(),
                                        "bfactor": atom.get_bfactor()
                                    })
                                
                                residue_key = (res_name, res_num)
                                residue_structures[residue_key] = {
                                    "residue_name": res_name,
                                    "residue_number": res_num,
                                    "chain_id": chain_id,
                                    "atoms": atoms_in_residue,
                                    "num_atoms": len(atoms_in_residue)
                                }
            except Exception as e:
                print(f"⚠️  Error extracting residue structures with BioPython: {e}")
                # Fallback to simple extraction from PDB file
                residue_structures = self._extract_residue_structures_simple(pdb_path, catalytic_atoms)
        else:
            # Fallback to simple extraction
            residue_structures = self._extract_residue_structures_simple(pdb_path, catalytic_atoms)

        return residue_structures

    def _extract_residue_structures_simple(
        self,
        pdb_path: str,
        catalytic_atoms: List[Dict]
    ) -> Dict[Tuple[str, int], Dict]:
        """Simple residue structure extraction from PDB file without BioPython."""
        residue_structures = {}
        
        # Get unique residues
        unique_residues = {}
        for atom in catalytic_atoms:
            key = (atom["residue_name"], atom["residue_number"])
            if key not in unique_residues:
                unique_residues[key] = atom

        # Read PDB file and extract atoms for each residue
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    chain_id = line[21].strip()
                    
                    key = (res_name, res_num)
                    if key in unique_residues:
                        if key not in residue_structures:
                            residue_structures[key] = {
                                "residue_name": res_name,
                                "residue_number": res_num,
                                "chain_id": chain_id,
                                "atoms": []
                            }
                        
                        atom_info = {
                            "atom_name": line[12:16].strip(),
                            "element": line[76:78].strip() if len(line) > 76 else "",
                            "coordinates": [
                                float(line[30:38].strip()),
                                float(line[38:46].strip()),
                                float(line[46:54].strip())
                            ],
                            "occupancy": float(line[54:60].strip()) if len(line) > 60 else 1.0,
                            "bfactor": float(line[60:66].strip()) if len(line) > 66 else 0.0
                        }
                        residue_structures[key]["atoms"].append(atom_info)

        # Add atom count
        for key in residue_structures:
            residue_structures[key]["num_atoms"] = len(residue_structures[key]["atoms"])

        return residue_structures


    def _infer_functional_role(self, residue_name: str, atom_name: str) -> str:
        """
        Infer functional role from residue type and atom name.

        Args:
            residue_name: Three-letter residue code
            atom_name: PDB atom name

        Returns:
            Inferred functional role
        """
        # Common functional role patterns
        role_patterns = {
            "HIS": {
                "NE2": "general base / proton acceptor",
                "ND1": "general base / proton acceptor",
                "general": "histidine (general acid/base)"
            },
            "SER": {
                "OG": "nucleophile",
                "general": "serine (nucleophile)"
            },
            "CYS": {
                "SG": "nucleophile / thiol",
                "general": "cysteine (nucleophile)"
            },
            "ASP": {
                "OD1": "general acid / proton donor",
                "OD2": "general acid / proton donor",
                "general": "aspartate (general acid)"
            },
            "GLU": {
                "OE1": "general acid / proton donor",
                "OE2": "general acid / proton donor",
                "general": "glutamate (general acid)"
            },
            "TYR": {
                "OH": "general acid / proton donor",
                "general": "tyrosine (general acid)"
            },
            "ARG": {
                "NH1": "electrostatic stabilizer",
                "NH2": "electrostatic stabilizer",
                "general": "arginine (electrostatic stabilizer)"
            },
            "LYS": {
                "NZ": "general base / nucleophile",
                "general": "lysine (general base)"
            },
            "ASN": {
                "OD1": "hydrogen bond acceptor",
                "ND2": "hydrogen bond donor",
                "general": "asparagine (hydrogen bonding)"
            },
            "GLN": {
                "OE1": "hydrogen bond acceptor",
                "NE2": "hydrogen bond donor",
                "general": "glutamine (hydrogen bonding)"
            },
            "SEC": {
                "SE": "nucleophile / selenocysteine",
                "general": "selenocysteine (nucleophile)"
            }
        }

        if residue_name in role_patterns:
            patterns = role_patterns[residue_name]
            if atom_name in patterns:
                return patterns[atom_name]
            else:
                return patterns.get("general", "")

        return ""

    def _generate_2d_structure(
        self,
        motif: CatalyticMotif,
        residue_structures: Dict[Tuple[str, int], Dict]
    ) -> str:
        """
        Generate 2D chemical structure SVG for the catalytic motif.

        Args:
            motif: CatalyticMotif object
            residue_structures: Dictionary of residue structures

        Returns:
            SVG string of the 2D structure
        """
        if not RDKIT_AVAILABLE:
            return ""

        # Amino acid SMILES mapping (simplified side chain representations)
        AA_SMILES = {
            'HIS': 'CC1=CNC=N1',  # Histidine side chain
            'SER': 'CO',           # Serine side chain
            'CYS': 'CS',           # Cysteine side chain
            'ASP': 'CC(=O)O',      # Aspartate side chain
            'GLU': 'CCC(=O)O',     # Glutamate side chain
            'TYR': 'CC1=CC=C(C=C1)O',  # Tyrosine side chain
            'ARG': 'CCCCNC(=N)N',  # Arginine side chain
            'LYS': 'CCCCN',        # Lysine side chain
            'ASN': 'CC(=O)N',      # Asparagine side chain
            'GLN': 'CCC(=O)N',     # Glutamine side chain
            'SEC': 'C[SeH]',       # Selenocysteine side chain
        }

        try:
            # Collect unique residues
            unique_residues = {}
            for atom in motif.anchor_atoms:
                key = (atom.residue_name, atom.residue_number)
                if key not in unique_residues:
                    unique_residues[key] = atom

            # Generate SMILES for each residue side chain
            mols = []
            labels = []
            
            for (res_name, res_num), anchor in unique_residues.items():
                if res_name in AA_SMILES:
                    smiles = AA_SMILES[res_name]
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        # Highlight the key atom
                        highlight_atoms = []
                        if res_name == "HIS" and anchor.atom_name in ["NE2", "ND1"]:
                            # Find N atoms in histidine
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "N":
                                    highlight_atoms.append(atom.GetIdx())
                        elif res_name == "SER" and anchor.atom_name == "OG":
                            # Find O atom in serine
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "O":
                                    highlight_atoms.append(atom.GetIdx())
                        elif res_name == "CYS" and anchor.atom_name == "SG":
                            # Find S atom in cysteine
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "S":
                                    highlight_atoms.append(atom.GetIdx())
                        # Add more patterns as needed...
                        
                        mols.append(mol)
                        role = anchor.role if anchor.role else f"{res_name}{res_num}"
                        labels.append(f"{res_name}{res_num}\n({role})")

            if not mols:
                return ""

            # Generate combined SVG
            from rdkit.Chem import rdDepictor
            
            # Calculate grid layout
            n_mols = len(mols)
            cols = min(3, n_mols)
            rows = (n_mols + cols - 1) // cols
            
            width = 300 * cols
            height = 200 * rows
            
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
            opts = drawer.drawOptions()
            opts.backgroundColour = (1.0, 1.0, 1.0)  # White background (RGB 0-1 range, no alpha)
            opts.addStereoAnnotation = True
            opts.annotationFontScale = 0.8
            
            # Draw each molecule in grid
            for idx, mol in enumerate(mols):
                row = idx // cols
                col = idx % cols
                
                # Compute 2D coordinates
                rdDepictor.Compute2DCoords(mol)
                
                # Calculate position
                x_offset = col * 300
                y_offset = row * 200
                
                # Draw molecule
                drawer.SetOffset(x_offset, y_offset)
                drawer.DrawMolecule(mol)

            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            # Add labels directly to SVG (more compatible than DrawAnnotation)
            if labels:
                import re
                # Find the closing </svg> tag
                svg_end = svg.rfind('</svg>')
                if svg_end != -1:
                    label_elements = []
                    for idx, label in enumerate(labels):
                        row = idx // cols
                        col = idx % cols
                        x_pos = col * 300 + 150
                        y_pos = row * 200 + 180
                        # Split label by newline for multi-line text
                        label_lines = label.split('\n')
                        for line_idx, line in enumerate(label_lines):
                            y_offset_line = y_pos + line_idx * 15
                            label_elements.append(
                                f'<text x="{x_pos}" y="{y_offset_line}" '
                                f'font-family="Arial" font-size="12" '
                                f'text-anchor="middle" fill="black">{line}</text>'
                            )
                    # Insert labels before closing </svg> tag
                    svg = svg[:svg_end] + '\n'.join(label_elements) + '\n' + svg[svg_end:]
            
            # Clean up SVG
            svg = svg.replace('white', 'none')
            
            return svg

        except Exception as e:
            print(f"⚠️  Error generating 2D structure: {e}")
            return ""
