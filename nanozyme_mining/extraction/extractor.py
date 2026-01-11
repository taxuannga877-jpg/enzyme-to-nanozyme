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
from typing import Dict, List, Optional, Tuple, Any

from .motif import CatalyticMotif, AnchorAtom, GeometryConstraint

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


# Catalytic residue definitions for common enzyme types
CATALYTIC_RESIDUES = {
    "POD": ["HIS", "ARG", "ASN"],      # Peroxidase
    "Peroxidase": ["HIS", "ARG", "ASN"],
    "SOD": ["HIS", "ASP", "CYS"],      # Superoxide dismutase
    "Superoxide Dismutase": ["HIS", "ASP", "CYS"],
    "CAT": ["HIS", "ASN", "TYR"],      # Catalase
    "Catalase": ["HIS", "ASN", "TYR"],
    "GSH": ["SEC", "CYS", "GLN"],      # Glutathione peroxidase
    "Glutathione Peroxidase": ["SEC", "CYS", "GLN"],
    "OXD": ["HIS", "CYS", "TYR"],      # Oxidase
    "Oxidase": ["HIS", "CYS", "TYR"],
    "LAC": ["HIS", "CYS", "ASP"],      # Laccase
    "Laccase": ["HIS", "CYS", "ASP"],
}

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
        nanozyme_type: str,
        active_site_indices: Optional[List[int]] = None
    ) -> List[Dict]:
        """
        Find catalytic residues in the structure.

        Args:
            atoms: List of atom dictionaries
            nanozyme_type: Type of nanozyme
            active_site_indices: Known active site residue indices

        Returns:
            List of catalytic residue atoms
        """
        target_residues = CATALYTIC_RESIDUES.get(nanozyme_type, [])
        catalytic_atoms = []

        for atom in atoms:
            res_name = atom["residue_name"]
            res_num = atom["residue_number"]

            # Check if residue is catalytic type
            if res_name not in target_residues:
                continue

            # Check if in known active sites
            if active_site_indices and res_num not in active_site_indices:
                continue

            # Check if atom is a donor atom
            if res_name in DONOR_ATOMS:
                if atom["atom_name"] in DONOR_ATOMS[res_name]:
                    catalytic_atoms.append(atom)

        return catalytic_atoms

    def extract_motif(
        self,
        pdb_path: str,
        uniprot_id: str,
        ec_number: str,
        nanozyme_type: str,
        active_site_indices: Optional[List[int]] = None,
        functional_roles: Optional[Dict[Tuple[str, int], str]] = None,
        mcsa_query: Optional[Any] = None
    ) -> Optional[CatalyticMotif]:
        """
        Extract catalytic motif from PDB structure with full residue information.

        Args:
            pdb_path: Path to PDB file
            uniprot_id: UniProt ID
            ec_number: EC number
            nanozyme_type: Nanozyme type
            active_site_indices: Known active site residue indices
            functional_roles: Dict mapping (residue_name, residue_number) to functional role
            mcsa_query: M-CSA query instance for getting functional roles

        Returns:
            CatalyticMotif object or None
        """
        atoms = self.parse_pdb(pdb_path)
        if not atoms:
            return None

        # Get functional roles from M-CSA if available
        if mcsa_query and not functional_roles:
            functional_roles = self._get_functional_roles_from_mcsa(mcsa_query, ec_number, atoms)

        catalytic_atoms = self.find_catalytic_residues(
            atoms, nanozyme_type, active_site_indices
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
            extraction_method="rule_based_enhanced"
        )

        # Add residue structure information
        motif.residue_structures = residue_structures

        # Generate 2D structure if RDKit is available
        if RDKIT_AVAILABLE:
            motif.structure_2d_svg = self._generate_2d_structure(motif, residue_structures)

        return motif

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

    def _get_functional_roles_from_mcsa(
        self,
        mcsa_query: Any,
        ec_number: str,
        atoms: List[Dict]
    ) -> Dict[Tuple[str, int], str]:
        """
        Get functional roles from M-CSA database.

        Args:
            mcsa_query: M-CSA query instance
            ec_number: EC number
            atoms: List of atoms to match

        Returns:
            Dictionary mapping (residue_name, residue_number) to functional role
        """
        functional_roles = {}

        try:
            # Query M-CSA for catalytic residues
            metal_sites = mcsa_query.get_metal_sites(ec_number)
            catalytic_residues = metal_sites.get("catalytic_residues", [])

            # Create mapping from M-CSA residues
            for residue in catalytic_residues:
                res_name = residue.get("residue_type", "")
                res_num = residue.get("residue_number")
                roles = residue.get("roles", [])
                roles_summary = residue.get("roles_summary", "")
                
                if res_name and res_num:
                    key = (res_name, res_num)
                    # Combine roles into a single string
                    if roles_summary:
                        functional_roles[key] = roles_summary
                    elif roles:
                        functional_roles[key] = ", ".join(roles)
        except Exception as e:
            print(f"⚠️  Error getting functional roles from M-CSA: {e}")

        return functional_roles

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
            opts.backgroundColour = (255, 255, 255, 0)  # Transparent background
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
                
                # Add label
                if idx < len(labels):
                    drawer.SetFontSize(12)
                    drawer.DrawAnnotation((x_offset + 150, y_offset + 180), labels[idx])

            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            # Clean up SVG
            svg = svg.replace('white', 'none')
            
            return svg

        except Exception as e:
            print(f"⚠️  Error generating 2D structure: {e}")
            return ""
