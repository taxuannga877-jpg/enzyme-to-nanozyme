"""
Environment Analyzer - Calculate Residue Environment Properties
==============================================================

Calculates various environment properties for residues:
- Solvent Accessible Surface Area (SASA)
- Residue depth (distance from surface)
- Hydrophobicity classification
- Interaction network (hydrogen bonds, salt bridges, hydrophobic interactions)
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict


# Hydrophobicity scale (Kyte-Doolittle)
HYDROPHOBICITY_SCALE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5,
    "CYS": 2.5, "GLN": -3.5, "GLU": -3.5, "GLY": -0.4,
    "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9,
    "MET": 1.9, "PHE": 2.8, "PRO": -1.6, "SER": -0.8,
    "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
    "SEC": 2.5  # Selenocysteine
}

# Hydrogen bond donor/acceptor atoms
H_BOND_DONORS = {
    "HIS": ["NE2", "ND1"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "CYS": ["SG"],
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2", "NE"],
    "ASN": ["ND2"],
    "GLN": ["NE2"],
    "TRP": ["NE1"]
}

H_BOND_ACCEPTORS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "ASN": ["OD1"],
    "GLN": ["OE1"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "HIS": ["NE2", "ND1"],
    "MET": ["SD"]
}

# Charged residues for salt bridges
POSITIVE_CHARGED = {"LYS", "ARG", "HIS"}
NEGATIVE_CHARGED = {"ASP", "GLU"}

# Hydrophobic residues
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}


class EnvironmentAnalyzer:
    """Analyze residue environment properties."""
    
    def __init__(self, probe_radius: float = 1.4):
        """
        Initialize the environment analyzer.
        
        Args:
            probe_radius: Radius of water probe for SASA calculation (default: 1.4 Å)
        """
        self.probe_radius = probe_radius
    
    def analyze_residue_environment(
        self,
        pdb_path: Path,
        target_residues: Optional[List[Tuple[str, int]]] = None
    ) -> Dict[Tuple[str, int], Dict]:
        """
        Analyze environment for target residues.
        
        Args:
            pdb_path: Path to PDB file
            target_residues: List of (residue_name, residue_number) tuples to analyze.
                           If None, analyzes all residues.
        
        Returns:
            Dictionary mapping (res_name, res_num) to environment properties
        """
        # Parse PDB file to get all atoms
        atoms = self._parse_atoms(pdb_path)
        
        # Group atoms by residue
        residues = self._group_atoms_by_residue(atoms)
        
        # Determine which residues to analyze
        if target_residues is None:
            target_residues = list(residues.keys())
        
        results = {}
        
        for res_key in target_residues:
            if res_key not in residues:
                continue
            
            res_atoms = residues[res_key]
            res_name = res_key[0]
            
            # Calculate properties
            sasa = self._calculate_sasa_simple(res_atoms, residues)
            depth = self._calculate_residue_depth(res_atoms, residues)
            hydrophobicity = self._classify_hydrophobicity(res_name)
            interactions = self._build_interaction_network(
                res_key, res_atoms, residues
            )
            
            results[res_key] = {
                "residue_name": res_name,
                "residue_number": res_key[1],
                "sasa": round(sasa, 2),
                "depth": round(depth, 2),
                "hydrophobicity": hydrophobicity,
                "hydrophobicity_score": HYDROPHOBICITY_SCALE.get(res_name, 0.0),
                "interactions": interactions
            }
        
        return results
    
    def _parse_atoms(self, pdb_path: Path) -> List[Dict]:
        """Parse atoms from PDB file."""
        atoms = []
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        atom = {
                            "atom_name": line[12:16].strip(),
                            "residue_name": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "residue_number": int(line[22:26].strip()),
                            "x": float(line[30:38].strip()),
                            "y": float(line[38:46].strip()),
                            "z": float(line[46:54].strip()),
                            "element": line[76:78].strip() if len(line) > 76 else ""
                        }
                        atoms.append(atom)
                    except (ValueError, IndexError):
                        continue
        
        return atoms
    
    def _group_atoms_by_residue(
        self,
        atoms: List[Dict]
    ) -> Dict[Tuple[str, int], List[Dict]]:
        """Group atoms by residue."""
        residues = defaultdict(list)
        
        for atom in atoms:
            key = (atom["residue_name"], atom["residue_number"])
            residues[key].append(atom)
        
        return dict(residues)
    
    def _calculate_sasa_simple(
        self,
        res_atoms: List[Dict],
        all_residues: Dict[Tuple[str, int], List[Dict]]
    ) -> float:
        """
        Calculate simplified SASA (Solvent Accessible Surface Area).
        
        This is a simplified version that estimates SASA based on
        atom exposure to solvent.
        
        Args:
            res_atoms: Atoms in the residue
            all_residues: All residues in the structure
        
        Returns:
            Estimated SASA in Å²
        """
        if not res_atoms:
            return 0.0
        
        # Get coordinates of all atoms in residue
        coords = np.array([
            [atom["x"], atom["y"], atom["z"]] for atom in res_atoms
        ])
        
        # Get coordinates of all other atoms
        other_coords = []
        for res_key, atoms in all_residues.items():
            if atoms != res_atoms:  # Exclude current residue
                for atom in atoms:
                    other_coords.append([atom["x"], atom["y"], atom["z"]])
        
        if not other_coords:
            # If no other atoms, all atoms are exposed
            # Rough estimate: 4πr² per atom
            return len(res_atoms) * 50.0  # Rough estimate
        
        other_coords = np.array(other_coords)
        
        # Simple SASA estimation: count atoms that are not buried
        # An atom is considered exposed if it's not within 3.5 Å of many other atoms
        exposed_count = 0
        for coord in coords:
            distances = np.linalg.norm(other_coords - coord, axis=1)
            # Count nearby atoms (within 5 Å)
            nearby = np.sum(distances < 5.0)
            # If fewer than 5 nearby atoms, consider it exposed
            if nearby < 5:
                exposed_count += 1
        
        # Rough SASA estimate: exposed atoms contribute ~50 Å²
        sasa = exposed_count * 50.0
        
        return sasa
    
    def _calculate_residue_depth(
        self,
        res_atoms: List[Dict],
        all_residues: Dict[Tuple[str, int], List[Dict]]
    ) -> float:
        """
        Calculate residue depth (distance from protein surface).
        
        Args:
            res_atoms: Atoms in the residue
            all_residues: All residues in the structure
        
        Returns:
            Average depth in Å
        """
        if not res_atoms:
            return 0.0
        
        # Get center of residue
        coords = np.array([
            [atom["x"], atom["y"], atom["z"]] for atom in res_atoms
        ])
        center = np.mean(coords, axis=0)
        
        # Get all other atom coordinates
        other_coords = []
        for res_key, atoms in all_residues.items():
            if atoms != res_atoms:
                for atom in atoms:
                    other_coords.append([atom["x"], atom["y"], atom["z"]])
        
        if not other_coords:
            return 0.0  # Surface residue
        
        other_coords = np.array(other_coords)
        
        # Find nearest neighbor distance (approximate surface distance)
        distances = np.linalg.norm(other_coords - center, axis=1)
        min_distance = float(np.min(distances))
        
        # Depth is roughly the distance to nearest neighbor
        # For buried residues, this will be small
        # For surface residues, this will be larger
        depth = max(0.0, 5.0 - min_distance)  # Invert: larger min_distance = smaller depth
        
        return float(depth)
    
    def _classify_hydrophobicity(self, residue_name: str) -> str:
        """
        Classify residue hydrophobicity.
        
        Args:
            residue_name: Three-letter residue code
        
        Returns:
            "hydrophobic", "hydrophilic", or "neutral"
        """
        score = HYDROPHOBICITY_SCALE.get(residue_name, 0.0)
        
        if score > 1.0:
            return "hydrophobic"
        elif score < -1.0:
            return "hydrophilic"
        else:
            return "neutral"
    
    def _build_interaction_network(
        self,
        res_key: Tuple[str, int],
        res_atoms: List[Dict],
        all_residues: Dict[Tuple[str, int], List[Dict]]
    ) -> Dict[str, List[Dict]]:
        """
        Build interaction network for a residue.
        
        Args:
            res_key: (residue_name, residue_number) tuple
            res_atoms: Atoms in the residue
            all_residues: All residues in the structure
        
        Returns:
            Dictionary with interaction types as keys
        """
        interactions = {
            "hydrogen_bonds": [],
            "salt_bridges": [],
            "hydrophobic_interactions": []
        }
        
        res_name = res_key[0]
        
        # Check interactions with other residues
        for other_key, other_atoms in all_residues.items():
            if other_key == res_key:
                continue
            
            other_name = other_key[0]
            
            # Check hydrogen bonds
            h_bonds = self._find_hydrogen_bonds(
                res_name, res_atoms, other_name, other_atoms
            )
            if h_bonds:
                interactions["hydrogen_bonds"].extend(h_bonds)
            
            # Check salt bridges
            salt_bridge = self._find_salt_bridge(
                res_name, res_atoms, other_name, other_atoms
            )
            if salt_bridge:
                interactions["salt_bridges"].append(salt_bridge)
            
            # Check hydrophobic interactions
            if res_name in HYDROPHOBIC_RESIDUES and other_name in HYDROPHOBIC_RESIDUES:
                hydrophobic = self._find_hydrophobic_interaction(
                    res_atoms, other_atoms
                )
                if hydrophobic:
                    interactions["hydrophobic_interactions"].append(hydrophobic)
        
        return interactions
    
    def _find_hydrogen_bonds(
        self,
        res1_name: str,
        res1_atoms: List[Dict],
        res2_name: str,
        res2_atoms: List[Dict]
    ) -> List[Dict]:
        """Find hydrogen bonds between two residues."""
        h_bonds = []
        
        # Get donor and acceptor atoms
        donors = []
        acceptors = []
        
        if res1_name in H_BOND_DONORS:
            donor_atoms = H_BOND_DONORS[res1_name]
            for atom in res1_atoms:
                if atom["atom_name"] in donor_atoms:
                    donors.append(atom)
        
        if res2_name in H_BOND_ACCEPTORS:
            acceptor_atoms = H_BOND_ACCEPTORS[res2_name]
            for atom in res2_atoms:
                if atom["atom_name"] in acceptor_atoms:
                    acceptors.append(atom)
        
        # Check distances (H-bond: 2.5-3.5 Å)
        for donor in donors:
            donor_coords = np.array([donor["x"], donor["y"], donor["z"]])
            for acceptor in acceptors:
                acceptor_coords = np.array([acceptor["x"], acceptor["y"], acceptor["z"]])
                distance = np.linalg.norm(donor_coords - acceptor_coords)
                
                if 2.5 <= distance <= 3.5:
                    h_bonds.append({
                        "donor_residue": res1_name,
                        "donor_atom": donor["atom_name"],
                        "acceptor_residue": res2_name,
                        "acceptor_atom": acceptor["atom_name"],
                        "distance": round(distance, 2)
                    })
        
        # Also check reverse (res2 as donor, res1 as acceptor)
        donors2 = []
        acceptors2 = []
        
        if res2_name in H_BOND_DONORS:
            donor_atoms = H_BOND_DONORS[res2_name]
            for atom in res2_atoms:
                if atom["atom_name"] in donor_atoms:
                    donors2.append(atom)
        
        if res1_name in H_BOND_ACCEPTORS:
            acceptor_atoms = H_BOND_ACCEPTORS[res1_name]
            for atom in res1_atoms:
                if atom["atom_name"] in acceptor_atoms:
                    acceptors2.append(atom)
        
        for donor in donors2:
            donor_coords = np.array([donor["x"], donor["y"], donor["z"]])
            for acceptor in acceptors2:
                acceptor_coords = np.array([acceptor["x"], acceptor["y"], acceptor["z"]])
                distance = np.linalg.norm(donor_coords - acceptor_coords)
                
                if 2.5 <= distance <= 3.5:
                    h_bonds.append({
                        "donor_residue": res2_name,
                        "donor_atom": donor["atom_name"],
                        "acceptor_residue": res1_name,
                        "acceptor_atom": acceptor["atom_name"],
                        "distance": round(distance, 2)
                    })
        
        return h_bonds
    
    def _find_salt_bridge(
        self,
        res1_name: str,
        res1_atoms: List[Dict],
        res2_name: str,
        res2_atoms: List[Dict]
    ) -> Optional[Dict]:
        """Find salt bridge between two residues."""
        # Salt bridge: positive and negative charged residues within 4-6 Å
        is_pos_neg = (
            (res1_name in POSITIVE_CHARGED and res2_name in NEGATIVE_CHARGED) or
            (res1_name in NEGATIVE_CHARGED and res2_name in POSITIVE_CHARGED)
        )
        
        if not is_pos_neg:
            return None
        
        # Find charged atoms
        res1_charged = []
        res2_charged = []
        
        if res1_name in POSITIVE_CHARGED:
            # Find NZ (LYS), NH1/NH2/NE (ARG), NE2/ND1 (HIS)
            for atom in res1_atoms:
                if res1_name == "LYS" and atom["atom_name"] == "NZ":
                    res1_charged.append(atom)
                elif res1_name == "ARG" and atom["atom_name"] in ["NH1", "NH2", "NE"]:
                    res1_charged.append(atom)
                elif res1_name == "HIS" and atom["atom_name"] in ["NE2", "ND1"]:
                    res1_charged.append(atom)
        elif res1_name in NEGATIVE_CHARGED:
            # Find OD1/OD2 (ASP), OE1/OE2 (GLU)
            for atom in res1_atoms:
                if res1_name == "ASP" and atom["atom_name"] in ["OD1", "OD2"]:
                    res1_charged.append(atom)
                elif res1_name == "GLU" and atom["atom_name"] in ["OE1", "OE2"]:
                    res1_charged.append(atom)
        
        if res2_name in POSITIVE_CHARGED:
            for atom in res2_atoms:
                if res2_name == "LYS" and atom["atom_name"] == "NZ":
                    res2_charged.append(atom)
                elif res2_name == "ARG" and atom["atom_name"] in ["NH1", "NH2", "NE"]:
                    res2_charged.append(atom)
                elif res2_name == "HIS" and atom["atom_name"] in ["NE2", "ND1"]:
                    res2_charged.append(atom)
        elif res2_name in NEGATIVE_CHARGED:
            for atom in res2_atoms:
                if res2_name == "ASP" and atom["atom_name"] in ["OD1", "OD2"]:
                    res2_charged.append(atom)
                elif res2_name == "GLU" and atom["atom_name"] in ["OE1", "OE2"]:
                    res2_charged.append(atom)
        
        # Check distances
        for atom1 in res1_charged:
            coords1 = np.array([atom1["x"], atom1["y"], atom1["z"]])
            for atom2 in res2_charged:
                coords2 = np.array([atom2["x"], atom2["y"], atom2["z"]])
                distance = np.linalg.norm(coords1 - coords2)
                
                if 4.0 <= distance <= 6.0:
                    return {
                        "residue1": res1_name,
                        "atom1": atom1["atom_name"],
                        "residue2": res2_name,
                        "atom2": atom2["atom_name"],
                        "distance": round(distance, 2)
                    }
        
        return None
    
    def _find_hydrophobic_interaction(
        self,
        res1_atoms: List[Dict],
        res2_atoms: List[Dict]
    ) -> Optional[Dict]:
        """Find hydrophobic interaction between two residues."""
        # Hydrophobic interaction: non-polar residues within 5 Å
        min_distance = float('inf')
        closest_pair = None
        
        for atom1 in res1_atoms:
            coords1 = np.array([atom1["x"], atom1["y"], atom1["z"]])
            for atom2 in res2_atoms:
                coords2 = np.array([atom2["x"], atom2["y"], atom2["z"]])
                distance = np.linalg.norm(coords1 - coords2)
                
                if distance < min_distance:
                    min_distance = distance
                    closest_pair = (atom1, atom2)
        
        if min_distance <= 5.0 and closest_pair:
            atom1, atom2 = closest_pair
            return {
                "atom1": atom1["atom_name"],
                "atom2": atom2["atom_name"],
                "distance": round(min_distance, 2)
            }
        
        return None


