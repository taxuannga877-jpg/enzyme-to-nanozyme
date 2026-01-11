"""
Nanozyme Structure Representation
==================================

Represents assembled nanozyme structures with:
- 3D atomic coordinates
- Atom types and properties
- Connectivity
- Catalytic sites
- Material properties
"""

from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass, field
import numpy as np
import json
from pathlib import Path

from .motif_enhanced import MaterialType, AnchorAtom, NanozymeMotif


@dataclass
class Atom:
    """Single atom in nanozyme structure"""
    element: str
    coordinates: np.ndarray
    index: int
    
    # Chemical properties
    partial_charge: float = 0.0
    oxidation_state: Optional[int] = None
    
    # Bonding
    bonded_to: List[int] = field(default_factory=list)
    
    # Role
    is_active_site: bool = False
    motif_id: Optional[str] = None
    
    def distance_to(self, other: "Atom") -> float:
        """Calculate distance to another atom"""
        return float(np.linalg.norm(self.coordinates - other.coordinates))
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'element': self.element,
            'coordinates': self.coordinates.tolist(),
            'index': self.index,
            'partial_charge': self.partial_charge,
            'oxidation_state': self.oxidation_state,
            'bonded_to': self.bonded_to,
            'is_active_site': self.is_active_site,
            'motif_id': self.motif_id,
        }


@dataclass
class Bond:
    """Bond between two atoms"""
    atom1_idx: int
    atom2_idx: int
    bond_order: float = 1.0
    bond_length: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'atom1_idx': self.atom1_idx,
            'atom2_idx': self.atom2_idx,
            'bond_order': self.bond_order,
            'bond_length': self.bond_length,
        }


class NanozymeStructure:
    """
    Complete nanozyme structure
    
    Contains:
    - All atoms with coordinates
    - Bonding information
    - Active site locations
    - Material properties
    - Motif mappings
    """
    
    def __init__(
        self,
        atoms: Optional[List[Atom]] = None,
        bonds: Optional[List[Bond]] = None,
        nanozyme_type: str = "",
        material_type: MaterialType = MaterialType.METAL_OXIDE,
        metadata: Optional[Dict] = None,
    ):
        self.atoms = atoms or []
        self.bonds = bonds or []
        self.nanozyme_type = nanozyme_type
        self.material_type = material_type
        self.metadata = metadata or {}
        
        # Motif tracking
        self.motif_origins: Dict[str, NanozymeMotif] = {}
        self.motif_atom_mapping: Dict[str, List[int]] = {}
    
    def add_atom(
        self,
        element: str,
        coordinates: np.ndarray,
        **kwargs
    ) -> int:
        """Add an atom and return its index"""
        idx = len(self.atoms)
        atom = Atom(
            element=element,
            coordinates=np.array(coordinates),
            index=idx,
            **kwargs
        )
        self.atoms.append(atom)
        return idx
    
    def add_bond(
        self,
        atom1_idx: int,
        atom2_idx: int,
        bond_order: float = 1.0,
    ):
        """Add a bond between two atoms"""
        bond_length = self.atoms[atom1_idx].distance_to(self.atoms[atom2_idx])
        bond = Bond(atom1_idx, atom2_idx, bond_order, bond_length)
        self.bonds.append(bond)
        
        # Update atom bonding lists
        self.atoms[atom1_idx].bonded_to.append(atom2_idx)
        self.atoms[atom2_idx].bonded_to.append(atom1_idx)
    
    def add_motif(
        self,
        motif: NanozymeMotif,
        atom_indices: List[int],
    ):
        """Add motif and track which atoms belong to it"""
        self.motif_origins[motif.motif_id] = motif
        self.motif_atom_mapping[motif.motif_id] = atom_indices
        
        # Mark atoms as active site
        for idx in atom_indices:
            if idx < len(self.atoms):
                self.atoms[idx].is_active_site = True
                self.atoms[idx].motif_id = motif.motif_id
    
    @property
    def num_atoms(self) -> int:
        return len(self.atoms)
    
    @property
    def num_bonds(self) -> int:
        return len(self.bonds)
    
    def get_coordinates(self) -> np.ndarray:
        """Get all atom coordinates as numpy array"""
        return np.array([atom.coordinates for atom in self.atoms])
    
    def get_elements(self) -> List[str]:
        """Get list of all elements"""
        return [atom.element for atom in self.atoms]
    
    def get_active_site_atoms(self) -> List[Atom]:
        """Get all atoms in active sites"""
        return [atom for atom in self.atoms if atom.is_active_site]
    
    def get_metal_atoms(self) -> List[Atom]:
        """Get all metal atoms"""
        from .motif_enhanced import AtomType
        return [atom for atom in self.atoms if AtomType.is_metal(atom.element)]
    
    def get_centroid(self) -> np.ndarray:
        """Get geometric center"""
        coords = self.get_coordinates()
        return coords.mean(axis=0)
    
    def translate(self, vector: np.ndarray):
        """Translate entire structure"""
        for atom in self.atoms:
            atom.coordinates += np.array(vector)
    
    def rotate(self, rotation_matrix: np.ndarray, center: Optional[np.ndarray] = None):
        """Rotate structure around center"""
        if center is None:
            center = self.get_centroid()
        
        for atom in self.atoms:
            coords_centered = atom.coordinates - center
            coords_rotated = rotation_matrix @ coords_centered
            atom.coordinates = coords_rotated + center
    
    def compute_composition(self) -> Dict[str, int]:
        """Compute chemical composition"""
        composition = {}
        for atom in self.atoms:
            composition[atom.element] = composition.get(atom.element, 0) + 1
        return composition
    
    def compute_formula(self) -> str:
        """Compute chemical formula"""
        composition = self.compute_composition()
        formula_parts = []
        for element in sorted(composition.keys()):
            count = composition[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
        return "".join(formula_parts)
    
    def to_xyz(self, filepath: str, comment: str = ""):
        """Save as XYZ file"""
        with open(filepath, 'w') as f:
            f.write(f"{self.num_atoms}\n")
            f.write(f"{comment}\n")
            for atom in self.atoms:
                x, y, z = atom.coordinates
                f.write(f"{atom.element:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    
    @classmethod
    def from_xyz(cls, filepath: str) -> "NanozymeStructure":
        """Load from XYZ file"""
        structure = cls()
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        num_atoms = int(lines[0].strip())
        for i in range(2, 2 + num_atoms):
            parts = lines[i].split()
            element = parts[0]
            coords = np.array([float(x) for x in parts[1:4]])
            structure.add_atom(element, coords)
        
        return structure
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'atoms': [atom.to_dict() for atom in self.atoms],
            'bonds': [bond.to_dict() for bond in self.bonds],
            'nanozyme_type': self.nanozyme_type,
            'material_type': self.material_type.value,
            'metadata': self.metadata,
            'motif_origins': {k: v.to_dict() for k, v in self.motif_origins.items()},
            'motif_atom_mapping': self.motif_atom_mapping,
        }
    
    def to_json(self, filepath: str):
        """Save as JSON"""
        # Need to convert np.ndarray to list first
        data = self.to_dict()
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
    
    @classmethod
    def from_json(cls, filepath: str) -> "NanozymeStructure":
        """Load from JSON"""
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        # Reconstruct structure
        structure = cls(
            nanozyme_type=data.get('nanozyme_type', ''),
            material_type=MaterialType(data.get('material_type', 'metal_oxide')),
            metadata=data.get('metadata', {}),
        )
        
        # Reconstruct atoms
        for atom_data in data.get('atoms', []):
            atom = Atom(
                element=atom_data['element'],
                coordinates=np.array(atom_data['coordinates']),
                index=atom_data['index'],
                partial_charge=atom_data.get('partial_charge', 0.0),
                oxidation_state=atom_data.get('oxidation_state'),
                bonded_to=atom_data.get('bonded_to', []),
                is_active_site=atom_data.get('is_active_site', False),
                motif_id=atom_data.get('motif_id'),
            )
            structure.atoms.append(atom)
        
        # Reconstruct bonds
        for bond_data in data.get('bonds', []):
            bond = Bond(
                atom1_idx=bond_data['atom1_idx'],
                atom2_idx=bond_data['atom2_idx'],
                bond_order=bond_data.get('bond_order', 1.0),
                bond_length=bond_data.get('bond_length', 0.0),
            )
            structure.bonds.append(bond)
        
        # Reconstruct motif mappings
        from .motif_enhanced import NanozymeMotif
        for motif_id, motif_data in data.get('motif_origins', {}).items():
            structure.motif_origins[motif_id] = NanozymeMotif.from_dict(motif_data)
        structure.motif_atom_mapping = data.get('motif_atom_mapping', {})
        
        return structure
    
    def visualize(self, output_path: Optional[str] = None):
        """Visualize structure (requires py3Dmol or similar)"""
        try:
            import py3Dmol
        except ImportError:
            print("py3Dmol not installed. Install with: pip install py3Dmol")
            return
        
        # Create temporary XYZ
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_path = f.name
            self.to_xyz(temp_path, comment=f"{self.nanozyme_type} nanozyme")
        
        # Visualize
        with open(temp_path, 'r') as f:
            xyz_str = f.read()
        
        view = py3Dmol.view(width=800, height=600)
        view.addModel(xyz_str, 'xyz')
        view.setStyle({'sphere': {'radius': 0.3}, 'stick': {'radius': 0.15}})
        view.zoomTo()
        
        if output_path:
            view.png(output_path)
        
        return view
    
    def __repr__(self) -> str:
        formula = self.compute_formula()
        return f"NanozymeStructure({formula}, {self.num_atoms} atoms, {self.num_bonds} bonds, type={self.nanozyme_type})"

