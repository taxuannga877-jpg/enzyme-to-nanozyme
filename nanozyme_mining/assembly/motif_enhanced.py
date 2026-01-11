"""
Enhanced Motif Representation for Nanozymes
===========================================

Extended from basic CatalyticMotif to support nanozyme materials:
- Metal oxides (Fe3O4, CeO2, MnO2, etc.)
- Metal-organic frameworks (MOFs)
- Single-atom catalysts (SACs)
- Carbon-based nanozymes
- Metal clusters

IMPORTANT: Nanozymes are NOT proteins!
They are artificial nanomaterials with enzyme-like activities.
"""

from enum import Enum
from typing import List, Dict, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
import json
import numpy as np


class MaterialType(Enum):
    """Types of nanozyme materials"""
    METAL_OXIDE = "metal_oxide"          # Fe3O4, CeO2, MnO2等
    CARBON_BASED = "carbon_based"        # Graphene, CNT, C-dots
    MOF = "metal_organic_framework"      # MOF, ZIF等
    METAL_CLUSTER = "metal_cluster"      # Pt, Au, Ag团簇
    SAC = "single_atom_catalyst"         # Fe-N4-C等单原子
    METAL_COMPLEX = "metal_complex"      # 金属配合物
    HYBRID = "hybrid"                    # 复合材料


class CoordinationType(Enum):
    """Metal coordination geometries"""
    OCTAHEDRAL = "octahedral"            # 八面体 (CN=6)
    TETRAHEDRAL = "tetrahedral"          # 四面体 (CN=4)
    SQUARE_PLANAR = "square_planar"      # 平面正方形 (CN=4)
    SQUARE_PYRAMIDAL = "square_pyramidal"  # 正方锥 (CN=5)
    TRIGONAL_BIPYRAMIDAL = "trigonal_bipyramidal"  # 三角双锥 (CN=5)
    LINEAR = "linear"                    # 线性 (CN=2)
    BENT = "bent"                        # 弯曲
    IRREGULAR = "irregular"              # 不规则
    UNKNOWN = "unknown"


class AtomType(Enum):
    """Extended atom types for nanozymes"""
    # Metals (过渡金属)
    FE = "Fe"
    CU = "Cu"
    MN = "Mn"
    CO = "Co"
    NI = "Ni"
    ZN = "Zn"
    CE = "Ce"
    PT = "Pt"
    AU = "Au"
    AG = "Ag"
    PD = "Pd"
    RU = "Ru"
    V = "V"
    TI = "Ti"
    
    # Non-metals
    C = "C"
    N = "N"
    O = "O"
    S = "S"
    P = "P"
    H = "H"
    F = "F"
    CL = "Cl"
    BR = "Br"
    I = "I"
    
    @classmethod
    def is_metal(cls, atom_symbol: str) -> bool:
        """Check if atom is a metal"""
        metals = {'Fe', 'Cu', 'Mn', 'Co', 'Ni', 'Zn', 'Ce', 'Pt', 'Au', 'Ag', 'Pd', 'Ru', 'V', 'Ti'}
        return atom_symbol in metals
    
    @classmethod
    def get_atomic_number(cls, atom_symbol: str) -> int:
        """Get atomic number"""
        atomic_numbers = {
            'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
            'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53,
            'Ti': 22, 'V': 23, 'Mn': 25, 'Fe': 26, 'Co': 27,
            'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ru': 44, 'Pd': 46,
            'Ag': 47, 'Pt': 78, 'Au': 79, 'Ce': 58,
        }
        return atomic_numbers.get(atom_symbol, 0)


@dataclass
class AnchorAtom:
    """
    Enhanced anchor atom for nanozymes
    
    Can represent:
    - Metal centers
    - Coordinating atoms (N, O, S in ligands)
    - Surface atoms on nanoparticles
    """
    atom_name: str
    element: str
    coordinates: List[float]
    
    # For protein-derived motifs
    residue_name: Optional[str] = None
    residue_number: Optional[int] = None
    chain_id: Optional[str] = None
    
    # Nanozyme-specific properties
    is_metal_center: bool = False
    is_surface_atom: bool = False
    oxidation_state: Optional[int] = None
    coordination_number: Optional[int] = None
    
    # Chemical properties
    partial_charge: Optional[float] = None
    is_donor: bool = False
    is_acceptor: bool = False
    
    # Role in catalysis
    role: str = ""  # "metal_center", "coordinating_atom", "substrate_binding"等
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnchorAtom":
        return cls(**data)
    
    def is_metal(self) -> bool:
        """Check if this atom is a metal"""
        return AtomType.is_metal(self.element)


@dataclass
class GeometryConstraint:
    """
    Geometric constraint for nanozyme assembly
    
    Types:
    - distance: atom-atom distance
    - angle: three-atom angle
    - dihedral: four-atom dihedral
    - coordination: metal coordination constraint
    """
    constraint_type: str
    atom_indices: List[int]
    value: float
    tolerance: float = 0.5
    unit: str = "angstrom"
    
    # For coordination constraints
    coordination_type: Optional[CoordinationType] = None
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        if self.coordination_type:
            data['coordination_type'] = self.coordination_type.value
        return data
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GeometryConstraint":
        if 'coordination_type' in data and data['coordination_type']:
            data['coordination_type'] = CoordinationType(data['coordination_type'])
        return cls(**data)


@dataclass
class MetalProperties:
    """Properties of metal center in nanozyme"""
    element: str
    oxidation_state: int
    coordination_number: int
    coordination_geometry: CoordinationType
    d_electron_count: int
    spin_state: Optional[str] = None  # "high_spin", "low_spin"
    ligand_field_strength: Optional[str] = None  # "weak", "medium", "strong"
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['coordination_geometry'] = self.coordination_geometry.value
        return data
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MetalProperties":
        if 'coordination_geometry' in data:
            data['coordination_geometry'] = CoordinationType(data['coordination_geometry'])
        return cls(**data)


@dataclass
class NanozymeMotif:
    """
    Enhanced catalytic motif for nanozymes
    
    This extends the basic CatalyticMotif to support:
    1. Various nanozyme materials (not just proteins)
    2. Metal centers and coordination chemistry
    3. Surface catalysis
    4. Material-specific properties
    """
    # Basic identification
    motif_id: str
    nanozyme_type: str  # POD/CAT/SOD等
    source_uniprot_id: Optional[str] = None
    source_ec_number: Optional[str] = None
    
    # Structural components
    anchor_atoms: List[AnchorAtom] = field(default_factory=list)
    geometry_constraints: List[GeometryConstraint] = field(default_factory=list)
    
    # Material properties
    material_type: MaterialType = MaterialType.METAL_OXIDE
    
    # Metal properties (if applicable)
    metal_centers: List[MetalProperties] = field(default_factory=list)
    
    # Chemical information
    chemistry_tag: str = ""
    reaction_smiles: str = ""
    reaction_template: str = ""
    
    # Quality metrics
    confidence_score: float = 0.0
    extraction_method: str = "rule_based"
    
    # Metadata
    notes: str = ""
    references: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'motif_id': self.motif_id,
            'nanozyme_type': self.nanozyme_type,
            'source_uniprot_id': self.source_uniprot_id,
            'source_ec_number': self.source_ec_number,
            'anchor_atoms': [atom.to_dict() for atom in self.anchor_atoms],
            'geometry_constraints': [c.to_dict() for c in self.geometry_constraints],
            'material_type': self.material_type.value,
            'metal_centers': [m.to_dict() for m in self.metal_centers],
            'chemistry_tag': self.chemistry_tag,
            'reaction_smiles': self.reaction_smiles,
            'reaction_template': self.reaction_template,
            'confidence_score': self.confidence_score,
            'extraction_method': self.extraction_method,
            'notes': self.notes,
            'references': self.references,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "NanozymeMotif":
        """Load from dictionary"""
        data = data.copy()
        data['anchor_atoms'] = [AnchorAtom.from_dict(a) for a in data.get('anchor_atoms', [])]
        data['geometry_constraints'] = [GeometryConstraint.from_dict(c) for c in data.get('geometry_constraints', [])]
        if 'material_type' in data:
            data['material_type'] = MaterialType(data['material_type'])
        data['metal_centers'] = [MetalProperties.from_dict(m) for m in data.get('metal_centers', [])]
        return cls(**data)
    
    def save(self, filepath: str):
        """Save motif to JSON file"""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    @classmethod
    def load(cls, filepath: str) -> "NanozymeMotif":
        """Load motif from JSON file"""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
    
    def get_metal_atoms(self) -> List[AnchorAtom]:
        """Get all metal center atoms"""
        return [atom for atom in self.anchor_atoms if atom.is_metal_center]
    
    def get_coordinating_atoms(self) -> List[AnchorAtom]:
        """Get coordinating atoms around metals"""
        return [atom for atom in self.anchor_atoms if atom.is_donor or atom.is_acceptor]
    
    def get_centroid(self) -> np.ndarray:
        """Get geometric center of motif"""
        if not self.anchor_atoms:
            return np.zeros(3)
        coords = np.array([atom.coordinates for atom in self.anchor_atoms])
        return coords.mean(axis=0)
    
    def translate(self, vector: np.ndarray):
        """Translate motif by vector"""
        for atom in self.anchor_atoms:
            atom.coordinates = (np.array(atom.coordinates) + vector).tolist()
    
    def rotate(self, rotation_matrix: np.ndarray, center: Optional[np.ndarray] = None):
        """Rotate motif around center"""
        if center is None:
            center = self.get_centroid()
        
        for atom in self.anchor_atoms:
            coords = np.array(atom.coordinates)
            coords_centered = coords - center
            coords_rotated = rotation_matrix @ coords_centered
            atom.coordinates = (coords_rotated + center).tolist()


def convert_basic_to_nanozyme_motif(
    basic_motif_path: str,
    material_type: MaterialType = MaterialType.METAL_OXIDE,
) -> NanozymeMotif:
    """
    Convert basic CatalyticMotif to NanozymeMotif
    
    This is a helper function to migrate existing motifs
    """
    with open(basic_motif_path, 'r') as f:
        data = json.load(f)
    
    # Convert anchor atoms
    anchor_atoms = []
    for atom_data in data.get('anchor_atoms', []):
        # Detect if metal
        element = atom_data.get('element', 'C')
        is_metal = AtomType.is_metal(element)
        
        anchor_atom = AnchorAtom(
            atom_name=atom_data.get('atom_name', ''),
            element=element,
            coordinates=atom_data.get('coordinates', [0, 0, 0]),
            residue_name=atom_data.get('residue_name'),
            residue_number=atom_data.get('residue_number'),
            chain_id=atom_data.get('chain_id'),
            is_metal_center=is_metal,
            is_donor=atom_data.get('is_donor', False),
            role=atom_data.get('role', ''),
        )
        anchor_atoms.append(anchor_atom)
    
    # Convert geometry constraints
    geometry_constraints = []
    for constraint_data in data.get('geometry_constraints', []):
        constraint = GeometryConstraint.from_dict(constraint_data)
        geometry_constraints.append(constraint)
    
    # Detect metal centers
    metal_centers = []
    for atom in anchor_atoms:
        if atom.is_metal_center:
            # Infer coordination from constraints
            coord_number = sum(
                1 for c in geometry_constraints
                if c.constraint_type == 'distance' and atom in [anchor_atoms[i] for i in c.atom_indices]
            )
            
            metal_prop = MetalProperties(
                element=atom.element,
                oxidation_state=atom.oxidation_state or 2,  # Default +2
                coordination_number=coord_number,
                coordination_geometry=CoordinationType.OCTAHEDRAL if coord_number == 6 else CoordinationType.TETRAHEDRAL,
                d_electron_count=0,  # Need to calculate
            )
            metal_centers.append(metal_prop)
    
    # Create NanozymeMotif
    nanozyme_motif = NanozymeMotif(
        motif_id=data.get('motif_id', ''),
        nanozyme_type=data.get('nanozyme_type', ''),
        source_uniprot_id=data.get('source_uniprot_id'),
        source_ec_number=data.get('source_ec_number'),
        anchor_atoms=anchor_atoms,
        geometry_constraints=geometry_constraints,
        material_type=material_type,
        metal_centers=metal_centers,
        chemistry_tag=data.get('chemistry_tag', ''),
        reaction_smiles=data.get('reaction_smiles', ''),
        reaction_template=data.get('reaction_template', ''),
        confidence_score=data.get('confidence_score', 0.0),
        extraction_method=data.get('extraction_method', 'rule_based'),
        notes=data.get('notes', ''),
    )
    
    return nanozyme_motif

