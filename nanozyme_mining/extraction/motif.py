"""
Catalytic Motif Data Structures - Stage 2 Core Module
======================================================

Defines the CatalyticMotif data structure for representing
extractable catalytic fragments from enzyme structures.

Based on user requirements:
- AnchorAtoms: Key donor atoms / metal or cofactor centers
- Geometry: Distances, angles, dihedral angles
- ChemistryTag: Reaction type / activity classification
- ReactionContext: Reaction SMILES / templates
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict
from enum import Enum


def convert_numpy_types(obj: Any) -> Any:
    """
    Recursively convert numpy types to Python native types for JSON serialization.
    
    Args:
        obj: Object that may contain numpy types
        
    Returns:
        Object with numpy types converted to Python native types
    """
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()  # Convert numpy scalar to Python native type
    elif isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert numpy array to list
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj


class AtomType(Enum):
    """Atom types for anchor atoms."""
    CARBON = "C"
    NITROGEN = "N"
    OXYGEN = "O"
    SULFUR = "S"
    METAL = "Metal"
    COFACTOR = "Cofactor"


@dataclass
class AnchorAtom:
    """
    Represents a key anchor atom in the catalytic motif.

    Attributes:
        atom_name: PDB atom name (e.g., "CA", "NE2", "OG")
        residue_name: Three-letter residue code (e.g., "HIS", "SER")
        residue_number: Residue sequence number
        chain_id: Chain identifier
        element: Element symbol
        coordinates: 3D coordinates [x, y, z]
        atom_type: Classification of atom type
        is_donor: Whether this is a donor atom
        role: Functional role (e.g., "nucleophile", "base", "acid")
    """
    atom_name: str
    residue_name: str
    residue_number: int
    chain_id: str = "A"
    element: str = ""
    coordinates: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    atom_type: str = "CARBON"
    is_donor: bool = False
    role: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnchorAtom":
        return cls(**data)


@dataclass
class GeometryConstraint:
    """
    Represents geometric constraints between anchor atoms.

    Attributes:
        constraint_type: Type of constraint (distance/angle/dihedral)
        atom_indices: Indices of atoms involved
        value: Measured value
        tolerance: Acceptable tolerance
        unit: Unit of measurement
    """
    constraint_type: str  # "distance", "angle", "dihedral"
    atom_indices: List[int] = field(default_factory=list)
    value: float = 0.0
    tolerance: float = 0.5
    unit: str = "angstrom"  # or "degree"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GeometryConstraint":
        return cls(**data)


@dataclass
class CatalyticMotif:
    """
    Complete catalytic motif representation.

    Core fields as specified:
    - AnchorAtoms: Key donor atoms / metal or cofactor centers
    - Geometry: Distances, angles, dihedral angles
    - ChemistryTag: Reaction type / activity classification
    - ReactionContext: Reaction SMILES / templates
    """
    motif_id: str
    source_uniprot_id: str
    source_ec_number: str
    nanozyme_type: str

    # Core components
    anchor_atoms: List[AnchorAtom] = field(default_factory=list)
    geometry_constraints: List[GeometryConstraint] = field(default_factory=list)

    # Chemistry and reaction context
    chemistry_tag: str = ""
    reaction_smiles: str = ""
    reaction_template: str = ""

    # Metadata
    confidence_score: float = 0.0
    extraction_method: str = ""
    notes: str = ""
    
    # Enhanced fields
    residue_structures: Dict[Tuple[str, int], Dict] = field(default_factory=dict)
    structure_2d_svg: str = ""
    
    # Extended PDB information fields
    metal_sites: List[Dict] = field(default_factory=list)
    chemical_bonds: Dict[str, List[Dict]] = field(default_factory=dict)  # SSBOND, LINK, CONECT
    ligands_and_cofactors: List[Dict] = field(default_factory=list)
    site_annotations: List[Dict] = field(default_factory=list)
    residue_environment: Dict[str, Dict] = field(default_factory=dict)  # Key: (res_name, res_num)
    secondary_structure: Dict[str, List[Dict]] = field(default_factory=dict)  # helix, sheet

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        # Convert residue_structures keys from tuples to strings for JSON serialization
        residue_structures_dict = {}
        for key, value in self.residue_structures.items():
            key_str = f"{key[0]}_{key[1]}"
            residue_structures_dict[key_str] = value
        
        result = {
            "motif_id": self.motif_id,
            "source_uniprot_id": self.source_uniprot_id,
            "source_ec_number": self.source_ec_number,
            "nanozyme_type": self.nanozyme_type,
            "anchor_atoms": [a.to_dict() for a in self.anchor_atoms],
            "geometry_constraints": [g.to_dict() for g in self.geometry_constraints],
            "chemistry_tag": self.chemistry_tag,
            "reaction_smiles": self.reaction_smiles,
            "reaction_template": self.reaction_template,
            "confidence_score": self.confidence_score,
            "extraction_method": self.extraction_method,
            "notes": self.notes,
            "residue_structures": residue_structures_dict,
            "structure_2d_svg": self.structure_2d_svg,
            "metal_sites": self.metal_sites,
            "chemical_bonds": self.chemical_bonds,
            "ligands_and_cofactors": self.ligands_and_cofactors,
            "site_annotations": self.site_annotations,
            "residue_environment": {f"{k[0]}_{k[1]}": v for k, v in self.residue_environment.items()},
            "secondary_structure": self.secondary_structure
        }
        
        # Convert numpy types to Python native types for JSON serialization
        return convert_numpy_types(result)

    def to_json(self, filepath: str) -> None:
        """Export to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def to_xyz(self, filepath: str) -> None:
        """Export anchor atoms to XYZ format."""
        lines = [str(len(self.anchor_atoms)), self.motif_id]
        for atom in self.anchor_atoms:
            x, y, z = atom.coordinates
            lines.append(f"{atom.element:2s} {x:10.5f} {y:10.5f} {z:10.5f}")
        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CatalyticMotif":
        """Create from dictionary."""
        anchor_atoms = [AnchorAtom.from_dict(a) for a in data.get("anchor_atoms", [])]
        geometry = [GeometryConstraint.from_dict(g) for g in data.get("geometry_constraints", [])]
        
        # Convert residue_structures back from string keys to tuple keys
        residue_structures = {}
        residue_structures_dict = data.get("residue_structures", {})
        for key_str, value in residue_structures_dict.items():
            parts = key_str.split("_")
            if len(parts) >= 2:
                res_name = parts[0]
                try:
                    res_num = int(parts[1])
                    residue_structures[(res_name, res_num)] = value
                except ValueError:
                    pass
        
        # Convert residue_environment back from string keys to tuple keys
        residue_environment = {}
        residue_environment_dict = data.get("residue_environment", {})
        for key_str, value in residue_environment_dict.items():
            parts = key_str.split("_")
            if len(parts) >= 2:
                res_name = parts[0]
                try:
                    res_num = int(parts[1])
                    residue_environment[(res_name, res_num)] = value
                except ValueError:
                    pass
        
        return cls(
            motif_id=data["motif_id"],
            source_uniprot_id=data["source_uniprot_id"],
            source_ec_number=data["source_ec_number"],
            nanozyme_type=data["nanozyme_type"],
            anchor_atoms=anchor_atoms,
            geometry_constraints=geometry,
            chemistry_tag=data.get("chemistry_tag", ""),
            reaction_smiles=data.get("reaction_smiles", ""),
            reaction_template=data.get("reaction_template", ""),
            confidence_score=data.get("confidence_score", 0.0),
            extraction_method=data.get("extraction_method", ""),
            notes=data.get("notes", ""),
            residue_structures=residue_structures,
            structure_2d_svg=data.get("structure_2d_svg", ""),
            metal_sites=data.get("metal_sites", []),
            chemical_bonds=data.get("chemical_bonds", {}),
            ligands_and_cofactors=data.get("ligands_and_cofactors", []),
            site_annotations=data.get("site_annotations", []),
            residue_environment=residue_environment,
            secondary_structure=data.get("secondary_structure", {})
        )
