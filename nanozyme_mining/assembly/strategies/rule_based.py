"""
Rule-Based Nanozyme Assembly
============================

Assemble nanozymes using chemical rules and heuristics.
Suitable for:
- Simple metal oxide nanozymes
- Single-atom catalysts
- Well-defined coordination complexes

Strategy:
1. Extract active site from motif
2. Build coordination environment
3. Add scaffolding/support
4. Connect multiple motifs if needed
"""

from typing import List, Dict, Optional, Tuple
import numpy as np
from scipy.spatial.transform import Rotation

from ..motif_enhanced import (
    NanozymeMotif,
    MaterialType,
    CoordinationType,
    AnchorAtom,
    AtomType,
)
from ..structure import NanozymeStructure, Atom, Bond


class RuleBasedAssembler:
    """
    Rule-based nanozyme assembler
    
    Uses chemical rules and templates to assemble nanozymes
    """
    
    # Typical bond lengths (Angstroms)
    BOND_LENGTHS = {
        ('Fe', 'O'): 2.0,
        ('Fe', 'N'): 2.1,
        ('Fe', 'S'): 2.3,
        ('Fe', 'C'): 2.0,
        ('Cu', 'O'): 1.95,
        ('Cu', 'N'): 2.0,
        ('Mn', 'O'): 1.9,
        ('Co', 'N'): 2.0,
        ('Ni', 'N'): 2.0,
        ('C', 'C'): 1.54,
        ('C', 'N'): 1.47,
        ('C', 'O'): 1.43,
        ('N', 'N'): 1.45,
    }
    
    # Coordination geometries and ideal angles
    COORDINATION_GEOMETRIES = {
        CoordinationType.OCTAHEDRAL: {
            'angles': [90, 90, 90, 90, 180],  # Ideal angles
            'positions': [  # Normalized positions
                [1, 0, 0], [-1, 0, 0],
                [0, 1, 0], [0, -1, 0],
                [0, 0, 1], [0, 0, -1],
            ]
        },
        CoordinationType.TETRAHEDRAL: {
            'angles': [109.47, 109.47, 109.47],
            'positions': [
                [1, 1, 1], [-1, -1, 1],
                [-1, 1, -1], [1, -1, -1],
            ]
        },
        CoordinationType.SQUARE_PLANAR: {
            'angles': [90, 90, 90, 90],
            'positions': [
                [1, 0, 0], [-1, 0, 0],
                [0, 1, 0], [0, -1, 0],
            ]
        },
    }
    
    def __init__(
        self,
        material_type: MaterialType = MaterialType.METAL_OXIDE,
        output_dir: str = "./assembled",
    ):
        self.material_type = material_type
        self.output_dir = output_dir
    
    def assemble(
        self,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble nanozyme from motif(s)
        
        Args:
            motifs: List of catalytic motifs
            **kwargs: Additional parameters
                - num_metal_centers: Number of metal centers
                - scaffold_size: Size of support scaffold
                - spacing: Distance between motifs
        """
        num_metal_centers = kwargs.get('num_metal_centers', len(motifs))
        scaffold_size = kwargs.get('scaffold_size', 20)
        spacing = kwargs.get('spacing', 5.0)
        
        # Create structure
        structure = NanozymeStructure(
            nanozyme_type=motifs[0].nanozyme_type if motifs else "Unknown",
            material_type=self.material_type,
        )
        
        # Strategy depends on material type
        if self.material_type == MaterialType.METAL_OXIDE:
            return self._assemble_metal_oxide(structure, motifs, **kwargs)
        elif self.material_type == MaterialType.SAC:
            return self._assemble_single_atom(structure, motifs, **kwargs)
        elif self.material_type == MaterialType.METAL_COMPLEX:
            return self._assemble_metal_complex(structure, motifs, **kwargs)
        elif self.material_type == MaterialType.MOF:
            return self._assemble_mof(structure, motifs, **kwargs)
        else:
            # Generic assembly
            return self._assemble_generic(structure, motifs, **kwargs)
    
    def _assemble_metal_oxide(
        self,
        structure: NanozymeStructure,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble metal oxide nanozyme (e.g., Fe3O4, CeO2)
        
        Strategy:
        1. Create metal oxide core/cluster
        2. Place active sites on surface
        3. Add coordination environment
        """
        size = kwargs.get('size', 20)  # Number of atoms
        
        # For each motif, create a metal center with coordination
        for i, motif in enumerate(motifs):
            # Get metal atoms from motif
            metal_atoms = motif.get_metal_atoms()
            if not metal_atoms:
                # If no metal specified, infer from nanozyme type
                metal_atoms = self._infer_metal_atoms(motif)
            
            # Place metal center
            position = self._get_motif_position(i, len(motifs), **kwargs)
            
            for metal_atom in metal_atoms:
                # Add metal
                metal_idx = structure.add_atom(
                    element=metal_atom.element,
                    coordinates=position,
                    is_active_site=True,
                    motif_id=motif.motif_id,
                    oxidation_state=metal_atom.oxidation_state,
                )
                
                # Add coordinating atoms (typically O for metal oxides)
                coord_atoms = motif.get_coordinating_atoms()
                if not coord_atoms:
                    # Generate default coordination sphere
                    coord_atoms = self._generate_coordination_sphere(
                        metal_atom,
                        coordination_type=CoordinationType.OCTAHEDRAL,
                    )
                
                # Place coordinating atoms
                for j, coord_atom in enumerate(coord_atoms):
                    # Get position relative to metal
                    offset = self._get_coordination_position(
                        metal_atom.element,
                        coord_atom.element,
                        j,
                        len(coord_atoms),
                    )
                    coord_pos = position + offset
                    
                    coord_idx = structure.add_atom(
                        element=coord_atom.element,
                        coordinates=coord_pos,
                        is_active_site=True,
                        motif_id=motif.motif_id,
                    )
                    
                    # Add bond
                    structure.add_bond(metal_idx, coord_idx)
            
            # Track motif
            structure.add_motif(motif, list(range(structure.num_atoms)))
        
        return structure
    
    def _assemble_single_atom(
        self,
        structure: NanozymeStructure,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble single-atom catalyst (e.g., Fe-N4-C)
        
        Strategy:
        1. Place metal center
        2. Add N4 coordination
        3. Add carbon support
        """
        for i, motif in enumerate(motifs):
            metal_atoms = motif.get_metal_atoms()
            if not metal_atoms:
                print(f"Warning: No metal atoms in motif {motif.motif_id}")
                continue
            
            metal_atom = metal_atoms[0]  # Single atom
            position = self._get_motif_position(i, len(motifs), **kwargs)
            
            # Add metal center
            metal_idx = structure.add_atom(
                element=metal_atom.element,
                coordinates=position,
                is_active_site=True,
                motif_id=motif.motif_id,
                oxidation_state=metal_atom.oxidation_state,
            )
            
            # Add N4 coordination (square planar)
            bond_length = self.BOND_LENGTHS.get((metal_atom.element, 'N'), 2.0)
            positions = self.COORDINATION_GEOMETRIES[CoordinationType.SQUARE_PLANAR]['positions']
            
            for j, rel_pos in enumerate(positions):
                n_pos = position + np.array(rel_pos) * bond_length
                n_idx = structure.add_atom(
                    element='N',
                    coordinates=n_pos,
                    is_active_site=True,
                    motif_id=motif.motif_id,
                )
                structure.add_bond(metal_idx, n_idx)
                
                # Add carbon connected to N (pyridinic-N or pyrrolic-N)
                c_offset = np.array(rel_pos) * (bond_length + 1.4)  # C-N bond ~1.4Å
                c_pos = position + c_offset
                c_idx = structure.add_atom(
                    element='C',
                    coordinates=c_pos,
                    motif_id=motif.motif_id,
                )
                structure.add_bond(n_idx, c_idx)
            
            structure.add_motif(motif, list(range(structure.num_atoms)))
        
        return structure
    
    def _assemble_metal_complex(
        self,
        structure: NanozymeStructure,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble metal complex nanozyme
        
        Strategy:
        1. Place metal center(s)
        2. Add ligands from motif
        3. Respect coordination geometry
        """
        for i, motif in enumerate(motifs):
            # Extract metal properties
            if not motif.metal_centers:
                print(f"Warning: No metal centers defined in motif {motif.motif_id}")
                continue
            
            metal_prop = motif.metal_centers[0]
            position = self._get_motif_position(i, len(motifs), **kwargs)
            
            # Add metal
            metal_idx = structure.add_atom(
                element=metal_prop.element,
                coordinates=position,
                is_active_site=True,
                oxidation_state=metal_prop.oxidation_state,
                motif_id=motif.motif_id,
            )
            
            # Add ligands according to coordination geometry
            coord_geometry = metal_prop.coordination_geometry
            coord_number = metal_prop.coordination_number
            
            if coord_geometry in self.COORDINATION_GEOMETRIES:
                positions = self.COORDINATION_GEOMETRIES[coord_geometry]['positions']
            else:
                # Generate positions for irregular coordination
                positions = self._generate_irregular_positions(coord_number)
            
            # Place ligands
            coordinating_atoms = motif.get_coordinating_atoms()
            for j in range(min(coord_number, len(positions))):
                if j < len(coordinating_atoms):
                    element = coordinating_atoms[j].element
                else:
                    element = 'O'  # Default to oxygen
                
                bond_length = self.BOND_LENGTHS.get((metal_prop.element, element), 2.0)
                ligand_pos = position + np.array(positions[j]) * bond_length
                
                ligand_idx = structure.add_atom(
                    element=element,
                    coordinates=ligand_pos,
                    is_active_site=True,
                    motif_id=motif.motif_id,
                )
                structure.add_bond(metal_idx, ligand_idx)
            
            structure.add_motif(motif, list(range(structure.num_atoms)))
        
        return structure
    
    def _assemble_mof(
        self,
        structure: NanozymeStructure,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble MOF nanozyme
        
        Strategy:
        1. Extract metal node from motif
        2. Add organic linkers
        3. Create periodic structure (simplified)
        """
        # For simplicity, just create metal nodes with carboxylate linkers
        for i, motif in enumerate(motifs):
            position = self._get_motif_position(i, len(motifs), **kwargs)
            
            # Metal node (e.g., Fe3-O cluster)
            metal_atoms = motif.get_metal_atoms()
            if not metal_atoms:
                continue
            
            # Create tri-metallic node (simplified)
            metal_element = metal_atoms[0].element
            for j, angle in enumerate([0, 120, 240]):
                angle_rad = np.radians(angle)
                metal_pos = position + np.array([
                    2.0 * np.cos(angle_rad),
                    2.0 * np.sin(angle_rad),
                    0
                ])
                
                metal_idx = structure.add_atom(
                    element=metal_element,
                    coordinates=metal_pos,
                    is_active_site=True,
                    motif_id=motif.motif_id,
                )
                
                # Add bridging O
                if j > 0:
                    prev_metal_idx = metal_idx - 1
                    mid_point = (structure.atoms[metal_idx].coordinates + 
                                structure.atoms[prev_metal_idx].coordinates) / 2
                    o_idx = structure.add_atom(
                        element='O',
                        coordinates=mid_point,
                        motif_id=motif.motif_id,
                    )
                    structure.add_bond(metal_idx, o_idx)
                    structure.add_bond(prev_metal_idx, o_idx)
            
            structure.add_motif(motif, list(range(structure.num_atoms)))
        
        return structure
    
    def _assemble_generic(
        self,
        structure: NanozymeStructure,
        motifs: List[NanozymeMotif],
        **kwargs
    ) -> NanozymeStructure:
        """Generic assembly: just place atoms from motifs"""
        for i, motif in enumerate(motifs):
            position = self._get_motif_position(i, len(motifs), **kwargs)
            
            atom_indices = []
            for anchor_atom in motif.anchor_atoms:
                # Place atom
                atom_pos = position + np.array(anchor_atom.coordinates)
                idx = structure.add_atom(
                    element=anchor_atom.element,
                    coordinates=atom_pos,
                    is_active_site=True,
                    motif_id=motif.motif_id,
                    oxidation_state=anchor_atom.oxidation_state,
                )
                atom_indices.append(idx)
            
            # Add bonds based on geometry constraints
            for constraint in motif.geometry_constraints:
                if constraint.constraint_type == 'distance' and len(constraint.atom_indices) == 2:
                    local_idx1, local_idx2 = constraint.atom_indices
                    if local_idx1 < len(atom_indices) and local_idx2 < len(atom_indices):
                        structure.add_bond(
                            atom_indices[local_idx1],
                            atom_indices[local_idx2],
                        )
            
            structure.add_motif(motif, atom_indices)
        
        return structure
    
    # Helper methods
    
    def _get_motif_position(
        self,
        index: int,
        total: int,
        **kwargs
    ) -> np.ndarray:
        """Get position for motif placement"""
        spacing = kwargs.get('spacing', 5.0)
        
        if total == 1:
            return np.zeros(3)
        elif total == 2:
            return np.array([index * spacing, 0, 0])
        else:
            # Arrange in a grid or circle
            angle = 2 * np.pi * index / total
            radius = spacing
            return np.array([
                radius * np.cos(angle),
                radius * np.sin(angle),
                0
            ])
    
    def _get_coordination_position(
        self,
        metal: str,
        ligand: str,
        index: int,
        total: int,
    ) -> np.ndarray:
        """Get position of coordinating atom relative to metal"""
        bond_length = self.BOND_LENGTHS.get((metal, ligand), 2.0)
        
        # Octahedral by default
        if total <= 6:
            positions = self.COORDINATION_GEOMETRIES[CoordinationType.OCTAHEDRAL]['positions']
            if index < len(positions):
                return np.array(positions[index]) * bond_length
        
        # Random positioning
        random_dir = np.random.randn(3)
        random_dir /= np.linalg.norm(random_dir)
        return random_dir * bond_length
    
    def _infer_metal_atoms(self, motif: NanozymeMotif) -> List[AnchorAtom]:
        """Infer metal atoms based on nanozyme type"""
        metal_map = {
            'POD': 'Fe',
            'CAT': 'Fe',
            'SOD': 'Cu',
            'GSH': 'Se',
            'OXD': 'Cu',
            'LAC': 'Cu',
            'GOX': 'Fe',
        }
        
        metal = metal_map.get(motif.nanozyme_type, 'Fe')
        return [AnchorAtom(
            atom_name=f"{metal}1",
            element=metal,
            coordinates=[0, 0, 0],
            is_metal_center=True,
        )]
    
    def _generate_coordination_sphere(
        self,
        metal_atom: AnchorAtom,
        coordination_type: CoordinationType = CoordinationType.OCTAHEDRAL,
    ) -> List[AnchorAtom]:
        """Generate default coordination sphere"""
        coord_atoms = []
        positions = self.COORDINATION_GEOMETRIES[coordination_type]['positions']
        bond_length = self.BOND_LENGTHS.get((metal_atom.element, 'O'), 2.0)
        
        for i, pos in enumerate(positions):
            coord = (np.array(metal_atom.coordinates) + 
                    np.array(pos) * bond_length).tolist()
            coord_atoms.append(AnchorAtom(
                atom_name=f"O{i+1}",
                element='O',
                coordinates=coord,
                is_donor=True,
            ))
        
        return coord_atoms
    
    def _generate_irregular_positions(self, num_positions: int) -> List[List[float]]:
        """Generate positions for irregular coordination"""
        positions = []
        for i in range(num_positions):
            # Distribute on sphere
            theta = np.arccos(1 - 2 * (i + 0.5) / num_positions)
            phi = np.pi * (1 + 5**0.5) * i
            
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            
            positions.append([x, y, z])
        
        return positions

