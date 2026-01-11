"""
Nanozyme Structure Validator
============================

Validates assembled nanozyme structures for:
1. Chemical validity (bond lengths, angles, coordination)
2. Physical plausibility (no overlapping atoms, reasonable geometry)
3. Catalytic competence (active site accessibility, substrate binding)
"""

from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field
import numpy as np

from .structure import NanozymeStructure, Atom
from .motif_enhanced import AtomType, CoordinationType


@dataclass
class ValidationResult:
    """Result of structure validation"""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    scores: Dict[str, float] = field(default_factory=dict)
    
    def add_error(self, message: str):
        """Add an error"""
        self.errors.append(message)
        self.is_valid = False
    
    def add_warning(self, message: str):
        """Add a warning"""
        self.warnings.append(message)
    
    def __repr__(self) -> str:
        status = "VALID" if self.is_valid else "INVALID"
        return f"ValidationResult({status}, {len(self.errors)} errors, {len(self.warnings)} warnings)"


class NanozymeValidator:
    """
    Validator for nanozyme structures
    
    Checks chemical and physical validity of assembled nanozymes
    """
    
    # Typical bond length ranges (Angstroms)
    BOND_LENGTH_RANGES = {
        ('Fe', 'O'): (1.7, 2.3),
        ('Fe', 'N'): (1.8, 2.4),
        ('Fe', 'S'): (2.0, 2.6),
        ('Cu', 'O'): (1.7, 2.2),
        ('Cu', 'N'): (1.8, 2.3),
        ('Mn', 'O'): (1.6, 2.2),
        ('Co', 'N'): (1.8, 2.3),
        ('Ni', 'N'): (1.8, 2.3),
        ('C', 'C'): (1.2, 1.7),
        ('C', 'N'): (1.2, 1.6),
        ('C', 'O'): (1.2, 1.5),
        ('N', 'O'): (1.2, 1.5),
    }
    
    # Minimum separation between non-bonded atoms
    MIN_SEPARATION = 1.5  # Angstroms
    
    # Typical coordination numbers
    COORDINATION_NUMBERS = {
        'Fe': [4, 5, 6],
        'Cu': [4, 5, 6],
        'Mn': [4, 5, 6],
        'Co': [4, 5, 6],
        'Ni': [4, 5, 6],
        'Zn': [4, 5, 6],
    }
    
    def __init__(
        self,
        strict: bool = False,
        check_coordination: bool = True,
        check_geometry: bool = True,
        check_accessibility: bool = True,
    ):
        self.strict = strict
        self.check_coordination = check_coordination
        self.check_geometry = check_geometry
        self.check_accessibility = check_accessibility
    
    def validate(self, structure: NanozymeStructure) -> ValidationResult:
        """
        Validate nanozyme structure
        
        Args:
            structure: Nanozyme structure to validate
        
        Returns:
            ValidationResult with errors, warnings, and scores
        """
        result = ValidationResult(is_valid=True)
        
        # Basic checks
        if structure.num_atoms == 0:
            result.add_error("Structure has no atoms")
            return result
        
        # Check bond lengths
        if self.check_geometry:
            self._check_bond_lengths(structure, result)
            self._check_atomic_overlaps(structure, result)
        
        # Check metal coordination
        if self.check_coordination:
            self._check_metal_coordination(structure, result)
        
        # Check active site accessibility
        if self.check_accessibility:
            self._check_active_site_accessibility(structure, result)
        
        # Compute scores
        result.scores['completeness'] = self._compute_completeness(structure)
        result.scores['geometry_quality'] = self._compute_geometry_quality(structure)
        result.scores['coordination_quality'] = self._compute_coordination_quality(structure)
        
        return result
    
    def _check_bond_lengths(self, structure: NanozymeStructure, result: ValidationResult):
        """Check if bond lengths are reasonable"""
        for bond in structure.bonds:
            atom1 = structure.atoms[bond.atom1_idx]
            atom2 = structure.atoms[bond.atom2_idx]
            
            length = bond.bond_length
            element_pair = tuple(sorted([atom1.element, atom2.element]))
            
            if element_pair in self.BOND_LENGTH_RANGES:
                min_len, max_len = self.BOND_LENGTH_RANGES[element_pair]
                if length < min_len or length > max_len:
                    msg = f"Bond {atom1.element}-{atom2.element} length {length:.2f}Å out of range [{min_len}, {max_len}]"
                    if self.strict:
                        result.add_error(msg)
                    else:
                        result.add_warning(msg)
    
    def _check_atomic_overlaps(self, structure: NanozymeStructure, result: ValidationResult):
        """Check for overlapping atoms"""
        coords = structure.get_coordinates()
        n_atoms = len(coords)
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                distance = np.linalg.norm(coords[i] - coords[j])
                if distance < self.MIN_SEPARATION:
                    # Check if they're bonded
                    is_bonded = j in structure.atoms[i].bonded_to
                    if not is_bonded:
                        msg = f"Atoms {i} ({structure.atoms[i].element}) and {j} ({structure.atoms[j].element}) too close: {distance:.2f}Å"
                        result.add_error(msg)
    
    def _check_metal_coordination(self, structure: NanozymeStructure, result: ValidationResult):
        """Check if metal coordination is reasonable"""
        metal_atoms = structure.get_metal_atoms()
        
        for metal in metal_atoms:
            coord_number = len(metal.bonded_to)
            
            # Check if coordination number is typical
            if metal.element in self.COORDINATION_NUMBERS:
                typical_coords = self.COORDINATION_NUMBERS[metal.element]
                if coord_number not in typical_coords:
                    msg = f"Metal {metal.element} (atom {metal.index}) has unusual coordination number {coord_number} (typical: {typical_coords})"
                    if self.strict:
                        result.add_error(msg)
                    else:
                        result.add_warning(msg)
            
            # Check coordination geometry
            if coord_number > 0:
                geometry_score = self._check_coordination_geometry(metal, structure)
                if geometry_score < 0.5:
                    msg = f"Metal {metal.element} (atom {metal.index}) has poor coordination geometry (score: {geometry_score:.2f})"
                    result.add_warning(msg)
    
    def _check_coordination_geometry(self, metal: Atom, structure: NanozymeStructure) -> float:
        """
        Check quality of coordination geometry
        
        Returns score 0-1 (1 = perfect geometry)
        """
        if len(metal.bonded_to) == 0:
            return 0.0
        
        # Get coordinating atom positions
        coord_positions = []
        for bonded_idx in metal.bonded_to:
            bonded_atom = structure.atoms[bonded_idx]
            # Vector from metal to ligand
            vec = bonded_atom.coordinates - metal.coordinates
            vec_norm = vec / np.linalg.norm(vec)
            coord_positions.append(vec_norm)
        
        # For octahedral (6-coordinate), check if angles are close to 90° or 180°
        if len(coord_positions) == 6:
            ideal_angles = [90, 90, 90, 90, 180]
            return self._compute_angle_score(coord_positions, ideal_angles)
        
        # For tetrahedral (4-coordinate), check 109.47°
        elif len(coord_positions) == 4:
            ideal_angles = [109.47] * 6
            return self._compute_angle_score(coord_positions, ideal_angles)
        
        # For square planar (4-coordinate)
        elif len(coord_positions) == 4:
            # Check if all in same plane
            planarity = self._check_planarity(coord_positions)
            return planarity
        
        return 0.5  # Default for unknown coordination
    
    def _compute_angle_score(self, positions: List[np.ndarray], ideal_angles: List[float]) -> float:
        """Compute how close angles are to ideal"""
        if len(positions) < 2:
            return 1.0
        
        # Compute all pairwise angles
        angles = []
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                cos_angle = np.dot(positions[i], positions[j])
                cos_angle = np.clip(cos_angle, -1, 1)
                angle = np.degrees(np.arccos(cos_angle))
                angles.append(angle)
        
        # Compare with ideal angles
        if len(angles) != len(ideal_angles):
            return 0.5
        
        # Compute RMSD from ideal
        angles_sorted = sorted(angles)
        ideal_sorted = sorted(ideal_angles)
        rmsd = np.sqrt(np.mean([(a - i)**2 for a, i in zip(angles_sorted, ideal_sorted)]))
        
        # Convert to score (0 RMSD = 1.0, 30° RMSD = 0.0)
        score = max(0, 1 - rmsd / 30)
        return score
    
    def _check_planarity(self, positions: List[np.ndarray]) -> float:
        """Check if points are coplanar (for square planar)"""
        if len(positions) < 4:
            return 1.0
        
        # Fit plane
        centroid = np.mean(positions, axis=0)
        centered = positions - centroid
        
        # SVD
        U, S, Vt = np.linalg.svd(centered)
        normal = Vt[-1]
        
        # Check distances from plane
        distances = [abs(np.dot(pos - centroid, normal)) for pos in positions]
        max_dist = max(distances)
        
        # Score: 0 distance = 1.0, 1Å distance = 0.0
        score = max(0, 1 - max_dist)
        return score
    
    def _check_active_site_accessibility(self, structure: NanozymeStructure, result: ValidationResult):
        """Check if active sites are accessible"""
        active_atoms = structure.get_active_site_atoms()
        
        if len(active_atoms) == 0:
            result.add_warning("No active site atoms found")
            return
        
        # Check if active atoms are buried
        for active_atom in active_atoms:
            neighbors_close = self._count_neighbors_within(
                active_atom,
                structure,
                radius=3.0,
            )
            
            if neighbors_close > 12:  # Highly coordinated = buried
                result.add_warning(
                    f"Active site atom {active_atom.index} ({active_atom.element}) "
                    f"may be buried ({neighbors_close} neighbors within 3Å)"
                )
    
    def _count_neighbors_within(self, atom: Atom, structure: NanozymeStructure, radius: float) -> int:
        """Count number of atoms within radius"""
        count = 0
        for other in structure.atoms:
            if other.index != atom.index:
                distance = np.linalg.norm(atom.coordinates - other.coordinates)
                if distance < radius:
                    count += 1
        return count
    
    def _compute_completeness(self, structure: NanozymeStructure) -> float:
        """Compute structure completeness score"""
        # Check if structure has:
        # - Atoms
        # - Bonds
        # - Active sites
        # - Metal centers
        
        score = 0.0
        
        if structure.num_atoms > 0:
            score += 0.25
        
        if structure.num_bonds > 0:
            score += 0.25
        
        if len(structure.get_active_site_atoms()) > 0:
            score += 0.25
        
        if len(structure.get_metal_atoms()) > 0:
            score += 0.25
        
        return score
    
    def _compute_geometry_quality(self, structure: NanozymeStructure) -> float:
        """Compute overall geometry quality"""
        if structure.num_bonds == 0:
            return 0.0
        
        # Count bonds with reasonable lengths
        good_bonds = 0
        for bond in structure.bonds:
            atom1 = structure.atoms[bond.atom1_idx]
            atom2 = structure.atoms[bond.atom2_idx]
            element_pair = tuple(sorted([atom1.element, atom2.element]))
            
            if element_pair in self.BOND_LENGTH_RANGES:
                min_len, max_len = self.BOND_LENGTH_RANGES[element_pair]
                if min_len <= bond.bond_length <= max_len:
                    good_bonds += 1
        
        return good_bonds / structure.num_bonds
    
    def _compute_coordination_quality(self, structure: NanozymeStructure) -> float:
        """Compute metal coordination quality"""
        metal_atoms = structure.get_metal_atoms()
        
        if len(metal_atoms) == 0:
            return 1.0  # No metals to check
        
        scores = []
        for metal in metal_atoms:
            score = self._check_coordination_geometry(metal, structure)
            scores.append(score)
        
        return np.mean(scores)

