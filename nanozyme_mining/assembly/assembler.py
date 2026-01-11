"""
Nanozyme Assembler - Main Assembly Engine
==========================================

Assembles nanozyme structures from catalytic motifs using different strategies:
1. Rule-based: Using chemical rules and templates
2. Diffusion-based: Using diffusion models (like DiffLinker)
3. Template-based: Using predefined templates (like stk)

Inspired by:
- DiffLinker: fragment linking with diffusion models
- LigandDiff: metal complex generation
- stk: template-based molecular assembly
"""

from typing import List, Dict, Optional, Union
from pathlib import Path
import numpy as np

from .motif_enhanced import NanozymeMotif, MaterialType
from .structure import NanozymeStructure
from .strategies.rule_based import RuleBasedAssembler
# from .strategies.diffusion_based import DiffusionBasedAssembler  # TODO
# from .strategies.template_based import TemplateBasedAssembler  # TODO


class NanozymeAssembler:
    """
    Main nanozyme assembly engine
    
    Coordinates different assembly strategies and provides
    a unified interface for nanozyme structure generation.
    """
    
    def __init__(
        self,
        strategy: str = "rule",  # "rule", "diffusion", "template"
        material_type: Optional[MaterialType] = None,
        device: str = "cpu",
        output_dir: Optional[str] = None,
    ):
        """
        Initialize assembler
        
        Args:
            strategy: Assembly strategy to use
            material_type: Default material type for nanozymes
            device: Device for computation (cpu/cuda)
            output_dir: Directory to save outputs
        """
        self.strategy = strategy
        self.material_type = material_type or MaterialType.METAL_OXIDE
        self.device = device
        self.output_dir = Path(output_dir) if output_dir else Path("./assembled_nanozymes")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize strategy-specific assembler
        self._init_strategy()
    
    def _init_strategy(self):
        """Initialize the selected assembly strategy"""
        if self.strategy == "rule":
            self.assembler = RuleBasedAssembler(
                material_type=self.material_type,
                output_dir=str(self.output_dir),
            )
        elif self.strategy == "diffusion":
            # TODO: Implement diffusion-based assembler
            raise NotImplementedError("Diffusion-based assembly not yet implemented")
        elif self.strategy == "template":
            # TODO: Implement template-based assembler
            raise NotImplementedError("Template-based assembly not yet implemented")
        else:
            raise ValueError(f"Unknown strategy: {self.strategy}")
    
    def assemble(
        self,
        motifs: Union[NanozymeMotif, List[NanozymeMotif]],
        **kwargs
    ) -> NanozymeStructure:
        """
        Assemble nanozyme from motif(s)
        
        Args:
            motifs: Single motif or list of motifs
            **kwargs: Strategy-specific parameters
        
        Returns:
            Assembled nanozyme structure
        """
        # Ensure motifs is a list
        if isinstance(motifs, NanozymeMotif):
            motifs = [motifs]
        
        # Validate motifs
        self._validate_motifs(motifs)
        
        # Delegate to strategy-specific assembler
        structure = self.assembler.assemble(motifs, **kwargs)
        
        # Post-processing
        structure = self._post_process(structure)
        
        return structure
    
    def assemble_batch(
        self,
        motif_groups: List[List[NanozymeMotif]],
        **kwargs
    ) -> List[NanozymeStructure]:
        """
        Assemble multiple nanozymes in batch
        
        Args:
            motif_groups: List of motif lists
            **kwargs: Strategy-specific parameters
        
        Returns:
            List of assembled structures
        """
        structures = []
        for i, motifs in enumerate(motif_groups):
            print(f"Assembling nanozyme {i+1}/{len(motif_groups)}...")
            structure = self.assemble(motifs, **kwargs)
            structures.append(structure)
        return structures
    
    def _validate_motifs(self, motifs: List[NanozymeMotif]):
        """Validate input motifs"""
        if not motifs:
            raise ValueError("No motifs provided")
        
        # Check if all motifs have the same nanozyme type
        nanozyme_types = set(m.nanozyme_type for m in motifs)
        if len(nanozyme_types) > 1:
            print(f"Warning: Multiple nanozyme types detected: {nanozyme_types}")
        
        # Check if motifs have anchor atoms
        for motif in motifs:
            if not motif.anchor_atoms:
                print(f"Warning: Motif {motif.motif_id} has no anchor atoms")
    
    def _post_process(self, structure: NanozymeStructure) -> NanozymeStructure:
        """Post-process assembled structure"""
        # Add metadata
        structure.metadata['assembly_strategy'] = self.strategy
        structure.metadata['material_type'] = self.material_type.value
        
        # Compute properties
        structure.metadata['composition'] = structure.compute_composition()
        structure.metadata['formula'] = structure.compute_formula()
        structure.metadata['num_active_sites'] = len(structure.get_active_site_atoms())
        structure.metadata['num_metals'] = len(structure.get_metal_atoms())
        
        return structure
    
    def save_structure(
        self,
        structure: NanozymeStructure,
        name: str,
        formats: List[str] = ['xyz', 'json'],
    ):
        """
        Save structure in multiple formats
        
        Args:
            structure: Nanozyme structure
            name: Base filename
            formats: List of formats ('xyz', 'json', 'pdb')
        """
        for fmt in formats:
            if fmt == 'xyz':
                filepath = self.output_dir / f"{name}.xyz"
                structure.to_xyz(str(filepath), comment=f"{structure.nanozyme_type} nanozyme")
            elif fmt == 'json':
                filepath = self.output_dir / f"{name}.json"
                structure.to_json(str(filepath))
            elif fmt == 'pdb':
                # TODO: Implement PDB export
                print(f"Warning: PDB export not yet implemented")
            else:
                print(f"Warning: Unknown format {fmt}")
        
        print(f"Saved structure to {self.output_dir}/{name}")


class MotifLibrary:
    """
    Library of catalytic motifs
    
    Manages loading and querying motifs for assembly
    """
    
    def __init__(self, library_dir: str = "./motif_library"):
        self.library_dir = Path(library_dir)
        self.motifs: Dict[str, NanozymeMotif] = {}
        self._load_library()
    
    def _load_library(self):
        """Load all motifs from library directory"""
        if not self.library_dir.exists():
            print(f"Warning: Motif library directory not found: {self.library_dir}")
            return
        
        # Load motifs from each nanozyme type directory
        for nanozyme_dir in self.library_dir.iterdir():
            if nanozyme_dir.is_dir():
                self._load_motifs_from_dir(nanozyme_dir)
        
        print(f"Loaded {len(self.motifs)} motifs from library")
    
    def _load_motifs_from_dir(self, directory: Path):
        """Load motifs from a specific directory"""
        for motif_file in directory.glob("*.json"):
            try:
                motif = NanozymeMotif.load(str(motif_file))
                self.motifs[motif.motif_id] = motif
            except Exception as e:
                print(f"Warning: Failed to load {motif_file}: {e}")
    
    def get_by_id(self, motif_id: str) -> Optional[NanozymeMotif]:
        """Get motif by ID"""
        return self.motifs.get(motif_id)
    
    def get_by_type(self, nanozyme_type: str) -> List[NanozymeMotif]:
        """Get all motifs of a specific nanozyme type"""
        return [m for m in self.motifs.values() if m.nanozyme_type == nanozyme_type]
    
    def get_by_material(self, material_type: MaterialType) -> List[NanozymeMotif]:
        """Get all motifs of a specific material type"""
        return [m for m in self.motifs.values() if m.material_type == material_type]
    
    def search(
        self,
        nanozyme_type: Optional[str] = None,
        material_type: Optional[MaterialType] = None,
        min_confidence: float = 0.0,
    ) -> List[NanozymeMotif]:
        """Search motifs by criteria"""
        results = list(self.motifs.values())
        
        if nanozyme_type:
            results = [m for m in results if m.nanozyme_type == nanozyme_type]
        
        if material_type:
            results = [m for m in results if m.material_type == material_type]
        
        if min_confidence > 0:
            results = [m for m in results if m.confidence_score >= min_confidence]
        
        return results
    
    def __len__(self) -> int:
        return len(self.motifs)
    
    def __repr__(self) -> str:
        types = {}
        for motif in self.motifs.values():
            types[motif.nanozyme_type] = types.get(motif.nanozyme_type, 0) + 1
        type_str = ", ".join(f"{k}: {v}" for k, v in types.items())
        return f"MotifLibrary({len(self)} motifs, types: {type_str})"

