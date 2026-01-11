"""
Nanozyme Assembly Module
========================

Nanozyme structure assembly from catalytic motifs.
Inspired by DiffLinker, LigandDiff, and stk.

Key modules:
- motif_enhanced: Enhanced motif representation for nanozymes
- assembler: Main assembly engine
- strategies: Different assembly strategies (rule-based, diffusion, template)
- validator: Chemical validity checking
- predictor: Catalytic activity prediction
"""

from .motif_enhanced import (
    MaterialType,
    CoordinationType,
    NanozymeMotif,
    AnchorAtom,
    GeometryConstraint,
    MetalProperties,
    AtomType,
    convert_basic_to_nanozyme_motif,
)
from .assembler import NanozymeAssembler, MotifLibrary
from .validator import NanozymeValidator, ValidationResult
from .structure import NanozymeStructure, Atom, Bond

__all__ = [
    "MaterialType",
    "CoordinationType",
    "NanozymeMotif",
    "AnchorAtom",
    "GeometryConstraint",
    "MetalProperties",
    "AtomType",
    "convert_basic_to_nanozyme_motif",
    "NanozymeAssembler",
    "MotifLibrary",
    "NanozymeValidator",
    "ValidationResult",
    "NanozymeStructure",
    "Atom",
    "Bond",
]

