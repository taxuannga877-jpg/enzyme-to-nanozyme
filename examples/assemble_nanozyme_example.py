"""
Nanozyme Assembly Example
=========================

Demonstrates how to assemble nanozymes from catalytic motifs.

Usage:
    python examples/assemble_nanozyme_example.py
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from nanozyme_mining.assembly import (
    MaterialType,
    CoordinationType,
    NanozymeMotif,
    AnchorAtom,
    MetalProperties,
    GeometryConstraint,
    NanozymeAssembler,
    MotifLibrary,
    NanozymeValidator,
)


def create_example_motif() -> NanozymeMotif:
    """
    Create an example POD-like nanozyme motif
    
    Represents a Fe-based peroxidase active site
    """
    # Define anchor atoms
    anchor_atoms = [
        # Metal center (Fe)
        AnchorAtom(
            atom_name="FE1",
            element="Fe",
            coordinates=[0.0, 0.0, 0.0],
            is_metal_center=True,
            oxidation_state=3,
            coordination_number=6,
            role="metal_center",
        ),
        # Coordinating atoms (His-His-Asn triad-like)
        AnchorAtom(
            atom_name="N1",
            element="N",
            coordinates=[2.0, 0.0, 0.0],
            is_donor=True,
            role="coordinating_atom",
        ),
        AnchorAtom(
            atom_name="N2",
            element="N",
            coordinates=[0.0, 2.0, 0.0],
            is_donor=True,
            role="coordinating_atom",
        ),
        AnchorAtom(
            atom_name="O1",
            element="O",
            coordinates=[0.0, 0.0, 2.0],
            is_donor=True,
            role="coordinating_atom",
        ),
    ]
    
    # Define geometry constraints
    geometry_constraints = [
        # Fe-N distances
        GeometryConstraint(
            constraint_type="distance",
            atom_indices=[0, 1],
            value=2.0,
            tolerance=0.2,
        ),
        GeometryConstraint(
            constraint_type="distance",
            atom_indices=[0, 2],
            value=2.0,
            tolerance=0.2,
        ),
        # Fe-O distance
        GeometryConstraint(
            constraint_type="distance",
            atom_indices=[0, 3],
            value=2.0,
            tolerance=0.2,
        ),
    ]
    
    # Define metal properties
    metal_centers = [
        MetalProperties(
            element="Fe",
            oxidation_state=3,
            coordination_number=6,
            coordination_geometry=CoordinationType.OCTAHEDRAL,
            d_electron_count=5,
            spin_state="high_spin",
        )
    ]
    
    # Create motif
    motif = NanozymeMotif(
        motif_id="FePOD_example_001",
        nanozyme_type="POD",
        source_ec_number="1.11.1.7",
        anchor_atoms=anchor_atoms,
        geometry_constraints=geometry_constraints,
        material_type=MaterialType.METAL_OXIDE,
        metal_centers=metal_centers,
        chemistry_tag="Fe-based peroxidase mimic",
        reaction_smiles="ROOH>>ROH.O",
        confidence_score=0.9,
        extraction_method="manual",
        notes="Example Fe-POD motif for demonstration",
    )
    
    return motif


def example_1_single_atom_catalyst():
    """Example 1: Assemble a single-atom catalyst (Fe-N4-C)"""
    print("=" * 60)
    print("Example 1: Single-Atom Catalyst (Fe-N4-C)")
    print("=" * 60)
    
    # Create motif
    motif = create_example_motif()
    
    # Initialize assembler for single-atom catalyst
    assembler = NanozymeAssembler(
        strategy="rule",
        material_type=MaterialType.SAC,
        output_dir="./examples/output/sac",
    )
    
    # Assemble
    print("\nAssembling Fe-N4-C single-atom catalyst...")
    nanozyme = assembler.assemble(
        motifs=motif,
        num_metal_centers=1,
    )
    
    # Print info
    print(f"\n{nanozyme}")
    print(f"Composition: {nanozyme.compute_composition()}")
    print(f"Formula: {nanozyme.compute_formula()}")
    print(f"Active site atoms: {len(nanozyme.get_active_site_atoms())}")
    print(f"Metal atoms: {len(nanozyme.get_metal_atoms())}")
    
    # Validate
    validator = NanozymeValidator()
    result = validator.validate(nanozyme)
    print(f"\nValidation: {result}")
    if result.errors:
        print("Errors:")
        for error in result.errors:
            print(f"  - {error}")
    if result.warnings:
        print("Warnings:")
        for warning in result.warnings:
            print(f"  - {warning}")
    print(f"Scores: {result.scores}")
    
    # Save
    assembler.save_structure(nanozyme, "Fe_N4_C_sac", formats=['xyz', 'json'])
    print(f"\nSaved to: {assembler.output_dir}")
    
    return nanozyme


def example_2_metal_oxide_nanozyme():
    """Example 2: Assemble a metal oxide nanozyme"""
    print("\n" + "=" * 60)
    print("Example 2: Metal Oxide Nanozyme (Fe3O4-like)")
    print("=" * 60)
    
    # Create motif
    motif = create_example_motif()
    
    # Initialize assembler for metal oxide
    assembler = NanozymeAssembler(
        strategy="rule",
        material_type=MaterialType.METAL_OXIDE,
        output_dir="./examples/output/metal_oxide",
    )
    
    # Assemble
    print("\nAssembling Fe3O4-like nanozyme...")
    nanozyme = assembler.assemble(
        motifs=motif,
        num_metal_centers=3,
        size=20,
    )
    
    # Print info
    print(f"\n{nanozyme}")
    print(f"Composition: {nanozyme.compute_composition()}")
    print(f"Formula: {nanozyme.compute_formula()}")
    
    # Validate
    validator = NanozymeValidator()
    result = validator.validate(nanozyme)
    print(f"\nValidation: {result}")
    print(f"Scores: {result.scores}")
    
    # Save
    assembler.save_structure(nanozyme, "Fe3O4_nanozyme", formats=['xyz', 'json'])
    print(f"\nSaved to: {assembler.output_dir}")
    
    return nanozyme


def example_3_metal_complex():
    """Example 3: Assemble a metal complex nanozyme"""
    print("\n" + "=" * 60)
    print("Example 3: Metal Complex Nanozyme")
    print("=" * 60)
    
    # Create motif
    motif = create_example_motif()
    
    # Initialize assembler for metal complex
    assembler = NanozymeAssembler(
        strategy="rule",
        material_type=MaterialType.METAL_COMPLEX,
        output_dir="./examples/output/metal_complex",
    )
    
    # Assemble
    print("\nAssembling Fe-complex nanozyme...")
    nanozyme = assembler.assemble(
        motifs=motif,
    )
    
    # Print info
    print(f"\n{nanozyme}")
    print(f"Composition: {nanozyme.compute_composition()}")
    print(f"Formula: {nanozyme.compute_formula()}")
    
    # Validate
    validator = NanozymeValidator(strict=False)
    result = validator.validate(nanozyme)
    print(f"\nValidation: {result}")
    print(f"Scores: {result.scores}")
    
    # Save
    assembler.save_structure(nanozyme, "Fe_complex_nanozyme", formats=['xyz', 'json'])
    print(f"\nSaved to: {assembler.output_dir}")
    
    return nanozyme


def example_4_multiple_motifs():
    """Example 4: Assemble from multiple motifs"""
    print("\n" + "=" * 60)
    print("Example 4: Multi-site Nanozyme (Multiple Motifs)")
    print("=" * 60)
    
    # Create multiple motifs
    motif1 = create_example_motif()
    motif1.motif_id = "FePOD_site1"
    
    motif2 = create_example_motif()
    motif2.motif_id = "FePOD_site2"
    # Slightly different coordinates
    for atom in motif2.anchor_atoms:
        atom.coordinates = [c + 0.5 for c in atom.coordinates]
    
    motifs = [motif1, motif2]
    
    # Initialize assembler
    assembler = NanozymeAssembler(
        strategy="rule",
        material_type=MaterialType.METAL_OXIDE,
        output_dir="./examples/output/multi_site",
    )
    
    # Assemble
    print(f"\nAssembling nanozyme from {len(motifs)} motifs...")
    nanozyme = assembler.assemble(
        motifs=motifs,
        spacing=10.0,  # Distance between sites
    )
    
    # Print info
    print(f"\n{nanozyme}")
    print(f"Composition: {nanozyme.compute_composition()}")
    print(f"Formula: {nanozyme.compute_formula()}")
    print(f"Number of motifs: {len(nanozyme.motif_origins)}")
    
    # Validate
    validator = NanozymeValidator()
    result = validator.validate(nanozyme)
    print(f"\nValidation: {result}")
    print(f"Scores: {result.scores}")
    
    # Save
    assembler.save_structure(nanozyme, "multi_site_nanozyme", formats=['xyz', 'json'])
    print(f"\nSaved to: {assembler.output_dir}")
    
    return nanozyme


def example_5_from_library():
    """Example 5: Assemble from motif library"""
    print("\n" + "=" * 60)
    print("Example 5: Assemble from Motif Library")
    print("=" * 60)
    
    # Try to load motif library
    library = MotifLibrary(library_dir="./motif_library")
    
    if len(library) == 0:
        print("\nMotif library is empty. Creating example motif...")
        motif = create_example_motif()
        # Save to library
        library_dir = Path("./motif_library/POD")
        library_dir.mkdir(parents=True, exist_ok=True)
        motif.save(str(library_dir / f"{motif.motif_id}.json"))
        print(f"Saved example motif to library")
        motifs = [motif]
    else:
        print(f"\nLoaded motif library: {library}")
        # Get POD motifs
        motifs = library.get_by_type("POD")
        if not motifs:
            print("No POD motifs found in library")
            return
        print(f"Found {len(motifs)} POD motifs")
    
    # Assemble from first motif
    motif = motifs[0]
    print(f"\nUsing motif: {motif.motif_id}")
    
    assembler = NanozymeAssembler(
        strategy="rule",
        material_type=MaterialType.METAL_OXIDE,
        output_dir="./examples/output/from_library",
    )
    
    nanozyme = assembler.assemble(motifs=motif)
    
    print(f"\n{nanozyme}")
    print(f"Formula: {nanozyme.compute_formula()}")
    
    assembler.save_structure(nanozyme, "library_nanozyme", formats=['xyz', 'json'])
    print(f"\nSaved to: {assembler.output_dir}")
    
    return nanozyme


def main():
    """Run all examples"""
    print("\n" + "=" * 60)
    print("Nanozyme Assembly Examples")
    print("=" * 60)
    print("\nThese examples demonstrate how to assemble nanozyme structures")
    print("from catalytic motifs using different strategies and material types.")
    print("\nIMPORTANT: Nanozymes are NOT proteins!")
    print("They are artificial nanomaterials with enzyme-like activities.")
    print("=" * 60)
    
    # Run examples
    try:
        example_1_single_atom_catalyst()
        example_2_metal_oxide_nanozyme()
        example_3_metal_complex()
        example_4_multiple_motifs()
        example_5_from_library()
        
        print("\n" + "=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        print("\nCheck the ./examples/output/ directory for generated structures.")
        
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

