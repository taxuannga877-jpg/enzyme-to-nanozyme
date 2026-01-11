"""
Test Motif Extraction - Simplified Version
==========================================

Test script to verify motif extraction functionality with a single EC number.
"""

import os
import sys
import json
from pathlib import Path

# Add ChemEnzyRetroPlanner to path
CHEMENZYRP_PATH = "/home/tangboshi/.111tangboshi/参考项目代码库/ChemEnzyRetroPlanner-main"
sys.path.insert(0, CHEMENZYRP_PATH)
sys.path.insert(0, os.path.join(CHEMENZYRP_PATH, "webapp"))
sys.path.insert(0, os.path.join(CHEMENZYRP_PATH, "retro_planner"))

# Import ChemEnzyRetroPlanner modules
try:
    from retro_planner.packages.easifa.easifa.interface.easifa_inference_api import EasIFAInferenceAPI
    from retro_planner.packages.uniprot_parser.uniprot_parser import UniProtParserEC
    print("✓ Successfully imported ChemEnzyRetroPlanner modules")
except ImportError as e:
    print(f"✗ Error importing modules: {e}")
    sys.exit(1)

# Configuration
EASIFA_CHECKPOINT = "/home/tangboshi/.111tangboshi/data/models/models/easifa/checkpoints/enzyme_site_type_predition_model/train_in_uniprot_ecreact_cluster_split_merge_dataset_all_limit_3_at_2023-12-19-16-06-42/global_step_284000"
UNIPROT_DATA = "/home/tangboshi/.111tangboshi/data/uniprot"
PDB_CACHE = "/home/tangboshi/.111tangboshi/data/pdb_cache"

# Test EC number (Catalase)
TEST_EC = "1.11.1.6"
TEST_REACTION = "OO>>O=O.O"  # H2O2 decomposition

def test_easifa_initialization():
    """Test EasIFA model initialization."""
    print("\n" + "="*60)
    print("TEST 1: EasIFA Model Initialization")
    print("="*60)

    try:
        print(f"Loading model from: {EASIFA_CHECKPOINT}")
        easifa = EasIFAInferenceAPI(
            model_checkpoint_path=EASIFA_CHECKPOINT,
            device="cpu"
        )
        print("✓ EasIFA model initialized successfully")
        return easifa
    except Exception as e:
        print(f"✗ Error initializing EasIFA: {e}")
        return None

def test_uniprot_parser():
    """Test UniProt parser initialization."""
    print("\n" + "="*60)
    print("TEST 2: UniProt Parser Initialization")
    print("="*60)

    # Create directories if they don't exist
    os.makedirs(os.path.join(UNIPROT_DATA, "json"), exist_ok=True)
    os.makedirs(os.path.join(UNIPROT_DATA, "csv"), exist_ok=True)
    os.makedirs(PDB_CACHE, exist_ok=True)

    try:
        print(f"UniProt data path: {UNIPROT_DATA}")
        print(f"PDB cache path: {PDB_CACHE}")

        parser = UniProtParserEC(
            json_folder=os.path.join(UNIPROT_DATA, "json"),
            csv_folder=os.path.join(UNIPROT_DATA, "csv"),
            alphafolddb_folder=PDB_CACHE,
            chebi_path=None,
            rxn_folder=None
        )
        print("✓ UniProt parser initialized successfully")
        return parser
    except Exception as e:
        print(f"✗ Error initializing UniProt parser: {e}")
        return None

def test_query_uniprot(parser):
    """Test querying UniProt for a specific EC number."""
    print("\n" + "="*60)
    print("TEST 3: Query UniProt Database")
    print("="*60)

    try:
        print(f"Querying EC number: {TEST_EC}")

        # Try to query UniProt
        df = parser.query_enzyme_pdb_by_ec(
            ec_number=TEST_EC,
            topk=3
        )

        if df.empty:
            print(f"⚠️  No results found for EC {TEST_EC}")
            print("This might be because:")
            print("  1. UniProt data needs to be downloaded")
            print("  2. The EC number is not in the database")
            return None

        print(f"✓ Found {len(df)} enzymes")
        print("\nResults:")
        for idx, row in df.iterrows():
            print(f"  - {row.get('Entry', 'Unknown')}: {row.get('AlphaFoldDB', 'Unknown')}")

        return df
    except Exception as e:
        print(f"✗ Error querying UniProt: {e}")
        return None

def test_predict_active_sites(easifa, pdb_path):
    """Test active site prediction."""
    print("\n" + "="*60)
    print("TEST 4: Predict Active Sites")
    print("="*60)

    try:
        print(f"PDB file: {pdb_path}")
        print(f"Reaction: {TEST_REACTION}")

        if not os.path.exists(pdb_path):
            print(f"⚠️  PDB file not found: {pdb_path}")
            return None

        # Predict active sites
        labels = easifa.inference(
            rxn=TEST_REACTION,
            enzyme_structure_path=pdb_path
        )

        # Count site types
        num_binding = sum(1 for l in labels if l == 1)
        num_catalytic = sum(1 for l in labels if l == 2)
        num_other = sum(1 for l in labels if l == 3)

        print(f"✓ Prediction complete")
        print(f"  Total residues: {len(labels)}")
        print(f"  Binding sites: {num_binding}")
        print(f"  Catalytic sites: {num_catalytic}")
        print(f"  Other sites: {num_other}")

        return labels
    except Exception as e:
        print(f"✗ Error predicting active sites: {e}")
        return None

def main():
    """Run all tests."""
    print("\n" + "="*80)
    print("MOTIF EXTRACTION TEST SUITE")
    print("="*80)

    # Test 1: Initialize EasIFA
    easifa = test_easifa_initialization()
    if not easifa:
        print("\n✗ Test suite failed: Cannot initialize EasIFA")
        return

    # Test 2: Initialize UniProt parser
    parser = test_uniprot_parser()
    if not parser:
        print("\n✗ Test suite failed: Cannot initialize UniProt parser")
        return

    # Test 3: Query UniProt
    df = test_query_uniprot(parser)
    if df is None or df.empty:
        print("\n⚠️  Cannot proceed with prediction test (no PDB data)")
        print("\nTo download UniProt data, you may need to:")
        print("  1. Run the UniProt data download script")
        print("  2. Or use the ChemEnzyRetroPlanner web interface to cache data")
        return

    # Test 4: Predict active sites
    pdb_path = df.iloc[0].get("pdb_fpath")
    if pdb_path:
        labels = test_predict_active_sites(easifa, pdb_path)
        if labels:
            print("\n" + "="*80)
            print("✓ ALL TESTS PASSED!")
            print("="*80)
        else:
            print("\n⚠️  Prediction test failed")
    else:
        print("\n⚠️  No PDB path found in results")

    print("\nTest suite complete.")

if __name__ == "__main__":
    main()
