"""
Extract Catalytic Motifs from Nanozyme Enzymes using ChemEnzyRetroPlanner
==========================================================================

This script uses ChemEnzyRetroPlanner's EasIFA model to predict active sites
for nanozyme-related enzymes and extract catalytic motif information.

Usage:
    python extract_nanozyme_motifs.py --output_dir ./motif_results
"""

import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime

# Import default constants
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
try:
    from nanozyme_mining.utils.constants import DEFAULT_TOPK_ENZYMES
except ImportError:
    DEFAULT_TOPK_ENZYMES = 50  # Fallback default

# Add ChemEnzyRetroPlanner to path
CHEMENZYRP_PATH = "/home/tangboshi/.111tangboshi/参考项目代码库/ChemEnzyRetroPlanner-main"
sys.path.insert(0, CHEMENZYRP_PATH)
sys.path.insert(0, os.path.join(CHEMENZYRP_PATH, "webapp"))
sys.path.insert(0, os.path.join(CHEMENZYRP_PATH, "retro_planner"))

# Import ChemEnzyRetroPlanner modules
from retro_planner.api import dirpath
from retro_planner.packages.easifa.easifa.interface.easifa_inference_api import EasIFAInferenceAPI
from retro_planner.packages.uniprot_parser.uniprot_parser import UniProtParserEC

# Import visualization utilities
from images import get_structure_html_and_active_data


# Nanozyme EC numbers mapping
NANOZYME_EC_NUMBERS = {
    "POD": {
        "name": "Peroxidase",
        "ec_numbers": ["1.11.1.7", "1.11.1.11", "1.11.1.21"],
        "typical_reaction": "ROOH>>ROH.O"  # Peroxide reduction
    },
    "CAT": {
        "name": "Catalase",
        "ec_numbers": ["1.11.1.6"],
        "typical_reaction": "OO>>O=O.O"  # H2O2 decomposition
    },
    "SOD": {
        "name": "Superoxide Dismutase",
        "ec_numbers": ["1.15.1.1"],
        "typical_reaction": "[O-][O]>>[O-][O-]"  # Superoxide dismutation
    },
    "GSH": {
        "name": "Glutathione Peroxidase",
        "ec_numbers": ["1.11.1.9", "1.11.1.12"],
        "typical_reaction": "ROOH>>ROH.O"  # Similar to POD
    },
    "GOX": {
        "name": "Glucose Oxidase",
        "ec_numbers": ["1.1.3.4"],
        "typical_reaction": "OCC(O)C(O)C(O)C(O)C=O>>OCC(O)C(O)C(O)C(=O)C=O"  # Glucose oxidation
    },
    "LAC": {
        "name": "Laccase",
        "ec_numbers": ["1.10.3.2"],
        "typical_reaction": "Oc1ccccc1>>O=C1C=CC=CC1"  # Phenol oxidation
    },
    "OXD": {
        "name": "Oxidase",
        "ec_numbers": ["1.4.3.4", "1.3.3.4"],
        "typical_reaction": "CC(N)C(=O)O>>CC(=O)C(=O)O"  # Amine oxidation
    },
    "PHOS": {
        "name": "Phosphatase",
        "ec_numbers": ["3.1.3.1"],
        "typical_reaction": "OP(=O)(O)O>>O"  # Phosphate ester hydrolysis
    },
    "DNASE": {
        "name": "DNase",
        "ec_numbers": ["3.1.21.1"],
        "typical_reaction": "[DNA]>>[nucleotides]"  # DNA hydrolysis
    }
}


class NanozymeMotifExtractor:
    """Extract catalytic motifs from nanozyme enzymes."""

    def __init__(self,
                 easifa_checkpoint_path: str,
                 uniprot_data_path: str,
                 pdb_cache_path: str,
                 device: str = "cpu"):
        """
        Initialize the motif extractor.

        Args:
            easifa_checkpoint_path: Path to EasIFA model checkpoint
            uniprot_data_path: Path to UniProt data directory
            pdb_cache_path: Path to PDB cache directory
            device: Device to run model on (cpu or cuda)
        """
        print(f"Initializing NanozymeMotifExtractor...")
        print(f"  Device: {device}")
        print(f"  EasIFA checkpoint: {easifa_checkpoint_path}")

        # Initialize EasIFA predictor
        self.easifa_annotator = EasIFAInferenceAPI(
            model_checkpoint_path=easifa_checkpoint_path,
            device=device
        )

        # Initialize UniProt parser
        self.uniprot_parser = UniProtParserEC(
            json_folder=os.path.join(uniprot_data_path, "json"),
            csv_folder=os.path.join(uniprot_data_path, "csv"),
            alphafolddb_folder=pdb_cache_path,
            chebi_path=None,  # Not needed for motif extraction
            rxn_folder=None   # Not needed for motif extraction
        )

        print("Initialization complete!")

    def extract_motif_for_ec(self,
                            ec_number: str,
                            reaction_smiles: str = None,
                            topk: int = DEFAULT_TOPK_ENZYMES) -> List[Dict]:
        """
        Extract catalytic motif for a specific EC number.

        Args:
            ec_number: EC number (e.g., "1.11.1.7")
            reaction_smiles: Optional reaction SMILES for context
            topk: Number of top enzymes to retrieve (default: 50)

        Returns:
            List of motif data dictionaries
        """
        print(f"\n{'='*60}")
        print(f"Processing EC {ec_number}")
        print(f"{'='*60}")

        results = []

        try:
            # Query UniProt for PDB structures
            if reaction_smiles:
                uniprot_df = self.uniprot_parser.query_enzyme_pdb_by_ec_with_rxn_ranking(
                    ec_number=ec_number,
                    rxn_smiles=reaction_smiles,
                    topk=topk
                )
            else:
                uniprot_df = self.uniprot_parser.query_enzyme_pdb_by_ec(
                    ec_number=ec_number,
                    topk=topk
                )

            if uniprot_df.empty:
                print(f"  ⚠️  No PDB structures found for EC {ec_number}")
                return results

            print(f"  ✓ Found {len(uniprot_df)} PDB structures")

            # Process each enzyme
            for idx, row in uniprot_df.iterrows():
                alphafolddb_id = row.get("AlphaFoldDB", "Unknown")
                uniprot_id = row.get("Entry", "Unknown")
                pdb_fpath = row.get("pdb_fpath", None)

                if not pdb_fpath or not os.path.exists(pdb_fpath):
                    print(f"    ⚠️  PDB file not found: {pdb_fpath}")
                    continue

                print(f"    Processing {uniprot_id} ({alphafolddb_id})...")

                try:
                    # Predict active sites using EasIFA
                    if reaction_smiles:
                        predicted_labels = self.easifa_annotator.inference(
                            rxn=reaction_smiles,
                            enzyme_structure_path=pdb_fpath
                        )
                    else:
                        # Use structure-only prediction
                        predicted_labels = self.easifa_annotator.inference(
                            rxn="",  # Empty reaction
                            enzyme_structure_path=pdb_fpath
                        )

                    # Extract active site residues
                    active_sites = self._extract_active_sites(
                        pdb_fpath=pdb_fpath,
                        site_labels=predicted_labels
                    )

                    # Store results
                    motif_data = {
                        "ec_number": ec_number,
                        "uniprot_id": uniprot_id,
                        "alphafolddb_id": alphafolddb_id,
                        "pdb_path": pdb_fpath,
                        "active_sites": active_sites,
                        "num_binding_sites": sum(1 for s in active_sites if s["type"] == "Binding"),
                        "num_catalytic_sites": sum(1 for s in active_sites if s["type"] == "Catalytic"),
                        "timestamp": datetime.now().isoformat()
                    }

                    results.append(motif_data)
                    print(f"      ✓ Found {len(active_sites)} active sites")

                except Exception as e:
                    print(f"      ✗ Error predicting active sites: {e}")
                    continue

        except Exception as e:
            print(f"  ✗ Error processing EC {ec_number}: {e}")

        return results

    def _extract_active_sites(self,
                             pdb_fpath: str,
                             site_labels: List[int]) -> List[Dict]:
        """
        Extract active site information from PDB file and labels.

        Args:
            pdb_fpath: Path to PDB file
            site_labels: List of site type labels (0=None, 1=Binding, 2=Catalytic, 3=Other)

        Returns:
            List of active site dictionaries
        """
        active_sites = []

        # Label type mapping
        LABEL_TO_TYPE = {
            0: None,
            1: "Binding",
            2: "Catalytic",
            3: "Other"
        }

        try:
            with open(pdb_fpath, 'r') as f:
                lines = f.readlines()

            res_idx = None
            first_res_idx = None

            for line in lines:
                if not line.startswith("ATOM"):
                    continue

                # Parse PDB line
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                residue_number = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Track residue index
                if first_res_idx is None:
                    first_res_idx = residue_number
                res_idx = residue_number - first_res_idx

                # Only process CA atoms for active sites
                if atom_name == "CA" and res_idx < len(site_labels):
                    site_type_label = site_labels[res_idx]
                    site_type = LABEL_TO_TYPE.get(site_type_label)

                    if site_type:  # Skip None types
                        active_sites.append({
                            "residue_name": residue_name,
                            "residue_number": residue_number,
                            "chain_id": chain_id,
                            "type": site_type,
                            "coordinates": {"x": x, "y": y, "z": z}
                        })

        except Exception as e:
            print(f"      ✗ Error extracting active sites: {e}")

        return active_sites

    def extract_all_nanozyme_motifs(self,
                                    output_dir: str,
                                    topk: int = DEFAULT_TOPK_ENZYMES) -> Dict:
        """
        Extract motifs for all nanozyme types.

        Args:
            output_dir: Directory to save results
            topk: Number of top enzymes per EC number (default: 50)

        Returns:
            Dictionary of all results
        """
        print("\n" + "="*80)
        print("NANOZYME MOTIF EXTRACTION")
        print("="*80)

        os.makedirs(output_dir, exist_ok=True)
        all_results = {}

        for nanozyme_type, info in NANOZYME_EC_NUMBERS.items():
            print(f"\n{'#'*80}")
            print(f"# {nanozyme_type}: {info['name']}")
            print(f"{'#'*80}")

            type_results = []

            for ec_number in info["ec_numbers"]:
                motifs = self.extract_motif_for_ec(
                    ec_number=ec_number,
                    reaction_smiles=info.get("typical_reaction"),
                    topk=topk
                )
                type_results.extend(motifs)

            all_results[nanozyme_type] = {
                "name": info["name"],
                "ec_numbers": info["ec_numbers"],
                "motifs": type_results,
                "total_enzymes": len(type_results)
            }

            # Save individual type results
            type_output_file = os.path.join(output_dir, f"{nanozyme_type}_motifs.json")
            with open(type_output_file, 'w') as f:
                json.dump(all_results[nanozyme_type], f, indent=2)
            print(f"\n  ✓ Saved {nanozyme_type} results to {type_output_file}")

        # Save combined results
        combined_output_file = os.path.join(output_dir, "all_nanozyme_motifs.json")
        with open(combined_output_file, 'w') as f:
            json.dump(all_results, f, indent=2)

        # Generate summary
        self._generate_summary(all_results, output_dir)

        print("\n" + "="*80)
        print("EXTRACTION COMPLETE!")
        print(f"Results saved to: {output_dir}")
        print("="*80)

        return all_results

    def _generate_summary(self, results: Dict, output_dir: str):
        """Generate summary statistics."""
        summary_data = []

        for nanozyme_type, data in results.items():
            summary_data.append({
                "Nanozyme Type": nanozyme_type,
                "Name": data["name"],
                "EC Numbers": ", ".join(data["ec_numbers"]),
                "Total Enzymes": data["total_enzymes"],
                "Total Binding Sites": sum(m["num_binding_sites"] for m in data["motifs"]),
                "Total Catalytic Sites": sum(m["num_catalytic_sites"] for m in data["motifs"])
            })

        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(output_dir, "summary.csv")
        summary_df.to_csv(summary_file, index=False)

        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(summary_df.to_string(index=False))
        print(f"\nSummary saved to: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description="Extract nanozyme catalytic motifs")
    parser.add_argument("--output_dir", type=str, default="./motif_results",
                       help="Output directory for results")
    parser.add_argument("--topk", type=int, default=DEFAULT_TOPK_ENZYMES,
                       help=f"Number of top enzymes per EC number (default: {DEFAULT_TOPK_ENZYMES})")
    parser.add_argument("--device", type=str, default="cpu",
                       help="Device to run model on (cpu or cuda)")

    args = parser.parse_args()

    # Paths configuration
    EASIFA_CHECKPOINT = os.path.join(
        CHEMENZYRP_PATH,
        "retro_planner/packages/easifa/checkpoints/enzyme_site_type_predition_model"
    )
    UNIPROT_DATA = os.path.join(CHEMENZYRP_PATH, "data/uniprot")
    PDB_CACHE = os.path.join(CHEMENZYRP_PATH, "data/pdb_cache")

    # Initialize extractor
    extractor = NanozymeMotifExtractor(
        easifa_checkpoint_path=EASIFA_CHECKPOINT,
        uniprot_data_path=UNIPROT_DATA,
        pdb_cache_path=PDB_CACHE,
        device=args.device
    )

    # Extract motifs
    results = extractor.extract_all_nanozyme_motifs(
        output_dir=args.output_dir,
        topk=args.topk
    )


if __name__ == "__main__":
    main()
