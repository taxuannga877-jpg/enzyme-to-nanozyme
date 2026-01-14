# Enzyme-to-Nanozyme Mining System

A comprehensive system for mining catalytic motifs from enzymes and mapping them to nanozyme design. This project extracts catalytic sites from enzyme structures based on EC numbers, predicts active sites using the EasIFA model, and organizes the data for nanozyme design applications.

## Features

- **EC Number to Nanozyme Type Mapping**: Automatically maps Enzyme Commission (EC) numbers to nanozyme function types (Peroxidase, Catalase, SOD, etc.)
- **Dual-Track Processing**: 
  - **Annotated Data**: Directly uses UniProt/M-CSA active site annotations
  - **Unannotated Data**: Uses EasIFA model to predict active sites
- **Catalytic Motif Extraction**: Extracts catalytic motifs from PDB structures with full residue information
- **Web Interface**: Interactive web viewer for enzyme structures and catalytic sites
- **Organized Data Storage**: Structured data repositories for PDB structures and extracted motifs

## Installation

### System Requirements

- Python 3.8 or higher
- Linux, macOS, or Windows (WSL2 recommended for Windows)

### Install Dependencies

```bash
# Clone the repository
git clone https://github.com/taxuannga877-jpg/enzyme-to-nanozyme.git
cd enzyme-to-nanozyme

# Install core dependencies
pip install -r requirements.txt

# For EasIFA model support (optional, for active site prediction)
pip install -r requirements.txt  # Includes PyTorch, DGL, dgllife, fair-esm
```

### Optional: Install in Development Mode

```bash
pip install -e .
```

## Quick Start

### 1. Process an EC Number

Process a single EC number to fetch enzyme data, download PDB structures, and extract motifs:

```bash
python scripts/run_pipeline.py --ec 1.11.1.7 --max_results 100
```

### 2. Process All Supported EC Numbers

```bash
python scripts/run_pipeline.py --all --max_results 50
```

### 3. Start Web Interface

```bash
cd enzyme_viewer
python app.py
```

Then open your browser to `http://localhost:5000`

### 4. Extract Motifs from Cache

If you already have cached data, extract motifs directly:

```bash
python scripts/extract_motifs_from_cache.py --ec 1.11.1.7
```

## Data Repositories

The system uses two main data directories that are generated at runtime:

### `pdb_library/`

PDB structure library organized by EC numbers. Contains:
- PDB files (AlphaFold structures and experimental PDBs)
- JSON index files (`{EC_number}_sites.json`) with enzyme metadata and active site annotations

**Structure:**
```
pdb_library/
├── library_index.json
└── {EC_number}/
    ├── AF-{uniprot_id}-F1-model_v6.pdb
    ├── {pdb_id}.pdb  (experimental structures)
    └── {EC_number}_sites.json
```

**Note**: This directory is generated at runtime and should be added to `.gitignore`. It can be large (several GB depending on the number of structures downloaded).

### `motif_library/`

Motif repository organized by EC numbers and motif types. Contains extracted catalytic motifs classified by type.

**Structure:**
```
motif_library/
└── {EC_number}/
    ├── catalytic_sites/
    │   └── {uniprot_id}_{site_id}.json
    ├── binding_sites/
    │   └── {uniprot_id}_{site_id}.json
    ├── metal_sites/
    │   └── {uniprot_id}_{site_id}.json
    └── other/
        └── {uniprot_id}_{site_id}.json
```

Each motif JSON file contains:
- Residue information (name, number, coordinates)
- Atom coordinates
- Functional roles
- Extraction source (UniProt, M-CSA, or EasIFA prediction)

**Note**: This directory is also generated at runtime and should be added to `.gitignore`.

## Project Structure

```
enzyme-to-nanozyme/
├── nanozyme_mining/          # Core mining module
│   ├── database/             # Database layer (UniProt, M-CSA)
│   ├── extraction/           # Motif extraction
│   ├── prediction/          # EasIFA active site prediction
│   ├── structure/           # PDB parsing and structure handling
│   ├── core/                # Dual-track processor
│   └── utils/               # Utilities and constants
├── enzyme_viewer/           # Web interface
│   ├── app.py              # Flask application
│   ├── templates/          # HTML templates
│   └── motif_db.py        # Motif database indexing
├── scripts/                # Utility scripts
│   ├── run_pipeline.py    # Main pipeline script
│   ├── extract_motifs_from_cache.py
│   └── init_local_data.py
├── pdb_library/            # PDB structure library (generated)
├── motif_library/         # Motif repository (generated)
├── cache/                 # Temporary cache (generated)
├── requirements.txt       # Python dependencies
├── setup.py              # Package setup
└── README.md             # This file
```

## Usage

### Command-Line Tools

#### Process EC Numbers

```bash
# Process a single EC number
python scripts/run_pipeline.py --ec 1.11.1.7 --max_results 100

# Process all supported EC numbers
python scripts/run_pipeline.py --all --max_results 50

# Skip EasIFA prediction (only process annotated data)
python scripts/run_pipeline.py --ec 1.11.1.7 --skip-prediction
```

#### Extract Motifs

```bash
# Extract motifs for a specific EC number
python scripts/extract_motifs_from_cache.py --ec 1.11.1.7

# Extract all motifs and rebuild index
python scripts/extract_motifs_from_cache.py --clear
```

### Python API

```python
from nanozyme_mining import (
    UniProtFetcher,
    EasIFAPredictor,
    MotifExtractor,
    DualTrackProcessor
)

# Fetch enzyme data
fetcher = UniProtFetcher(cache_dir="./cache")
annotated, unannotated = fetcher.fetch_and_classify(
    ec_number="1.11.1.7",
    nanozyme_type="POD",
    max_results=100
)

# Predict active sites for unannotated data
processor = DualTrackProcessor(output_dir="./processed")
results = processor.predict_unannotated_batch(
    unannotated,
    reaction_smiles="C>>C"
)

# Extract motifs
extractor = MotifExtractor(output_dir="./motifs")
motif = extractor.extract_motif(
    pdb_path="path/to/structure.pdb",
    uniprot_id="P12345",
    ec_number="1.11.1.7",
    nanozyme_type="POD",
    active_site_indices=[50, 51, 52]
)
```

## Supported EC Numbers

The system currently supports the following nanozyme types and their corresponding EC numbers:

- **Peroxidases (POD)**: 1.11.1.7, 1.11.1.11, 1.11.1.21
- **Catalases (CAT)**: 1.11.1.6, 1.11.1.21
- **Superoxide Dismutase (SOD)**: 1.15.1.1
- **Glutathione Peroxidases (GSH)**: 1.11.1.9, 1.11.1.12
- **Oxidases (OXD)**: 1.4.3.4, 1.3.3.4
- **Laccase (LAC)**: 1.10.3.2
- **Glucose Oxidase (GOX)**: 1.1.3.4
- **Phosphatase (PHOS)**: 3.1.3.1
- **DNase (DNASE)**: 3.1.21.1

## Configuration

### EasIFA Model Path

The EasIFA model checkpoint path is configured in `nanozyme_mining/prediction/easifa_predictor.py`. By default, it looks for:

```
data/models/models/easifa/checkpoints/enzyme_site_type_predition_model/
train_in_uniprot_ecreact_cluster_split_merge_dataset_all_limit_3_at_2023-12-19-16-06-42/
global_step_284000
```

You can override this by setting the `ENZYME_MODEL_PATH` environment variable or modifying the code.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines and contribution instructions.

### Running Tests

```bash
pytest tests/
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **ChemEnzyRetroPlanner**: Architecture inspiration and EasIFA model
- **UniProt**: Enzyme data and annotations
- **M-CSA**: Catalytic Site Atlas for active site annotations
- **AlphaFold Database**: Protein structure predictions
- **DGL**: Deep Graph Library for graph neural networks
- **RDKit**: Chemical informatics toolkit

## Citation

If you use this project in your research, please cite:

```bibtex
@software{enzyme_to_nanozyme,
  title = {Enzyme-to-Nanozyme Mining System},
  author = {Nanozyme Design Team},
  year = {2024},
  url = {https://github.com/taxuannga877-jpg/enzyme-to-nanozyme}
}
```

## Contact

For questions, issues, or contributions, please open an issue on GitHub.

