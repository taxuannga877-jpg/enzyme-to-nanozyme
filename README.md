# Nanozyme Mining System

<div align="center">

![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)

**A comprehensive system for mining catalytic motifs from enzyme structures and designing nanozyme materials**

[Features](#-features) • [Installation](#-installation) • [Quick Start](#-quick-start) • [Documentation](#-documentation) • [Contributing](#-contributing)

</div>

---

## 📖 Overview

**Nanozyme Mining System** is a three-stage pipeline for extracting catalytic motifs from natural enzyme structures and using them to design artificial nanozyme materials. The system bridges the gap between enzyme catalysis and nanozyme design by:

1. **Database Layer**: Mapping EC numbers to nanozyme function types and fetching enzyme data
2. **Extraction Layer**: Extracting catalytic motifs from PDB structures
3. **Assembly Layer**: Assembling nanozyme structures from extracted motifs

### Key Concept

> **Nanozymes are NOT proteins!** They are artificial nanomaterials (metal oxides, MOFs, single-atom catalysts, etc.) with enzyme-like activities. This system extracts catalytic patterns from natural enzymes and adapts them for nanozyme design.

---

## ✨ Features

### 🔍 Data Acquisition
- **UniProt Integration**: Automatic enzyme data fetching by EC number
- **AlphaFold Database**: Download high-quality protein structures
- **M-CSA Database**: Query metal sites and catalytic residue information
- **Dual-Track Processing**: Handle both annotated and unannotated data

### 🧬 Motif Extraction
- **Active Site Detection**: 
  - Direct annotation from UniProt/M-CSA
  - ML-based prediction using EasIFA model for unannotated data
- **Catalytic Motif Extraction**: Extract anchor atoms, geometry constraints, and chemical tags
- **Multi-format Support**: JSON, PDB, and custom motif formats

### 🏗️ Nanozyme Assembly
- **Multiple Material Types**:
  - Metal oxides (Fe₃O₄, CeO₂, MnO₂)
  - Single-atom catalysts (Fe-N₄-C, Cu-N₄)
  - Metal-organic frameworks (MOFs)
  - Carbon-based nanomaterials
  - Metal clusters
- **Assembly Strategies**: Rule-based, template-based (diffusion-based in development)
- **Structure Validation**: Chemical validity checking and geometry validation

### 🌐 Web Interface
- **Interactive Visualization**: 3D structure viewer with py3Dmol
- **EC Number Browser**: Search and explore enzymes by EC number
- **Motif Library**: Browse and manage extracted motifs
- **M-CSA Integration**: Visualize metal sites and coordination

### 🎯 Supported Nanozyme Types
- **POD** (Peroxidase-like, EC 1.11.1.7)
- **CAT** (Catalase-like, EC 1.11.1.6)
- **SOD** (Superoxide dismutase-like, EC 1.15.1.1)
- **GPx** (Glutathione peroxidase-like, EC 1.11.1.9)
- **Phosphatase** (EC 3.1.3.1)
- **DNase** (EC 3.1.21.1)

---

## 🚀 Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Basic Installation

```bash
# Clone the repository
git clone git@github.com:taxuannga877-jpg/enzyme-to-nanozyme.git
cd enzyme-to-nanozyme

# Install core dependencies
pip install -r requirements.txt
```

### Full Installation (with ML models)

For full functionality including EasIFA prediction:

```bash
# Install with all dependencies
pip install -r requirements.txt

# Install PyTorch (if not already installed)
# For CPU:
pip install torch torchvision torchaudio

# For GPU (CUDA 11.8):
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

### Optional Dependencies

```bash
# For advanced visualization
pip install matplotlib seaborn

# For molecular structure generation
pip install rdkit-pypi  # or conda install -c conda-forge rdkit
```

---

## 📦 Quick Start

### 1. Basic Usage: Extract Motifs from EC Number

```python
from nanozyme_mining.database import UniProtFetcher
from nanozyme_mining.extraction import MotifExtractor

# Initialize fetcher
fetcher = UniProtFetcher(cache_dir="./cache")

# Fetch enzyme data for POD (EC 1.11.1.7)
annotated, unannotated = fetcher.fetch_and_classify(
    ec_number="1.11.1.7",
    nanozyme_type="POD",
    max_results=10
)

# Extract motifs
extractor = MotifExtractor(output_dir="./motifs")
for entry in annotated:
    motif = extractor.extract_motif(
        pdb_path=entry["pdb_path"],
        uniprot_id=entry["uniprot_id"],
        ec_number="1.11.1.7",
        nanozyme_type="POD",
        active_site_indices=entry["active_site_indices"]
    )
    print(f"Extracted motif: {motif.motif_id}")
```

### 2. Using Dual-Track Processor

```python
from nanozyme_mining.core import DualTrackProcessor

# Initialize processor (handles both annotated and unannotated data)
processor = DualTrackProcessor(
    output_dir="./processed",
    device="cpu"  # or "cuda" for GPU
)

# Process enzyme entries
results = processor.process_batch(
    entries=enzyme_entries,
    reaction_smiles="C>>C"
)
```

### 3. Assembling Nanozyme Structures

```python
from nanozyme_mining.assembly import (
    NanozymeAssembler,
    MaterialType,
    convert_basic_to_nanozyme_motif
)

# Convert basic motif to nanozyme motif
nanozyme_motif = convert_basic_to_nanozyme_motif(
    basic_motif=motif,
    material_type=MaterialType.METAL_OXIDE
)

# Initialize assembler
assembler = NanozymeAssembler(
    strategy="rule",
    material_type=MaterialType.METAL_OXIDE,
    output_dir="./output"
)

# Assemble nanozyme
nanozyme = assembler.assemble(
    motifs=nanozyme_motif,
    num_metal_centers=3
)

# Save structure
assembler.save_structure(nanozyme, "Fe3O4_nanozyme", formats=['xyz', 'json'])
```

### 4. Running the Web Interface

```bash
# Navigate to enzyme_viewer directory
cd enzyme_viewer

# Run the Flask app
python app.py

# Or use the provided script
bash run.sh
```

Then open your browser to `http://localhost:5000`

---

## 📚 Documentation

### Project Structure

```
nanozyme_mining/
├── nanozyme_mining/          # Core package
│   ├── database/             # Data acquisition (UniProt, M-CSA)
│   ├── extraction/           # Motif extraction
│   ├── prediction/           # EasIFA model integration
│   ├── core/                 # Dual-track processor
│   ├── assembly/             # Nanozyme assembly
│   └── utils/                # Utilities and constants
├── enzyme_viewer/            # Web interface
├── scripts/                  # Batch processing scripts
├── docs/                     # Documentation
│   ├── ASSEMBLY_GUIDE.md    # Assembly guide
│   ├── MIGRATION_PLAN.md    # Migration from reference projects
│   └── IMPLEMENTATION_SUMMARY.md
├── examples/                 # Example scripts
├── cache/                    # Cached data
├── pdb_library/             # PDB structure library
└── motif_library/           # Extracted motif library
```

### Key Modules

#### Database Layer
- **`UniProtFetcher`**: Fetches enzyme data from UniProt API
- **`MCSAFetcher`**: Queries M-CSA database for metal sites
- **`NanozymeDatabase`**: Manages enzyme entry database

#### Extraction Layer
- **`MotifExtractor`**: Extracts catalytic motifs from PDB structures
- **`CatalyticMotif`**: Data structure for motifs (anchor atoms, geometry constraints)

#### Assembly Layer
- **`NanozymeAssembler`**: Main assembly engine
- **`NanozymeValidator`**: Structure validation
- **`NanozymeMotif`**: Enhanced motif representation for nanozymes

### Detailed Documentation

- **[Assembly Guide](docs/ASSEMBLY_GUIDE.md)**: Comprehensive guide for nanozyme assembly
- **[Migration Plan](docs/MIGRATION_PLAN.md)**: Migration from reference projects (DiffLinker, LigandDiff, stk)
- **[Implementation Summary](docs/IMPLEMENTATION_SUMMARY.md)**: Detailed implementation notes

---

## 🔧 Usage Examples

### Example 1: Batch Processing EC Numbers

```bash
# Process a single EC number
python scripts/run_pipeline.py --ec 1.11.1.7 --max_results 100

# Process all supported EC numbers
python scripts/run_pipeline.py --all --max_results 50
```

### Example 2: Extract Motifs from Cache

```bash
# Extract motifs from cached JSON data
python scripts/extract_motifs_from_cache.py --ec 1.11.1.7
```

### Example 3: Predict Active Sites for Unannotated Data

```bash
# Use EasIFA to predict active sites
python scripts/predict_unannotated_pdb.py --ec 1.11.1.7
```

See the `examples/` directory for more detailed examples.

---

## 🧪 Testing

```bash
# Run tests
python -m pytest tests/

# Run specific test module
python -m pytest tests/test_extraction.py
```

---

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Setup

```bash
# Clone your fork
git clone git@github.com:taxuannga877-jpg/enzyme-to-nanozyme.git
cd enzyme-to-nanozyme

# Install in development mode
pip install -e .

# Install development dependencies
pip install -r requirements-dev.txt  # if available
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

### Reference Projects

This project is inspired by and adapts concepts from:

- **[DiffLinker](https://github.com/gcorso/DiffLinker)**: Fragment linking with diffusion models
- **[LigandDiff](https://github.com/THUNLP-MT/LigandDiff)**: Metal complex generation
- **[stk](https://github.com/lukasturcani/stk)**: Template-based molecular assembly
- **[ChemEnzyRetroPlanner](https://github.com/yourusername/ChemEnzyRetroPlanner)**: Enzyme retro-planning framework

### Data Sources

- **UniProt**: Enzyme sequence and annotation data
- **AlphaFold Database**: Protein structure predictions
- **M-CSA**: Metal-containing active sites database

### Models

- **EasIFA**: Active site prediction model (integrated from ChemEnzyRetroPlanner)

---

## 📧 Contact

For questions, issues, or suggestions:

- **Issues**: [GitHub Issues](https://github.com/taxuannga877-jpg/enzyme-to-nanozyme/issues)

---

## 🗺️ Roadmap

- [ ] Diffusion-based assembly strategy
- [ ] Template-based assembly from MOF databases
- [ ] Catalytic activity prediction
- [ ] Multi-motif assembly
- [ ] Web-based assembly interface
- [ ] Database integration for motif sharing

---

## 📊 Citation

If you use this project in your research, please cite:

```bibtex
@software{nanozyme_mining,
  title = {Nanozyme Mining System: From Enzyme Structures to Nanozyme Design},
  author = {Nanozyme Design Team},
  year = {2024},
  url = {https://github.com/taxuannga877-jpg/enzyme-to-nanozyme},
  version = {0.2.0}
}
```

---

<div align="center">

**Made with ❤️ by the Nanozyme Design Team**

⭐ Star this repo if you find it useful!

</div>

