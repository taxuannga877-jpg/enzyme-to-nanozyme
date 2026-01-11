"""
EC Number Mappings for Nanozyme Types
======================================

Based on literature and IUBMB nomenclature.
"""

from typing import Dict, List
from .constants import NanozymeType


# EC_PATTERNS: NanozymeType -> List of EC numbers
# Used by database module for batch fetching
EC_PATTERNS: Dict[NanozymeType, List[str]] = {
    NanozymeType.POD: ["1.11.1.7", "1.11.1.11", "1.11.1.21"],  # Peroxidases
    NanozymeType.CAT: ["1.11.1.6", "1.11.1.21"],               # Catalases
    NanozymeType.SOD: ["1.15.1.1"],                            # Superoxide dismutase
    NanozymeType.GSH: ["1.11.1.9", "1.11.1.12"],               # Glutathione peroxidases
    NanozymeType.OXD: ["1.4.3.4", "1.3.3.4"],                  # Oxidases
    NanozymeType.LAC: ["1.10.3.2"],                            # Laccase
    NanozymeType.GOX: ["1.1.3.4"],                             # Glucose oxidase
    NanozymeType.HRP: ["1.11.1.7"],                            # Horseradish peroxidase
    NanozymeType.PHOS: ["3.1.3.1"],                            # Phosphatase
    NanozymeType.DNASE: ["3.1.21.1"],                          # DNase
}


# Core EC mappings for nanozyme database
NANOZYME_EC_MAPPINGS: Dict[str, Dict] = {
    "POD": {
        "name": "POD-like (过氧化物酶样)",
        "primary_ec": "1.11.1.7",
        "description": "经典 TMB/ABTS/OPD + H₂O₂ 体系",
        "reaction": "donor + H2O2 → oxidized donor + H2O",
    },
    "CAT": {
        "name": "CAT-like (过氧化氢酶样)",
        "primary_ec": "1.11.1.6",
        "description": "2 H₂O₂ → O₂ + 2 H₂O，ROS scavenger",
        "reaction": "2 H2O2 → O2 + 2 H2O",
    },
    "KatG": {
        "name": "CAT-POD双功能 (KatG-like)",
        "primary_ec": "1.11.1.21",
        "description": "同一活性中心同时表现强CAT和POD特征",
        "reaction": "Bifunctional CAT + POD",
    },
    "SOD": {
        "name": "SOD-like (超氧化物歧化酶样)",
        "primary_ec": "1.15.1.1",
        "description": "O₂•⁻ 的歧化反应，抗氧化nanozyme",
        "reaction": "2 O2•- + 2 H+ → O2 + H2O2",
    },
    "GPx": {
        "name": "GPx-like (谷胱甘肽过氧化物酶样)",
        "primary_ec": "1.11.1.9",
        "description": "以GSH为还原底物、清除H₂O₂/ROOH",
        "reaction": "2 GSH + H2O2 → GSSG + 2 H2O",
    },
    "Phosphatase": {
        "name": "Phosphatase-like (磷酸酶样)",
        "primary_ec": "3.1.3.1",
        "description": "水解磷酸单酯，释放无机磷酸",
        "reaction": "phosphate monoester + H2O → alcohol + phosphate",
    },
    "DNase": {
        "name": "DNase-like (脱氧核酸酶样)",
        "primary_ec": "3.1.21.1",
        "description": "水解DNA链中的磷酸二酯键，产生单核苷酸或寡核苷酸",
        "reaction": "DNA + H2O → nucleotides/oligonucleotides",
    },
}

# Quick lookup: EC -> nanozyme type
EC_TO_TYPE: Dict[str, str] = {
    v["primary_ec"]: k for k, v in NANOZYME_EC_MAPPINGS.items()
}
