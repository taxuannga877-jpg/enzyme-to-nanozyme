"""
Constants for Nanozyme Mining System
=====================================

Defines EC number mappings to nanozyme function types and other constants.
"""

from typing import Dict, List, Tuple
from enum import Enum

# ============================================================================
# Download Configuration Constants
# ============================================================================

# Default download limits (increased from original small values)
DEFAULT_UNIPROT_QUERY_SIZE = 500  # Default number of entries to query from UniProt
DEFAULT_MAX_ENTRIES = 500  # Default max entries for fetch_and_populate
DEFAULT_TOPK_ENZYMES = 50  # Default top-k enzymes for motif extraction

# Use -1 to fetch all available results
FETCH_ALL_RESULTS = -1

# UniProt API pagination settings
UNIPROT_PAGE_SIZE = 500  # Maximum page size supported by UniProt API


class NanozymeType(Enum):
    """Nanozyme function types based on catalytic activity."""
    POD = "Peroxidase"           # 过氧化物酶
    SOD = "Superoxide Dismutase" # 超氧化物歧化酶
    OXD = "Oxidase"              # 氧化酶
    CAT = "Catalase"             # 过氧化氢酶
    GSH = "Glutathione Peroxidase"  # 谷胱甘肽过氧化物酶
    LAC = "Laccase"              # 漆酶
    HRP = "Horseradish Peroxidase"  # 辣根过氧化物酶
    GOX = "Glucose Oxidase"      # 葡萄糖氧化酶
    PHOS = "Phosphatase"         # 磷酸酶
    DNASE = "DNase"              # 脱氧核酸酶
    UNKNOWN = "Unknown"


class ActiveSiteType(Enum):
    """Active site classification types."""
    NONE = 0           # 非活性位点
    BINDING = 1        # 结合位点
    CATALYTIC = 2      # 催化位点（活性中心）
    OTHER = 3          # 其他位点


# EC Number to Nanozyme Type Mapping
# Based on enzyme classification system
EC_TO_NANOZYME_TYPE: Dict[str, NanozymeType] = {
    # Peroxidase (POD) - EC 1.11.1.x
    "1.11.1.7": NanozymeType.POD,   # Peroxidase
    "1.11.1.11": NanozymeType.POD,  # L-ascorbate peroxidase
    "1.11.1.21": NanozymeType.POD,  # Catalase-peroxidase (KatG)

    # Catalase (CAT) - EC 1.11.1.6
    "1.11.1.6": NanozymeType.CAT,   # Catalase

    # Superoxide Dismutase (SOD) - EC 1.15.1.1
    "1.15.1.1": NanozymeType.SOD,   # Superoxide dismutase

    # Glutathione Peroxidase (GSH) - EC 1.11.1.9
    "1.11.1.9": NanozymeType.GSH,   # Glutathione peroxidase
    "1.11.1.12": NanozymeType.GSH,  # Phospholipid-hydroperoxide GPx

    # Oxidase (OXD) - EC 1.1.3.x, 1.4.3.x
    "1.1.3.4": NanozymeType.GOX,    # Glucose oxidase
    "1.10.3.2": NanozymeType.LAC,   # Laccase
    "1.4.3.4": NanozymeType.OXD,    # Monoamine oxidase
    "1.3.3.4": NanozymeType.OXD,    # Protoporphyrinogen oxidase

    # Phosphatase (PHOS) - EC 3.1.3.1
    "3.1.3.1": NanozymeType.PHOS,   # Phosphatase

    # DNase (DNASE) - EC 3.1.21.1
    "3.1.21.1": NanozymeType.DNASE, # DNase
}
