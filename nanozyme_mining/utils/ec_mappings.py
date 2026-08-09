"""
EC Number Mappings for Nanozyme Types
======================================

Based on literature and IUBMB nomenclature.

This module is the **single source of truth** for nanozyme EC mappings.
Downstream consumers (enzyme_viewer/app.py, scripts/*) should import from here
rather than maintaining parallel tables.
"""

from collections import defaultdict
from typing import Dict, List

from .constants import NanozymeType


# ============================================================
# EC_PATTERNS: NanozymeType -> List of EC numbers
#
# Used by database module for batch fetching. Single key per NanozymeType
# (Python dict semantics would silently drop duplicate keys; see PR0-1 fix).
# OXD entries from line 19 and line 22 of the legacy version are merged here.
# ============================================================
EC_PATTERNS: Dict[NanozymeType, List[str]] = {
    NanozymeType.POD: ["1.11.1.7", "1.11.1.11", "1.11.1.21"],            # Peroxidases (incl. KatG bifunctional)
    NanozymeType.CAT: ["1.11.1.6", "1.11.1.21"],                          # Catalases (incl. KatG bifunctional)
    NanozymeType.SOD: ["1.15.1.1"],                                       # Superoxide dismutase
    NanozymeType.GSH: ["1.11.1.9", "1.11.1.12"],                          # Glutathione peroxidases
    NanozymeType.OXD: ["1.4.3.4", "1.3.3.4", "1.1.3.5", "1.1.3.9"],       # Oxidases (merged superset; was duplicated)
    NanozymeType.LAC: ["1.10.3.2", "1.10.3.3", "1.10.3.4"],               # Laccase family
    NanozymeType.GOX: ["1.1.3.4"],                                        # Glucose oxidase
    NanozymeType.HRP: ["1.11.1.7"],                                       # Horseradish peroxidase (subtype of POD)
    NanozymeType.PHOS: ["3.1.3.1"],                                       # Phosphatase
    NanozymeType.DNASE: ["3.1.21.1"],                                     # DNase
    NanozymeType.URE: ["3.5.1.5"],                                        # Urease
}


# ============================================================
# NANOZYME_EC_MAPPINGS: short-code -> { name, primary_ec, ... }
#
# Display-oriented metadata for the web UI. Distinct short-codes (POD, KatG,
# GPx, ...) can legitimately share a primary_ec (e.g. POD + HRP both = 1.11.1.7);
# downstream lookups MUST therefore not assume a single-valued EC -> code map
# (see EC_TO_TYPE below).
# ============================================================
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


# ============================================================
# EC_TO_TYPE: EC number -> List[short-code]
#
# CRITICAL: This is now Dict[str, List[str]] (was Dict[str, str] which silently
# dropped duplicates). If someone adds {"HRP": {"primary_ec": "1.11.1.7"}} to
# NANOZYME_EC_MAPPINGS, the previous List[str] form would lose POD; this form
# preserves both:  EC_TO_TYPE["1.11.1.7"] == ["POD", "HRP"]
# ============================================================
def _build_ec_to_type() -> Dict[str, List[str]]:
    result: Dict[str, List[str]] = defaultdict(list)
    for short_code, meta in NANOZYME_EC_MAPPINGS.items():
        ec = meta.get("primary_ec")
        if ec:
            result[ec].append(short_code)
    return dict(result)


EC_TO_TYPE: Dict[str, List[str]] = _build_ec_to_type()


def lookup_types(ec_number: str) -> List[str]:
    """Return all nanozyme short-codes for an EC number (may be empty)."""
    if not ec_number:
        return []
    return EC_TO_TYPE.get(ec_number, [])


# ============================================================
# EC_ACTIVITY_LABELS: prefix/exact EC -> short label
#
# Migrated from enzyme_viewer/app.py to make this module the single source.
# Mixed prefix + exact matching is preserved for backward compatibility, but
# all entries explicitly named in NANOZYME_EC_MAPPINGS take precedence via the
# exact match in get_ec_activity_label() below.
#
# Prefix entries cover EC families that aren't enumerated above (e.g. esterase
# 3.1.1.x, protease 3.4.x.x).
# ============================================================
EC_ACTIVITY_LABELS: Dict[str, str] = {
    # Prefix entries (legacy from app.py — used as fallback)
    "1.11.1": "POD",   # Peroxidase family
    "1.11.2": "POD",   # Peroxidase family
    "1.1.3":  "OXD",   # Oxidase family
    "1.4.3":  "OXD",   # Amino-acid oxidase family
    "1.3.3":  "OXD",   # CH-CH oxidase family
    "1.10.3": "LAC",   # Laccase family
    "1.13":   "OXD",   # Oxygenase family
    "1.14":   "OXD",   # Monooxygenase family
    "3.1.1":  "EST",   # Esterase
    "3.1.3":  "PHO",   # Phosphatase
    "3.4":    "PRO",   # Protease
    "2.7.1":  "KIN",   # Kinase
    # Exact entries (override the prefixes above when EC matches exactly)
    "1.15.1.1":  "SOD",  # Superoxide dismutase
    "1.11.1.6":  "CAT",  # Catalase
    "1.11.1.21": "CAT/POD",  # KatG bifunctional — CRITICAL: prefix "1.11.1" would mis-label this as POD only
    "3.5.1.5":   "URE",  # Urease
    "3.1.21.1":  "DNase",
    "1.1.3.4":   "GOX",  # Glucose oxidase
}


def get_ec_activity_label(ec_number: str) -> str:
    """
    Return the short activity label for an EC number.

    Resolution order:
      1. Exact match in EC_ACTIVITY_LABELS (most specific)
      2. Exact match via NANOZYME_EC_MAPPINGS / EC_TO_TYPE
      3. Prefix match in EC_ACTIVITY_LABELS (longest prefix wins)
      4. Empty string
    """
    if not ec_number:
        return ""
    # 1. Exact match in label table
    if ec_number in EC_ACTIVITY_LABELS:
        return EC_ACTIVITY_LABELS[ec_number]
    # 2. Exact match via the canonical type table (e.g. 1.11.1.7 -> "POD")
    types = lookup_types(ec_number)
    if types:
        return types[0]
    # 3. Longest prefix match (so "1.11.1.x" hits "1.11.1" before "1.11")
    best_prefix = ""
    best_label = ""
    for prefix, label in EC_ACTIVITY_LABELS.items():
        if ec_number.startswith(prefix) and len(prefix) > len(best_prefix):
            best_prefix = prefix
            best_label = label
    return best_label
