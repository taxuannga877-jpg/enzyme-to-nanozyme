"""
Database module for Nanozyme Mining System
==========================================

Stage 1: EC -> Nanozyme Function Type Mapping Database
"""

from .nanozyme_db import NanozymeDatabase, EnzymeEntry
from .uniprot_fetcher import UniProtFetcher

__all__ = [
    "NanozymeDatabase",
    "EnzymeEntry",
    "UniProtFetcher",
]
