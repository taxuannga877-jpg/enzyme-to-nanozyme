__version__ = "0.2.0"

from .utils.exceptions import ModelNotFoundError

try:
    from .database import NanozymeDatabase, UniProtFetcher, EnzymeEntry
except Exception:
    NanozymeDatabase = UniProtFetcher = EnzymeEntry = None

try:
    from .extraction import MotifExtractor, CatalyticMotif
except Exception:
    MotifExtractor = CatalyticMotif = None

__all__ = [
    "NanozymeDatabase",
    "UniProtFetcher",
    "EnzymeEntry",
    "MotifExtractor",
    "CatalyticMotif",
    "ModelNotFoundError",
]
