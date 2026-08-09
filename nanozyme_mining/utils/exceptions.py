"""
Custom exceptions for Nanozyme Mining System
=============================================

This module defines custom exception classes used throughout the system
for better error handling and clearer error messages.
"""


class NanozymeMiningError(Exception):
    """Base exception for all Nanozyme Mining errors."""
    pass


class ModelNotFoundError(NanozymeMiningError):
    """
    Raised when a required model is not found or unavailable.

    This exception is raised during initialization when:
    - DiffLinker model checkpoints are not found
    - MetalloGen model is unavailable
    """

    def __init__(self, model_name: str, path: str = "", details: str = ""):
        """
        Initialize ModelNotFoundError.

        Args:
            model_name: Name of the model that was not found
            path: Path where the model was expected
            details: Additional details about the error
        """
        self.model_name = model_name
        self.path = path
        self.details = details

        message = f"Model '{model_name}' not found"
        if path:
            message += f" at path: {path}"
        if details:
            message += f". {details}"

        super().__init__(message)


class LinkerGenerationError(NanozymeMiningError):
    """Raised when linker generation fails."""
    pass


class MetalNotFoundError(NanozymeMiningError):
    """Raised when no metal atom is found in a structure."""
    pass


class AssemblyError(NanozymeMiningError):
    """Raised when nanozyme assembly fails."""
    pass


class ParseError(NanozymeMiningError):
    """Raised when parsing a file fails."""
    pass
