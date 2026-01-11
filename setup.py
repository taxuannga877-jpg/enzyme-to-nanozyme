"""
Nanozyme Mining System Setup
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="nanozyme_mining",
    version="0.2.0",
    author="Nanozyme Design Team",
    author_email="",
    description="EC-based Nanozyme Database and Catalytic Motif Extraction System",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/taxuannga877-jpg/enzyme-to-nanozyme",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "requests>=2.25.0",
    ],
    extras_require={
        "full": [
            "torch>=2.0.0",
            "dgl>=2.0.0",
            "dgllife>=0.3.2",
            "rdkit>=2022.9.0",
            "fair-esm>=2.0.0",
        ],
    },
)
