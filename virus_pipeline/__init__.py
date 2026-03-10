# virus_pipeline/__init__.py

"""
Virus Variant Calling Pipeline

This package provides modular scripts for processing dengue virus sequencing data.
It includes steps for:
- Sample sheet generation
- Read mapping
- SAM/BAM conversion
- Consensus sequence generation
- SnpEff database creation
- Variant calling and annotation
- Result summarization

Each script is designed to be imported as a module or run independently.

Author: Rajindra Napit
License: MIT
"""

__version__ = "1.0.0"

from .create_samplesheet import main as create_samplesheet
from .map_reads import main as map_reads
from .samtobamdenv import main as samtobamdenv
from .create_snpeff_database import main as create_snpeff_database
from .variant_calling_consensus import main as variant_calling_consensus
from .summarize_result import main as summarize_result
from .summarize_snpEff import main as summarize_snpEff