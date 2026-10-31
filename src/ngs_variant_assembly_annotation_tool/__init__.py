"""NGS Variant Assembly Annotation Tool package."""

from .core import (
    compute_fastq_quality_summary,
    parse_vcf,
    design_primers_for_variants,
)

__all__ = [
    "compute_fastq_quality_summary",
    "parse_vcf",
    "design_primers_for_variants",
]
