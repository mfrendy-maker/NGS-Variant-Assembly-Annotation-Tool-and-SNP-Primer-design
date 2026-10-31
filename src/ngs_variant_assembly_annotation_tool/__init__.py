"""NGS Variant Assembly Annotation Tool package."""

from .core import (
    annotate_variant,
    assemble_contigs,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    map_reads_to_reference,
    lookup_reference_gene,
    parse_vcf,
    summarize_reference_alignment,
)

__all__ = [
    "compute_fastq_quality_summary",
    "parse_vcf",
    "assemble_contigs",
    "map_reads_to_reference",
    "annotate_variant",
    "lookup_reference_gene",
    "design_primers_for_variants",
    "summarize_reference_alignment",
]
