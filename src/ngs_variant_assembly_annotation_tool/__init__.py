"""NGS Variant Assembly Annotation Tool package."""

from .core import (
    align_sequence_to_reference,
    annotate_variant,
    assemble_contigs,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    load_reference_database,
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
    "align_sequence_to_reference",
    "load_reference_database",
    "annotate_variant",
    "lookup_reference_gene",
    "design_primers_for_variants",
    "summarize_reference_alignment",
]
