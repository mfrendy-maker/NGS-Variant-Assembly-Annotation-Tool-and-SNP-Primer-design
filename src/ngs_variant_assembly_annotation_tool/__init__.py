"""NGS Variant Assembly Annotation Tool package."""

from .core import (
    align_sequence_to_reference,
    annotate_variant,
    assemble_contigs,
    call_variants_from_alignment,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    load_reference_database,
    map_reads_to_reference,
    lookup_reference_gene,
    parse_fasta,
    parse_vcf,
    summarize_reference_alignment,
    summarize_variant_samples,
)

__all__ = [
    "compute_fastq_quality_summary",
    "parse_vcf",
    "parse_fasta",
    "assemble_contigs",
    "map_reads_to_reference",
    "align_sequence_to_reference",
    "call_variants_from_alignment",
    "load_reference_database",
    "annotate_variant",
    "lookup_reference_gene",
    "design_primers_for_variants",
    "summarize_reference_alignment",
    "summarize_variant_samples",
]
