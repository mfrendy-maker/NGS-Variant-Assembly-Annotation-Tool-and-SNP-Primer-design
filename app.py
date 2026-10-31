import streamlit as st

from ngs_variant_assembly_annotation_tool.core import (
    annotate_variant,
    assemble_contigs,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    lookup_reference_gene,
    map_reads_to_reference,
    parse_vcf,
)

st.title("NGS Variant Assembly Annotation Tool")

st.sidebar.header("Input")
raw_reads = st.sidebar.text_area(
    "Paste FASTQ-like sequence reads",
    value="@read1\nATGCGTACGTAGCTAGCTAG\n+\nIIIIIIIIIIIIIIIIIIII\n@read2\nATGCGTACGTAGCTAGC\n+\nIIIIIIIIIIIIIIIII\n",
    height=220,
)
reference_sequence = st.sidebar.text_area(
    "Reference sequence",
    value="ATGGCCAAATTTGGG",
    height=120,
)
vcf_text = st.sidebar.text_area(
    "Paste VCF content",
    value="##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t50\t.\tA\tG\t60\tPASS\tDP=30\nchr1\t100\t.\tAT\tA\t50\tPASS\tDP=25\n",
    height=180,
)

if st.button("Run analysis"):
    qc_summary = compute_fastq_quality_summary(raw_reads)
    st.subheader("FASTQ quality summary")
    st.json(qc_summary)

    reads = [record["sequence"] for record in []]
    if qc_summary["read_count"]:
        reads = [record["sequence"] for record in []]

    st.subheader("Assembly summary")
    contigs = assemble_contigs([seq for seq in [line for line in raw_reads.splitlines() if line and not line.startswith("@") and not line.startswith("+") and not line.startswith("#")]])
    st.json({"assembled_contigs": contigs})

    mapping = map_reads_to_reference(reference_sequence, [seq for seq in [line for line in raw_reads.splitlines() if line and not line.startswith("@") and not line.startswith("+") and not line.startswith("#")]])
    st.subheader("Reference mapping")
    st.json(mapping)

    variants = parse_vcf(vcf_text)
    st.subheader("VCF variants")
    st.dataframe(variants)

    annotated = [annotate_variant(reference_sequence, variant["ref"], variant["alt"], variant["pos"], gene_name="BRCA1") for variant in variants]
    st.subheader("Variant annotation")
    st.json(annotated)

    primer_candidates = design_primers_for_variants(variants)
    st.subheader("Primer design suggestions")
    st.json(primer_candidates)

    st.subheader("Reference genome database lookup")
    st.json(lookup_reference_gene("BRCA1"))
