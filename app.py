import streamlit as st

from ngs_variant_assembly_annotation_tool.core import (
    compute_fastq_quality_summary,
    parse_vcf,
    design_primers_for_variants,
)

st.title("NGS Variant Assembly Annotation Tool")

st.sidebar.header("Input")
raw_reads = st.sidebar.text_area(
    "Paste FASTQ-like sequence reads",
    value="@read1\nATGCGTACGTAGCTAGCTAG\n+\nIIIIIIIIIIIIIIIIIIII\n@read2\nATGCGTACGTAGCTAGC\n+\nIIIIIIIIIIIIIIIII\n",
    height=220,
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

    variants = parse_vcf(vcf_text)
    st.subheader("VCF variants")
    st.dataframe(variants)

    primer_candidates = design_primers_for_variants(variants)
    st.subheader("Primer design suggestions")
    st.json(primer_candidates)
