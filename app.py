import streamlit as st

from ngs_variant_assembly_annotation_tool.core import (
    align_sequence_to_reference,
    annotate_variant,
    assemble_contigs,
    call_variants_from_alignment,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    load_reference_database,
    lookup_reference_gene,
    map_reads_to_reference,
    parse_fasta,
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

uploaded_fasta = st.sidebar.file_uploader("Load reference FASTA", type=["fasta", "fa", "txt"])
uploaded_vcf = st.sidebar.file_uploader("Load VCF file", type=["vcf", "txt"])

if uploaded_fasta is not None:
    fasta_records = parse_fasta(uploaded_fasta.read().decode("utf-8"))
    if fasta_records:
        reference_sequence = next(iter(fasta_records.values()))
if uploaded_vcf is not None:
    vcf_text = uploaded_vcf.read().decode("utf-8")

if st.button("Run analysis"):
    qc_summary = compute_fastq_quality_summary(raw_reads)
    st.subheader("FASTQ quality summary")
    st.json(qc_summary)

    read_sequences = [
        line.strip()
        for line in raw_reads.splitlines()
        if line.strip() and not line.startswith("@") and not line.startswith("+")
    ]

    st.subheader("Assembly summary")
    contigs = assemble_contigs(read_sequences)
    st.json({"assembled_contigs": contigs})

    mapping = map_reads_to_reference(reference_sequence, read_sequences)
    st.subheader("Reference mapping")
    st.json(mapping)

    if read_sequences:
        variant_candidates = call_variants_from_alignment(reference_sequence, read_sequences)
        st.subheader("Alignment-derived variants")
        st.json(variant_candidates[:10])

    variants = parse_vcf(vcf_text)
    st.subheader("VCF variants")
    st.dataframe(variants)

    annotated = [
        annotate_variant(reference_sequence, variant["ref"], variant["alt"], variant["pos"], gene_name="BRCA1")
        for variant in variants
    ]
    st.subheader("Variant annotation")
    st.json(annotated)

    primer_candidates = design_primers_for_variants(variants)
    st.subheader("Primer design suggestions")
    st.json(primer_candidates)

    st.subheader("Reference genome database lookup")
    st.json(lookup_reference_gene("BRCA1"))

    st.subheader("Sequence alignment")
    st.json(align_sequence_to_reference(reference_sequence, reference_sequence))

    st.subheader("Database loader")
    st.json(load_reference_database({"BRCA1": lookup_reference_gene("BRCA1")}))
