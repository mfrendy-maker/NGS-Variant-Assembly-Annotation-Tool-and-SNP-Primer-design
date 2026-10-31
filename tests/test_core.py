import json
import sys

from ngs_variant_assembly_annotation_tool.core import (
    compute_fastq_quality_summary,
    parse_vcf,
    design_primers_for_variants,
)


def test_compute_fastq_quality_summary():
    reads = "@read1\nATGC\n+\nIIII\n@read2\nAAAA\n+\nIIII\n"
    result = compute_fastq_quality_summary(reads)

    assert result["read_count"] == 2
    assert result["average_length"] == 4
    assert result["gc_content"] > 0


def test_parse_vcf():
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=30\nchr1\t150\t.\tAT\tA\t50\tPASS\tDP=25\n"
    variants = parse_vcf(vcf)

    assert len(variants) == 2
    assert variants[0]["chrom"] == "chr1"
    assert variants[0]["ref"] == "A"
    assert variants[0]["alt"] == "G"


def test_design_primers_for_variants():
    variants = [
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "type": "SNP"},
        {"chrom": "chr1", "pos": 150, "ref": "AT", "alt": "A", "type": "INDEL"},
    ]
    primers = design_primers_for_variants(variants)

    assert len(primers) == 2
    assert "forward_primer" in primers[0]
    assert "reverse_primer" in primers[0]
    assert primers[0]["variant_type"] in {"SNP", "INDEL"}


def test_assemble_contigs_simple_overlap():
    from ngs_variant_assembly_annotation_tool.core import assemble_contigs

    contigs = assemble_contigs(["ATGCG", "GCGTA", "GTAGC"])

    assert len(contigs) == 1
    assert contigs[0].startswith("ATGCG")
    assert len(contigs[0]) >= 10


def test_map_reads_to_reference_and_annotate_variant():
    from ngs_variant_assembly_annotation_tool.core import (
        map_reads_to_reference,
        annotate_variant,
    )

    reference = "ATGGCCAAATTTGGG"
    reads = ["ATGGCCAAATTTGGG", "ATGGCCAAATTTGGG"]
    mapping = map_reads_to_reference(reference, reads)
    annotation = annotate_variant(reference, "A", "G", 5)

    assert mapping["mapped_reads"] == 2
    assert mapping["coverage"] > 0
    assert annotation["effect"] in {"missense", "synonymous", "nonsynonymous", "in-frame"}


def test_lookup_reference_gene_and_snp_primer_length():
    from ngs_variant_assembly_annotation_tool.core import lookup_reference_gene, design_primers_for_variants

    gene = lookup_reference_gene("BRCA1")
    primers = design_primers_for_variants([
        {"chrom": "chr17", "pos": 430, "ref": "C", "alt": "T", "type": "SNP"}
    ], target_length=24)

    assert gene["gene_name"] == "BRCA1"
    assert len(primers[0]["forward_primer"]) == 24
    assert len(primers[0]["reverse_primer"]) == 24


def test_load_reference_database_and_align_read():
    from ngs_variant_assembly_annotation_tool.core import (
        align_sequence_to_reference,
        load_reference_database,
    )

    catalog = load_reference_database({
        "BRCA1": {"gene_name": "BRCA1", "chromosome": "chr17", "sequence": "ATGCGTACGA"}
    })
    alignment = align_sequence_to_reference("ATGCGTACGA", "ATGCGTACGA")

    assert catalog["BRCA1"]["gene_name"] == "BRCA1"
    assert alignment["score"] > 0
    assert alignment["identity"] >= 0.95


def test_parse_fasta_and_call_variants_from_alignment():
    from ngs_variant_assembly_annotation_tool.core import (
        call_variants_from_alignment,
        parse_fasta,
    )

    fasta = ">chr1\nATGCGTACGA\n"
    records = parse_fasta(fasta)
    variants = call_variants_from_alignment("ATGCGTACGA", ["ATGCGTACGA", "ATGCGTACCA"]) 

    assert records["chr1"].startswith("ATGCGT")
    assert len(variants) >= 1
    assert variants[0]["ref"] in {"A", "C", "G", "T"}


def test_summarize_variant_samples():
    from ngs_variant_assembly_annotation_tool.core import summarize_variant_samples

    sample_variants = {
        "sample_1": [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "type": "SNP"},
            {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "type": "SNP"},
        ],
        "sample_2": [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "type": "SNP"},
        ],
    }

    summary = summarize_variant_samples(sample_variants)

    assert summary["sample_count"] == 2
    assert summary["total_variants"] == 3
    assert "chr1:10" in summary["shared_variants"]
    assert summary["variant_types"]["SNP"] == 3


def test_cli_sample_summary(tmp_path, monkeypatch, capsys):
    from ngs_variant_assembly_annotation_tool import cli

    sample_variants = {
        "sample_1": [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "type": "SNP"},
            {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "type": "SNP"},
        ],
        "sample_2": [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "type": "SNP"},
        ],
    }

    sample_file = tmp_path / "sample_variants.json"
    sample_file.write_text(json.dumps(sample_variants), encoding="utf-8")

    monkeypatch.setattr(sys, "argv", ["ngs-variant-tool", "sample-summary", str(sample_file)])
    cli.main()

    captured = capsys.readouterr().out
    assert "sample_count" in captured
    assert "shared_variants" in captured
    assert "chr1:10" in captured
