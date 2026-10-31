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
