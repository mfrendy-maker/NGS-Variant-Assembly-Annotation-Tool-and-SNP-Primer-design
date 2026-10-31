import argparse

from .core import compute_fastq_quality_summary, design_primers_for_variants, parse_vcf


def main() -> None:
    parser = argparse.ArgumentParser(description="NGS Variant Assembly Annotation Tool")
    subparsers = parser.add_subparsers(dest="command")

    qc_parser = subparsers.add_parser("qc", help="Summarize FASTQ quality")
    qc_parser.add_argument("fastq_path")

    vcf_parser = subparsers.add_parser("vcf-summary", help="Summarize VCF content")
    vcf_parser.add_argument("vcf_path")

    primer_parser = subparsers.add_parser("design-primer", help="Generate a primer design suggestion")
    primer_parser.add_argument("ref")
    primer_parser.add_argument("alt")
    primer_parser.add_argument("target_length", type=int, default=20)

    args = parser.parse_args()

    if args.command == "qc":
        with open(args.fastq_path, "r", encoding="utf-8") as handle:
            summary = compute_fastq_quality_summary(handle.read())
        print(summary)

    elif args.command == "vcf-summary":
        with open(args.vcf_path, "r", encoding="utf-8") as handle:
            variants = parse_vcf(handle.read())
        print(variants)

    elif args.command == "design-primer":
        variant = {"chrom": "chr1", "pos": 100, "ref": args.ref, "alt": args.alt, "type": "SNP"}
        print(design_primers_for_variants([variant], args.target_length)[0])

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
