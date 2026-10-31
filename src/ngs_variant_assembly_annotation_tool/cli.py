import argparse
import json

from .core import (
    align_sequence_to_reference,
    compute_fastq_quality_summary,
    design_primers_for_variants,
    load_reference_database,
    parse_fasta,
    parse_vcf,
    summarize_variant_samples,
)


def _summarize_fasta_file(fasta_path: str) -> dict:
    with open(fasta_path, "r", encoding="utf-8") as handle:
        records = parse_fasta(handle.read())

    return {
        "record_count": len(records),
        "total_bases": sum(len(sequence) for sequence in records.values()),
        "records": {header: len(sequence) for header, sequence in records.items()},
        "sequence_names": list(records.keys()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="NGS Variant Assembly Annotation Tool")
    subparsers = parser.add_subparsers(dest="command")

    qc_parser = subparsers.add_parser("qc", help="Summarize FASTQ quality")
    qc_parser.add_argument("fastq_path")

    fasta_parser = subparsers.add_parser("fasta-summary", help="Summarize FASTA content")
    fasta_parser.add_argument("fasta_path")

    vcf_parser = subparsers.add_parser("vcf-summary", help="Summarize VCF content")
    vcf_parser.add_argument("vcf_path")

    sample_parser = subparsers.add_parser("sample-summary", help="Summarize variants across multiple samples")
    sample_parser.add_argument("sample_json_path")

    primer_parser = subparsers.add_parser("design-primer", help="Generate a primer design suggestion")
    primer_parser.add_argument("ref")
    primer_parser.add_argument("alt")
    primer_parser.add_argument("target_length", type=int, default=20)

    align_parser = subparsers.add_parser("align", help="Align a read to a reference sequence")
    align_parser.add_argument("read")
    align_parser.add_argument("reference")

    db_parser = subparsers.add_parser("load-db", help="Load a reference gene database")
    db_parser.add_argument("database_path")

    args = parser.parse_args()

    if args.command == "qc":
        with open(args.fastq_path, "r", encoding="utf-8") as handle:
            summary = compute_fastq_quality_summary(handle.read())
        print(summary)

    elif args.command == "fasta-summary":
        print(_summarize_fasta_file(args.fasta_path))

    elif args.command == "vcf-summary":
        with open(args.vcf_path, "r", encoding="utf-8") as handle:
            variants = parse_vcf(handle.read())
        print(variants)

    elif args.command == "sample-summary":
        with open(args.sample_json_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        print(summarize_variant_samples(payload))

    elif args.command == "design-primer":
        variant = {"chrom": "chr1", "pos": 100, "ref": args.ref, "alt": args.alt, "type": "SNP"}
        print(design_primers_for_variants([variant], args.target_length)[0])

    elif args.command == "align":
        print(align_sequence_to_reference(args.read, args.reference))

    elif args.command == "load-db":
        print(load_reference_database(args.database_path))

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
