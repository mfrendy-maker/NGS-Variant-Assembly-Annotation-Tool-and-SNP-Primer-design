from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List


def _parse_fastq_records(fastq_text: str) -> List[Dict[str, str]]:
    lines = [line.strip() for line in fastq_text.splitlines() if line.strip()]
    records: List[Dict[str, str]] = []

    for index in range(0, len(lines) - 3, 4):
        if lines[index].startswith("@"):
            records.append(
                {
                    "header": lines[index],
                    "sequence": lines[index + 1],
                    "quality_header": lines[index + 2],
                    "quality": lines[index + 3],
                }
            )

    return records


def compute_fastq_quality_summary(fastq_text: str) -> Dict[str, Any]:
    records = _parse_fastq_records(fastq_text)

    if not records:
        return {
            "read_count": 0,
            "average_length": 0,
            "gc_content": 0.0,
            "base_counts": {},
            "reads": [],
        }

    sequences = [record["sequence"].upper() for record in records]
    total_length = sum(len(sequence) for sequence in sequences)
    total_bases = "".join(sequences)
    gc_count = total_bases.count("G") + total_bases.count("C")
    gc_content = (gc_count / len(total_bases)) * 100 if total_bases else 0.0

    base_counts = dict(Counter(total_bases))

    return {
        "read_count": len(records),
        "average_length": total_length / len(records),
        "gc_content": round(gc_content, 2),
        "base_counts": {base: base_counts.get(base, 0) for base in sorted(base_counts)},
        "reads": [
            {"header": record["header"], "length": len(record["sequence"])} for record in records
        ],
    }


def parse_vcf(vcf_text: str) -> List[Dict[str, Any]]:
    variants: List[Dict[str, Any]] = []

    for raw_line in vcf_text.splitlines():
        if not raw_line or raw_line.startswith("#"):
            continue

        parts = raw_line.strip().split("\t")
        if len(parts) < 8:
            continue

        chrom, pos, _identifier, ref, alt, qual, filter_value, info = parts[:8]

        if not ref or not alt:
            continue

        variant_type = "INDEL" if len(ref) != len(alt) else "SNP"

        variants.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "quality": float(qual) if qual not in {".", ""} else 0.0,
                "filter": filter_value,
                "info": info,
                "type": variant_type,
            }
        )

    return variants


def design_primers_for_variants(variants: List[Dict[str, Any]], target_length: int = 20) -> List[Dict[str, str]]:
    designed: List[Dict[str, str]] = []

    for variant in variants:
        variant_type = variant.get("type", "SNP")
        pos = variant.get("pos", 0)
        ref = str(variant.get("ref", "A")).upper()
        alt = str(variant.get("alt", "G")).upper()

        sequence_seed = ref + alt + "GATCGATCGATCGATCGATC"
        forward_primer = (sequence_seed[:target_length]).upper()
        reverse_primer = (sequence_seed[-target_length:][::-1]).upper()

        designed.append(
            {
                "chrom": str(variant.get("chrom", "chr1")),
                "position": pos,
                "variant_type": variant_type,
                "forward_primer": forward_primer,
                "reverse_primer": reverse_primer,
                "target_context": f"{ref}>{alt}",
            }
        )

    return designed


def summarize_reference_alignment(reference_name: str, alignment_score: float, mapped_reads: int) -> Dict[str, Any]:
    return {
        "reference_name": reference_name,
        "alignment_score": alignment_score,
        "mapped_reads": mapped_reads,
        "status": "mapped" if mapped_reads > 0 else "not_mapped",
    }
