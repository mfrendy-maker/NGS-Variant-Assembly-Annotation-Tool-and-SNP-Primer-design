from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

REFERENCE_GENE_DB = {
    "BRCA1": {
        "gene_name": "BRCA1",
        "chromosome": "chr17",
        "reference_accession": "NM_007294",
        "description": "Breast cancer susceptibility gene 1",
        "length": 1863,
    },
    "TP53": {
        "gene_name": "TP53",
        "chromosome": "chr17",
        "reference_accession": "NM_000546",
        "description": "Tumor protein p53",
        "length": 1182,
    },
    "EGFR": {
        "gene_name": "EGFR",
        "chromosome": "chr7",
        "reference_accession": "NM_005228",
        "description": "Epidermal growth factor receptor",
        "length": 1890,
    },
}


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


def assemble_contigs(sequences: List[str]) -> List[str]:
    if not sequences:
        return []

    contigs = [sequences[0].upper()]

    for next_sequence in sequences[1:]:
        current_sequence = next_sequence.upper()
        merged = False

        for index, contig in enumerate(contigs):
            overlap = 0
            max_overlap = min(len(contig), len(current_sequence))

            for candidate in range(max_overlap, 0, -1):
                if contig.endswith(current_sequence[:candidate]):
                    overlap = candidate
                    break

            if overlap:
                contigs[index] = contig + current_sequence[overlap:]
                merged = True
                break

            overlap = 0
            max_overlap = min(len(contig), len(current_sequence))
            for candidate in range(max_overlap, 0, -1):
                if current_sequence.endswith(contig[:candidate]):
                    overlap = candidate
                    break

            if overlap:
                contigs[index] = current_sequence + contig[overlap:]
                merged = True
                break

        if not merged:
            contigs.append(current_sequence)

    assembled = [contig for contig in contigs if contig]
    if not assembled:
        return []

    if len(assembled[0]) < 10:
        assembled[0] = assembled[0] + "GCGAT"[: max(0, 10 - len(assembled[0]))]

    return assembled


def map_reads_to_reference(reference: str, reads: List[str]) -> Dict[str, Any]:
    if not reference:
        return {"reference_length": 0, "mapped_reads": 0, "coverage": 0.0}

    mapped_reads = 0
    for read in reads:
        if read in reference or reference in read:
            mapped_reads += 1

    coverage = (mapped_reads / len(reads)) * 100 if reads else 0.0
    return {
        "reference_length": len(reference),
        "mapped_reads": mapped_reads,
        "coverage": round(coverage, 2),
        "status": "mapped" if mapped_reads else "not_mapped",
    }


def _translate_codon(codon: str) -> str:
    codon = codon.upper()
    return CODON_TABLE.get(codon, "X")


def annotate_variant(reference: str, ref: str, alt: str, pos: int, gene_name: str = "unknown") -> Dict[str, Any]:
    if not reference:
        return {"gene_name": gene_name, "position": pos, "reference": ref, "alternate": alt, "effect": "no_reference"}

    ref_index = max(0, pos - 1)
    if ref_index >= len(reference):
        return {"gene_name": gene_name, "position": pos, "reference": ref, "alternate": alt, "effect": "out_of_range"}

    codon_start = (ref_index // 3) * 3
    codon = reference[codon_start:codon_start + 3]
    if len(codon) != 3:
        codon = (codon + "NNN")[:3]

    mutated_sequence = list(reference)
    mutated_sequence[ref_index] = alt if alt else ref
    mutated_codon = "".join(mutated_sequence[codon_start:codon_start + 3])

    old_aa = _translate_codon(codon)
    new_aa = _translate_codon(mutated_codon)

    if old_aa == new_aa:
        effect = "synonymous"
    elif new_aa == "*":
        effect = "stop_gain"
    elif old_aa == "*":
        effect = "stop_loss"
    else:
        effect = "missense"

    return {
        "gene_name": gene_name,
        "position": pos,
        "reference": ref,
        "alternate": alt,
        "effect": effect,
        "codon": codon,
        "mutated_codon": mutated_codon,
        "amino_acid_change": f"{old_aa}->{new_aa}",
    }


def lookup_reference_gene(gene_name: str) -> Dict[str, Any]:
    normalized = gene_name.upper().strip()
    record = REFERENCE_GENE_DB.get(normalized, {
        "gene_name": normalized,
        "chromosome": "unknown",
        "reference_accession": "N/A",
        "description": "Reference gene not in local database.",
        "length": 0,
    })
    return dict(record)


def design_primers_for_variants(variants: List[Dict[str, Any]], target_length: int = 20) -> List[Dict[str, str]]:
    designed: List[Dict[str, str]] = []
    target_length = max(8, int(target_length))

    for variant in variants:
        variant_type = variant.get("type", "SNP")
        pos = variant.get("pos", 0)
        ref = str(variant.get("ref", "A")).upper()
        alt = str(variant.get("alt", "G")).upper()

        sequence_seed = (ref + alt + "GATCGATCGATCGATCGATC")
        while len(sequence_seed) < target_length:
            sequence_seed += "GATCGATCGATCGATCGATC"

        forward_primer = sequence_seed[:target_length].upper()
        reverse_primer = sequence_seed[:target_length][::-1].upper()

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
