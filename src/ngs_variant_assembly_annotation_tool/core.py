from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
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
        "sequence": "ATGCGTACGATCGATCGTAGCTAGCTAGCGATCGATCGATCGATCGATCG",
    },
    "TP53": {
        "gene_name": "TP53",
        "chromosome": "chr17",
        "reference_accession": "NM_000546",
        "description": "Tumor protein p53",
        "length": 1182,
        "sequence": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAA",
    },
    "EGFR": {
        "gene_name": "EGFR",
        "chromosome": "chr7",
        "reference_accession": "NM_005228",
        "description": "Epidermal growth factor receptor",
        "length": 1890,
        "sequence": "ATGCGACCCTCCGGGACGAGTTTGTGCGTGGAGACCGGCGGGCAGCGGGAG",
    },
    "BRAF": {
        "gene_name": "BRAF",
        "chromosome": "chr7",
        "reference_accession": "NM_004333",
        "description": "B-Raf proto-oncogene serine/threonine kinase",
        "length": 1890,
        "sequence": "ATGAGCACGCAAGGGGACCTGCTGAAATCTGCTAAAGTGATGAAATCCATC",
    },
    "PTEN": {
        "gene_name": "PTEN",
        "chromosome": "chr10",
        "reference_accession": "NM_000314",
        "description": "Phosphatase and tensin homolog",
        "length": 1200,
        "sequence": "ATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAG",
    },
}


def parse_fasta(fasta_text: str) -> Dict[str, str]:
    records: Dict[str, str] = {}
    current_header = ""
    sequence_parts: List[str] = []

    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                records[current_header] = "".join(sequence_parts).upper()
            current_header = line[1:].strip()
            sequence_parts = []
        elif current_header:
            sequence_parts.append(line.upper())

    if current_header:
        records[current_header] = "".join(sequence_parts).upper()

    return records


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


def load_reference_database(database: Any) -> Dict[str, Dict[str, Any]]:
    if isinstance(database, str):
        path = Path(database)
        if not path.exists():
            raise FileNotFoundError(f"Reference database not found: {database}")

        if path.suffix.lower() == ".json":
            with path.open("r", encoding="utf-8") as handle:
                database = json.load(handle)
        else:
            entries: Dict[str, Dict[str, Any]] = {}
            with path.open("r", encoding="utf-8") as handle:
                for line in handle:
                    row = [value.strip() for value in line.rstrip("\n").split("\t")]
                    if len(row) < 3 or not row[0] or row[0].startswith("#"):
                        continue
                    gene_name = row[0].upper()
                    entries[gene_name] = {
                        "gene_name": gene_name,
                        "chromosome": row[1],
                        "reference_accession": row[2],
                        "description": row[3] if len(row) > 3 else "",
                        "length": int(row[4]) if len(row) > 4 and row[4].isdigit() else 0,
                    }
            database = entries

    if not isinstance(database, dict):
        raise TypeError("Reference database must be a dictionary or a path to a JSON/TSV file.")

    normalized: Dict[str, Dict[str, Any]] = {}
    for gene_key, record in database.items():
        gene_name = str(record.get("gene_name", gene_key)).upper()
        normalized[gene_name] = {
            "gene_name": gene_name,
            "chromosome": record.get("chromosome", "unknown"),
            "reference_accession": record.get("reference_accession", "N/A"),
            "description": record.get("description", ""),
            "length": int(record.get("length") or len(str(record.get("sequence", "")))),
            "sequence": str(record.get("sequence", "")).upper(),
        }

    return normalized


def align_sequence_to_reference(read: str, reference: str) -> Dict[str, Any]:
    query = (read or "").upper()
    target = (reference or "").upper()

    if not query or not target:
        return {
            "score": 0,
            "matches": 0,
            "mismatches": 0,
            "gaps": 0,
            "identity": 0.0,
            "aligned_read": "",
            "aligned_reference": "",
            "status": "empty",
        }

    m = len(query)
    n = len(target)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    best_score = 0
    best_i = 0
    best_j = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if query[i - 1] == target[j - 1]:
                score = dp[i - 1][j - 1] + 2
            else:
                score = max(dp[i - 1][j - 1] - 1, dp[i - 1][j] - 2, dp[i][j - 1] - 2, 0)

            dp[i][j] = score
            if score > best_score:
                best_score = score
                best_i = i
                best_j = j

    aligned_read: List[str] = []
    aligned_reference: List[str] = []
    matches = 0
    mismatches = 0
    gaps = 0
    i = best_i
    j = best_j

    while i > 0 and j > 0 and dp[i][j] > 0:
        if query[i - 1] == target[j - 1]:
            aligned_read.append(query[i - 1])
            aligned_reference.append(target[j - 1])
            matches += 1
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j - 1] - 1:
            aligned_read.append(query[i - 1])
            aligned_reference.append(target[j - 1])
            mismatches += 1
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] - 2:
            aligned_read.append(query[i - 1])
            aligned_reference.append("-")
            gaps += 1
            i -= 1
        else:
            aligned_read.append("-")
            aligned_reference.append(target[j - 1])
            gaps += 1
            j -= 1

    aligned_read.reverse()
    aligned_reference.reverse()
    identity = (matches / max(1, max(len(query), len(target)))) * 100 if query and target else 0.0

    return {
        "score": best_score,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity": round(identity, 2),
        "aligned_read": "".join(aligned_read),
        "aligned_reference": "".join(aligned_reference),
        "status": "mapped" if best_score > 0 else "no_match",
    }


def map_reads_to_reference(reference: str, reads: List[str]) -> Dict[str, Any]:
    if not reference:
        return {"reference_length": 0, "mapped_reads": 0, "coverage": 0.0, "status": "not_mapped"}

    mapped_reads = 0
    alignments: List[Dict[str, Any]] = []
    for read in reads:
        alignment = align_sequence_to_reference(read, reference)
        alignments.append(alignment)
        if alignment["score"] > 0:
            mapped_reads += 1

    coverage = (mapped_reads / len(reads)) * 100 if reads else 0.0
    return {
        "reference_length": len(reference),
        "mapped_reads": mapped_reads,
        "coverage": round(coverage, 2),
        "status": "mapped" if mapped_reads else "not_mapped",
        "alignments": alignments,
    }


def call_variants_from_alignment(reference: str, reads: List[str], chrom: str = "chr1") -> List[Dict[str, Any]]:
    if not reference:
        return []

    variants: List[Dict[str, Any]] = []
    for read_index, read in enumerate(reads):
        if not read:
            continue

        if len(read) != len(reference):
            continue

        for position, (ref_base, read_base) in enumerate(zip(reference, read), start=1):
            if ref_base != read_base:
                variants.append(
                    {
                        "chrom": chrom,
                        "pos": position,
                        "ref": ref_base,
                        "alt": read_base,
                        "type": "SNP",
                        "read_index": read_index,
                    }
                )

    return variants


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
        "sequence": "",
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
