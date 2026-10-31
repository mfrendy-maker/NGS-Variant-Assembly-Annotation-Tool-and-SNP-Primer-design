# NGS Variant Assembly Annotation Tool

End-to-end NGS workflow for FASTQ processing, read assembly, reference genome alignment, variant calling, functional annotation, mutant identification, and primer design for detected mutations.

## Quick start

This repository includes example data under the examples folder so you can run the tool immediately.

```bash
python -m pip install -e .
ngs-variant-tool qc examples/sample.fastq
ngs-variant-tool vcf-summary examples/variants.vcf
ngs-variant-tool design-primer A G 24
streamlit run app.py
```

- Example FASTQ file: `examples/sample.fastq`
- Example VCF file: `examples/variants.vcf`

## Overview

This project is designed for next-generation sequencing (NGS) analysis workflows, from raw FASTQ data to biological interpretation and experimental validation. It integrates:

- FASTQ quality control and summary statistics
- simple read assembly and contig reconstruction
- alignment and mapping to a reference genome
- variant calling and VCF parsing
- reference genome database lookups for gene context
- annotation of mutation effects
- mutant prioritization for downstream studies
- primer design for mutation validation and PCR screening

## Features

- FASTQ quality summary for read length, count, and base composition
- basic contig assembly from overlapping reads
- mapping summary against a reference sequence
- VCF parsing and variant extraction
- mutation annotation and effect prediction
- local reference gene lookup for BRCA1, TP53, and EGFR
- primer design suggestions for SNP and indel targets
- lightweight CLI interface for common tasks
- browser-based dashboard for exploration

## Example workflow

1. Input FASTQ files
2. Run quality checks
3. Assemble overlapping reads into contigs
4. Map reads against a reference genome
5. Call variants and produce VCF files
6. Annotate mutation effects against reference context
7. Design primers for mutation validation

## Scope

This project is intended for gene-scale and small-to-medium NGS analyses, with a focus on practical mutation discovery, reference-based interpretation, and downstream validation workflows.

## Installation

```bash
python -m pip install -e .
```

## Quick start

### Install the project

```bash
python -m pip install -e .
```

### CLI commands

After installation, you can run the command-line tools directly from the terminal.

#### 1. FASTQ quality check

```bash
ngs-variant-tool qc sample.fastq
```

This prints a summary with:
- read count
- average read length
- GC content
- base composition

#### 2. VCF summary

```bash
ngs-variant-tool vcf-summary variants.vcf
```

This parses a VCF file and prints variant records with chromosome, position, REF/ALT, quality, and type.

#### 3. Primer design for a mutation

```bash
ngs-variant-tool design-primer A T 24
```

This generates a primer pair for a variant context using the reference and alternative allele and the target primer length.

#### 4. Align a read to a reference sequence

```bash
ngs-variant-tool align ACGTACGTTT ACGTACGTTT
```

This performs a simple sequence alignment and reports:
- alignment score
- matches
- mismatches
- gaps
- percent identity

#### 5. Load a gene reference database

```bash
ngs-variant-tool load-db reference_database.json
```

This loads a JSON or TSV-style gene database and returns gene entries with chromosome, accession, and sequence metadata.

### Streamlit dashboard

```bash
streamlit run app.py
```

The dashboard accepts:
- FASTQ-like reads
- FASTA reference inputs
- VCF content
- alignment and summary visual outputs

### Example workflow

```bash
ngs-variant-tool qc sample.fastq
ngs-variant-tool align ATGCGTACGA ATGCGTACGA
ngs-variant-tool vcf-summary variants.vcf
ngs-variant-tool design-primer A G 24
```

## License

MIT
