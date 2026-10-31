<<<<<<< HEAD
# NGS Variant Assembly Annotation Tool and SNP Primer Design

A lightweight genome-analysis toolkit for reading FASTQ, FASTA, and VCF data, assembling reads, aligning sequences to a reference, summarizing variants, annotating mutations, and designing primers for SNP targets.

## Quick start

```bash
python -m pip install -e .
ngs-variant-tool qc examples/sample.fastq
ngs-variant-tool vcf-summary examples/variants.vcf
streamlit run app.py
```

## Features

- FASTQ quality summaries
- FASTA sequence parsing and validation
- contig assembly from overlapping reads
- reference-based alignment and mapping
- VCF parsing and variant extraction
- alignment-derived variant calling
- mutation annotation and effect summary
- local reference gene lookup
- SNP and indel primer design
- Streamlit browser dashboard

## Typical workflow

1. Run quality control on raw FASTQ reads.
2. Assemble reads into contigs.
3. Align reads or sequences to a reference.
4. Parse variants from VCF or alignment output.
5. Annotate candidate mutations.
6. Design PCR primers for validation.

## Command-line usage

### FASTQ QC

```bash
ngs-variant-tool qc sample.fastq
```

### VCF summary

```bash
ngs-variant-tool vcf-summary variants.vcf
```

### Alignment

```bash
ngs-variant-tool align ATGCGTACGA ATGCGTACGA
```

### Primer design

```bash
ngs-variant-tool design-primer A G 24
```

### Load a reference gene database

```bash
ngs-variant-tool load-db reference_database.json
```

## Project scope

This project is intended for gene-scale and small-to-medium NGS analysis workflows, with a focus on practical variant review, reference-based interpretation, and experimental validation.
=======
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

- read quality control
- genome assembly
- alignment to reference genomes
- variant calling to produce VCF files
- functional annotation of mutations
- reference genome database integration
- mutant identification and prioritization
- primer design for detected variants

## Features

- FASTQ quality summary for read length, count, and base composition
- reference-guided assembly summaries
- VCF parsing and variant extraction
- mutation classification into SNP or indel categories
- primer design suggestions for target variants
- lightweight CLI interface for common tasks
- browser-based dashboard for exploration

## Project status

This is an early-stage NGS toolkit scaffold intended for extension into a practical mutation-discovery workflow.

## Installation

```bash
python -m pip install -e .
```

## Quick start

### CLI

```bash
ngs-variant-tool qc sample.fastq
ngs-variant-tool vcf-summary variants.vcf
ngs-variant-tool design-primer A T G C 20
```

### Streamlit dashboard

```bash
streamlit run app.py
```

## Example workflow

1. Input FASTQ files
2. Run quality checks
3. Map reads to a reference genome
4. Call variants and produce VCF files
5. Annotate and filter mutations
6. Design primers for candidate mutant sites

## Scope

This project is intended for gene-scale and small-to-medium NGS analyses, with a focus on practical mutation discovery and downstream validation workflows.
>>>>>>> 68b6327 (Initialize NGS variant assembly annotation project scaffold)

## License

MIT
