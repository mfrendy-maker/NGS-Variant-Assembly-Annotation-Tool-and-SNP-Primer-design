# NGS Variant Assembly Annotation Tool

End-to-end NGS workflow for FASTQ processing, read assembly, reference genome alignment, variant calling, functional annotation, mutant identification, and primer design for detected mutations.

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

### CLI

```bash
ngs-variant-tool qc sample.fastq
ngs-variant-tool vcf-summary variants.vcf
ngs-variant-tool design-primer A T 24
```

### Streamlit dashboard

```bash
streamlit run app.py
```

## License

MIT
