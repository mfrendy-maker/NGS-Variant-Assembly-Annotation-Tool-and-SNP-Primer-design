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

### Primer design method and parameter defaults

This project uses a lightweight custom algorithm for SNP primer generation rather than an external primer-design engine such as Primer3. The implementation is in `src/ngs_variant_assembly_annotation_tool/core.py` and is exposed through the CLI in `src/ngs_variant_assembly_annotation_tool/cli.py`.

The primer-generation function is:

```python
design_primers_for_variants(variants, target_length=20)
```

Parameter behavior:

- `variants`: list of variant dictionaries, each typically containing `chrom`, `pos`, `ref`, `alt`, and `type`.
- `target_length`: target primer length in nucleotides.
- default `target_length`: `20`
- minimum enforced primer length: `8`

The function builds a seed sequence from the variant context using the reference and alternate bases and then pads with a repeated `GATC` motif until it reaches the target length:

```python
sequence_seed = (ref + alt + "GATCGATCGATCGATCGATC")
while len(sequence_seed) < target_length:
    sequence_seed += "GATCGATCGATCGATCGATC"
```

This produces both a forward and reverse primer:

```python
forward_primer = sequence_seed[:target_length].upper()
reverse_primer = sequence_seed[:target_length][::-1].upper()
```

Important note: this is a demonstration-level primer generator meant for teaching, portfolio work, and prototype workflows. It does not include thermodynamic optimization, Tm checking, hairpin/dimer screening, or PCR-specific constraints used in production primer design workflows.

### Load a reference gene database

```bash
ngs-variant-tool load-db reference_database.json
```

## Project scope

This project is intended for gene-scale and small-to-medium NGS analysis workflows, with a focus on practical variant review, reference-based interpretation, and experimental validation.

## License

MIT
