# Organelle_RNA_sequencing_analysis

Bioinformatics pipeline and sequencing data supporting the manuscript:  
**Hua, Z. (2025).** *Rapid and Cost-Effective Digital Quantification of RNA Editing and Maturation in Organelle Transcripts*. *The Plant Journal* (in review).

This repository contains all raw sequence data and custom scripts used in the analysis of nanopore long-read sequencing for RNA editing and intron retention in chloroplast transcripts.

---

## 📁 Contents

### Appendices
- `Appendix_S1`: `raw_FASTQs/` — Raw FASTQ files from nanopore barcoded amplicon sequencing  
- `Appendix_S2`: `pseudo_genome_ndhBD.fa` — Synthetic pseudo-genome reference for alignment 
- `Appendix_S3`: `ndhB_intron.fa` — Reference sequence of the Group II intron in `ndhB`

### Methods
- `Method_S1`: `filter_barcoded_fastq_seqs.pl` — Perl script for barcode filtering and strand correction  
- `Method_S2`: `run_minimap_alignments.sh` — Bash script for `minimap2` alignment and `samtools` processing  
- `Method_S3`: `count_groupII_inserts_and_extract.py` — Python script for detecting unspliced reads and exporting FASTA files  
- `Method_S4`: `local_align_inserts_vs_intron_blast.pl` — Perl script for `BLASTN`-based local alignment of unspliced reads  

### Guide
- `USAGE_GUIDE.md`: `Step-by-Step Usage Guide`  — Reproducing the full organelle RNA sequencing analysis pipeline
  
---

## 🧪 System & Dependencies

**Tested on:** macOS (ARM64)

### Required tools:
- **Perl** ≥ 5.40.1 with **BioPerl** ≥ 1.7.8  
- **Python** ≥ 3.11 with:
  - `biopython` ≥ 1.81  
  - `pandas` ≥ 2.2.2  
- **minimap2** ≥ v2.29  
- **samtools** ≥ v1.21  
- **NCBI BLAST+** ≥ v2.16.0  

---

## 📬 Contact
**Zhihua Hua, Ph.D.**  
Associate Professor  
Department of Environmental & Plant Biology  
Ohio University  
📧 hua@ohio.edu  
🔗 ORCID: [0000-0003-1177-1612](https://orcid.org/0000-0003-1177-1612)
