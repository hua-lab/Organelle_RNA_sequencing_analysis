# Organelle_RNA_sequencing_analysis

Bioinformatics pipeline and sequencing data supporting the manuscript:  
**Hua, Z. (2025).** *Rapid and Cost-Effective Digital Quantification of RNA Editing and Maturation in Organelle Transcripts*. *The Plant Journal* (in review).

This repository includes all raw sequence data and custom scripts used to analyze Oxford Nanopore long-read sequencing of RNA editing and intron retention in chloroplast transcripts.

The supplementary files listed below match those submitted with the manuscript to *The Plant Journal*.

---

## 📁 Contents

### Appendices  
- `Appendix S1.zip`: `raw_FASTQs/` — Raw FASTQ files from nanopore barcoded amplicon sequencing  
- `Appendix S2.txt`: `pseudo_genome_ndhBD.fa` — Synthetic pseudo-genome reference for alignment  
- `Appendix S3.txt`: `ndhB_intron.fa` — Reference sequence of the Group II intron in `ndhB`

### Methods  
- `Method S1.txt`: `filter_barcoded_fastq_seqs.pl` — Perl script for barcode filtering and strand correction  
- `Method S2.txt`: `run_minimap_alignments.sh` — Bash script for `minimap2` alignment and `samtools` processing  
- `Method S3.txt`: `count_groupII_inserts_and_extract.py` — Python script for detecting unspliced reads and exporting FASTA files  
- `Method S4.txt`: `local_align_inserts_vs_intron_blast.pl` — Perl script for `BLASTN`-based local alignment of unspliced reads  

### Guide  
- For detailed step-by-step instructions for reproducing the full organelle RNA sequencing analysis pipeline, see [USAGE_GUIDE.md](USAGE_GUIDE.md).


---

## 🧪 System & Dependencies

**Tested on:** macOS (ARM64)

### Required Tools  
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
