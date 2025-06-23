# Step-by-Step Usage Guide

This guide walks you through reproducing the full organelle RNA sequencing analysis pipeline described in:

**Hua, Z. (2025).** _Rapid and Cost-Effective Digital Quantification of RNA Editing and Maturation in Organelle Transcripts_. *The Plant Journal* (in review).

All scripts and data files referenced are provided in this repository.

---

## 📁 1. Clone the Repository

```bash
git clone https://github.com/hua-lab/Organelle_RNA_sequencing_analysis.git
cd Organelle_RNA_sequencing_analysis


## ⚙️ 2. Run Method S1 – Demultiplex and Strand-Correct FASTQ Files
```bash
mv *.* ./raw_FASTQs
cd raw_FASTQs
perl filter_barcoded_fastq_seqs.pl

## 🧬 3. Run Method S2 – Align Strand-Corrected Reads to the Pseudo-Genome

chmod +x ./run_minimap_alignments.sh
./run_minimap_alignments.sh

## 🧪 4. Run Method S3 – Extract and Quantify Intron-Retaining Transcripts

# First, move the necessary files into the minimap_alignments folder:

mv count_groupII_inserts_and_extract.py ./minimap_alignments
mv ndhB_intron.fasta ./minimap_alignments
mv local_align_inserts_vs_intron_blast.pl ./minimap_alignments
cd minimap_alignments

# Then, set up and activate a virtual Python environment:

python3 -m venv venv
source venv/bin/activate
pip install biopython pandas

# Now run the script:

python count_groupII_inserts_and_extract.py

## 🔍 5. Run Method S4 – BLAST Analysis of Unspliced Reads

perl local_align_inserts_vs_intron_blast.pl

## ✅ Outputs

ndhB_groupII_insertion_count_chr1-3.tsv: Table summarizing intron-retention frequencies
insert_fastas/: Per-sample FASTA files containing unspliced reads
merged_fastas/: Merged genotype-wise FASTA files for downstream alignment or BLAST
ndhB_insert_vs_intron_blast.tsv: Summary of local alignment results

