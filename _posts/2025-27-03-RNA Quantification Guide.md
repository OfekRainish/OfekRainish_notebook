# RNA-Seq Quantification Guide

## Overview
This guide assumes that your libraries have already passed quality checks and trimming. The next step will be to associate (map) reads to regions in the genome and quantify how many reads are associated with each gene.

## Step A – Mapping with Bowtie
Since we are working with bacteria (which lack alternative splicing), we use **Bowtie** for mapping. Bowtie assigns each read to a genomic location by providing:
- The corresponding sequence in the genome
- The starting and ending nucleotide positions
- Quality information on the alignment

### Inputs:
- RNA-seq reads (FASTQ format)
- Indexed reference genome

### Preparing the Genome Index
Bowtie requires an **indexed reference genome** instead of a raw FASTA file. Indexing compresses and structures the FASTA file into `.ebwt` files.

```bash
#!/bin/bash
set -e

# Define paths
GENOM_DIR="/path/to/FASTA/directory"
GENOME_FILE="/path/to/FASTA/file.fna"
INDEX_PREFIX="/path/to/save/genome/index"
READS_DIR="/path/to/reads/directory"
OUTPUT_DIR="/path/to/bowtie_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "=========== Step 1: Building Bowtie Index ==========="

# Create a Bowtie index from the FASTA genome file
bowtie-build "$GENOME_FILE" "$INDEX_PREFIX"

echo "Indexing completed successfully!"
echo "--------------------------------------------------------------------------------"
echo "=========== Step 2: Aligning RNA-seq Reads ==========="

# Loop through all .fq files in the reads directory
for FILE in "$READS_DIR"/*.fq; do 
    BASENAME=$(basename "$FILE" .fq) # Extract filename without extension 
    OUTPUT_SAM="$OUTPUT_DIR/${BASENAME}.sam" # Define output SAM filename 
    echo "Processing $FILE..." 

    # Run Bowtie for single-end reads 
    bowtie -x "$INDEX_PREFIX" -q "$FILE" -S "$OUTPUT_SAM" 

    echo "Alignment completed for $FILE -> Output: $OUTPUT_SAM"
done

echo "--------------------------------------------------------------------------------"
echo "All RNA-seq files have been aligned successfully!"
echo "========== Bowtie Workflow Completed! =============="
```

The output of Bowtie is `.sam` files stored in `OUTPUT_DIR`.

---

## Step B - Quantification with FeatureCounts
To quantify gene expression, we use **featureCounts**, which maps each BAM file to genes or genomic regions using a GTF annotation file.

### Inputs:
- Mapped reads (`.bam` files)
- Gene annotation file (`.gtf`)

### Preparing BAM Files
FeatureCounts does not accept `.sam` files directly. They must be converted, sorted, and indexed:

```bash
#!/bin/bash

# Define directories
SAM_DIR="/path/to/bowtie_output" # Folder A: SAM files (Bowtie output)
BAM_DIR="/path/to/bam_files" # Folder B: Stores BAM files and indexes
RESULTS_DIR="/path/to/featureCounts_results" # Folder C: Stores count results
GTF_FILE="/path/to/gtf/file.gtf" # Replace with actual path
THREADS=4

# Create necessary directories
mkdir -p "$BAM_DIR"
mkdir -p "$RESULTS_DIR"

# Loop over all SAM files
for SAM_FILE in "$SAM_DIR"/*.sam; do 
    BASENAME=$(basename "$SAM_FILE" .sam) 

    # Convert SAM to BAM 
    echo "Converting $SAM_FILE to BAM..." 
    samtools view -b -o "$BAM_DIR/${BASENAME}.bam" "$SAM_FILE" 

    # Sort BAM file 
    echo "Sorting $BAM_DIR/${BASENAME}.bam..." 
    samtools sort -o "$BAM_DIR/${BASENAME}.sorted.bam" "$BAM_DIR/${BASENAME}.bam" 
    rm "$BAM_DIR/${BASENAME}.bam" # Remove unsorted BAM

    # Index BAM file 
    echo "Indexing $BAM_DIR/${BASENAME}.sorted.bam..." 
    samtools index "$BAM_DIR/${BASENAME}.sorted.bam"
done

# ===============================
# Step 1: Run FeatureCounts
# ===============================
echo "Running FeatureCounts..."
featureCounts -T "$THREADS" -t gene -g gene_id -F GTF -a "$GTF_FILE" -o "$RESULTS_DIR/counts.txt" "$BAM_DIR"/*.sorted.bam

if [ $? -eq 0 ]; then
    echo "FeatureCounts analysis complete. Results saved in $RESULTS_DIR/counts.txt"
else
    echo "Error: FeatureCounts failed."
    exit 1
fi
```

### Output:
FeatureCounts generates a text file (`counts.txt`) containing gene expression data with:
- A **GeneID** column
- **n sample columns**, where n is the number of BAM files

---

## Step C – Adding Annotations
To annotate the output, merge the **featureCounts** result with the original **GTF/GBFF file** by matching:
- Contig name
- Start and end nucleotide positions

> **Note:** The GTF file contains both gene and non-gene regions, so filtering may be necessary to align the annotations correctly.

**Example Table (Final Merged Output):**
| GeneID | Count_Sample1 | Count_Sample2 | Annotation |
|--------|--------------|--------------|------------|
| GeneA  | 105          | 98           | Metabolic Enzyme |
| GeneB  | 50           | 72           | Transcription Factor |

---

## Summary
1. **Index the genome** (`bowtie-build`)
2. **Map reads to genome** (`bowtie` → `.sam`)
3. **Convert, sort, and index BAM files** (`samtools`)
4. **Quantify gene expression** (`featureCounts`)
5. **Annotate results** (`merge with GTF`)

This workflow provides a robust pipeline for RNA-Seq quantification in bacteria using Bowtie and featureCounts.

---

📌 **Notes:**
- Ensure FASTA and GTF files are from the same source (e.g., NCBI) to avoid contig mismatches.
- Adjust paths and thread counts (`-T`) based on your system capabilities.

🛠 Happy RNA-seq analysis! 🚀