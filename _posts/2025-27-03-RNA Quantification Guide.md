# RNA-Seq Quantification Guide

## Introduction
This guide assumes that your libraries have already passed quality checks and trimming. The next step will be to associate (map) regions in the genome and quantify how many reads are associated with each region (gene).

## Step A – Mapping

The mapping, since we are dealing with bacteria without alternative splicing, is done with a tool called **Bowtie**. What Bowtie does is go through our read sequences and give each read an ID card: which sequence it corresponds to in the genome, the starting nucleotide number, the ending nucleotide number, and information about the quality of the alignment.

To do this, it needs a reference genome. Therefore, one of the inputs that needs to be put into the Bowtie code is the genome of the organism from which the RNA sequences were taken. Bowtie does not work directly with the genomic FASTA file since it is cumbersome and too long. Instead, the FASTA file needs to go through a process called **indexing**, in which it is divided and compressed. The results of the indexing are **EBWT** files that cannot be opened manually (but Bowtie can process them).

Now, Bowtie has everything it needs: RNA reads and an indexed reference genome. You can run the following code (just make sure to change the file paths):

```
#!/bin/bash
set -e

# Define paths
GENOM_DIR="/pathway/to/FASTA/directory"
GENOME_FILE="/pathway/to/FASTA/file.fna"
INDEX_PREFIX="/pathway/to/save/genome/index"
READS_DIR="/pathway/to/reads/directory"
OUTPUT_DIR="/pathway/to/bowtie_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "=========== Step 1: Building Bowtie Index ==========="

# Create a Bowtie index from the FASTA genome file
bowtie-build "$GENOME_FILE" "$INDEX_PREFIX"

echo "Indexing completed successfully! (.ebwt files created in $OUTPUT_DIR)"
echo "--------------------------------------------------------------------------------"
echo "=========== Step 2: Aligning RNA-seq Reads ===========" 

# Loop through all .fq files in the trimmed reads directory and align them
for FILE in "$READS_DIR"/*.fq; do
    BASENAME=$(basename "$FILE" .fq) # Extract the filename without extension
    OUTPUT_SAM="$OUTPUT_DIR/${BASENAME}.sam" # Define output SAM filename
    echo "Processing $FILE..."

    # Run Bowtie for single-end reads
    bowtie -x "$INDEX_PREFIX" -q "$FILE" -S "$OUTPUT_SAM"

    echo "Alignment completed for $FILE -> Output: $OUTPUT_SAM"
done

echo "--------------------------------------------------------------------------------"
echo "All RNA-seq files have been aligned successfully!"
echo "========== Bowtie Workflow Completed! ============="

```

The code will save the output of Bowtie in the `OUTPUT_DIR` folder according to the path you entered. The output of Bowtie consists of `.sam` files.

---
## Step B - Quantification

The tool used for quantification is **bedcov**, which counts how many nucleotides are mapped (via Bowtie) to each defined region in the genome. To use bedcov efficiently, we need to prepare two inputs:

### 1. Preparing the BED File
The BED file contains three columns:
- **Contig name** (in a closed genome, there is usually only one value here)
- **Starting nucleotide number**
- **Ending nucleotide number**

Example:
![](../images/rna_bioinformatics/mapping&quantification/Screenshot%202025-03-27%20142927.png)

To create a BED file, you can download a **GBFF/GTF** file from NCBI, which provides information about gene locations and annotations. It is recommended to download the **genomic FASTA** (for Bowtie) and the **GBFF** from the same source to ensure matching contig names between the Bowtie results and the GBFF/GTF file.

Once you have the file, open it in Excel and remove irrelevant columns. You should be left with only:
- Contig name
- Start nucleotide number
- End nucleotide number

Next, convert the Excel file into BED format using the following Python script:

```
import pandas as pd

# Specify the path to your Excel file
excel_file_path = 'pathway/to/excel/file.xlsx'  # Replace with the actual path

# Read the Excel file
df = pd.read_excel(excel_file_path)

# Ensure correct column names
required_columns = ["contig", "start", "end"]
if not all(col in df.columns for col in required_columns): 
    raise ValueError(f"Excel file must contain columns: {required_columns}")

# Convert DataFrame to a BED file format
with open("genes_annotation.bed", "w") as outfile: 
    for _, row in df.iterrows(): 
        contig = row["contig"] 
        start = row["start"] 
        end = row["end"] 

        # BED format: contig start end (tab-separated) 
        bed_line = f"{contig}\t{start}\t{end}\n" 
        outfile.write(bed_line)

print("BED file 'genes_annotation.bed' has been created successfully.")

```

### 2. Preparing the Read Files
bedcov cannot accept the `.sam` files that are the output of Bowtie. These need to go through a three-step process:
1. **Conversion to BAM format** – BAM files are binary and more compact, making them easier to process.
2. **Sorting** – Sorting arranges sequences in genomic order so that overlapping or nearby sequences are grouped together. This improves efficiency when bedcov processes them.
3. **Indexing** – Indexing creates a reference table that allows tools to quickly access parts of the BAM file without scanning the entire file.

Run the following code to process your `.sam` files:

```
#!/bin/bash

# Define directories
SAM_DIR="/pathway" # Folder A: SAM files (Bowtie output)
QUANTIFICATION_DIR="/pathway"
BAM_DIR="/pathway" # Folder B: Will store BAM files and indexes
RESULTS_DIR="/pathway" # Folder C: Will store coverage results
BED_FILE="/pathway/genes_annotation_gbff.bed" # BED file with gene coordinates

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
    samtools sort -o "$BAM_DIR/${BASENAME}_sorted.bam" "$BAM_DIR/${BASENAME}.bam"
    mv "$BAM_DIR/${BASENAME}_sorted.bam" "$BAM_DIR/${BASENAME}.bam"

    # Index BAM file
    echo "Indexing $BAM_DIR/${BASENAME}.bam..."
    samtools index "$BAM_DIR/${BASENAME}.bam"

    # Run bedcov to quantify coverage
    echo "Running samtools bedcov for $BAM_DIR/${BASENAME}.bam..."
    samtools bedcov "$BED_FILE" "$BAM_DIR/${BASENAME}.bam" > "$RESULTS_DIR/${BASENAME}_coverage.txt"

    echo "Finished processing $BASENAME"
done

echo "All files processed. Results are saved in $RESULTS_DIR."

```

The output of bedcov is a text file containing four columns:
- Contig name
- Starting nucleotide number
- Ending nucleotide number
- The total number of nucleotides mapped to that region of the genome

**Note:** This is **not** the number of reads.

---
## Step C – Calculating the Number of Reads & Adding Annotations

### 1. Calculating the Read Count
Open the text file in Excel and add two new columns:
- **Region length** (calculated as `end - start`)
- **Average cover (or sum all reads)** – This represents the average number of nucleotides per position in the region, which is a proxy for read count.

Formula:
```
Average cover = (total nucleotides mapped) / (region length)
```

Example output:

![](../images/rna_bioinformatics/mapping&quantification/Screenshot%202025-03-27%20143318.png)

This is not a precise read count, but it is sufficient for comparing libraries and treatments.

### 2. Adding Annotations
To add gene annotations, merge this table with the original GBFF/GTF table. Since the contig name, start, and end positions are the same in both tables, they can be easily matched.

I tested 12 libraries and merged them into one final table, which you can see [here](../exel%20files/quantification/overall_quantification_summary.csv).



## Summary
By following these steps, you can:
1. Map RNA reads to the genome using Bowtie.
2. Prepare the BED file and process `.sam` files for quantification.
3. Use bedcov to quantify nucleotide counts in genomic regions.
4. Estimate the number of reads per region and add annotations.

This process enables efficient RNA-Seq quantification for bacterial datasets.

