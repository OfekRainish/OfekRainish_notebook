# 🧬 RNA-seq Pipeline (Eukaryotes, Single-End 100 bp)

This pipeline processes **single-end RNA-seq data (100 bp reads)** from raw FASTQ files to gene-level counts.

---

# 🔧 STEP 1 — Install Tools

```bash
conda create -n rnaseq_env -y
conda activate rnaseq_env

conda install -c bioconda fastqc multiqc trim-galore star salmon samtools subread -y
```

---

# 📊 STEP 2 — Initial QC

```bash
#!/bin/bash

input_dir="/path/to/fastq_files"
fastqc_output_dir="/path/to/qc_raw/fastqc"
multiqc_output_dir="/path/to/qc_raw/multiqc"

mkdir -p "$fastqc_output_dir"
mkdir -p "$multiqc_output_dir"

fastqc "$input_dir"/*.fastq.gz -o "$fastqc_output_dir"
multiqc "$fastqc_output_dir" -o "$multiqc_output_dir"

echo "Initial QC completed"
```

---

# ✂️ STEP 3 — Trimming (Single-End)

```bash
#!/bin/bash

input_dir="/path/to/fastq_files"
trimmed_dir="/path/to/trimmed_reads"

mkdir -p "$trimmed_dir"

for file in "$input_dir"/*.fastq.gz; do
    trim_galore \
        --length 30 \
        --nextseq 20 \
        "$file" \
        -o "$trimmed_dir"
done

echo "Trimming completed"
```

🧠 Notes:

* `--length 30` keeps reads ≥30 bp (appropriate for 100 bp reads)
* `--nextseq 20` removes poly-G artifacts

---

# 🔁 STEP 4 — QC After Trimming

```bash
#!/bin/bash

trimmed_dir="/path/to/trimmed_reads"
fastqc_output_dir="/path/to/qc_trimmed/fastqc"
multiqc_output_dir="/path/to/qc_trimmed/multiqc"

mkdir -p "$fastqc_output_dir"
mkdir -p "$multiqc_output_dir"

fastqc "$trimmed_dir"/*.fq.gz -o "$fastqc_output_dir"
multiqc "$fastqc_output_dir" -o "$multiqc_output_dir"

echo "Post-trimming QC completed"
```

---

# 🧬 STEP 5 — Mapping (STAR, Single-End)

```bash
#!/bin/bash

genome_fasta="/path/to/genome.fa"
gtf_file="/path/to/annotation.gtf"
trimmed_dir="/path/to/trimmed_reads"
star_index="/path/to/star_index"
mapped_dir="/path/to/post_mapping"

mkdir -p "$star_index"
mkdir -p "$mapped_dir"

# Build index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$star_index" \
     --genomeFastaFiles "$genome_fasta" \
     --sjdbGTFfile "$gtf_file"

# Map reads
for file in "$trimmed_dir"/*_trimmed.fq.gz; do
    base=$(basename "$file" _trimmed.fq.gz)

    STAR --genomeDir "$star_index" \
         --readFilesIn "$file" \
         --readFilesCommand zcat \
         --runThreadN 8 \
         --outFileNamePrefix "$mapped_dir/${base}_" \
         --outSAMtype BAM SortedByCoordinate
done

echo "Mapping completed"
```

---

# 📊 STEP 6 — Quantification (featureCounts, Single-End)

```bash
#!/bin/bash

mapped_dir="/path/to/post_mapping"
counts_dir="/path/to/post_quantification"
gtf_file="/path/to/annotation.gtf"

mkdir -p "$counts_dir"

featureCounts -T 8 \
    -t exon \
    -g gene_id \
    -a "$gtf_file" \
    -o "$counts_dir/counts.txt" \
    "$mapped_dir"/*Aligned.sortedByCoord.out.bam

echo "Quantification completed"
```

---

# ⚡ OPTIONAL — Salmon (Single-End alternative)

```bash
# Index
salmon index -t transcripts.fa -i salmon_index

# Quantification
for file in *.fastq.gz; do
    base=$(basename "$file" .fastq.gz)

    salmon quant -i salmon_index \
        -l A \
        -r "$file" \
        -p 8 \
        -o "salmon_${base}"
done
```

---

# 🎯 Summary

Pipeline flow:

1. QC
2. Trimming
3. QC again
4. STAR mapping
5. featureCounts quantification

OR:

👉 Salmon replaces steps 4–5

---

# ⚠️ Final Notes

* Make sure FASTA and GTF match
* Check strandedness (`-s` in featureCounts)
* Use absolute paths (`pwd`)
* Keep directory structure consistent

---

✅ This pipeline is:

* Correct for single-end RNA-seq
* Clean and reproducible
* Ready for DESeq2

---
