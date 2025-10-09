# 📘 qPCR Reference Gene Selection from RNA-Seq Data

This R script identifies **stable genes** from RNA-Seq (DESeq2) results that can serve as **reference (normalization) genes** for qPCR validation.  
The analysis was performed on *Paenibacillus dendritiformis* grown **with and without surfactin treatment**.

---

## 🧩 Overview

**Goal:**  
To find genes with stable expression (no significant changes between conditions) and adequate expression levels, suitable for qPCR normalization.

**Main Steps:**
1. Load RNA-seq counts and metadata.
2. Prepare data for DESeq2.
3. Filter out low-expression genes.
4. Run DESeq2 differential expression analysis.
5. Select stable genes (non-significant and low fold-change).
6. Annotate genes using the genome GTF file.
7. Rank stability and identify top candidates.
8. Compare against known housekeeping genes.

---
the setup code is the one used for all rna seq analysis:
```r
# --- Setup ---
setwd("~/home/oreinish/RNA_seq/deseq2/dseq2_input/")

# Load counts
f_counts <- "/home/oreinish/RNA_seq/deseq2/dseq2_input/Counts.csv"
counts1 <- read.csv(f_counts, sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# Load metadata
f_design <- "/home/oreinish/RNA_seq/deseq2/dseq2_input/metadata.csv"
design1 <- read.csv(f_design, sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)


# removefirst column in counts. if your geneID column isnt the first one, adust the code

counts3 = counts1[, -1]
colnames(counts1) #check column names of counts1
rownames(design1) #check row names of design1 (the metadata)
colnames(counts1) <- rownames(design1) # make colnames as rownames
identical(colnames(counts1), rownames(design1)) # check if identical

#should return TRUE if both match and are in the same order

#design1->my metadata, counts1-> my counts
all(rownames(design1) == colnames(counts1)) # check if identical

# Convert directly to factors
colData <- design1
colData$Treatment <- factor(colData$Treatment)
colData$TimePoint <- factor(colData$TimePoint)


# Create model matrix for experimental design
modelMatrixTest <- as.data.frame(model.matrix(~ Treatment + TimePoint, data = colData))
# View model matrix

#View(modelMatrixTest)

# Filter genes with low counts
keep = rowSums(counts1 > 5) >= 3 & rowSums(counts1) >100
counts3 = counts1[keep, ]

# Load DESeq2 (and make sure its installed)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts3,
  colData = colData,
  design = ~ TimePoint * Treatment
)

# Check current levels of 'treatment'. should be the name of the chosen column in the metadata file (colData).
levels(dds$Treatment)

# Set 'control' as the reference level
dds$treatment <- relevel(dds$Treatment, ref = "control")

# Stats
vcd <- vst(dds, blind = FALSE) # VST normalization
plotPCA(vcd, intgroup = c("TimePoint","Treatment")) # PCA


# Run DESeq2 analysis
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)

resultsNames(dds) #shows you the comparisons deseq made
compare1 <- resultsNames(dds)[2] # Variable assignment of TimePoint comparison
compare2 <- resultsNames(dds)[3] #Variable assignment of Surfactin treatment comparison

```

the next thing we want to do is to select the genes which expression dont significantly change when adding the surfactin:
* Selects non-significant genes (padj > 0.05).
* Restricts to genes with minimal expression changes (|LFC| < 0.3).
* Ensures moderate expression (baseMean > 100).

These are the most stable genes, ideal for qPCR normalization.
```R
res <- results(dds, name = compare2)
res <- na.omit(as.data.frame(res))
res$gene <- rownames(res)

stable_genes <- res %>%
  dplyr::filter(padj > 0.05 | is.na(padj)) %>%
  dplyr::filter(abs(log2FoldChange) < 0.3) %>%
  dplyr::filter(baseMean > 100) %>%
  dplyr::arrange(abs(log2FoldChange))
```

this next segment is cosnetic, just to get a "clean" gene name:
```R
library(stringr)
stable_genes$gene_clean <- str_extract(stable_genes$gene, "PDENDC454_[0-9]+")
stable_genes$gene_clean[is.na(stable_genes$gene_clean)] <- stable_genes$gene[is.na(stable_genes$gene_clean)]
stable_genes$gene <- stable_genes$gene_clean
```
the following is for annotation.
```R
library(rtracklayer)
gtf <- import("genomic.gtf")
gtf_df <- as.data.frame(gtf)

gtf_slim <- gtf_df %>%
  dplyr::select(seqnames, start, end, strand, gene_id, locus_tag, product) %>%
  distinct(gene_id, .keep_all = TRUE)

annotated_stable <- stable_genes %>%
  left_join(gtf_slim, by = c("gene" = "gene_id"))
```
finally we rank the most stable genes according to the categories mentiond above (fold change, number if reads, p-val), and ask for the top stable genes.
note, this dosnt mean these genes are suitable to be reference genes for qpcr. we still need to make sure these are houskeeping genes.
```R
annotated_stable$stability_score <- (1 / (1 + abs(annotated_stable$log2FoldChange))) * log10(annotated_stable$baseMean + 1)

top_stable <- annotated_stable %>%
  arrange(desc(stability_score)) %>%
  select(gene, locus_tag, baseMean, log2FoldChange, padj, stability_score) %>%
  head(1000)
```
another thing i did is goung the other way around. meaning i got a list of houskeeping genes used priviuosly as reference genes in Paenibacillus. then i looked for them in the rna seq results to see whether they are stable.
the initial list was provided by chat GPT and litriture. 
```R
candidate_genes <- c("PDENDC454_19328", "PDENDC454_09860",
                     "PDENDC454_23409", "PDENDC454_05621",
                     "PDENDC454_20892", "PDENDC454_09875")

candidate_table <- annotated_stable[annotated_stable$gene %in% candidate_genes,
                                    c("gene", "baseMean", "log2FoldChange", "padj", "product")]

candidate_table$stable_expr <- with(candidate_table,
                                    padj > 0.05 & abs(log2FoldChange) < 0.3 & baseMean > 100)

```