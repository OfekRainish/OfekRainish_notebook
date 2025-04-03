# RNA-Seq Differential Expression Analysis with DESeq2
For this stage i used R studios.
## 1. Set Working Directory and Load Data

First, set the working directory and load the count and metadata files.

```r
# Set working directory
setwd("~/path/to/deseq2/input/directory")

# Load count data (Counts.csv)
f_counts <- "/path/to/deseq2/input/directory/Counts.csv"
counts1 <- read.csv(f_counts, sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# Load metadata (metadata.csv)
f_design <- "/path/to/deseq2/input/directory/metadata.csv"
design1 <- read.csv(f_design, sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
```
## 2. Remove "GeneID" as a column from Counts
you need to make sure the rows of the factors in the metadata file (TimePint & Treatment) match the columns in the counts file. these will be the names of your samples. additionally, they should appear in the **same order**.

in order for the to match, the GeneID columns must be removed from the counts file (i also think you can not give this column a name from the getgo but i havent checked it)

```r
# removefirst column in counts. if your geneID column isnt the first one, adust the code

counts3 = counts3[, -1]
```


## 3. Ensure Columns Match Rows and Convert Strings to Factors
Ensure that the columns in the count data and rows in the metadata match. Convert metadata variables to factors, as DESeq2 requires factors, not strings, for experimental design.

```r
#should return TRUE if both match and are in the same order
#design1->my metadata, counts1-> my counts
all(rownames(design1) == colnames(counts1))
# Convert metadata (in my case - design1) to factors
colData <- as.data.frame(apply(design1, c(1, 2), as.factor))
```
## 4. Create Model Matrix (optional)
Sometimes, people create a model matrix to check how their design is structured before running DESeq2, i.e. how we divide the data (in this case by time points and with and without surfactant treatment). In our case, this step is just to make sure that the system knows our factors, and is **not mandatory**.

```r
# Create model matrix for experimental design
modelMatrixTest <- as.data.frame(model.matrix(~ treatment + TimePoint, data = colData))

# View model matrix
View(modelMatrixTest)

```

## 5. Filter Genes with Low Counts
Filter out genes with low counts, retaining only those with a count greater than 5 in at least 2 samples. It will go through each row (which represents a gene, in this case) and if at least two samples (2/12) have a read count higher than 5, the gene will be kept. If there are not at least two such samples - the gene will be filtered out.

```r
# Filter genes with low counts (at least 2 samples with count > 5)
keep = rowSums(counts1 > 5) >= 2
counts3 = counts1[keep, ]
```

## 6. Create DESeqDataSet
Create a DESeqDataSet object that contains the count data, the metadata, and the experimental design. tihs dataset combines everything deseq2 needs, and we ill run deseq2 on this dataset.

```r
# Load DESeq2 (and make sure its installed)
library(DESeq2)

# Create DESeqDataSet
dds = DESeqDataSetFromMatrix(
  countData = counts3,   # Use filtered count data
  colData = colData,     # Metadata
  design = ~ TimePoint + Treatment   # Experimental design formula (your factors)
)
```
## 7. Set Reference Level for Treatment
DESeq2 compares gene expression between groups, and it needs to know which group is the "baseline" for comparisons. This is called the reference level.
in my case, the 'control' (no surfactin treatment) will be that baseline.

```r
# Check current levels of 'treatment'. should be the name of the chosen column in the metadata file (colData).
levels(dds$Treatment) 

# Set 'control' as the reference level
dds$treatment <- relevel(dds$Treatment, ref = "control")
```
 ## 8. Initial PCA to Check Data
Perform an initial Principal Component Analysis (PCA) to visualize variation in the data and check for batch effects or outliers.

✅ Good sign: Samples cluster by treatment (e.g., "Social" and "Isolated" separate).                                                             
❌ Bad sign: Samples cluster by something unrelated (e.g., sequencing date).

PCA helps us spot issues early before running DESeq2.

```r
vcd <- vst(dds, blind = FALSE)
PlotPCA(vcd, intgroup = C("TimePoint","Treatment")) #your factors
```


![](../images/rna_bioinformatics/deseq2/initial%20pca%20plot%20-%20before%20running%20deseq.png)



## 9. Run DESeq2
Run DESeq2 to perform differential expression analysis.
```r
# Run DESeq2 analysis
dds <- DESeq(dds)
```
