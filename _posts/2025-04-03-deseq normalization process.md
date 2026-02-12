# **DESeq2 Normalization Process**

1. **Geometric Mean Calculation**: For each gene, calculate the **geometric mean** of its raw counts across all samples. This provides a "typical" expression level for each gene.

2. **Calculate Ratios**: For each sample, calculate the **ratio** of the gene's raw count to the gene's geometric mean. This gives an idea of how much the gene's expression deviates from the typical expression across all samples.

3. **Compute Size Factors**: For each sample, calculate the **median ratio** of all genes. This median ratio is called the **size factor** for that sample, and it represents the relative sequencing depth of the sample.

4. **Normalize Counts**: Normalize the raw counts for each gene by dividing the raw counts by the size factor for that sample. This makes the counts across samples comparable, correcting for differences in sequencing depth.

The normalization allows for more accurate comparison of gene expression across samples, ensuring that differences are biological, not due to technical factors like varying sequencing depths.

---

### **Example of DESeq2 Normalization**:

Suppose we have the following raw counts for three genes across three samples:

| Gene | Sample 1 | Sample 2 | Sample 3 |
|------|----------|----------|----------|
| Gene A | 100      | 150      | 120      |
| Gene B | 200      | 300      | 250      |
| Gene C | 50       | 75       | 60       |

1. **Geometric Mean Calculation**: 

   For each gene, calculate the geometric mean of its counts across all samples:

   - **Gene A**: Geometric mean = sqrt(100 * 150 * 120) = 120.0
   - **Gene B**: Geometric mean = sqrt(200 * 300 * 250) = 250.0
   - **Gene C**: Geometric mean = sqrt(50 * 75 * 60) = 61.2

2. **Calculate Ratios**: 

   For each sample, calculate the ratio of the gene's raw count to the gene's geometric mean:

   - **Sample 1**:
     - Gene A ratio = 100 / 120.0 = 0.833
     - Gene B ratio = 200 / 250.0 = 0.8
     - Gene C ratio = 50 / 61.2 = 0.816
     
   - **Sample 2**:
     - Gene A ratio = 150 / 120.0 = 1.25
     - Gene B ratio = 300 / 250.0 = 1.2
     - Gene C ratio = 75 / 61.2 = 1.22
     
   - **Sample 3**:
     - Gene A ratio = 120 / 120.0 = 1.0
     - Gene B ratio = 250 / 250.0 = 1.0
     - Gene C ratio = 60 / 61.2 = 0.980

3. **Compute Size Factors**: 

   For each sample, calculate the **median ratio** across all genes:

   - **Sample 1 size factor** = median(0.833, 0.8, 0.816) = 0.816
   - **Sample 2 size factor** = median(1.25, 1.2, 1.22) = 1.22
   - **Sample 3 size factor** = median(1.0, 1.0, 0.980) = 1.0

4. **Normalize Counts**:

   Finally, normalize the counts by dividing each gene's raw count by the size factor for the corresponding sample:

   - **Normalized counts for Sample 1**:
     - Gene A = 100 / 0.816 = 122.45
     - Gene B = 200 / 0.816 = 245.1
     - Gene C = 50 / 0.816 = 61.2
     
   - **Normalized counts for Sample 2**:
     - Gene A = 150 / 1.22 = 122.95
     - Gene B = 300 / 1.22 = 245.9
     - Gene C = 75 / 1.22 = 61.48
     
   - **Normalized counts for Sample 3**:
     - Gene A = 120 / 1.0 = 120.0
     - Gene B = 250 / 1.0 = 250.0
     - Gene C = 60 / 1.0 = 60.0

---

This normalization process ensures that the counts from all samples are on the same scale, allowing for accurate comparison of gene expression across samples, regardless of differences in sequencing depth.
