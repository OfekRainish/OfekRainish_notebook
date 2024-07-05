# Protocol for Analyzing qPCR Results from Exel file

## Overview

This protocol outlines the steps for analyzing quantitative PCR (qPCR) data to determine the relative gene expression levels between control and treated samples. The process involves calculating ΔCt, ΔΔCt, and fold change values using ubiquitin as a reference gene.

## Materials

- Excel table with Ct (Cycle Threshold) values


## Data Description

You are provided with an Excel table containing Ct values, where:
- Columns represent different genes.
- Rows correspond to different samples (control group without treatment and treatment group with an inhibitor).
- Ubiquitin is used as a reference gene due to its stable expression under the tested conditions.

## Steps

### Step 1: Calculate ΔCt Values

**Objective**: Normalize the Ct values of the target genes relative to the reference gene (ubiquitin in this case). This normalization corrects for sample-to-sample variation in qPCR efficiency and total RNA amount.

**Procedure**:
1. Create a new table to record ΔCt values.
2. For each gene and each sample, subtract the Ct value of the reference gene (ubiquitin) from the Ct value of the target gene.

   $$\Delta Ct_{\text{gene}} = Ct_{\text{gene}} - Ct_{\text{reference}}$$

   **Example**: For gene A after treatment:
   
   $$\Delta Ct_{A, \text{treatment}} = Ct_{A, \text{treatment}} - Ct_{\text{ubiquitin, treatment}}$$

**Note**: This normalization is critical because it accounts for variations in RNA input and reverse transcription efficiency.

### Step 2: Calculate ΔΔCt Values

**Objective**: Determine the relative expression change of each gene before and after treatment.

**Procedure**:
1. Create another table to record ΔΔCt values.
2. For each gene, calculate the difference between the ΔCt value of the treated sample and the ΔCt value of the control sample.

   $$\Delta \Delta Ct_{\text{gene}} = \Delta Ct_{\text{gene, treatment}} - \Delta Ct_{\text{gene, control}}$$

   **Example**: For gene A:
   
   $$\Delta \Delta Ct_A = \Delta Ct_{A, \text{treatment}} - \Delta Ct_{A, \text{control}}$$

### Step 3: Calculate Fold Change in Gene Expression

**Objective**: Quantify the relative change in gene expression between treated and control samples using the fold change method.

**Procedure**:
1. Use the following formula to convert ΔΔCt values to fold change values:

   $$\text{Fold Change} = 2^{-\Delta \Delta Ct}$$

   **Example**: For gene A:
   
   $$\text{Fold Change}_A = 2^{-\Delta \Delta Ct_A}$$

**Explanation**: The fold change represents how many times the gene expression has increased or decreased in the treated sample relative to the control. The base of 2 is used because qPCR amplification is exponential and typically, each cycle doubles the amount of the target DNA.

![alt text](../images/qpcr%20exel%20analysis.png)

### Step 4: Graphical Representation of Results

**Objective**: Visualize the fold change in gene expression for each gene.

**Procedure**:
1. Using the fold change values (step 3), create a table suitable for graphing.
2. Select a bar graph format to display the data.
3. Adjust the y-axis settings:
   - Set the minimum value to 1.0 to effectively represent fold change on a logarithmic scale, where 1.0 indicates no change in expression.
4. Label the graph appropriately with gene names on the x-axis and fold change on the y-axis.

![alt text](../images/fold%20change%20graph%20qpcr.png)

## Conclusion

By following these steps, you can accurately analyze and visualize the relative gene expression levels in control and treated samples using qPCR data. This analysis provides insights into the effects of the treatment on the expression of specific genes.
