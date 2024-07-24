# Protocol for Loading qPCR Plates

This protocol outlines the steps to prepare stocks and load a qPCR plate for testing gene expression using SYBR Green and specific primers.

---

## Preparation of Stocks

### 1. Primer Stock for One Gene

Each primer stock corresponds to a specific gene and includes both the forward and reverse primers, diluted in water, combined with SYBR Green.

**Preparation Steps:**

1. **Prepare Primer Solution:**
   - Mix 10 µL of forward primer and 10 µL of reverse primer from their original 100 µM stocks.
   - Add 810 µL of nuclease-free water to dilute the primers into a working solution.

2. **Calculate Primer Volume for qPCR:**
   - Determine the number of wells needed for this primer pair (gene).
   - For each well, you will need 3 µL of the diluted primer solution.
   - Multiply 3 µL by the total number of wells to get the required volume of primer solution.

3. **Add SYBR Green:**
   - You will need to add 3 µL of SYBR Green for each well.
   - Multiply 3 µL by the total number of wells to get the required volume of SYBR Green.

4. **Prepare Final Primer Stock:**
   - Combine the calculated volumes of the diluted primer solution and SYBR Green to form the primer stock for one gene.

**Note:** Repeat these steps to prepare separate primer stocks for each gene you are testing.

### 2. cDNA Stock for One Sample

Prepare a cDNA stock for each sample by combining the cDNA solution with SYBR Green.

**Preparation Steps:**

1. **Calculate cDNA Volume for qPCR:**
   - Determine the number of wells needed for each cDNA sample.
   - For each well, you will need 3 µL of the cDNA solution.
   - Multiply 3 µL by the total number of wells to get the required volume of cDNA solution.

2. **Add SYBR Green:**
   -You will need to add another 3 µL of SYBR Green for each well.
   - Multiply 3 µL by the total number of wells to get the required volume of SYBR Green.

3. **Prepare Final cDNA Stock:**
   - Combine the calculated volumes of the cDNA solution and SYBR Green to form the cDNA stock for one sample.

**Note:** Repeat these steps to prepare separate cDNA stocks for each sample.

---

## Loading the qPCR Plate

Each well of the qPCR plate will contain 10 µL of the combined stock solutions.

1. **Plate Configuration:**
   - Use a 96-well plate (8 rows x 12 columns).
   - Assign rows to represent different genes (primers).
   - Assign columns to represent different samples (cDNA).

2. **Loading Steps:**
   - First, pipette 5 µL of the primer stock into each well.
   - Then, add 5 µL of the cDNA stock into each well.

3. **Control Wells:**
   - Include negative controls by adding water instead of cDNA.
   - Include RNA controls to check for genomic DNA contamination.

4. **Technical Repeats:**
   - Perform triplicate tests for each gene-sample combination to ensure accuracy.

### Example Plate Setup

![alt text](../images/qpcr%20loading.png)

---

## Important Considerations

- **Contamination Checks:**
  - Include water controls to ensure no contamination.
  - Include RNA controls to verify the absence of genomic DNA.

- **Plate Setup:**
  - Clearly label wells to identify gene-sample combinations.
  - Ensure precise pipetting to maintain the consistency of 10 µL total volume per well.

---

## Conclusion

By following this protocol, you can systematically prepare and load qPCR plates, ensuring accurate and reproducible results for your gene expression studies. Adjust as needed for different plate sizes or experimental designs.

