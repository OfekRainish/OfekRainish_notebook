# RNA Extraction and Library Preparation for RNA Sequencing

## 1. Sample Collection
Bacterial samples were collected at two different time points for this experiment:
- 5 samples after 20 hours with surfactin
- 5 samples after 20 hours without surfactin
- 5 samples after 44 hours with surfactin
- 5 samples after 44 hours without surfactin

A total of 20 samples were collected. RNA extraction and sequencing will be performed on triplicate samples from each group, selecting the one with the highest RNA quality for sequencing.

## 2. RNA Extraction
RNA will be extracted from the bacterial samples using the NucleoSpin® RNA Stool Kit (according to the manufacturer’s protocol). The procedure can be accessed in detail [here](../_posts/2024-07-24-RNA%20Extraction%20Protocol.md).

## 3. RNA Quality Testing
To ensure the RNA is of sufficient quality for sequencing, we will assess the RNA samples using two methods in our laboratory:
1. **NanoDrop Spectrophotometer**: To assess RNA concentration and purity.
2. **TapeStation 4200**: To determine RNA integrity and size distribution.

The samples selected for sequencing will are highlighted in green.

![image1](../images/rna%20extract/lab%20QC.png)

[raw data](../exel%20files/)

## 4. Library Preparation
After selecting the RNA samples for sequencing, we will send the samples to the Azrieli Technion Genome Center for library preparation. RNA sequencing libraries will be constructed simultaneously using the NEBNext UltraExpress RNA Library Prep Kit for Illumina (NEB, Cat. No. E3330) according to the manufacturer's protocol.

The library preparation steps include:

### Step 1 – rRNA Depletion
The majority of RNA in the sample is ribosomal RNA (rRNA). To remove rRNA, we will use the **Ribo-Zero Plus rRNA Depletion Kit**. After rRNA removal, we will purify the RNA, retaining 200 ng for the subsequent steps.

### Step 2 – Reverse Transcription to cDNA
To convert the RNA into a format suitable for sequencing, we will use reverse transcriptase and other enzymes to synthesize complementary DNA (cDNA) from the RNA template.

### Step 3 – Adapter Ligation
Adapters will be added to the cDNA fragments. These adapters serve several functions:
- **Barcode**: To uniquely identify samples.
- **Flow cell attachment**: To allow the cDNA fragments to bind to the sequencing flow cell.
- **Enable PCR and sequencing**: The adapters allow for both amplification (PCR) and sequencing of the cDNA fragments.

### Step 4 – PCR Amplification
The cDNA library will undergo PCR amplification to generate sufficient quantities of the library for sequencing.

### Step 5 – Library Quality Testing
To ensure that the library meets the necessary quality standards, the following tests will be performed:
- **Library concentration** will be measured using the Qubit dsDNA HS Assay Kit.
- **Library size distribution** will be analyzed using the TapeStation 4200 with the High Sensitivity D1000 kit.

Results of the quality testing will be shown below:

*(Leave space for image here)*

### Step 6 – Pooling Libraries
Once the library quality is verified, we will pool all the libraries in equal molar ratios to create a balanced sample pool. This ensures that all samples are represented accurately during sequencing.

### Step 7 – Sequencing
Sequencing will be performed on the Illumina NextSeq2000 using the P2 100-cycle kit (Illumina, Cat. No. 20046811). The sequencing run will follow the protocol for the P2 100-cycle kit, which includes:
- **Read 1**: 100 cycles
- **Index 1**: 8 cycles
- **Index 2**: 8 cycles
- **Read 2**: 0 cycles

Sequencing will provide the data necessary for downstream analysis.
