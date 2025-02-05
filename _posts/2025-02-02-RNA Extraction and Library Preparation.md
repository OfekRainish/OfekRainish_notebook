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


## 4. Library Preparation
After selecting the RNA samples for sequencing, we send the samples (12) to the [Azrieli Technion Genome Center](https://tgc.net.technion.ac.il/) for library preparation. RNA sequencing libraries will be constructed simultaneously using the [NEBNext UltraExpress RNA Library Prep Kit for Illumina (NEB, Cat. No. E3330)](../pdf%20protocols/PB30.11-UltraScript-cDNA-Synthesis-Kit-Manual.pdf) according to the manufacturer's protocol.

The library preparation steps include:

### Step 1 – rRNA Depletion
The majority of RNA in the sample is ribosomal RNA (rRNA). To remove rRNA from the starting material (200 ng), we will use the [**Ribo-Zero Plus rRNA Depletion Kit**](../pdf%20protocols/illumina-stranded-total-rna-prep-data-sheet-m-gl-02148.pdf)
### Step 2 – Reverse Transcription to cDNA
To convert the RNA into a format suitable for sequencing, we will use reverse transcriptase and other enzymes to synthesize complementary DNA (cDNA) from the RNA template. 2.5 microliter from each sample. (find the whole protocol [here](../pdf%20protocols/PB30.11-UltraScript-cDNA-Synthesis-Kit-Manual.pdf))

### Step 3 – Adapter Ligation
Adapters will be added to the cDNA fragments. These adapters serve several functions:
- **Barcode**: To uniquely identify samples.
- **Flow cell attachment**: To allow the cDNA fragments to bind to the sequencing flow cell.
- **Enable PCR and sequencing**: The adapters allow for both amplification (PCR) and sequencing of the cDNA fragments.

### Step 4 – PCR Amplification
The cDNA library will undergo PCR amplification to generate sufficient quantities of the library for sequencing.

### Step 5 – Library Quality Testing
To ensure that the library meets the necessary quality standards, the following tests will be performed:
- **Library concentration** will be measured using [the Qubit dsDNA HS Assay Kit](../pdf%20protocols/Qubit_dsDNA_HS_Assay_UG.pdf).
- **Library size distribution** will be analyzed using the [TapeStation 4200 with the High Sensitivity D1000 kit](../pdf%20protocols/ScreenTape_HSD1000_QG.pdf).

Results of the quality testing will be shown below:

- [library QC results1](../pdf%20files/qc%20rna%20library/2024-12-25%20-%20TalLuzato_RNA_libs_2,4.pdf)

- [library QC results2](../pdf%20files/qc%20rna%20library/2024-12-26%20-TalLuzato_RNA_libs.pdf)

### Step 6 – Pooling Libraries
Once the library quality is verified,a pool of all the libraries in equal molar ratios was created, in order to create a balanced sample pool. This ensures that all samples are represented accurately during sequencing.

### Step 7 – Sequencing
Sequencing will be performed on the Illumina NextSeq2000 using the [P2 100-cycle kit (Illumina, Cat. No. 20046811)](../pdf%20protocols/nextseq-1000-2000-spec-sheet-m-na-00008.pdf).

Sequencing will provide the data necessary for downstream analysis.
