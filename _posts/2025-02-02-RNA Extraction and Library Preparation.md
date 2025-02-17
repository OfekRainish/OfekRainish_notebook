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
After selecting the RNA samples for sequencing, we send the samples (12) to the [Azrieli Technion Genome Center](https://tgc.net.technion.ac.il/) for library preparation. RNA sequencing libraries will be constructed simultaneously using the [NEBNext UltraExpress RNA Library Prep Kit for Illumina (NEB, Cat. No. E3330)](../pdf%20protocols/manualE3330%20protocol%20-library%20preparation.pdf) according to the manufacturer's protocol.

The library preparation steps include:

### Step 1 – rRNA Depletion
The majority of RNA in the sample is ribosomal RNA (rRNA). To remove rRNA from the starting material (200 ng), we will use the [**Ribo-Zero Plus rRNA Depletion Kit**](../pdf%20protocols/illumina-stranded-total-rna-prep-data-sheet-m-gl-02148.pdf)

### Step 2 - RNA Fragmentation
RNA molecules are often too long for efficient sequencing. The sequencing machine reads short fragments (~200 bp), so RNA must be broken into smaller pieces.

In this protocol, fragmentation occurs after mRNA isolation using the NEBNext UltraExpress RNA Fragmentation Mix, which is heat-activated (94°C for 15 minutes). read more in the [full protocol](../pdf%20protocols/manualE3330%20protocol%20-library%20preparation.pdf)



### Step 3 – Reverse Transcription to cDNA (chek with tecnion what wad dne here)
To convert the RNA into a format suitable for sequencing, we will use reverse transcriptase and other enzymes to synthesize complementary DNA (cDNA) from the RNA template. 2.5 microliter from each sample. (find the whole protocol [here](../pdf%20protocols/PB30.11-UltraScript-cDNA-Synthesis-Kit-Manual.pdf))

### Step 4 – Adapter Ligation
Now that we have double-stranded cDNA, we need to prepare it for sequencing by adding adapters. Adapters are short DNA sequences that serve several important functions:

* They allow cDNA to attach to the sequencing flow cell.
* They contain barcodes (indexes) for identifying   different samples.
* They enable PCR amplification.

The process involves:

1. End Repair & dA-Tailing – The dsDNA ends are modified to create sticky ends that help the adapters attach. This is done using the NEBNext UltraExpress End Prep Enzyme Mix.
2. Adapter Ligation – The NEBNext Adaptor for Illumina is ligated (attached) to the cDNA fragments using the NEBNext UltraExpress Ligation Master Mix.
3. USER Enzyme Treatment – This removes unwanted parts of the adapter to ensure clean ligation.

At this stage, the cDNA is fully prepared with adapters and ready for amplification.

### Step 5 – PCR Amplification (Library Enrichment)
Because the amount of cDNA at this point is very low, we need to make many copies of it using PCR (Polymerase Chain Reaction).

The PCR reaction:

1. Uses the NEBNext MSTC High Yield Master Mix, which contains DNA polymerase and primers that recognize the adapter sequences.
2. Runs 12 cycles of PCR, doubling the number of DNA copies each cycle.

This step ensures that we have enough DNA for sequencing and helps eliminate any unligated adapters.

### Step 6 – Library Quality Testing
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
