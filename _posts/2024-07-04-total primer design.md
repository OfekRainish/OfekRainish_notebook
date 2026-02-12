# Protocol for Designing Primers for **Paenibacillus dendritiformis C454** Genes

## Purpose
This protocol outlines the steps for designing primers for genes in **Paenibacillus dendritiformis**, focusing on genes within biosynthetic gene clusters (BGCs) and physiological genes involved in various cellular functions.

---

## Step 1: Getting the Gene Sequences

### A. Genes from Biosynthetic Gene Clusters (BGCs)
1. **Identification**: Obtain sequences of core genes from BGCs using antiSMASH software.
2. **Selection**: Choose the core gene from each BGC for primer design.

### B. Physiological Genes
1. **Gene Identification**:
   - Identify physiological genes responsible for essential cellular actions (e.g., cell division, movement, DNA replication).
   - Utilize resources like ChatGPT to list known genes in **Paenibacillus dendritiformis**.
   
2. **Sequence Retrieval**:
   - Search for the desired gene on NCBI with the query "Paenibacillus dendritiformis + gene name".
   - Prefer sequences from strain J2TS7 if strain C454 sequences are not available directly.
   
3. **Sequence Validation**:
   - Click on the FASTA sequence for the gene.
   - Select "RUN BLAST" to compare against a database.
   - Choose 'Whole genome shotgun contigs (wgs)' with "Paenibacillus dendritiformis str. C454 (taxid: 1131935)" as the organism.
   
4. **Alignment Check**:
   - Verify that the query sequence aligns with the C454 strain sequence.
   - If matched, copy the sequence from GenBank or directly from the alignment.

5. **Sequence Comparison**:
   - If needed, compare sequences using BioEdit software to ensure similarity.

---

## Step 2: Primer Design

### A. Using Primer3 Software
1. **Input Sequence**: Enter the validated gene sequence into Primer3.
2. **Adjust Parameters**:
   - For qPCR primers, set the product size to 130-160 bp.
   - Generate five pairs of primer candidates.

### B. Primer Evaluation
1. **PCR Primer Stats Check**:
   - Run each primer through PCR Primer Stats software.
   - Ensure the primer passes the basic criteria. A warning for temperatures above 58°C is acceptable, with a normal range even up to 62°C.

2. **Hairpin Formation (Entropy) Test**:
   - Use IDT's OligoAnalyzer to check the ΔG value for hairpin formation.
   - Select primers with the lowest ΔG value, ideally lower than -1.

3. **Specificity Test using PrimerBLAST**:
   - Run PrimerBLAST with the following settings:
     - Organism: Paenibacillus dendritiformis c454 (taxid 1131935).
     - Database: Custom -> Full genome FASTA from NCBI.
   - Choose primer pairs with a single expected match, ensuring the amplicon length matches expectations and the match is complete.

---

## Step 3: Primer List

Prepare a list of primers for each gene, formatted as follows:
<<<<<<< HEAD
=======

>>>>>>> 455a2a65068423c4fa075ae2284bc802f46b75a9
| Serial Number | Gene Name               | Gene Type                        | Forward Primer          | Reverse Primer         |
|---------------|--------------------------|----------------------------------|-------------------------|------------------------|
| 1             | Asparagine synthase       | Natural product                  | TTATTTCTGGGCAGCCTCCA     | GCTGGTCTGCTAATTCGTCC   |
| 2             | Bacitracin synthetase 3   | Natural product                  | GAATGTCGGTTGGAGTACGC     | TTCCTCCTCCGTGAGCATTT   |
| 3             | Isochorismate synthase    | Natural product                  | TCTCGACCAACATAACCGGG     | ATAGTACGCCCGATCGAAGG   |
| 4             | L-ectoine synthase        | Natural product                  | ACAAACATCATGTGGAGGCG     | TCTGGCTCTTTCCTCTCAGC   |
| 5             | dnaA                     | Physiological - DNA replication   | TGGTTCAAAGCCACTCAAGC     | TCCACCTGATTGCCCGTAAT   |
| 6             | dnaG                     | Physiological - DNA replication   | GGCATCTTGGTGAACGGTTT     | TATGGAGGCAGCGACTTTCT   |
| 7             | dnaN                     | Physiological - DNA replication   | TCCAGACCTTTCTCCGTTCC     | GTGGAAATGGCAAAGACCGT   |
| 8             | gyrA                     | Physiological - DNA replication   | CTATGACGGGGAAGAGACGG     | CGCCTGAATGCCATCAATGA   |
| 9             | gyrB                     | Physiological - DNA replication   | TTGAAGTCAGCTCTCTCCCG     | AGATTTTCCCCTTGAGCGGA   |
| 10            | ftsA                     | Physiological - cell division     | AGCTACTTCGACCTTGCCAA     | GCGGACGACTTTGAACACTT   |
| 11            | ftsZ                     | Physiological - cell division     | GGTAAATACGGATGCGCAGG     | AGTTCACGGGACTCTTCAGC   |
| 12            | ftsI                     | Physiological - cell division     | CGGAATCCTGCAATGTCGTC     | GGCATAACGGGTTGAAGCAT   |
| 13            | minC                     | Physiological - cell division     | GCTTGACGATCAGTGCGAAT     | TCTGTCTTTTGCTCATCGGC   |
| 14            | minD                     | Physiological - cell division     | CGTCAAGGATAAGCGGTTCG     | CGGCAGGGCAGTCGATAATA   |
| 15            | motA                     | Physiological - movement          | GATGACTTCCTTCGCAACGG     | CCTGGGAAAAGATAAGCGCG   |
| 16            | motB                     | Physiological - movement          | GCGGGACAAGATGAACGAAT     | GACTTGAGGCTTGCTTTCCC   |
| 17            | mreB                     | Physiological - movement          | ACGCTGTCCATCTCTTCGAA     | GTGAGAACAATGCCCCGATC   |
| 18            | recA                     | Potential reference gene          | GGAATCTCCCATCTCGCCTT     | ATCGACGAACTGCTTCTGTC   |
| 19            | gyrB                     | Potential reference gene          | CAACTTGAGCGGGGATGATG     | GAACAAGGACTCGACGATGC   |
| 20            | dacF                     | Protease                         | CCCACTTCGAGAACAGCAAC     | TGGCGCAGGTAATCTTCGTA   |
<<<<<<< HEAD
=======

>>>>>>> 455a2a65068423c4fa075ae2284bc802f46b75a9


---

## Conclusion

This protocol provides a detailed guide for designing and validating primers for various genes in **Paenibacillus dendritiformis**. Following these steps ensures accurate and efficient primer design for further genetic and functional studies.
