## Primer Design for qPCR and Specificity Testing Using Gel Electrophoresis

Primers were designed to target four natural substances produced by the bacterium *Paenibacillus dendritiformis* strain C454. The biosynthetic gene clusters (BGCs) of interest were identified using the antiSMASH tool. The primers were design to target core genes within the BGCs.

### 1. Primer Design

Primers were designed using the [IDT PrimerQuest Tool](https://www.idtdna.com/Primerquest/Home/Index). The following criteria were applied during the design process:

- **Target Regions:** Selected within the core genes of the BGCs.
- **Melting Temperature (Tm):** Suitable for consistent annealing.
- **GC Content:** Between 40-60% to enhance primer stability.
- **Amplification Product Size:** Ranging between 75-150 base pairs (bp) for optimal qPCR efficiency.

### 2. Primer Validation

After design, primers were validated using the [PCR Primer Stats Tool](http://biotools.nubic.northwestern.edu/OligoCalc.html) to ensure they met standard qPCR parameters. This includes checking for:
- Primer length.
- Melting temperature.
- Lack of secondary structures (hairpins, dimers).

Subsequently, specificity was checked using [Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) against the *Paenibacillus dendritiformis* C454 genome to confirm that each primer pair binds specifically to the intended target sequence without off-target binding.
- Data base: custom -> FASTA (from NCBI)
- Organism: *Paenibacillus dendritiformis c454 (taxid 1131935)*

### 3. Primer Specificity and Gel Electrophoresis Testing

To test the primers for specificity and functionality, PCR amplification was performed, followed by gel electrophoresis. DNA was extracted using two different kits (Kit1 and Kit2) to compare their effectiveness. PCR was conducted with each primer set, including controls for validation.

#### Primers Designed for Specific Genes

| Gene Name                     | Forward Primer (5' → 3')            | Reverse Primer (5' → 3')          |
|-------------------------------|-------------------------------------|-----------------------------------|
| Bacitracin synthetase 3 (C454)| GAACGCTGATCGAGCATAAGA               | GAAGCAGAACGAGTGGAACA              |
| Asparagine synthase (C454)    | CCGTTCAAGCTCGGTTAAGA                | GAGGCTTGTTGTTGGCTTTC              |
| Isochorismate synthase (C454) | CGGAGCCGGAGAGATTTATG                | CGACAGGACCACCTTCTGA               |
| L-ectoine synthase (C454)     | ATTCTCGGCACGGAACAGGA                | AAACCAACTCCGTCCTTCTTG             |

### Workflow for Primer Testing

1. **DNA Extraction:**
   - DNA was extracted from *P. dendritiformis* C454 using two different kits to assess their extraction efficiency and the quality of the DNA;
   - kit1 - PureLink microbiom DNA Purification Kit invitrogen.
   - kit2 - NucleoSpin DNA stool
    

2. **PCR Amplification:**
   - Each primer pair was used to amplify its target gene from the extracted DNA.
   - The 16S rRNA gene was used as a control to validate the PCR process with primers previously tested successfully.

3. **Gel Electrophoresis:**
   - PCR products were analyzed using agarose gel electrophoresis to verify the presence and size of the amplified fragments.
   - Initial testing involved running a 1% agarose gel with 4 µL of SafeRed. This test included a negative control (no DNA) and a molecular weight ladder at different volumes (5 µL and 15 µL) to assess distribution accuracy.

   ![alt text](../images/trail%202%20-%20gel.png)

   - In the second gel electrophoresis run, DNA from Kit1 was tested for each target gene, along with negative controls and the 16S control. This run used a 2% agarose gel with 7 µL of SafeRed, and a molecular weight ladder at volumes of 20 µL and 25 µL. Additionally, a recA gene with validated primers was included to further validate the PCR setup.

   ![Gel Electrophoresis Run 2](../images/trail%202%20-%20gel.png)

The specificity and quality of the PCR amplification were confirmed by the presence of distinct bands corresponding to the expected product sizes, validating the designed primers for their intended use in qPCR analysis.

The negative control of 16S in the second run appears to show a band. This can be explained by human error. In any case, the primers seem to be normal, showing a single band each of the appropriate size of around 100bp.


# PCR Experiment Protocol

## Primers

The primers used in this experiment were designed in the laboratory ahead of time and were intended for several genes with different functions of *Paenibacillus dendritiformis* strain C454.

| Serial Number | Gene Name               | Gene Type     | Forward Primer             | Reverse Primer             |
|---------------|--------------------------|---------------|----------------------------|----------------------------|
| 1             | PgpA peptidase           | Protease      | TTACAAGTATGGCGGAACG        | GACTTCGAATTGGGAAGAGG       |
| 2             | M34 peptidases           | Protease      | ATTTGAAGGGCGTTACCC         | GTCCCTTCTTGCTGTAACC        |
| 3             | XcnG peptidase           | Protease      | GAAATCACATCAGGGTGATAGG     | GGTGCTCTATCATCGTTTGG       |
| 4             | S12 peptidases (MER0354617)           | Protease      | TGTTCCTGAACAGCTTGC         | CCGTAGGCGTAATACAAAGG       |
| 5             | C40 peptidases           | Protease      | GCTCTTGGTACAAAGTGAAGG      | CTTGTCAGATGGCTGATTGG       |
| 6             | Carboxypep DacF          | Protease      | AGCCTTGCTTCAATACCC         | CTTGTTCGTGTTGACAAGC        |
| 7             | S26A peptidases          | Protease      | CGATGTCGTCATTATCGATCC      | CGTTCACCTTCCAGCTTCC        |
| 8             | S12 peptidases (MER0353974)          | Protease      | GGACGTAATCGGCAATCC         | GAATAGACCCGCAAGTATCG       |
| 9             | S12 peptidases (MER1004512)           | Protease      | CTTCGTAACCGAGAAGAATGG      | GCCGCATAATGGGATATGG        |
| 10            | S12 peptidases  (MER0355148)         | Protease      | AGCTCGATACCGTTCCG          | TATAAGCCTCGTAGGACATACC     |
| 11            | dnaH                     | Physiological | GACAGAGCAAACTGATATAGGG     | ACTTCCATAGTGGACATTGC       |
| 12            | dnaA                     | Physiological | AATAAAGCCGTCGACTTCC        | TCCTGGGTTTGTTCTTTCC        |



## PCR Procedure

1. **Prepare the Reaction Mix:**
   - Prepare a master mix for each gene.
   - For each reaction, use the following components:
     - **Forward Primer:** 0.5 µL
     - **Reverse Primer:** 0.5 µL
     - **DNA Template:** 0.5 µL (except for the negative control, use UPW instead)
     - **PCR Master Mix:** 12.5 µL
     - **UPW:** 11 µL (or adjusted to final volume of 25 µL)

2. **Set up PCR Reaction:**
   - Label each tube appropriately for each gene.
   - Add the master mix to each tube.
   - Include negative controls (UPW instead of DNA) for each gene.
   - Include a positive control with Isochorismate synthase (C454):
     - **Forward Primer:** CGGAGCCGGAGAGATTTATG
     - **Reverse Primer:** CGACAGGACCACCTTCTGA

3. **PCR Cycling Conditions:**
   - Initial denaturation: 94°C for 1 minutes
   - 35 cycles of:
     - Denaturation: 94°C for 30 seconds
     - Annealing: 55°C for 30 seconds
     - Extension: 72°C for 30 minute
   - Final extension: 72°C for 5 minutes
   - Hold: 4°C

4. **Analyze PCR Products:**
   - Run the PCR products on an agarose gel (2%)
   - Visualize bands under UV light.

   ![alt text](../images/existing%20primers%20trial%201kb%20ladder%20(2%20microL)%202%25%20agar.png)

5. **Notes**

- As can be seen, the protease genes showed a single band.
- The scale did not really help in determining the size of the replicon,and it should be replaced in the following experiments.
- It should be noted that we know the size of the band of the positive control (slightly smaller than 100bp, the image is shown above - Isochorismate synthase), so we can compare it and determine with some certainty that the replicons are the correct size and we can continue working with these primers.
- The DNA replication genes (dnaA and dnaH) did not show good results, so I will repeat the experiment (perhaps already at the cDNA stage with primers I designed myself).




