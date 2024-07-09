# Finding Sequences for Primer Design for qPCR

## Part 1 - Natural Products (4 sequences)
- Obtain BGCs sequences using antiSMASH software on Paenibacillus dendritiformis.
- Choose the core gene from each BGC for primer design.

### genes selected:
      1.Asparagine Synthase.

      2.Bacitracin Synthetase 3.

      3.Isochorismate Synthase  

      4.L-Ectoine Synthase


---

## Part 2 – "Physiological Genes" (13 sequences)

### Acquisition Process
1. Search for interesting genes (e.g., cell division, cell movement, DNA division, etc.) using ChatGPT for known genes in Paenibacillus dendritiformis (strain C454).
2. Search for the desired gene in NCBI: "Paenibacillus dendritiformis + gene name".
3. Retrieve the sequence in FASTA format and verify by running BLAST against 'whole genome shotgun contigs (wgs)' database with Paenibacillus dendritiformis str. C454 as the organism.
4. Ensure sequence accuracy by comparing with GenBank and using alignment tools like BioEdit.

- List of Genes:
  1. dnaA
  2. DnaG
  3. DnaN
  4. gyrA
  5. GyrB
  6. ftsA
  7. ftsZ
  8. ftsI
  9. motA
  10. motB
  11. minC
  12. minD
  13. mreB

---

## Part 3A – Proteases from Article on Secretome of Paenibacillus (5 sequences)

### Methodology
1. Extract protein sequences of proteases from the article **"Secretome of Paenibacillus sp. S 12 provides an insight about its survival and possible pathogenicity"**.
2. Convert protein sequences to amino acid sequences.
3. Blast amino acid sequences in UniProt with strain C454.
4. Verify translation accuracy using ExPASy Translate Tool.
5. Confirm sequence fidelity by comparing with original protein sequence using alignment tools like BioEdit.

- List of Proteases:
  1. K5ABD2 COG0793 Periplasmic protease 
  2. S9SXA2 ATP-dependent protease ATP-binding subunit ClpX
  3. A0A383R5G4 - ATP-dependent Clp endopeptidase proteolytic subunit ClpP
  4. A0A7Y4IB32 - ATP-dependent zinc metalloprotease FtsH
  5. A0A7Y4I7I5 - Trypsin-like serine protease 

---

## Part 3B – Proteases from MEROPS (11 sequences)

### Procedure
1. Visit MEROPS site and select 'Paenibacillus dendritiformis' as the organism.
2. Retrieve protease sequences and convert them to FASTA format using ChatGPT.
3. BLAST protein sequences against Paenibacillus dendritiformis in NCBI to obtain sequences in FASTA format.
4. Use SignalP 6.0 to predict secretion pathways (Sec/SPI or LIPO/SPII).
5. Confirm secretion type by comparing with UniProt and verifying translation accuracy with ExPASy Translate Tool.

- List of Proteases:
  1. MER1004512 - Family S12 unassigned peptidases [S12.UPW] (Paenibacillus dendritiformis)
  2. MER0354617 - Family S12 unassigned peptidases [S12.UPW] (Paenibacillus dendritiformis)
  3. MER0355148 - Family S12 unassigned peptidases [S12.UPW] (Paenibacillus dendritiformis)
  4. MER0297385 - D-Ala-D-Ala carboxypeptidase DacF [S11.005] (Paenibacillus dendritiformis)
  5. MER0353802 - XcnG peptidase ({Xenorhabdus nematophila}) [S12.011] (Paenibacillus dendritiformis)
  6. MER0452375 - PgpA peptidase [C40.006] (Paenibacillus dendritiformis)
  7. MER0454194 - Family C40 unassigned peptidases [C40.UPW] (Paenibacillus dendritiformis)
  8. MER0295257 - Family M34 unassigned peptidases [M34.UPW] (Paenibacillus dendritiformis)
  9. MER0353974 - Family S12 unassigned peptidases [S12.UPW] (Paenibacillus dendritiformis)
  10. MER0373982 - Subfamily S26A unassigned peptidases [S26.UPA] (Paenibacillus dendritiformis)
  11. M15 family metallopeptidase [Paenibacillus dendritiformis]
