# Protocol for cDNA Examination by PCR

## Materials
- cDNA of **Paenibacillus dendritiformis** and **Bacillus subtilis** (grown in liquid medium without special treatment)
- RNA used for making cDNA (negative control without reverse transcriptase enzymes)
- UPW (Ultra Pure Water)
- PCR reagents
- Primers (listed in the table below)

## cDNA and RNA Concentrations
Both cDNA and RNA extracts were at a concentration of 20 µg/µL.

## Dilution for qPCR
we tried 3 different concentrations: 4,2,0.5 µg/µL.
- Use the formula  C1*V1 = C2*V2 (where C2= C1/dilution factor)

  In our case, we did a series of dilutions with the amounts we decided on, but any other way that reaches these concentrations is fine:
  - for c=4 µg/µL: 6 µL cDNA extraction (of 20 µg/µL) + 24 µL UPW.
  - for c=2 µg/µL: 12 µL cDNA extraction (of 4 µg/µL) + 12 µL UPW.
  - for c=0.5 µg/µL: 6 µL cDNA extraction (of 2 µg/µL) + 18 µL UPW.

  

## Genes and Primers

| Gene Name           | Type             | Forward Primer (F)        | Reverse Primer (R)          |
|---------------------|------------------|---------------------------|-----------------------------|
| Asparagine synthase | Natural Products | CCGTTCAAGCTCGGTTAAGA       | GAGGCTTGTTGTTGGCTTTC         |
| 16s RNA             | Control          | ACACTGACGACATGGTTCTACAGTGYCAGCMGCCGCGGTAA       | TACGGTAGCAGAGACTTGGTCTCCGYCAATTYMTTTRAGTTT-        |


## PCR Reaction
1. **Prepare PCR Mix:**
   - Mix the cDNA or RNA extracts with the respective primers and PCR reagents.

2. **PCR Tests:**


## Results
1.   in the first image We did 12 tests. With each test of cDNA we also did a negative control (NC) in which the transition from RNA to cDNA was without the reverse transcriptase enzymes. Each of the two genes I chose, 16S and Asparagine synthase, we tested on two cDNA extracts of two bacteria: Paenibacillus dendritiformis and Bacillus subtilis, and at three different concentrations 4,2,0.5 µg/µL.

      ![results](../images/cdna%20validation.png)

      As you can see the transition to cDNA was successful and you can see it in the upper part of the gel, of the 16S. In the asparagine synthase gene no bands were seen at all, which indicates that the gene was not expressed during RNA extraction. This makes sense since Asparagine synthase catalyzes the biosynthesis of the amino acid asparagine. In LB medium, cells can efficiently obtain asparagine from the medium.

      It can also be seen that the most prominent band in the concentration of 4  µg/µL. **However**, it is important to note that we ran 28 cycles in PCR and not 35. A concentration of 4 µg/µL is **very high** for a qPCR process, so it may be better to work with the concentration of 0.5 µg/µL.



2. In the first image We did 24 tests. For each of the four genes recA, 16S, Carboxypep DacF
and XcnG peptidase We performed 6 tests:

   1. Versus genomic DNA (DNA)
   2. Against cDNA of Pd. (cDNA Pd)
   3. Against RNA Pd in ​​which no reverse transcriptase enzyme was inserted (RNA control Pd)
   4. against cDNA of Bs. (cDNA Bs)
   5. Against RNA Bs where no reverse transcriptase enzyme was inserted. (RNA control Bs)
   6. water control (W.C)

| Gene Name         | Type                     | Forward Primer (F)                          | Reverse Primer (R)                          |
|-------------------|--------------------------|---------------------------------------------|---------------------------------------------|
| recA              | Potential normalizing gene| GAG ATT GAA GGC GAG ATG G                   | GTC TTG GAC TTG CTG ATC G                   |
| 16S               | Control                  | ACACTGACGACATGGTTCTACAGTGGYCAGCMGCCGCGGTAA- | TACGGTAGCAGAGACTTGGTCTCCGYCAATTYMTTTRAGTTT- |
| carboxypep DacF   | Protease                 | AGC CTT GCT TCA ATA CCC                     | CTT GTT CGT GTT GAC AAG C                   |
| XcnG peptidase    | Protease                 | GAA ATC ACA TCA GGG TGA TAG G               | GGT GCT CTA TCA TCG TTT GG                  |




   In each row there are two PCRBIO ladder 3 (pcr biosystems) ladders, once diluted 5 times and once diluted 10 times, in this order.

   The gel is 2% agar and 4µL safeRed, and runs for 50 minutes at 120 volts.

PCR Program

**Initial Denaturation:**
- 94°C for 1 minute

**35 Cycles:**
- Denaturation: 94°C for 30 seconds
- Annealing: 53°C for 30 seconds
- Extension: 72°C for 30 seconds

**Final Extension:**
- 72°C for 5 minutes


   ![results](../images/cdna%202.png)

   3. For each of the four genes recA, 16S, Carboxypep DacF
and gyrB we performed 5 tests:

   1. Versus genomic DNA of P.d (PdG).
   2. Against cDNA of P.d (Pd+).
   3. Against RNA Pd in ​​which no reverse transcriptase enzyme was inserted (Pd-).
   4. against genomic DNA of Bs (BsG).
   5. water control (W.C).

**Notice we used different primers for RecA and Carboxypep DacF**

| Gene Name          | Gene Type                  | Forward Primer                                               | Reverse Primer                                               |
|--------------------|----------------------------|--------------------------------------------------------------|--------------------------------------------------------------|
| 16S                | Control                    | ACACTGACGACATGGTTCTACAGTGGYCAGCMGCCGCGGTAA                    | TACGGTAGCAGAGACTTGGTCTCCGYCAATTYMTTTRAGTTT                    |
| gyrB               | Potential normalizing gene | CAACTTGAGCGGGGATGATG                                          | GAACAAGGACTCGACGATGC                                          |
| RecA               | Potential normalizing gene | GGAATCTCCCATCTCGCCTT                                          | ATCGACGAACTGCTTCTGTC                                          |
| carboxypep DacF    | Protease                   | CCCACTTCGAGAACAGCAAC                                          | TGGCGCAGGTAATCTTCGTA                                          |



   In each row there are two PCRBIO ladder 3 (pcr biosystems) ladders, once diluted 5 times and once diluted 10 times, in this order.

   The gel is 2% agar and 4µL safeRed, and runs for 50 minutes at 120 volts.

PCR Program

**Initial Denaturation:**
- 94°C for 1 minute

**35 Cycles:**
- Denaturation: 94°C for 30 seconds
- Annealing: 53°C for 30 seconds
- Extension: 72°C for 30 seconds

**Final Extension:**
- 72°C for 5 minutes

![](../images/gel/16s%20gyrB%20dacF.png)

As observed, there are smears and unexpected bands appearing in locations other than anticipated. After evaluation, it was decided to repeat the PCR reaction at a higher annealing temperature of 63°C instead of 53°C. The second run focused on the recA, gyrB, and DacF genes.

![](../images/gel/recA%20gyrB%20DacF%207.10.png)