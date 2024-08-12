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

   We did 12 tests. With each test of cDNA we also did a negative control (NC) in which the transition from RNA to cDNA was without the reverse transcriptase enzymes. Each of the two genes I chose, 16S and Asparagine synthase, we tested on two cDNA extracts of two bacteria: Paenibacillus dendritiformis and Bacillus subtilis, and at three different concentrations 4,2,0.5 µg/µL.
## Results
![results](../images/cdna%20validation.png)

As you can see the transition to cDNA was successful and you can see it in the upper part of the gel, of the 16S. In the asparagine synthase gene no bands were seen at all, which indicates that the gene was not expressed during RNA extraction. This makes sense since Asparagine synthase catalyzes the biosynthesis of the amino acid asparagine. In LB medium, cells can efficiently obtain asparagine from the medium.

It can also be seen that the most prominent band in the concentration of 4  µg/µL. **However**, it is important to note that we ran 28 cycles in PCR and not 35. A concentration of 4 µg/µL is **very high** for a qPCR process, so it may be better to work with the concentration of 0.5 µg/µL.
