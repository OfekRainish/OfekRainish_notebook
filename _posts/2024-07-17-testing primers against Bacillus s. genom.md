

# Testing the primers against the genome of Bacillus subtilis

These are the primers of Benibetselos that have undergone a successful PCR test. Since I want to proceed with qPCR and part of the test is to examine Peanibacillus together with Bacillus subtilis, it is necessary to make sure that the primers do not amplify any segment in the genome of the bacillus. For this I used primerBLAST and ran these primers against the genome of NCIB 3610 (taxid:1423), the strain I work with.

| Serial Number | Gene Name                       | Gene Type     | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Tested Against a B.s Genome |
|---------------|----------------------------------|---------------|--------------------------|--------------------------|-----------------------------|
| 1             | PgpA peptidase                  | Protease      | TTACAAGTATGGCGGAACG      | GACTTCGAATTGGGAAGAGG      | pass                        |
| 2             | M34 peptidases                  | Protease      | ATTTGAAGGGCGTTACCC       | GTCCCTTCTTGCTGTAACC       | pass                        |
| 3             | XcnG peptidase                  | Protease      | GAAATCACATCAGGGTGATAGG   | GGTGCTCTATCATCGTTTGG      | pass                       |
| 4             | S12 peptidases (MER0354617)     | Protease      | TGTTCCTGAACAGCTTGC       | CCGTAGGCGTAATACAAAGG      | Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 5             | C40 peptidases                  | Protease      | GCTCTTGGTACAAAGTGAAGG    | CTTGTCAGATGGCTGATTGG      | Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 6             | Carboxypep DacF                 | Protease      | AGCCTTGCTTCAATCCC        | CTTGTTCGTGTTGACAAGC       | Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 7             | S26A peptidases                 | Protease      | CGATGTCGTCATTATCGATCC    | CGTTCACCTTCCAGCTTCC       | Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 8             | S12 peptidases (MER0353974)     | Protease      | GGACGTAATCGGGCAATCC      | GAATAGACCCGCAAGTATCG      | pass                        |
| 9             | S12 peptidases (MER1004512)     | Protease      | CTTCGTAACCGAGAAGAATGG    | GCCGCATAATGGGGATATGG      | pass                        |
| 10            | S12 peptidases (MER0355148)     | Protease      | AGCTCGATACCGTTCCG        | TATAAGCCTCGTAGGACATACC    | Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 11            | Bacitracin synthetase 3 (C454)  | Natural Products | GAACGCTGATCGAGCATAAGA   | GAAGCAGAACGAGTGGAACA     | pass                        |
| 12            | Asparagine synthase (C454)      | Natural Products | CCGTTCAAGCTCGGTTAAGA    | GAGGCTTGTTGTTGGCTTTC     | pass                        |
| 13            | Isochorismate synthase (C454)   | Natural Products | CGGAGCCGGAGAGATTTATG    | CGACAGGACCACCTTCTGA      |  Bacillus subtilis subsp. subtilis str. 168 complete genome                        |
| 14            | L-ectoine synthase (C454)       | Natural Products | ATTCTCGGCACGGAACAGGA    | AAACCAACTCCGTCCTTCTTG    | pass                        |

Note that even primers that had matches, in most of them the replicons are very long and therefore it can be resolved by the elongation time. None of the matches are perfect with matches of five or six nucleotides.