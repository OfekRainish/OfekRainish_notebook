# Primer selection for qPCR Experiment

## Objective
To investigate the effect of light and radiation on Thiamine production in the red alga *Gracilaria gracilis*.

## Reference Gene
**Glyceraldehyde-3-phosphate dehydrogenase (GAPDH)** is selected as the reference gene due to its consistent expression under various conditions. GAPDH has been successfully used in previous studies as a stable reference gene [1], making it a reliable choice for normalization in qPCR experiments.

### Background on GAPDH
GAPDH is a key enzyme in the glycolytic pathway, facilitating the conversion of glyceraldehyde-3-phosphate to 1,3-bisphosphoglycerate. Due to its essential role in basic cellular metabolism, its expression levels remain relatively constant across different conditions. 

However, since Glyceraldehyde-3-phosphate dehydrogenase (GAPDH) is involved in sugar metabolism, its expression might be influenced by varying light levels due to the plant's photosynthetic activity and subsequent sugar production. Therefore, it may not be the most suitable reference gene for a qPCR experiment studying exposure to different light levels. Since no additional reference genes were found in the literature and genes that are usually used as a reference gene in algae (tubulin and actin) were not found in Gracilaria gracilis, this gene is worth checking.

## Test Gene
I was not able to find genes related to the biosynthesis of carotenoids in Gracilaria gracilis, so I decided to test another gene that I think will also be affected by exposing the organism to different intensities of light and radiation.

**Thiamine biosynthesis protein G (thiG)** is selected for this experiment. Thiamine (vitamin B1) is crucial for cellular functions and stress responses. ThiG is involved in the biosynthesis of thiamine, which plays an important role in the stress response of organisms. 

### Background on thiG
Thiamine is essential for metabolic processes and serves as a cofactor for several enzymes. It has been shown to aid in oxidative stress resistance in various organisms [2]. Studies have indicated a positive correlation between thiamine levels and resistance to oxidative stress, including in the context of alcohol consumption in humans [3] and different types of stress in phytoplankton [4]. Given that both high light intensities and radiation can cause stress and oxidative damage [5, 6], it is hypothesized that thiG expression will change under these conditions.

## Primers
The following primers will be used for the qPCR experiment:

| Primer name     | Primer sequence         | Sequence length | GC content (%) | Molecular weight (Daltons) | nmol/A260 | micrograms/A260 | Basic Tm (°C) | Salt adjusted Tm (°C) | Nearest neighbor Tm (°C) |
|-----------------|-------------------------|-----------------|----------------|----------------------------|-----------|------------------|---------------|-----------------------|--------------------------|
| fwd GAPDH (reference)       | ATGGCTATCCCGAAAGTTGG    | 20              | 50.00          | 6157.07                    | 5.10      | 31.41            | 52            | 47                    | 62.70                    |
| rev GAPDH (reference)       | GACCGTGCGTGGAGTCATAC    | 20              | 60.00          | 6158.05                    | 5.15      | 31.71            | 56            | 51                    | 65.66                    |
| fwd thiG (tested)     | GGTTCTGGCCAAGGTTTACA    | 20              | 50.00          | 6148.06                    | 5.22      | 32.07            | 52            | 47                    | 63.07                    |
| rev thiG (tested)       | CAGATGCGCCATTTCCATA     | 19              | 47.37          | 5747.80                    | 5.55      | 31.90            | 49            | 44                    | 60.52                    |







## References
1. Chang, L., Sui, Z., Fu, F. et al. Relationship between gene expression of UDP-glucose pyrophosphorylase and agar yield in *Gracilariopsis lemaneiformis* (Rhodophyta). J Appl Phycol 26, 2435–2441 (2014). https://doi.org/10.1007/s10811-014-0277-7
2. Kartal, B., & Palabiyik, B. (2019). Thiamine leads to oxidative stress resistance via regulation of the glucose metabolism. Cellular and Molecular Biology, 65(1), 73-77
3. de Carvalho Gonçalves, Á., Soldi, L. R., & Portari, G. V. (2020). Thiamine, oxidative stress, and ethanol. In Molecular Nutrition (pp. 207-223). Academic Press.
4. Sylvander, P., Häubner, N. & Snoeijs, P. The Thiamine Content of Phytoplankton Cells Is Affected by Abiotic Stress and Growth Rate. Microb Ecol 65, 566–577 (2013). https://doi.org/10.1007/s00248-012-0156-1
5. Xue, S., Zang, Y., Chen, J., Shang, S., Gao, L., & Tang, X. (2022). Ultraviolet-B radiation stress triggers reactive oxygen species and regulates the antioxidant defense and photosynthesis systems of intertidal red algae *Neoporphyra haitanensis*. Frontiers in Marine Science, 9, 1043462
6. Coulombier, N., Nicolau, E., Le Déan, L., Antheaume, C., Jauffrais, T., & Lebouvier, N. (2020). Impact of light intensity on antioxidant activity of tropical microalgae. Marine drugs, 18(2), 122
