# genes of interest

The overall goal of this research is to investigate how surfactin affects gene expression in the bacterium *Paenibacillus dendritiformis* (Pd).

We had two main points of interest beforehand:

1. **Effect on the production of natural products**
2. **Identification of proteases that break down surfactin**, which we believe are secreted outside the cell

However, we also aimed to explore additional trends. To help focus on the "genes of interest" for follow-up work, we performed **RNA-seq analysis**, comparing Pd samples exposed to surfactin with untreated controls at two time points:

- **20 hours (logarithmic growth phase)**
- **44 hours (stationary phase)**

Following RNA-seq analysis, we obtained hundreds of genes and narrowed the list down to a subset of "genes of interest." Here are our main observations:

---

## 1. Stress Response
- There appears to be a stress response to surfactin, which makes the use of surfactin degradation products even more intriguing, as the bacterium seems willing to “sacrifice” itself to utilize surfactin.  
- Genes like **RNA polymerase sigma factor SigY** and **carbon starvation protein CstA** increase in expression.  
- Interestingly, we did **not** see classic stress responses (e.g., SOS response). Some stress-related genes actually **decrease** (e.g., **ECF subfamily RNA polymerase sigma24 subunit**).  
- This suggests the stress is likely **targeted at the cell membrane/envelope**, not DNA damage.


** In order to ensure that no classic SOS response genes were significantly activated,  
I compiled a list of canonical SOS genes and matched them to their corresponding names in *P. dendritiformis*:

| Gene ID | Annotation |
|----------|-------------|
| PDENDC454_09755 | radA |
| PDENDC454_11340 | recN |
| PDENDC454_11345 | recN |
| PDENDC454_19328 | recA |
| PDENDC454_04681 | SOS response |
| PDENDC454_04686 | SOS response |
| PDENDC454_26503 | ruvB |
| PDENDC454_26493 | revC |
| PDENDC454_26498 | revC |
| PDENDC454_22519 | DinB-related |
| PDENDC454_22524 | DinB-related |
| PDENDC454_22539 | DinB-related |
| PDENDC454_20987 | DinB-related |
| PDENDC454_20992 | DinB-related |
| PDENDC454_02050 | DinB-related |
| PDENDC454_02055 | DinB-related |
| PDENDC454_05381 | recF pathway |
| PDENDC454_05386 | recF pathway |
| PDENDC454_03929 | – |
| PDENDC454_03934 | – |

Apart from one gene — **PDENDC454_05386**, which was upregulated after 44 hours of surfactin exposure (*padj = 0.03*) —  
none of these were classified as differentially expressed genes (DEGs).



#### Python Code Used to Compare the SOS Gene List with RNA-seq Results

```python
import pandas as pd

# Load the two Excel files
# Replace the paths with your actual filenames
df1 = pd.read_excel(
    "C:/Users/USER/OneDrive - University of Haifa/Documents/HAIFA/research/data analyzing/sos_genes/sos_genes.xlsx"
)
df2 = pd.read_excel(
    "C:/Users/USER/OneDrive - University of Haifa/Documents/HAIFA/research/data analyzing/mapping/deseq2/output/treatment - downup regulation/signalp0.6 res/full_table_with_secretion_predictions_anti17.xlsx"
)

# Extract the relevant columns
genes1 = df1["Gene_ID"].dropna().astype(str).str.strip()
genes2 = df2["GeneID_clean"].dropna().astype(str).str.strip()

# Find the intersection
common_genes = set(genes1).intersection(set(genes2))

# Print results
if common_genes:
    print("Values that appear in both tables:")
    for gene in common_genes:
        print(gene)
else:
    print("No common values found.")

```

## 2. ABC Transporters & Efflux Systems
- These genes encode proteins that transport substances into and out of the cell.  
- We see both **increases and decreases** in their expression. This is common as the cell adapts to a new environment, balancing import and export needs.  
- **ABC transporters** can import nutrients (amino acids, sugars, ions) or export toxic compounds, antibiotics, and signals.  
- **Efflux pumps** mainly export harmful substances, often contributing to **antibiotic resistance** and detoxification systems.

---

## 3. Prophage Genes
- Increased expression of **prophage genes** suggests a phage may be transitioning from a **lysogenic** to a **lytic** state.  
- This transition likely results from stress, meaning the **cell membrane may be under attack from both directions**:  
  - **Outside:** surfactin  
  - **Inside:** phage toxins  
- Alternatively, the internal attack may occur only after surfactin is degraded.

*Implication:* Genes related to **membrane stress responses** and **cell wall remodeling** are worth deeper exploration.

---

## 4. Natural Product Biosynthetic Gene Clusters (BGCs)
- Several BGCs identified by **antiSMASH** show **decreased expression** in response to surfactin.  
- We will focus on their **core genes**.  
- This decrease fits the hypothesis that surfactin is a **public good** used to reduce metabolic costs.  
- Question: Is this cost-saving alone sufficient to justify the stress response in Pd?

---

## 5. Proteases
- Several genes encoding **secreted proteases** (using signalP 6.0) show increased expression with surfactin exposure.  
- These are prime candidates for being involved in **surfactin degradation**.

---

## 6. Chemotaxis and Motility
- We observe increased expression of **chemotaxis-related genes**, consistent with the phenotype where Pd is **attracted to and moves toward surfactin**.

---

## 7. Antimicrobial Production
- Genes involved in **antimicrobial substance production** decrease in expression.  
- Possible explanations:  
  - pd utilizes surfactin to produce antimicrobial natural products and "skips" synthesis steps. Therefore, we see a decrease in the expression of key enzymes in the synthesis of these substances. This supports the hypothesis that pd uses surfactin from the environment to save the metabolic cost of producing natural products. 
  - Pd may avoid harming *Bacillus subtilis*, the surfactin producer, since Pd benefits from surfactin.  
- Indeed, Pd is attracted to and surrounds *B. subtilis* colonies, but *B. subtilis* does not die as a result of this interaction.

---

## Bottom Line
Even though Pd can degrade surfactin and shows a clear attraction to it, the bacterium still **suffers stress from its presence**.  
This raises the key question:

**What justifies this sacrifice?**  
- Is it simply the metabolic cost of producing natural products?  
- Or does surfactin provide additional benefits that outweigh the stress?

---

## Next Steps
- Further reduce the list of candidate genes.  
- Validate their expression in real time using **qPCR**.  
- The reduced gene list is attached [here](../exel%20files/deseq2/genes_of%20interest1.xlsx).
