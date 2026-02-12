# Guide: Annotating Metabolic Pathways Using blastKOALA and KEGG Mapper (Color)

This guide explains how to obtain **KO terms** for your organism using **blastKOALA**, prepare the data for **KEGG Mapper (Color)**, and visualize metabolic pathways with expression-based color coding.

---

## 1. Obtaining KO Terms Using blastKOALA

1. Go to the [blastKOALA](https://www.kegg.jp/blastkoala/) website:  
   
2. Paste your **amino acid sequences** into the input window.

3. Under *taxonomy group*, select **Prokaryotes**.

4. Enter your **university email address**.

5. Submit the job.

![](../images/koala/koalablast.png)

6. You will receive an email from KEGG.  
   Click **Submit** to confirm.

7. When the analysis is complete, you will receive a second email.  
   Click the link to access your results.

8. Download the results file and open it in Excel.

**Important:**  
Send **all genes** from your genome/proteome, *not just DEGs*.  
KEGG Mapper requires a complete set of KO annotations to properly mark pathways.

---

## 2. Preparing the Input Table for KEGG Mapper (Color)

1. Open the blastKOALA results file in Excel.

2. Prepare a table with **three columns**:

   - **Column 1:** Gene name  
   - **Column 2:** KO term (format: `Kxxxxx`)  
   - **Column 3:** Color code

3. Choose colors based on your gene expression results.  
   For example, for the 20-hour comparison:

   - **Red** – Upregulated  
   - **Blue** – Downregulated  
   - **Yellow** – Genes that have a KO term but are not differentially expressed

![](../images/koala/excelTable.png)

4. Save the file or keep it open for copying.

---

## 3. Running KEGG Mapper (Color)

1. Search for **"KEGG Mapper"** in Google.

2. Click [KEGG Mapper – Color](https://www.genome.jp/kegg/mapper/color.html):  


3. Copy **Columns 2 and 3** (KO term + color) from your Excel file.

4. Paste them into the input box in KEGG Mapper (Color).

5. Check the box:  
   **"Use uncolored diagrams"**

6. Click **EXEC**.

![](../images/koala/keggMapper.png)

7. KEGG Mapper will generate a list of metabolic pathways.  
   Each pathway diagram will display your selected **color codes** directly on the enzymes.

---

## 4. Interpreting the Pathways

- Look for **clusters**, not isolated colored enzymes.  
- Patterns such as:
  - Entire pathways turning blue or red  
  - Sequential genes in a pathway showing coordinated expression  
  - A group of connected reactions with consistent color trends  

These patterns usually indicate a biologically meaningful response.

*(Place for pictures)*

---

## Notes and Tips

- **blastKOALA assigns only one KO per gene.**  
  A gene may participate in multiple pathways, but KOALA assigns the *primary* function.  
  For expanded mapping, consider tools like **Galaxy** or **eggNOG mapper**.

- **This method is not KEGG enrichment.**  
  We are not calculating statistical enrichment.  
  Instead, we are overlaying expression data onto KEGG maps to visually identify:
  - Upregulated modules  
  - Downregulated modules  
  - Pathway-level expression trends

This visual inspection can reveal important biological insights even without enrichment statistics.

---

If you'd like, I can also create a **flowchart**, **image placeholders**, or **example Excel templates** to accompany the guide.

## Resalts 
* reminder - these results are for the 20hr time point
### ABC tansporters
![](../images/koala/pathways/ABC%20transporters%20-%20Reference%20pathway.png)

A clear decrease in proteins that make up the iron and oligopeptide transporter can be seen, alongside an increase in proteins that make up the ABC transporters of alduronate and glucose/mannose.

possibly, Surfactin pushes *P. dendritiformis* to temporarily shut down nutrient systems it doesn’t need (iron, peptides) and activate transporters that help it import energy-rich sugars and environmental carbohydrates. it is unlikely that the peptidic product of surfactin is transported through Opp sinse its bulky and branched.

### Argenine biosynthis
![](../images/koala/pathways/Arginine%20biosynthesis%20-%20Reference%20pathway.png)

We see an upregulation of genes that constetute the arganine biosynthesis process. seems that argininosuccinate accumulates rather than being converted to arginine/fumarate.

### Chemotaxix & Flagellar assembly
![](../images/koala/pathways/Bacterial%20chemotaxis%20-%20Reference%20pathway.png)
![](../images/koala/pathways/Flagellar%20assembly%20-%20Reference%20pathway.png)

motA & motB Encode parts of the flagellar motor stator, FliC Encodes flagellin, the main structural protein of the flagellar filament, and MCP (Methyl-accepting chemotaxis protein) Acts as a sensor that detects chemical gradients in the environment.

### L histadine --> L glutmate slowdown ?
![](../images/koala/pathways/Histidine%20metabolism%20-%20Reference%20pathway.png)
This pathway is the main route for converting histidine into glutamate.

### Increse in pyruvat metabolism
![](../images/koala/pathways/Lipoic%20acid%20metabolism%20-%20Reference%20pathway.png)
Pd is increasing the flux from pyruvate to acetyl-CoA, potentially enhancing energy production.

### Peptidoglican biosynthasis
![](../images/koala/pathways/Peptidoglycan%20biosynthesis%20-%20Reference%20pathway.png)
The downregulation of these DD-carboxypeptidases suggests the bacterium is modulating cell wall remodeling

### Quorum sensing
![](../images/koala/pathways/Quorum%20sensing%20-%20Reference%20pathway.png)
A decrease in expression of genes encoding oligopeptide transport system substrate-binding proteins, which bind a variety of autoinducers, can be seen. It is possible that Pd reduces its sensitivity to reading QS signals and begins to act according to the signals it receives. the downstream effect might lead to a decrease amount of biofilm (can add a biofilm assay?)


