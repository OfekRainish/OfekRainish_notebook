# Identifying Genes with Specific Functional Roles from RNA-seq Data

This document describes the workflow used to identify genes associated with a specific biological function (for example **membrane activity**) from RNA-seq differential expression data.
In addition to gene annotations provided in the GTF file (NCBI), two complementary approaches were used to locate candidate genes among the differentially expressed genes (DEGs):

1. **Keyword-based search using InterPro annotations**
2. **Gene Ontology (GO) term mapping**

Because annotation coverage is often incomplete, combining multiple approaches increases the likelihood of identifying relevant genes.

---

# 1. Keyword-Based Search Using InterPro

The **InterPro** database predicts protein family membership, domains, and functional activity based on amino acid sequence. All DEG protein sequences were submitted to the InterPro web server. The resulting output file can be downloaded in **TSV format**.

The InterPro results were then searched for functional keywords using Python.

Example workflow:

```python
import pandas as pd

# Define input file path
input_file = "C:/path/to/DEGs_interpro_output.tsv"

# Load the TSV file (tab-separated values)
df = pd.read_csv(input_file, sep='\t', header=None, comment='#')

# InterPro format: annotation information is typically found in columns 11–12
annotation_columns = [11, 12]

# Define keywords related to the activity of interest
keywords = ['protease', 'hydrolase', 'peptidase', 'proteinase', 'lipase', 'amidase']

# Search for keywords in annotation columns
matches = df[df[annotation_columns].apply(
    lambda row: any(any(k in str(cell).lower() for k in keywords) for cell in row),
    axis=1
)]

# Extract unique gene IDs (usually column 0)
gene_ids = sorted(matches[0].unique())

# Save gene list
with open("cutting_genes.txt", "w") as f:
    for gene in gene_ids:
        f.write(gene + "\n")

# Print results
print("Genes with cutting activity:")
for gene in gene_ids:
    print(gene)
```

This script scans the InterPro annotations for specified keywords and extracts genes whose predicted protein functions match the search terms.

Note that **InterPro cannot assign predictions to every gene**, therefore this method alone may miss relevant candidates.

---

# 2. GO Term–Based Identification

A second strategy uses **Gene Ontology (GO)** annotations to identify genes involved in a specific biological process or molecular function.

This approach consists of **two steps**.

---

## Step A – Identify Genes in the Organism with Relevant GO Terms

First, genes in the organism of interest (*Paenibacillus dendritiformis*) that contain GO terms associated with the desired function were identified.

The script below uses:

* an **Excel file containing all genes with GO annotations**
* the **GO ontology file (`go-basic.obo`)**, which describes relationships between GO terms

The code expands selected **root GO terms** to include all their descendant terms in the GO hierarchy, and then identifies genes annotated with those terms.

```python
import pandas as pd

# ---------- FILE PATHS ----------
excel_file = "C:/Users/goTerms.xlsx"
obo_file = "C:/Users/go-basic.obo"

# ---------- LOAD EXCEL ----------
df = pd.read_excel(excel_file)

# ---------- DEFINE ROOT GO TERMS ----------
parent_go_terms = {
    "GO:0061024": "membrane organization",
    "GO:0099500": "vesicle fusion to plasma membrane",
    "GO:0110165": "cellular anatomical structure",
    "GO:0006886": "intracellular protein transport",
    "GO:0007009": "plasma membrane organization",
    "GO:0071711": "basement membrane organization"
}

# ---------- PARSE OBO FILE ----------
go_children = {}
current_term = None
parents = []

with open(obo_file, "r") as f:
    for line in f:
        line = line.strip()

        if line == "[Term]":
            current_term = None
            parents = []

        elif line.startswith("id: GO:"):
            current_term = line.split("id: ")[1]

        elif line.startswith("is_a: GO:"):
            parent = line.split("is_a: ")[1].split(" !")[0]
            parents.append(parent)

        elif line.startswith("relationship: part_of GO:"):
            parent = line.split("relationship: part_of ")[1].split(" !")[0]
            parents.append(parent)

        elif line == "" and current_term:
            for p in parents:
                go_children.setdefault(p, set()).add(current_term)

# ---------- FUNCTION TO GET ALL DESCENDANTS ----------
def get_all_children(term, children_dict):
    result = set()
    stack = [term]

    while stack:
        t = stack.pop()

        for child in children_dict.get(t, []):
            if child not in result:
                result.add(child)
                stack.append(child)

    return result

# ---------- BUILD ROOT → DESCENDANT MAP ----------
root_to_terms = {}

for root in parent_go_terms:
    expanded = get_all_children(root, go_children)
    expanded.add(root)
    root_to_terms[root] = expanded

# ---------- MATCH GENES ----------
results = []

for _, row in df.iterrows():

    gene = row["Gene ID"]
    go_terms = str(row["GO Term"]).split(";")
    go_terms = [term.strip() for term in go_terms]

    for root, term_set in root_to_terms.items():
        for go_term in go_terms:

            if go_term in term_set:
                results.append({
                    "Gene ID": gene,
                    "GO Term": go_term,
                    "Root GO": root,
                    "Root Description": parent_go_terms[root]
                })

# Convert to DataFrame
results_df = pd.DataFrame(results)

# Remove duplicates
results_df = results_df.drop_duplicates()

# Save detailed table
results_df.to_excel("envelope_related_genes_detailed.xlsx", index=False)

# Save unique gene list
unique_genes = sorted(results_df["Gene ID"].unique())

with open("envelope_related_genes.txt", "w") as f:
    for gene in unique_genes:
        f.write(gene + "\n")

print(f"Total unique envelope-related genes: {len(unique_genes)}")
print("Results saved to:")
print(" - envelope_related_genes_detailed.xlsx")
print(" - envelope_related_genes.txt")

print(results_df["Root Description"].value_counts())
```

The output includes:

* A **detailed table** linking genes to GO terms and root categories
* A **unique gene list** of all genes associated with the function of interest

---

## Step B – Identify Which of These Genes Are Differentially Expressed

The next step is to determine which of the function-related genes identified in Step A are also present in the **DEG list**.

This is done by matching gene IDs between the DEG table and the GO-filtered gene list.

```python
import pandas as pd

# ==============================
# Input files
# ==============================

degs_file = "C:/Users/path/to/all_degs.xlsx"
go_file = "C:/Users/path/to/membrane_related_genes_detailed.xlsx"

# Load DEG list
degs_df = pd.read_excel(degs_file)
deg_genes = set(degs_df["GeneID_clean"].unique())

# Load GO results
go_df = pd.read_excel(go_file)

# ==============================
# Filter GO table for DEGs
# ==============================

filtered_df = go_df[go_df["Gene ID"].isin(deg_genes)].copy()

# ==============================
# Export result
# ==============================

filtered_df.to_excel("DEG_GO_filtered_membrane_related.xlsx", index=False)

# ==============================
# Summary
# ==============================

print(f"Total DEGs: {len(deg_genes)}")
print(f"Genes with GO terms: {filtered_df['Gene ID'].nunique()}")
print("Filtered table saved as: DEG_GO_filtered_membrane_related.xlsx")
```

This script produces a final table containing **only genes that are both differentially expressed and associated with the desired GO categories**.

---

# Limitations

A limitation of the GO-based approach is that **many genes in bacterial genomes lack GO annotations**. This is common in newly sequenced or less-studied organisms.

Therefore, relying on GO terms alone may overlook relevant genes.

---

# Final Candidate Gene List

By combining:

* **InterPro keyword searches**
* **GO term–based filtering**

we generated a **refined list of candidate genes** associated with the biological function of interest. These genes were then examined further for their potential role in the biological response observed in the RNA-seq experiment.
# Identifying Genes with Specific Functional Roles from RNA-seq Data

This document describes the workflow used to identify genes associated with a specific biological function (for example **membrane activity**) from RNA-seq differential expression data.
In addition to gene annotations provided in the GTF file (NCBI), two complementary approaches were used to locate candidate genes among the differentially expressed genes (DEGs):

1. **Keyword-based search using InterPro annotations**
2. **Gene Ontology (GO) term mapping**

Because annotation coverage is often incomplete, combining multiple approaches increases the likelihood of identifying relevant genes.

---

# 1. Keyword-Based Search Using InterPro

The **InterPro** database predicts protein family membership, domains, and functional activity based on amino acid sequence. All DEG protein sequences were submitted to the InterPro web server. The resulting output file can be downloaded in **TSV format**.

The InterPro results were then searched for functional keywords using Python.

Example workflow:

```python
import pandas as pd

# Define input file path
input_file = "C:/path/to/DEGs_interpro_output.tsv"

# Load the TSV file (tab-separated values)
df = pd.read_csv(input_file, sep='\t', header=None, comment='#')

# InterPro format: annotation information is typically found in columns 11–12
annotation_columns = [11, 12]

# Define keywords related to the activity of interest
keywords = ['protease', 'hydrolase', 'peptidase', 'proteinase', 'lipase', 'amidase']

# Search for keywords in annotation columns
matches = df[df[annotation_columns].apply(
    lambda row: any(any(k in str(cell).lower() for k in keywords) for cell in row),
    axis=1
)]

# Extract unique gene IDs (usually column 0)
gene_ids = sorted(matches[0].unique())

# Save gene list
with open("cutting_genes.txt", "w") as f:
    for gene in gene_ids:
        f.write(gene + "\n")

# Print results
print("Genes with cutting activity:")
for gene in gene_ids:
    print(gene)
```

This script scans the InterPro annotations for specified keywords and extracts genes whose predicted protein functions match the search terms.

Note that **InterPro cannot assign predictions to every gene**, therefore this method alone may miss relevant candidates.

---

# 2. GO Term–Based Identification

A second strategy uses **Gene Ontology (GO)** annotations to identify genes involved in a specific biological process or molecular function.

This approach consists of **two steps**.

---

## Step A – Identify Genes in the Organism with Relevant GO Terms

First, genes in the organism of interest (*Paenibacillus dendritiformis*) that contain GO terms associated with the desired function were identified.

The script below uses:

* an **Excel file containing all genes with GO annotations**
* the **GO ontology file (`go-basic.obo`)**, which describes relationships between GO terms

The code expands selected **root GO terms** to include all their descendant terms in the GO hierarchy, and then identifies genes annotated with those terms.

```python
import pandas as pd

# ---------- FILE PATHS ----------
excel_file = "C:/Users/goTerms.xlsx"
obo_file = "C:/Users/go-basic.obo"

# ---------- LOAD EXCEL ----------
df = pd.read_excel(excel_file)

# ---------- DEFINE ROOT GO TERMS ----------
parent_go_terms = {
    "GO:0061024": "membrane organization",
    "GO:0099500": "vesicle fusion to plasma membrane",
    "GO:0110165": "cellular anatomical structure",
    "GO:0006886": "intracellular protein transport",
    "GO:0007009": "plasma membrane organization",
    "GO:0071711": "basement membrane organization"
}

# ---------- PARSE OBO FILE ----------
go_children = {}
current_term = None
parents = []

with open(obo_file, "r") as f:
    for line in f:
        line = line.strip()

        if line == "[Term]":
            current_term = None
            parents = []

        elif line.startswith("id: GO:"):
            current_term = line.split("id: ")[1]

        elif line.startswith("is_a: GO:"):
            parent = line.split("is_a: ")[1].split(" !")[0]
            parents.append(parent)

        elif line.startswith("relationship: part_of GO:"):
            parent = line.split("relationship: part_of ")[1].split(" !")[0]
            parents.append(parent)

        elif line == "" and current_term:
            for p in parents:
                go_children.setdefault(p, set()).add(current_term)

# ---------- FUNCTION TO GET ALL DESCENDANTS ----------
def get_all_children(term, children_dict):
    result = set()
    stack = [term]

    while stack:
        t = stack.pop()

        for child in children_dict.get(t, []):
            if child not in result:
                result.add(child)
                stack.append(child)

    return result

# ---------- BUILD ROOT → DESCENDANT MAP ----------
root_to_terms = {}

for root in parent_go_terms:
    expanded = get_all_children(root, go_children)
    expanded.add(root)
    root_to_terms[root] = expanded

# ---------- MATCH GENES ----------
results = []

for _, row in df.iterrows():

    gene = row["Gene ID"]
    go_terms = str(row["GO Term"]).split(";")
    go_terms = [term.strip() for term in go_terms]

    for root, term_set in root_to_terms.items():
        for go_term in go_terms:

            if go_term in term_set:
                results.append({
                    "Gene ID": gene,
                    "GO Term": go_term,
                    "Root GO": root,
                    "Root Description": parent_go_terms[root]
                })

# Convert to DataFrame
results_df = pd.DataFrame(results)

# Remove duplicates
results_df = results_df.drop_duplicates()

# Save detailed table
results_df.to_excel("envelope_related_genes_detailed.xlsx", index=False)

# Save unique gene list
unique_genes = sorted(results_df["Gene ID"].unique())

with open("envelope_related_genes.txt", "w") as f:
    for gene in unique_genes:
        f.write(gene + "\n")

print(f"Total unique envelope-related genes: {len(unique_genes)}")
print("Results saved to:")
print(" - envelope_related_genes_detailed.xlsx")
print(" - envelope_related_genes.txt")

print(results_df["Root Description"].value_counts())
```

The output includes:

* A **detailed table** linking genes to GO terms and root categories
* A **unique gene list** of all genes associated with the function of interest

---

## Step B – Identify Which of These Genes Are Differentially Expressed

The next step is to determine which of the function-related genes identified in Step A are also present in the **DEG list**.

This is done by matching gene IDs between the DEG table and the GO-filtered gene list.

```python
import pandas as pd

# ==============================
# Input files
# ==============================

degs_file = "C:/Users/path/to/all_degs.xlsx"
go_file = "C:/Users/path/to/membrane_related_genes_detailed.xlsx"

# Load DEG list
degs_df = pd.read_excel(degs_file)
deg_genes = set(degs_df["GeneID_clean"].unique())

# Load GO results
go_df = pd.read_excel(go_file)

# ==============================
# Filter GO table for DEGs
# ==============================

filtered_df = go_df[go_df["Gene ID"].isin(deg_genes)].copy()

# ==============================
# Export result
# ==============================

filtered_df.to_excel("DEG_GO_filtered_membrane_related.xlsx", index=False)

# ==============================
# Summary
# ==============================

print(f"Total DEGs: {len(deg_genes)}")
print(f"Genes with GO terms: {filtered_df['Gene ID'].nunique()}")
print("Filtered table saved as: DEG_GO_filtered_membrane_related.xlsx")
```

This script produces a final table containing **only genes that are both differentially expressed and associated with the desired GO categories**.

---

# Limitations

A limitation of the GO-based approach is that **many genes in bacterial genomes lack GO annotations**. This is common in newly sequenced or less-studied organisms.

Therefore, relying on GO terms alone may overlook relevant genes.

---

# Final Candidate Gene List

By combining:

* **InterPro keyword searches**
* **GO term–based filtering**

we generated a **refined list of candidate genes** associated with the biological function of interest. These genes were then examined further for their potential role in the biological response observed in the RNA-seq experiment.
