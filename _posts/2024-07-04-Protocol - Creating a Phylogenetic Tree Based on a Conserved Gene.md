# Protocol: Creating a Phylogenetic Tree Based on a Conserved Gene

This protocol outlines the steps to create a phylogenetic tree using the conserved **rbcL** gene from various red algae species as an example.

## Step 1: Selecting the Conserved Gene and Species Group

1. **Choose a Species Group:**
   - in this case  **red algae**.

2. **Select a Conserved Gene:**
   - in this casee **rbcL**, which is conserved in red algae and plants.

## Step 2: Gathering Gene Sequences from NCBI

1. **Access NCBI:**
   - Visit [NCBI](https://www.ncbi.nlm.nih.gov/).

2. **Search for the Gene in Different Species of your selected group:**
   - Enter the search term: `red algae + rbcL`.


   ![alt text](../images/ncb1%20(1).png)
   
3. **Select Species and Download Sequences:**
    
   - Choose the species you want to include in the tree. Here are the selected species for this example:
     - *Gracilaria dura* strain BOR4
     - *Gracilaria gracilis* from United Kingdom
     - *Solieria* sp. isolate LUCXXXX
     - *Ptilophora* sp. voucher Ptilo4
     - *Chondria* sp. isolate PD1961

4. **Download the FASTA Sequences:**
   - Go to the gene’s NCBI page for each species (click on the page).
   - Click on `FASTA` to download the nucleotide sequence.
   - Save these sequences in a Word or Notepad file, including the labels starting with `>`.

   ![alt text](../images/ncb1%20(2).png)

   


## Step 3: Aligning Sequences in Mega11

1. **Open Mega11 Software:**
   - If not installed, [download Mega11](https://www.megasoftware.net/).

2. **Create a New Alignment:**
   - Navigate to `ALIGN` -> `Edit/Build Alignment` -> `Create a New Alignment`.
   - Select `DNA` (or `Protein` depending on your study).

   ![alt text](../images/mega11(1).png)

    ![alt text](../images/mega11(2).png)

    ![alt text](../images/mega11(3).png)


3. **Input Sequences:**
   - Copy the FASTA sequences, including titles starting with `>` for all your species together, and paste them into the alignment window in Mega11.

   ![alt text](../images/mega11(5).png)

4. **Align the Sequences:**
   - Select all sequences by holding `Shift` and clicking.
   - Click `Alignment` -> `Align by ClustalW` in the toolbar.

   ![alt text](../images/mega11(6).png)

5. **Save the Alignment:**
   - Save the aligned sequences to a known location on your computer.

## Step 4: Constructing the Phylogenetic Tree

1. **Open the Tree Construction Tool:**
   - In Mega11, go to `PHYLOGENY` -> `Construct/Test Neighbor-Joining Tree`.

![alt text](../images/mega11(7).png)

2. **Select the Alignment File:**
   - Choose the alignment file saved in the previous step (ensure you check "All files" to view all file types).

3. **View the Tree:**
   - To see the branch lengths, select the `Branch Length` option.

4. **Customize the Tree View:**
   - Use Mega11’s tools to adjust the appearance and labels of the tree as required.
![alt text](../images/mega11(8).png)
---

This protocol provides a step-by-step guide for constructing a phylogenetic tree using the conserved **rbcL** gene from red algae species. By following these steps, you can visualize the evolutionary relationships between the chosen species based on their gene sequences.
