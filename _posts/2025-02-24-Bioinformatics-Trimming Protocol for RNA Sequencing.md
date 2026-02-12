# Trimming Protocol for RNA Sequencing

This protocol explains how to use the **Trim Galore** tool to clean your sequences by removing adapter sequences, low-quality bases, and poly-G chains. **Trim Galore** also includes additional functions that you can explore by running the following command in the Linux terminal:

```bash
trim_galore --help
```
## Prerequisites
1. Download Trim Galore in your working enviroment.
Follow the installation guide [here](../_posts/2025-02-05-Bioinformatics%20-%20Installing%20with%20Miniconda.md) to download Trim Galore.

2. Make Trim Galore Executable
Use the following command to make Trim Galore executable:
```bash
chmod +x trim_galore
```
## Step-by-Step Procedure
1. Make sure you are in the right working enviroment, where trim galore is on.
2. Open **Input Files** Folder.
Ensure you have the folder containing the input files, which should end with *.fastq*.

3. Create a Directory for **Code Files**.
If you don’t have one, create a folder to store the script used for the process done in the enviroment. This will make it easier to locate and manage your scripts.
```bash
#create a directory
mkdir code_scripts
```
4. Create the Script File
Navigate to the newly created folder (code_scripts) and open a text file to input your commands for trim galore:

```bash
#create a directory
cd code_scripts #go into the code_scripts directory
nano trim_galore_code.sh #open a txt file for the trim galore code
```
This will open the nano text editor. After typing in the code, save the file by pressing Ctrl+X, then Y to confirm saving, and Enter to finalize.

5. Input the Code.

Paste the following script into the file:

```bash
#!/bin/bash

# Input and output directories - use 'pwd' to get the absolute path of each directory
input_dir="/complete/path/to/fastq_files"
output_dir="/complete/path/to/new/directory/trim_res"
output_dir_2="/complete/path/to/new/directory/trim_res_polyG"

# Create output directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$output_dir_2"

# First run: Trim adapters + filter reads <70bp
for file in "$input_dir"/*.fastq; do
  trim_galore --length 70 "$file" -o "$output_dir"
done

# Second run: Trim poly-G sequences
for file in "$output_dir"/*_trimmed.fq; do
  trim_galore --nextseq 20 "$file" -o "$output_dir_2"
done

echo "Trimming complete!"

```
6. Explanation of the Code

    * **Input and Output Directories:**

        The code starts by specifying the paths to the input files (input_dir) and two output directories (output_dir and output_dir_2) where the results will be saved.
        The **input_dir** is the folder with your .fastq files, while **output_dir** and **output_dir_2** are for storing results after the initial and secondary processing steps (output directories do not exist yet, and will be opend when running the code).

    * **Initial Processing (Trim Adapters + Filter Short Reads)**

        The first for loop processes each .fastq file in the input_dir folder. The command *trim_galore --length 70* is applied to each file to trim adapter sequences, low-quality bases, and remove reads shorter than 70 base pairs. The processed files are saved in the output_dir folder.

    * **Secondary Processing (Trim Poly-G Sequences)**

        The second for loop processes the files from the output_dir folder (the output from the first processing step). The command *trim_galore --nextseq 20* is applied to remove poly-G sequences, and the processed files are saved in output_dir_2.


                     
7. Results

After the process completes, you will have two directories:

**trim_res**: Contains files processed for adapter trimming and short-read filtering.

**trim_res_polyG**: Contains files processed to remove poly-G sequences.
You can now proceed with further analysis or sequencing steps.


