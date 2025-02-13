# Quality Check Protocol for Tiling Libraries - FastQC and MultiQC

## General Actions:

### 1. Activate Conda
Activate your Conda environment by running the following command:

```bash
source /usr/local/miniconda3/bin/activate
```
### 2. Setting the Path
```bash
conda activate your_environment_name
```
```bash
which toolName 
#for example 'fastqc'
#This will display the path to the fastqc command.
```
```bash
conda info --envs
#optional
#This will list all available Conda environments.
```
```bash 
export PATH=/home/user/miniconda3/envs/bioinf/bin:$PATH                     
#Set the PATH temporarily
```
```bash
echo 'export PATH=/home/user/miniconda3/envs/bioinf/bin:$PATH' >> ~/.bashrc
#Set the PATH permanently
```
## 3. FastQC
```bash 
cd /path/to/your/files
#Change to the directory where your files are stored.
```



```bash 
find /path/to/root_directory -name "*.fastq.gz" -exec fastqc {} \;
#run the tool
#This will generate relevant output files (e.g., .zip, .gz, and .html files) next to each raw data file.
```
### 4. Running MultiQC on FastQC Reports:

```bash 
cd /path/to/your/fastqc/reports
#Navigate to the directory containing the FastQC output files
```
```bash 
multiqc .
#run multiqc on current directory
```
```bash 
multiqc /path/to/your/fastqc/reports
#Alternatively, you can specify the path to the FastQC reports. this is just another way to do it.
```
### Results
![](../images/rna_bioinformatics/libraryQC/fastqc_adapter_content_plot.png)

![](../images/rna_bioinformatics/libraryQC/fastqc_per_base_sequence_quality_plot.png)

![](../images/rna_bioinformatics/libraryQC/fastqc_per_sequence_quality_scores_plot.png)

![](../images/rna_bioinformatics/libraryQC/fastqc_sequence_counts_plot.png)

** all results can be shown on my user in the station: home>RNA_seq>results