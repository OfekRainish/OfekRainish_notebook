# Protocol: Installing software with Miniconda

This protocol outlines the steps for installing softwarw using **Miniconda** on a Linux workstation where you don't have administrative write permissions. As an example of this protocol we will download **FastQC** and **MultiQC**. 

### Prerequisites
- **Miniconda** is already installed on the system.
- You have **access** to your home directory where you can create environments and install packages.

---

## Steps
** All comand are ritten in the command line (terminal)

### 1. **Activate Miniconda**

Activate the Miniconda environment by running the following command:

```bash
 source  /usr/local/miniconda3/bin/activate

# This will set up the necessary environment for using conda
```

### 2. Create a Conda Environment in Your Home Directory

```bash
 conda create --prefix ~/environmentName python=3.9

# This command creates a conda environment called environmentName in your home directory (~/environmentName) with Python 3.9.
```
### 3. Activate Your New Environment

```bash
 conda activate ~/environmentName

# Note: Your terminal prompt will change to reflect the active environment, showing (environmentName) instead of (base)
```
### 4. Install 
```bash
 conda install -c bioconda fastqc multiqc


# This installs tools from the bioconda channel, which is a bioinformatics package repository. Here we "installed" fastqc & multiqc, but it can be other tools
```
### 5. Verify Installation
```bash
fastqc --version
multiqc --version

# If the versions are displayed correctly, the tools were successfully installed.
```
______________________________
### Deactivating the Environment (optional)
```bash
conda deactivate

# This will deactivate the current environment, returning you to the base environment (or no environment).
```
### Reactivate
```bash
conda activate ~/environmentName

# similar to step 3
```