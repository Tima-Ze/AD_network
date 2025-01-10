# Nextflow Workflow for Running R Script (Calls_wTO.R)

This repository contains a Nextflow workflow (main.nf), that includes all scripts (Rscripts) to construct the wTO networks, check the network topology in each step, filter the links of each network based on cutoff, and finally to integrate networks to make a Consensus networks.

## Cloning the Repository and Setting Up

First, clone the git repository to your local machine:

```bash
git clone https://github.com/Tima-Ze/AD_network.git
```

After cloning, make sure to make the scripts in the `bin` directory executable by running:

```bash
chmod +x bin/*
```

The workflow is configured to run on a **SLURM** cluster and uses a **Conda environment** for managing R dependencies. For the most smooth usage, create a Conda environment from the `R.yml` file by running the following command:

```bash
conda env create -f R.yml
```

This will set up the necessary environment with all required dependencies specified in the `R.yml` file.

---

Before runing the workflow, set the config file for these parameters:

### Setting Up Conda for Nextflow

To run this Nextflow pipeline on your personal computer, you need to specify the path to your Conda installation.

**Steps to configure:**

1. **Install Anaconda**: If you haven't already, install Anaconda on your system. You can download it from the [Anaconda website](https://www.anaconda.com/products/distribution). It is also recommended to install Mamba for faster package management. You can install Mamba by running:

   ```bash
   conda install mamba -n base -c conda-forge
   ```

2. **Find your Conda path**: Determine the path where Conda is installed on your system. This is typically something like `/home/yourusername/miniconda3` or `/home/yourusername/mambaforge`.

3. **Update nextflow.config**: Open the `nextflow.config` file in a text editor and update the `conda` line with your Conda installation path. For example:
   ```nextflow-config
   conda = '/home/yourusername/miniconda3'
   ```

By setting this path, Nextflow will be able to use Conda to manage the software dependencies for the pipeline.

### Setting Bootstrap Iterations

Please refer to the wTO R package [tutorial](https://deisygysi.github.io/rpackages/wto/) for detailed instructions on setting bootstrap iterations.

### Adjusting wTO Process Parameters

You may customize the following parameters for the wTO process according to your requirements:

- **Memory**: Adjust the memory allocation.
- **CPU**: Modify the number of CPU cores to be used.

---

## Workflow Overview

### 1- Network construction

- Ensure your input data (gene expression data) is placed in the `Data` directory. The tables should be formatted with samples in columns and genes in rows, with gene names specified in the row names.
- The wTO networks constructed from each input will be stored in the `Results/Raw_wTO` directory.
- The process (1) wTO calculates the co-expression correlations by these parameters, `method='p'` (pearson coefficient), `sign='sign'`, and `delta=0.2`. For more details please refer to wTO R package [tutorial] (https://deisygysi.github.io/rpackages/wto/). You are able to determin the number of bootstrap through `nextflow.config` file.

## Directory Structure

```plaintext
.
├── main.nf                   # The main Nextflow workflow file
├── nextflow.config            # The configuration file for SLURM and Conda
├── Results/Raw_wTO/           # Directory where the output files will be stored
├── Calls_wTO.R                # The R script that will be executed
└── README.md                  # This README file
```

## Running the Nextflow Workflow

nextflow run -c nextflow.config main.nf
