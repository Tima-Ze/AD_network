#!/bin/bash 
#SBATCH -J cns
#SBATCH -D /data/scratch2/timaz/wto/script
#SBATCH -o check.out 
#SBATCH --partition=nowick
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=5-00:00:00
#SBATCH --mail-type=end 
#SBATCH --mail-user=timaz@zedat.fu-berlin.de
hostname
date
eval "$(/home/timaz/mambaforge-pypy3/bin/conda shell.bash hook)"
conda activate R
Rscript ./check.R

