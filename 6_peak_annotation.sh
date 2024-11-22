#!/bin/bash
#SBATCH --job-name=6_peak_annotation
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/6_peak_annotation.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/6_peak_annotation.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Run R script
Rscript scripts/6_peak_annotation.R