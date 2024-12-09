#!/bin/bash
#SBATCH --job-name=snakemake_mrg
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/snakemake_merged_%j.err"
#SBATCH --output="logs/snakemake_merged_%j.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create logs directory
mkdir -p logs/slurm

# Run snakemake with cluster configuration
snakemake \
    --snakefile Snakefile.merged \
    --profile .snakemake/slurm \
    2>&1 | tee logs/snakemake_merged.log