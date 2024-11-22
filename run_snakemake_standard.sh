#!/bin/bash
#SBATCH --job-name=snakemake_std
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/snakemake_standard_%j.err"
#SBATCH --output="logs/snakemake_standard_%j.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create logs directory
mkdir -p logs/slurm

# Run snakemake with cluster configuration
# snakemake \
#     --snakefile Snakefile.standard \
#     --profile .snakemake/slurm \
#     2>&1 | tee logs/snakemake_standard.log 

snakemake -j 100 --cluster-config config/cluster.yaml \
    --cluster "sbatch --parsable \
        --partition={resources.partition} \
        --account={resources.account} \
        --mem={resources.mem_mb} \
        --time={resources.time} \
        --cpus-per-task={threads} \
        --job-name=smk-{rule} \
        --output=logs/slurm/{rule}-%j.out \
        --error=logs/slurm/{rule}-%j.err"