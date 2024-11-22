#!/bin/bash
#SBATCH --job-name=SRF_H2AK119Ub
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/snakemake_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/snakemake_%j.out"

# Load necessary modules or activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb/bin/activate

# Create logs directory
mkdir -p logs/slurm

# Set temporary directory to beegfs
export TMPDIR=/beegfs/scratch/ric.broccoli/kubacki.michal/tmp
mkdir -p $TMPDIR

# Run Snakemake
snakemake \
    --latency-wait 60 \
    --printshellcmds \
    --keep-going \
    --rerun-incomplete \
    --jobs 1 \
    --verbose