#!/bin/bash
#SBATCH --job-name=5_differential_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/5_differential_binding_%A.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/5_differential_binding_%A.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Install/update BiocManager and DiffBind
# R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DiffBind', force = TRUE)"

# Change to project directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub

# Print current directory and list key files
echo "Current directory: $(pwd)"
echo "Contents of analysis/peaks:"
ls -l analysis/peaks/
echo "Contents of analysis/aligned:"
ls -l analysis/aligned/

# Run R script
Rscript scripts/5_differential_binding.R