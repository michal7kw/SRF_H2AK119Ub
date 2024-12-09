#!/bin/bash
#SBATCH --job-name=5_differential_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding_%A.err"
#SBATCH --output="logs/5_differential_binding_%A.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Print current directory and list key files
echo "Current directory: $(pwd)"
echo "Contents of analysis/peaks:"
ls -l analysis/peaks/
echo "Contents of analysis/aligned:"
ls -l analysis/aligned/

# Install required R packages if needed
# R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('DiffBind', 'ChIPseeker', 'clusterProfiler', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), force = TRUE)"

# Run R script
Rscript ../scripts/5_differential_binding.R