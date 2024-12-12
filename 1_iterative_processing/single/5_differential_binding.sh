#!/bin/bash
#SBATCH --job-name=5_differential_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding.err"
#SBATCH --output="logs/5_differential_binding.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

# Change to the correct directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Create log directory if it doesn't exist
mkdir -p logs

# Print current directory and check for required directories
echo "Current directory: $(pwd)"
echo "Checking required directories..."
for dir in analysis/peaks analysis/aligned; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory $dir does not exist"
        exit 1
    fi
    echo "Contents of $dir:"
    ls -l "$dir"
done

# Install required R packages if needed
# R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('DiffBind', 'ChIPseeker', 'clusterProfiler', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), force = TRUE)"

# Run R script
Rscript scripts/5_differential_binding.R