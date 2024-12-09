#!/bin/bash
#SBATCH --job-name=5_diff_bind_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding_merged.err"
#SBATCH --output="logs/5_differential_binding_merged.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Print current directory and list key files
echo "Current directory: $(pwd)"
echo "Contents of analysis/peaks_merged:"
ls -l analysis/peaks_merged/
echo "Contents of analysis/aligned_merged:"
ls -l analysis/aligned_merged/

# Verify input files exist
if [ ! -d "analysis/peaks_merged" ] || [ ! -d "analysis/aligned_merged" ]; then
    echo "ERROR: Required input directories not found"
    exit 1
fi

# Install required R packages if needed
# R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('DiffBind', 'ChIPseeker', 'clusterProfiler', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), force = TRUE)"

# Run R script
echo "Starting differential binding analysis..."
Rscript scripts/5_differential_binding_merged.R

# Check if output files were created
if [ ! -f "analysis/diffbind_merged/all_peaks.rds" ] || [ ! -f "analysis/diffbind_merged/significant_peaks.rds" ]; then
    echo "ERROR: Differential binding analysis failed to generate output files"
    exit 1
fi

echo "Differential binding analysis completed successfully" 