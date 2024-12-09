#!/bin/bash
#SBATCH --job-name=6_peak_anno_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/6_peak_annotation_merged.err"
#SBATCH --output="logs/6_peak_annotation_merged.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Verify input files exist
if [ ! -f "analysis/diffbind_merged/all_peaks.rds" ] || [ ! -f "analysis/diffbind_merged/significant_peaks.rds" ]; then
    echo "ERROR: Required input files from differential binding analysis not found"
    exit 1
fi

# Install required R packages if needed
R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('ChIPseeker', 'clusterProfiler', 'DOSE', 'enrichplot', 'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene'), force = TRUE)"

# Run R script
echo "Starting peak annotation and pathway analysis..."
Rscript scripts/6_peak_annotation_merged.R

# Check if output files were created
if [ ! -f "analysis/annotation_merged/tables/peak_annotation_full.txt" ]; then
    echo "ERROR: Peak annotation failed to generate output files"
    exit 1
fi

echo "Peak annotation and pathway analysis completed successfully" 