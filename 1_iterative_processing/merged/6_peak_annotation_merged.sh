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

# Add version logging
echo "Software versions:"
R --quiet -e "sessionInfo()"

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

# Add more specific output checks
for file in "analysis/annotation_merged/tables/peak_annotation_full.txt" \
            "analysis/annotation_merged/figures/annotation_plots.pdf"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Expected output file not found: $file"
        exit 1
    fi
done

echo "Peak annotation and pathway analysis completed successfully" 