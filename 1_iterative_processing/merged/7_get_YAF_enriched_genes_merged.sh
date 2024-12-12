#!/bin/bash
#SBATCH --job-name=7_get_YAF_enriched_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_get_YAF_enriched_genes.err"
#SBATCH --output="logs/7_get_YAF_enriched_genes.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Add version logging
echo "Software versions:"
R --quiet -e "sessionInfo()"

# Create output directories
mkdir -p analysis/plots_merged
mkdir -p analysis/gene_lists_merged

# Add input validation
if [ ! -f "analysis/diffbind_merged/significant_peaks.rds" ]; then
    echo "ERROR: Required input file not found"
    exit 1
fi

echo "Generating YAF-enriched gene lists and GO analysis..."
if Rscript scripts/get_YAF_enriched_genes_merged.R; then
    echo "YAF-enriched gene lists and GO analysis completed successfully"
else
    echo "Error: R script failed with exit code $?"
    exit 1
fi

# Add output validation
for file in "analysis/gene_lists_merged/YAF_enriched_H2AK119Ub_genes.txt" \
            "analysis/gene_lists_merged/YAF_enriched_GO_terms.txt"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Expected output file not found: $file"
        exit 1
    fi
done