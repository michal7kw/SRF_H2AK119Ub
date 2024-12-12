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

# Create output directories
mkdir -p analysis/plots
mkdir -p analysis/gene_lists

echo "Generating YAF-enriched gene lists and GO analysis..."
if Rscript scripts/get_YAF_enriched_genes.R; then
    echo "YAF-enriched gene lists and GO analysis completed successfully"
else
    echo "Error: R script failed with exit code $?"
    exit 1
fi