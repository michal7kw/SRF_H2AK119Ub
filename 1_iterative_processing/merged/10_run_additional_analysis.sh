#!/bin/bash
#SBATCH --job-name=additional_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/additional_analysis.err"
#SBATCH --output="logs/additional_analysis.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Create output directories
mkdir -p analysis/plots
mkdir -p analysis/gene_lists

# Run R scripts
echo "Creating YAF vs GFP plots..."
Rscript scripts/plot_YAF_vs_GFP.R

echo "Generating YAF-enriched gene lists and GO analysis..."
Rscript scripts/get_YAF_enriched_genes.R

echo "Additional analysis completed" 