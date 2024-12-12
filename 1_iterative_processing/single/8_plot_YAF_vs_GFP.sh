#!/bin/bash
#SBATCH --job-name=8_plot_YAF_vs_GFP
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/8_plot_YAF_vs_GFP.err"
#SBATCH --output="logs/8_plot_YAF_vs_GFP.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Create output directories
mkdir -p analysis/plots
mkdir -p analysis/gene_lists

# Run R scripts
echo "Creating YAF vs GFP plots..."
if Rscript scripts/plot_YAF_vs_GFP.R; then
    echo "plots created successfully"
else
    echo "Error: R script failed with exit code $?"
    exit 1
fi 