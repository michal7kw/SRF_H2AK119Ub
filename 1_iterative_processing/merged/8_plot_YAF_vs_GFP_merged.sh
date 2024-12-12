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
mkdir -p analysis/plots_merged
mkdir -p analysis/gene_lists_merged

# Add version logging
echo "Software versions:"
R --quiet -e "sessionInfo()"

# Add input validation
if [ ! -f "analysis/diffbind_merged/significant_peaks.rds" ]; then
    echo "ERROR: Required input file not found"
    exit 1
fi

# Run R scripts
echo "Creating YAF vs GFP plots..."
Rscript scripts/plot_YAF_vs_GFP_merged.R

# Add output validation
for file in "analysis/plots_merged/MA_plot_YAF_vs_GFP.pdf" \
           "analysis/plots_merged/volcano_plot_YAF_vs_GFP.pdf"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Expected output file not found: $file"
        exit 1
    fi
done

echo "Additional analysis completed" 