#!/bin/bash
#SBATCH --job-name=9_advanced_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/9_advanced_analysis.err"
#SBATCH --output="logs/9_advanced_analysis.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

echo "Generating advanced analysis..."
if Rscript scripts/10_advanced_analysis.R; then
    echo "advanced analysis completed successfully"
else
    echo "Error: R script failed with exit code $?"
    exit 1
fi

