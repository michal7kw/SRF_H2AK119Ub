#!/bin/bash
#SBATCH --job-name=run_all
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/run_all.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/run_all.out"

# Run all analysis steps
bash 1_genome_preparation.sh
bash 2_quality_control.sh
bash 3_alignment.sh
bash 4_peak_calling.sh
bash 5_differential_binding.sh
bash 6_peak_annotation.sh
bash 7_motif_analysis.sh
# bash scripts/8_sox2_integration.sh
# Rscript scripts/9_visualization.R

echo "Complete analysis pipeline finished"