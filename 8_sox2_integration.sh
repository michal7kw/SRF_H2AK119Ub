#!/bin/bash
#SBATCH --job-name=8_sox2_integration
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/8_sox2_integration.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/8_sox2_integration.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create necessary directories
mkdir -p results/integration
mkdir -p logs/integration

# NOTE: !!!!!!!!!! Replace with the actual path !!!!!!!!!!
SOX2_PEAKS="analysis/peaks/SOX2_peaks.bed"

# Find overlaps between H2AK119ub and SOX2 peaks
echo "Finding overlaps with SOX2 peaks..."
bedtools intersect -a analysis/peaks/merged_peaks.bed \
    -b ${SOX2_PEAKS} \
    -wa -wb > results/integration/h2ak119ub_sox2_overlap.bed

# Generate overlap statistics
bedtools intersect -a analysis/peaks/merged_peaks.bed \
    -b ${SOX2_PEAKS} \
    -c > results/integration/overlap_counts.bed

echo "SOX2 integration analysis completed"