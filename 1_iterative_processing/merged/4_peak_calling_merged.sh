#!/bin/bash
#SBATCH --job-name=4_peak_calling_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5  # One job per merged sample
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged/logs/4_peak_calling_merged_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged/logs/4_peak_calling_merged_%a.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Load sample configuration
source config/samples.conf

# Get current merged sample
merged_sample=${MERGED_SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing merged sample: ${merged_sample}"

# Create output directories
mkdir -p analysis/peaks_merged
mkdir -p analysis/visualization_merged
mkdir -p logs/peaks_merged

# Call peaks with MACS2 (CUT&Tag optimized parameters)
echo "Calling peaks for ${merged_sample}..."
macs2 callpeak \
    -t analysis/aligned_merged/tech_reps/${merged_sample}.dedup.bam \
    -f BAMPE \
    -g hs \
    --broad \
    --broad-cutoff 0.1 \
    --nomodel \
    --extsize 147 \
    --shift -75 \
    --keep-dup all \
    -q 0.1 \
    --min-length 200 \
    --outdir analysis/peaks_merged \
    -n ${merged_sample} \
    2> logs/peaks_merged/${merged_sample}.macs2.log

# Generate bigWig file for visualization
echo "Generating bigWig file..."
bamCoverage \
    --bam analysis/aligned_merged/tech_reps/${merged_sample}.dedup.bam \
    --outFileName analysis/visualization_merged/${merged_sample}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --smoothLength 30 \
    --extendReads \
    --centerReads \
    --numberOfProcessors 16

# Generate peak statistics
echo "Generating peak statistics..."
peak_count=$(wc -l < analysis/peaks_merged/${merged_sample}_peaks.broadPeak)
peak_stats=$(awk '{sum+=$3-$2; sumsq+=($3-$2)*($3-$2)} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' \
    analysis/peaks_merged/${merged_sample}_peaks.broadPeak)

echo "${merged_sample} peak statistics:" > analysis/peaks_merged/${merged_sample}_peak_stats.txt
echo "Total peaks: ${peak_count}" >> analysis/peaks_merged/${merged_sample}_peak_stats.txt
echo "Mean peak width: $(echo $peak_stats | cut -d' ' -f1)" >> analysis/peaks_merged/${merged_sample}_peak_stats.txt
echo "Peak width std: $(echo $peak_stats | cut -d' ' -f2)" >> analysis/peaks_merged/${merged_sample}_peak_stats.txt

# Calculate FRiP score
total_reads=$(samtools view -c analysis/aligned_merged/tech_reps/${merged_sample}.dedup.bam)
reads_in_peaks=$(bedtools intersect -a analysis/aligned_merged/tech_reps/${merged_sample}.dedup.bam \
    -b analysis/peaks_merged/${merged_sample}_peaks.broadPeak -u | samtools view -c)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)

echo "FRiP score: ${frip}" >> analysis/peaks_merged/${merged_sample}_peak_stats.txt

echo "Peak calling completed for ${merged_sample}" 

set -e
set -u

# Add version logging
echo "Software versions:"
macs2 --version
bamCoverage --version

# Add input validation
if [ ! -f "analysis/aligned_merged/tech_reps/${merged_sample}.dedup.bam" ]; then
    echo "ERROR: Input BAM file not found for ${merged_sample}"
    exit 1
fi

# Add completion check
if [ ! -f "analysis/peaks_merged/${merged_sample}_peaks.broadPeak" ]; then
    echo "ERROR: Peak calling failed for ${merged_sample}"
    exit 1
fi 