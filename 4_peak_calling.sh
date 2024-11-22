#!/bin/bash
#SBATCH --job-name=4_peak_calling
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/4_peak_calling_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/4_peak_calling_%a.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

mkdir -p analysis/peaks
mkdir -p logs/peaks
mkdir -p analysis/visualization
mkdir -p analysis/qc

# Define samples
samples=(YAF_1 YAF_2 YAF_3 GFP_1 GFP_2 GFP_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Generate bigwig for visualization
echo "Generating bigWig file for ${sample}..."
bamCoverage -b analysis/aligned/${sample}.dedup.bam \
    -o analysis/visualization/${sample}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --smoothLength 30 \
    --extendReads \
    --centerReads \
    --numberOfProcessors 16

# Calculate library complexity metrics
echo "Calculating library complexity for ${sample}..."
picard EstimateLibraryComplexity \
    I=analysis/aligned/${sample}.dedup.bam \
    O=analysis/qc/${sample}_complexity.txt

# Call peaks with optimized parameters for H2AK119Ub CUT&Tag
echo "Calling peaks for ${sample}..."
macs2 callpeak \
    -t analysis/aligned/${sample}.dedup.bam \
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
    --max-gap 500 \
    -n ${sample} \
    --outdir analysis/peaks \
    2> logs/peaks/${sample}.log

# Calculate and store QC metrics
echo "Calculating QC metrics for ${sample}..."
# FRiP score
total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
reads_in_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
    -b analysis/peaks/${sample}_peaks.broadPeak -u | samtools view -c)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)

# Peak statistics
peak_count=$(wc -l < analysis/peaks/${sample}_peaks.broadPeak)
mean_peak_length=$(awk '{sum += $3-$2} END {print sum/NR}' analysis/peaks/${sample}_peaks.broadPeak)

# Save QC metrics
echo "${sample},${total_reads},${reads_in_peaks},${frip},${peak_count},${mean_peak_length}" >> analysis/qc/qc_metrics.csv