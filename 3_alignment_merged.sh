#!/bin/bash
#SBATCH --job-name=3_alignment_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5  # Adjusted for all technical replicates
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3_alignment_merged_%a.err"
#SBATCH --output="logs/3_alignment_merged_%a.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Load sample configuration
source config/samples.conf

# Get current sample
sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: ${sample}"

# Create output directories
mkdir -p analysis/aligned_merged/tech_reps
mkdir -p analysis/qc_merged/tech_reps
mkdir -p logs/aligned_merged

# CUT&Tag-optimized alignment with Bowtie2
echo "Aligning ${sample}..."
bowtie2 -p 16 \
    -x genome/GRCh38 \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    --no-overlap --no-dovetail \
    -I 10 -X 700 \
    -1 analysis/trimmed_merged/${sample}_R1_paired.fastq.gz \
    -2 analysis/trimmed_merged/${sample}_R2_paired.fastq.gz \
    2> logs/aligned_merged/${sample}.align.log | \
    samtools view -@ 16 -bS - > analysis/aligned_merged/tech_reps/${sample}.bam

# Sort BAM
echo "Sorting ${sample}..."
samtools sort -@ 16 \
    analysis/aligned_merged/tech_reps/${sample}.bam \
    -o analysis/aligned_merged/tech_reps/${sample}.sorted.bam

# Index BAM
echo "Indexing ${sample}..."
samtools index analysis/aligned_merged/tech_reps/${sample}.sorted.bam

# Mark duplicates
echo "Marking duplicates for ${sample}..."
picard MarkDuplicates \
    I=analysis/aligned_merged/tech_reps/${sample}.sorted.bam \
    O=analysis/aligned_merged/tech_reps/${sample}.dedup.bam \
    M=analysis/qc_merged/tech_reps/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# Index final BAM
samtools index analysis/aligned_merged/tech_reps/${sample}.dedup.bam

# Generate alignment statistics
echo "Generating alignment statistics..."
samtools flagstat analysis/aligned_merged/tech_reps/${sample}.dedup.bam > \
    analysis/qc_merged/tech_reps/${sample}_flagstat.txt

# Generate fragment size distribution
echo "Generating fragment size distribution..."
picard CollectInsertSizeMetrics \
    I=analysis/aligned_merged/tech_reps/${sample}.dedup.bam \
    O=analysis/qc_merged/tech_reps/${sample}_insert_size_metrics.txt \
    H=analysis/qc_merged/tech_reps/${sample}_insert_size_histogram.pdf \
    M=0.5

# Clean up intermediate files
rm analysis/aligned_merged/tech_reps/${sample}.bam
rm analysis/aligned_merged/tech_reps/${sample}.sorted.bam
rm analysis/aligned_merged/tech_reps/${sample}.sorted.bam.bai

echo "Alignment completed for ${sample}" 