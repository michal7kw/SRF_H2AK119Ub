#!/bin/bash
#SBATCH --job-name=2_quality_control_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_quality_control_merged_%a.err"
#SBATCH --output="logs/2_quality_control_merged_%a.out"
#SBATCH --array=0-11  # Adjusted for all technical replicates

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Load sample configuration
source config/samples.conf

# Create separate directories for merged pipeline
mkdir -p analysis/fastqc_merged/pre_trim
mkdir -p analysis/fastqc_merged/post_trim
mkdir -p analysis/trimmed_merged
mkdir -p logs/trimming_merged

# Get current sample
sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: $sample"

# Initial FastQC
echo "Running initial FastQC for ${sample}..."
fastqc -o analysis/fastqc_merged/pre_trim -t 16 \
    DATA/fastq/${sample}_R1_001.fastq.gz \
    DATA/fastq/${sample}_R2_001.fastq.gz

# CUT&Tag-specific trimming parameters
echo "Trimming ${sample}..."
trimmomatic PE -threads 16 \
    DATA/fastq/${sample}_R1_001.fastq.gz \
    DATA/fastq/${sample}_R2_001.fastq.gz \
    analysis/trimmed_merged/${sample}_R1_paired.fastq.gz \
    analysis/trimmed_merged/${sample}_R1_unpaired.fastq.gz \
    analysis/trimmed_merged/${sample}_R2_paired.fastq.gz \
    analysis/trimmed_merged/${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25 \
    2> logs/trimming_merged/${sample}.log

# Post-trimming FastQC
echo "Running post-trim FastQC for ${sample}..."
fastqc -o analysis/fastqc_merged/post_trim -t 16 \
    analysis/trimmed_merged/${sample}_R1_paired.fastq.gz \
    analysis/trimmed_merged/${sample}_R2_paired.fastq.gz

echo "Quality control completed for ${sample}" 