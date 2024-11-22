#!/bin/bash
#SBATCH --job-name=3_alignment
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/3_alignment_%A_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/3_alignment_%A_%a.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create necessary directories
mkdir -p analysis/aligned
mkdir -p logs/aligned

# Define sample array
samples=(GFP-1 GFP-2 GFP-3 YAF-1 YAF-2 YAF-3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: ${sample}"

# Align with Bowtie2
echo "Aligning ${sample}..."
bowtie2 -p 16 -x genome/GRCh38 \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    -1 analysis/trimmed/${sample}_R1_paired.fastq.gz \
    -2 analysis/trimmed/${sample}_R2_paired.fastq.gz \
    2> logs/aligned/${sample}.align.log | \
    samtools view -@ 16 -bS - > analysis/aligned/${sample}.bam

# Sort BAM
echo "Sorting ${sample}..."
samtools sort -@ 16 analysis/aligned/${sample}.bam \
    -o analysis/aligned/${sample}.sorted.bam

# Index BAM
echo "Indexing ${sample}..."
samtools index analysis/aligned/${sample}.sorted.bam

# Remove duplicates
echo "Removing duplicates for ${sample}..."
picard MarkDuplicates \
    I=analysis/aligned/${sample}.sorted.bam \
    O=analysis/aligned/${sample}.dedup.bam \
    M=logs/aligned/${sample}.metrics.txt \
    REMOVE_DUPLICATES=true

# Index final BAM
samtools index analysis/aligned/${sample}.dedup.bam

# Generate mapping statistics
samtools flagstat analysis/aligned/${sample}.dedup.bam \
    > logs/aligned/${sample}_flagstat.txt

echo "Alignment completed for ${sample}"


####### Serial #######
# # Activate conda environment
# source /opt/common/tools/ric.cosr/miniconda3/bin/activate
# conda activate jupyter_nb

# # Create necessary directories
# mkdir -p analysis/aligned
# mkdir -p logs/aligned

# # Alignment and processing
# for sample in GFP-{1,2,3} YAF-{1,2,3}
# do
#     echo "Processing sample: ${sample}"
    
#     # Align with Bowtie2 (increased threads)
#     echo "Aligning ${sample}..."
#     bowtie2 -p 16 -x genome/GRCh38 \
#         --local --very-sensitive-local \
#         --no-mixed --no-discordant \
#         -1 analysis/trimmed/${sample}_R1_paired.fastq.gz \
#         -2 analysis/trimmed/${sample}_R2_paired.fastq.gz \
#         2> logs/aligned/${sample}.align.log | \
#         samtools view -@ 16 -bS - > analysis/aligned/${sample}.bam
    
#     # Sort BAM (increased threads)
#     echo "Sorting ${sample}..."
#     samtools sort -@ 16 analysis/aligned/${sample}.bam \
#         -o analysis/aligned/${sample}.sorted.bam
    
#     # Index BAM
#     echo "Indexing ${sample}..."
#     samtools index analysis/aligned/${sample}.sorted.bam
    
#     # Remove duplicates
#     echo "Removing duplicates for ${sample}..."
#     picard MarkDuplicates \
#         I=analysis/aligned/${sample}.sorted.bam \
#         O=analysis/aligned/${sample}.dedup.bam \
#         M=logs/aligned/${sample}.metrics.txt \
#         REMOVE_DUPLICATES=true
    
#     # Index final BAM
#     samtools index analysis/aligned/${sample}.dedup.bam
    
#     # Generate mapping statistics
#     samtools flagstat analysis/aligned/${sample}.dedup.bam \
#         > logs/aligned/${sample}_flagstat.txt
# done

# echo "Alignment completed"