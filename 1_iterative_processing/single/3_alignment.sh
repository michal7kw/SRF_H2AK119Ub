#!/bin/bash
#SBATCH --job-name=3_alignment
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3_alignment_%a.err"
#SBATCH --output="logs/3_alignment_%a.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate smarcb1_analysis

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Define samples
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: ${sample}"

# # CUT&Tag-optimized alignment with Bowtie2
# echo "Aligning ${sample}..."
# bowtie2 -p 16 \
#     -x ../../genome/GRCh38 \
#     --local --very-sensitive-local \
#     --no-mixed --no-discordant \
#     --no-overlap --no-dovetail \
#     -I 10 -X 700 \
#     -1 analysis/trimmed/${sample}_R1_paired.fastq.gz \
#     -2 analysis/trimmed/${sample}_R2_paired.fastq.gz \
#     2> logs/aligned/${sample}.align.log | \
#     samtools view -@ 16 -bS - > analysis/aligned/${sample}.bam

# # Sort BAM
# echo "Sorting ${sample}..."
# samtools sort -@ 16 \
#     analysis/aligned/${sample}.bam \
#     -o analysis/aligned/${sample}.sorted.bam

# # Index BAM
# echo "Indexing ${sample}..."
# samtools index analysis/aligned/${sample}.sorted.bam

# Mark duplicates
echo "Marking duplicates for ${sample}..."
picard MarkDuplicates \
    I=analysis/aligned/${sample}.sorted.bam \
    O=analysis/aligned/${sample}.dedup.bam \
    M=analysis/qc/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# Index final BAM
samtools index analysis/aligned/${sample}.dedup.bam

# Generate alignment statistics
echo "Generating alignment statistics..."
samtools flagstat analysis/aligned/${sample}.dedup.bam > \
    analysis/qc/${sample}_flagstat.txt

# Generate fragment size distribution
echo "Generating fragment size distribution..."
picard CollectInsertSizeMetrics \
    I=analysis/aligned/${sample}.dedup.bam \
    O=analysis/qc/${sample}_insert_size_metrics.txt \
    H=analysis/qc/${sample}_insert_size_histogram.pdf \
    M=0.5

echo "Alignment completed for ${sample}"