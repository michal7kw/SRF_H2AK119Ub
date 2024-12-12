#!/bin/bash
#SBATCH --job-name=3.5_merge_reps
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5  # One job per merged sample
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3.5_merge_replicates_%a.err"
#SBATCH --output="logs/3.5_merge_replicates_%a.out"

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
mkdir -p analysis/aligned_merged/merged
mkdir -p analysis/qc_merged/merged
mkdir -p logs/merged

# Get technical replicates for this merged sample
tech_reps=(${TECH_REPS[$merged_sample]})

# Add input validation
for rep in ${tech_reps[@]}; do
    if [ ! -f "analysis/aligned_merged/tech_reps/${rep}.dedup.bam" ]; then
        echo "ERROR: Input BAM file not found for replicate ${rep}"
        exit 1
    fi
done

# Merge BAM files
echo "Merging technical replicates for ${merged_sample}..."
samtools merge -@ 16 \
    analysis/aligned_merged/merged/${merged_sample}.merged.bam \
    $(for rep in ${tech_reps[@]}; do echo analysis/aligned_merged/tech_reps/${rep}.dedup.bam; done)

# Sort merged BAM
echo "Sorting merged BAM..."
samtools sort -@ 16 \
    analysis/aligned_merged/merged/${merged_sample}.merged.bam \
    -o analysis/aligned_merged/merged/${merged_sample}.merged.sorted.bam

# Index merged BAM
echo "Indexing merged BAM..."
samtools index analysis/aligned_merged/merged/${merged_sample}.merged.sorted.bam

# Generate merged statistics
echo "Generating merged statistics..."
samtools flagstat analysis/aligned_merged/merged/${merged_sample}.merged.sorted.bam \
    > analysis/qc_merged/merged/${merged_sample}.merged.flagstat

# Generate fragment size distribution for merged BAM
echo "Generating fragment size distribution..."
picard CollectInsertSizeMetrics \
    I=analysis/aligned_merged/merged/${merged_sample}.merged.sorted.bam \
    O=analysis/qc_merged/merged/${merged_sample}.merged_insert_size_metrics.txt \
    H=analysis/qc_merged/merged/${merged_sample}.merged_insert_size_histogram.pdf \
    M=0.5

# Clean up intermediate files
rm analysis/aligned_merged/merged/${merged_sample}.merged.bam

echo "Technical replicate merging completed for ${merged_sample}"

# Add completion check
if [ ! -f "analysis/aligned_merged/merged/${merged_sample}.merged.sorted.bam" ]; then
    echo "ERROR: Merging failed for ${merged_sample}"
    exit 1
fi 