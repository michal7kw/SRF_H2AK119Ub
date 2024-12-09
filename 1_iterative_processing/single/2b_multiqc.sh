#!/bin/bash
#SBATCH --job-name=2_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_multiqc.err"
#SBATCH --output="logs/2_multiqc.out"

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Create directory for MultiQC reports
mkdir -p analysis/multiqc
mkdir -p analysis/qc

# Run MultiQC for pre-trimming data
echo "Creating pre-trimming MultiQC report..."
multiqc \
    --filename "multiqc_pre_trimming" \
    --outdir analysis/multiqc \
    --title "Pre-trimming QC Report" \
    analysis/fastqc/pre_trim/*_fastqc*

# Run MultiQC for post-trimming data
echo "Creating post-trimming MultiQC report..."
multiqc \
    --filename "multiqc_post_trimming" \
    --outdir analysis/multiqc \
    --title "Post-trimming QC Report" \
    analysis/fastqc/post_trim/*_paired_fastqc*

# Generate QC summary
echo "Sample,RawReads,CleanReads,DuplicationRate" > analysis/qc/qc_summary.csv
for sample in GFP_{1,2,3} YAF_{1,2,3}; do
    raw_reads=$(zcat DATA/fastq/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))
    clean_reads=$(zcat analysis/trimmed/${sample}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
    dup_rate=$(grep "PERCENT_DUPLICATION" analysis/qc/${sample}_dup_metrics.txt | cut -f9)
    echo "${sample},${raw_reads},${clean_reads},${dup_rate}" >> analysis/qc/qc_summary.csv
done

echo "MultiQC report generation completed"