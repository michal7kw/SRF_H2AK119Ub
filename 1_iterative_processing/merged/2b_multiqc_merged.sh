#!/bin/bash
#SBATCH --job-name=2_multiqc_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_multiqc_merged.err"
#SBATCH --output="logs/2_multiqc_merged.out"

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Load sample configuration
source config/samples.conf

# Create directories
mkdir -p analysis/multiqc_merged
mkdir -p analysis/qc_merged

# Run MultiQC for pre-trimming data
echo "Creating pre-trimming MultiQC report..."
multiqc \
    --filename "multiqc_pre_trimming_merged" \
    --outdir analysis/multiqc_merged \
    --title "Pre-trimming QC Report (Technical Replicates)" \
    analysis/fastqc_merged/pre_trim/*_fastqc*

# Run MultiQC for post-trimming data
echo "Creating post-trimming MultiQC report..."
multiqc \
    --filename "multiqc_post_trimming_merged" \
    --outdir analysis/multiqc_merged \
    --title "Post-trimming QC Report (Technical Replicates)" \
    analysis/fastqc_merged/post_trim/*_paired_fastqc*

# Generate QC summary for technical replicates
echo "Sample,RawReads,CleanReads,DuplicationRate" > analysis/qc_merged/qc_summary_tech_reps.csv
for sample in "${SAMPLES[@]}"; do
    raw_reads=$(zcat DATA/fastq/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))
    clean_reads=$(zcat analysis/trimmed_merged/${sample}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
    dup_rate=$(grep "PERCENT_DUPLICATION" analysis/qc_merged/${sample}_dup_metrics.txt | cut -f9)
    echo "${sample},${raw_reads},${clean_reads},${dup_rate}" >> analysis/qc_merged/qc_summary_tech_reps.csv
done

# Generate QC summary for merged samples
echo "Sample,TotalRawReads,TotalCleanReads,MeanDuplicationRate" > analysis/qc_merged/qc_summary_merged.csv
for merged_sample in "${MERGED_SAMPLES[@]}"; do
    tech_reps=(${TECH_REPS[$merged_sample]})
    total_raw=0
    total_clean=0
    dup_rates=()
    
    for tech_rep in ${tech_reps[@]}; do
        raw_reads=$(zcat ../../DATA/fastq/${tech_rep}_R1_001.fastq.gz | echo $((`wc -l`/4)))
        clean_reads=$(zcat analysis/trimmed_merged/${tech_rep}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
        dup_rate=$(grep "PERCENT_DUPLICATION" analysis/qc_merged/${tech_rep}_dup_metrics.txt | cut -f9)
        
        total_raw=$((total_raw + raw_reads))
        total_clean=$((total_clean + clean_reads))
        dup_rates+=(${dup_rate})
    done
    
    # Calculate mean duplication rate
    mean_dup=$(echo "${dup_rates[@]}" | awk '{sum=0; for(i=1;i<=NF;i++){sum+=$i}; print sum/NF}')
    
    echo "${merged_sample},${total_raw},${total_clean},${mean_dup}" >> analysis/qc_merged/qc_summary_merged.csv
done

echo "MultiQC report generation completed" 