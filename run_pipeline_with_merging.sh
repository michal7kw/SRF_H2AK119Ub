#!/bin/bash

# Create necessary directories
mkdir -p logs

# Function to submit job and get job ID
submit_job() {
    local job_script=$1
    local dependency=$2
    local job_id
    
    if [ -z "$dependency" ]; then
        job_id=$(sbatch "$job_script" | grep -oP '\d+')
    else
        job_id=$(sbatch --dependency=afterok:$dependency "$job_script" | grep -oP '\d+')
    fi
    
    # Verify job submission
    if [ -z "$job_id" ]; then
        echo "Error: Failed to submit job $job_script" >&2
        return 1
    fi
    
    echo "Submitted job $job_script with ID: $job_id" >&2
    echo "$job_id"
}

echo "Starting CUT&Tag pipeline with technical replicate merging..."

# Clear any existing jobs (optional)
scancel -u $USER

# 1. Genome preparation
echo "1. Submitting genome preparation..."
# genome_job=$(submit_job 1_genome_preparation.sh)
# echo "Genome preparation job ID: $genome_job"

# 2. Quality control for technical replicates
echo "2. Submitting quality control for technical replicates..."
# qc_job=$(submit_job 2_quality_control_merged.sh $genome_job)
qc_job=$(submit_job 2_quality_control_merged.sh)

echo "Quality control job ID: $qc_job"

# 2.1 MultiQC
echo "2.1 Submitting MultiQC..."
multiqc_job=$(submit_job 2_multiqc_merged.sh $qc_job)
echo "MultiQC job ID: $multiqc_job"

# 3. Alignment of technical replicates
echo "3. Submitting alignment of technical replicates..."
align_job=$(submit_job 3_alignment_merged.sh $multiqc_job)
echo "Alignment job ID: $align_job"

# 3.5 Merge technical replicates
echo "3.5 Submitting technical replicate merging..."
merge_job=$(submit_job 3.5_merge_technical_replicates.sh $align_job)
echo "Merge job ID: $merge_job"

# 4. Peak calling on merged data
echo "4. Submitting peak calling on merged data..."
peak_job=$(submit_job 4_peak_calling_merged.sh $merge_job)
echo "Peak calling job ID: $peak_job"

# 5. Differential binding on merged data
echo "5. Submitting differential binding analysis..."
diff_job=$(submit_job 5_differential_binding_merged.sh $peak_job)
echo "Differential binding job ID: $diff_job"

# 6. Peak annotation on merged data
echo "6. Submitting peak annotation..."
anno_job=$(submit_job 6_peak_annotation_merged.sh $diff_job)
echo "Peak annotation job ID: $anno_job"

# 7. Motif analysis on merged data
echo "7. Submitting motif analysis..."
motif_job=$(submit_job 7_motif_analysis_merged.sh $anno_job)
echo "Motif analysis job ID: $motif_job"

# Print job dependency chain
echo -e "\nJob dependency chain:"
echo "genome_job ($genome_job) -> qc_job ($qc_job) -> multiqc_job ($multiqc_job) -> align_job ($align_job) -> merge_job ($merge_job) -> peak_job ($peak_job) -> diff_job ($diff_job) -> anno_job ($anno_job) -> motif_job ($motif_job)"

echo -e "\nPipeline submitted successfully!"
echo "Monitor progress with: squeue -u \$USER"
echo "Check logs in: logs/" 