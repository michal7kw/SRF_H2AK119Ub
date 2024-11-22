#!/bin/bash

# Create necessary directories
mkdir -p logs

# Function to submit job and get job ID
submit_job() {
    local job_script=$1
    local dependency=$2
    
    if [ -z "$dependency" ]; then
        echo "$(sbatch $job_script | awk '{print $4}')"
    else
        echo "$(sbatch --dependency=afterok:$dependency $job_script | awk '{print $4}')"
    fi
}

echo "Starting CUT&Tag pipeline (standard version)..."

# 1. Genome preparation
echo "1. Submitting genome preparation..."
# genome_job=$(submit_job 1_genome_preparation.sh)

# 2. Quality control
echo "2. Submitting quality control..."
# qc_job=$(submit_job 2_quality_control.sh $genome_job)
qc_job=$(submit_job 2_quality_control.sh)

# 2.1 MultiQC
echo "2.1 Submitting MultiQC..."
multiqc_job=$(submit_job 2_multiqc.sh $qc_job)

# 3. Alignment
echo "3. Submitting alignment..."
align_job=$(submit_job 3_alignment.sh $multiqc_job)

# 4. Peak calling
echo "4. Submitting peak calling..."
peak_job=$(submit_job 4_peak_calling.sh $align_job)

# 5. Differential binding
echo "5. Submitting differential binding analysis..."
diff_job=$(submit_job 5_differential_binding.sh $peak_job)

# 6. Peak annotation
echo "6. Submitting peak annotation..."
anno_job=$(submit_job 6_peak_annotation.sh $diff_job)

# 7. Motif analysis
echo "7. Submitting motif analysis..."
motif_job=$(submit_job 7_motif_analysis.sh $anno_job)

echo "Pipeline submitted successfully!"
echo "Monitor progress with: squeue -u \$USER"
echo "Check logs in: logs/" 