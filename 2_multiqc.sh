#!/bin/bash
#SBATCH --job-name=2_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_multiqc.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_multiqc.out"

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create directory for MultiQC reports
mkdir -p analysis/multiqc

# Run MultiQC only on post-trimmed paired samples
echo "Creating MultiQC report..."
multiqc \
    --filename "multiqc_report_final" \
    --outdir analysis/multiqc \
    --title "SRF_H2AK119Ub Final Quality Control Report" \
    --comment "FastQC results of paired post-trimmed samples" \
    analysis/fastqc/post_trim/*_paired_fastqc* \
    --ignore "*unpaired*"

echo "MultiQC report generation completed"