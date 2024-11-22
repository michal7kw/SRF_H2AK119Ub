#!/bin/bash
#SBATCH --job-name=7_motif_merged
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_motif_analysis_merged.err"
#SBATCH --output="logs/7_motif_analysis_merged.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create necessary directories
mkdir -p analysis/motifs_merged
mkdir -p logs/motifs_merged

# Verify input files exist
if [ ! -f "analysis/diffbind_merged/YAF_enriched.bed" ] || [ ! -f "analysis/diffbind_merged/GFP_enriched.bed" ]; then
    echo "ERROR: Required BED files not found"
    exit 1
fi

# Function to run HOMER motif analysis
run_homer_analysis() {
    local bed_file=$1
    local output_dir=$2
    local genome_file="genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    
    echo "Running HOMER on ${bed_file}..."
    findMotifsGenome.pl \
        "$bed_file" \
        "$genome_file" \
        "$output_dir" \
        -size 200 \
        -mask \
        -p 32 \
        -nomotif \
        -preparsedDir analysis/motifs_merged/preparsed
}

# Run motif analysis for each condition using merged peaks
echo "Running motif analysis for YAF-enriched regions..."
run_homer_analysis \
    "analysis/diffbind_merged/YAF_enriched.bed" \
    "analysis/motifs_merged/YAF"

echo "Running motif analysis for GFP-enriched regions..."
run_homer_analysis \
    "analysis/diffbind_merged/GFP_enriched.bed" \
    "analysis/motifs_merged/GFP"

# Create summary report
echo "Creating motif analysis summary..."
{
    echo "Motif Analysis Summary (Merged Analysis)"
    echo "----------------------------------------"
    echo "YAF-enriched regions:"
    cat analysis/motifs_merged/YAF/knownResults.txt | head -n 10
    echo
    echo "GFP-enriched regions:"
    cat analysis/motifs_merged/GFP/knownResults.txt | head -n 10
} > analysis/motifs_merged/motif_summary.txt

# Check if output files were created
if [ ! -f "analysis/motifs_merged/motif_summary.txt" ]; then
    echo "ERROR: Motif analysis failed to generate output files"
    exit 1
fi

echo "Motif analysis completed successfully" 