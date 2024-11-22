#!/bin/bash
#SBATCH --job-name=2_quality_control
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB           # Increased memory
#SBATCH --time=INFINITE
#SBATCH --nodes=1           # each array task gets 1 node
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_quality_control_%A_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_quality_control_%A_%a.out"
#SBATCH --array=0-5         # This will create 6 separate jobs, each with 1 node

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Function to check if file is valid gzipped file
check_gzip_integrity() {
    local file=$1
    if ! gzip -t "$file" 2>/dev/null; then
        echo "ERROR: $file is corrupted or incomplete"
        return 1
    fi
    return 0
}

# Function to check disk space
check_disk_space() {
    local dir=$1
    local required_space=20  # Required space in GB
    local available_space=$(df -BG "$dir" | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "$available_space" -lt "$required_space" ]; then
        echo "ERROR: Insufficient disk space in $dir. Available: ${available_space}GB, Required: ${required_space}GB"
        return 1
    fi
    return 0
}

# Function to check if FastQC result exists and is valid
check_fastqc_exists() {
    local fastq=$1
    local output_dir=$2
    local basename=$(basename "$fastq" .fastq.gz)
    if [ -f "${output_dir}/${basename}_fastqc.html" ] && [ -f "${output_dir}/${basename}_fastqc.zip" ]; then
        # Verify zip file integrity
        if unzip -t "${output_dir}/${basename}_fastqc.zip" >/dev/null 2>&1; then
            return 0
        else
            echo "WARNING: Existing FastQC results for $basename are corrupted, will rerun"
            rm -f "${output_dir}/${basename}_fastqc.html" "${output_dir}/${basename}_fastqc.zip"
            return 1
        fi
    fi
    return 1
}

# Function to check if trimmed files exist and are valid
check_trimmed_exists() {
    local sample=$1
    local all_files_valid=true
    
    for file in "analysis/trimmed/${sample}_R1_paired.fastq.gz" \
                "analysis/trimmed/${sample}_R2_paired.fastq.gz" \
                "analysis/trimmed/${sample}_R1_unpaired.fastq.gz" \
                "analysis/trimmed/${sample}_R2_unpaired.fastq.gz"; do
        if [ ! -f "$file" ] || ! check_gzip_integrity "$file"; then
            all_files_valid=false
            # Remove incomplete/corrupted files
            [ -f "$file" ] && rm -f "$file"
        fi
    done
    
    $all_files_valid && return 0 || return 1
}

# Create necessary directories
mkdir -p analysis/fastqc/pre_trim
mkdir -p analysis/fastqc/post_trim
mkdir -p analysis/trimmed
mkdir -p logs/trimming

# Define samples array
declare -a SAMPLES=(
    "GFP-1"
    "GFP-2"
    "GFP-3"
    "YAF-1"
    "YAF-2"
    "YAF-3"
)

# Calculate which samples to process on this node
SAMPLES_PER_NODE=$((${#SAMPLES[@]} / 6))
START_IDX=$((SLURM_ARRAY_TASK_ID * SAMPLES_PER_NODE))
END_IDX=$((START_IDX + SAMPLES_PER_NODE))

echo "Node ${SLURM_ARRAY_TASK_ID} processing samples ${START_IDX} to $((END_IDX-1))"

# Process assigned samples
for ((i=START_IDX; i<END_IDX; i++)); do
    sample=${SAMPLES[$i]}
    echo "Processing sample: $sample"
    
    # Check disk space before processing
    if ! check_disk_space "analysis"; then
        echo "ERROR: Insufficient disk space. Exiting."
        exit 1
    fi

    # Verify input files exist and are valid
    for read in R1 R2; do
        input_fastq="DATA/fastq/${sample}_${read}_001.fastq.gz"
        if [ ! -f "$input_fastq" ]; then
            echo "ERROR: Input file $input_fastq not found"
            exit 1
        fi
        if ! check_gzip_integrity "$input_fastq"; then
            echo "ERROR: Input file $input_fastq is corrupted"
            exit 1
        fi
    done

    # Initial FastQC
    for read in R1 R2; do
        input_fastq="DATA/fastq/${sample}_${read}_001.fastq.gz"
        if ! check_fastqc_exists "$input_fastq" "analysis/fastqc/pre_trim"; then
            echo "Running initial FastQC for ${input_fastq}..."
            if ! fastqc -o analysis/fastqc/pre_trim -t 16 "$input_fastq"; then
                echo "ERROR: FastQC failed for ${input_fastq}"
                exit 1
            fi
        else
            echo "Initial FastQC already exists for ${input_fastq}, skipping..."
        fi
    done

    # Trimming
    if ! check_trimmed_exists "$sample"; then
        echo "Trimming ${sample}..."
        if ! trimmomatic PE -threads 16 \
            DATA/fastq/${sample}_R1_001.fastq.gz \
            DATA/fastq/${sample}_R2_001.fastq.gz \
            analysis/trimmed/${sample}_R1_paired.fastq.gz \
            analysis/trimmed/${sample}_R1_unpaired.fastq.gz \
            analysis/trimmed/${sample}_R2_paired.fastq.gz \
            analysis/trimmed/${sample}_R2_unpaired.fastq.gz \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20 \
            2> logs/trimming/${sample}.log; then
            echo "ERROR: Trimming failed for ${sample}"
            exit 1
        fi
        
        # Verify trimmed files
        if ! check_trimmed_exists "$sample"; then
            echo "ERROR: Trimmed files verification failed for ${sample}"
            exit 1
        fi
    else
        echo "Trimmed files already exist for ${sample}, skipping trimming..."
    fi

    # Post-trimming FastQC
    for type in paired unpaired; do
        for read in R1 R2; do
            trimmed_fastq="analysis/trimmed/${sample}_${read}_${type}.fastq.gz"
            if ! check_fastqc_exists "$trimmed_fastq" "analysis/fastqc/post_trim"; then
                echo "Running post-trim FastQC for ${trimmed_fastq}..."
                if ! fastqc -o analysis/fastqc/post_trim -t 16 "$trimmed_fastq"; then
                    echo "ERROR: Post-trim FastQC failed for ${trimmed_fastq}"
                    exit 1
                fi
            else
                echo "Post-trim FastQC already exists for ${trimmed_fastq}, skipping..."
            fi
        done
    done
done

echo "Quality control completed for node ${SLURM_ARRAY_TASK_ID}"



############### Serial ###############
# #!/bin/bash
# #SBATCH --job-name=2_quality_control
# #SBATCH --account=kubacki.michal
# #SBATCH --mem=128GB
# #SBATCH --time=INFINITE
# #SBATCH --nodes=1
# #SBATCH --ntasks=32
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=kubacki.michal@hsr.it
# #SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_quality_control.err"
# #SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/2_quality_control.out"

# source /opt/common/tools/ric.cosr/miniconda3/bin/activate
# conda activate jupyter_nb

# mkdir -p analysis/fastqc/pre_trim
# mkdir -p analysis/fastqc/post_trim
# mkdir -p analysis/trimmed
# mkdir -p logs/trimming

# # Initial FastQC
# echo "Running initial FastQC..."
# for fastq in DATA/fastq/*fastq.gz
# do
#     fastqc -o analysis/fastqc/pre_trim -t 8 $fastq
# done

# # Modified trimming parameters for CUT&TAG
# echo "Trimming reads..."
# for sample in GFP-{1,2,3} YAF-{1,2,3}
# do
#     trimmomatic PE -threads 8 \
#         DATA/fastq/${sample}_R1_001.fastq.gz \
#         DATA/fastq/${sample}_R2_001.fastq.gz \
#         analysis/trimmed/${sample}_R1_paired.fastq.gz \
#         analysis/trimmed/${sample}_R1_unpaired.fastq.gz \
#         analysis/trimmed/${sample}_R2_paired.fastq.gz \
#         analysis/trimmed/${sample}_R2_unpaired.fastq.gz \
#         ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
#         LEADING:3 TRAILING:3 MINLEN:25 \
#         2> logs/trimming/${sample}.log
# done

# # Post-trimming FastQC
# echo "Running post-trim FastQC..."
# for fastq in analysis/trimmed/*paired.fastq.gz
# do
#     fastqc -o analysis/fastqc/post_trim -t 8 $fastq
# done

# echo "Quality control completed"