#!/bin/bash
#SBATCH --job-name=1_genome_preparation
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/1_genome_preparation.err"
#SBATCH --output="logs/1_genome_preparation.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/merged

# Create necessary directories
mkdir -p genome
mkdir -p logs
mkdir -p analysis/fastqc/pre_trim
mkdir -p analysis/fastqc/post_trim
mkdir -p analysis/trimmed
mkdir -p analysis/aligned
mkdir -p analysis/peaks
mkdir -p analysis/visualization
mkdir -p analysis/qc
mkdir -p analysis/annotation

# Download and prepare reference genome
cd ../genome

# Download primary assembly to avoid alternative contigs
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build bowtie2 indices
echo "Building human genome index..."
bowtie2-build --threads 32 Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38


# # Generate chromosome sizes file for downstream analysis
# samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
# cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > hg38.chrom.sizes

cd ..

echo "Genome preparation completed"