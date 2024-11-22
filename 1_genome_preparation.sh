#!/bin/bash
#SBATCH --job-name=1_genome_preparation
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/1_genome_preparation.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/logs/1_genome_preparation.out"


# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Download and prepare reference genome
cd genome
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Run bowtie2-build with parallel processing
bowtie2-build --threads 32 Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38
cd ..

echo "Genome preparation completed"