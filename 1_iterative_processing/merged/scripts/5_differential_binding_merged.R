# Load required libraries
library(DiffBind)
library(tidyverse)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(clusterProfiler)

# Create necessary directories
dir.create("analysis/diffbind_merged", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/annotation_merged", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/enrichment_merged", recursive = TRUE, showWarnings = FALSE)

# Load sample configuration
source("config/samples.conf")

# Function to check if files exist and have content
check_files <- function(samples_df) {
  missing_files <- character()
  empty_files <- character()
  
  for(i in 1:nrow(samples_df)) {
    bam_file <- samples_df$bamReads[i]
    peak_file <- samples_df$Peaks[i]
    
    # Check BAM files
    if(!file.exists(bam_file)) {
      missing_files <- c(missing_files, bam_file)
    }
    
    # Check peak files
    if(!file.exists(peak_file)) {
      missing_files <- c(missing_files, peak_file)
    } else if(file.info(peak_file)$size == 0) {
      empty_files <- c(empty_files, peak_file)
    }
  }
  
  if(length(missing_files) > 0) {
    cat("Missing files:\n")
    cat(paste(missing_files, collapse="\n"))
    cat("\n")
    stop("Required input files are missing. Please check the paths.")
  }
  
  if(length(empty_files) > 0) {
    cat("Empty files:\n")
    cat(paste(empty_files, collapse="\n"))
    cat("\n")
    stop("Some peak files are empty. Please check the files.")
  }
}

# Create sample sheet for DiffBind
samples <- data.frame(
    SampleID = MERGED_SAMPLES,
    Factor = rep("H2AK119Ub", length(MERGED_SAMPLES)),
    Condition = sub("_[0-9]+$", "", MERGED_SAMPLES),
    Replicate = as.numeric(sub("^.*_", "", MERGED_SAMPLES)),
    bamReads = file.path("analysis/aligned_merged", 
                        paste0(sub("_", "-", MERGED_SAMPLES), ".dedup.bam")),
    Peaks = file.path("analysis/peaks_merged", 
                     paste0(sub("_", "-", MERGED_SAMPLES), "_peaks.broadPeak")),
    PeakCaller = rep("broad", length(MERGED_SAMPLES))
)

# Print current working directory and sample sheet for debugging
print("Current working directory:")
print(getwd())
print("\nSample sheet:")
print(samples)

# Check if all required files exist and have content
check_files(samples)

# Create DiffBind object with more verbose output
print("Creating DiffBind object...")
dba_data <- tryCatch({
    dba(sampleSheet=samples,
        minOverlap=2,
        peakFormat="bed",
        peakCaller="raw",
        config=data.frame(AnalysisMethod="max",
                         fragmentSize=150,
                         doBlacklist=TRUE))
}, error = function(e) {
    print("Error creating DiffBind object:")
    print(e)
    stop(e)
})

# Count reads in peaks
print("Counting reads in peaks...")
tryCatch({
    dba_data <- dba.count(dba_data,
                          bParallel=TRUE,
                          summits=FALSE,
                          score=DBA_SCORE_READS,
                          bRemoveDuplicates=FALSE,
                          fragmentSize=150)
    print("Successfully counted reads in peaks")
}, error = function(e) {
    print("Error in dba.count:")
    print(e)
    print("Checking if BAM files exist and are readable:")
    for(bam in samples$bamReads) {
        print(paste("BAM file:", bam, "exists:", file.exists(bam)))
    }
    stop(e)
})

# Generate QC plots
print("Generating QC plots...")
pdf("analysis/diffbind_merged/qc_plots.pdf")
# Sample correlation heatmap
dba.plotHeatmap(dba_data, correlations=TRUE,
                main="Sample Correlations (Merged Analysis)")
# PCA plot
dba.plotPCA(dba_data, 
            attributes=DBA_CONDITION,
            label=DBA_ID,
            labelSize=0.8,
            main="PCA Plot (Merged Analysis)")
# Peak overlap rate
dba.plotVenn(dba_data, c(1:3), main="YAF Replicates Overlap (Merged)")
dba.plotVenn(dba_data, c(4:6), main="GFP Replicates Overlap (Merged)")
dev.off()

# Perform differential binding analysis
print("Performing differential binding analysis...")
dba_data <- dba.normalize(dba_data, method=DBA_DESEQ2)
dba_data <- dba.contrast(dba_data, 
                        categories=DBA_CONDITION,
                        minMembers=2)
dba_data <- dba.analyze(dba_data, method=DBA_DESEQ2)

# Extract results
print("Extracting differential binding results...")
res_all <- dba.report(dba_data)
res_significant <- dba.report(dba_data, th=0.05)

# Save results
print("Saving results...")
saveRDS(res_all, "analysis/diffbind_merged/all_peaks.rds")
saveRDS(res_significant, "analysis/diffbind_merged/significant_peaks.rds")

# Generate analysis plots
pdf("analysis/diffbind_merged/analysis_plots.pdf")
# MA plot
dba.plotMA(dba_data, main="MA Plot (Merged Analysis)")
# Volcano plot
dba.plotVolcano(dba_data, main="Volcano Plot (Merged Analysis)")
# Binding patterns heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n=13)
dba.plotHeatmap(dba_data, correlations=FALSE,
                scale="row", 
                colScheme=hmap,
                main="Binding Patterns (Merged Analysis)")
dev.off()

# Create BED files for visualization
print("Creating BED files for visualization...")
df_all <- as.data.frame(res_significant)
yaf_enriched <- df_all %>% 
    filter(FDR < 0.05 & Fold > 0) %>% 
    select(seqnames, start, end)
gfp_enriched <- df_all %>% 
    filter(FDR < 0.05 & Fold < 0) %>% 
    select(seqnames, start, end)

write.table(yaf_enriched, 
            file="analysis/diffbind_merged/YAF_enriched.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)
write.table(gfp_enriched, 
            file="analysis/diffbind_merged/GFP_enriched.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)

print("Differential binding analysis completed")