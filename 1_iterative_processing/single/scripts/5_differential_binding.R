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
dir.create("analysis/diffbind", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/annotation", recursive = TRUE, showWarnings = FALSE)

# Create sample sheet for DiffBind
samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Factor = rep("H2AK119Ub", 6),
    Condition = rep(c("YAF", "GFP"), each=3),
    Replicate = rep(1:3, 2),
    bamReads = file.path("analysis/aligned", 
                        c(paste0("YAF_", 1:3, ".dedup.bam"),
                          paste0("GFP_", 1:3, ".dedup.bam"))),
    Peaks = file.path("analysis/peaks",
                     c(paste0("YAF_", 1:3, "_peaks.broadPeak"),
                       paste0("GFP_", 1:3, "_peaks.broadPeak"))),
    PeakCaller = rep("broad", 6)
)

# Add error checking for file existence
check_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
}

# Check if all files exist
print("Checking input files...")
check_files(samples$bamReads)
check_files(samples$Peaks)

# Create DiffBind object
print("Creating DiffBind object...")
dba_data <- dba(sampleSheet=samples,
                minOverlap=2,
                peakFormat="bed",
                peakCaller="broad",
                config=data.frame(AnalysisMethod="max",
                                fragmentSize=150,
                                doBlacklist=TRUE))

# Count reads in peaks
print("Counting reads in peaks...")
dba_data <- dba.count(dba_data,
                      bParallel=TRUE,
                      summits=FALSE,  # Disable summits for broad peaks
                      score=DBA_SCORE_READS,
                      bRemoveDuplicates=FALSE,
                      fragmentSize=150)

# Generate QC plots
print("Generating QC plots...")
pdf("../analysis/diffbind/qc_plots.pdf")

# Sample correlation heatmap
dba.plotHeatmap(dba_data, correlations=TRUE,
                main="Sample Correlations")
# PCA plot
dba.plotPCA(dba_data, 
            attributes=DBA_CONDITION,
            label=DBA_ID,
            labelSize=0.8)
# Peak overlap rate
dba.plotVenn(dba_data, c(1:3), main="YAF Replicates Overlap")
dba.plotVenn(dba_data, c(4:6), main="GFP Replicates Overlap")
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
saveRDS(res_all, "analysis/diffbind/all_peaks.rds")
saveRDS(res_significant, "analysis/diffbind/significant_peaks.rds")

# Generate analysis plots
pdf("analysis/diffbind/analysis_plots.pdf")
# MA plot
dba.plotMA(dba_data)
# Volcano plot
dba.plotVolcano(dba_data)
# Binding patterns heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n=13)
dba.plotHeatmap(dba_data, correlations=FALSE,
                scale="row", 
                colScheme=hmap)
dev.off()

# Create BED files for visualization
print("Creating BED files for visualization...")
df_all <- as.data.frame(res_significant)

# Update column names to match expected format
names(df_all)[names(df_all) == "seqnames"] <- "Chr"
names(df_all)[names(df_all) == "start"] <- "Start"
names(df_all)[names(df_all) == "end"] <- "End"

# Create separate files for YAF and GFP enriched regions
yaf_enriched <- df_all %>% 
    filter(FDR < 0.05 & Fold > 0) %>% 
    select(Chr, Start, End, Fold, FDR)  # Include Fold and FDR for downstream analysis

gfp_enriched <- df_all %>% 
    filter(FDR < 0.05 & Fold < 0) %>% 
    select(Chr, Start, End, Fold, FDR)

# Save BED files with consistent paths
write.table(yaf_enriched, 
            file="analysis/diffbind/YAF_enriched.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=TRUE)  # Include headers for downstream analysis

write.table(gfp_enriched, 
            file="analysis/diffbind/GFP_enriched.bed", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=TRUE)

# Save the complete results table for downstream analysis
write.table(df_all,
            file="analysis/diffbind/all_differential_peaks.txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE)

print("Differential binding analysis completed")