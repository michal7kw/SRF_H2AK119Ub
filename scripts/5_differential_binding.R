# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffBind")

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
dir.create("analysis/enrichment", recursive = TRUE, showWarnings = FALSE)

# Read QC metrics
qc_metrics <- read.csv("analysis/qc/qc_metrics.csv", 
                      header=FALSE,
                      col.names=c("Sample", "TotalReads", "ReadsInPeaks", 
                                "FRiP", "PeakCount", "MeanPeakLength"))

# Filter QC metrics to match your samples
qc_metrics_filtered <- qc_metrics[qc_metrics$Sample %in% c(
    paste0("YAF-", 1:3),
    paste0("GFP-", 1:3)
), ]

# Add these prints after reading QC metrics
print("All samples in QC metrics:")
print(qc_metrics$Sample)

print("Filtered samples:")
print(qc_metrics_filtered$Sample)

# Create sample sheet with filtered QC information
samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Factor = rep("H2AK119Ub", 6),
    Condition = rep(c("YAF", "GFP"), each=3),
    Replicate = rep(1:3, 2),
    bamReads = c(
        paste0("analysis/aligned/YAF-", 1:3, ".dedup.bam"),
        paste0("analysis/aligned/GFP-", 1:3, ".dedup.bam")
    ),
    Peaks = c(
        paste0("analysis/peaks/YAF-", 1:3, "_peaks.broadPeak"),
        paste0("analysis/peaks/GFP-", 1:3, "_peaks.broadPeak")
    ),
    PeakCaller = rep("broad", 6),
    FRiP = qc_metrics_filtered$FRiP,
    stringsAsFactors = FALSE
)

# Add these prints after creating the sample sheet
print("Final sample sheet:")
print(samples)

# Quality checks
if (any(qc_metrics_filtered$FRiP < 0.01, na.rm = TRUE)) {
    warning("Low FRiP scores detected (<1%)")
}
if (any(qc_metrics_filtered$PeakCount < 1000, na.rm = TRUE)) {
    warning("Low peak counts detected (<1000)")
}

# Verify peak files
print("Verifying peak files...")
test_peaks <- read.table(samples$Peaks[1], header=FALSE)
print("First few lines of peak file:")
print(head(test_peaks))
print("Peak file structure:")
print(str(test_peaks))

# Add after reading peaks but before DiffBind analysis
print("Performing QC checks...")

# Check peak counts
peak_counts <- sapply(samples$Peaks, function(x) {
  peaks <- read.table(x)
  nrow(peaks)
})
print("Peak counts per sample:")
print(data.frame(Sample=samples$SampleID, Peaks=peak_counts))

# Check peak overlaps
peaks_list <- lapply(samples$Peaks, function(x) {
  read.table(x)[,1:3]
})
names(peaks_list) <- samples$SampleID

# Print overlap statistics
for(i in 1:length(peaks_list)) {
  for(j in i:length(peaks_list)) {
    if(i != j) {
      overlap <- GenomicRanges::countOverlaps(
        makeGRangesFromDataFrame(data.frame(
          chr=peaks_list[[i]][,1],
          start=peaks_list[[i]][,2],
          end=peaks_list[[i]][,3]
        )),
        makeGRangesFromDataFrame(data.frame(
          chr=peaks_list[[j]][,1],
          start=peaks_list[[j]][,2],
          end=peaks_list[[j]][,3]
        ))
      )
      print(sprintf("Overlap between %s and %s: %d peaks",
                   names(peaks_list)[i],
                   names(peaks_list)[j],
                   sum(overlap > 0)))
    }
  }
}

# Check if peak files exist and have content
for (peak_file in samples$Peaks) {
  if (!file.exists(peak_file)) {
    stop(paste("Peak file does not exist:", peak_file))
  }
  peak_content <- read.table(peak_file)
  print(paste(peak_file, "has", nrow(peak_content), "peaks"))
}

# Check format of first peak file
test_peak <- read.table(samples$Peaks[1])
print("Peak file columns:")
print(names(test_peak))
print("First few rows:")
print(head(test_peak))

# Add these checks right before creating the DBA object
print("Checking peak files...")
for (peak_file in samples$Peaks) {
    if (!file.exists(peak_file)) {
        stop(paste("Peak file missing:", peak_file))
    }
    
    peaks <- try(read.table(peak_file))
    if (inherits(peaks, "try-error")) {
        stop(paste("Cannot read peak file:", peak_file))
    }
    
    if (nrow(peaks) == 0) {
        stop(paste("Empty peak file:", peak_file))
    }
    
    print(paste(peak_file, "contains", nrow(peaks), "peaks"))
    print("First few rows:")
    print(head(peaks))
}

# Create DBA object with modified parameters
dba_data <- dba(sampleSheet=samples,
                minOverlap=1,  # Reduced from default
                peakFormat="bed",  # Changed from broadPeak
                peakCaller="raw",  # Changed from broad
                config=data.frame(AnalysisMethod="max",
                                fragmentSize=150,
                                doBlacklist=FALSE,
                                doQC=FALSE))

# Modified counting parameters
dba_data <- dba.count(dba_data,
                      bParallel=FALSE,  # Changed to FALSE to debug
                      summits=FALSE,    # Disabled summits for broad peaks
                      score=DBA_SCORE_READS,
                      bRemoveDuplicates=FALSE,
                      bScaleControl=FALSE,
                      filter=0)

# Add more detailed diagnostics
print("DBA object summary:")
print(dba_data)
print("Peak counts per sample:")
print(dba.show(dba_data)$Intervals)

# Count reads with optimized parameters
print("Counting reads in peaks...")
dba_data <- dba.count(dba_data,
                      score=DBA_SCORE_READS,
                      bParallel=TRUE,
                      summits=FALSE,  # Disable summits calculation since these are broad peaks
                      fragmentSize=150,
                      filter=0,
                      minCount=0,
                      bRemoveDuplicates=FALSE,
                      mapQCth=10,
                      bScaleControl=FALSE,
                      readFormat="bam")

# Print count summary
print("Count summary:")
print(dba.show(dba_data))

# Enhanced QC analysis
print("Performing QC analysis...")
pdf("analysis/diffbind/qc_plots.pdf")
# Sample correlation heatmap
dba.plotHeatmap(dba_data, correlations=TRUE,
                main="Sample Correlations")
# PCA with variance explained
dba.plotPCA(dba_data, 
            attributes = DBA_CONDITION,  # Only plot conditions
            label = DBA_ID,             # Label with sample IDs
            main = "PCA of Sample Data", # Single main title
            sub = "",                    # Empty subtitle
            labelSize = 0.8,            # Adjust label size if needed
            vColors = c("red", "blue"),  # Colors for different conditions
            dotSize = 2)                # Size of points
# Peak overlap rate
dba.plotVenn(dba_data, c(1:3), main="YAF Replicates Overlap")
dba.plotVenn(dba_data, c(4:6), main="GFP Replicates Overlap")
# Add FRiP score plot
barplot(qc_metrics_filtered$FRiP, names.arg=qc_metrics_filtered$Sample,
        main="FRiP Scores", las=2)
dev.off()

# Normalize and analyze
print("Performing differential binding analysis...")
dba_data <- dba.normalize(dba_data, method=DBA_DESEQ2, normalize=DBA_NORM_LIB)
dba_data <- dba.contrast(dba_data, 
                        categories=DBA_CONDITION,
                        minMembers=2)
dba_data <- dba.analyze(dba_data, method=DBA_DESEQ2)

# Extract analysis with multiple thresholds
print("Extracting analysis...")
dba_analysis <- dba.report(dba_data)
significant_analysis <- dba_analysis[dba_analysis$FDR < 0.05 & abs(dba_analysis$Fold) > 1,]

# Save analysis with different stringency levels
write.csv(as.data.frame(dba_analysis), 
          file="analysis/diffbind/all_peaks.csv",
          row.names=FALSE)
write.csv(as.data.frame(significant_analysis), 
          file="analysis/diffbind/significant_peaks.csv",
          row.names=FALSE)

# Enhanced peak annotation
print("Annotating peaks...")
peakAnno <- annotatePeak(significant_analysis,
                         TxDb=txdb,
                         tssRegion=c(-3000, 3000),
                         verbose=FALSE,
                         annoDb="org.Hs.eg.db")

# Generate comprehensive annotation plots
pdf("analysis/annotation/peak_annotation.pdf")
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)
plotAnnoBar(peakAnno)
dev.off()

# Perform pathway enrichment
genes <- unique(as.data.frame(peakAnno)$SYMBOL)
ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

pdf("analysis/enrichment/pathway_analysis.pdf")
dotplot(ego, showCategory=20)
dev.off()

# Generate analysis plots
pdf("analysis/diffbind/analysis_plots.pdf")
# MA plot with highlighted significant peaks
dba.plotMA(dba_data, bXY=TRUE, sub="FDR < 0.05 & |FC| > 1")
# Volcano plot with thresholds
dba.plotVolcano(dba_data)
# Binding patterns
dba.plotBox(dba_data)
dev.off()

print("Analysis completed successfully")

# Print peak statistics
for (i in 1:nrow(samples)) {
  peaks <- read.table(samples$Peaks[i])
  cat(sprintf("Sample %s: %d peaks\n", samples$SampleID[i], nrow(peaks)))
  cat(sprintf("Peak width statistics:\n"))
  print(summary(peaks$V3 - peaks$V2))
}