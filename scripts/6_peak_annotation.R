# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)

# Add these functions at the start
standardize_chromosomes <- function(gr) {
  # First, check if chromosomes have 'chr' prefix
  current_chroms <- seqlevels(gr)
  
  # Add 'chr' prefix if missing
  if (!all(grepl("^chr", current_chroms))) {
    seqlevels(gr) <- paste0("chr", seqlevels(gr))
  }
  
  # Define standard chromosomes
  standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
  
  # Keep only standard chromosomes
  gr <- gr[seqnames(gr) %in% standard_chroms]
  
  # Set standard chromosome levels
  seqlevels(gr) <- standard_chroms[standard_chroms %in% seqlevels(gr)]
  
  # Ensure genome info is set
  genome(gr) <- "hg38"
  
  return(gr)
}

# List files and print them for verification
gfp_files <- list.files("./analysis/peaks", pattern="GFP-.*_peaks.broadPeak", full.names=TRUE)
yaf_files <- list.files("./analysis/peaks", pattern="YAF-.*_peaks.broadPeak", full.names=TRUE)

print("GFP files found:")
print(gfp_files)
print("YAF files found:")
print(yaf_files)

# Custom function to read broadPeak files
read_broadPeak <- function(file) {
  # Read the file
  df <- read.table(file, header=FALSE, stringsAsFactors=FALSE)
  
  # Add 'chr' prefix if missing
  if (!all(grepl("^chr", df$V1))) {
    df$V1 <- paste0("chr", df$V1)
  }
  
  # Create GRanges object
  peaks <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3),
    name = df$V4,
    score = round(as.numeric(df$V5)),
    strand = "*",
    signalValue = as.numeric(df$V7),
    pValue = as.numeric(df$V8),
    qValue = as.numeric(df$V9)
  )
  
  # Set genome info
  genome(peaks) <- "hg38"
  
  return(peaks)
}

# Import peaks with error checking
import_peaks_safely <- function(file) {
  tryCatch({
    peaks <- read_broadPeak(file)
    if (length(peaks) == 0) {
      warning(paste("No peaks found in file:", file))
      return(NULL)
    }
    # Standardize chromosome names
    peaks <- standardize_chromosomes(peaks)
    return(peaks)
  }, error = function(e) {
    warning(paste("Error reading file:", file, "\nError:", e$message))
    return(NULL)
  })
}

# Import and check peaks for each condition
gfp_peaks_list <- lapply(gfp_files, import_peaks_safely)
yaf_peaks_list <- lapply(yaf_files, import_peaks_safely)

# Remove any NULL entries
gfp_peaks_list <- Filter(Negate(is.null), gfp_peaks_list)
yaf_peaks_list <- Filter(Negate(is.null), yaf_peaks_list)

# Check if we have valid peaks before merging
if (length(gfp_peaks_list) > 0) {
  gfp_peaks <- reduce(do.call(c, gfp_peaks_list))
  mcols(gfp_peaks)$condition <- "GFP"
} else {
  stop("No valid GFP peaks found")
}

if (length(yaf_peaks_list) > 0) {
  yaf_peaks <- reduce(do.call(c, yaf_peaks_list))
  mcols(yaf_peaks)$condition <- "YAF"
} else {
  stop("No valid YAF peaks found")
}

# Combine all peaks
all_peaks <- c(gfp_peaks, yaf_peaks)

# Get TxDb object
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks
peakAnno <- annotatePeak(all_peaks, 
                        tssRegion=c(-3000, 3000),
                        TxDb=txdb, 
                        annoDb="org.Hs.eg.db")

# Save detailed annotation results
write.table(as.data.frame(peakAnno), 
            file="./analysis/annotation/peak_annotation.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# Generate visualization plots with modified settings
pdf("./analysis/annotation/annotation_plots.pdf")
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
upsetplot(peakAnno)
plotAnnoPie(peakAnno)
dev.off()

# Save gene lists with symbols instead of Entrez IDs
genes_df <- as.data.frame(peakAnno)

# Create a data frame with both IDs and symbols
gene_mapping <- data.frame(
    EntrezID = genes_df$geneId,
    Symbol = mapIds(org.Hs.eg.db,
                   keys = genes_df$geneId,
                   column = "SYMBOL",
                   keytype = "ENTREZID",
                   multiVals = "first")
)

# Remove rows with NA values
gene_mapping <- na.omit(gene_mapping)
# Remove duplicates
gene_mapping <- unique(gene_mapping)

# Save the mapping file
write.table(gene_mapping,
            file = "./analysis/annotation/target_genes_with_ids.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save just the unique gene symbols
write.table(unique(gene_mapping$Symbol),
            file = "./analysis/annotation/target_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
