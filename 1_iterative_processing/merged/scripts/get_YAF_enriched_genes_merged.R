# Add suppressPackageStartupMessages and enhanced error handling
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(dplyr)
})

# Create output directories
dir.create("analysis/gene_lists_merged", recursive = TRUE, showWarnings = FALSE)

# Read DiffBind results with error checking
diff_peaks_file <- "analysis/diffbind_merged/all_differential_peaks.txt"
if (!file.exists(diff_peaks_file)) {
    stop("Input file not found: ", diff_peaks_file)
}

# Add debugging information
message("Dimensions of input diff_peaks: ", paste(dim(diff_peaks), collapse = " x "))

# Read the DiffBind results
diff_peaks <- readRDS("analysis/diffbind_merged/significant_peaks.rds")

# Filter for YAF-enriched peaks (positive fold change and significant)
yaf_enriched <- diff_peaks[diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05]

# Convert to GRanges object
yaf_peaks_gr <- GRanges(
  seqnames = yaf_enriched$Chr,
  ranges = IRanges(start = yaf_enriched$Start, end = yaf_enriched$End),
  strand = "*",
  fold_change = yaf_enriched$Fold,
  FDR = yaf_enriched$FDR
)

# Annotate peaks with nearby genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(yaf_peaks_gr, 
                        TxDb = txdb,
                        tssRegion = c(-3000, 3000),
                        annoDb = "org.Hs.eg.db")

# Convert annotation to data frame
anno_df <- as.data.frame(peakAnno)

# Enhanced error handling for gene mapping
anno_df$gene_symbol <- tryCatch({
    mapIds(org.Hs.eg.db,
           keys = anno_df$geneId,
           column = "SYMBOL",
           keytype = "ENTREZID",
           multiVals = "first")
}, error = function(e) {
    message("Error mapping gene symbols: ", e$message)
    return(NA)
})

# Create output directory if it doesn't exist
dir.create("analysis/gene_lists", showWarnings = FALSE)

# Write full results to file
write.table(anno_df,
            "analysis/gene_lists/YAF_enriched_H2AK119Ub_full_annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create simplified gene list with key information
simplified_list <- anno_df %>%
  select(gene_symbol, annotation, fold_change = fold_change, FDR,
         distance_to_TSS = distanceToTSS) %>%
  arrange(desc(fold_change))

write.table(simplified_list,
            "analysis/gene_lists/YAF_enriched_H2AK119Ub_genes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Perform GO enrichment analysis on YAF-enriched genes
gene_list <- unique(anno_df$geneId)
go_bp <- enrichGO(gene = gene_list,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

# Save GO results
write.table(as.data.frame(go_bp),
            "analysis/gene_lists/YAF_enriched_GO_terms.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create GO dotplot
pdf("analysis/plots/GO_enrichment_dotplot.pdf", width = 10, height = 8)
dotplot(go_bp, showCategory = 20)
dev.off()

# Print summary statistics
cat("Summary of YAF-enriched H2AK119Ub regions:\n",
    "Total significant peaks:", nrow(yaf_enriched), "\n",
    "Total annotated genes:", length(unique(anno_df$gene_symbol)), "\n",
    "Number of significant GO terms:", nrow(go_bp), "\n",
    file = "analysis/gene_lists/analysis_summary.txt") 
