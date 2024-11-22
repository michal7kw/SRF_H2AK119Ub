# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# Function to standardize chromosome names
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

# Create output directories
dir.create("analysis/annotation/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/annotation/tables", recursive = TRUE, showWarnings = FALSE)

# Load differential binding results
print("Loading differential binding results...")
all_peaks <- readRDS("analysis/diffbind/all_peaks.rds")
sig_peaks <- readRDS("analysis/diffbind/significant_peaks.rds")

# Standardize chromosome names
print("Standardizing chromosome names...")
all_peaks <- standardize_chromosomes(all_peaks)
sig_peaks <- standardize_chromosomes(sig_peaks)

# Get TxDb object
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks with enhanced parameters for CUT&Tag
print("Annotating peaks...")
tryCatch({
    # Check if input data is valid
    print(paste("Number of significant peaks:", length(sig_peaks)))
    print("Peak data structure:")
    print(str(sig_peaks))
    
    peakAnno <- annotatePeak(sig_peaks, 
                            tssRegion=c(-3000, 3000),
                            TxDb=txdb, 
                            annoDb="org.Hs.eg.db",
                            level="gene")
    
    # Verify annotation results
    print("Annotation complete. Results summary:")
    print(summary(peakAnno))
    
}, error = function(e) {
    print(paste("Error in peak annotation:", e))
    stop("Peak annotation failed. See error message above.")
})

# Generate comprehensive annotation visualizations
print("Generating annotation plots...")
tryCatch({
    pdf("analysis/annotation/figures/annotation_plots.pdf", width=10, height=8)
    
    # Check if plotting is possible
    if (length(peakAnno) == 0) {
        stop("No peaks were annotated - cannot generate plots")
    }
    
    # Add individual try-catch blocks for each plot
    tryCatch({
        plotDistToTSS(peakAnno, title="Peak Distribution Relative to TSS")
    }, error = function(e) print(paste("Error in plotDistToTSS:", e)))
    
    tryCatch({
        plotAnnoBar(peakAnno, title="Peak Annotation Distribution")
    }, error = function(e) print(paste("Error in plotAnnoBar:", e)))
    
    tryCatch({
        plotAnnoPie(peakAnno, title="Genomic Feature Distribution")
    }, error = function(e) print(paste("Error in plotAnnoPie:", e)))
    
    tryCatch({
        upsetplot(peakAnno, vennpie=TRUE, title="Overlap of Different Genomic Features")
    }, error = function(e) print(paste("Error in upsetplot:", e)))
    
    dev.off()
}, error = function(e) {
    print(paste("Error generating plots:", e))
    if (dev.cur() > 1) dev.off()  # Close PDF device if open
})

# Save detailed annotation results
print("Saving annotation results...")
anno_df <- as.data.frame(peakAnno)
write.table(anno_df, 
            file="analysis/annotation/tables/peak_annotation_full.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# Create gene lists with proper symbols
print("Creating gene lists...")
# First, check the EntrezIDs
print("Number of EntrezIDs before mapping:", length(unique(anno_df$geneId)))

# Add error checking and verbose mapping
genes_df <- tryCatch({
    # Try mapping with more detailed error checking
    symbols <- mapIds(org.Hs.eg.db,
                     keys = anno_df$geneId,
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multiVals = "first")
    
    # Print diagnostic information
    print("Number of mapped symbols:", length(symbols[!is.na(symbols)]))
    
    # Create dataframe only with valid mappings
    valid_indices <- !is.na(symbols)
    data.frame(
        EntrezID = anno_df$geneId[valid_indices],
        Symbol = symbols[valid_indices],
        Type = ifelse(anno_df$fold[valid_indices] > 0, "YAF_enriched", "GFP_enriched"),
        stringsAsFactors = FALSE
    )
}, error = function(e) {
    print(paste("Error in gene mapping:", e))
    return(NULL)
})

# Check if genes_df was created successfully
if (is.null(genes_df) || nrow(genes_df) == 0) {
    warning("No genes could be mapped successfully. Stopping analysis.")
    quit(status = 1)
}

# Remove duplicates
genes_df <- unique(genes_df)
print(paste("Final number of mapped genes:", nrow(genes_df)))

# Save gene lists
write.table(genes_df,
            file = "analysis/annotation/tables/target_genes_with_ids.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save separate lists for up/down regulated genes
write.table(unique(subset(genes_df, Type == "YAF_enriched")$Symbol),
            file = "analysis/annotation/tables/YAF_enriched_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(unique(subset(genes_df, Type == "GFP_enriched")$Symbol),
            file = "analysis/annotation/tables/GFP_enriched_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Perform pathway enrichment analysis
print("Performing pathway enrichment analysis...")
# For YAF-enriched genes
yaf_genes <- subset(genes_df, Type == "YAF_enriched")$EntrezID
if(length(yaf_genes) > 0) {
    ego_yaf <- enrichGO(gene = yaf_genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
    
    if(nrow(ego_yaf) > 0) {
        pdf("analysis/annotation/figures/YAF_pathway_analysis.pdf", width=12, height=8)
        print(dotplot(ego_yaf, showCategory=20, title="YAF-enriched Pathways"))
        print(emapplot(ego_yaf, showCategory=30))
        dev.off()
        
        write.csv(as.data.frame(ego_yaf), 
                 "analysis/annotation/tables/YAF_pathway_analysis.csv",
                 row.names=FALSE)
    }
}

# For GFP-enriched genes
gfp_genes <- subset(genes_df, Type == "GFP_enriched")$EntrezID
if(length(gfp_genes) > 0) {
    ego_gfp <- enrichGO(gene = gfp_genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
    
    if(nrow(ego_gfp) > 0) {
        pdf("analysis/annotation/figures/GFP_pathway_analysis.pdf", width=12, height=8)
        print(dotplot(ego_gfp, showCategory=20, title="GFP-enriched Pathways"))
        print(emapplot(ego_gfp, showCategory=30))
        dev.off()
        
        write.csv(as.data.frame(ego_gfp), 
                 "analysis/annotation/tables/GFP_pathway_analysis.csv",
                 row.names=FALSE)
    }
}

print("Peak annotation and pathway analysis completed")
