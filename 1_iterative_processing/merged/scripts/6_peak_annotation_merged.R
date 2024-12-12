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
dir.create("analysis/annotation_merged/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/annotation_merged/tables", recursive = TRUE, showWarnings = FALSE)

# Load differential binding results
print("Loading differential binding results...")
all_peaks <- readRDS("analysis/diffbind_merged/all_peaks.rds")
sig_peaks <- readRDS("analysis/diffbind_merged/significant_peaks.rds")

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
pdf("analysis/annotation_merged/figures/annotation_plots.pdf", width=10, height=8)

# Distribution of peaks relative to TSS
plotDistToTSS(peakAnno, 
              title="Peak Distribution Relative to TSS (Merged Analysis)")

# Annotation bar plot
plotAnnoBar(peakAnno,
            title="Peak Annotation Distribution (Merged Analysis)")

# Pie chart of genomic feature distribution
plotAnnoPie(peakAnno,
            title="Genomic Feature Distribution (Merged Analysis)")

# Venn diagram of genomic features
upsetplot(peakAnno, 
          vennpie=TRUE,
          title="Overlap of Different Genomic Features (Merged Analysis)")

dev.off()

# Save detailed annotation results
print("Saving annotation results...")
anno_df <- as.data.frame(peakAnno)
write.table(anno_df, 
            file="analysis/annotation_merged/tables/peak_annotation_full.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# Create gene lists with proper symbols
print("Creating gene lists...")
genes_df <- data.frame(
    EntrezID = anno_df$geneId,
    Symbol = mapIds(org.Hs.eg.db,
                   keys = anno_df$geneId,
                   column = "SYMBOL",
                   keytype = "ENTREZID",
                   multiVals = "first"),
    Type = ifelse(anno_df$fold > 0, "YAF_enriched", "GFP_enriched")
)

# Remove NAs and duplicates
genes_df <- na.omit(genes_df)
genes_df <- unique(genes_df)

# Save gene lists
write.table(genes_df,
            file = "analysis/annotation_merged/tables/target_genes_with_ids.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save separate lists for up/down regulated genes
write.table(unique(subset(genes_df, Type == "YAF_enriched")$Symbol),
            file = "analysis/annotation_merged/tables/YAF_enriched_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(unique(subset(genes_df, Type == "GFP_enriched")$Symbol),
            file = "analysis/annotation_merged/tables/GFP_enriched_genes.txt",
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
        pdf("analysis/annotation_merged/figures/YAF_pathway_analysis.pdf", width=12, height=8)
        print(dotplot(ego_yaf, showCategory=20, title="YAF-enriched Pathways (Merged Analysis)"))
        print(emapplot(ego_yaf, showCategory=30))
        print(cnetplot(ego_yaf, showCategory=10, circular=TRUE))
        dev.off()
        
        write.csv(as.data.frame(ego_yaf), 
                 "analysis/annotation_merged/tables/YAF_pathway_analysis.csv",
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
        pdf("analysis/annotation_merged/figures/GFP_pathway_analysis.pdf", width=12, height=8)
        print(dotplot(ego_gfp, showCategory=20, title="GFP-enriched Pathways (Merged Analysis)"))
        print(emapplot(ego_gfp, showCategory=30))
        print(cnetplot(ego_gfp, showCategory=10, circular=TRUE))
        dev.off()
        
        write.csv(as.data.frame(ego_gfp), 
                 "analysis/annotation_merged/tables/GFP_pathway_analysis.csv",
                 row.names=FALSE)
    }
}

# Compare pathways between conditions
if(exists("ego_yaf") && exists("ego_gfp")) {
    pdf("analysis/annotation_merged/figures/pathway_comparison.pdf", width=12, height=8)
    compareCluster(list(YAF=yaf_genes, GFP=gfp_genes),
                  fun="enrichGO",
                  OrgDb=org.Hs.eg.db,
                  ont="BP") %>%
        dotplot(title="Pathway Comparison Between Conditions (Merged Analysis)")
    dev.off()
}

print("Peak annotation and pathway analysis completed") 