# Add suppressPackageStartupMessages
suppressPackageStartupMessages({
    library(ggplot2)
    library(GenomicRanges)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)
})

# Add error handling for reading files
diff_peaks <- tryCatch({
    readRDS("analysis/diffbind_merged/significant_peaks.rds")
}, error = function(e) {
    stop("Error reading significant peaks: ", e$message)
})

# Enhanced volcano plot
pdf("results/visualization/volcano_plot_merged.pdf", width = 8, height = 7)
ggplot(as.data.frame(all_peaks), aes(x = Fold, y = -log10(FDR))) +
    geom_point(aes(color = FDR < 0.05 & abs(Fold) > 1), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("gray70", "red3")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change (YAF/GFP)",
         y = "-log10(FDR)",
         title = "H2AK119ub Differential Binding (Merged)",
         subtitle = paste0("Significant peaks: ", sum(diff_peaks$FDR < 0.05)),
         color = "Significant") +
    theme(legend.position = "bottom")
dev.off()

# Read overlap data
overlap_data <- read.table("results/integration/overlap_counts.bed")

# Create overlap distribution plot
pdf("results/visualization/overlap_distribution.pdf")
ggplot(overlap_data, aes(x=V4)) + 
    geom_histogram(bins=30, fill="steelblue", color="black") +
    theme_minimal() +
    labs(x="Number of SOX2 overlaps", 
         y="Frequency",
         title="Distribution of H2AK119ub-SOX2 overlaps")
dev.off()

# Create summary statistics
summary_stats <- data.frame(
    total_peaks = nrow(overlap_data),
    peaks_with_sox2 = sum(overlap_data$V4 > 0),
    mean_overlaps = mean(overlap_data$V4),
    median_overlaps = median(overlap_data$V4)
)

write.csv(summary_stats, 
          "results/visualization/overlap_statistics.csv")

print("Visualization completed")