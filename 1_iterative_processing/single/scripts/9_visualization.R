# Load required libraries
library(ggplot2)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)

# Read overlap data
overlap_data <- read.table("../results/integration/overlap_counts.bed")

# Create overlap distribution plot
pdf("../results/visualization/overlap_distribution.pdf")
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
          "../results/visualization/overlap_statistics.csv")

print("Visualization completed")