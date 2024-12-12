# Add suppressPackageStartupMessages and enhanced error handling
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(GenomicRanges)
    library(rtracklayer)
    library(DiffBind)
})

# Create output directories
dir.create("analysis/plots_merged", recursive = TRUE, showWarnings = FALSE)

# Read DiffBind results with error handling
tryCatch({
    diff_peaks <- readRDS("analysis/diffbind_merged/significant_peaks.rds")
    all_peaks <- readRDS("analysis/diffbind_merged/all_peaks.rds")
}, error = function(e) {
    stop("Error reading DiffBind results: ", e$message)
})

# Enhanced MA plot
ma_plot <- ggplot(as.data.frame(diff_peaks), aes(x = log2(Conc), y = Fold)) +
    geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
    scale_color_manual(values = c("grey50", "red3")) +
    scale_size_continuous(range = c(0.5, 2)) +
    labs(x = "log2 Mean Concentration",
         y = "log2 Fold Change (YAF/GFP)",
         title = "MA Plot: YAF vs GFP H2AK119Ub Peaks (Merged)",
         subtitle = paste0("Significant peaks: ", sum(diff_peaks$FDR < 0.05))) +
    theme_bw()

# Create volcano plot
volcano_plot <- ggplot(as.data.frame(diff_peaks), 
                      aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey50", "red")) +
  labs(x = "log2 Fold Change (YAF/GFP)",
       y = "-log10(FDR)",
       title = "Volcano Plot: YAF vs GFP H2AK119Ub Peaks") +
  theme_bw() +
  theme(legend.position = "none")

# Save volcano plot
ggsave("analysis/plots/volcano_plot_YAF_vs_GFP.pdf", volcano_plot, width = 8, height = 6)

# Create profile plot using deepTools output if available
if(file.exists("analysis/visualization_merged/YAF_vs_GFP_matrix.gz")) {
  profile_data <- read.table("analysis/visualization_merged/YAF_vs_GFP_matrix.gz", 
                            header = TRUE)
  
  profile_plot <- ggplot(profile_data, aes(x = Position, y = Signal, color = Condition)) +
    geom_line() +
    geom_ribbon(aes(ymin = Signal - SEM, ymax = Signal + SEM, fill = Condition), 
                alpha = 0.2) +
    labs(x = "Distance from peak center (bp)",
         y = "Average H2AK119Ub signal",
         title = "H2AK119Ub Signal Profile") +
    theme_bw()
  
  ggsave("analysis/plots/profile_plot_YAF_vs_GFP.pdf", profile_plot, width = 8, height = 6)
} 
