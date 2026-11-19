# 从原 scRNA_analysis_151530_all.R 的热图和比例图分节筛选；只补充组图导出。
# 在本目录运行，先运行 calculate_time_subtypes.R。
library(ggplot2)
library(pheatmap)
library(cowplot)

output_dir <- "results"
figure_width <- 12
figure_height <- 7
data <- readRDS(file.path(output_dir, "time_subtype_data.rds"))
subtypes <- c("IA", "ISM", "ISS", "IE", "IR")
cells <- c("B cells", "CAFs", "TECs", "TAMs", "CTL", "Th", "Treg", "Cancer cells")
colors <- setNames(c("#bf5b17", "#66c2a5", "#abd9e9", "#8da0cb",
                     "#d7191c", "#fc8d62", "#a6d854", "#386cb0"), cells)
annotation <- data.frame(TIME_subtype = factor(data$subtype, levels = subtypes),
                         row.names = colnames(data$heatmap))
heatmap <- data$heatmap
rownames(heatmap) <- paste0(rownames(heatmap), "-genes")
groups <- rle(data$subtype)$lengths
gaps <- head(cumsum(groups), -1)
heatmap_plot <- pheatmap(heatmap, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_col = gaps, annotation_col = annotation, annotation_names_col = FALSE,
         annotation_colors = list(TIME_subtype = c(IA = "#E41A1C", ISM = "#377EB8",
                                                  ISS = "#4DAF4A", IE = "#984EA3", IR = "#FF7F00")),
         border_color = NA, fontsize = 10, cellwidth = 10, cellheight = 24,
         main = "Heatmap of TIME Subtype Signatures Across Samples",
         silent = TRUE)
composition_plots <- list()
for (level in c("sample", "subtype")) {
  values <- data[[paste0("by_", level)]]
  values$cell_type <- factor(values$cell_type, levels = cells)
  order <- if (level == "sample") sort(unique(as.character(values$group))) else subtypes
  values$group <- factor(values$group, levels = order)
  plot <- ggplot(values, aes(group, fraction * 100, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(values = colors, drop = FALSE) +
    theme_minimal(base_size = 10) +
    labs(x = if (level == "sample") "Sample ID" else "TIME subtype",
         y = "Proportion (%)", fill = "Cell type",
         title = if (level == "sample") "Proportion of main cell types across samples" else
           "Proportion of main cell types\nacross TIME subtypes") +
    theme(axis.text.x = element_text(angle = if (level == "sample") 45 else 0, hjust = 1),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 8), plot.margin = margin(18, 5, 5, 5))
  composition_plots[[level]] <- plot
}
bottom <- plot_grid(plotlist = composition_plots, nrow = 1, rel_widths = c(1.6, 1),
                    labels = c("(B)", "(C)"), label_size = 12)
figure <- plot_grid(heatmap_plot$gtable, bottom, ncol = 1, rel_heights = c(1, 1),
                    labels = c("(A)", ""), label_size = 12)
ggsave(file.path(output_dir, "time_subtypes.pdf"), figure,
       width = figure_width, height = figure_height, device = cairo_pdf)
ggsave(file.path(output_dir, "time_subtypes.png"), figure,
       width = figure_width, height = figure_height, dpi = 300)
