# 从原脚本的 netVisual_circle 和 A_Function.R 的散点图分节筛选；只补充组图导出。
# 在本目录运行，先运行 calculate_communication.R。
library(CellChat)
library(ggplot2)
library(ggrepel)
library(cowplot)

output_dir <- "results"
figure_width <- 12
figure_height <- 10
data <- readRDS(file.path(output_dir, "communication_data.rds"))
network_plots <- lapply(names(data$networks), function(name) {
  net <- data$networks[[name]]
  as_grob(function() {
    par(mar = c(1.5, 0, 1.5, 0), xpd = TRUE)
    netVisual_circle(net$weight, vertex.weight = net$cells, weight.scale = TRUE,
                     label.edge = FALSE, title.name = paste("Interaction Strength of", name))
  })
})
colors <- setNames(scPalette(length(data$cell_types)), data$cell_types)
role_plots <- lapply(c("TNF", "IFN-II", "TGFb"), function(pathway) {
  values <- data$roles[data$roles$pathway == pathway, , drop = FALSE]
  values$cell_type <- factor(values$cell_type, levels = data$cell_types)
  plot <- ggplot(values, aes(outgoing, incoming)) +
    geom_point(aes(size = Count, colour = cell_type, fill = cell_type)) +
    CellChat_theme_opts() +
    scale_colour_manual(values = colors, drop = FALSE) +
    scale_fill_manual(values = alpha(colors, 0.6), drop = FALSE) +
    scale_size_continuous(limits = data$size_limits, range = c(2, 6)) +
    geom_text_repel(aes(label = cell_type, colour = cell_type), size = 3,
                    segment.size = 0.2, segment.alpha = 0.5, seed = 123, show.legend = FALSE) +
    guides(colour = "none", fill = "none") +
    labs(title = paste("S004 -", pathway), x = "Outgoing interaction strength",
         y = "Incoming interaction strength") +
    theme(text = element_text(size = 10), plot.title = element_text(size = 10, face = "plain"),
          legend.key.height = grid::unit(0.15, "in"))
  plot
})
top <- plot_grid(plotlist = network_plots, ncol = 2,
                 labels = c("(A)", "(B)", "(C)", "(D)"), label_size = 12)
bottom <- plot_grid(plotlist = role_plots, nrow = 1,
                    labels = c("(E)", "(F)", "(G)"), label_size = 12)
figure <- plot_grid(top, bottom, ncol = 1, rel_heights = c(2, 1))
ggsave(file.path(output_dir, "cell_communication.pdf"), figure,
       width = figure_width, height = figure_height, device = cairo_pdf)
ggsave(file.path(output_dir, "cell_communication.png"), figure,
       width = figure_width, height = figure_height, dpi = 300)
