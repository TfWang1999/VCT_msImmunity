# GO enrichment plot from the saved numerical table.
sc_go_panel <- function(plot_data, config = NULL) {
  plot_data <- plot_data %>% arrange(GeneRatio)
  plot_data$ID <- seq_len(nrow(plot_data))
  max_p <- max(abs(plot_data$signed_logP))

  # Reserve extra space for GO terms that wrap to a second line.
  labels <- str_wrap(plot_data$Description, width = 48)
  row_heights <- 1 + 0.8 * str_count(labels, "\n")
  plot_data_y <- cumsum(row_heights) - row_heights / 2
  p_final <- ggplot(plot_data, aes(x = plot_data_y, y = GeneRatio)) +
    geom_point(aes(color = Group), alpha = 0) +
    geom_col(aes(fill = signed_logP), width = 0.6) +
    geom_text(aes(y = 0, label = labels, hjust = ifelse(GeneRatio > 0, 1, 0)),
              nudge_y = ifelse(plot_data$GeneRatio > 0, -0.001, 0.001),
              family = "Arial", size = 2.85, lineheight = 1) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    coord_flip() +
    scale_x_continuous("Description", breaks = NULL, expand = expansion(add = 0.15)) +
    scale_y_continuous("Gene Ratio", labels = abs, limits = c(-0.135, NA),
                       expand = expansion(mult = c(0.015, 0.03))) +
    scale_fill_gradientn(
      name = "-log10(P.adjust)", colors = c("#2166AC", "#A6D1E6", "#FDDBC7", "#F46D43"),
      limits = c(-max_p, max_p), values = c(0, 0.45, 0.55, 1),
      breaks = seq(-20, 20, 5), labels = abs,
      guide = guide_colorbar(order = 1, direction = "horizontal",
                             theme = theme(legend.title.position = "left",
                                           legend.key.width = grid::unit(2.8, "cm"),
                                           legend.key.height = grid::unit(0.26, "cm")))
    ) +
    scale_color_manual(
      name = NULL, values = c("High scImmunity" = "#F46D43", "Low scImmunity" = "#2166AC"),
      labels = c("High scImmunity" = "Enriched in high\nscImmunity group",
                 "Low scImmunity" = "Enriched in low\nscImmunity group"),
      guide = guide_legend(order = 2, nrow = 1, override.aes = list(alpha = 1, size = 3, shape = 15))
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      legend.position = "top", legend.box = "horizontal", legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
      legend.spacing.x = grid::unit(5, "pt"), legend.key.size = grid::unit(10, "pt"),
      legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(3, 5, 3, 3)
    )
  if (!is.null(config)) {
    directory <- config$go$output
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    sc_write_csv(plot_data, file.path(directory, "go_enrichment.csv"))
    sc_save_plot(p_final, directory, "go_enrichment", 5.3, 4.5, config)
  }
  list(plot = p_final, data = plot_data)
}
