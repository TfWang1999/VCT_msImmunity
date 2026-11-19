# Clinical comparisons and score-variant comparison plots.
# Colors, geometry and typography retain the manuscript style.

sc_crc_panels <- function(plot_data_clean) {
  title_select <- c( "scImmunity ","scImmunity ")
  cols_1 <- c("#1b9e77", "#d95f02", "#7570b3")
  cols_2 <- c("#1b9e77", "#d95f02")
  cols_3 <- c("#377eb8", "#e34a33")
  my_theme <- (
    theme_classic(base_size = 12) +
      theme(
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold",vjust = 0.6),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
  crc_response_plot <-  ggplot(plot_data_clean, aes(x =  Response_Status, y = Index)) +
    geom_boxplot(aes(fill = Response_Status), lwd = 0.6, width = 0.5,
                 outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 2.5, shape = 21,
                fill = "white", color = "black", stroke = 0.8, alpha = 0.8) +
    scale_fill_manual(values = cols_3) +
    labs(y = title_select[[2]], x = "Response Status",title = "CRC cohort") +
    my_theme+
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      )
    )+
    theme(legend.position = "none")
  crc_outcome_plot <- ggplot(plot_data_clean, aes(x =  Outcome, y = Index)) +
    geom_boxplot(aes(fill = Outcome), lwd = 0.6, width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 3, shape = 21, fill = "white", stroke = 0.8) +
    geom_smooth(method = "lm", se=FALSE, aes(group=1), color="grey", linetype="dashed") +
    scale_fill_manual(values = cols_1) +
    labs(y = title_select[[2]], x = "Clinical outcome",title = "CRC cohort") +
    my_theme+
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 7)
      )
    )
  crc_trg_plot <- ggplot(plot_data_clean, aes(x =  TRG, y = Index)) +
    geom_boxplot(aes(fill = TRG), lwd = 0.6, width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 3, shape = 21, fill = "white", stroke = 0.8) +
    geom_smooth(method = "lm", se=FALSE, aes(group=1), color="grey", linetype="dashed") +
    scale_fill_manual(values = cols_1) +
    labs(y = title_select[[2]],x = "Tumor Regression Grade",title = "CRC cohort") +
    my_theme+
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 7)
      )
    )
  list(crc_response = crc_response_plot, crc_outcome = crc_outcome_plot, crc_trg = crc_trg_plot)
}

sc_nsclc_panels <- function(plot_data) {
  my_theme <- (
    theme_classic(base_size = 12) +
      theme(
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold",vjust = 0.6),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
  cols_1 <- c("#1b9e77", "#d95f02", "#7570b3","#377eb8")
  cols_2 <- c("#1b9e77", "#d95f02")
  cols_3 <- c("#377eb8","#e34a33")
  plot_data$Response <- factor(plot_data$Response, levels = c("Responder", "Non-responder"))
  nsclc_response_plot <- ggplot(plot_data, aes(x = Response, y = Index_list)) +
    geom_violin(aes(fill = Response),
                trim = FALSE,
                alpha = 0.6,
                lwd = 0.8) +
    geom_boxplot(width = 0.1,
                 fill = "white",
                 outlier.shape = NA,
                 alpha = 0.9) +
    geom_jitter(width = 0.15, size = 2.5, shape = 21, fill = "white", stroke = 0.5,alpha = 0.9) +
    scale_fill_manual(values = c("#377eb8","#e34a33")) +
    labs(title = "NSCLC_1 cohort", x = "Radiological response", y = "scImmunity") +
    my_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 25)
      )
    )
  nsclc_histology_plot <- ggplot(plot_data, aes(x = cancer_type, y = Index_list)) +
    geom_violin(aes(fill = cancer_type),
                trim = FALSE,
                alpha = 0.5,
                lwd = 0.8) +
    geom_boxplot(width = 0.1,
                 fill = "white",
                 outlier.shape = NA,
                 alpha = 1) +
    geom_jitter(width = 0.15, size = 2.5, shape = 21, fill = "white", stroke = 0.5,alpha = 0.9) +
    scale_fill_manual(values = cols_3) +
    labs(title = "NSCLC_1 cohort", x = "Cancer type", y = "scImmunity") +
    my_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 25)
      )
    )
  pathologic <- plot_data[!is.na(plot_data$Pathological_group), ]
  pathologic$Pathological_group <- factor(pathologic$Pathological_group, levels = c("pCR", "non-pCR"))
  nsclc_pathology_plot <- ggplot(pathologic, aes(Pathological_group, Index_list)) +
    geom_violin(aes(fill = Pathological_group), trim = FALSE, alpha = 0.5, scale = "width", linewidth = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 2.5, shape = 21, fill = "white", stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = c("pCR" = "#1b9e77", "non-pCR" = "#d95f02")) +
    scale_x_discrete(labels = c("pCR" = "pCR", "non-pCR" = "non-pCR")) +
    labs(title = "NSCLC_1 cohort", x = "Pathological response", y = "scImmunity") + nsclc_response_plot$theme
  list(nsclc_radiological_response = nsclc_response_plot, nsclc_pathological_response = nsclc_pathology_plot, nsclc_histology = nsclc_histology_plot)
}

sc_comparison_panel <- function(all_data) {
  all_data <- subset(all_data, Comparison != "CR vs PR")
  my_order <- c(
    "R vs NR \n dMMR/MSI group",
    "TRG",
    "PR vs SD",
    "CR vs SD",
    "Stage I-II vs III-IV",
    "LUAD vs LUSC",
    "MPR vs Non-MPR",
    "pCR vs Non-pCR","R vs NR"
  )
  plot_data_final <- all_data %>%
    pivot_longer(cols = c(P_Full, P_Prop, P_Gene), names_to = "Model", values_to = "P_value") %>%
    mutate(
      logP = -log10(P_value),
      Model = case_when(
        Model == "P_Full" ~ "scImmunity",
        Model == "P_Prop" ~ "scImmunity (Only Proportion)",
        Model == "P_Gene" ~ "scImmunity (Only Expression)"
      ),
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      Dataset = factor(Dataset, levels = c("CRC", "NSCLC_1", "NSCLC_2", "BCC")),
      Model = factor(Model, levels = c(
        "scImmunity (Only Expression)",
        "scImmunity (Only Proportion)",
        "scImmunity"
      )),
      Comparison = factor(Comparison, levels = unique(my_order))
    )
  p_grand <- ggplot(plot_data_final, aes(x = Comparison, y = logP, fill = Model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.4) +
    geom_text(aes(label = Significance),
              position = position_dodge(width = 0.8),
              hjust = -0.2, vjust = 0.8, size = 3.5, fontface = "bold") +
    facet_grid(Dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +
    coord_flip() +
    scale_fill_manual(values = c(
      "scImmunity"                   = "#F28E2B",
      "scImmunity (Only Proportion)" = "#4E79A7",
      "scImmunity (Only Expression)"       = "#59A14F"
    ),
    breaks = c("scImmunity", "scImmunity (Only Proportion)", "scImmunity (Only Expression)") ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "", subtitle = "", y = "-log10 (P-value)", x = "") +
    theme_classic() +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text.y.left = element_text(
        angle = 0, family = "Arial", size = 10, face = "bold", color = "black"
      ),
      axis.text.y = element_text(size = 12, color = "black", family = "Arial"),
      axis.text.x = element_text(size = 12, family = "Arial"),
      axis.title = element_text(family = "Arial"),
      legend.position = c(0.7, 0.90),
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold", family = "Arial"),
      legend.key.size = unit(0.5, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5, family = "Arial")
    )
  list(plot = p_grand, data = plot_data_final)
}

# scImmunity clinical comparisons and composite layout.
sc_plot_immunity <- function(config, report, export = TRUE) {
  directory <- config$sc_immunity_output
  clinical <- report$clinical
  statistics <- report$statistics
  clinical_comparisons <- statistics[nzchar(statistics$Clinical_plot), ]
  panels <- c(sc_crc_panels(clinical$CRC), sc_nsclc_panels(clinical$NSCLC_1))
  add_verified_labels <- function(plot, rows) {
    set.seed(config$plot_seed)
    built <- suppressMessages(ggplot_build(plot))
    limits <- range(unlist(lapply(built$data, function(layer) {
      fields <- intersect(c("y", "ymin", "ymax", "lower", "upper"), names(layer))
      unlist(layer[fields])
    })), finite = TRUE)
    span <- diff(limits)
    if (span <= 0) span <- max(abs(limits), 0.1)
    y_positions <- limits[2] + span * seq(0.10, by = 0.14, length.out = nrow(rows))
    format_p <- function(p) format(signif(p, 2), scientific = p < 0.001, trim = TRUE)
    if (nrow(rows) == 1L) {
      plot <- plot + annotate("text", x = 1.5, y = y_positions, label = paste0("p = ", format_p(rows$P_value)), family = config$font, size = 5)
    } else {
      # CR/TRG0 versus SD/TRG3 forms the long bracket; the second forms the short bracket.
      for (i in seq_len(nrow(rows))) {
        x1 <- if (i == 1L) 1 else 2
        x2 <- 3
        y <- y_positions[i]
        plot <- plot + annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.4) +
          annotate("segment", x = c(x1, x2), xend = c(x1, x2), y = y, yend = y - span * 0.025, linewidth = 0.4) +
          annotate("text", x = (x1 + x2) / 2, y = y + span * 0.04, label = format_p(rows$P_value[i]), family = config$font, size = 4)
      }
    }
    plot + scale_y_continuous(expand = expansion(mult = c(0.04, 0.04))) +
      coord_cartesian(ylim = c(limits[1], max(y_positions) + span * 0.10))
  }
  for (panel in names(panels)) panels[[panel]] <- add_verified_labels(panels[[panel]], clinical_comparisons[clinical_comparisons$Clinical_plot == panel, ])

  # Compare computed P values for all three score variants.
  # Pathological response is compared as pCR versus all known non-pCR.
  comparison_statistics <- statistics[!statistics$Comparison %in% c("TRG0 vs TRG3", "TRG2 vs TRG3"), ]
  comparison_table <- comparison_statistics[, c("Dataset", "Comparison", "Score", "P_value")]
  comparison_table$Comparison[comparison_table$Comparison == "TRG (overall)"] <- "TRG"
  comparison_table$Comparison[comparison_table$Comparison == "R vs NR (dMMR/MSI)"] <- "R vs NR \n dMMR/MSI group"
  comparison_table$Comparison[comparison_table$Comparison == "pCR vs non-pCR (MPR + non-MPR)"] <- "pCR vs Non-pCR"
  comparison_table$Score <- c(Full = "P_Full", Proportion = "P_Prop", Expression = "P_Gene")[comparison_table$Score]
  comparison_table <- tidyr::pivot_wider(comparison_table, names_from = Score, values_from = P_value)
  stopifnot(nrow(comparison_table) == 11L, nrow(comparison_statistics) == 33L)
  comparison <- sc_comparison_panel(comparison_table)
  panels$score_comparison <- comparison$plot +
    scale_x_discrete(labels = function(x) ifelse(x == "pCR vs Non-pCR", "pCR vs non-pCR", x)) +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.955),
          legend.justification = c(1, 1), legend.text = element_text(size = 10, face = "bold", family = config$font),
          legend.key.size = grid::unit(0.35, "cm"), legend.margin = margin(0, 0, 0, 0)) +
    guides(fill = guide_legend(ncol = 1))
  left <- (panels$crc_response | panels$crc_outcome | panels$crc_trg) / patchwork::plot_spacer() / (panels$nsclc_radiological_response | panels$nsclc_pathological_response | panels$nsclc_histology) +
    patchwork::plot_layout(heights = c(10, 0.5, 10))
  composite <- (left | panels$score_comparison) + patchwork::plot_layout(widths = c(0.7, 0.3)) +
    patchwork::plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag.position = c(0.08, 0.98), plot.tag = element_text(size = 20, face = "bold", family = config$font))
  if (export) {
    sc_save_plot(composite, directory, "sc_immunity_overview", 17, 9.6, config)
  }
  message("scImmunity exported to ", directory)
  invisible(list(panels = panels, composite = composite, statistics = statistics))
}
