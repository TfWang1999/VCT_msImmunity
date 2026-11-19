# CRC differential expression and GO enrichment using the original analysis settings.

sc_go_subset <- function(seurat_obj,
                         high_patients = c("P21", "P26", "P11", "P20", "P09"),
                         low_patients = c("P12", "P02", "P25", "P16", "P17")) {
  suppressPackageStartupMessages(library(Seurat))
  stopifnot(!anyDuplicated(c(high_patients, low_patients)),
            all(c(high_patients, low_patients) %in% as.character(seurat_obj$Patient)))
  seu_subset <- subset(seurat_obj, subset = Patient %in% c(high_patients, low_patients))
  seu_subset$group <- ifelse(seu_subset$Patient %in% high_patients, "High", "Low")
  Idents(seu_subset) <- "group"
  seu_subset
}

sc_compute_go <- function(seu_subset) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  degs <- FindMarkers(seu_subset,
                      ident.1 = "High",
                      ident.2 = "Low",
                      logfc.threshold = 0.25,
                      min.pct = 0.1)
  up_degs <- degs %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.05)
  down_degs <- degs %>% filter(avg_log2FC < -0.25 & p_val_adj < 0.05)
  top_up_genes <- up_degs %>% slice_max(order_by = avg_log2FC, n = 200) %>% rownames()
  ids_up <- bitr(top_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego_high_raw <- enrichGO(gene          = ids_up$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENTREZID',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable      = TRUE)
  ego_high_sim <- simplify(ego_high_raw, cutoff = 0.7, by = "p.adjust", select_fun = min)
  immune_keywords <- "T cell|lymphocyte|leukocyte|antigen|immune|cytokine|interferon|inflammatory|MHC|chemokine|cell-cell adhesion"
  low_keywords <- "stress|negative regulation|myeloid|detoxification|hypoxia|exhaustion|tolerance|suppression|inhibitory|anergy|negative regulation of immune|apoptotic|PD-L1|TGF"
  universal_keywords <- paste0(immune_keywords, "|", low_keywords)
  high_plot_data <- ego_high_sim@result %>%
    filter(sapply(GeneRatio, function(x) eval(parse(text=x))) > 0.02) %>%
  filter(grepl(universal_keywords, Description, ignore.case = TRUE)) %>%
    slice_head(n = 10)
  top_down_genes <- down_degs %>% slice_max(order_by = abs(avg_log2FC), n = 200) %>% rownames()
  ids_down <- bitr(top_down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego_low_raw <- enrichGO(gene          = ids_down$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENTREZID',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          readable      = TRUE)
  ego_low_sim <- simplify(ego_low_raw, cutoff = 0.7, by = "p.adjust", select_fun = min)
  low_plot_data <- ego_low_sim@result %>%
    filter(sapply(GeneRatio, function(x) eval(parse(text=x))) > 0.02) %>%
    filter(grepl(universal_keywords, Description, ignore.case = TRUE)) %>%
    slice_head(n = 10)
  df_high <- sc_go_term_data(ego_high_sim, high_plot_data$Description)
  df_high$Group <- "High scImmunity"
  df_low <- sc_go_term_data(ego_low_sim, low_plot_data$Description)
  df_low$Group <- "Low scImmunity"
  df_low$GeneRatio <- -df_low$GeneRatio
  plot_data <- rbind(df_high, df_low)
  plot_data$logP <- -log10(plot_data$p.adjust)
  plot_data$signed_logP <- ifelse(plot_data$Group == "Low scImmunity", -plot_data$logP, plot_data$logP)
  plot_data$Description <- factor(plot_data$Description,
                                  levels = unique(plot_data$Description[order(plot_data$GeneRatio)]))
  plot_data$label_y_pos <- ifelse(plot_data$Group == "High scImmunity", -0.002, 0.002)
  plot_data$label_hjust <- ifelse(plot_data$Group == "High scImmunity", 1, 0)
  list(degs = degs, high_raw = ego_high_raw, low_raw = ego_low_raw,
       high = ego_high_sim, low = ego_low_sim, plot_data = plot_data)
}

sc_run_go <- function(config, score_file) {
  if (!file.exists(score_file)) stop("Select an existing scImmunity score table.")
  scores <- read.csv(score_file)
  crc <- scores[scores$Cohort == "CRC" & is.finite(scores$Index_list), c("Sample", "Index_list")]
  crc <- crc[order(-crc$Index_list, crc$Sample), ]
  count <- config$go$group_size
  stopifnot(nrow(crc) >= 2L * count, !anyDuplicated(crc$Sample))
  crc$Rank <- seq_len(nrow(crc))
  crc$Group <- "Unselected"
  crc$Group[seq_len(count)] <- "High"
  crc$Group[seq.int(nrow(crc) - count + 1L, nrow(crc))] <- "Low"
  high_patients <- crc$Sample[crc$Group == "High"]
  low_patients <- crc$Sample[crc$Group == "Low"]
  message("GO High: ", paste(high_patients, collapse = ", "))
  message("GO Low: ", paste(low_patients, collapse = ", "))

  suppressPackageStartupMessages(library(Seurat))
  seurat_obj <- sc_read(config$go$input)
  stopifnot(DefaultAssay(seurat_obj) == "RNA",
            "data" %in% Layers(seurat_obj[["RNA"]]))
  if ("Treatment" %in% colnames(seurat_obj@meta.data)) {
    stopifnot(all(as.character(seurat_obj$Treatment) == "I"))
  }
  counts <- table(as.character(seurat_obj$Patient))
  crc$Cell_count <- as.integer(counts[crc$Sample])
  seu_subset <- sc_go_subset(seurat_obj, high_patients, low_patients)
  rm(seurat_obj)
  gc(verbose = FALSE)
  message("GO input: ", ncol(seu_subset), " cells; High = ", sum(seu_subset$group == "High"),
          "; Low = ", sum(seu_subset$group == "Low"))
  result <- sc_compute_go(seu_subset)
  directory <- config$go$output
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  sc_write_csv(crc, file.path(directory, "GO_input_groups.csv"))
  saveRDS(result, file.path(directory, "GO_results.rds"))
  message("GO numerical results saved in ", directory)
  invisible(list(results = result, groups = crc))
}


# Match the numerical data and factor ordering used by enrichplot::dotplot.
sc_go_term_data <- function(enrichment, descriptions) {
  data <- ggplot2::fortify(enrichment, showCategory = descriptions)
  data$x <- data$GeneRatio
  index <- order(data$x, decreasing = TRUE)
  data$Description <- factor(data$Description, levels = rev(unique(data$Description[index])))
  data
}
