sc_read_cellchat <- function(config) {
  suppressPackageStartupMessages(library(CellChat))
  folder <- config$cellchat$input
  lapply(config$cellchat$samples, function(sample) {
    sc_read(file.path(folder, paste0("cellchat_", sample, ".rds")))
  })
}

sc_cellchat_merge <- function(objects) {
  CellChat::mergeCellChat(list(low = objects$low, high = objects$high),
                         add.names = c("low", "high"), cell.prefix = TRUE)
}

sc_information_flow_data <- function(merged, config) {
  net_stats <- CellChat::subsetCommunication(merged, slot.name = "netP")
  all_pathways <- dplyr::bind_rows(net_stats, .id = "sample") %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(total_prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(total_prob)) %>%
    dplyr::slice_head(n = config$cellchat$top_pathways)
  list(pathways = all_pathways, networks = net_stats)
}

sc_cellchat_differential <- function(merged, config) {
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("The original CellChat default requires presto. Run the environment installation script.")
  }
  merged <- CellChat::identifyOverExpressedGenes(
    merged, group.dataset = "datasets", pos.dataset = "high",
    features.name = "high.merged", only.pos = FALSE,
    thresh.pc = config$cellchat$de$thresh.pc,
    thresh.fc = config$cellchat$de$thresh.fc,
    thresh.p = config$cellchat$de$thresh.p,
    group.DE.combined = FALSE
  )
  net <- CellChat::netMappingDEG(merged, features.name = "high.merged", variable.all = TRUE)
  up <- CellChat::subsetCommunication(merged, net = net, datasets = "high",
                                     ligand.logFC = config$cellchat$up_logfc,
                                     receptor.logFC = config$cellchat$up_logfc)
  down <- CellChat::subsetCommunication(merged, net = net, datasets = "low",
                                       ligand.logFC = config$cellchat$down_logfc,
                                       receptor.logFC = config$cellchat$down_logfc)
  list(merged = merged, net = net, up = up, down = down)
}


sc_calculate_cellchat <- function(config) {
  directory <- config$cellchat$output
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  objects <- sc_read_cellchat(config)
  merged <- sc_cellchat_merge(objects)
  flow <- sc_information_flow_data(merged, config)
  differential <- sc_cellchat_differential(merged, config)
  # Plotting needs the differential tables, not another copy of expression data.
  comparison <- list(flow = flow, differential = differential[c("net", "up", "down")])
  saveRDS(comparison, file.path(directory, "cellchat_comparison.rds"))
  sc_write_csv(flow$pathways, file.path(directory, "pathway_order.csv"))
  for (group in names(flow$networks)) {
    sc_write_csv(flow$networks[[group]], file.path(directory, paste0(group, "_network.csv")))
  }
  pathways <- list(high = c("MHC-II", "IFN-II"), low = c("MK", "GDF"))
  for (group in names(pathways)) for (pathway in pathways[[group]]) {
    net <- CellChat::subsetCommunication(objects[[group]], slot.name = "netP", signaling = pathway)
    sc_write_csv(net, file.path(directory, paste0("pathway_network_", group, "_", pathway, ".csv")))
  }
  sc_write_csv(differential$up, file.path(directory, "upregulated_LR_pairs.csv"))
  sc_write_csv(differential$down, file.path(directory, "downregulated_LR_pairs.csv"))
  sc_write_csv(differential$net, file.path(directory, "differential_network.csv"))
  invisible(comparison)
}
