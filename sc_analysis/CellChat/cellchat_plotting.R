sc_plot_information_flow <- function(merged, pathways, config) {
  plot <- CellChat::rankNet(merged, signaling = pathways$pathway_name,
                           mode = "comparison", comparison = c(2, 1),
                           stacked = TRUE, font.size = 10) +
    theme(axis.text.y = element_text(size = 8)) +
    scale_fill_manual(values = c(high = "#F46D43", low = "#2166AC"),
                      labels = c(high = "High-scImmunity", low = "Low-scImmunity"), name = NULL) +
    theme(axis.text.y = element_text(size = 7.5, family = config$font),
          axis.text.x = element_text(size = 9, family = config$font),
          legend.position = "top", legend.text = element_text(size = 8))
  plot
}

sc_chord <- function(object, links, title, config, small_gap = 1.5) {
  if (!nrow(links)) stop("No ligand-receptor pairs for ", title)
  types <- config$cellchat$cell_types
  object@idents <- factor(object@idents, levels = types)
  types <- types[types %in% levels(object@idents)]
  sector_group <- NULL
  if (small_gap > 1.5) {
    # Shared ligand/receptor names require an explicit group for gaps to apply.
    for (endpoint in list(c("ligand", "source"), c("receptor", "target"))) {
      gene <- links[[endpoint[[1]]]]
      cell <- links[[endpoint[[2]]]]
      sector <- vapply(seq_along(gene), function(i) {
        suffix <- match(cell[[i]], unique(cell[gene == gene[[i]]])) - 1L
        paste0(gene[[i]], strrep(" ", suffix))
      }, character(1))
      index <- order(match(cell, types), -links$prob)
      sector_group <- c(sector_group, stats::setNames(rep("links", length(index)), sector[index]))
    }
    sector_group <- sector_group[!duplicated(names(sector_group))]
  }
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  CellChat::netVisual_chord_gene(
    object, slot.name = "net", net = links, sources.use = types, targets.use = types,
    lab.cex = 0.8, small.gap = small_gap,
    big.gap = if (is.null(sector_group)) 10 else small_gap, group = sector_group,
    annotationTrackHeight = c(0.03),
    legend.pos.x = 100, show.legend = FALSE, title.name = title
  )
}


sc_plot_cellchat <- function(config, comparison) {
  directory <- config$cellchat$output
  objects <- sc_read_cellchat(config)
  merged <- sc_cellchat_merge(objects)
  flow_plot <- sc_plot_information_flow(merged, comparison$flow$pathways, config)
  sc_save_plot(flow_plot, directory, "information_flow", 3.7, 4.5, config)
  pathways <- list(high = c("MHC-II", "IFN-II"), low = c("MK", "GDF"))
  for (group in names(pathways)) for (pathway in pathways[[group]]) {
    object <- objects[[group]]
    if (!pathway %in% object@netP$pathways) stop("Missing pathway: ", group, "/", pathway)
    draw <- function() CellChat::netVisual_aggregate(object, signaling = pathway,
                                                   layout = "circle", weight.scale = TRUE)
    sc_save_base(draw, directory, paste0("pathway_network_", group, "_", pathway), 4, 4, config)
  }
  for (direction in c("up", "down")) {
    object <- if (direction == "up") objects$high else objects$low
    links <- comparison$differential[[direction]]
    title <- if (direction == "up") "Up-regulated (High)" else "Down-regulated (High)"
    stem <- if (direction == "up") "upregulated_chord" else "downregulated_chord"
    small_gap <- if (direction == "up") config$cellchat$up_small_gap else 1.5
    sc_save_base(function() sc_chord(object, links, title, config, small_gap),
                 directory, stem, 5.5, 5.5, config)
  }
  invisible(flow_plot)
}
