# Run from this directory with R 4.4.1.
library(Seurat)
library(GSVA)

input_dir <- "../input_data"
output_dir <- "results"
required <- file.path(input_dir, c("seurat_all.rds", "Sce_filter_new.rds"))
if (!all(file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))

# Heatmap: aggregate the original, unfiltered object by sample, as in the original script.
sce <- readRDS(required[1])
stopifnot("S_ID" %in% colnames(sce[[]]))
signatures <- read.csv("../TME_type.csv", fileEncoding = "UTF-8-BOM")
subtypes <- c("IA", "ISM", "ISS", "IE", "IR")
genes <- rownames(sce)
sets <- list(all_related = genes[toupper(genes) %in% signatures$gene])
for (type in subtypes) {
  sets[[type]] <- genes[toupper(genes) %in% signatures$gene[signatures$TIME_subtype == paste0("TIME-", type)]]
}
stopifnot(all(lengths(sets) > 0))
bulk <- AggregateExpression(sce, assays = "RNA", group.by = "S_ID", return.seurat = TRUE)
counts <- GetAssayData(bulk, assay = "RNA", layer = "counts")
scores <- t(gsva(ssgseaParam(exprData = counts, geneSets = sets, normalize = TRUE), verbose = FALSE))
scaled <- t(scale(t(scores[, subtypes, drop = FALSE])))
stopifnot(all(is.finite(scaled)))
assignment <- subtypes[max.col(scaled, ties.method = "first")]
sample_order <- unlist(lapply(subtypes, function(type) {
  ids <- which(assignment == type)
  ids[order(scaled[ids, type])]
}), use.names = FALSE)
heatmap <- t(scaled[sample_order, , drop = FALSE])
samples <- data.frame(sample = rownames(scaled), subtype = assignment)

# Composition: use the saved filtered object and its existing cell annotations/subtypes.
filtered <- readRDS(required[2])
meta <- filtered[[]]
stopifnot(all(c("S_ID", "singleR_anno_new", "TME_subtype") %in% names(meta)))
cell_types <- c("B cells", "CAFs", "TECs", "TAMs", "CTL", "Th", "Treg", "Cancer cells")
stopifnot(!anyNA(meta[, c("S_ID", "singleR_anno_new", "TME_subtype")]),
          all(meta$singleR_anno_new %in% cell_types), all(meta$TME_subtype %in% subtypes))
composition <- function(group) {
  result <- as.data.frame(prop.table(table(meta$singleR_anno_new, group), margin = 2))
  names(result) <- c("cell_type", "group", "fraction")
  result
}
by_sample <- composition(meta$S_ID)
by_subtype <- composition(meta$TME_subtype)
# Keep both sample sets visible for checking the manuscript's sample-count statement.
sample_check <- merge(samples, unique(meta[, c("S_ID", "TME_subtype")]),
                      by.x = "sample", by.y = "S_ID", all = TRUE)
names(sample_check)[2:3] <- c("heatmap_subtype", "saved_composition_subtype")
dir.create(output_dir, showWarnings = FALSE)
write.csv(scores, file.path(output_dir, "ssgsea_scores.csv"))
write.csv(heatmap, file.path(output_dir, "heatmap_values.csv"))
write.csv(sample_check, file.path(output_dir, "sample_membership.csv"), row.names = FALSE)
write.csv(by_sample, file.path(output_dir, "composition_by_sample.csv"), row.names = FALSE)
write.csv(by_subtype, file.path(output_dir, "composition_by_subtype.csv"), row.names = FALSE)
saveRDS(list(heatmap = heatmap, subtype = assignment[sample_order],
             by_sample = by_sample, by_subtype = by_subtype),
        file.path(output_dir, "time_subtype_data.rds"))
message("Heatmap samples: ", ncol(heatmap), "; composition samples: ", length(unique(meta$S_ID)))
