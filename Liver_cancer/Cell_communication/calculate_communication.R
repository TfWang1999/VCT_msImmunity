# Run from this directory. Reuse saved inference; do not rerun communication probabilities.
library(CellChat)

input_dir <- "../input_data"
output_dir <- "results"
# Secreted Signaling objects used for the manuscript scatter panels.
files <- c("IA-S007" = "CellChat_S007_IA.rds", "ISM-S011" = "CellChat_S011_ISM.rds",
           "IE-S006" = "CellChat_S006_IE.rds", "IR-S004" = "CellChat_S004_IR.rds")
missing <- file.path(input_dir, files)[!file.exists(file.path(input_dir, files))]
if (length(missing)) stop("Missing input: ", paste(missing, collapse = ", "))
objects <- lapply(file.path(input_dir, files), readRDS)
names(objects) <- names(files)
for (name in names(objects)) {
  expected <- sub(".*-", "", name)
  actual <- unique(as.character(objects[[name]]@meta$S_ID))
  if (length(actual) != 1L || is.na(actual) || actual != expected) {
    stop(files[[name]], " contains S_ID = ", paste(actual, collapse = ", "),
         "; expected ", expected, ". Check the saved object before plotting.")
  }
}
networks <- lapply(objects, function(object) {
  list(weight = object@net$weight, count = object@net$count,
       cells = as.numeric(table(object@idents)))
})
# Original custom scatter uses pathway-specific link counts for dot size.
object <- netAnalysis_computeCentrality(objects[["IR-S004"]], slot.name = "netP")
pathways <- c("TNF", "IFN-II", "TGFb")
stopifnot(all(pathways %in% object@netP$pathways))
roles <- do.call(rbind, lapply(pathways, function(pathway) {
  centrality <- object@netP$centr[[pathway]]
  count <- aggregateNet(object, signaling = pathway, return.object = FALSE, remove.isolate = FALSE)$count
  cells <- levels(object@idents)
  stopifnot(identical(rownames(count), cells), identical(colnames(count), cells))
  data.frame(pathway = pathway, cell_type = cells,
             outgoing = centrality$outdeg, incoming = centrality$indeg,
             Count = rowSums(count) + colSums(count) - diag(count), row.names = NULL)
}))
size_limits <- range(unlist(lapply(networks, function(net) {
  rowSums(net$count) + colSums(net$count) - diag(net$count)
})))
dir.create(output_dir, showWarnings = FALSE)
write.csv(roles, file.path(output_dir, "signaling_roles_S004.csv"), row.names = FALSE)
for (name in names(networks)) {
  write.csv(networks[[name]]$weight, file.path(output_dir, paste0("interaction_strength_", name, ".csv")))
}
saveRDS(list(networks = networks, roles = roles, size_limits = size_limits,
             cell_types = levels(object@idents)), file.path(output_dir, "communication_data.rds"))
