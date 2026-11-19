# Input paths and environment for scImmunity and GO/CellChat.
sc_config <- function(project) {
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  sc_directory <- file.path(project, "scImmunity")
  go_directory <- file.path(project, "GO_enrichment")
  cellchat_directory <- file.path(project, "CellChat")
  list(
    project = project,
    sc_immunity_input = file.path(sc_directory, "input_data"),
    sc_immunity_directory = sc_directory,
    sc_immunity_output = file.path(sc_directory, "results"),
    font = "Arial",
    plot_seed = 123L,
    go = list(
      directory = go_directory,
      output = file.path(go_directory, "results"),
      input = file.path(go_directory, "input_data", "CRC_Treatment_I.rds"),
      group_size = 5L
    ),
    cohorts = list(
      CRC = list(folder = "CRC", h = 0.15, macrophages = TRUE),
      NSCLC_1 = list(folder = "NSCLC_1", h = 0.15, macrophages = TRUE),
      NSCLC_2 = list(folder = "NSCLC_2", h = 0.15, macrophages = FALSE),
      BCC = list(folder = "BCC", h = 0.15, macrophages = FALSE)
    ),
    cellchat = list(
      directory = cellchat_directory,
      output = file.path(cellchat_directory, "results"),
      input = file.path(cellchat_directory, "input_data"),
      samples = c(high = "High_tumor_P21", low = "Low_tumor_P12"),
      # Chord diagrams use these five classes as both sources and targets.
      cell_types = c("CD4T", "CD8T", "TAM", "Treg", "tumor"),
      top_pathways = 50L,
      de = list(thresh.pc = 0.05, thresh.fc = 0.01, thresh.p = 0.05),
      # Upregulated links use 0.7 for both ligand and receptor logFC.
      up_logfc = 0.7,
      up_small_gap = 4,
      down_logfc = -0.05
    )
  )
}

sc_load_environment <- function(project) {
  if (getRversion() != "4.4.1") stop("Use R 4.4.1 for this project.")
  config <- sc_config(project)
  if (.Platform$OS.type == "windows") {
    grDevices::windowsFonts(Arial = grDevices::windowsFont(paste("TT", config$font)))
  }
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(config$project)
  source(".Rprofile", local = new.env(parent = globalenv()))
  for (package in c("dplyr", "tidyr", "ggplot2", "ggpubr", "patchwork", "stringr")) {
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  }
  sys.source(file.path(config$project, "R", "io.R"), envir = globalenv())
  config
}
