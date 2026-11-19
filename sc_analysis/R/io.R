sc_read <- function(path, object = NULL) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  container <- new.env(parent = baseenv())
  names <- load(path, envir = container)
  if (is.null(object)) {
    if (length(names) != 1L) stop("Select the RData object explicitly: ", path)
    object <- names[[1]]
  }
  if (!exists(object, envir = container, inherits = FALSE)) stop("Missing object: ", object)
  container[[object]]
}

sc_read_cohort <- function(config, cohort) {
  folder <- file.path(config$sc_immunity_input, config$cohorts[[cohort]]$folder)
  files <- if (cohort == "CRC") {
    c(all = "df_all.RData", cd4 = "df_all_CD4.RData", cd8 = "df_all_CD8.RData", tam = "df_all_TAM.RData")
  } else {
    c(all = "df_all.rds", cd4 = "df_CD4.RData", cd8 = "df_CD8.RData")
  }
  if (cohort == "NSCLC_1") files <- c(files, tam = "df_TAM.RData")
  data <- lapply(files, function(file) sc_read(file.path(folder, file)))
  ids <- names(data$all)
  for (component in names(data)) {
    if (!setequal(ids, names(data[[component]]))) stop("Sample IDs differ: ", cohort, " / ", component)
    data[[component]] <- data[[component]][ids]
  }
  data
}

sc_write_csv <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

sc_save_plot <- function(plot, directory, filename, width, height, config) {
  for (extension in c("pdf", "png")) {
    set.seed(config$plot_seed) # Only point jitter and graphics; score calculations do not use this seed.
    path <- file.path(directory, paste0(filename, ".", extension))
    if (extension == "pdf") {
      ggplot2::ggsave(path, plot, width = width, height = height,
                     device = grDevices::cairo_pdf, family = config$font, bg = "white")
    } else {
      ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
    }
  }
}

sc_save_base <- function(draw, directory, filename, width, height, config) {
  for (extension in c("pdf", "png")) {
    path <- file.path(directory, paste0(filename, ".", extension))
    if (extension == "pdf") {
      grDevices::cairo_pdf(path, width = width, height = height, family = config$font)
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = 300, type = "cairo")
    }
    tryCatch({
      set.seed(config$plot_seed)
      draw()
    }, finally = grDevices::dev.off())
  }
}
