# Run with R 4.4.1, or source this file in the RStudio project.
.sc_entry <- local({
  files <- vapply(rev(sys.frames()), function(frame) if (is.null(frame$ofile)) "" else frame$ofile, character(1))
  files <- files[nzchar(files)]
  if (length(files)) normalizePath(files[[1]], winslash = "/", mustWork = TRUE)
  else {
    argument <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(argument)) stop("Source this script or run it with Rscript.")
    normalizePath(sub("^--file=", "", argument[[1]]), winslash = "/", mustWork = TRUE)
  }
})
.sc_project <- dirname(dirname(.sc_entry))
source(file.path(.sc_project, "config.R"))
config <- sc_load_environment(.sc_project)

source(file.path(config$go$directory, "go_plotting.R"))
plot_data <- readRDS(file.path(config$go$output, "GO_results.rds"))$plot_data
invisible(sc_go_panel(plot_data, config))
