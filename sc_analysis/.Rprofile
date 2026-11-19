if (getRversion() != "4.4.1") {
  stop("This project uses R 4.4.1. Select R 4.4.1 in RStudio or use its Rscript executable.")
}
local({
  project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  cache <- file.path(dirname(dirname(project)), "R_env_cache", "renv")
  Sys.setenv(
    RENV_PATHS_ROOT = cache,
    RENV_PATHS_CACHE = file.path(cache, "cache"),
    RENV_PATHS_SANDBOX = file.path(cache, "sandbox"),
    R_MAKEVARS_USER = file.path(dirname(project), "config", "sc_environment", "Makevars.sc_compat.win")
  )
})
source("renv/activate.R")
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(install.opts = "--preclean")
