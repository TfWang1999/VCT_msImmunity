# Windows setup: R 4.4.1 and Rtools44 must already be installed.
# R: https://cran.r-project.org/bin/windows/base/old/4.4.1/
# Rtools: https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html
# In the sc_analysis RStudio project, run:
# source("../config/sc_environment/install_sc_environment.R")
# Or use R 4.4.1 Rscript --vanilla with this script's absolute path.
# Existing packages with the correct versions are left in place.

install_sc_environment <- function() {
  if (.Platform$OS.type != "windows" || getRversion() != "4.4.1") {
    stop("Use Windows R 4.4.1 to install this environment.")
  }

  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) {
    script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
    project <- file.path(dirname(dirname(dirname(script_file))), "sc_analysis")
  } else {
    candidates <- c(getwd(), file.path(getwd(), "sc_analysis"),
                    file.path(getwd(), "VCT_msImmunity_release", "sc_analysis"))
    candidates <- candidates[file.exists(file.path(candidates, "renv.lock"))]
    if (!length(candidates)) stop("Open the sc_analysis project before sourcing this script.")
    project <- candidates[[1]]
  }
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  config <- file.path(dirname(project), "config", "sc_environment")
  cache <- file.path(dirname(dirname(project)), "R_env_cache")
  Sys.setenv(
    LANG = "en_US.UTF-8", LC_ALL = "en_US.UTF-8",
    RENV_PATHS_ROOT = file.path(cache, "renv"),
    RENV_PATHS_CACHE = file.path(cache, "renv", "cache"),
    RENV_PATHS_SANDBOX = file.path(cache, "renv", "sandbox"),
    RENV_CONFIG_AUTO_SNAPSHOT = "FALSE",
    R_MAKEVARS_USER = file.path(config, "Makevars.sc_compat.win")
  )
  options(timeout = max(600, getOption("timeout")),
          repos = c(CRAN = "https://cloud.r-project.org"), install.opts = "--preclean")
  sources <- read.csv(file.path(config, "package_sources.csv"), stringsAsFactors = FALSE)

  # Downloads are temporary; the project library and renv cache are persistent.
  staging <- tempfile("sc-install-")
  dir.create(staging)
  staging <- normalizePath(staging, winslash = "/", mustWork = TRUE)
  temporary_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  stopifnot(startsWith(tolower(staging), paste0(tolower(temporary_root), "/")))
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  download_package <- function(row) {
    extension <- if (row$Type == "binary") ".zip" else ".tar.gz"
    archive <- file.path(staging, paste0(row$Package, "_", row$Version, extension))
    urls <- row$URL
    if (row$Type == "source" && startsWith(row$URL, "https://cloud.r-project.org/")) {
      urls <- unique(c(urls, paste0("https://cloud.r-project.org/src/contrib/Archive/",
                                  row$Package, "/", row$Package, "_", row$Version, ".tar.gz")))
    }
    downloaded <- FALSE
    for (url in urls) {
      downloaded <- tryCatch(
        isTRUE(utils::download.file(url, archive, mode = "wb", method = "libcurl") == 0),
        error = function(error) FALSE
      )
      if (downloaded) break
    }
    if (!downloaded) stop("Download failed for ", row$Package, "; rerun to resume installation.")
    if (nzchar(row$MD5) && !identical(unname(tools::md5sum(archive)), row$MD5)) {
      stop("Download checksum mismatch: ", row$Package)
    }
    archive
  }

  bootstrap <- file.path(cache, "bootstrap", "4.4")
  dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(bootstrap, .libPaths()))
  renv_version <- suppressWarnings(utils::packageDescription("renv", fields = "Version"))
  if (!identical(renv_version, "1.3.0")) {
    if ("renv" %in% loadedNamespaces()) stop("Restart R before installing renv 1.3.0.")
    archive <- download_package(sources[sources$Package == "renv", ])
    utils::install.packages(archive, repos = NULL, type = "source", lib = bootstrap)
    unlink(archive)
  }
  if (!requireNamespace("renv", quietly = TRUE) || packageVersion("renv") != "1.3.0") {
    stop("renv 1.3.0 installation failed.")
  }

  renv::load(project = project)
  project_library <- renv::paths$library(project = project)
  lock <- renv::lockfile_read(file.path(project, "renv.lock"))
  records <- lock$Packages
  installed_version <- function(name) {
    files <- file.path(.libPaths(), name, "DESCRIPTION")
    files <- files[file.exists(files)]
    if (!length(files)) return(NA_character_)
    read.dcf(files[[1]], fields = "Version")[[1]]
  }
  needs_install <- function(name) {
    if (!identical(installed_version(name), records[[name]]$Version)) return(TRUE)
    record <- records[[name]]
    if (identical(record$Source, "GitHub")) {
      descriptions <- file.path(.libPaths(), name, "DESCRIPTION")
      descriptions <- descriptions[file.exists(descriptions)]
      description <- read.dcf(descriptions[[1]])
      return(!"RemoteSha" %in% colnames(description) ||
               !identical(unname(description[1, "RemoteSha"]), record$RemoteSha))
    }
    FALSE
  }
  stopifnot(!anyDuplicated(sources$Package),
            all(sources$Package %in% names(records)),
            all(sources$Version == vapply(records[sources$Package], `[[`, character(1), "Version")))

  # Order dependencies using the lockfile, without the original server export.
  ordered <- visited <- character()
  visit <- function(name) {
    if (name %in% visited) return(invisible(NULL))
    visited <<- c(visited, name)
    fields <- unlist(records[[name]][c("Depends", "Imports", "LinkingTo")], use.names = FALSE)
    dependencies <- trimws(sub("[[:space:]]*\\(.*$", "", fields))
    for (dependency in intersect(dependencies, names(records))) visit(dependency)
    ordered <<- c(ordered, name)
  }
  for (name in names(records)) visit(name)
  pending <- ordered[vapply(ordered, needs_install, logical(1))]
  unavailable <- setdiff(pending, sources$Package)
  if (length(unavailable)) {
    stop("R 4.4.1 bundled package versions differ: ", paste(unavailable, collapse = ", "),
         ". Check the R installation.")
  }

  if (any(sources$Type[match(pending, sources$Package)] == "source")) {
    rtools <- Sys.getenv("RTOOLS44_HOME", unset = "C:/rtools44")
    Sys.setenv(PATH = paste(file.path(rtools, "usr", "bin"),
                           file.path(rtools, "x86_64-w64-mingw32.static.posix", "bin"),
                           Sys.getenv("PATH"), sep = .Platform$path.sep))
    if (!nzchar(Sys.which("make"))) stop("Install Rtools44 before compiling source packages.")
  }

  for (name in pending) {
    row <- sources[sources$Package == name, ]
    cat("INSTALL:", name, row$Version, "\n")
    archive <- download_package(row)
    utils::install.packages(archive, repos = NULL, lib = project_library,
                            type = if (row$Type == "binary") "win.binary" else "source",
                            INSTALL_opts = "--preclean")
    if (!identical(installed_version(name), records[[name]]$Version)) stop("Installation failed: ", name)
    if (identical(records[[name]]$Source, "GitHub")) {
      # Preserve the commit identity when installing a GitHub tarball with base R.
      description <- file.path(project_library, name, "DESCRIPTION")
      metadata <- as.list(read.dcf(description)[1, ])
      remote <- records[[name]]
      for (field in c("RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteRef", "RemoteSha")) {
        if (!is.null(remote[[field]])) metadata[[field]] <- remote[[field]]
      }
      metadata$GithubSHA1 <- remote$RemoteSha
      write.dcf(as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE), description)
      # packageDescription() reads the installed metadata cache as well.
      metadata_file <- file.path(project_library, name, "Meta", "package.rds")
      package_metadata <- readRDS(metadata_file)
      package_metadata$DESCRIPTION <- unlist(metadata, use.names = TRUE)
      saveRDS(package_metadata, metadata_file)
    }
    if (needs_install(name)) stop("Installed package differs from its locked source: ", name)
    unlink(archive)
  }

  actual <- vapply(names(records), installed_version, character(1))
  expected <- vapply(records, `[[`, character(1), "Version")
  if (any(is.na(actual) | actual != expected)) stop("Package version verification failed.")
  for (name in c("Seurat", "CellChat", "ggplot2", "patchwork", "ComplexHeatmap")) {
    if (!requireNamespace(name, quietly = TRUE)) stop("Package cannot load: ", name)
  }
  cat("Verified", length(records), "locked packages under", R.version.string, "\n")
  renv::status(project = project)
  invisible(project_library)
}

install_sc_environment()
