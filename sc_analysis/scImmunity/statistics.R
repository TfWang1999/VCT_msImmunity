# Clinical groups and statistical analysis for scImmunity.
sc_immunity_clinical <- function(config, scores) {
  clinical <- list(
    CRC = sc_read(file.path(config$sc_immunity_input, config$cohorts$CRC$folder, "plot_data_clean.RData")),
    NSCLC_1 = sc_read(file.path(config$sc_immunity_input, config$cohorts$NSCLC_1$folder, "plot_data.RData")),
    NSCLC_2 = read.csv(file.path(config$sc_immunity_input, config$cohorts$NSCLC_2$folder, "NSCLC_2_clinical.csv")),
    BCC = read.csv(file.path(config$sc_immunity_input, config$cohorts$BCC$folder, "BCC_clinical.csv"))
  )
  key <- paste(scores$Cohort, scores$Sample, sep = "|")
  stopifnot(!anyDuplicated(key))
  for (cohort in names(clinical)) {
    meta <- clinical[[cohort]]
    ids <- if (cohort %in% c("CRC", "NSCLC_1")) meta$Patient else if (cohort == "NSCLC_2") meta$sample else meta$patient
    matched <- match(paste(cohort, ids, sep = "|"), key)
    stopifnot(!anyNA(matched))
    for (field in c("Index_list", "Index_gene", "Index_type", "Index")) meta[[field]] <- scores[[field]][matched]
    clinical[[cohort]] <- meta
  }
  clinical$NSCLC_1$Response <- ifelse(clinical$NSCLC_1$radiological_response %in% c("CR", "PR"), "Responder", "Non-responder")
  clinical$NSCLC_1$Pathological_group <- ifelse(clinical$NSCLC_1$pathological_response %in% "pCR", "pCR",
    ifelse(clinical$NSCLC_1$pathological_response %in% c("MPR", "non-MPR"), "non-pCR", NA_character_))
  stopifnot(identical(as.character(clinical$CRC$Response_Status), as.character(clinical$CRC$Response)))
  clinical
}

sc_immunity_statistics <- function(clinical) {
  fields <- c(Full = "Index_list", Proportion = "Index_type", Expression = "Index_gene")
  rows <- list()
  compare <- function(cohort, label, groups, first = "", second = "", panel = "",
                      mask = rep(TRUE, length(groups)), omnibus = FALSE,
                      directional = TRUE, variants = names(fields)) {
    meta <- clinical[[cohort]]
    for (variant in variants) {
      values <- meta[[fields[[variant]]]]
      keep <- !is.na(mask) & mask & !is.na(groups) & is.finite(values)
      one <- NA_real_
      if (omnibus) {
        two <- kruskal.test(values[keep], factor(groups[keep]))$p.value
        n1 <- sum(keep); n2 <- NA_integer_
        group1 <- paste(sort(unique(groups[keep])), collapse = "/")
        test <- "Kruskal-Wallis"
      } else {
        x <- values[keep & groups %in% first]
        y <- values[keep & groups %in% second]
        stopifnot(length(x) > 0L, length(y) > 0L)
        two <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided")$p.value)
        if (directional) one <- suppressWarnings(wilcox.test(x, y, alternative = "greater")$p.value)
        n1 <- length(x); n2 <- length(y); group1 <- first
        test <- if (cohort %in% c("NSCLC_2", "BCC")) "Wilcoxon one-sided: Group1 > Group2" else "Wilcoxon two-sided"
      }
      rows[[length(rows) + 1L]] <<- data.frame(
        Model = "scImmunity", Dataset = cohort, Comparison = label, Score = variant,
        Clinical_plot = if (variant == "Full") panel else "", Group1 = group1, Group2 = second,
        N1 = n1, N2 = n2, Original_test = test,
        P_value = if (!omnibus && cohort %in% c("NSCLC_2", "BCC")) one else two,
        P_two_sided = two, P_one_sided = one,
        One_sided_direction = if (!omnibus && directional) paste(first, ">", second) else "")
    }
  }
  c <- clinical$CRC
  compare("CRC", "R vs NR", c$Response, "Responder", "Non-responder", "crc_response")
  compare("CRC", "CR vs SD", c$Outcome, "CR", "SD", "crc_outcome")
  compare("CRC", "PR vs SD", c$Outcome, "PR", "SD", "crc_outcome")
  compare("CRC", "TRG (overall)", c$TRG, omnibus = TRUE)
  compare("CRC", "TRG0 vs TRG3", c$TRG, "TRG0", "TRG3", "crc_trg", variants = "Full")
  compare("CRC", "TRG2 vs TRG3", c$TRG, "TRG2", "TRG3", "crc_trg", variants = "Full")
  compare("CRC", "R vs NR (dMMR/MSI)", c$Response, "Responder", "Non-responder",
          mask = c$MMR_Group == "dMMR" | c$MSI_Group == "MSI")
  c <- clinical$NSCLC_1
  compare("NSCLC_1", "R vs NR", c$Response, "Responder", "Non-responder", "nsclc_radiological_response")
  compare("NSCLC_1", "pCR vs non-pCR (MPR + non-MPR)", c$Pathological_group, "pCR", "non-pCR", "nsclc_pathological_response")
  compare("NSCLC_1", "LUAD vs LUSC", c$cancer_type, "LUAD", "LUSC", "nsclc_histology", directional = FALSE)
  stage <- rep(NA_character_, nrow(c))
  stage[c$pre_treatment_staging %in% c("IA2", "IB", "IIA", "IIB")] <- "Stage I-II"
  stage[c$pre_treatment_staging %in% c("III", "IIIA", "IIIB", "IIIC", "IVA", "IVB")] <- "Stage III-IV"
  compare("NSCLC_1", "Stage I-II vs III-IV", stage, "Stage I-II", "Stage III-IV")
  for (cohort in c("NSCLC_2", "BCC")) compare(cohort, "R vs NR", clinical[[cohort]]$response, "Yes", "No")
  statistics <- do.call(rbind, rows)
  stopifnot(nrow(statistics) == 35L, sum(nzchar(statistics$Clinical_plot)) == 8L)
  statistics
}

sc_run_immunity_statistics <- function(config, scores) {
  clinical <- sc_immunity_clinical(config, scores)
  statistics <- sc_immunity_statistics(clinical)
  out <- config$sc_immunity_output
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  sc_write_csv(statistics, file.path(out, "clinical_statistics.csv"))
  invisible(list(clinical = clinical, statistics = statistics))
}
