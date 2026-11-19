# Unchanged filtering, proportions and module means from the processed df_data.
# This file extracts inputs only; the fixed scImmunity formula is in scoring.R.
sc_immunity_patient_features <- function(data, patient, cohort) {
  df <- data$all[[patient]]
  df4 <- data$cd4[[patient]]
  df8 <- data$cd8[[patient]]
  keep <- c("CD4_Th", "Treg", "CD8_em", "CD8_ex")
  if (cohort == "CRC") keep <- c(keep, "TAM")
  if (cohort == "NSCLC_1") keep <- c(keep, "M1", "M2")
  if (cohort == "CRC") {
    final_df <- df[df$Treatment %in% "I" & df$Type %in% keep, , drop = FALSE]
    df4 <- df4[df4$Treatment %in% "I" & df4$Type %in% keep, , drop = FALSE]
    df8 <- df8[df8$Treatment %in% "I" & df8$Type %in% keep, , drop = FALSE]
  } else {
    final_df <- df[df$Type %in% keep, , drop = FALSE]
    df4 <- df4[df4$Type %in% keep, , drop = FALSE]
    df8 <- df8[df8$Type %in% keep, , drop = FALSE]
  }
  proportions <- prop.table(table(final_df$Type))
  Th_pro <- as.numeric(proportions["CD4_Th"])
  Treg_pro <- as.numeric(proportions["Treg"])
  T8em_pro <- as.numeric(proportions["CD8_em"])
  T8ex_pro <- as.numeric(proportions["CD8_ex"])
  mean_score <- function(table, column, type) mean(table[[column]][table$Type %in% type], na.rm = TRUE)
  Act_4 <- mean_score(df4, "Activation/Effector function", "CD4_Th")
  Costim_4 <- mean_score(df4, "Costimulatory molecules", "CD4_Th")
  IFN_4 <- mean_score(df4, "IFN response", "CD4_Th")
  etah <- mean_score(df4, "Cytotoxicity", "CD4_Th")
  Treg_sig <- mean_score(df4, "Treg signature", "Treg")
  Act8 <- mean_score(df8, "Activation:Effector function", "CD8_em")
  eta8 <- mean_score(df8, "Cytotoxicity", "CD8_em")
  IFN_8 <- mean_score(df8, "IFN Response", "CD8_em")
  Exh_ex <- mean_score(df8, "Exhaustion", "CD8_ex")
  eta8_ex <- mean_score(df8, "Cytotoxicity", "CD8_ex")
  S_8 <- (1 + mean_score(df8, "Anti-apoptosis", "CD8_em")) / (1 + mean_score(df8, "Pro-apoptosis", "CD8_em"))
  S_ex <- (1 + mean_score(df8, "Anti-apoptosis", "CD8_ex")) / (1 + mean_score(df8, "Pro-apoptosis", "CD8_ex"))
  S_4 <- (1 + mean_score(df4, "Anti-apoptosis", "CD4_Th")) / (1 + mean_score(df4, "Pro-apoptosis", "CD4_Th"))
  S_Treg <- (1 + mean_score(df4, "Anti-apoptosis", "Treg")) / (1 + mean_score(df4, "Pro-apoptosis", "Treg"))
  M1_pro <- M2_pro <- m1_sig <- m2_sig <- NA_real_
  if (cohort %in% c("CRC", "NSCLC_1")) {
    dfM <- data$tam[[patient]]
    if (cohort == "CRC") {
      sub_proportions <- prop.table(table(final_df$SubCellType))
      m1_types <- c("c61_Mph_FCGR3A", "c62_Mph_S100A8", "c63_Mph_CCL20")
      m2_types <- c("c64_Mph_SPP1", "c65_Mph_APOE")
      value <- function(type) ifelse(is.na(sub_proportions[type]), 0, as.numeric(sub_proportions[type]))
      M1_pro <- value(m1_types[1]) + value(m1_types[2]) + value(m1_types[3])
      M2_pro <- value(m2_types[1]) + value(m2_types[2])
      m1_sig <- mean(dfM$M1_SCORE[dfM$SubCellType %in% m1_types], na.rm = TRUE)
      m2_sig <- mean(dfM$M2_SCORE[dfM$SubCellType %in% m2_types], na.rm = TRUE)
    } else {
      M1_pro <- as.numeric(proportions["M1"])
      M2_pro <- as.numeric(proportions["M2"])
      m1_sig <- mean_score(dfM, "M1_SCORE", "M1")
      m2_sig <- mean_score(dfM, "M2_SCORE", "M2")
    }
  }
  fields <- c("Th_pro", "Treg_pro", "T8em_pro", "T8ex_pro", "M1_pro", "M2_pro",
              "Act_4", "Costim_4", "IFN_4", "etah", "Treg_sig", "S_4", "S_Treg",
              "Act8", "eta8", "IFN_8", "S_8", "Exh_ex", "eta8_ex", "S_ex", "m1_sig", "m2_sig")
  values <- lapply(mget(fields, envir = environment(), inherits = FALSE), unname)
  data.frame(Sample = patient, Cohort = cohort, as.data.frame(values), row.names = NULL)
}

sc_immunity_extract_features <- function(config) {
  tables <- lapply(names(config$cohorts), function(cohort) {
    data <- sc_read_cohort(config, cohort)
    features <- do.call(rbind, lapply(names(data$all), function(patient) sc_immunity_patient_features(data, patient, cohort)))
    message(cohort, ": ", nrow(features), " input feature rows extracted from df_data.")
    features
  })
  features <- do.call(rbind, tables)
  rownames(features) <- NULL
  features
}
