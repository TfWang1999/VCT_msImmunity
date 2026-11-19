# scImmunity scoring with kp = 0.15 and kf = 0.5 for all four cohorts.
sc_immunity_scores <- function(features) {
  kp <- 0.15
  kf <- 0.5
  Hp <- function(x) x / (x + kf)
  Hi <- function(x) kf / (x + kf)
  A <- function(p) p / (p + kp)
  B <- function(p) kp / (p + kp)
  f <- features
  Fh <- Hp(1 + f$etah * Hp(1 + f$Act_4) * Hp(1 + f$Costim_4) * Hp(1 + f$IFN_4)* Hp(1 + f$S_4)) #
  Fr <- Hi(1 + f$Treg_sig * Hp(1+ f$S_Treg)) #
  Fe <- Hp(1 + f$eta8 * Hp(1 + f$Act8) * Hp(1 + f$IFN_8)  * Hp(1 + f$S_8) ) #
  Fx <- Hp(1 + f$eta8_ex / (1 + f$Exh_ex) * Hp(1 + f$S_ex))# 
  F1 <- Hp(1 + f$m1_sig)
  F2 <- Hi(1 + f$m2_sig)
  full <- A(f$Th_pro) * Fh * B(f$Treg_pro) * Fr + A(f$T8em_pro) * Fe + A(f$T8ex_pro) * Fx
  expression <- Fh * Fr + Fe + Fx
  proportion <- A(f$Th_pro) * B(f$Treg_pro) + A(f$T8em_pro) + A(f$T8ex_pro)
  tam <- f$Cohort %in% c("CRC", "NSCLC_1")
  # NSCLC_2/BCC are T-cell-only because no TAM inputs were supplied.
  full[tam] <- full[tam] + (A(f$M1_pro) * F1 * B(f$M2_pro) * F2)[tam]
  expression[tam] <- expression[tam] + (F1 * F2)[tam]
  proportion[tam] <- proportion[tam] + (A(f$M1_pro) * B(f$M2_pro))[tam]
  data.frame(Sample = f$Sample, Cohort = f$Cohort, Index_list = full,
             Index_gene = expression, Index_type = proportion, Index = full)
}

sc_run_immunity_scores <- function(config) {
  directory <- config$sc_immunity_output
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  features <- sc_immunity_extract_features(config)
  stopifnot(setequal(unique(features$Cohort), names(config$cohorts)),
            !anyDuplicated(paste(features$Cohort, features$Sample)))
  scores <- sc_immunity_scores(features)
  sc_write_csv(scores, file.path(directory, "scImmunity_scores.csv"))
  message("scImmunity: ", nrow(scores), " scores saved to ", directory)
  invisible(scores)
}
