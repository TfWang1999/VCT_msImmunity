# Upstream df_data calculations, extracted without changing the original operations.
# Pass the same annotated/filtered Seurat checkpoint used by the original script.
# Marker lists are explicit arguments; no global environment is cleared.

sc_write_df_data <- function(data, cohort, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  df_all <- data$all
  df_all_CD4 <- data$cd4
  df_all_CD8 <- data$cd8
  if (cohort == "CRC") {
    df_all_TAM <- data$tam
    save(df_all, file = file.path(directory, "df_all.RData"))
    save(df_all_CD4, file = file.path(directory, "df_all_CD4.RData"))
    save(df_all_CD8, file = file.path(directory, "df_all_CD8.RData"))
    save(df_all_TAM, file = file.path(directory, "df_all_TAM.RData"))
  } else {
    saveRDS(df_all, file.path(directory, "df_all.rds"))
    save(df_all_CD4, file = file.path(directory, "df_CD4.RData"))
    save(df_all_CD8, file = file.path(directory, "df_CD8.RData"))
    if (!is.null(data$tam)) {
      df_all_TAM <- data$tam
      save(df_all_TAM, file = file.path(directory, "df_TAM.RData"))
    }
    if (cohort == "NSCLC_1") {
      fields <- c("pathological_response_rate", "radiological_response", "pathological_response",
                  "pre_treatment_staging", "age", "targeted_therapy", "gender",
                  "cancer_type", "smoking_history")
      plot_data <- do.call(rbind, lapply(names(df_all), function(id) {
        data.frame(Patient = id, df_all[[id]][1, fields, drop = FALSE], row.names = NULL)
      }))
      names(plot_data)[names(plot_data) == "age"] <- "Age"
      names(plot_data)[names(plot_data) == "gender"] <- "Gender"
      save(plot_data, file = file.path(directory, "plot_data.RData"))
    }
  }
  invisible(directory)
}

sc_annotate_crc <- function(seurat_obj_all) {
  suppressPackageStartupMessages(library(Seurat))
  {
  CD4_Tn <- c("c01_CD4_Tn_CCR7","c02_CD4_Tn_SELL","c03_CD4_Tn_NR4A2")
  CD4_Th <- c("c04_CD4_Tcm_ANXA1","c05_CD4_Tcm_GPR183","c06_CD4_Trm_HSPA1A",
              "c07_CD4_Th17_CTSH","c08_CD4_Tfh_CXCL13_IL6ST","c09_CD4_Th1_CXCL13_HAVCR2",
              "c10_CD4_Temra_GZMB")
  Treg <- c("c11_CD4_Treg_FOXP3","c12_CD4_Treg_KLRB1","c13_CD4_Treg_TNFRSF9")
  CD4_Special <- c("c14_CD4_MT")
  CD8_Tn <- c("c15_CD8_Tn_CCR7", "c16_CD8_Tn_SELL")
  CD8_em <- c("c17_CD8_Tcm_GPR183", "c18_CD8_Tcm_ANXA1","c19_CD8_Tem_CMC1", "c20_CD8_Tem_GZMK",
    "c21_CD8_Trm_XCL1", "c22_CD8_Trm_HSPA1B","c24_CD8_Temra_CX3CR1", "c25_CD8_Temra_TYROBP","c28_CD8_IEL_CD160")
  CD8_ex <- c("c23_CD8_Tex_LAYN")
  MAIT <- c("c26_CD8_MAIT_KLRB1", "c27_CD8_MAIT_SLC4A10")
  T_pro <- c("c29_T_MKI67")
  T_Special <- c("c30_gdT_TRDV2","c31_DNT","c32_DPT")
  ILC <- c("c33_ILC_GZMK","c34_ILC_NR4A2","c35_ILC_FCGR3A","c36_ILC_SYNE2","c37_ILC_SLC4A10","c38_ILC_MKI67")
  B_cell <- c("c39_NaiveB_TCL1A", "c40_NaiveB_IGHD","c41_MemB_CD27", "c42_MemB_GPR183","c43_GCB_LRMP",
             "c44_GCB_MKI67","c45_PlasmaB_TXNDC5", "c46_PlasmaB_IGHA1", "c47_PlasmaB_IGHG1", "c48_PlasmaB_MKI67")
  Mono <- c("c49_Mono_CD14", "c50_Mono_FCGR3A", "c51_Mono_NFKBIZ")
  Mast<- c("c52_Mast_TPSAB1")
  Neu <- c("c53_Neu_CXCR2")
  pDC <- c("c54_pDC_LILRA4", "c55_pDC_GZMB")
  DC <- c("c56_tDC_SIGLEC6","c57_cDC_CLEC9A",
          "c58_cDC_CD1C", "c59_cDC_IL1B", "c60_cDC_LAMP3")
  TAM <- c("c61_Mph_FCGR3A", "c62_Mph_S100A8", "c63_Mph_CCL20", "c64_Mph_SPP1", "c65_Mph_APOE")
  Endo <- c("c66_Endo_LYVE1", "c67_Endo_ACKR1", "c68_Endo_FABP5",
            "c69_Endo_CXCL12", "c70_Endo_COL4A1", "c71_Endo_MKI67")
  Fibro <- c("c72_Fibro_SOX6", "c73_Fibro_ADAMDEC1", "c74_Fibro_C7",
             "c75_Fibro_CCL19", "c76_Fibro_PI16", "c77_Fibro_ACTA2",
             "c78_Fibro_DES", "c79_Fibro_FAP", "c80_Fibro_MKI67")
  Other <- c("c81_Pericyte_RGS5", "c82_SMC_MYH11","c83_Glia_PLP1",
             "c84_Coloncyte_CA2", "c85_Coloncyte_SLC26A3", "c86_Coloncyte_BEST4",
             "c87_Goblet_MUC2", "c88_EEC_PCSK1N", "c89_Tuft_PTGS1")
  Epi <- c("c90_Epi_MKI67", "c91_Epi_Tumor")
  Bcell <- list(
    NaiveB = c("c39_NaiveB_TCL1A", "c40_NaiveB_IGHD"),
    MemB = c("c41_MemB_CD27", "c42_MemB_GPR183"),
    GCB = c("c43_GCB_LRMP", "c44_GCB_MKI67"),
    PlasmaB = c("c45_PlasmaB_TXNDC5", "c46_PlasmaB_IGHA1", "c47_PlasmaB_IGHG1", "c48_PlasmaB_MKI67"))
  type_list <- list( CD4_Tn = CD4_Tn,  CD4_Th = CD4_Th,
    Treg = Treg,CD4_Special = CD4_Special, CD8_Tn = CD8_Tn,CD8_em = CD8_em,
    CD8_ex = CD8_ex,MAIT = MAIT, T_pro = T_pro, T_Special = T_Special,
    ILC = ILC, B_cell = B_cell,Mono = Mono,Mast = Mast,Neu = Neu, pDC = pDC,
    DC = DC,TAM = TAM, Endo = Endo,Fibro = Fibro,Other = Other,Epi = Epi)
  }
  seurat_obj_all$Type <- NA
  for (type_name in names(type_list)) {
    cells <- type_list[[type_name]]
    seurat_obj_all$Type[seurat_obj_all$SubCellType %in% cells] <- type_name
  }
  seurat_obj_all
}

sc_filter_crc <- function(seurat_obj_all) {
  suppressPackageStartupMessages(library(Seurat))
  keep_types <- c("CD4_Tn", "CD4_Th", "Treg",
                  "CD8_Tn", "CD8_em", "CD8_ex",
                  "DC", "TAM", "Epi", "Endo","Fibro")
  seurat_obj <- subset(seurat_obj_all,
                       subset = Type %in% keep_types & Treatment == "I")
  P_list <- unique(seurat_obj$Patient)
  seurat_obj
}

sc_prepare_crc <- function(seurat_obj, CD4_marker, CD8_marker) {
  suppressPackageStartupMessages(library(Seurat))
  norm_minmax <- function(x) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  df <- seurat_obj@meta.data
  df_all <- split(df, f = df$Patient)
  df_all_CD4 <- list()
  df_all_CD8 <- list()
  gene_sets <- CD4_marker
  seu_obj = subset(seurat_obj,
                   subset = Type %in% c("CD4_Th", "Treg"))
  for (module_name in names(gene_sets)) {
      genes <- gene_sets[[module_name]]
      valid_genes <- intersect(genes, rownames(seurat_obj))
       seu_obj <- AddModuleScore(
       object = seu_obj,
        features = list(valid_genes),
        name = paste0(module_name))
        score_col <- paste0(module_name, "1")
        colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd4_sig_names <- names(CD4_marker)
  df[, cd4_sig_names] <- lapply(df[, cd4_sig_names], norm_minmax)
  df_all_CD4 <- split(df, f = df$Patient)
  seu_obj = subset(seurat_obj, subset = Type %in% c("CD8_em", "CD8_ex"))
  gene_sets <- CD8_marker
  for (module_name in names(gene_sets)) {
    genes <- gene_sets[[module_name]]
    valid_genes <- intersect(genes, rownames(seurat_obj))
    seu_obj <- AddModuleScore(
      object = seu_obj,
      features = list(valid_genes),
      name = paste0(module_name))
    score_col <- paste0(module_name, "1")
    colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd8_sig_names <- names(CD8_marker)
  df[, cd8_sig_names] <- lapply(df[, cd8_sig_names], norm_minmax)
  df_all_CD8 <- split(df, f = df$Patient)
  seu_obj <- subset(seurat_obj,  Type %in% c("TAM") )
  M1_genes <- c(
    "CXCL9", "CXCL10", "CXCL11",
    "CD80", "CD86", "HLA-DRA",
    "TNF", "IL1B", "NOS2",
    "FCGR1A"
  )
  M2_genes <- c(
    "IL10", "TGFB1", "ARG1", "VSIG4",
    "VEGFA", "SPP1", "FN1",
    "CD163", "MRC1", "MS4A4A"
  )
  all_genes <- rownames(seu_obj)
  M1_list_clean <- list(intersect(M1_genes, all_genes))
  M2_list_clean <- list(intersect(M2_genes, all_genes))
  M1_set <- M1_list_clean
  M2_set <- M2_list_clean
  survival_genes <- list(c("MCL1", "BCL2A1", "XIAP", "BIRC5", "NFKB1"))
  apoptosis_genes <- list(c("BAX", "BAK1", "CASP3", "CASP7", "FAS"))
  seu_obj <- AddModuleScore(
    object = seu_obj,
    features = survival_genes,
    name = "Anti-apoptosis"
  )
  seu_obj <- AddModuleScore(
    object = seu_obj,
    features = apoptosis_genes,
    name = "Pro-apoptosis"
  )
    seu_obj <- AddModuleScore(JoinLayers(seu_obj),
                                 features = M1_set,
                                 name = "M1_SCORE")
    seu_obj <- AddModuleScore(JoinLayers(seu_obj),
                                 features = M2_set,
                                 name = "M2_SCORE")
    df <- seu_obj@meta.data
    colnames(df)[colnames(df) == "M1_SCORE1"] <- "M1_SCORE"
    colnames(df)[colnames(df) == "M2_SCORE1"] <- "M2_SCORE"
    colnames(df)[colnames(df) == "Anti-apoptosis1"] <- "Anti-apoptosis"
    colnames(df)[colnames(df) == "Pro-apoptosis1"] <- "Pro-apoptosis"
    df$M1_SCORE <- norm_minmax(df$M1_SCORE)
    df$M2_SCORE <- norm_minmax(df$M2_SCORE)
    df$`Anti-apoptosis` <- norm_minmax(df$`Anti-apoptosis`)
    df$`Pro-apoptosis`<- norm_minmax(df$`Pro-apoptosis`)
    df_all_TAM <- split(df, f = df$Patient)
  list(all = df_all, cd4 = df_all_CD4, cd8 = df_all_CD8, tam = df_all_TAM)
}

sc_select_nsclc1_samples <- function(df_all) {
  populations <- c("CD4_Th", "Treg", "CD8_em", "CD8_ex", "M1", "M2")
  names(Filter(function(cells) {
    counts <- table(factor(cells$Type, levels = populations))
    all(counts >= 20L & counts > 0.02 * sum(counts))
  }, df_all))
}

sc_prepare_nsclc1 <- function(seurat_obj, CD4_marker, CD8_marker) {
  suppressPackageStartupMessages(library(Seurat))
  {
    CD4_Th <- c("CD4T_Tem_GZMA","CD4T_Tfh_CXCL13","CD4T_Th1-like_CXCL13","CD4T_Tm_ANXA1","CD4T_Tm_XCL1")
    Treg <- c("CD4T_Treg_CCR8","CD4T_Treg_FOXP3","CD4T_Treg_MKI67")
    CD8_ex <- c("CD8T_terminal_Tex_LAYN","CD8T_Tex_CXCL13")
    CD8_em <- c("CD8T_Tem_GZMK+GZMH+","CD8T_Tem_GZMK+NR4A1+")
    M1 <- c("Mφ_CXCL10","Mφ_CXCL2","Mφ_ISG15","Mφ_FCGR3A")
    M2 <- c("Mφ_FOLR2","Mφ_MARCO","Mφ_VCAN","Mφ_MMP9","Mφ_S100A10")
    type_list <- list( CD4_Th = CD4_Th,Treg = Treg, CD8_em = CD8_em,
                       CD8_ex = CD8_ex, M1=M1,M2=M2)
  }
  seurat_obj$Type <- NA
  for (type_name in names(type_list)) {
    cells <- type_list[[type_name]]
    seurat_obj$Type[seurat_obj$sub_cell_type %in% cells] <- type_name
  }
  seurat_obj <- subset(seurat_obj,
                       subset = Type %in% c("CD4_Th", "Treg","CD8_em","CD8_ex","M1","M2") &
                         radiological_response %in% c("CR", "PR", "PD", "SD") )
  df <- seurat_obj@meta.data
  df_all <- split(df, f = df$sampleID)
  # Score modules over samples with all six populations, then apply the final thresholds.
  df_with_six_populations <- Filter(function(x) length(unique(x$Type)) == 6L, df_all)
  seurat_obj <- subset(seurat_obj,
                       subset = sampleID %in% names(df_with_six_populations) )
  df <- seurat_obj@meta.data
  df_all <- split(df, f = df$sampleID)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  norm_minmax <- function(x) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  df_all_CD4 <- list()
  df_all_CD8 <- list()
  seu_obj = subset(seurat_obj,
                   subset = Type %in% c("CD4_Th", "Treg"))
  gene_sets <- CD4_marker
  for (module_name in names(gene_sets)) {
    genes <- gene_sets[[module_name]]
    valid_genes <- intersect(genes, rownames(seurat_obj))
    seu_obj <- AddModuleScore(
      object = seu_obj,
      features = list(valid_genes),
      name = paste0(module_name))
    score_col <- paste0(module_name, "1")
    colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd4_sig_names <- names(CD4_marker)
  df[, cd4_sig_names] <- lapply(df[, cd4_sig_names], norm_minmax)
  df_all_CD4 <- split(df, f = df$sampleID)
  seu_obj = subset(seurat_obj,
                   subset = Type %in% c("CD8_ex", "CD8_em"))
  gene_sets <- CD8_marker
  for (module_name in names(gene_sets)) {
    genes <- gene_sets[[module_name]]
    valid_genes <- intersect(genes, rownames(seurat_obj))
    seu_obj <- AddModuleScore(
      object = seu_obj,
      features = list(valid_genes),
      name = paste0(module_name))
    score_col <- paste0(module_name, "1")
    colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd8_sig_names <- names(CD8_marker)
  df[, cd8_sig_names] <- lapply(df[, cd8_sig_names], norm_minmax)
  df_all_CD8 <- split(df, f = df$sampleID)
  seu_obj <- subset(seurat_obj,  Type %in% c("M1","M2") )
  M1_genes <- c(
    "CXCL9", "CXCL10", "CXCL11",
    "CD80", "CD86", "HLA-DRA",
    "TNF", "IL1B", "NOS2",
    "FCGR1A"
  )
  M2_genes <- c(
    "IL10", "TGFB1", "ARG1", "VSIG4",
    "VEGFA", "SPP1", "FN1",
    "CD163", "MRC1", "MS4A4A"
  )
  all_genes <- rownames(seu_obj)
  M1_list_clean <- list(intersect(M1_genes, all_genes))
  M2_list_clean <- list(intersect(M2_genes, all_genes))
  M1_set <- M1_list_clean
  M2_set <- M2_list_clean
  seu_obj <- AddModuleScore(JoinLayers(seu_obj),
                            features = M1_set,
                            name = "M1_SCORE")
  seu_obj <- AddModuleScore(JoinLayers(seu_obj),
                            features = M2_set,
                            name = "M2_SCORE")
  df <- seu_obj@meta.data
  colnames(df)[colnames(df) == "M1_SCORE1"] <- "M1_SCORE"
  colnames(df)[colnames(df) == "M2_SCORE1"] <- "M2_SCORE"
  df$M1_SCORE <- norm_minmax(df$M1_SCORE)
  df$M2_SCORE <- norm_minmax(df$M2_SCORE)
  df_all_TAM <- split(df, f = df$sampleID)
  selected <- sc_select_nsclc1_samples(df_all)
  list(all = df_all[selected], cd4 = df_all_CD4[selected],
       cd8 = df_all_CD8[selected], tam = df_all_TAM[selected])
}

sc_prepare_nsclc2 <- function(seurat_obj, CD4_marker, CD8_marker) {
  suppressPackageStartupMessages(library(Seurat))
  {
    CD4_Tn <- c("CD4_C1-Naive")
    CD4_Th <- c("CD4_C2-Tcm","CD4_C3-Tem","CD4_C7-Th1-like")
    Treg <- c("CD4_C8-Treg")
    CD8_em <- c("Non-exhausted")
    CD8_ex <- c("Tex")
    type_list <- list( CD4_Tn = CD4_Tn, CD4_Th = CD4_Th,
                       Treg = Treg, CD8_em = CD8_em,
                       CD8_ex = CD8_ex)
  seurat_obj$Type <- NA
  for (type_name in names(type_list)) {
    cells <- type_list[[type_name]]
    seurat_obj$Type[seurat_obj$cluster %in% cells] <- type_name
  }
  seurat_obj <- subset(seurat_obj,
                       subset = Type %in% c("CD4_Th", "Treg","CD8_em","CD8_ex"))
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  norm_minmax <- function(x) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  df <- seurat_obj@meta.data
  df_all <- split(df, f = df$sample)
  df_all_CD4 <- list()
  df_all_CD8 <- list()
  seu_obj = subset(seurat_obj,
                   subset = Type %in% c("CD4_Th", "Treg"))
  gene_sets <- CD4_marker
  for (module_name in names(gene_sets)) {
    genes <- gene_sets[[module_name]]
    valid_genes <- intersect(genes, rownames(seurat_obj))
    seu_obj <- AddModuleScore(
      object = seu_obj,
      features = list(valid_genes),
      name = paste0(module_name))
    score_col <- paste0(module_name, "1")
    colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd4_sig_names <- names(CD4_marker)
  df[, cd4_sig_names] <- lapply(df[, cd4_sig_names], norm_minmax)
  df_all_CD4 <- split(df, f = df$sample)
  seu_obj = subset(seurat_obj,
                   subset = Type %in% c("CD8_ex", "CD8_em"))
  gene_sets <- CD8_marker
  for (module_name in names(gene_sets)) {
    genes <- gene_sets[[module_name]]
    valid_genes <- intersect(genes, rownames(seurat_obj))
    seu_obj <- AddModuleScore(
      object = seu_obj,
      features = list(valid_genes),
      name = paste0(module_name))
    score_col <- paste0(module_name, "1")
    colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
  }
  df <- seu_obj@meta.data
  cd8_sig_names <- names(CD8_marker)
  df[, cd8_sig_names] <- lapply(df[, cd8_sig_names], norm_minmax)
  df_all_CD8 <- split(df, f = df$sample)
  P_list <- names(df_all_CD4)
  }
  list(all = df_all, cd4 = df_all_CD4, cd8 = df_all_CD8)
}

sc_prepare_bcc <- function(seurat_obj_all, CD4_marker, CD8_marker) {
  suppressPackageStartupMessages(library(Seurat))
    CD4_Th <- c("CD4_T_cells")
    Treg <- c("Tregs")
    CD8_em <- c("CD8_act_T_cells","CD8_mem_T_cells")
    CD8_ex <- c("CD8_ex_T_cells")
    M1 <- c("M1-like TAM")
    M2 <- c("M2-like TAM")
    type_list <- list( CD4_Th = CD4_Th,
                       Treg = Treg, CD8_em = CD8_em,
                       CD8_ex = CD8_ex, M1 = M1, M2 = M2)
    seurat_obj_all$Type <- NA
    for (type_name in names(type_list)) {
      cells <- type_list[[type_name]]
      seurat_obj_all$Type[seurat_obj_all$final_annotation %in% cells] <- type_name
    }
    keep_types = c("CD4_Th","CD8_em","CD8_ex","Treg")
    seurat_obj <- subset(seurat_obj_all,
                         subset = Type %in% keep_types & treatment == "post")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    norm_minmax <- function(x) {
      return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    }
    df <- seurat_obj@meta.data
    df_all <- split(df, f = df$patient)
    {
    df_all_CD4 <- list()
    df_all_CD8 <- list()
    seu_obj = subset(seurat_obj,
                     subset = Type %in% c("CD4_Th", "Treg"))
    gene_sets <- CD4_marker
    for (module_name in names(gene_sets)) {
      genes <- gene_sets[[module_name]]
      valid_genes <- intersect(genes, rownames(seurat_obj))
      seu_obj <- AddModuleScore(
        object = seu_obj,
        features = list(valid_genes),
        name = paste0(module_name))
      score_col <- paste0(module_name, "1")
      colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
    }
    df <- seu_obj@meta.data
    cd4_sig_names <- names(CD4_marker)
    df[, cd4_sig_names] <- lapply(df[, cd4_sig_names], norm_minmax)
    df_all_CD4 <- split(df, f = df$patient)
    seu_obj = subset(seurat_obj,
                     subset = Type %in% c("CD8_ex", "CD8_em"))
    gene_sets <- CD8_marker
    for (module_name in names(gene_sets)) {
      genes <- gene_sets[[module_name]]
      valid_genes <- intersect(genes, rownames(seurat_obj))
      seu_obj <- AddModuleScore(
        object = seu_obj,
        features = list(valid_genes),
        name = paste0(module_name))
      score_col <- paste0(module_name, "1")
      colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == score_col] <- module_name
    }
    df <- seu_obj@meta.data
    cd8_sig_names <- names(CD8_marker)
    df[, cd8_sig_names] <- lapply(df[, cd8_sig_names], norm_minmax)
    df_all_CD8 <- split(df, f = df$patient)
    P_list <- names(df_all_CD4)
    }
  list(all = df_all, cd4 = df_all_CD4, cd8 = df_all_CD8)
}
