# Single-cell modules named by analysis content (2026-11-08)

- Renamed the single-cell modules to scImmunity, GO_enrichment and CellChat, with independent calculation and plotting entries.
- Replaced figure-number identifiers and output names with analysis names; original data and result contents are preserved.
- New calculations and plots use each module's reproduced directory; usage READMEs are shortened.

# Archive consolidation (2026-11-05)

- Moved the documentation-only root docs, data and outputs folders to the shared archive outside the release.
- Retained current data-processing notes, fixed-A references and fitting bounds in the corresponding analysis READMEs; moved the model preservation record into analysis/README_figures.md.
- Updated documentation links and removed references to archive folders that are no longer present; code, inputs, results and figures are unchanged.

# Single-cell data and figure paths (2026-11-03)

- Grouped existing scripts, inputs and outputs into Figure7_scImmunity, Figure8_GO and Figure8_CellChat.
- Moved cohort clinical CSV files beside the corresponding df_data inputs; GO reads the single shared scImmunity score table.
- Stored GO RDS results under results and CSV tables under source_data; the existing GO plotting function now owns Figure 8A exports for both callers.
- Preserved the manual Cellchat.pdf composite and all existing analysis settings; plotting uses saved GO results and CellChat objects without rerunning enrichment or communication inference.
- Updated entry paths and documentation; the MATLAB modules and R environment are unchanged.

# Model data and figure paths (2026-11-03)

- Moved existing CRC and mouse analysis drivers into their Figure 3/4 folders; moved virtual-cohort and scoring drivers into analysis/VCT_analysis.
- Separated original inputs, calculation MAT files, and exported tables into input_data, results, and source_data beside the corresponding analysis drivers.
- Moved existing model figures and plot-source tables beside their plotting scripts; Figure 5/6/S4/S5 share one set of VCT calculations.
- Removed output-directory overrides from plotting and cohort generation; plot options are now the fourth render_paper_figure argument.
- Preserved model equations, fitted parameters, existing cohort simulations, and single-cell paths.

# Mouse calibration with fixed anti-PD-1 concentration (2026-11-01)

- Fixed treatment A at 200 nM for all four mouse models, supported as an effective concentration scale by published mouse anti-PD-1 pharmacokinetic studies.
- Removed A from mouse optimization; rescaled initial p_2 to retain the original p_2*A products and initial trajectories.
- Added fixed-concentration metadata and explicit fitted/fixed status to mouse result exports; documented the literature and parameter interpretation (now in analysis/Figure4_mouse/README.md).
- Clarified the shared parameter comments to use A in nM and p_2 in nM^-1; retained the CRC and virtual-cohort parameter settings.

# Figure-code organization (2026-10-06)

- Grouped existing main Figure 3B/4/5/6 visualization code into individual folders with Chinese-commented, section-based one-click preview entries.
- Added Supplementary Figure S4/S5 entries by reusing the corresponding original main_quantify scatter plotting blocks.
- Preserved standalone msImmunity scatter and ROC views under the Figure 6 folder, retaining existing function/export names for compatibility.
- Kept default previews export-free and made CRC redraw use saved optimized trajectories without a model rerun; palettes, data, typography and layout are unchanged.
- Updated setup_paths and kept analysis/main.m as a section-based combined entry.
- Retained all CRC/mouse calibration entries, objective functions, optimizers, model code, source data and reference results.
- Added a figure-entry check for clean-path execution, figure/panel counts and no-export/no-refit behavior; the checks were moved to the sibling archive on 2026-10-31.

# Corrected release changes (2026-08-17)

- Removed the dual mapping workflow and retained only the active 76-parameter mapping.
- Corrected virtual-patient columns 9, 12, 15, and 63.
- Replaced independent random draws with randomized Latin-hypercube sampling using the original bounded beta marginals.
- Replaced repeated per-patient random-seed initialization with a stored patient-specific initial tumor fraction.
- Reclassified virtual-patient response using equivalent-diameter RECIST thresholds.
- Regenerated four 1,000-patient cohorts and Figures 5-6.
- Treated all-zero CRC volume triplets as missing placeholders; P02 day 12 is included.
- Recalibrated eight CRC cases and four mouse models with bounded multistart optimization.
- Regenerated Figures 3-4 from optimized parameters.
- Added vector PDF, 600-dpi PNG, CSV source data, and release regression tests.
- Restored the author's original MATLAB visualization settings for CRC panels, mouse validation, the msImmunity response plot, and ROC curves; only the corrected data interface was changed.

The ODE equations, multiscale coupling equations, and the three immunity-score definitions were not changed.
