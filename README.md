# VCT-msImmunity

`analysis/` 包含产生论文模型结果的 MATLAB 代码。进入相应目录运行计算和绘图入口。已有 `results/` 数据可直接用于绘图，`input_data/` 保存计算输入。

| 目录 | 结果 | 计算入口 | 绘图入口 |
| --- | --- | --- | --- |
| `CRC_fit` | 图 3，8 位 CRC 患者拟合 | `run_crc_analysis` | `plot_crc_fit` |
| `Mouse_fit` | 图 4，四个小鼠数据集拟合 | `run_mouse_analysis` | `plot_mouse_fit` |
| `VCT_analysis` | 图 5、6、S4、S5，虚拟患者轨迹、评分与疗效 | `run_virtual_cohort`、`run_vct_analysis` | `plot_tumor_response`、`plot_score_performance` |
| `LHS_PRCC` | 图 S7，末次肿瘤体积敏感性 | `calculate_prcc` | `plot_prcc_bars` |
| `Robustness_analysis` | 图 S8，丰度近似与 AUC 置信区间 | `calculate_abundance_robustness`、`calculate_auc_bootstrap` | `plot_robustness` |
| `src` | 共用的模型、拟合和工具函数 | | |

MATLAB 使用 R2025a，统计计算需要 Statistics and Machine Learning Toolbox。重新优化 CRC 或小鼠参数可运行 `optimize_crc_parameters` 或 `optimize_mouse_parameters`。拟合使用 Optimization Toolbox 中的 `lsqnonlin`，也支持已有的 `fminsearch` 实现。重积分支持 Parallel Computing Toolbox。

图 7、8 和 S6 见 [单细胞分析](sc_analysis/README.md)，图 S2、S3 见 [肝癌 TIME 分型与通讯分析](Liver_cancer/README.md)。

## 单独提供的输入

从 Zenodo 下载 [processed_single_cell_inputs.zip](https://zenodo.org/records/23828850/files/processed_single_cell_inputs.zip?download=1)（约 1.59 GB），解压到仓库根目录（与本 README 同级）并保留压缩包内的目录结构。 [Zenodo 数据记录](https://doi.org/10.5281/zenodo.23828850)

| 文件 | 放置路径（相对于仓库根目录） |
| --- | --- |
| `CRC_Treatment_I.rds` | `sc_analysis/GO_enrichment/input_data/CRC_Treatment_I.rds` |
| `seurat_all.rds` | `Liver_cancer/input_data/seurat_all.rds` |
| `Sce_filter_new.rds` | `Liver_cancer/input_data/Sce_filter_new.rds` |

CRC 拟合需要去标识化临床观测 `analysis/CRC_fit/input_data/patient_all.mat`，可向通讯作者申请。图 3 的 `results/crc_results.mat` 随 CRC 计算生成，其中含临床观测值。
