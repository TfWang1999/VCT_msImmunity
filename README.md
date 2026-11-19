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

| 本地路径 | 获取方式 |
| --- | --- |
| `analysis/CRC_fit/input_data/patient_all.mat` | 向通讯作者申请去标识化临床观测 |
| `sc_analysis/GO_enrichment/input_data/CRC_Treatment_I.rds` | 获取已处理的 GSE236581 CRC 表达对象后放到此路径 |
| `Liver_cancer/input_data/seurat_all.rds`、`Liver_cancer/input_data/Sce_filter_new.rds` | 获取已处理的肝癌表达和细胞注释对象后放到对应路径 |

以上文件用于从观测或表达数据重算。图 3 的 `results/crc_results.mat` 随 CRC 计算生成，其中含临床观测值。

