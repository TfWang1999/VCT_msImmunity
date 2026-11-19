# 肝癌 TIME 分型与细胞通讯

使用与 [单细胞分析](../sc_analysis/README.md) 相同的 R 4.4.1 环境。在打开 `sc_analysis.Rproj` 的 R 会话中，切换到对应分析目录后运行：

| 目录 | 计算 | 绘图 |
| --- | --- | --- |
| `TIME_subtypes` | `source("calculate_time_subtypes.R")` | `source("plot_time_subtypes.R")` |
| `Cell_communication` | `source("calculate_communication.R")` | `source("plot_communication.R")` |

表格、绘图数据和图保存在各目录的 `results/`。TIME 分型对应补充图 S2，通讯分析对应补充图 S3。

## 输入

| 文件 | 用途 |
| --- | --- |
| `input_data/seurat_all.rds` | 46 个样本的 ssGSEA 分型及热图 |
| `input_data/Sce_filter_new.rds` | 28 个保留样本的细胞构成 |
| `input_data/CellChat_S007_IA.rds`、`CellChat_S011_ISM.rds`、`CellChat_S006_IE.rds`、`CellChat_S004_IR.rds` | 四个代表性样本的分泌信号网络，S004 用于 TNF、IFN-II 和 TGF-β 散点图 |
| `TME_type.csv` | TIME 分型基因集 |

热图按样本聚合 RNA counts 后计算 ssGSEA，对每类评分做样本间标准化。比例图使用过滤后细胞注释，按样本和亚型汇总。通讯图使用保存的 CellChat 网络和通路特异性互作数。

