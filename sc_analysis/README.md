# 单细胞分析

使用 R 4.4.1。安装依赖的方法见 [R 环境说明](../config/sc_environment/README.md)。

| 分析 | 计算 | 绘图 |
| --- | --- | --- |
| scImmunity 与临床应答 | `scImmunity/calculate_scores.R` | `scImmunity/plot_sc_immunity.R` |
| CRC 差异表达与 GO 富集 | `GO_enrichment/calculate_go.R` | `GO_enrichment/plot_go.R` |
| CellChat 通讯比较 | `CellChat/calculate_cellchat.R` | `CellChat/plot_cellchat.R` |

在 `sc_analysis/` 运行对应脚本。例如：

```r
source("scImmunity/calculate_scores.R")
source("scImmunity/plot_sc_immunity.R")
```

各分析的输入在 `input_data/`，数值和图在 `results/`。GO 计算使用 scImmunity 的 CRC 评分。已有计算结果可直接用于绘图。评分公式、数据来源和统计见 [METHODS.md](METHODS.md)。

