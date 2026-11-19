# scImmunity 评分与临床比较

在本目录依次运行 `source("calculate_scores.R")` 和 `source("plot_sc_immunity.R")`。计算读取 `input_data/` 中的四队列模块评分及临床信息，生成三种评分和临床统计。NSCLC_1 的正式输入为六类细胞各占保留细胞 >2% 且各至少 20 个的 140 条样本记录。

`results/` 中的 `scImmunity_scores.csv` 也供 GO 分组使用，`clinical_statistics.rds` 和 `clinical_statistics.csv` 分别提供绘图数据和统计表。`sc_immunity_overview.pdf` 对应图 7。

上游细胞注释、模块评分和样本筛选见 `prepare_scores.R`。评分与统计方法见 [METHODS.md](../METHODS.md)。

