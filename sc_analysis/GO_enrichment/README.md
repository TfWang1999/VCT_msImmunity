# CRC 差异表达与 GO 富集

在本目录依次运行 `source("calculate_go.R")` 和 `source("plot_go.R")`。计算读取 `input_data/CRC_Treatment_I.rds` 和 `../scImmunity/results/scImmunity_scores.csv`，选取 scImmunity 最高和最低各 5 名 CRC 患者。

`results/GO_results.rds` 保存差异表达和富集结果，`GO_input_groups.csv` 和 `go_enrichment.csv` 分别记录分组和绘图数据。`go_enrichment.pdf` 与 PNG 对应图 8A。

