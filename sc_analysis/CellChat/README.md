# CellChat 通讯比较

在本目录依次运行 `source("calculate_cellchat.R")` 和 `source("plot_cellchat.R")`。计算读取 `input_data/` 中 P21 和 P12 两个 CellChat 对象，比较信息流和差异配体–受体对。`infer_cellchat.R` 给出从已注释 CRC RNA 对象推断通讯网络的上游步骤。

`results/` 保存比较结果、对应数据表和论文图。绘图脚本从保存的结果重画通路信息流、网络和弦图。图 8 的定稿组图为 `Cellchat.pdf`，补充图 S6 为 `Sup_cellchat.pdf`，两张组图由对应面板排版得到。

