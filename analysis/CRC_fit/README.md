# CRC 拟合

`run_crc_analysis` 读取 `input_data/` 中八位患者的临床观测、给药表及最终拟合参数。参数、指标和轨迹写入 `results/`。`plot_crc_fit` 读取这些结果重画图 3。脚本顶部的 `save_crc_fit` 设为 `true` 时导出 PDF 和 PNG。

重拟合可运行 `optimize_crc_parameters`。例如 `optimize_crc_parameters(1,3,120)` 从三个起点拟合 P01，每个起点最多迭代 120 次。优化另存为 `results/refit_*`。

输入、目标函数和拟合参数见 [METHODS.md](METHODS.md)。

