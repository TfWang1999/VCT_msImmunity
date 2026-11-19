# LHS 与 PRCC

`calculate_prcc` 使用 `input_data/` 中的 2,000 组参数和五种情境下的终末肿瘤体积，计算 PRCC 和置信区间。运行 `plot_prcc_bars` 绘制各情境绝对 PRCC 前五名的并集（图 S7），或运行 `plot_prcc_bars(true,'off')` 导出。数值和图保存在 `results/`。

用 `calculate_prcc(true)` 可先重积分同一设计，`calculate_prcc(true,8)` 使用 8 个并行进程。详见 [METHODS.md](METHODS.md)。

