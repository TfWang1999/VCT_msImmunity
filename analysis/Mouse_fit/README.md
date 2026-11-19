# 小鼠拟合

`run_mouse_analysis` 使用 `input_data/final_parameters.mat` 对 LLC、MC38、A375、CT26 重算分组肿瘤轨迹和拟合指标。数值写入 `results/`。`plot_mouse_fit` 读取这些结果重画图 4。脚本顶部的 `save_mouse_fit` 设为 `true` 时导出 PDF 和 PNG。

`optimize_mouse_parameters` 对各数据集内的对照组和治疗组联合拟合，结果写入 `results/refit_*`。输入及拟合方法见 [METHODS.md](METHODS.md)。

