# 评分稳健性

`calculate_abundance_robustness` 从共用虚拟队列计算五个丰度近似常数下的 msImmunity AUC。`calculate_auc_bootstrap` 对患者进行 5,000 次 bootstrap，生成三种评分的 AUC 及 95% 置信区间。

`plot_robustness` 读取 `results/` 中两张结果表绘制图 S8。`plot_robustness(true,'off')` 导出组图。详见 [METHODS.md](METHODS.md)。

