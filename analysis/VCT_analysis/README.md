# 虚拟临床试验

四种方案共用 `input_data/VP_parameter_matrix.mat` 中的 1,000 名患者及初始肿瘤比例。四个 `Therapy_*.mat` 保存给药方案。[METHODS.md](METHODS.md) 说明抽样、模拟和评分方法。

`results/` 中的 `vct_analysis_results.mat` 和 CSV 提供论文用的患者级结果和 AUC。直接运行绘图入口即可从这套结果绘图。`plot_tumor_response(true)` 导出图 5 的各面板。`plot_score_performance("tsImmunity",true)` 导出图 S4，改用 `"molsImmunity"` 可导出图 S5。

```matlab
plot_tumor_response;                     % Fig. 5
plot_score_performance;                  % Fig. 6
plot_score_performance("tsImmunity");    % Fig. S4
plot_score_performance("molsImmunity");  % Fig. S5
```

从保存的设计重算时，在本目录运行：

```matlab
run_virtual_cohort;
run_vct_analysis;
```

第一步向 `results/cohort/` 写入四种治疗的轨迹，第二步从同一目录读取并产生论文分析所用的汇总。可用 `run_virtual_cohort(8)` 并行重积分。

