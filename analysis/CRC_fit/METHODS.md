## 数据处理与拟合参数

路径以发布版根目录为起点。

`analysis/CRC_fit/input_data/patient_all.mat` contains eight matrices (`P1`-`P8`). Columns are study day, radiotherapy dose, capecitabine, oxaliplatin, anti-PD-1, mean tumor volume, and two interval limits.

A row is a valid observation when the three volume fields are finite and not all zero. All-zero triplets are treated as missing-value placeholders. Thus P01 excludes its terminal zero row, while P02 includes the nonzero day-12 measurement and excludes its terminal zero row. Interval limits are normalized with `min`/`max` before fitting.

The treatment matrix contains only the four columns consumed by `run_simulation`. A reported study day `d` selects trajectory element `d`, because element 1 represents study day 1 at simulation time 0.

拟合参数及边界以 `analysis/src/fitting/fit_crc_case.m` 和 `analysis/src/config/crc_cases.m` 为准：

```csv
analysis,parameter,role,default_lower_bound,default_upper_bound,notes
CRC,TSM.r_C,fitted,0.25,0.60,Tumor proliferation rate
CRC,par_DNA.d_0,fitted,0.10,0.35,DNA-damage-induced death rate
CRC,par_DNA.lambda_0,fitted,0.10,0.40,DNA-damage-associated growth inhibition
CRC,CY.p_1,fitted,1,10,PD-1 PD-L1 coupling parameter
CRC,CY.p_2,fitted for P01 P02 P03,1e3,1e6,Fixed at 1e6 when there is no anti-PD-1 administration
```

拟合配置为 `fixed_structure_pd1`。除上述拟合参数外，所有模型参数在八例间保持相同，包括 `C_max=1e11 cells`、`q_F=2.9928`、`q_P=0.50580 μg/mL`、`q_A=33.082 nM`、`lambda_F=0.5`、`lambda_P=1.5625`、`rho_F=0.1`、`a=0.03`、`b=0.003`、`eta_8=4e-8`、`eta_h=1e-8` 及原有放化疗相对损伤系数。三个药效耦合系数是共同校准值，不是逐患者拟合参数或独立实验测量值。病例的初始细胞数与实际给药日期分别保留；无抗 PD-1 给药时，`p_2` 不影响轨迹，不进行不可识别的优化。

初始癌细胞数由 `analysis/src/config/crc_cases.m` 给出。拟合最小化肿瘤体积残差平方和，Ki-67 用于对照。`results/crc_final_parameters.csv` 和 `results/crc_fit_summary.csv` 分别列出患者级拟合参数及指标。


## 固定药代与给药方案

各患者的拟合配置由 `analysis/src/config/crc_cases.m` 指定，PK 参数固定。CRC 日期来自 `Supplementary Data.xlsx` 的 `Table 5`：

| 病例 | 卡培他滨服药日（含首尾） | 奥沙利铂输注日 | 抗 PD-1 注射日 |
| --- | --- | --- | --- |
| P01 | 7–20 | 7 | 18 |
| P02 | 1–14、23–36 | 1、7、23 | 18、30 |
| P03 | 1–14、31–40 | 1、7、31 | 35 |
| P04 | 1–14、23–33 | 1、23 | 无 |
| P05 | 4–17、25–37 | 4、25 | 无 |
| P06 | 5–18、26–32 | 5、26 | 无 |
| P07 | 1–14、27–42 | 1、27 | 无 |
| P08 | 1–14、23–32 | 1、23 | 无 |

卡培他滨在表中标记的服药日内按每天两次输入；奥沙利铂仅在输注日进行两小时输入；抗 PD-1 仅在注射日增加浓度，此后继续清除。研究日 1 对应模拟时间 0；第 d 天的观测对应第 d 个轨迹点，随后进入第 d 天的治疗区间。PK 配置中的 21 天参考间隔用于从参考暴露推算单次剂量增量，不会自动安排注射。


图中阴影表示模拟体积上下 15% 的展示范围，不是统计置信区间。
