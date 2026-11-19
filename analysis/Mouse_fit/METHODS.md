## 数据处理与固定暴露尺度依据

路径以发布版根目录为起点。

四个数据集的观测值见 `analysis/src/config/mouse_cases.m`。LLC 的组均值和标准差由个体数据计算，其余数据集直接使用已汇总的观测值。

## 联合拟合

每个数据集内的全部实验组共享待估参数。仿真体积由 cm³ 转为 mm³，在观测时间点计算预测值与观测组均值的差，最小化所有组的未加权残差平方和。采用有界非线性最小二乘，使用三个起点，每个起点最多迭代 120 次。在相对坐标中优化各参数，`p_2` 使用以 10 为底的对数坐标。

LLC、MC38 和 A375 拟合 `r_C`、`C_max`、`p_1`、`p_2` 和 `epC`。CT26 的四组共享前四个参数，`epC` 在常规表达组设为 100，在低 PD-L1 表达组设为 1。

`PK.q_A` 表示每次抗 PD-1 给药引起的肿瘤部位有效浓度增量（nM），在 `analysis/src/config/mouse_cases.m` 中固定：LLC 和 MC38 为 44.7593，A375 为 104.4383，CT26 为 29.8395。采用 223 小时的清除半衰期。给药量由共同参考暴露 200 nM 按 `q_A = k_A × 参考给药间隔 × 200` 换算；各数据集的参考间隔依次为 3、3、7、2 天。

参考暴露的量级依据为 Cowles et al., *mAbs* (2022), [An affinity threshold for maximum efficacy in anti-PD-1 immunotherapy](https://doi.org/10.1080/19420862.2022.2088454), Fig. 3d 的肿瘤游离抗体模拟浓度。清除半衰期依据为 Boch et al., *Biomedicines* (2022), [Evaluating Antibody Pharmacokinetics as Prerequisite for Determining True Efficacy as Shown by Dual Targeting of PD-1 and CD96](https://doi.org/10.3390/biomedicines10092146), Fig. 2B 的 RMP1-14 小鼠测量值。

拟合参数边界以 `analysis/src/fitting/fit_mouse_case.m` 为准：

```csv
analysis,parameter,role,default_lower_bound,default_upper_bound,notes
Mouse,TSM.r_C,fitted,0.20,0.80,Joint control/intervention fit
Mouse,TSM.C_max,fitted,max(1.02x initial cancer count or 0.5x manual),2x manual,Joint control/intervention fit
Mouse,CY.p_1,fitted,0.1x manual,10x manual,PD-1/PD-L1 binding parameter
Mouse,CY.p_2,fitted,0.01x manual,100x manual,Effective association parameter in nM^-1 conditional on fixed q_A
Mouse,TSM.epC,fitted except CT26,0.1x manual,10x manual,CT26 arms explicitly fix epC
Mouse,PK.q_A,fixed,case-specific,case-specific,Effective tumor-site concentration increment per anti-PD-1 dose in nM
```


## 固定药代与给药方案

`q_A` 和药代时间常数固定，`p_2` 参与拟合。

LLC 第 9、12、15、18、21 天给药（Chen et al. 及对应 WO2025113697A1）；MC38 按 Zhang et al. 改为第 7、10、13 天三次给药；A375 按 Rosato et al. 每周给药；CT26 按 Liu et al. 每两天给药，起始日沿用输入模型的第 3 天。各方案在 `analysis/src/config/mouse_cases.m` 中列出。

运行 `run_mouse_analysis` 从 `input_data/final_parameters.mat` 重算曲线和拟合指标，输出到 `results/`。`plot_mouse_fit` 从 `results/mouse_results.mat` 读取轨迹绘图。
