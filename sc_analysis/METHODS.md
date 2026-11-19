# 单细胞分析方法

## scImmunity 公式和输出字段

四个队列均使用比例半饱和常数 `kp=0.15`、功能半饱和常数 `kf=0.5`。
比例的促进和抑制函数分别为 `A(p)=p/(p+kp)`、`B(p)=kp/(p+kp)`；
功能的促进和抑制函数分别为 `Hp(z)=z/(z+kf)`、`Hi(z)=kf/(z+kf)`。

四种 T 细胞均将存活促进项 `Hp(1+S)` 乘入功能输入，再通过外层功能
Hill 函数；Th、CD8_em、CD8_ex 使用促进函数，Treg 使用抑制函数。
当前代码的功能项为：

```r
Fh = Hp(1 + etah * Hp(1+Act_4) * Hp(1+Costim_4) * Hp(1+IFN_4) * Hp(1+S_4))
Fr = Hi(1 + Treg_sig * Hp(1+S_Treg))
Fe = Hp(1 + eta8 * Hp(1+Act8) * Hp(1+IFN_8) * Hp(1+S_8))
Fx = Hp(1 + eta8_ex / (1+Exh_ex) * Hp(1+S_ex))
F1 = Hp(1 + m1_sig)
F2 = Hi(1 + m2_sig)
```

完整 scImmunity 得分由三个 T 细胞贡献和适用队列的巨噬细胞贡献组成：

```r
T_cell_score = A(Th_pro) * Fh * B(Treg_pro) * Fr +
               A(T8em_pro) * Fe + A(T8ex_pro) * Fx
TAM_score = A(M1_pro) * F1 * B(M2_pro) * F2
```

CRC、NSCLC_1 的完整得分为 `T_cell_score + TAM_score`；NSCLC_2、BCC 的
完整得分为 `T_cell_score`。CD8_em 与 CD8_ex 等权相加。Th 的贡献包含
Treg 比例和功能的双重抑制。各亚群的存活比值为
`S=(1+mean(Anti-apoptosis))/(1+mean(Pro-apoptosis))`，均值在对应亚群内计算。

仅表达得分为 `Fh*Fr + Fe + Fx`，CRC、NSCLC_1 额外加 `F1*F2`；
仅比例得分为 `A(Th_pro)*B(Treg_pro) + A(T8em_pro) + A(T8ex_pro)`，
CRC、NSCLC_1 额外加 `A(M1_pro)*B(M2_pro)`。
完整实现见 `scImmunity/scoring.R`，输入特征的提取见 `scImmunity/features.R`。

## 临床统计与图形

CRC、NSCLC_1 的两组比较使用双侧 Wilcoxon 秩和检验；CRC 总体 TRG 比较使用
Kruskal–Wallis 检验。NSCLC_2、BCC 使用单侧“应答者评分更高”的 Wilcoxon 检验。
统计表同时提供双侧 p 值和具有明确方向的单侧 p 值；LUAD/LUSC 不指定单侧方向。
`P_value` 为绘图使用的值，`P_two_sided` 和 `P_one_sided` 为对应方向的结果。
所有 p 值为未校正的名义 p 值；NSCLC_2 原有重复患者/时间点条目保持原样。

| 图形内容 | 队列与比较 |
| --- | --- |
| CRC 应答 | CRC：Responder vs Non-responder。 |
| CRC 疗效 | CRC：CR vs SD、PR vs SD。 |
| CRC TRG | CRC：TRG0 vs TRG3、TRG2 vs TRG3。 |
| NSCLC 影像学应答 | NSCLC_1：影像学 Responder vs Non-responder。 |
| NSCLC 病理应答 | NSCLC_1：pCR vs 全部已知 non-pCR（MPR + non-MPR，41 vs 98），排除未知病理分组。 |
| NSCLC 组织学类型 | NSCLC_1：LUAD vs LUSC。 |
| 三种评分比较 | 四队列共 11 种临床比较，分别展示完整、仅比例和仅表达评分的 p 值；不包含 MPR vs non-MPR。 |

共有 35 项不同的统计结果。临床图使用其中 8 项完整评分比较；三种评分比较图使用
33 项结果，其中部分与临床图重合。绘图保持既有颜色和排版；抖动随机种子为 123，
PDF 为矢量格式，PNG 分辨率为 300 dpi。

## 输入与上游准备

四个队列分别为 CRC（GSE236581）、NSCLC_1（GSE243013）、NSCLC_2（GSE179994）和 BCC（GSE123813）。
现有 df_data 和临床映射位于 `scImmunity/input_data/<队列>/`，共 22、140、47、11 条样本记录；记录数不等同于独立患者数。
NSCLC_1 在具有影像学应答标签且六类细胞均存在的样本中计算和归一化模块分数，随后筛选六类细胞各占保留细胞 >2% 且各至少 20 个的样本。
这些输入足以复现评分和临床统计。`prepare_scores.R` 保留原始 Seurat 对象的细胞注释、功能模块计算与归一化函数；从原始 RNA 重建 df_data 还需要各队列原始检查点及 CD4/CD8 marker 列表，这些完整上游输入尚未齐备。

GO 使用 `GO_enrichment/input_data/CRC_Treatment_I.rds` 和固定基线评分，高分组为 P21、P26、P11、P20、P09，低分组为 P17、P16、P25、P02、P12。
`FindMarkers` 使用 `logfc.threshold=0.25`、`min.pct=0.1`；上下调基因分别要求调整后 p 值小于 0.05、绝对 log2FC 大于 0.25，按效应大小各取前 200 个并保留并列。
GO 使用 BP、BH 校正、`pvalueCutoff=0.05` 和语义去重阈值 0.7。显示条目的关键词和 GeneRatio 筛选保留在 `GO_enrichment/go_analysis.R`。

CellChat 使用 `CellChat/input_data/` 中已保存的 P21/P12 对象。比较时的差异表达阈值为 `thresh.pc=0.05`、`thresh.fc=0.01`、`thresh.p=0.05`；上下调配体和受体 logFC 阈值分别为 0.7 和 −0.05。
弦图发送端和接收端均限定为 CD4T、CD8T、TAM、Treg、tumor；上调弦图扇区间距为 4，仅影响布局。
`CellChat/infer_cellchat.R` 保留原始通讯概率推断函数，需要带细胞注释的 RNA 对象。现有比较入口直接读取已保存的 CellChat 对象。

环境由 `renv.lock` 固定：R 4.4.1、Seurat 5.3.0、CellChat 2.2.0、clusterProfiler 4.14.6、org.Hs.eg.db 3.20.0、enrichplot 1.26.6、GO.db 3.20.0；presto 1.0.0 使用提交 `7636b3d0465c468c35853f82f1717d3a64b3c8f6`。
