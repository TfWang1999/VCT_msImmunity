# R 环境

使用 Windows R 4.4.1；安装源代码包时需要 Rtools44。在 `sc_analysis` 目录运行：

```r
source("../config/sc_environment/install_sc_environment.R")
```

`sc_analysis/renv.lock` 固定包版本，`package_sources.csv` 给出下载地址、安装类型和压缩包 MD5。
安装入口沿用清单中的版本，已安装且版本相符的包不重复安装。

肝癌分析使用 GSVA 2.0.7（Bioconductor 3.20）和 pheatmap 1.0.13；两者及所需依赖已纳入上述清单。
肝癌脚本使用同一包库，运行方式见 `Liver_cancer/README.md`。
