#import "@local/modern-cug-report:0.1.3": *

#let Ic = none

#let length = none

#let Ips = none
#show: doc => template(doc, footer: "水文气象学", header: "")

= 操作步骤与存在的问题

== 1步骤

+ 使用AI通读 `component.R`，识别出四个函数：`RemoveSeasonal`、`LS_MultiBreaks`、`LS_OneBreak`、`lm_solve`。
+ 逐函数理解：去季节（按月日分组求均值再求距平）、分段回归（同斜率不同截距，构造指示变量设计矩阵）、最小二乘求解（`.lm.fit`）。
+ 对照 `docs` 下已有 `note_01.typ` 的模板与排版风格，保持标题、代码块、公式写法一致。
+ 为每个函数配公式 + R 代码 + 简要说明，整理成笔记。

== 2存在的问题（AI发现）

+ 文档与代码不一致：`LS_OneBreak` 的 `@param` 写 `T_anorm`，实参却是 `X`；注释 `T_anorm = T - mean(T)` 但代码未中心化。
+ 命名不统一：`Y/X` 大小写混用，`anorm` 与 `anomaly`、`T_anorm` 与 `X` 并存。
+ 拼写错误：`perid` 应为 `period`。
+ 数值稳健性：`t_mmdd == inds[i]` 用浮点相等比较，非整数或浮点编码时有精度风险，宜用 `match()` 或整数编码。
+ 全 `NA` 组未处理：某月日样本全为 `NA` 时 `mean(..., na.rm = TRUE)` 返回 `NaN`，距平随之失效。
+ `lm_solve` 未处理 `NA` 与秩亏设计矩阵，`.lm.fit` 遇奇异矩阵可能报错。
+ 性能：`RemoveSeasonal` 逐组循环，可改用 `tapply` / `aggregate` 向量化。

= component.R 笔记

== 1 RemoveSeasonal：去季节

给定观测值 $y$ 及其月日标签 $t_(m m d d)$，按月日升序求季节均值，再求距平。

$
  mu_d = 1/n_d sum_(t_(m m d d)=d) y,
  quad
  a_t = y_t - mu_(t_(m m d d))
$

其中 $d$ 为去重后的月日，$n_d$ 为该月日的样本数。月尺度最多 12 个、日尺度最多 366 个。

```R
RemoveSeasonal <- function(y, t_mmdd) {
  inds <- sort(unique(t_mmdd))
  mu <- numeric(length(inds))
  anomaly <- numeric(length(y))
  for (i in seq_along(inds)) {
    ind <- which(t_mmdd == inds[i])
    mu[i] <- mean(y[ind], na.rm = TRUE)
    anomaly[ind] <- y[ind] - mu[i]
  }
  list(mu = mu, anomaly = anomaly)
}
```

== 2 LS_MultiBreaks：多断点分段回归

同斜率、不同截距的分段线性回归。$Ips$ 为各段末尾位置（末位为序列长度）。

$
  y = a + b (x - overline(x)) + sum_(i=1)^(m-1) c_i bold(1)_((Ips_i, Ips_(i+1)])
$

其中 $m = length(Ips)$，段数为 $m-1$。设计矩阵第一列为 1，第二列为中心化自变量，其后每列对应一个分段指示变量。

```R
LS_MultiBreaks <- function(Y, X, Ips, ...) {
  n_breaks <- length(Ips) - 1
  design <- cbind(1, X - mean(X))
  for (i in seq_len(n_breaks)) {
    segment <- numeric(length(Y))
    segment[(Ips[i] + 1):Ips[i + 1]] <- 1
    design <- cbind(design, segment)
  }
  lm_solve(design, Y, ...)
}
```

== 3 LS_OneBreak：单断点分段回归

$Ic$ 为变点位置，$[1, Ic]$ 与 $[Ic+1, n]$ 为前后两段。

$
  y = a + b x + c bold(1)_({t > Ic})
$

$Ic = "NA"$ 时退化为普通一元线性回归。

```R
LS_OneBreak <- function(Y, X, Ic, ...) {
  n <- length(Y)
  design <- cbind(1, X)
  if (!is.na(Ic)) {
    design <- cbind(design, c(rep(0, Ic), rep(1, n - Ic)))
  }
  lm_solve(design, Y, ...)
}
```

== 4 lm_solve：最小二乘求解

用 R 内部例程 `.lm.fit` 直接求解，避免构造并求逆正规方程，更快更稳。

$
  beta = (X^T X)^(-1) X^T y
$

返回 `coef`（系数）、`fitted`（拟合值）、`resi`（残差）、`SSE`（残差平方和）；`only.SSE = TRUE` 时只返回 $S S E$。

```R
lm_solve <- function(design, response, only.SSE = FALSE, ...) {
  l <- .lm.fit(design, as.matrix(response))
  fitted <- design %*% l$coefficients
  SSE <- sum(l$residuals^2)
  if (!only.SSE) {
    list(
      coef = l$coefficients,
      fitted = as.vector(fitted),
      resi = l$residuals,
      SSE = SSE
    )
  } else {
    list(SSE = SSE)
  }
}
```
