# WENO5-JS 与 UW5 格式的近似色散关系分析

MATLAB 实现的准线性近似色散关系（ADR）分析，适用于 WENO5-JS 和五阶迎风（UW5）格式，复现了以下文献的结果：

> Pirozzoli, S. (2006). *On the spectral properties of shock-capturing schemes*. Journal of Computational Physics, 219(2), 489–497.

## 方法

对周期网格（共 `N` 个格点）中每个约化波数 `φ ∈ (0, π)`：

1. 初始化单模场：`v_j(0) = cos(j·φ)`
2. 以目标格式推进线性对流方程 `v_t + a·v_x = 0`，时间步长 `τ = σ·h/a`（`σ ≪ 1`）
3. 计算 DFT，提取波数 `φ` 处的复数振幅
4. 还原修正波数：

$$\tilde{\phi}(\varphi) = \frac{i}{\sigma} \log\!\left(\frac{\hat{v}(\varphi,\tau)}{\hat{v}(\varphi,0)}\right)$$

- **Re(Φ)** → 近似相速度（色散误差）
- **Im(Φ)** → 数值耗散（负值表示稳定）

时间积分采用 Shu & Osher 的三阶 TVD Runge-Kutta 格式。同时计算 UW5 的解析修正波数作为正确性验证。

## 格式说明

| 格式 | 描述 |
|------|------|
| **WENO5-JS** | 五阶 WENO 格式，使用 Jiang-Shu 光滑度指示子（1996） |
| **UW5** | 五阶线性迎风格式，等价于权重固定为理想权重的 WENO5 |

## 环境要求

- MATLAB R2016b 或更高版本（脚本中使用了局部函数特性）

## 使用方法

```matlab
% 所有参数在 main.m 顶部设置
N        = 1000;    % 周期网格格点数
sigma    = 1e-2;    % ADR 探测步的 CFL 数（保持 ≪ 1）
eps_weno = 1e-6;    % WENO5-JS 正则化参数（Jiang-Shu 默认值）

run('main.m')
```

脚本运行后将输出 UW5 数值与解析修正波数之间的误差（正确性检验），并生成双面板图：

- **上面板**：色散关系（Re(Φ) vs φ）
- **下面板**：耗散关系（Im(Φ) vs φ）

## 参考文献

```bibtex
@article{pirozzoli2006spectral,
  author  = {Pirozzoli, Sergio},
  title   = {On the spectral properties of shock-capturing schemes},
  journal = {Journal of Computational Physics},
  volume  = {219},
  number  = {2},
  pages   = {489--497},
  year    = {2006}
}
```
