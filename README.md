# Approximate Dispersion Relation of WENO5-JS and UW5

MATLAB implementation of the quasi-linear Approximate Dispersion Relation (ADR) analysis for the WENO5-JS and 5th-order upwind (UW5) schemes, reproducing the results of:

> Pirozzoli, S. (2006). *On the spectral properties of shock-capturing schemes*. Journal of Computational Physics, 219(2), 489–497.

## Method

For each reduced wavenumber `φ ∈ (0, π)` supported by a periodic grid of `N` points:

1. Initialize a single-mode field: `v_j(0) = cos(j·φ)`
2. Advance the linear advection equation `v_t + a·v_x = 0` for a small time `τ = σ·h/a` (with `σ ≪ 1`) using the target scheme
3. Compute the DFT and extract the complex amplitude at wavenumber `φ`
4. Recover the modified wavenumber:

$$\tilde{\phi}(\varphi) = \frac{i}{\sigma} \log\!\left(\frac{\hat{v}(\varphi,\tau)}{\hat{v}(\varphi,0)}\right)$$

- **Re(Φ)** → approximate phase speed (dispersion error)
- **Im(Φ)** → numerical dissipation (negative = stable)

The time integration uses the 3rd-order TVD Runge-Kutta scheme of Shu & Osher. The analytical modified wavenumber of UW5 is also computed as a sanity check.

## Schemes

| Scheme | Description |
|--------|-------------|
| **WENO5-JS** | 5th-order WENO with Jiang-Shu smoothness indicators (1996) |
| **UW5** | 5th-order linear upwind — equivalent to WENO5 with ideal weights frozen |

## Requirements

- MATLAB R2016b or later (uses local functions in scripts)

## Usage

```matlab
% All parameters are set at the top of main.m
N        = 1000;    % number of periodic grid points
sigma    = 1e-2;    % CFL number for ADR probe step (keep ≪ 1)
eps_weno = 1e-6;    % WENO5-JS regularization (Jiang-Shu default)

run('main.m')
```

The script prints a sanity-check error between the numerical and analytical UW5 modified wavenumbers, then produces a two-panel figure:

- **Top panel**: dispersion (Re(Φ) vs φ)
- **Bottom panel**: dissipation (Im(Φ) vs φ)

## Reference

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
