# ModalForecast 0.2.0

## New features

* Seasonal models. `fit_modal_arima()` gains a `seasonal` argument, with the
  same form as in `stats::arima()`, to fit Modal SARIMA(p,d,q)(P,D,Q)[s]
  models. The seasonal coefficients are reported as `sar1`, `sma1`, and so on.
  Existing non-seasonal calls are unchanged.
* `auto.modal.arima()` searches over seasonal orders (`max.P`, `max.Q`), chooses
  the seasonal differencing order `D` with `forecast::nsdiffs()`, and gains a
  `trace` argument.
* `forecast()` gains `point = c("joint", "marginal")`. The default, `"joint"`, is
  the most probable future trajectory, as in 0.1.0. `"marginal"` gives the
  conditional mode of each future value, which differs from the joint path for
  skewed errors and horizons greater than one. The returned object also
  contains `mode_path` and `mode_shift`.
* Asymptotic prediction intervals now use the exact distribution of the
  forecast error at the estimated parameters, computed by numerical
  convolution of the SKD densities.
* New `vcov()` and `nobs()` methods. `summary()` reports a test of symmetry,
  H0: gamma = 1.
* `diagnostics()` and `envelope()` support seasonal models; for these the
  residual ACF, PACF and Ljung-Box test cover two seasonal periods.
* `envelope()` gains `refit = TRUE` (the default): each replication simulates a
  series from the fitted model and re-estimates it, so the envelope accounts
  for parameter estimation. Without it, the smallest distances of a correct
  Skewed Laplace fit fall below the envelope, because its maximum likelihood
  fit pulls several residuals to zero. `refit = FALSE` keeps the faster
  envelope based on simulated innovations.
* Sergio Luis Mercado Londoño and Víctor Hugo Lachos join as contributors.

## Changes in parameterization

* The Skew-Normal now uses the parameterization of Galarza et al. (2017),
  `f(y) = 4p(1-p)/sqrt(2*pi*sigma^2) * exp(-2*rho_p((y-mu)/sigma)^2)`, the same
  as the Skewed Student-t and Skewed Laplace. `sigma` therefore has the same
  meaning for the three distributions. For Skew-Normal fits, the log-likelihood,
  the systematic coefficients and `gamma` are unchanged; the reported `sigma`
  equals the 0.1.0 value divided by `(gamma + 1/gamma)/2`.
* `BIC()` and `logLik()` use the number of observations after differencing.

## Bug fixes

* The analytical gradients of the Skewed Student-t and Skewed Laplace
  log-likelihoods had errors, so the optimizer could stop far from the maximum
  while reporting convergence. For example, a Skewed Laplace ARIMA(2,1,1) for
  `log10(lynx)` reached a log-likelihood of -189.99 in 0.1.0 and -0.14 now.
  Estimates from these two distributions may change substantially.
* The Skewed Laplace fit is refined with Nelder-Mead after BFGS, because its
  log-likelihood is not differentiable when a residual is zero, and its
  standard errors use the outer product of the scores instead of the Hessian.
* Standard errors of `sigma`, `gamma` and `nu` in `summary()` were those of
  their logarithms; they now use the delta method.
* The random generators used for bootstrap intervals and envelopes did not
  match the fitted Student-t and Laplace densities (the Laplace draws were half
  as wide even for `gamma = 1`). All draws now come from the fitted density.
* The Skewed Laplace CDF was wrong below the mode, which affected the lower
  asymptotic prediction bounds and the randomized quantile residuals.
* The theoretical quantiles of the Laplace envelope were half their correct
  value.
* Removed an unrelated script, `inst/extras/devTests.R`.

# ModalForecast 0.1.0

* Initial release.
* Added `fit_modal_arima()` for fitting Parametric Modal ARIMA models using the generalized SKD family.
* Implemented S3 methods: `print`, `coef`, `AIC`, `BIC`, `summary`, and `plot`.
* Added comprehensive simulation tests demonstrating parameter recovery (c, ar1, sigma, gamma).
