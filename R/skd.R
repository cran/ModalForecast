# ========== Internal SKD distribution functions ==========
#
# All three members use the parameterization of Galarza et al. (2017):
#   f(y | mu, sigma, gamma) = 4p(1-p)/sigma * g(2 * rho_p((y - mu)/sigma)),
# with p = 1/(1 + gamma^2), rho_p(u) = u * (p - I(u < 0)) and g a standard
# symmetric kernel (normal, Student-t or Laplace). The mode is mu, the mode is
# also the p-th quantile, and sigma and gamma have the same meaning for all
# three members.

.p_skew <- function(gamma) 1 / (1 + gamma^2)

.rho <- function(z, p) z * (p - (z < 0))

# Kernel g: log-density, CDF, quantile and derivative of the log-density.
.skd_kernel <- function(dist, nu = NULL) {
  switch(dist,
    normal = list(
      logd = function(x) stats::dnorm(x, log = TRUE),
      cdf  = stats::pnorm,
      qf   = stats::qnorm,
      dlog = function(x) -x),
    t = list(
      logd = function(x) stats::dt(x, df = nu, log = TRUE),
      cdf  = function(x) stats::pt(x, df = nu),
      qf   = function(u) stats::qt(u, df = nu),
      dlog = function(x) -(nu + 1) * x / (nu + x^2)),
    laplace = list(
      logd = function(x) -log(2) - abs(x),
      cdf  = function(x) ifelse(x < 0, 0.5 * exp(pmin(x, 0)), 1 - 0.5 * exp(-pmax(x, 0))),
      qf   = function(u) ifelse(u < 0.5, log(2 * pmin(u, 0.5)), -log(2 * (1 - pmax(u, 0.5)))),
      dlog = function(x) -sign(x)),
    stop("Unknown distribution: ", dist))
}

.dskd <- function(x, mu = 0, sigma = 1, gamma = 1, dist = "normal", nu = NULL, log = FALSE) {
  k <- .skd_kernel(dist, nu)
  p <- .p_skew(gamma)
  z <- (x - mu) / sigma
  ld <- log(4 * p * (1 - p)) - log(sigma) + k$logd(2 * .rho(z, p))
  if (log) ld else exp(ld)
}

.pskd <- function(x, mu = 0, sigma = 1, gamma = 1, dist = "normal", nu = NULL) {
  k <- .skd_kernel(dist, nu)
  p <- .p_skew(gamma)
  z <- (x - mu) / sigma
  ifelse(z < 0,
         2 * p * k$cdf(2 * (1 - p) * pmin(z, 0)),
         p + 2 * (1 - p) * (k$cdf(2 * p * pmax(z, 0)) - 0.5))
}

.qskd <- function(prob, mu = 0, sigma = 1, gamma = 1, dist = "normal", nu = NULL) {
  k <- .skd_kernel(dist, nu)
  p <- .p_skew(gamma)
  prob <- rep_len(prob, max(length(prob), length(sigma)))
  z <- numeric(length(prob))
  lo <- prob < p
  if (any(lo)) z[lo] <- k$qf(prob[lo] / (2 * p)) / (2 * (1 - p))
  if (any(!lo)) z[!lo] <- k$qf(0.5 + (prob[!lo] - p) / (2 * (1 - p))) / (2 * p)
  mu + sigma * z
}

.rskd <- function(n, mu = 0, sigma = 1, gamma = 1, dist = "normal", nu = NULL) {
  .qskd(stats::runif(n), mu = mu, sigma = sigma, gamma = gamma, dist = dist, nu = nu)
}

.dist_label <- function(dist) {
  switch(dist, normal = "Skew-Normal", t = "Skewed Student-t", laplace = "Skewed Laplace")
}
