# ========== Internal core for Modal (S)ARIMA models ==========
#
# The model on the differenced series w_t = (1-B)^d (1-B^s)^D y_t is
#   phi(B) Phi(B^s) w_t = c + theta(B) Theta(B^s) eps_t,  eps_t ~ SKD(0, sigma, gamma),
# so the conditional mode is mu_t = w_t - eps_t. The non-seasonal model is the
# special case P = D = Q = 0.
#
# Parameter vector used by the optimizer:
#   (c, phi_1..p, theta_1..q, Phi_1..P, Theta_1..Q, log sigma, log gamma[, log(nu - 2)])

.polymul <- function(a, b) {
  r <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) {
    idx <- i:(i + length(b) - 1)
    r[idx] <- r[idx] + a[i] * b
  }
  r
}

# Spread the coefficients of a polynomial in B^s over powers of B.
.seasonal_poly <- function(poly, s) {
  out <- numeric((length(poly) - 1) * s + 1)
  out[seq(1, by = s, length.out = length(poly))] <- poly
  out
}

# (1 - B)^d (1 - B^s)^D
.diff_poly <- function(d, D, s) {
  r <- 1
  for (i in seq_len(d)) r <- .polymul(r, c(1, -1))
  for (i in seq_len(D)) r <- .polymul(r, .seasonal_poly(c(1, -1), s))
  r
}

.difference <- function(y, d, D, s) {
  w <- as.numeric(y)
  if (D > 0) w <- diff(w, lag = s, differences = D)
  if (d > 0) w <- diff(w, differences = d)
  w
}

.parse_seasonal <- function(seasonal, y) {
  if (is.null(seasonal)) seasonal <- list(order = c(0L, 0L, 0L))
  if (is.numeric(seasonal)) seasonal <- list(order = seasonal)
  if (is.null(seasonal$order) || length(seasonal$order) != 3)
    stop("'seasonal' must be a list with component 'order' = c(P, D, Q)")
  period <- seasonal$period
  if (is.null(period) || is.na(period)) period <- stats::frequency(y)
  order <- as.integer(seasonal$order)
  if (any(order < 0)) stop("seasonal orders must be non-negative")
  if (any(order > 0) && period <= 1)
    stop("a seasonal model needs a seasonal period > 1: pass 'seasonal$period' or use a ts with frequency > 1")
  list(order = order, period = as.integer(period))
}

.make_spec <- function(order, seasonal, dist) {
  so <- seasonal$order
  list(p = order[1], d = order[2], q = order[3],
       P = so[1], D = so[2], Q = so[3], s = seasonal$period, dist = dist)
}

.coef_names <- function(spec) {
  nm <- "intercept"
  if (spec$p > 0) nm <- c(nm, paste0("ar", seq_len(spec$p)))
  if (spec$q > 0) nm <- c(nm, paste0("ma", seq_len(spec$q)))
  if (spec$P > 0) nm <- c(nm, paste0("sar", seq_len(spec$P)))
  if (spec$Q > 0) nm <- c(nm, paste0("sma", seq_len(spec$Q)))
  nm <- c(nm, "sigma", "gamma")
  if (spec$dist == "t") nm <- c(nm, "nu")
  nm
}

.nu_from_par <- function(x) min(exp(x), 200) + 2

# Split an optimizer vector into named components.
.split_par <- function(par, spec) {
  i <- 1
  take <- function(k) { v <- if (k > 0) par[i + seq_len(k)] else numeric(0); i <<- i + k; v }
  out <- list(c = par[1])
  out$phi <- take(spec$p); out$theta <- take(spec$q)
  out$Phi <- take(spec$P); out$Theta <- take(spec$Q)
  out$sigma <- exp(par[i + 1]); out$gamma <- exp(par[i + 2])
  out$nu <- if (spec$dist == "t") .nu_from_par(par[i + 3]) else NULL
  out
}

# Coefficients on the natural scale -> optimizer vector.
.par_from_coef <- function(coefs, spec) {
  nb <- 1 + spec$p + spec$q + spec$P + spec$Q
  par <- c(as.numeric(coefs[seq_len(nb)]), log(coefs["sigma"]), log(coefs["gamma"]))
  if (spec$dist == "t") par <- c(par, log(coefs["nu"] - 2))
  unname(par)
}

.lag_polys <- function(pp, spec) {
  list(ar = .polymul(c(1, -pp$phi), .seasonal_poly(c(1, -pp$Phi), spec$s)),
       ma = .polymul(c(1, pp$theta), .seasonal_poly(c(1, pp$Theta), spec$s)))
}

.roots_outside <- function(poly) {
  if (length(poly) <= 1 || all(poly[-1] == 0)) return(TRUE)
  all(Mod(polyroot(poly)) > 1.001)
}

.admissible <- function(pp) {
  .roots_outside(c(1, -pp$phi)) && .roots_outside(c(1, -pp$Phi)) &&
    .roots_outside(c(1, pp$theta)) && .roots_outside(c(1, pp$Theta))
}

# Apply a lag polynomial (coefficients in increasing powers of B) to x.
.apply_poly <- function(x, poly) {
  if (length(poly) == 1) return(poly * x)
  as.numeric(stats::filter(x, poly, method = "convolution", sides = 1))
}

# eps_t = beta(B)^{-1} x_t with zero pre-sample values.
.ma_inverse <- function(x, ma) {
  if (length(ma) == 1) return(x)
  if (is.matrix(x)) return(apply(x, 2, .ma_inverse, ma = ma))
  as.numeric(stats::filter(x, -ma[-1], method = "recursive"))
}

# Modal residual recursion. Pre-sample values of w are set to the sample mean of
# w and pre-sample residuals to zero.
.recursion <- function(pp, spec, w) {
  polys <- .lag_polys(pp, spec)
  m <- length(polys$ar) - 1
  n <- length(w)
  w_ext <- c(rep(mean(w), m), w)
  x <- .apply_poly(w_ext, polys$ar)[m + seq_len(n)] - pp$c
  eps <- .ma_inverse(x, polys$ma)
  list(eps = eps, mu = w - eps, polys = polys, w_ext = w_ext, m = m)
}

.neg_loglik <- function(par, spec, w) {
  pp <- .split_par(par, spec)
  if (!.admissible(pp)) return(1e10)
  rec <- .recursion(pp, spec, w)
  ans <- -sum(.dskd(rec$eps, 0, pp$sigma, pp$gamma, spec$dist, pp$nu, log = TRUE))
  if (!is.finite(ans)) ans <- 1e10
  ans
}

# Analytical gradient of the negative log-likelihood.
.grad_neg_loglik <- function(par, spec, w) {
  pp <- .split_par(par, spec)
  if (!.admissible(pp)) return(rep(0, length(par)))
  colSums(.scores(par, spec, w))
}

# Per-observation gradients of the negative log-likelihood (n x length(par)).
.scores <- function(par, spec, w) {
  pp <- .split_par(par, spec)
  rec <- .recursion(pp, spec, w)
  n <- length(w); m <- rec$m; s <- spec$s
  eps <- rec$eps
  r <- length(rec$polys$ma) - 1
  eps_ext <- c(rep(0, r), eps)

  # Direct derivatives Z of mu_t (residuals held fixed), one column per coefficient.
  lagged <- function(x_ext, poly, lag, pad) {
    u <- .apply_poly(x_ext, poly)
    u[pad + seq_len(n) - lag]
  }
  Z <- matrix(1, n, 1)
  sphi <- .seasonal_poly(c(1, -pp$Phi), s)
  stheta <- .seasonal_poly(c(1, pp$Theta), s)
  for (i in seq_len(spec$p)) Z <- cbind(Z, lagged(rec$w_ext, sphi, i, m))
  for (j in seq_len(spec$q)) Z <- cbind(Z, lagged(eps_ext, stheta, j, r))
  for (i in seq_len(spec$P)) Z <- cbind(Z, lagged(rec$w_ext, c(1, -pp$phi), i * s, m))
  for (j in seq_len(spec$Q)) Z <- cbind(Z, lagged(eps_ext, c(1, pp$theta), j * s, r))

  # V_t = d mu_t / d beta satisfies V_t = Z_t - sum_k beta_k V_{t-k}.
  V <- .ma_inverse(Z, rec$polys$ma)
  if (!is.matrix(V)) V <- matrix(V, ncol = 1)

  k <- .skd_kernel(spec$dist, pp$nu)
  p <- .p_skew(pp$gamma)
  z <- eps / pp$sigma
  x2 <- 2 * .rho(z, p)
  L1 <- k$dlog(x2)
  score <- L1 * 2 * (p - (z < 0)) / pp$sigma          # d log f / d eps

  dp_dlg <- -2 * pp$gamma^2 / (1 + pp$gamma^2)^2
  S <- cbind(score * V,
             1 + L1 * x2,
             -((1 - 2 * p) / (p * (1 - p)) + L1 * 2 * z) * dp_dlg)

  if (spec$dist == "t") {
    nu <- pp$nu
    dnu <- 0.5 * digamma((nu + 1) / 2) - 0.5 * digamma(nu / 2) - 1 / (2 * nu) -
      0.5 * log1p(x2^2 / nu) + (nu + 1) * x2^2 / (2 * nu^2 * (1 + x2^2 / nu))
    raw <- exp(par[length(par)])
    S <- cbind(S, -dnu * (if (raw < 200) raw else 0))
  }
  unname(S)
}

# Observed information on the optimization scale: the Hessian of the negative
# log-likelihood, or the outer product of the per-observation scores when the
# log-likelihood is not smooth (Skewed Laplace) or the Hessian is not positive
# definite.
.information <- function(par, spec, w, fn, gr) {
  if (spec$dist != "laplace") {
    H <- stats::optimHess(par, fn, gr)
    free <- abs(diag(H)) > 1e-8
    ev <- tryCatch(eigen(H[free, free, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) -1)
    if (all(is.finite(ev)) && min(ev) > 0) return(list(matrix = H, type = "hessian"))
  }
  S <- .scores(par, spec, w)
  list(matrix = crossprod(S), type = "opg")
}

# Simulate n values of the differenced series w from the model, after a burn-in.
.simulate_differenced <- function(pp, spec, n, burn = 200) {
  polys <- .lag_polys(pp, spec)
  r <- length(polys$ma) - 1
  N <- n + burn
  eps <- .rskd(N + r, 0, pp$sigma, pp$gamma, spec$dist, pp$nu)
  x <- pp$c + .apply_poly(eps, polys$ma)[r + seq_len(N)]
  w <- if (length(polys$ar) > 1)
    as.numeric(stats::filter(x, -polys$ar[-1], method = "recursive")) else x
  w[burn + seq_len(n)]
}

# Differenced series, residuals and modes for a fitted object.
.modal_filter <- function(object) {
  spec <- .object_spec(object)
  w <- .difference(object$y, spec$d, spec$D, spec$s)
  pp <- .split_par(.par_from_coef(object$coefficients, spec), spec)
  rec <- .recursion(pp, spec, w)
  list(w = w, eps = rec$eps, mu = rec$mu, pp = pp, spec = spec, polys = rec$polys)
}

# Specification of a fitted object; objects from version 0.1.0 had no seasonal part.
.object_spec <- function(object) {
  seasonal <- object$seasonal
  if (is.null(seasonal)) seasonal <- list(order = c(0L, 0L, 0L), period = 1L)
  dist <- if (!is.null(object$dist)) object$dist else "normal"
  .make_spec(object$order, seasonal, dist)
}

.is_seasonal <- function(spec) any(c(spec$P, spec$D, spec$Q) > 0)

.model_label <- function(object) {
  spec <- .object_spec(object)
  lab <- paste0("(", spec$p, ",", spec$d, ",", spec$q, ")")
  if (.is_seasonal(spec))
    paste0("Modal SARIMA", lab, "(", spec$P, ",", spec$D, ",", spec$Q, ")[", spec$s, "]")
  else paste0("Modal ARIMA", lab)
}
