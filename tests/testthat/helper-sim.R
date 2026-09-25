# Simulate a Modal SARIMA series: phi(B)Phi(B^s) w_t = c + theta(B)Theta(B^s) eps_t,
# then integrate with (1-B)^d (1-B^s)^D.
sim_modal_sarima <- function(n, c = 0, phi = numeric(0), theta = numeric(0),
                             Phi = numeric(0), Theta = numeric(0), s = 12,
                             d = 0, D = 0, sigma = 1, gamma = 1,
                             dist = "normal", nu = NULL, burn = 300) {
  pp <- list(phi = phi, theta = theta, Phi = Phi, Theta = Theta)
  spec <- list(s = s)
  polys <- .lag_polys(pp, spec)
  N <- n + burn
  eps <- .rskd(N, 0, sigma, gamma, dist, nu)
  x <- c + .apply_poly(c(rep(0, length(polys$ma) - 1), eps), polys$ma)[length(polys$ma) - 1 + seq_len(N)]
  w <- if (length(polys$ar) > 1)
    as.numeric(stats::filter(x, -polys$ar[-1], method = "recursive")) else x
  w <- w[burn + seq_len(n)]
  y <- w
  if (D > 0) for (i in seq_len(D)) {
    z <- numeric(length(y))
    for (t in seq_along(y)) z[t] <- y[t] + (if (t > s) z[t - s] else 0)
    y <- z
  }
  if (d > 0) for (i in seq_len(d)) y <- cumsum(y)
  stats::ts(y, frequency = s)
}

num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))
}
