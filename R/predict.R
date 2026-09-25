#' Modal Trajectory Prediction for Modal (S)ARIMA Models
#'
#' Returns the joint modal trajectory of the future values: the recursion is
#' iterated with future innovations set to their mode, zero. For the one-step
#' horizon this is the conditional mode of \eqn{y_{T+1}}; for longer horizons
#' it is the most probable path, see \code{\link{forecast.modal_arima}}.
#'
#' @param object A \code{modal_arima} object.
#' @param n.ahead The forecast horizon.
#' @param ... Additional arguments (unused).
#' @return A numeric vector with the modal trajectory on the scale of \code{y}.
#' @export
predict.modal_arima <- function(object, n.ahead = 10, ...) {
  st <- .modal_filter(object)
  .modal_path(st, object$y, n.ahead)
}

# Joint modal path: iterate the model with future innovations equal to zero.
.modal_path <- function(st, y, h) {
  w_future <- .simulate_w(st, matrix(0, h, 1))[, 1]
  .integrate_path(w_future, y, st$spec)
}

# Iterate the model forward from the end of the sample for each column of the
# matrix of future innovations E (h x npaths). Returns future w values.
.simulate_w <- function(st, E) {
  h <- nrow(E); npaths <- ncol(E)
  alpha <- -st$polys$ar[-1]; beta <- st$polys$ma[-1]
  m <- length(alpha); r <- length(beta)
  n <- length(st$w)
  w_hist <- c(rep(mean(st$w), m), st$w)
  e_hist <- c(rep(0, r), st$eps)
  W <- matrix(rep(w_hist, npaths), ncol = npaths)
  Eall <- rbind(matrix(rep(e_hist, npaths), ncol = npaths), E)
  W <- rbind(W, matrix(0, h, npaths))
  for (k in seq_len(h)) {
    iw <- m + n + k; ie <- r + n + k
    val <- st$pp$c + Eall[ie, ]
    for (j in seq_len(m)) if (alpha[j] != 0) val <- val + alpha[j] * W[iw - j, ]
    for (j in seq_len(r)) if (beta[j] != 0) val <- val + beta[j] * Eall[ie - j, ]
    W[iw, ] <- val
  }
  W[m + n + seq_len(h), , drop = FALSE]
}

# Undo (1-B)^d (1-B^s)^D: y_t = w_t - sum_k delta_k y_{t-k}.
.integrate_path <- function(w_future, y, spec) {
  dp <- .diff_poly(spec$d, spec$D, spec$s)
  if (length(dp) == 1) return(w_future)
  W <- as.matrix(w_future)
  k <- length(dp) - 1
  yh <- as.numeric(y)
  Y <- rbind(matrix(rep(utils::tail(yh, k), ncol(W)), ncol = ncol(W)), matrix(0, nrow(W), ncol(W)))
  for (t in seq_len(nrow(W))) {
    val <- W[t, ]
    for (j in seq_len(k)) if (dp[j + 1] != 0) val <- val - dp[j + 1] * Y[k + t - j, ]
    Y[k + t, ] <- val
  }
  out <- Y[k + seq_len(nrow(W)), , drop = FALSE]
  if (ncol(out) == 1) out[, 1] else out
}

# psi-weights of y_{T+h} with respect to future innovations (psi_0 = 1).
.psi_weights <- function(st, h) {
  if (h == 1) return(1)
  ar_full <- .polymul(st$polys$ar, .diff_poly(st$spec$d, st$spec$D, st$spec$s))
  c(1, stats::ARMAtoMA(ar = -ar_full[-1], ma = st$polys$ma[-1], lag.max = h - 1))
}

# Exact distribution of sum_j psi_j eps_j (eps_j iid SKD) by FFT convolution.
# Returns the mode and a quantile function of the sum.
.error_distribution <- function(psi, pp, dist, ngrid = 2048) {
  if (length(psi) == 1)
    return(list(mode = 0, quantile = function(u) psi * .qskd(u, 0, pp$sigma, pp$gamma, dist, pp$nu)))
  tail_p <- 1e-6
  lo <- .qskd(tail_p, 0, pp$sigma, pp$gamma, dist, pp$nu)
  hi <- .qskd(1 - tail_p, 0, pp$sigma, pp$gamma, dist, pp$nu)
  psi <- psi[abs(psi) > 1e-10]
  dx <- sum(abs(psi)) * (hi - lo) / ngrid
  dens <- NULL; start <- 0
  for (a in psi) {
    b <- sort(c(a * lo, a * hi))
    x <- seq(b[1], b[2], by = dx)
    if (length(x) < 3) next                  # negligible component
    f <- .dskd(x / a, 0, pp$sigma, pp$gamma, dist, pp$nu) / abs(a)
    f <- f / (sum(f) * dx)
    if (is.null(dens)) { dens <- f; start <- x[1] }
    else {
      dens <- pmax(stats::convolve(dens, rev(f), type = "open") * dx, 0)
      start <- start + x[1]
    }
  }
  grid <- start + dx * (seq_along(dens) - 1)
  cdf <- cumsum(dens) * dx
  cdf <- cdf / cdf[length(cdf)]
  keep <- !duplicated(cdf)
  list(mode = grid[which.max(dens)],
       quantile = function(u) stats::approx(cdf[keep], grid[keep], xout = u, rule = 2)$y)
}

#' @importFrom forecast forecast
#' @export
forecast::forecast

#' Forecasting with Modal (S)ARIMA Models
#'
#' Produces modal point forecasts and prediction intervals for a fitted
#' \code{modal_arima} object, with or without a seasonal component.
#'
#' @details
#' Because the mode is not a linear operator, two modal point forecasts are
#' available. \code{point = "joint"} (the default) returns the joint modal
#' trajectory: the most probable future path, obtained by iterating the model
#' with future innovations set to zero. \code{point = "marginal"} returns, for
#' each horizon \eqn{h}, the conditional mode of \eqn{y_{T+h}} itself. The two
#' coincide at \eqn{h = 1} and for symmetric errors (\eqn{\gamma = 1}); for
#' skewed errors and \eqn{h \ge 2} the marginal mode equals the joint path plus
#' the mode \eqn{\delta_h} of the forecast error
#' \eqn{\sum_{j=0}^{h-1}\psi_j\epsilon_{T+h-j}}, which can be substantial for
#' integrated series.
#'
#' With \code{interval = "asymptotic"}, the prediction intervals use the exact
#' distribution of the forecast error at the estimated parameters, computed by
#' numerical convolution of the SKD densities. With \code{interval = "bootstrap"},
#' they use empirical quantiles of \code{npaths} simulated future paths.
#'
#' @param object A modal_arima object.
#' @param h The forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param interval Method for computing prediction intervals ("asymptotic" or "bootstrap").
#' @param npaths Number of simulated paths for bootstrap intervals. Defaults to 1000.
#' @param point Modal point forecast: \code{"joint"} (most probable trajectory,
#'   default) or \code{"marginal"} (conditional mode at each horizon).
#' @param ... Additional arguments.
#' @return An object of class \code{forecast}. Besides the usual components, it
#'   contains \code{mode_path} (the joint modal trajectory) and \code{mode_shift}
#'   (the correction \eqn{\delta_h} from the joint path to the marginal mode).
#' @export
#'
#' @examples
#' library(forecast)
#'
#' # Non-seasonal
#' fit <- fit_modal_arima(log10(lynx), order = c(2, 0, 0))
#' pred <- forecast(fit, h = 5, level = c(80, 95))
#' autoplot(pred)
#' accuracy(pred)
#'
#' # Seasonal airline model: joint trajectory versus marginal mode
#' fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
#'                            seasonal = list(order = c(0, 1, 1), period = 12))
#' fc_joint <- forecast(fit_air, h = 24)
#' fc_marg <- forecast(fit_air, h = 24, point = "marginal")
#' round(cbind(joint = fc_joint$mean, marginal = fc_marg$mean)[c(1, 12, 24), ], 4)
#' fc_boot <- forecast(fit_air, h = 12, interval = "bootstrap", npaths = 200)
forecast.modal_arima <- function(object, h = 10, level = c(80, 95), interval = c("asymptotic", "bootstrap"),
                                 npaths = 1000, point = c("joint", "marginal"), ...) {
  interval <- match.arg(interval)
  point <- match.arg(point)
  st <- .modal_filter(object)
  pp <- st$pp; dist <- st$spec$dist

  path <- .modal_path(st, object$y, h)
  psi <- .psi_weights(st, h)
  err <- lapply(seq_len(h), function(k) .error_distribution(psi[seq_len(k)], pp, dist))
  shift <- vapply(err, function(e) e$mode, numeric(1))
  point_fc <- if (point == "joint") path else path + shift

  lower <- matrix(NA, nrow = h, ncol = length(level))
  upper <- matrix(NA, nrow = h, ncol = length(level))
  if (interval == "asymptotic") {
    for (i in seq_along(level)) {
      alpha <- 1 - level[i] / 100
      lower[, i] <- path + vapply(err, function(e) e$quantile(alpha / 2), numeric(1))
      upper[, i] <- path + vapply(err, function(e) e$quantile(1 - alpha / 2), numeric(1))
    }
  } else {
    E <- matrix(.rskd(h * npaths, 0, pp$sigma, pp$gamma, dist, pp$nu), h, npaths)
    sim_paths <- as.matrix(.integrate_path(.simulate_w(st, E), object$y, st$spec))
    for (i in seq_along(level)) {
      alpha <- 1 - level[i] / 100
      lower[, i] <- apply(sim_paths, 1, stats::quantile, probs = alpha / 2)
      upper[, i] <- apply(sim_paths, 1, stats::quantile, probs = 1 - alpha / 2)
    }
  }

  y_ts <- if (stats::is.ts(object$y)) object$y else stats::ts(object$y)
  tsp_y <- stats::tsp(y_ts)
  start_fc <- tsp_y[2] + 1 / tsp_y[3]
  as_fc_ts <- function(x) stats::ts(x, start = start_fc, frequency = tsp_y[3])
  colnames(lower) <- colnames(upper) <- paste0(level, "%")

  res <- list(mean = as_fc_ts(point_fc), lower = as_fc_ts(lower), upper = as_fc_ts(upper), level = level,
              method = paste0(.model_label(object), " [", .dist_label(dist), ", ", point, " mode]"),
              x = y_ts, model = object,
              fitted = stats::ts(as.numeric(object$fitted.values), start = tsp_y[1], frequency = tsp_y[3]),
              residuals = stats::ts(as.numeric(object$residuals), start = tsp_y[1], frequency = tsp_y[3]),
              mode_path = as_fc_ts(path), mode_shift = as_fc_ts(shift))
  class(res) <- "forecast"
  res
}
