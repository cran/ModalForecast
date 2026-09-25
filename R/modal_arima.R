#' Fit a Parametric Modal ARIMA or Seasonal ARIMA Model using the SKD Family
#'
#' Fits a Modal ARIMA or Modal SARIMA model, where the conditional mode follows
#' an (optionally seasonal) ARIMA recursion and the innovations follow a member
#' of the SKD (Skewed Distribution) family: Skew-Normal, Skewed Student-t or
#' Skewed Laplace.
#'
#' @details
#' Let \eqn{w_t = (1-B)^d (1-B^s)^D y_t}. The model is
#' \deqn{\phi(B)\Phi(B^s) w_t = c + \theta(B)\Theta(B^s)\epsilon_t,}
#' where \eqn{\epsilon_t} are independent \eqn{SKD(0, \sigma, \gamma)} errors
#' with mode zero, so that \eqn{\mu_t = w_t - \epsilon_t} is the conditional
#' mode of \eqn{w_t}. The non-seasonal Modal ARIMA model is the special case
#' \eqn{P = D = Q = 0}.
#'
#' All members share the parameterization of Galarza et al. (2017),
#' \deqn{f(y \mid \mu, \sigma, \gamma) = \frac{4p(1-p)}{\sigma} g\{2\rho_p((y-\mu)/\sigma)\},}
#' with \eqn{p = 1/(1+\gamma^2)} and \eqn{\rho_p(u) = u(p - I(u<0))}, where
#' \eqn{g} is the standard normal, Student-t or Laplace density. The mode
#' \eqn{\mu} is also the \eqn{p}-th quantile, and \eqn{\gamma > 1} gives a
#' heavier right tail.
#'
#' Parameters are estimated by maximum likelihood with BFGS and analytical
#' gradients, started from a conditional-sum-of-squares Gaussian (S)ARIMA fit.
#' The Skewed Laplace log-likelihood is not differentiable when a residual is
#' zero, so for \code{dist = "laplace"} the BFGS solution is refined with
#' Nelder-Mead. Stationarity and invertibility of the regular and seasonal
#' polynomials are enforced during optimization.
#'
#' @references
#' Galarza, C. E., Lachos, V. H., Cabral, C. R. B., and Castro, L. M. (2017).
#' Robust quantile regression using a generalized class of skewed distributions.
#' Stat, 6(1), 113-130.
#'
#' @seealso \code{\link{auto.modal.arima}}, \code{\link{forecast.modal_arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations.
#' @param order A specification of the non-seasonal part of the model: the
#'   three components (p, d, q) are the AR order, the degree of differencing,
#'   and the MA order.
#' @param seasonal A specification of the seasonal part of the model, as in
#'   \code{\link[stats]{arima}}: a list with components \code{order}
#'   (the seasonal orders (P, D, Q)) and \code{period} (the seasonal period
#'   \eqn{s}; defaults to \code{frequency(y)}). A numeric vector of length 3 is
#'   taken as the seasonal order. The default is a non-seasonal model.
#' @param dist Character string specifying the error distribution from the
#'   SKD family. One of \code{"normal"} (default) for the Skew-Normal,
#'   \code{"t"} for the Skewed Student-t (adds degrees-of-freedom parameter
#'   \code{nu}), or \code{"laplace"} for the Skewed Laplace distribution.
#'
#' @return An object of class \code{modal_arima} containing:
#'   \describe{
#'     \item{\code{y}}{The original time series.}
#'     \item{\code{order}}{The non-seasonal order \code{(p,d,q)}.}
#'     \item{\code{seasonal}}{List with the seasonal \code{order} \code{(P,D,Q)} and \code{period}.}
#'     \item{\code{coefficients}}{Named vector of estimated parameters.}
#'     \item{\code{vcov}}{Asymptotic covariance matrix of the coefficients, from the observed information (delta method for \code{sigma}, \code{gamma} and \code{nu}).}
#'     \item{\code{loglik}}{The maximized log-likelihood.}
#'     \item{\code{nobs}}{Number of observations used in the likelihood, \eqn{T - d - Ds}.}
#'     \item{\code{hessian}}{The observed information on the optimization scale: the Hessian of the negative log-likelihood, or the outer product of the per-observation scores (OPG) for the Skewed Laplace, whose log-likelihood is not twice differentiable, and whenever the Hessian is not positive definite.}
#'     \item{\code{information}}{Which estimator was used for \code{hessian}: \code{"hessian"} or \code{"opg"}.}
#'     \item{\code{convergence}}{Convergence code from \code{optim}.}
#'     \item{\code{dist}}{The distribution used (\code{"normal"}, \code{"t"}, or \code{"laplace"}).}
#'     \item{\code{fitted.values}, \code{residuals}}{Fitted modes and modal residuals on the scale of \code{y} (\code{NA} for the first \eqn{d + Ds} observations).}
#'   }
#' @importFrom stats arima optim optimHess qnorm runif dnorm pnorm coef AIC BIC logLik printCoefmat
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' library(forecast)
#'
#' # Non-seasonal: Lynx data
#' y <- log10(lynx)
#' fit_n <- fit_modal_arima(y, order = c(2, 0, 0), dist = "normal")
#' fit_t <- fit_modal_arima(y, order = c(2, 0, 0), dist = "t")
#' fit_l <- fit_modal_arima(y, order = c(2, 0, 0), dist = "laplace")
#' c(Normal = AIC(fit_n), Student = AIC(fit_t), Laplace = AIC(fit_l))
#' summary(fit_n)
#'
#' # Seasonal: the airline model for monthly air passengers
#' fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
#'                            seasonal = list(order = c(0, 1, 1), period = 12))
#' summary(fit_air)
#' fc <- forecast(fit_air, h = 24)
#' autoplot(fc)
fit_modal_arima <- function(y, order = c(1, 0, 0),
                            seasonal = list(order = c(0, 0, 0), period = NA),
                            dist = c("normal", "t", "laplace")) {
  dist <- match.arg(dist)
  if (length(order) != 3) stop("'order' must have length 3 (p, d, q)")
  if (any(order < 0)) stop("'order' must be non-negative")
  if (anyNA(y)) stop("Missing values (NAs) are not currently supported in ModalForecast.")
  order <- as.integer(order)
  seasonal <- .parse_seasonal(seasonal, y)
  spec <- .make_spec(order, seasonal, dist)

  w <- .difference(y, spec$d, spec$D, spec$s)
  n <- length(w)
  n_coef <- 1 + spec$p + spec$q + spec$P + spec$Q + 2 + (dist == "t")
  if (n <= n_coef + 1) stop("Not enough observations after differencing to fit this model.")

  fn <- function(par) .neg_loglik(par, spec, w)
  gr <- function(par) .grad_neg_loglik(par, spec, w)
  opt <- .mle(w, spec, .initial_values(w, spec))
  info <- .information(opt$par, spec, w, fn, gr)
  hess <- info$matrix

  pp <- .split_par(opt$par, spec)
  coefficients <- c(opt$par[seq_len(1 + spec$p + spec$q + spec$P + spec$Q)],
                    pp$sigma, pp$gamma, pp$nu)
  names(coefficients) <- .coef_names(spec)

  # Delta method: the optimizer works with log(sigma), log(gamma) and log(nu - 2).
  jac <- rep(1, length(coefficients))
  jac[names(coefficients) == "sigma"] <- pp$sigma
  jac[names(coefficients) == "gamma"] <- pp$gamma
  if (dist == "t") jac[names(coefficients) == "nu"] <- pp$nu - 2
  # A parameter on the boundary (nu at its upper limit) has no information;
  # invert the block of the remaining parameters.
  free <- abs(diag(hess)) > 1e-8
  vcov <- matrix(NA_real_, length(jac), length(jac))
  vcov[free, free] <- tryCatch(solve(hess[free, free, drop = FALSE]),
                               error = function(e) matrix(NA_real_, sum(free), sum(free)))
  vcov <- vcov * outer(jac, jac)
  dimnames(vcov) <- list(names(coefficients), names(coefficients))

  rec <- .recursion(pp, spec, w)
  n_lost <- length(y) - n
  eps_vals <- c(rep(NA, n_lost), rec$eps)
  fit_vals <- c(rep(NA, n_lost), as.numeric(y)[n_lost + seq_len(n)] - rec$eps)
  if (stats::is.ts(y)) {
    eps_vals <- stats::ts(eps_vals, start = stats::start(y), frequency = stats::frequency(y))
    fit_vals <- stats::ts(fit_vals, start = stats::start(y), frequency = stats::frequency(y))
  }

  out <- list(
    y = y,
    order = order,
    seasonal = seasonal,
    coefficients = coefficients,
    vcov = vcov,
    loglik = -opt$value,
    nobs = n,
    hessian = hess,
    information = info$type,
    convergence = opt$convergence,
    dist = dist,
    fitted.values = fit_vals,
    residuals = eps_vals
  )
  class(out) <- "modal_arima"
  out
}

# Maximize the log-likelihood of the differenced series w from 'init': BFGS with
# analytical gradients, refined by Nelder-Mead for the non-smooth Laplace case.
.mle <- function(w, spec, init, reltol = 1e-12) {
  fn <- function(par) .neg_loglik(par, spec, w)
  gr <- function(par) .grad_neg_loglik(par, spec, w)
  opt <- stats::optim(init, fn, gr, method = "BFGS", control = list(maxit = 1000))
  if (spec$dist == "laplace") {
    nm <- stats::optim(opt$par, fn, method = "Nelder-Mead",
                       control = list(maxit = 20000, reltol = reltol))
    if (nm$value < opt$value) opt <- nm
  }
  opt
}

# Starting values from a conditional-sum-of-squares Gaussian (S)ARIMA fit of w.
.initial_values <- function(w, spec) {
  seas <- if (spec$P + spec$Q > 0) list(order = c(spec$P, 0, spec$Q), period = spec$s)
          else list(order = c(0, 0, 0))
  init_fit <- tryCatch(
    suppressWarnings(stats::arima(w, order = c(spec$p, 0, spec$q), seasonal = seas, method = "CSS")),
    error = function(e) NULL)
  get <- function(prefix, k) {
    if (k == 0) return(numeric(0))
    v <- if (!is.null(init_fit)) init_fit$coef[paste0(prefix, seq_len(k))] else rep(0, k)
    v[!is.finite(v)] <- 0
    as.numeric(v)
  }
  phi <- get("ar", spec$p); theta <- get("ma", spec$q)
  Phi <- get("sar", spec$P); Theta <- get("sma", spec$Q)
  pp <- list(phi = phi, theta = theta, Phi = Phi, Theta = Theta)
  if (!.admissible(pp)) {
    phi[] <- 0; theta[] <- 0; Phi[] <- 0; Theta[] <- 0
    pp <- list(phi = phi, theta = theta, Phi = Phi, Theta = Theta)
  }
  # arima() reports the mean; the model constant is c = mean * phi(1) * Phi(1).
  mean_w <- if (!is.null(init_fit) && "intercept" %in% names(init_fit$coef))
    init_fit$coef[["intercept"]] else mean(w)
  c0 <- mean_w * sum(.lag_polys(pp, spec)$ar)
  sigma0 <- if (!is.null(init_fit)) sqrt(init_fit$sigma2) else stats::sd(w)
  if (!is.finite(sigma0) || sigma0 < 1e-8) sigma0 <- 1
  par <- c(c0, phi, theta, Phi, Theta, log(sigma0), 0)
  if (spec$dist == "t") par <- c(par, log(3))
  unname(par)
}
