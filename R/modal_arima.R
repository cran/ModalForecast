#' Fit a Parametric Modal ARIMA Model using the SKD Family
#'
#' Fits a Modal ARIMA model where the conditional mode follows an ARIMA recursion
#' and the innovations follow a member of the SKD (Skewed Distribution) family.
#' Supports the Skew-Normal, Skewed Student-t, and Skewed Laplace distributions.
#'
#' @references
#' Galarza, C. E., Lachos, V. H., Cabral, C. R. B., and Castro, L. M. (2017).
#' Robust quantile regression using a generalized class of skewed distributions.
#' Stat, 6(1), 113-130.
#'
#' @seealso \code{\link{auto.modal.arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations.
#' @param order A specification of the non-seasonal part of the ARIMA
#'   model: the three components (p, d, q) are the AR order, the
#'   degree of differencing, and the MA order.
#' @param dist Character string specifying the error distribution from the
#'   SKD family. One of \code{"normal"} (default) for the Skew-Normal,
#'   \code{"t"} for the Skewed Student-t (adds degrees-of-freedom parameter
#'   \code{nu}), or \code{"laplace"} for the Skewed Laplace distribution.
#'
#' @return An object of class \code{modal_arima} containing:
#'   \describe{
#'     \item{\code{y}}{The original time series.}
#'     \item{\code{order}}{The ARIMA order \code{(p,d,q)}.}
#'     \item{\code{coefficients}}{Named vector of estimated parameters.}
#'     \item{\code{loglik}}{The maximized log-likelihood.}
#'     \item{\code{hessian}}{The Hessian matrix at the optimum.}
#'     \item{\code{convergence}}{Convergence code from \code{optim}.}
#'     \item{\code{dist}}{The distribution used (\code{"normal"}, \code{"t"}, or \code{"laplace"}).}
#'   }
#' @importFrom stats arima optim qnorm runif dnorm pnorm coef AIC BIC logLik printCoefmat
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' library(forecast)
#' 
#' # 1. Load Empirical Data (Lynx)
#' data(lynx)
#' y <- log10(lynx)
#'
#' # 2. Find the best SKD Error Distribution (Normal vs T vs Laplace) 
#' fit_n <- fit_modal_arima(y, order = c(2, 0, 0), dist = "normal")
#' fit_t <- fit_modal_arima(y, order = c(2, 0, 0), dist = "t")
#' fit_l <- fit_modal_arima(y, order = c(2, 0, 0), dist = "laplace")
#' c(Normal = AIC(fit_n), Student = AIC(fit_t), Laplace = AIC(fit_l))
#'
#' # 3. Auto Model Selection globally on the winning distribution (Skew-Normal)
#' fit_auto <- auto.modal.arima(y, d=0, max.p=2, max.q=2, dist="normal")
#'
#' # 4. Summary & Inferences
#' summary(fit_auto)
#'
#' # 5. Run residual diagnostics and Envelopes
#' diagnostics(fit_auto)
#' envelope(fit_auto, B=10)
#'
#' # 6. Produce forecasts with multiple prediction bands (alphas)
#' pred <- forecast(fit_auto, h=5, level = c(80, 95))
#'
#' # 7. Native integration with 'forecast' ecosystem
#' autoplot(pred)    
#' accuracy(pred)
fit_modal_arima <- function(y, order = c(1, 0, 0), dist = c("normal", "t", "laplace")) {
  dist <- match.arg(dist)
  if (length(order) != 3) stop("'order' must have length 3 (p, d, q)")
  if (anyNA(y)) stop("Missing values (NAs) are not currently supported in ModalForecast.")

  p <- order[1]
  d <- order[2]
  q <- order[3]

  y_orig <- y
  if (d > 0) {
    y <- diff(y, differences = d)
  }

  n <- length(y)

  # ==== Helper: ARIMA recursion (shared) ====
  arima_recursion <- function(params, p, q, n, y) {
    c_mu <- params[1]
    phi <- if (p > 0) params[2:(p+1)] else numeric(0)
    theta <- if (q > 0) params[(p+2):(1+p+q)] else numeric(0)

    if (p > 0) {
      ar_roots <- polyroot(c(1, -phi))
      if (any(abs(ar_roots) <= 1.001)) return(NULL)
    }
    if (q > 0) {
      ma_roots <- polyroot(c(1, theta))
      if (any(abs(ma_roots) <= 1.001)) return(NULL)
    }

    mu_t <- numeric(n)
    eps <- numeric(n)
    mean_y <- mean(y)

    for (t in 1:n) {
      ar_term <- 0
      for (i in seq_len(p)) {
        if (t - i > 0) ar_term <- ar_term + phi[i] * y[t - i]
        else ar_term <- ar_term + phi[i] * mean_y
      }
      ma_term <- 0
      for (j in seq_len(q)) {
        if (t - j > 0) ma_term <- ma_term + theta[j] * eps[t - j]
      }
      mu_t[t] <- c_mu + ar_term + ma_term
      eps[t] <- y[t] - mu_t[t]
    }
    return(list(mu_t = mu_t, eps = eps, phi = phi, theta = theta))
  }

  # ========================================================================
  # NORMAL (Skew-Normal) log-likelihood
  # ========================================================================
  neg_log_lik_normal <- function(params) {
    log_sigma <- params[length(params) - 1]
    log_gamma <- params[length(params)]
    sigma <- exp(log_sigma)
    gamma <- exp(log_gamma)

    rec <- arima_recursion(params, p, q, n, y)
    if (is.null(rec)) return(1e10)

    z <- (y - rec$mu_t) / sigma
    log_f0 <- function(val) -0.5 * log(2 * pi) - 0.5 * val^2

    idx_ge <- which(y >= rec$mu_t)
    idx_lt <- which(y < rec$mu_t)

    log_lik <- numeric(n)
    if (length(idx_ge) > 0) log_lik[idx_ge] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_ge] / gamma)
    if (length(idx_lt) > 0) log_lik[idx_lt] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_lt] * gamma)

    ans <- -sum(log_lik)
    if (is.na(ans) || is.infinite(ans)) ans <- 1e10
    return(ans)
  }

  grad_normal <- function(params) {
    log_sigma <- params[length(params) - 1]
    log_gamma <- params[length(params)]
    sigma <- exp(log_sigma)
    gamma <- exp(log_gamma)
    num_beta <- 1 + p + q
    G <- numeric(length(params))
    mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y)
    V <- matrix(0, nrow = n, ncol = num_beta)
    phi <- if (p > 0) params[2:(p+1)] else numeric(0)
    theta <- if (q > 0) params[(p+2):(1+p+q)] else numeric(0)

    for (t in 1:n) {
      ar_term <- 0
      if (p > 0) for (i in seq_len(p)) {
        if (t-i > 0) ar_term <- ar_term + phi[i]*y[t-i]
        else ar_term <- ar_term + phi[i]*mean_y
      }
      ma_term <- 0
      if (q > 0) for (j in seq_len(q)) {
        if (t-j > 0) ma_term <- ma_term + theta[j]*eps[t-j]
      }
      mu_t[t] <- params[1] + ar_term + ma_term
      eps[t] <- y[t] - mu_t[t]

      Z_t <- numeric(num_beta); Z_t[1] <- 1
      if (p > 0) for (i in 1:p) Z_t[1+i] <- if (t-i > 0) y[t-i] else mean_y
      if (q > 0) for (j in 1:q) Z_t[1+p+j] <- if (t-j > 0) eps[t-j] else 0
      V_ma <- numeric(num_beta)
      if (q > 0) for (j in 1:q) if (t-j > 0) V_ma <- V_ma + theta[j]*V[t-j,]
      V[t,] <- Z_t - V_ma

      e_t <- eps[t]
      if (e_t >= 0) {
        W_t <- 1/(sigma^2*gamma^2); d_log_gamma <- -(gamma^2-1)/(gamma^2+1) + W_t*e_t^2
      } else {
        W_t <- gamma^2/sigma^2; d_log_gamma <- -(gamma^2-1)/(gamma^2+1) - W_t*e_t^2
      }
      G[1:num_beta] <- G[1:num_beta] - W_t*e_t*V[t,]
      G[length(params)-1] <- G[length(params)-1] + (1 - W_t*e_t^2)
      G[length(params)] <- G[length(params)] - d_log_gamma
    }
    return(G)
  }

  # ========================================================================
  # T (Skewed Student-t) log-likelihood
  # ========================================================================
  neg_log_lik_t <- function(params) {
    log_sigma <- params[length(params) - 2]
    log_gamma <- params[length(params) - 1]
    log_nu    <- params[length(params)]
    sigma <- exp(log_sigma); gamma <- exp(log_gamma)
    nu <- min(exp(log_nu), 200) + 2

    rec <- arima_recursion(params, p, q, n, y)
    if (is.null(rec)) return(1e10)

    p_skew <- 1 / (gamma^2 + 1)
    z <- (y - rec$mu_t) / sigma
    rho <- z * (p_skew - (z < 0))

    log_lik <- lgamma((nu+1)/2) - lgamma(nu/2) - log(sigma) - 0.5*log(pi*nu) +
      log(4*p_skew*(1-p_skew)) - ((nu+1)/2)*log(4*rho^2/nu + 1)

    ans <- -sum(log_lik)
    if (is.na(ans) || is.infinite(ans)) ans <- 1e10
    return(ans)
  }

  grad_t <- function(params) {
    log_sigma <- params[length(params)-2]
    log_gamma <- params[length(params)-1]
    log_nu    <- params[length(params)]
    sigma <- exp(log_sigma); gamma <- exp(log_gamma)
    nu <- min(exp(log_nu), 200) + 2
    num_beta <- 1+p+q; G <- numeric(length(params))
    mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y)
    V <- matrix(0, nrow=n, ncol=num_beta)
    phi <- if (p>0) params[2:(p+1)] else numeric(0)
    theta <- if (q>0) params[(p+2):(1+p+q)] else numeric(0)
    p_skew <- 1/(gamma^2+1)

    for (t in 1:n) {
      ar_term <- 0
      if (p>0) for (i in seq_len(p)) {
        if (t-i>0) ar_term <- ar_term+phi[i]*y[t-i] else ar_term <- ar_term+phi[i]*mean_y
      }
      ma_term <- 0
      if (q>0) for (j in seq_len(q)) if (t-j>0) ma_term <- ma_term+theta[j]*eps[t-j]
      mu_t[t] <- params[1]+ar_term+ma_term; eps[t] <- y[t]-mu_t[t]

      Z_t <- numeric(num_beta); Z_t[1] <- 1
      if (p>0) for (i in 1:p) Z_t[1+i] <- if (t-i>0) y[t-i] else mean_y
      if (q>0) for (j in 1:q) Z_t[1+p+j] <- if (t-j>0) eps[t-j] else 0
      V_ma <- numeric(num_beta)
      if (q>0) for (j in 1:q) if (t-j>0) V_ma <- V_ma+theta[j]*V[t-j,]
      V[t,] <- Z_t-V_ma

      e_t <- eps[t]; z_t <- e_t/sigma
      xi_t <- if (e_t>=0) p_skew else (1-p_skew)
      rho_t <- z_t*(p_skew-(z_t<0)); denom <- nu+4*rho_t^2
      w_t <- (nu+1)/denom

      if (e_t>=0) dmu <- -w_t*8*p_skew^2*z_t/(sigma*nu)
      else dmu <- -w_t*8*(1-p_skew)^2*z_t/(sigma*nu)
      G[1:num_beta] <- G[1:num_beta]+dmu*V[t,]
      G[length(params)-2] <- G[length(params)-2]+(1-w_t*8*xi_t^2*z_t^2/nu)

      dp_dg <- -2*gamma/(gamma^2+1)^2
      d_log_p <- dp_dg/p_skew+(-dp_dg)/(1-p_skew)
      drho_dp <- z_t
      d_kernel <- -w_t*8*rho_t*drho_dp*dp_dg/nu
      G[length(params)-1] <- G[length(params)-1]-(d_log_p*gamma+d_kernel*gamma)

      dnu <- exp(log_nu)
      d_lgamma_1 <- 0.5*digamma((nu+1)/2)-0.5*digamma(nu/2) - 0.5/nu
      d_kernel_nu <- -0.5*log(4*rho_t^2/nu+1)+(nu+1)*4*rho_t^2/(2*nu^2*(4*rho_t^2/nu+1))
      G[length(params)] <- G[length(params)]-(d_lgamma_1+d_kernel_nu)*dnu
    }
    return(G)
  }

  # ========================================================================
  # LAPLACE (Skewed Laplace) log-likelihood
  # f(y|mu,sigma,p) = 2*p*(1-p)/sigma * exp(-2*rho_p((y-mu)/sigma))
  # ========================================================================
  neg_log_lik_laplace <- function(params) {
    log_sigma <- params[length(params) - 1]
    log_gamma <- params[length(params)]
    sigma <- exp(log_sigma); gamma <- exp(log_gamma)

    rec <- arima_recursion(params, p, q, n, y)
    if (is.null(rec)) return(1e10)

    p_skew <- 1 / (gamma^2 + 1)
    z <- (y - rec$mu_t) / sigma
    rho <- z * (p_skew - (z < 0))

    log_lik <- log(2*p_skew*(1-p_skew)) - log(sigma) - 2*rho

    ans <- -sum(log_lik)
    if (is.na(ans) || is.infinite(ans)) ans <- 1e10
    return(ans)
  }

  grad_laplace <- function(params) {
    log_sigma <- params[length(params)-1]
    log_gamma <- params[length(params)]
    sigma <- exp(log_sigma); gamma <- exp(log_gamma)
    num_beta <- 1+p+q; G <- numeric(length(params))
    mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y)
    V <- matrix(0, nrow=n, ncol=num_beta)
    phi <- if (p>0) params[2:(p+1)] else numeric(0)
    theta <- if (q>0) params[(p+2):(1+p+q)] else numeric(0)
    p_skew <- 1/(gamma^2+1)

    for (t in 1:n) {
      ar_term <- 0
      if (p>0) for (i in seq_len(p)) {
        if (t-i>0) ar_term <- ar_term+phi[i]*y[t-i] else ar_term <- ar_term+phi[i]*mean_y
      }
      ma_term <- 0
      if (q>0) for (j in seq_len(q)) if (t-j>0) ma_term <- ma_term+theta[j]*eps[t-j]
      mu_t[t] <- params[1]+ar_term+ma_term; eps[t] <- y[t]-mu_t[t]

      Z_t <- numeric(num_beta); Z_t[1] <- 1
      if (p>0) for (i in 1:p) Z_t[1+i] <- if (t-i>0) y[t-i] else mean_y
      if (q>0) for (j in 1:q) Z_t[1+p+j] <- if (t-j>0) eps[t-j] else 0
      V_ma <- numeric(num_beta)
      if (q>0) for (j in 1:q) if (t-j>0) V_ma <- V_ma+theta[j]*V[t-j,]
      V[t,] <- Z_t-V_ma

      e_t <- eps[t]; z_t <- e_t/sigma
      # d/d_mu of -2*rho_p(z) where rho_p(z) = z*(p - I(z<0))
      # d(rho)/d(mu) = -(1/sigma)*(p - I(z<0))
      sign_contrib <- if (e_t >= 0) p_skew else -(1-p_skew)
      # grad w.r.t. beta: increases neg_log_lik
      G[1:num_beta] <- G[1:num_beta] - (2/sigma)*sign_contrib*V[t,]

      # grad w.r.t. log_sigma: d/d_log_sigma = 1 - 2*|rho_p(z)| (chain rule)
      rho_t <- z_t*(p_skew-(z_t<0))
      G[length(params)-1] <- G[length(params)-1] + (1 - 2*rho_t)

      # grad w.r.t. log_gamma
      dp_dg <- -2*gamma/(gamma^2+1)^2
      d_log_p <- dp_dg/p_skew + (-dp_dg)/(1-p_skew)
      drho_dp <- z_t
      G[length(params)] <- G[length(params)] - (d_log_p + 2*drho_dp*dp_dg)*gamma
    }
    return(G)
  }

  # ========================================================================
  # Initial values based on standard ARIMA
  # ========================================================================
  init_fit <- suppressWarnings(stats::arima(y, order = c(p, 0, q), method = "CSS"))
  init_c <- if ("intercept" %in% names(init_fit$coef)) init_fit$coef["intercept"] else 0
  init_phi <- if (p > 0) init_fit$coef[paste0("ar", 1:p)] else numeric(0)
  init_theta <- if (q > 0) init_fit$coef[paste0("ma", 1:q)] else numeric(0)
  init_sigma <- sqrt(init_fit$sigma2)
  if (init_sigma < 1e-4) init_sigma <- 1

  if (dist == "normal") {
    init_params <- c(init_c, init_phi, init_theta, log(init_sigma), 0)
    opt <- optim(par = init_params, fn = neg_log_lik_normal, gr = grad_normal,
                 method = "BFGS", hessian = TRUE)
    est_sigma <- exp(opt$par[length(opt$par)-1])
    est_gamma <- exp(opt$par[length(opt$par)])
    coef_names <- c("intercept")
    if (p > 0) coef_names <- c(coef_names, paste0("ar", 1:p))
    if (q > 0) coef_names <- c(coef_names, paste0("ma", 1:q))
    coef_names <- c(coef_names, "sigma", "gamma")
    coefficients <- c(opt$par[1],
                      if (p>0) opt$par[2:(p+1)] else NULL,
                      if (q>0) opt$par[(p+2):(1+p+q)] else NULL,
                      est_sigma, est_gamma)
    names(coefficients) <- coef_names

  } else if (dist == "t") {
    init_params <- c(init_c, init_phi, init_theta, log(init_sigma), 0, log(3))
    opt <- optim(par = init_params, fn = neg_log_lik_t, gr = grad_t,
                 method = "BFGS", hessian = TRUE)
    est_sigma <- exp(opt$par[length(opt$par)-2])
    est_gamma <- exp(opt$par[length(opt$par)-1])
    est_nu    <- min(exp(opt$par[length(opt$par)]), 200) + 2
    coef_names <- c("intercept")
    if (p > 0) coef_names <- c(coef_names, paste0("ar", 1:p))
    if (q > 0) coef_names <- c(coef_names, paste0("ma", 1:q))
    coef_names <- c(coef_names, "sigma", "gamma", "nu")
    coefficients <- c(opt$par[1],
                      if (p>0) opt$par[2:(p+1)] else NULL,
                      if (q>0) opt$par[(p+2):(1+p+q)] else NULL,
                      est_sigma, est_gamma, est_nu)
    names(coefficients) <- coef_names

  } else if (dist == "laplace") {
    init_params <- c(init_c, init_phi, init_theta, log(init_sigma), 0)
    opt <- optim(par = init_params, fn = neg_log_lik_laplace, gr = grad_laplace,
                 method = "BFGS", hessian = TRUE)
    est_sigma <- exp(opt$par[length(opt$par)-1])
    est_gamma <- exp(opt$par[length(opt$par)])
    coef_names <- c("intercept")
    if (p > 0) coef_names <- c(coef_names, paste0("ar", 1:p))
    if (q > 0) coef_names <- c(coef_names, paste0("ma", 1:q))
    coef_names <- c(coef_names, "sigma", "gamma")
    coefficients <- c(opt$par[1],
                      if (p>0) opt$par[2:(p+1)] else NULL,
                      if (q>0) opt$par[(p+2):(1+p+q)] else NULL,
                      est_sigma, est_gamma)
    names(coefficients) <- coef_names
  }

  rec <- arima_recursion(opt$par, p, q, n, y)
  eps_vals <- rep(NA, length(y_orig))
  fit_vals <- rep(NA, length(y_orig))
  
  if (!is.null(rec)) {
     if (d > 0) {
        eps_vals[(d+1):length(y_orig)] <- rec$eps
        fit_vals[(d+1):length(y_orig)] <- y_orig[(d+1):length(y_orig)] - rec$eps
     } else {
        eps_vals <- rec$eps
        fit_vals <- rec$mu_t
     }
  }
  if (stats::is.ts(y_orig)) {
     eps_vals <- stats::ts(eps_vals, start=stats::start(y_orig), frequency=stats::frequency(y_orig))
     fit_vals <- stats::ts(fit_vals, start=stats::start(y_orig), frequency=stats::frequency(y_orig))
  }

  out <- list(
    y = y_orig,
    order = order,
    coefficients = coefficients,
    loglik = -opt$value,
    hessian = opt$hessian,
    convergence = opt$convergence,
    dist = dist,
    fitted.values = fit_vals,
    residuals = eps_vals
  )
  class(out) <- "modal_arima"
  return(out)
}
