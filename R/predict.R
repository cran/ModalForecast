#' @export
predict.modal_arima <- function(object, n.ahead = 10, ...) {
  y_orig <- object$y; order <- object$order; coefs <- object$coefficients
  p <- order[1]; d <- order[2]; q <- order[3]
  y_diff <- y_orig; if (d > 0) y_diff <- diff(y_orig, differences = d); n <- length(y_diff)
  c_mu <- coefs["intercept"]
  phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) for (i in 1:p) ar_term <- ar_term + phi[i] * (if (t-i > 0) y_diff[t-i] else mean_y)
    ma_term <- 0
    if (q > 0) for (j in 1:q) ma_term <- ma_term + theta[j] * (if (t-j > 0) eps[t-j] else 0)
    mu_t[t] <- c_mu + ar_term + ma_term; eps[t] <- y_diff[t] - mu_t[t]
  }
  pred_diff <- numeric(n.ahead)
  for (h in 1:n.ahead) {
    ar_term <- 0
    if (p > 0) for (i in 1:p) { val <- if (h-i > 0) pred_diff[h-i] else y_diff[n+h-i]; ar_term <- ar_term + phi[i]*val }
    ma_term <- 0
    if (q > 0) for (j in 1:q) { val <- if (h-j > 0) 0 else eps[n+h-j]; ma_term <- ma_term + theta[j]*val }
    pred_diff[h] <- c_mu + ar_term + ma_term
  }
  if (d > 0) { last_y <- utils::tail(y_orig, d); pred <- stats::diffinv(pred_diff, differences = d, xi = last_y)[-(1:d)] }
  else pred <- pred_diff
  return(pred)
}

#' @importFrom forecast forecast
#' @export
forecast::forecast

# ========== Internal distribution functions ==========

# --- SKN (normal) ---
.qsn <- function(prob, mu = 0, sigma = 1, gamma = 1) {
  p_mu <- 1/(gamma^2+1); res <- numeric(length(prob))
  for (i in seq_along(prob)) {
    if (prob[i] < p_mu) z <- (1/gamma)*stats::qnorm(prob[i]*(1+gamma^2)/2)
    else { val <- (prob[i]*(1+gamma^2)+gamma^2-1)/(2*gamma^2); z <- gamma*stats::qnorm(val) }
    res[i] <- mu + sigma[min(i,length(sigma))]*z
  }; return(res)
}
.rsn <- function(n, mu = 0, sigma = 1, gamma = 1) {
  u <- stats::runif(n); p_mu <- 1/(gamma^2+1); z <- numeric(n)
  idx_l <- which(u < p_mu); idx_u <- which(u >= p_mu)
  if (length(idx_l) > 0) z[idx_l] <- (1/gamma)*stats::qnorm(u[idx_l]*(1+gamma^2)/2)
  if (length(idx_u) > 0) { val <- (u[idx_u]*(1+gamma^2)+gamma^2-1)/(2*gamma^2); z[idx_u] <- gamma*stats::qnorm(val) }
  return(mu + sigma*z)
}

# --- SKT (t) ---
.dskt <- function(x, mu = 0, sigma = 1, gamma = 1, nu = 5) {
  p_skew <- 1/(gamma^2+1); z <- (x-mu)/sigma; rho <- z*(p_skew-(z<0))
  exp(lgamma((nu+1)/2)-lgamma(nu/2)-log(sigma)-0.5*log(pi*nu))*4*p_skew*(1-p_skew)*(4*rho^2/nu+1)^(-(nu+1)/2)
}
.pskt <- function(x, mu = 0, sigma = 1, gamma = 1, nu = 5) {
  cdf <- numeric(length(x))
  for (i in seq_along(x)) {
    cdf[i] <- tryCatch(stats::integrate(function(y) .dskt(y, mu=mu, sigma=sigma, gamma=gamma, nu=nu), lower=-Inf, upper=x[i], rel.tol=1e-8)$value,
                        error=function(e) 0.5)
  }; pmax(0, pmin(1, cdf))
}
.qskt <- function(prob, mu = 0, sigma = 1, gamma = 1, nu = 5) {
  res <- numeric(length(prob))
  for (i in seq_along(prob)) {
    sig_i <- sigma[min(i,length(sigma))]
    f_root <- function(x) .pskt(x, mu=mu, sigma=sig_i, gamma=gamma, nu=nu) - prob[i]
    res[i] <- tryCatch(stats::uniroot(f_root, c(mu-10*sig_i, mu+10*sig_i), tol=1e-8)$root,
                        error=function(e) tryCatch(stats::uniroot(f_root, c(mu-50*sig_i, mu+50*sig_i), tol=1e-8)$root, error=function(e2) mu))
  }; return(res)
}
.rskt <- function(n, mu = 0, sigma = 1, gamma = 1, nu = 5) {
  u_vals <- stats::rgamma(n, shape=nu/2, rate=nu/2)
  mu + sigma*(1/sqrt(u_vals))*.rsn(n, mu=0, sigma=1, gamma=gamma)
}

# --- SKL (laplace) ---
# CDF of Skewed Laplace: F(y) = p*(1-exp(-2*(1-p)*(mu-y)/sigma)) for y <= mu
#                                1-(1-p)*exp(-2*p*(y-mu)/sigma) for y > mu
.dskl <- function(x, mu = 0, sigma = 1, gamma = 1) {
  p_skew <- 1/(gamma^2+1); z <- (x-mu)/sigma; rho <- z*(p_skew-(z<0))
  2*p_skew*(1-p_skew)/sigma*exp(-2*rho)
}
.pskl <- function(x, mu = 0, sigma = 1, gamma = 1) {
  p_skew <- 1/(gamma^2+1); cdf <- numeric(length(x))
  mu_v <- rep_len(mu, length(x)); sig_v <- rep_len(sigma, length(x))
  idx_le <- which(x <= mu_v); idx_gt <- which(x > mu_v)
  if (length(idx_le) > 0) cdf[idx_le] <- p_skew*(1-exp(-2*(1-p_skew)*(mu_v[idx_le]-x[idx_le])/sig_v[idx_le]))
  if (length(idx_gt) > 0) cdf[idx_gt] <- 1-(1-p_skew)*exp(-2*p_skew*(x[idx_gt]-mu_v[idx_gt])/sig_v[idx_gt])
  pmax(0, pmin(1, cdf))
}
.qskl <- function(prob, mu = 0, sigma = 1, gamma = 1) {
  p_skew <- 1/(gamma^2+1); res <- numeric(length(prob))
  for (i in seq_along(prob)) {
    sig_i <- sigma[min(i,length(sigma))]
    if (prob[i] <= p_skew) {
      # Invert: prob = p*(1-exp(-2*(1-p)*(mu-x)/sigma))
      res[i] <- mu + sig_i*log(1-prob[i]/p_skew)/(2*(1-p_skew))
    } else {
      # Invert: prob = 1-(1-p)*exp(-2*p*(x-mu)/sigma)
      res[i] <- mu - sig_i*log((1-prob[i])/(1-p_skew))/(2*p_skew)
    }
  }; return(res)
}
.rskl <- function(n, mu = 0, sigma = 1, gamma = 1) {
  # Stochastic representation: Y = mu + sigma*sqrt(U)*Z where U ~ Exp(2), Z ~ SKN(0,1,p)
  u_vals <- stats::rexp(n, rate = 2)
  mu + sigma*sqrt(u_vals)*.rsn(n, mu=0, sigma=1, gamma=gamma)
}

#' Forecast methodology for Modal ARIMA
#'
#' @param object A modal_arima object.
#' @param h The forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param interval Method for computing prediction intervals ("asymptotic" or "bootstrap").
#' @param npaths Number of simulated paths for bootstrap intervals. Defaults to 1000.
#' @param ... Additional arguments.
#' @return A forecast object.
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
forecast.modal_arima <- function(object, h = 10, level = c(80, 95), interval = c("asymptotic", "bootstrap"), npaths = 1000, ...) {
  interval <- match.arg(interval)
  dist <- if (!is.null(object$dist)) object$dist else "normal"

  if (length(object$y) <= 50 && interval == "asymptotic")
    message("Notice: The sample size is small (n <= 50). Consider using interval = 'bootstrap'.")

  pred <- predict(object, n.ahead = h)
  coefs <- object$coefficients
  p <- object$order[1]; d <- object$order[2]; q <- object$order[3]
  ar_coef <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  ma_coef <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  sigma <- coefs["sigma"]; gamma <- coefs["gamma"]
  nu <- if (dist == "t") coefs["nu"] else NULL
  c_mu <- coefs["intercept"]

  lower <- matrix(NA, nrow=h, ncol=length(level)); upper <- matrix(NA, nrow=h, ncol=length(level))

  if (interval == "asymptotic") {
    poly_diff <- c(1)
    if (d > 0) for (i in seq_len(d)) poly_diff <- stats::convolve(poly_diff, rev(c(1,-1)), type="open")
    ar_poly <- c(1, -ar_coef)
    total_ar <- stats::convolve(ar_poly, rev(poly_diff), type="open")
    total_ar_coef <- -total_ar[-1]; total_ar_coef <- total_ar_coef[total_ar_coef != 0]
    psi_weights <- if (h > 1) c(1, stats::ARMAtoMA(ar=total_ar_coef, ma=ma_coef, lag.max=h-1)) else c(1)
    sigma_h <- sigma*sqrt(cumsum(psi_weights^2))

    for (i in seq_along(level)) {
      alpha <- 1-level[i]/100; p_lower <- alpha/2; p_upper <- 1-alpha/2
      if (dist == "normal") {
        q_low <- .qsn(rep(p_lower,h), mu=0, sigma=sigma_h, gamma=gamma)
        q_up  <- .qsn(rep(p_upper,h), mu=0, sigma=sigma_h, gamma=gamma)
      } else if (dist == "t") {
        q_low <- .qskt(rep(p_lower,h), mu=0, sigma=sigma_h, gamma=gamma, nu=nu)
        q_up  <- .qskt(rep(p_upper,h), mu=0, sigma=sigma_h, gamma=gamma, nu=nu)
      } else if (dist == "laplace") {
        q_low <- .qskl(rep(p_lower,h), mu=0, sigma=sigma_h, gamma=gamma)
        q_up  <- .qskl(rep(p_upper,h), mu=0, sigma=sigma_h, gamma=gamma)
      }
      lower[,i] <- pred + q_low; upper[,i] <- pred + q_up
    }
  } else {
    y_orig <- object$y; y_diff <- y_orig
    if (d > 0) y_diff <- diff(y_orig, differences=d); n_hist <- length(y_diff)
    mean_y <- mean(y_diff); eps_hist <- numeric(n_hist)
    for (t in 1:n_hist) {
      ar_term <- 0; if (p > 0) for (i in seq_len(p)) ar_term <- ar_term+ar_coef[i]*(if (t-i>0) y_diff[t-i] else mean_y)
      ma_term <- 0; if (q > 0) for (j in seq_len(q)) ma_term <- ma_term+ma_coef[j]*(if (t-j>0) eps_hist[t-j] else 0)
      eps_hist[t] <- y_diff[t]-(c_mu+ar_term+ma_term)
    }
    sim_paths <- matrix(NA, nrow=h, ncol=npaths)
    for (b in 1:npaths) {
      yt_sim <- numeric(h)
      eps_sim <- if (dist == "normal") .rsn(h, mu=0, sigma=sigma, gamma=gamma)
                 else if (dist == "t") .rskt(h, mu=0, sigma=sigma, gamma=gamma, nu=nu)
                 else .rskl(h, mu=0, sigma=sigma, gamma=gamma)
      for (t in 1:h) {
        ar_term <- 0; if (p > 0) for (i in 1:p) { val <- if (t-i>0) yt_sim[t-i] else y_diff[n_hist+t-i]; ar_term <- ar_term+ar_coef[i]*val }
        ma_term <- 0; if (q > 0) for (j in 1:q) { val <- if (t-j>0) eps_sim[t-j] else eps_hist[n_hist+t-j]; ma_term <- ma_term+ma_coef[j]*val }
        yt_sim[t] <- c_mu+ar_term+ma_term+eps_sim[t]
      }
      if (d > 0) { last_y <- utils::tail(y_orig,d); path <- stats::diffinv(yt_sim, differences=d, xi=last_y)[-(1:d)] }
      else path <- yt_sim
      sim_paths[,b] <- path
    }
    for (i in seq_along(level)) {
      alpha <- 1-level[i]/100
      lower[,i] <- apply(sim_paths, 1, stats::quantile, probs=alpha/2)
      upper[,i] <- apply(sim_paths, 1, stats::quantile, probs=1-alpha/2)
    }
  }

  dist_label <- switch(dist, "normal"="Skew-Normal", "t"="Skewed Student-t", "laplace"="Skewed Laplace")
  y_ts <- if(stats::is.ts(object$y)) object$y else stats::ts(object$y)
  tsp_y <- stats::tsp(y_ts)
  if (is.null(tsp_y)) tsp_y <- c(1, length(y_ts), 1)
  
  mean_ts <- stats::ts(pred, start = tsp_y[2] + 1/tsp_y[3], frequency = tsp_y[3])
  lower_ts <- stats::ts(lower, start = tsp_y[2] + 1/tsp_y[3], frequency = tsp_y[3])
  upper_ts <- stats::ts(upper, start = tsp_y[2] + 1/tsp_y[3], frequency = tsp_y[3])
  fitted_ts <- stats::ts(object$fitted.values, start = tsp_y[1], frequency = tsp_y[3])
  resid_ts <- stats::ts(object$residuals, start = tsp_y[1], frequency = tsp_y[3])
  
  res <- list(mean=mean_ts, lower=lower_ts, upper=upper_ts, level=level,
              method=paste0("Modal ARIMA(", p, ",", d, ",", q, ") [", dist_label, "]"),
              x=y_ts, model=object, fitted=fitted_ts, residuals=resid_ts)
  class(res) <- "forecast"; return(res)
}
