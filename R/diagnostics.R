#' Diagnostic Plots for Modal ARIMA Models
#'
#' @description
#' Provides visual and statistical diagnostics for the residuals of a fitted modal ARIMA model.
#' Produces a comprehensive diagnostic panel including time series plot with
#' fitted modes, ACF/PACF of residuals, QQ-plot for normality, histogram of
#' residuals, and Ljung-Box p-values, all implemented using ggplot2.
#'
#' @param object An object of class \code{modal_arima}.
#' @param ... Additional arguments (unused).
#' @return A list of ggplot objects (invisibly) and draws the panel.
#' @import ggplot2
#' @importFrom stats acf pacf shapiro.test Box.test qqnorm qqline predict is.ts ts start frequency dnorm qnorm pnorm time residuals rgamma integrate uniroot
#' @importFrom scales comma_format
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra grid.arrange
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
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.modal_arima <- function(object, ...) {
  # Global variable hack for ggplot2
  Time <- Observed <- Fitted <- Lag <- ACF <- PACF <- Residuals <- Pval <- density <- NULL

  y_orig <- object$y
  order <- object$order
  coefs <- object$coefficients
  dist <- if (!is.null(object$dist)) object$dist else "normal"

  p <- order[1]
  d <- order[2]
  q <- order[3]

  y_diff <- y_orig
  if (d > 0) y_diff <- diff(y_orig, differences = d)
  n <- length(y_diff)

  c_mu <- coefs["intercept"]
  phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  sigma_hat <- coefs["sigma"]
  gamma_hat <- coefs["gamma"]
  nu_hat <- if (dist == "t") coefs["nu"] else NULL

  # Reconstruct fitted modes and residuals
  mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0; if (p>0) for (i in 1:p) ar_term <- ar_term+phi[i]*(if (t-i>0) y_diff[t-i] else mean_y)
    ma_term <- 0; if (q>0) for (j in 1:q) ma_term <- ma_term+theta[j]*(if (t-j>0) eps[t-j] else 0)
    mu_t[t] <- c_mu+ar_term+ma_term; eps[t] <- y_diff[t]-mu_t[t]
  }

  # Randomized Quantile Residuals via CDF
  if (dist == "normal") {
    u <- p_skewnorm(y_diff, mu_t, sigma_hat, gamma_hat)
  } else if (dist == "t") {
    u <- numeric(n); for (i in 1:n) u[i] <- .pskt(y_diff[i], mu=mu_t[i], sigma=sigma_hat, gamma=gamma_hat, nu=nu_hat)
  } else if (dist == "laplace") {
    u <- .pskl(y_diff, mu=mu_t, sigma=sigma_hat, gamma=gamma_hat)
  }
  u <- pmax(1e-7, pmin(1-1e-7, u)); rqr <- stats::qnorm(u)

  time_idx <- if(stats::is.ts(y_diff)) as.numeric(stats::time(y_diff)) else 1:n
  df_fits <- data.frame(Time=time_idx, Observed=as.numeric(y_diff), Fitted=as.numeric(mu_t))

  p1 <- ggplot(df_fits, aes(x=Time)) + geom_line(aes(y=Observed, color="Observed"), alpha=0.6) +
    geom_line(aes(y=Fitted, color="Fitted Mode"), linewidth=0.75) +
    scale_color_manual(values=c("Observed"="gray40", "Fitted Mode"="blue")) +
    labs(title="Series & Fitted Modes", y="Value", x="Time", color=NULL) + theme_minimal(base_size=7) + theme(legend.position="top")

  bacf <- stats::acf(rqr, plot=FALSE, lag.max=20)
  df_acf <- data.frame(Lag=as.numeric(bacf$lag[-1]), ACF=as.numeric(bacf$acf[-1]))
  p2 <- ggplot(df_acf, aes(x=Lag, y=ACF)) + geom_bar(stat="identity", width=0.1, fill="blue") +
    geom_hline(yintercept=c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype="dashed", color="red") + labs(title="ACF of Residuals (RQR)") + theme_minimal(base_size=7)

  bpacf <- stats::pacf(rqr, plot=FALSE, lag.max=20)
  df_pacf <- data.frame(Lag=as.numeric(bpacf$lag), PACF=as.numeric(bpacf$acf))
  p3 <- ggplot(df_pacf, aes(x=Lag, y=PACF)) + geom_bar(stat="identity", width=0.1, fill="blue") +
    geom_hline(yintercept=c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype="dashed", color="red") + labs(title="PACF of Residuals (RQR)") + theme_minimal(base_size=7)

  df_rqr <- data.frame(Residuals=rqr)
  p4 <- ggplot(df_rqr, aes(sample=Residuals)) + stat_qq(color="steelblue", alpha=0.6) + stat_qq_line(color="red", linewidth=0.8) +
    labs(title="Normal QQ-Plot (RQR)") + theme_minimal(base_size=7)

  p5 <- ggplot(df_rqr, aes(x=Residuals)) + geom_histogram(aes(y=after_stat(density)), bins=20, fill="steelblue", alpha=0.4, color="white") +
    stat_function(fun=stats::dnorm, color="red", linewidth=0.8) + labs(title="Histogram (RQR)") + theme_minimal(base_size=7)

  max_lag <- min(20, n-1)
  lb_pvals <- sapply(1:max_lag, function(lag) stats::Box.test(rqr, lag=lag, type="Ljung-Box")$p.value)
  df_lb <- data.frame(Lag=1:max_lag, Pval=lb_pvals)
  p6 <- ggplot(df_lb, aes(x=Lag, y=Pval)) + geom_point(color="blue", size=2) +
    geom_hline(yintercept=0.05, linetype="dashed", color="red", linewidth=0.8) + ylim(0,1) + labs(title="Ljung-Box p-values", y="p-value") + theme_minimal(base_size=7)

  dist_label <- switch(dist, "normal"="[Normal]", "t"="[Student-t]", "laplace"="[Laplace]")
  title_grob <- grid::textGrob(paste("Modal ARIMA", dist_label, "Diagnostic Panel"), gp=grid::gpar(fontsize=10, fontface="bold"))
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, top=title_grob)

  cat("\n=== Diagnostic Tests ===\n")
  cat(sprintf("Distribution: %s\n", dist))
  sw <- stats::shapiro.test(rqr)
  cat(sprintf("Shapiro-Wilk Normality Test: W = %.4f, p-value = %.4f\n", sw$statistic, sw$p.value))
  if (sw$p.value > 0.05) cat("  -> Normality assumption met (Expected/Good fit).\n") else cat("  -> Normality assumption rejected (Requires investigation).\n")
  
  lb <- stats::Box.test(rqr, lag = 10, type = "Ljung-Box")
  cat(sprintf("Ljung-Box Test (lag=10):     X2 = %.4f, p-value = %.4f\n", lb$statistic, lb$p.value))
  if (lb$p.value > 0.05) cat("  -> Residuals are independent / No Autocorrelation (Expected/Good fit).\n") else cat("  -> Residuals are autocorrelated (Investigate lags).\n")
  cat(sprintf("\nEstimated gamma (skewness):  %.4f\n", gamma_hat))
  gamma_se <- tryCatch(sqrt(diag(solve(object$hessian))[length(object$coefficients) - if(dist=="t") 1 else 0]), error=function(e) NA)
  if (!is.na(gamma_se)) {
      p_gamma <- 2 * (1 - stats::pnorm(abs((gamma_hat - 1) / gamma_se)))
      if (p_gamma > 0.05) cat("  -> Symmetry confirmed (gamma = 1) (Expected if true series is symmetric).\n\n") else cat("  -> Symmetry rejected / SKD Skewness present (Expected by model).\n\n")
  }
  if (dist == "t") cat(sprintf("Estimated nu (d.f.):        %.4f\n\n", nu_hat))

  invisible(list(p1, p2, p3, p4, p5, p6))
}

p_skewnorm <- function(y, mu, sigma, gamma) {
  z <- (y-mu)/sigma; p_thr <- 1/(gamma^2+1); cdf <- numeric(length(y))
  idx_lt <- which(y < mu); if (length(idx_lt) > 0) cdf[idx_lt] <- 2*p_thr*stats::pnorm(z[idx_lt]*gamma)
  idx_ge <- which(y >= mu); if (length(idx_ge) > 0) cdf[idx_ge] <- p_thr+2*(1-p_thr)*(stats::pnorm(z[idx_ge]/gamma)-0.5)
  return(cdf)
}

#' Simulation Envelope Diagnostics for Modal ARIMA Models
#'
#' Constructs simulation envelopes based on the theoretical distance
#' distributions from the SKD family (Galarza et al., 2017).
#'
#' @param object An object of class \code{modal_arima}.
#' @param B Number of Monte Carlo replications for envelope construction. Default is 100.
#' @param ... Additional arguments (unused).
#' @return A ggplot object (invisibly) and draws the envelope plot.
#' @import ggplot2
#' @export
envelope <- function(object, ...) { UseMethod("envelope") }

#' @rdname envelope
#' @export
envelope.modal_arima <- function(object, B = 100, ...) {
  Theoretical <- Observed <- Lower <- Upper <- Median <- NULL # ggplot2 hack
  coefs <- object$coefficients; dist <- if (!is.null(object$dist)) object$dist else "normal"
  y_orig <- object$y; order <- object$order; p <- order[1]; d <- order[2]; q <- order[3]
  y_diff <- y_orig; if (d > 0) y_diff <- diff(y_orig, differences = d); n <- length(y_diff)
  c_mu <- coefs["intercept"]; phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0); theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  sigma_hat <- coefs["sigma"]; gamma_hat <- coefs["gamma"]; nu_hat <- if (dist == "t") coefs["nu"] else NULL
  p_skew <- 1/(gamma_hat^2+1)

  mu_t <- numeric(n); eps <- numeric(n); mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0; if (p>0) for (i in 1:p) ar_term <- ar_term+phi[i]*(if (t-i>0) y_diff[t-i] else mean_y)
    ma_term <- 0; if (q>0) for (j in 1:q) ma_term <- ma_term+theta[j]*(if (t-j>0) eps[t-j] else 0)
    mu_t[t] <- c_mu+ar_term+ma_term; eps[t] <- y_diff[t]-mu_t[t]
  }

  z <- eps/sigma_hat; d_obs <- z*(p_skew-(z<0)); d_obs_sorted <- sort(d_obs)

  envelope_mat <- matrix(NA, nrow=n, ncol=B)
  for (b in 1:B) {
    sim_eps <- if (dist == "normal") .rsn(n, 0, sigma_hat, gamma_hat)
               else if (dist == "t") .rskt(n, 0, sigma_hat, gamma_hat, nu_hat)
               else .rskl(n, 0, sigma_hat, gamma_hat)
    z_sim <- sim_eps/sigma_hat; d_sim <- z_sim*(p_skew-(z_sim<0)); envelope_mat[,b] <- sort(d_sim)
  }

  env_lower <- apply(envelope_mat, 1, min); env_upper <- apply(envelope_mat, 1, max); env_median <- apply(envelope_mat, 1, stats::median)
  probs <- (1:n-0.5)/n
  if (dist == "normal") theo_q <- stats::qnorm((1+probs)/2)/2
  else if (dist == "t") theo_q <- stats::qt((1+probs)/2, df=nu_hat)/2
  else theo_q <- -log(1-probs)/4 # Exponential(2) distance for SKL scaled by 1/2. CDF is 1-exp(-2*(2x)). So x = -ln(1-u)/4

  df_env <- data.frame(Theoretical=theo_q, Observed=d_obs_sorted, Lower=env_lower, Upper=env_upper, Median=env_median)
  dist_label <- switch(dist, "normal"="Skew-Normal", "t"="Skewed Student-t", "laplace"="Skewed Laplace")
  pl <- ggplot(df_env, aes(x=Theoretical)) +
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="lightblue", alpha=0.5) +
    geom_line(aes(y=Median), color="blue", linetype="dashed", linewidth=0.6) +
    geom_point(aes(y=Observed), color="black", size=1.5, alpha=0.7) +
    geom_abline(intercept=0, slope=1, color="red", linewidth=0.8) +
    labs(title=paste0("Envelope (", dist_label, ")"), x="Theoretical Quantiles", y="Observed Distances") + theme_minimal(base_size=10)

  print(pl)
  outside <- sum(d_obs_sorted < env_lower | d_obs_sorted > env_upper)
  cat(sprintf("\n=== Envelope Diagnostics [%s] ===\n", dist_label))
  cat(sprintf("Points outside envelope: %d / %d (%.1f%%)\n", outside, n, 100*outside/n))
  invisible(pl)
}

