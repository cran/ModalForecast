#' Diagnostic Plots for Modal ARIMA and SARIMA Models
#'
#' @description
#' Provides visual and statistical diagnostics for the residuals of a fitted
#' Modal (S)ARIMA model. Produces a panel with the differenced series and
#' fitted modes, the ACF and PACF of the randomized quantile residuals (RQR),
#' a normal QQ-plot, a histogram and Ljung-Box p-values. For seasonal models
#' the ACF, PACF and Ljung-Box statistics cover at least two seasonal periods.
#'
#' @param object An object of class \code{modal_arima}.
#' @param ... Additional arguments (unused).
#' @return A list of ggplot objects (invisibly) and draws the panel.
#' @import ggplot2
#' @importFrom stats acf pacf shapiro.test Box.test qqnorm qqline predict is.ts ts start frequency dnorm qnorm pnorm time residuals
#' @importFrom scales comma_format
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' # Non-seasonal
#' fit <- fit_modal_arima(log10(lynx), order = c(2, 0, 0))
#' diagnostics(fit)
#'
#' # Seasonal airline model
#' fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
#'                            seasonal = list(order = c(0, 1, 1), period = 12))
#' diagnostics(fit_air)
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.modal_arima <- function(object, ...) {
  # Global variable hack for ggplot2
  Time <- Observed <- Fitted <- Lag <- ACF <- PACF <- Residuals <- Pval <- density <- NULL

  st <- .modal_filter(object)
  dist <- st$spec$dist
  y_diff <- st$w; mu_t <- st$mu; n <- length(y_diff)
  sigma_hat <- st$pp$sigma; gamma_hat <- st$pp$gamma; nu_hat <- st$pp$nu

  # Randomized Quantile Residuals via the SKD CDF
  u <- .pskd(y_diff, mu_t, sigma_hat, gamma_hat, dist, nu_hat)
  u <- pmax(1e-7, pmin(1-1e-7, u)); rqr <- stats::qnorm(u)

  time_idx <- if (stats::is.ts(object$y)) utils::tail(as.numeric(stats::time(object$y)), n) else seq_len(n)
  df_fits <- data.frame(Time=time_idx, Observed=as.numeric(y_diff), Fitted=as.numeric(mu_t))

  p1 <- ggplot(df_fits, aes(x=Time)) + geom_line(aes(y=Observed, color="Observed"), alpha=0.6) +
    geom_line(aes(y=Fitted, color="Fitted Mode"), linewidth=0.75) +
    scale_color_manual(values=c("Observed"="gray40", "Fitted Mode"="blue")) +
    labs(title="Series & Fitted Modes", y="Value", x="Time", color=NULL) + theme_minimal(base_size=7) + theme(legend.position="top")

  max_lag <- min(max(20, 2 * st$spec$s + 1), n - 1)
  bacf <- stats::acf(rqr, plot=FALSE, lag.max=max_lag)
  df_acf <- data.frame(Lag=as.numeric(bacf$lag[-1]), ACF=as.numeric(bacf$acf[-1]))
  p2 <- ggplot(df_acf, aes(x=Lag, y=ACF)) + geom_bar(stat="identity", width=0.1, fill="blue") +
    geom_hline(yintercept=c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype="dashed", color="red") + labs(title="ACF of Residuals (RQR)") + theme_minimal(base_size=7)

  bpacf <- stats::pacf(rqr, plot=FALSE, lag.max=max_lag)
  df_pacf <- data.frame(Lag=as.numeric(bpacf$lag), PACF=as.numeric(bpacf$acf))
  p3 <- ggplot(df_pacf, aes(x=Lag, y=PACF)) + geom_bar(stat="identity", width=0.1, fill="blue") +
    geom_hline(yintercept=c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype="dashed", color="red") + labs(title="PACF of Residuals (RQR)") + theme_minimal(base_size=7)

  df_rqr <- data.frame(Residuals=rqr)
  p4 <- ggplot(df_rqr, aes(sample=Residuals)) + stat_qq(color="steelblue", alpha=0.6) + stat_qq_line(color="red", linewidth=0.8) +
    labs(title="Normal QQ-Plot (RQR)") + theme_minimal(base_size=7)

  p5 <- ggplot(df_rqr, aes(x=Residuals)) + geom_histogram(aes(y=after_stat(density)), bins=20, fill="steelblue", alpha=0.4, color="white") +
    stat_function(fun=stats::dnorm, color="red", linewidth=0.8) + labs(title="Histogram (RQR)") + theme_minimal(base_size=7)

  lb_pvals <- sapply(1:max_lag, function(lag) stats::Box.test(rqr, lag=lag, type="Ljung-Box")$p.value)
  df_lb <- data.frame(Lag=1:max_lag, Pval=lb_pvals)
  p6 <- ggplot(df_lb, aes(x=Lag, y=Pval)) + geom_point(color="blue", size=2) +
    geom_hline(yintercept=0.05, linetype="dashed", color="red", linewidth=0.8) + ylim(0,1) + labs(title="Ljung-Box p-values", y="p-value") + theme_minimal(base_size=7)

  dist_label <- paste0("[", .dist_label(dist), "]")
  title_grob <- grid::textGrob(paste(.model_label(object), dist_label, "Diagnostic Panel"), gp=grid::gpar(fontsize=10, fontface="bold"))
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, top=title_grob)

  cat("\n=== Diagnostic Tests ===\n")
  cat(sprintf("Distribution: %s\n", dist))
  sw <- stats::shapiro.test(rqr)
  cat(sprintf("Shapiro-Wilk Normality Test: W = %.4f, p-value = %.4f\n", sw$statistic, sw$p.value))
  if (sw$p.value > 0.05) cat("  -> Normality assumption met (Expected/Good fit).\n") else cat("  -> Normality assumption rejected (Requires investigation).\n")
  
  lb_lag <- if (st$spec$s > 1) min(2 * st$spec$s, n - 1) else min(10, n - 1)
  lb <- stats::Box.test(rqr, lag = lb_lag, type = "Ljung-Box")
  cat(sprintf("Ljung-Box Test (lag=%d):     X2 = %.4f, p-value = %.4f\n", lb_lag, lb$statistic, lb$p.value))
  if (lb$p.value > 0.05) cat("  -> Residuals are independent / No Autocorrelation (Expected/Good fit).\n") else cat("  -> Residuals are autocorrelated (Investigate lags).\n")
  cat(sprintf("\nEstimated gamma (skewness):  %.4f\n", gamma_hat))
  gamma_se <- suppressWarnings(sqrt(vcov(object)["gamma", "gamma"]))
  if (is.finite(gamma_se)) {
      p_gamma <- 2 * (1 - stats::pnorm(abs((gamma_hat - 1) / gamma_se)))
      if (p_gamma > 0.05) cat("  -> Symmetry confirmed (gamma = 1) (Expected if true series is symmetric).\n\n") else cat("  -> Symmetry rejected / SKD Skewness present (Expected by model).\n\n")
  }
  if (dist == "t") cat(sprintf("Estimated nu (d.f.):        %.4f\n\n", nu_hat))

  invisible(list(p1, p2, p3, p4, p5, p6))
}

#' Simulation Envelope Diagnostics for Modal ARIMA and SARIMA Models
#'
#' Constructs simulation envelopes for the distances
#' \eqn{\rho_p(\epsilon_t/\sigma)} of the modal residuals, based on their
#' theoretical distribution under the SKD family (Galarza et al., 2017):
#' twice the distance follows the half-normal, half-t or standard exponential
#' distribution for the Skew-Normal, Skewed Student-t and Skewed Laplace
#' members, respectively.
#'
#' @details
#' With \code{refit = TRUE} (the default), each replication simulates a series
#' from the fitted model, re-estimates the model and computes the distances of
#' the re-estimated residuals, so the envelope accounts for parameter
#' estimation. This matters most for the Skewed Laplace, whose maximum
#' likelihood fit, like least absolute deviations, pulls several residuals to
#' zero: without refitting, the smallest distances fall below the envelope even
#' when the model is correct. With \code{refit = FALSE}, the envelope is built
#' from simulated innovations at the estimated parameters, which is faster.
#'
#' @param object An object of class \code{modal_arima}.
#' @param B Number of Monte Carlo replications for envelope construction. Default is 100.
#' @param refit Logical. If \code{TRUE} (default), re-estimate the model on
#'   each simulated series.
#' @param ... Additional arguments (unused).
#' @return A ggplot object (invisibly) and draws the envelope plot.
#' @import ggplot2
#' @export
#'
#' @examples
#' fit <- fit_modal_arima(log10(lynx), order = c(2, 0, 0))
#' envelope(fit, B = 10)
#'
#' # Seasonal model; refit = FALSE is faster, refit = TRUE accounts for estimation
#' fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
#'                            seasonal = list(order = c(0, 1, 1), period = 12))
#' envelope(fit_air, B = 10, refit = FALSE)
envelope <- function(object, ...) { UseMethod("envelope") }

#' @rdname envelope
#' @export
envelope.modal_arima <- function(object, B = 100, refit = TRUE, ...) {
  Theoretical <- Observed <- Lower <- Upper <- Median <- NULL # ggplot2 hack
  st <- .modal_filter(object)
  dist <- st$spec$dist
  sigma_hat <- st$pp$sigma; gamma_hat <- st$pp$gamma; nu_hat <- st$pp$nu
  p_skew <- .p_skew(gamma_hat)
  eps <- st$eps; n <- length(eps)

  # Distances d = rho_p(eps / sigma); 2d is distributed as |X| with X ~ g.
  z <- eps/sigma_hat; d_obs <- .rho(z, p_skew); d_obs_sorted <- sort(d_obs)

  # Refits work on the differenced series, so the model has d = D = 0.
  spec_w <- st$spec; spec_w$d <- 0; spec_w$D <- 0
  par_hat <- .par_from_coef(object$coefficients, st$spec)
  envelope_mat <- matrix(NA, nrow=n, ncol=B)
  for (b in 1:B) {
    if (refit) {
      w_sim <- .simulate_differenced(st$pp, spec_w, n)
      opt <- tryCatch(.mle(w_sim, spec_w, par_hat, reltol = 1e-8), error = function(e) NULL)
      if (is.null(opt)) next
      pp_b <- .split_par(opt$par, spec_w)
      e_b <- .recursion(pp_b, spec_w, w_sim)$eps
      envelope_mat[,b] <- sort(.rho(e_b/pp_b$sigma, .p_skew(pp_b$gamma)))
    } else {
      sim_eps <- .rskd(n, 0, sigma_hat, gamma_hat, dist, nu_hat)
      envelope_mat[,b] <- sort(.rho(sim_eps/sigma_hat, p_skew))
    }
  }
  envelope_mat <- envelope_mat[, colSums(is.na(envelope_mat)) == 0, drop = FALSE]

  env_lower <- apply(envelope_mat, 1, min); env_upper <- apply(envelope_mat, 1, max); env_median <- apply(envelope_mat, 1, stats::median)
  probs <- (1:n-0.5)/n
  theo_q <- if (dist == "normal") stats::qnorm((1+probs)/2)/2
            else if (dist == "t") stats::qt((1+probs)/2, df=nu_hat)/2
            else -log(1-probs)/2  # |X| ~ Exp(1) for the standard Laplace kernel

  df_env <- data.frame(Theoretical=theo_q, Observed=d_obs_sorted, Lower=env_lower, Upper=env_upper, Median=env_median)
  dist_label <- switch(dist, "normal"="Skew-Normal", "t"="Skewed Student-t", "laplace"="Skewed Laplace")
  pl <- ggplot(df_env, aes(x=Theoretical)) +
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="lightblue", alpha=0.5) +
    geom_line(aes(y=Median), color="blue", linetype="dashed", linewidth=0.6) +
    geom_point(aes(y=Observed), color="black", size=1.5, alpha=0.7) +
    geom_abline(intercept=0, slope=1, color="red", linewidth=0.8) +
    labs(title=paste0("Envelope (", dist_label, ")"), subtitle=.model_label(object), x="Theoretical Quantiles", y="Observed Distances") + theme_minimal(base_size=10)

  print(pl)
  outside <- sum(d_obs_sorted < env_lower | d_obs_sorted > env_upper)
  cat(sprintf("\n=== Envelope Diagnostics [%s] ===\n", dist_label))
  cat(sprintf("Points outside envelope: %d / %d (%.1f%%)\n", outside, n, 100*outside/n))
  if (refit)
    cat(sprintf("Expected under a correct model: about %.1f%% (each point falls outside a %d-replication envelope with probability 2/(B+1))\n",
                200 / (ncol(envelope_mat) + 1), ncol(envelope_mat)))
  invisible(pl)
}

