#' @export
print.modal_arima <- function(x, ...) {
  dist <- if (!is.null(x$dist)) x$dist else "normal"
  cat("\nCall:\nfit_modal_arima(order = c(", paste(x$order, collapse=", "), "), dist = \"", dist, "\")\n", sep="")
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nLog-likelihood:", round(x$loglik, 2), "\n")
  invisible(x)
}

#' @export
coef.modal_arima <- function(object, ...) {
  object$coefficients
}

#' @export
fitted.modal_arima <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.modal_arima <- function(object, ...) {
  object$residuals
}

#' @export
logLik.modal_arima <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients), nobs = length(object$y), class = "logLik")
}

#' @export
AIC.modal_arima <- function(object, ..., k = 2) {
  -2 * object$loglik + k * length(object$coefficients)
}

#' @export
BIC.modal_arima <- function(object, ...) {
  -2 * object$loglik + log(length(object$y)) * length(object$coefficients)
}

#' @export
summary.modal_arima <- function(object, ...) {
  dist <- if (!is.null(object$dist)) object$dist else "normal"

  se <- suppressWarnings(tryCatch(sqrt(diag(solve(object$hessian))), error = function(e) rep(NA, length(object$coefficients))))
  z <- object$coefficients / se
  pval <- 2 * (1 - pnorm(abs(z)))

  res <- cbind(Estimate = object$coefficients, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = pval)
  
  # Separate systematic parameters from SKD (distribution) parameters
  sys_names <- grep("^(intercept|ar[0-9]+|ma[0-9]+)$", rownames(res), value = TRUE)
  skd_names <- grep("^(sigma|gamma|nu)$", rownames(res), value = TRUE)
  
  sys_mat <- res[sys_names, , drop = FALSE]
  skd_mat <- res[skd_names, , drop = FALSE]

  dist_label <- switch(dist, "normal"="Skew-Normal", "t"="Skewed Student-t", "laplace"="Skewed Laplace")
  cat("\nModal ARIMA(", paste(object$order, collapse=","), ") Model [", dist_label, "]\n", sep="")
  cat("====================================================\n")
  
  cat("\nSystematic Component:\n")
  printCoefmat(sys_mat)
  
  cat("\nSKD Family Scale & Skewness:\n")
  printCoefmat(skd_mat)
  
  cat("---\n")
  cat(sprintf("Log-likelihood: %.2f   AIC: %.2f   BIC: %.2f\n", object$loglik, AIC(object), BIC(object)))
  
  status <- if(object$convergence == 0) "Success" else paste("Code", object$convergence)
  cat(sprintf("Convergence: %s\n", status))
}

#' @export
plot.modal_arima <- function(x, ...) {
  Time <- Value <- Type <- NULL # Hack for ggplot2 notes
  dist <- if (!is.null(x$dist)) x$dist else "SKN"
  
  time_vec <- if(stats::is.ts(x$y)) as.numeric(stats::time(x$y)) else 1:length(x$y)
  
  df_orig <- data.frame(Time = time_vec, Value = as.numeric(x$y), Type = "Observed")
  df_fit <- data.frame(Time = time_vec, Value = as.numeric(x$fitted.values), Type = "Fitted Mode")
  df <- rbind(df_orig, df_fit)
  
  # Remove NAs from fitted for plotting
  df <- df[!is.na(df$Value), ]
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Value, color = Type, linetype = Type)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Observed" = "gray50", "Fitted Mode" = "blue3")) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid", "Fitted Mode" = "solid")) +
    ggplot2::labs(title = paste0("Modal ARIMA(", paste(x$order, collapse=","), ") [", dist, "] Fit"),
                  x = "Time", y = "Value", color = "", linetype = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
    
  print(p)
  invisible(p)
}
