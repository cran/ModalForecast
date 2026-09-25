#' @export
print.modal_arima <- function(x, ...) {
  dist <- if (!is.null(x$dist)) x$dist else "normal"
  spec <- .object_spec(x)
  cat("\n", .model_label(x), " [", .dist_label(dist), "]\n", sep = "")
  cat("\nCall:\nfit_modal_arima(order = c(", paste(x$order, collapse = ", "), ")", sep = "")
  if (.is_seasonal(spec))
    cat(", seasonal = list(order = c(", spec$P, ", ", spec$D, ", ", spec$Q, "), period = ", spec$s, ")", sep = "")
  cat(", dist = \"", dist, "\")\n", sep = "")
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
vcov.modal_arima <- function(object, ...) {
  if (!is.null(object$vcov)) return(object$vcov)
  # Objects created with version 0.1.0 store only the Hessian on the optimization scale.
  v <- tryCatch(solve(object$hessian), error = function(e) NULL)
  if (is.null(v)) return(matrix(NA, length(object$coefficients), length(object$coefficients)))
  cf <- object$coefficients
  jac <- ifelse(names(cf) %in% c("sigma", "gamma"), cf, ifelse(names(cf) == "nu", cf - 2, 1))
  v <- v * outer(jac, jac)
  dimnames(v) <- list(names(cf), names(cf))
  v
}

#' @export
fitted.modal_arima <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.modal_arima <- function(object, ...) {
  object$residuals
}

.nobs_modal <- function(object) {
  if (!is.null(object$nobs)) object$nobs else length(object$y)
}

#' @export
nobs.modal_arima <- function(object, ...) {
  .nobs_modal(object)
}

#' @export
logLik.modal_arima <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients), nobs = .nobs_modal(object), class = "logLik")
}

#' @export
AIC.modal_arima <- function(object, ..., k = 2) {
  -2 * object$loglik + k * length(object$coefficients)
}

#' @export
BIC.modal_arima <- function(object, ...) {
  -2 * object$loglik + log(.nobs_modal(object)) * length(object$coefficients)
}

#' @export
summary.modal_arima <- function(object, ...) {
  dist <- if (!is.null(object$dist)) object$dist else "normal"
  cf <- object$coefficients
  se <- suppressWarnings(sqrt(diag(vcov(object))))

  sys_names <- grep("^(intercept|ar[0-9]+|ma[0-9]+|sar[0-9]+|sma[0-9]+)$", names(cf), value = TRUE)
  z <- cf[sys_names] / se[sys_names]
  sys_mat <- cbind(Estimate = cf[sys_names], `Std. Error` = se[sys_names],
                   `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))

  skd_names <- grep("^(sigma|gamma|nu)$", names(cf), value = TRUE)
  skd_mat <- cbind(Estimate = cf[skd_names], `Std. Error` = se[skd_names])

  cat("\n", .model_label(object), " Model [", .dist_label(dist), "]\n", sep = "")
  cat("====================================================\n")

  cat("\nSystematic Component:\n")
  printCoefmat(sys_mat)

  cat("\nSKD Family Scale & Skewness:\n")
  print(round(skd_mat, 5))
  z_gamma <- (cf[["gamma"]] - 1) / se[["gamma"]]
  cat(sprintf("Symmetry test H0: gamma = 1:  z = %.3f, p-value = %.4f\n",
              z_gamma, 2 * stats::pnorm(-abs(z_gamma))))

  cat("---\n")
  cat(sprintf("Log-likelihood: %.2f   AIC: %.2f   BIC: %.2f   n = %d\n",
              object$loglik, AIC(object), BIC(object), .nobs_modal(object)))

  status <- if (object$convergence == 0) "Success" else paste("Code", object$convergence)
  cat(sprintf("Convergence: %s\n", status))
  invisible(object)
}

#' @export
plot.modal_arima <- function(x, ...) {
  Time <- Value <- Type <- NULL # Hack for ggplot2 notes
  dist <- if (!is.null(x$dist)) x$dist else "normal"

  time_vec <- if (stats::is.ts(x$y)) as.numeric(stats::time(x$y)) else seq_along(x$y)

  df_orig <- data.frame(Time = time_vec, Value = as.numeric(x$y), Type = "Observed")
  df_fit <- data.frame(Time = time_vec, Value = as.numeric(x$fitted.values), Type = "Fitted Mode")
  df <- rbind(df_orig, df_fit)
  df <- df[!is.na(df$Value), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Value, color = Type, linetype = Type)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Observed" = "gray50", "Fitted Mode" = "blue3")) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid", "Fitted Mode" = "solid")) +
    ggplot2::labs(title = paste0(.model_label(x), " [", .dist_label(dist), "] Fit"),
                  x = "Time", y = "Value", color = "", linetype = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  print(p)
  invisible(p)
}
