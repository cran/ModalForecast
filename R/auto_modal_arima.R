#' Automatic selection of Parametric Modal ARIMA and SARIMA models
#'
#' @details
#' This function performs a grid search over AR and MA orders, and over
#' seasonal AR and MA orders when the series is seasonal, to find the Modal
#' (S)ARIMA model with the smallest information criterion (AIC or BIC). As in
#' \code{forecast::auto.arima}, the seasonal differencing order \code{D} is
#' chosen with \code{forecast::nsdiffs} and then the regular order \code{d} with
#' \code{forecast::ndiffs} on the seasonally differenced series, unless they are
#' supplied.
#'
#' @references
#' Galarza, C. E., Lachos, V. H., Cabral, C. R. B., and Castro, L. M. (2017).
#' Robust quantile regression using a generalized class of skewed distributions.
#' Stat, 6(1), 113-130.
#'
#' @seealso \code{\link{fit_modal_arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations
#' @param d Integer, degree of differencing. If NA, it's determined automatically.
#' @param D Integer, degree of seasonal differencing. If NA, it's determined
#'   automatically. Ignored for non-seasonal models.
#' @param max.p Maximum AR order
#' @param max.q Maximum MA order
#' @param max.P Maximum seasonal AR order
#' @param max.Q Maximum seasonal MA order
#' @param seasonal Logical. If \code{TRUE} (default) and the seasonal period is
#'   greater than one, seasonal models are considered.
#' @param period Seasonal period. Defaults to \code{frequency(y)}.
#' @param ic Information criterion to be used in model selection ("aic", "bic")
#' @param dist Character string specifying the error distribution.
#'   \code{"normal"} (default) for Skew-Normal, \code{"t"} for Skewed Student-t,
#'   \code{"laplace"} for Skewed Laplace.
#' @param trace Logical. If \code{TRUE}, prints each model and its criterion.
#'
#' @return An object of class \code{modal_arima}.
#' @importFrom stats AIC BIC
#' @importFrom forecast ndiffs nsdiffs
#' @export
#'
#' @examples
#' library(forecast)
#'
#' # Non-seasonal
#' fit_auto <- auto.modal.arima(log10(lynx), d = 0, max.p = 2, max.q = 2)
#' summary(fit_auto)
#'
#' # Seasonal, with a small search grid
#' fit_sauto <- auto.modal.arima(log(AirPassengers), max.p = 1, max.q = 1,
#'                               max.P = 1, max.Q = 1)
#' fit_sauto
auto.modal.arima <- function(y, d = NA, D = NA, max.p = 5, max.q = 5,
                             max.P = 1, max.Q = 1, seasonal = TRUE,
                             period = stats::frequency(y),
                             ic = c("aic", "bic"),
                             dist = c("normal", "t", "laplace"),
                             trace = FALSE) {
  ic <- match.arg(ic)
  dist <- match.arg(dist)
  seasonal <- seasonal && period > 1
  if (!seasonal) { D <- 0; max.P <- 0; max.Q <- 0 }

  if (is.na(D)) {
    D <- if (length(y) >= 2 * period + 1)
      forecast::nsdiffs(stats::ts(as.numeric(y), frequency = period)) else 0
  }
  if (is.na(d)) {
    ys <- if (D > 0) diff(as.numeric(y), lag = period, differences = D) else as.numeric(y)
    d <- forecast::ndiffs(ys)
  }

  best_ic <- Inf
  best_mod <- NULL
  for (P in 0:max.P) for (Q in 0:max.Q) for (p in 0:max.p) for (q in 0:max.q) {
    if (p + q + P + Q == 0) next
    mod <- tryCatch(
      fit_modal_arima(y, order = c(p, d, q),
                      seasonal = list(order = c(P, D, Q), period = period), dist = dist),
      error = function(e) NULL)
    if (!is.null(mod) && mod$convergence == 0) {
      current_ic <- if (ic == "aic") AIC(mod) else BIC(mod)
      if (trace) cat(sprintf("%-40s %s = %.3f\n", .model_label(mod), toupper(ic), current_ic))
      if (current_ic < best_ic) {
        best_ic <- current_ic
        best_mod <- mod
      }
    }
  }

  if (is.null(best_mod)) {
    best_mod <- fit_modal_arima(y, order = c(0, d, 0),
                                seasonal = list(order = c(0, D, 0), period = period), dist = dist)
  }
  best_mod
}
