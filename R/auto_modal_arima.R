#' Automatic selection of Parametric Modal ARIMA model
#'
#' @details
#' This function performs a grid search over AR and MA orders to find the
#' optimal Modal ARIMA model based on the selected information criterion (AIC or BIC).
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
#' @param max.p Maximum AR order
#' @param max.q Maximum MA order
#' @param ic Information criterion to be used in model selection ("aic", "bic")
#' @param dist Character string specifying the error distribution.
#'   \code{"normal"} (default) for Skew-Normal, \code{"t"} for Skewed Student-t,
#'   \code{"laplace"} for Skewed Laplace.
#'
#' @return An object of class \code{modal_arima}.
#' @importFrom stats AIC BIC
#' @importFrom forecast ndiffs
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
auto.modal.arima <- function(y, d = NA, max.p = 5, max.q = 5,
                              ic = c("aic", "bic"),
                              dist = c("normal", "t", "laplace")) {
  ic <- match.arg(ic)
  dist <- match.arg(dist)

  if (is.na(d)) {
    if (!requireNamespace("forecast", quietly = TRUE)) {
      stop("Package 'forecast' is needed to automatically determine 'd'. Please install it or specify 'd'.")
    }
    d <- forecast::ndiffs(y)
  }

  best_ic <- Inf
  best_mod <- NULL

  for (p in 0:max.p) {
    for (q in 0:max.q) {
      if (p == 0 && q == 0) next
      mod <- tryCatch({
        fit_modal_arima(y, order = c(p, d, q), dist = dist)
      }, error = function(e) NULL)

      if (!is.null(mod) && mod$convergence == 0) {
        current_ic <- if (ic == "aic") AIC(mod) else BIC(mod)
        if (current_ic < best_ic) {
          best_ic <- current_ic
          best_mod <- mod
        }
      }
    }
  }

  if (is.null(best_mod)) {
    best_mod <- fit_modal_arima(y, order = c(0, d, 0), dist = dist)
  }

  return(best_mod)
}
