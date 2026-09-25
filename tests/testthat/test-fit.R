test_that("non-seasonal fits keep the 0.1.0 interface", {
  set.seed(123)
  y <- stats::arima.sim(n = 50, list(ar = 0.5)) + 10
  fit <- fit_modal_arima(as.numeric(y), order = c(1, 0, 0))
  expect_s3_class(fit, "modal_arima")
  expect_equal(fit$convergence, 0)
  expect_named(coef(fit), c("intercept", "ar1", "sigma", "gamma"))
  expect_equal(fit$seasonal$order, c(0L, 0L, 0L))
})

test_that("auto.modal.arima works for non-seasonal series", {
  set.seed(123)
  y <- stats::arima.sim(n = 50, list(ar = 0.5)) + 10
  fit <- auto.modal.arima(as.numeric(y), max.p = 1, max.q = 0)
  expect_s3_class(fit, "modal_arima")
  expect_true(fit$convergence == 0)
})

test_that("the airline model fits for every distribution", {
  y <- log(AirPassengers)
  ref <- stats::arima(y, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  for (dist in c("normal", "t", "laplace")) {
    fit <- fit_modal_arima(y, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12), dist = dist)
    expect_equal(fit$convergence, 0)
    expect_equal(fit$nobs, length(y) - 13)
    expect_true(all(c("ma1", "sma1") %in% names(coef(fit))))
    expect_equal(unname(coef(fit)[c("ma1", "sma1")]), unname(coef(ref)), tolerance = 0.25)
    expect_true(all(is.na(fit$residuals[1:13])))
    expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  }
})

test_that("seasonal specification accepts a vector and uses frequency(y)", {
  y <- log(AirPassengers)
  a <- fit_modal_arima(y, order = c(0, 1, 1), seasonal = c(0, 1, 1))
  b <- fit_modal_arima(y, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  expect_equal(a$loglik, b$loglik)
  expect_error(fit_modal_arima(as.numeric(y), order = c(0, 1, 1), seasonal = c(0, 1, 1)), "period")
})

test_that("seasonal parameters are recovered from simulated data", {
  set.seed(2026)
  est <- replicate(4, {
    y <- sim_modal_sarima(800, c = 0.1, phi = 0.5, theta = 0.3, Phi = 0.6, Theta = -0.3, s = 4,
                          d = 1, D = 1, sigma = 1, gamma = 2, dist = "laplace")
    coef(fit_modal_arima(y, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 1), period = 4),
                         dist = "laplace"))
  })
  truth <- c(ar1 = 0.5, ma1 = 0.3, sar1 = 0.6, sma1 = -0.3, sigma = 1, gamma = 2)
  expect_equal(rowMeans(est)[names(truth)], truth, tolerance = 0.15)
})

test_that("the Skew-Normal log-likelihood equals the 0.1.0 Fernandez-Steel one", {
  # Same data, same maximum: the reparameterization only rescales sigma.
  y <- log10(lynx)
  fit <- fit_modal_arima(y, order = c(2, 0, 0), dist = "normal")
  g <- coef(fit)[["gamma"]]; s_fs <- coef(fit)[["sigma"]] * (g + 1 / g) / 2
  e <- stats::na.omit(as.numeric(residuals(fit)))
  ll_fs <- sum(log(2) - log(s_fs) - log(g + 1 / g) +
                 stats::dnorm(ifelse(e >= 0, e / (s_fs * g), e * g / s_fs), log = TRUE))
  expect_equal(fit$loglik, ll_fs, tolerance = 1e-8)
})

test_that("standard errors exist when nu reaches its upper bound", {
  fit <- fit_modal_arima(USAccDeaths, order = c(0, 1, 1), seasonal = c(0, 1, 1), dist = "t")
  expect_gt(coef(fit)[["nu"]], 150)
  se <- sqrt(diag(vcov(fit)))
  expect_true(all(is.finite(se[names(se) != "nu"])))
})
