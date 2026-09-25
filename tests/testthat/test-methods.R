fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
                           seasonal = list(order = c(0, 1, 1), period = 12), dist = "laplace")

test_that("print and summary label the seasonal model", {
  expect_output(print(fit_air), "Modal SARIMA\\(0,1,1\\)\\(0,1,1\\)\\[12\\]")
  expect_output(summary(fit_air), "sma1")
  expect_output(summary(fit_air), "H0: gamma = 1")
})

test_that("information criteria use the effective sample size", {
  expect_equal(nobs(fit_air), length(AirPassengers) - 13)
  expect_equal(BIC(fit_air), -2 * fit_air$loglik + log(nobs(fit_air)) * length(coef(fit_air)))
  expect_equal(attr(logLik(fit_air), "nobs"), nobs(fit_air))
})

test_that("standard errors of sigma and gamma use the delta method", {
  v <- vcov(fit_air)
  expect_equal(v["gamma", "gamma"], solve(fit_air$hessian)[5, 5] * coef(fit_air)[["gamma"]]^2)
})

test_that("diagnostics, envelope and plot run on seasonal fits", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (dist in c("normal", "t", "laplace")) {
    fit <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1), seasonal = c(0, 1, 1), dist = dist)
    expect_output(res <- diagnostics(fit), "Ljung-Box Test \\(lag=24\\)")
    expect_length(res, 6)
    expect_output(pl <- envelope(fit, B = 20), "Points outside envelope")
    expect_s3_class(pl, "ggplot")
  }
  expect_s3_class(plot(fit_air), "ggplot")
})

test_that("envelope theoretical quantiles match the simulated distances", {
  set.seed(3)
  for (dist in c("normal", "t", "laplace")) {
    nu <- if (dist == "t") 5 else NULL
    e <- .rskd(20000, 0, 1, 2, dist, nu)
    d <- .rho(e, .p_skew(2))
    u <- c(0.25, 0.5, 0.9)
    theo <- switch(dist, normal = stats::qnorm((1 + u) / 2) / 2,
                   t = stats::qt((1 + u) / 2, df = nu) / 2, laplace = -log(1 - u) / 2)
    expect_equal(unname(stats::quantile(d, u)), theo, tolerance = 0.03)
  }
})

test_that("auto.modal.arima selects a seasonal model", {
  fit <- auto.modal.arima(log(AirPassengers), max.p = 1, max.q = 1, max.P = 1, max.Q = 1)
  expect_s3_class(fit, "modal_arima")
  expect_equal(fit$seasonal$period, 12L)
  expect_equal(fit$seasonal$order[2], 1L)
  expect_true(sum(fit$seasonal$order[c(1, 3)]) >= 1)
})

test_that("objects saved with version 0.1.0 still work", {
  old <- fit_modal_arima(log10(lynx), order = c(2, 0, 0))
  old$seasonal <- NULL; old$vcov <- NULL; old$nobs <- NULL
  expect_output(summary(old), "Modal ARIMA\\(2,0,0\\)")
  expect_length(predict(old, n.ahead = 3), 3)
  expect_s3_class(forecast(old, h = 3), "forecast")
})

test_that("envelope works with and without refitting", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1), seasonal = c(0, 1, 1))
  expect_output(envelope(fit, B = 10), "Expected under a correct model: about 18.2%")
  out <- capture.output(envelope(fit, B = 10, refit = FALSE))
  expect_false(any(grepl("Expected under", out)))
})

test_that("simulated differenced series follow the fitted model", {
  set.seed(4)
  pp <- list(c = 0.2, phi = 0.5, theta = numeric(0), Phi = numeric(0), Theta = numeric(0),
             sigma = 1, gamma = 1, nu = NULL)
  spec <- list(p = 1, d = 0, q = 0, P = 0, D = 0, Q = 0, s = 1, dist = "normal")
  w <- .simulate_differenced(pp, spec, 5000)
  # sd of the mean of this AR(1) is about 0.028
  expect_lt(abs(mean(w) - 0.2 / 0.5), 0.12)
  expect_equal(stats::acf(w, plot = FALSE)$acf[2], 0.5, tolerance = 0.05)
})
