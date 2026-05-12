test_that("auto.modal_arima works", {
  set.seed(123)
  y <- stats::arima.sim(n = 50, list(ar = 0.5)) + 10
  fit <- auto.modal.arima(as.numeric(y), max.p = 1, max.q = 0)
  expect_s3_class(fit, "modal_arima")
  expect_true(fit$convergence == 0)
})

test_that("predict and forecast work", {
  set.seed(123)
  y <- stats::arima.sim(n = 50, list(ar = 0.5)) + 10
  fit <- fit_modal_arima(as.numeric(y), order = c(1, 0, 0))
  
  pred <- predict(fit, n.ahead = 5)
  expect_length(pred, 5)
  
  fc <- forecast(fit, h = 5)
  expect_s3_class(fc, "forecast")
  expect_length(fc$mean, 5)
})
