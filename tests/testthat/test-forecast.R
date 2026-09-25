fit_air <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1),
                           seasonal = list(order = c(0, 1, 1), period = 12))

test_that("seasonal forecasts have the right shape and time index", {
  fc <- forecast(fit_air, h = 24, level = c(80, 95))
  expect_s3_class(fc, "forecast")
  expect_length(fc$mean, 24)
  expect_equal(stats::start(fc$mean), c(1961, 1))
  expect_equal(stats::frequency(fc$mean), 12)
  expect_true(all(fc$lower[, "95%"] < fc$lower[, "80%"]))
  expect_true(all(fc$lower[, "80%"] < fc$mean & fc$mean < fc$upper[, "80%"]))
  expect_true(all(fc$upper[, "80%"] < fc$upper[, "95%"]))
})

test_that("predict undoes regular and seasonal differencing", {
  st <- .modal_filter(fit_air)
  y <- as.numeric(fit_air$y); h <- 15
  w_fc <- .simulate_w(st, matrix(0, h, 1))[, 1]
  y_ext <- c(y, rep(NA, h)); n <- length(y)
  for (k in seq_len(h)) {
    t <- n + k
    y_ext[t] <- w_fc[k] + y_ext[t - 1] + y_ext[t - 12] - y_ext[t - 13]
  }
  expect_equal(predict(fit_air, n.ahead = h), y_ext[n + seq_len(h)])
})

test_that("joint and marginal modes agree at h = 1 and differ later", {
  fc_j <- forecast(fit_air, h = 12)
  fc_m <- forecast(fit_air, h = 12, point = "marginal")
  expect_equal(as.numeric(fc_j$mean), as.numeric(fc_j$mode_path))
  expect_equal(fc_m$mean[1], fc_j$mean[1])
  expect_equal(as.numeric(fc_m$mean - fc_j$mean), as.numeric(fc_j$mode_shift))
})

test_that("the marginal mode correction vanishes for symmetric errors", {
  st <- .modal_filter(fit_air)
  st$pp$gamma <- 1
  psi <- .psi_weights(st, 10)
  expect_lt(abs(.error_distribution(psi, st$pp, "normal")$mode), 0.01 * st$pp$sigma)
})

test_that("the marginal mode matches Monte Carlo simulation", {
  set.seed(11)
  pp <- list(sigma = 1, gamma = 2, nu = NULL)
  psi <- c(1, 0.8, 0.8, 0.5)
  ed <- .error_distribution(psi, pp, "laplace")
  E <- matrix(.rskd(4e5 * length(psi), 0, 1, 2, "laplace"), ncol = length(psi))
  S <- as.numeric(E %*% psi)
  h <- graphics::hist(S[S > -3 & S < 8], breaks = seq(-3, 8, by = 0.05), plot = FALSE)
  sm <- stats::filter(h$counts, rep(1 / 7, 7))
  mc_mode <- h$mids[which.max(sm)]
  expect_equal(ed$mode, mc_mode, tolerance = 0.15)
  expect_equal(ed$quantile(c(0.1, 0.5, 0.9)), unname(stats::quantile(S, c(0.1, 0.5, 0.9))), tolerance = 0.03)
})

test_that("bootstrap and exact prediction intervals agree", {
  set.seed(5)
  fc_a <- forecast(fit_air, h = 12, level = 90)
  fc_b <- forecast(fit_air, h = 12, level = 90, interval = "bootstrap", npaths = 4000)
  expect_equal(as.numeric(fc_b$lower), as.numeric(fc_a$lower), tolerance = 0.01)
  expect_equal(as.numeric(fc_b$upper), as.numeric(fc_a$upper), tolerance = 0.01)
})

test_that("forecasts work for Student-t and Laplace seasonal fits", {
  for (dist in c("t", "laplace")) {
    fit <- fit_modal_arima(log(AirPassengers), order = c(0, 1, 1), seasonal = c(0, 1, 1), dist = dist)
    fc <- forecast(fit, h = 6, point = "marginal")
    expect_true(all(is.finite(fc$mean)))
    expect_true(all(fc$lower[, 2] < fc$upper[, 2]))
  }
})
