test_that("analytical gradients match finite differences for seasonal models", {
  set.seed(7)
  w <- as.numeric(stats::arima.sim(list(ar = 0.5, ma = 0.3), 150)) + 0.3
  cases <- list(
    list(spec = list(p = 1, d = 0, q = 1, P = 1, D = 0, Q = 1, s = 12),
         beta = c(0.2, 0.4, 0.25, 0.3, -0.35)),
    list(spec = list(p = 2, d = 0, q = 2, P = 2, D = 0, Q = 1, s = 4),
         beta = c(0.1, 0.3, -0.2, 0.2, 0.1, 0.3, 0.2, -0.4)),
    list(spec = list(p = 1, d = 0, q = 0, P = 0, D = 0, Q = 0, s = 1),
         beta = c(0.1, 0.6)))
  for (cs in cases) for (dist in c("normal", "t", "laplace")) {
    spec <- c(cs$spec, dist = dist)
    par <- c(cs$beta, log(0.9), log(1.6), if (dist == "t") log(3))
    analytic <- .grad_neg_loglik(par, spec, w)
    numeric <- num_grad(function(x) .neg_loglik(x, spec, w), par)
    expect_equal(analytic, numeric, tolerance = 1e-5)
  }
})

test_that("non-admissible parameters are penalized", {
  spec <- list(p = 1, d = 0, q = 0, P = 1, D = 0, Q = 0, s = 12, dist = "normal")
  w <- rnorm(100)
  expect_equal(.neg_loglik(c(0, 0.5, 1.2, 0, 0), spec, w), 1e10)
})
