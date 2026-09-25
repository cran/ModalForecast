dists <- c("normal", "t", "laplace")

test_that("SKD densities integrate to one and have P(Y < mu) = p", {
  for (dist in dists) for (g in c(0.5, 1, 2.5)) {
    nu <- if (dist == "t") 4 else NULL
    f <- function(x) .dskd(x, 0.3, 1.2, g, dist, nu)
    left <- stats::integrate(f, -Inf, 0.3)$value
    right <- stats::integrate(f, 0.3, Inf)$value
    expect_equal(left + right, 1, tolerance = 1e-6)
    expect_equal(left, 1 / (1 + g^2), tolerance = 1e-6)
    expect_equal(.pskd(0.3, 0.3, 1.2, g, dist, nu), 1 / (1 + g^2))
  }
})

test_that("the mode of the SKD density is mu", {
  for (dist in dists) {
    nu <- if (dist == "t") 4 else NULL
    opt <- stats::optimize(function(x) .dskd(x, 0.7, 1, 2, dist, nu), c(-3, 3), maximum = TRUE, tol = 1e-10)
    expect_equal(opt$maximum, 0.7, tolerance = 1e-6)
  }
})

test_that("the Skew-Normal is the SKN of Galarza et al. (2017)", {
  x <- seq(-4, 4, by = 0.25); g <- 2; p <- 1 / (1 + g^2); s <- 1.3
  galarza <- 4 * p * (1 - p) / sqrt(2 * pi * s^2) * exp(-2 * (x / s * (p - (x < 0)))^2)
  expect_equal(.dskd(x, 0, s, g, "normal"), galarza)
})

test_that("CDF, quantile and density are consistent", {
  for (dist in dists) for (g in c(0.6, 2)) {
    nu <- if (dist == "t") 5 else NULL
    u <- c(0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999)
    expect_equal(.pskd(.qskd(u, 1, 2, g, dist, nu), 1, 2, g, dist, nu), u, tolerance = 1e-10)
    q <- c(-2, -0.5, 0.4, 3)
    num <- vapply(q, function(v) stats::integrate(function(x) .dskd(x, 0, 1, g, dist, nu), -Inf, v)$value, 1)
    expect_equal(.pskd(q, 0, 1, g, dist, nu), num, tolerance = 1e-6)
  }
})

test_that("random generation follows the SKD distribution", {
  set.seed(42)
  for (dist in dists) {
    nu <- if (dist == "t") 4 else NULL
    x <- .rskd(4000, 0, 1.5, 2, dist, nu)
    expect_gt(stats::ks.test(x, function(q) .pskd(q, 0, 1.5, 2, dist, nu))$p.value, 0.001)
  }
})
