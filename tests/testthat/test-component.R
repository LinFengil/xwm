tol <- 1e-10

test_that("RemoveSeasonal returns seasonal means and anomalies", {
  y <- c(1, 2, 3, 4)
  t <- c(1, 1, 2, 2)
  out <- RemoveSeasonal(y, t)
  expect_equal(out$mu, c(1.5, 3.5))
  expect_equal(out$anomaly, c(-0.5, 0.5, -0.5, 0.5))
})

test_that("lm_solve fits simple linear regression", {
  X <- 1:4
  Y <- 1 + 2 * X
  out <- lm_solve(cbind(1, X), Y)
  expect_equal(as.numeric(out$coef), c(1, 2), tolerance = tol)
  expect_equal(out$fitted, Y, tolerance = tol)
  expect_equal(out$SSE, 0, tolerance = tol)
})

test_that("LS_OneBreak with no break equals linear fit", {
  X <- 1:4
  Y <- 1 + 2 * X
  out <- LS_OneBreak(Y, X, Ic = NA)
  expect_equal(as.numeric(out$coef), c(1, 2), tolerance = tol)
  expect_equal(out$SSE, 0, tolerance = tol)
})

test_that("LS_OneBreak recovers step change", {
  X <- 1:4
  Y <- c(0, 0, 3, 3)
  out <- LS_OneBreak(Y, X, Ic = 2)
  expect_equal(as.numeric(out$coef), c(0, 0, 3), tolerance = tol)
  expect_equal(out$fitted, Y, tolerance = tol)
  expect_equal(out$SSE, 0, tolerance = tol)
})

test_that("LS_MultiBreaks recovers step change", {
  X <- 1:4
  Y <- c(1, 1, 2, 2)
  out <- LS_MultiBreaks(Y, X, Ips = c(2, 4))
  expect_equal(as.numeric(out$coef), c(1, 0, 1), tolerance = tol)
  expect_equal(out$fitted, Y, tolerance = tol)
  expect_equal(out$SSE, 0, tolerance = tol)
})
