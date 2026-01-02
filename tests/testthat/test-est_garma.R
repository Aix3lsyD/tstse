test_that("est_garma returns correct structure", {
  x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
                 plot = FALSE, seed = 123)

  # Coarse grid for speed
  fit <- est_garma(x, u_range = c(0.6, 0.9, 0.1),
                   lambda_range = c(0.1, 0.4, 0.1), p_max = 3)

  expect_s3_class(fit, "est_garma")
  expect_named(fit, c("u", "lambda", "phi", "p", "var_a", "aic", "table"))
  expect_true(is.numeric(fit$u))
  expect_true(is.numeric(fit$lambda))
  expect_true(is.numeric(fit$phi))
  expect_true(is.numeric(fit$var_a))
  expect_true(is.numeric(fit$aic))
  expect_s3_class(fit$table, "data.frame")
})

test_that("est_garma finds parameters close to true values", {
  # Generate data with known parameters
  x <- gen_garma(n = 300, u = 0.8, lambda = 0.3, phi = 0.5,
                 plot = FALSE, seed = 123)

  # Estimate with grid centered on true values
  fit <- est_garma(x, u_range = c(0.6, 1.0, 0.1),
                   lambda_range = c(0.1, 0.5, 0.1), p_max = 3)

  # Should find parameters reasonably close to true
  expect_true(abs(fit$u - 0.8) <= 0.2)
  expect_true(abs(fit$lambda - 0.3) <= 0.2)
})

test_that("est_garma grid search produces valid results", {
  x <- gen_garma(n = 150, u = 0.5, lambda = 0.2, plot = FALSE, seed = 42)

  fit <- est_garma(x, u_range = c(0.3, 0.7, 0.1),
                   lambda_range = c(0.1, 0.4, 0.1), p_max = 2)

  # Table should have all grid combinations
  u_vals <- seq(0.3, 0.7, by = 0.1)
  lambda_vals <- seq(0.1, 0.4, by = 0.1)
  expected_rows <- length(u_vals) * length(lambda_vals)
  expect_equal(nrow(fit$table), expected_rows)

  # All table values should be finite
  expect_true(all(is.finite(fit$table$u)))
  expect_true(all(is.finite(fit$table$lambda)))
  expect_true(all(is.finite(fit$table$garma_aic)))
})

test_that("est_garma validates inputs", {
  x <- rnorm(100)

  expect_error(est_garma(1:10, u_range = c(0, 1, 0.1), lambda_range = c(0, 0.5, 0.1)),
               "at least 20 observations")
  expect_error(est_garma(x, u_range = c(0.5, 0.3, 0.1), lambda_range = c(0, 0.5, 0.1)),
               "must be less than")
  expect_error(est_garma(x, u_range = c(0, 1, -0.1), lambda_range = c(0, 0.5, 0.1)),
               "must be positive")
  expect_error(est_garma(x, u_range = c(-2, 0, 0.1), lambda_range = c(0, 0.5, 0.1)),
               "must be in")
  expect_error(est_garma(x, u_range = c(0, 1, 0.1), lambda_range = c(-0.1, 0.5, 0.1)),
               "must be non-negative")
  expect_error(est_garma(x, u_range = c(0, 1, 0.1), lambda_range = c(0, 0.5, 0.1),
                         p_max = -1),
               "must be a non-negative integer")
})

test_that("est_garma print method works", {
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.2, plot = FALSE, seed = 123)
  fit <- est_garma(x, u_range = c(0.3, 0.7, 0.2),
                   lambda_range = c(0.1, 0.3, 0.1), p_max = 2)

  output <- capture.output(print(fit))
  expect_true(any(grepl("GARMA", output)))
  expect_true(any(grepl("Gegenbauer frequency", output)))
  expect_true(any(grepl("lambda", output)))
})

test_that("est_garma handles p_max = 0", {
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.2, plot = FALSE, seed = 123)

  # Only fit Gegenbauer, no AR
  fit <- est_garma(x, u_range = c(0.3, 0.7, 0.2),
                   lambda_range = c(0.1, 0.3, 0.1), p_max = 0)

  expect_equal(fit$p, 0)
  expect_equal(fit$phi, 0)
})

test_that("est_garma matches tswge::est.garma.wge", {
  skip_if_not_installed("tswge")

  x <- gen_garma(n = 200, u = 0.7, lambda = 0.25, phi = 0.4,
                 plot = FALSE, seed = 789)

  # Use same grid
  old <- suppressWarnings(tswge::est.garma.wge(x,
    low.u = 0.5, high.u = 0.9, inc.u = 0.1,
    low.lambda = 0.1, high.lambda = 0.4, inc.lambda = 0.1,
    p.max = 3, nback = 500))

  new <- est_garma(x, u_range = c(0.5, 0.9, 0.1),
                   lambda_range = c(0.1, 0.4, 0.1),
                   p_max = 3, n_back = 500)

  # Should find same or very close parameters
  expect_equal(new$u, old$u, tolerance = 0.11)
  expect_equal(new$lambda, old$lambda, tolerance = 0.11)
})
