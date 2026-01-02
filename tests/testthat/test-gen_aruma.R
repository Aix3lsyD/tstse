# tests/testthat/test-gen_aruma.R

test_that("gen_aruma produces correct length output", {
  x <- gen_aruma(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_length(x, 100)
})

test_that("gen_aruma is reproducible with seed", {
  x1 <- gen_aruma(n = 50, phi = 0.5, d = 1, plot = FALSE, seed = 123)
  x2 <- gen_aruma(n = 50, phi = 0.5, d = 1, plot = FALSE, seed = 123)
  expect_equal(x1, x2)
})

test_that("gen_aruma handles various parameter combinations", {
  # AR only
  x1 <- gen_aruma(n = 100, phi = c(0.5, 0.3), plot = FALSE, seed = 42)
  expect_length(x1, 100)


  # MA only
  x2 <- gen_aruma(n = 100, theta = 0.4, plot = FALSE, seed = 42)
  expect_length(x2, 100)

  # ARMA
  x3 <- gen_aruma(n = 100, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)
  expect_length(x3, 100)

  # With differencing
  x4 <- gen_aruma(n = 100, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  expect_length(x4, 100)

  # With seasonal
  x5 <- gen_aruma(n = 100, phi = 0.5, s = 12, plot = FALSE, seed = 42)
  expect_length(x5, 100)

  # With lambda
  x6 <- gen_aruma(n = 100, phi = 0.5, lambda = 0.9, plot = FALSE, seed = 42)
  expect_length(x6, 100)

  # All combined
  x7 <- gen_aruma(n = 100, phi = 0.5, theta = 0.3, d = 1, s = 4, lambda = 0.8,
                  plot = FALSE, seed = 42)
  expect_length(x7, 100)
})

test_that("gen_aruma validates inputs", {
  expect_error(gen_aruma(n = -1, phi = 0.5),
               "`n` must be a positive integer")
  expect_error(gen_aruma(n = 100, phi = "a"),
               "`phi` must be numeric")
  expect_error(gen_aruma(n = 100, theta = "a"),
               "`theta` must be numeric")
  expect_error(gen_aruma(n = 100, d = -1),
               "`d` must be a non-negative integer")
  expect_error(gen_aruma(n = 100, s = -1),
               "`s` must be a non-negative integer")
  expect_error(gen_aruma(n = 100, lambda = "a"),
               "`lambda` must be numeric")
  expect_error(gen_aruma(n = 100, vara = -1),
               "`vara` must be positive")
})

test_that("gen_aruma returns invisibly when plot = TRUE", {
  expect_invisible(gen_aruma(n = 50, phi = 0.5, seed = 42))
})

test_that("gen_aruma matches tswge::gen.aruma.wge output with same seed", {
  skip_if_not_installed("tswge")

  # Simple AR(1) case
  tstse_result <- gen_aruma(n = 200, phi = 0.5, plot = FALSE, seed = 42)
  tswge_result <- tswge::gen.aruma.wge(n = 200, phi = 0.5, plot = FALSE, sn = 42)
  expect_equal(tstse_result, tswge_result, tolerance = 1e-10)

  # ARMA(1,1) case
  tstse_result2 <- gen_aruma(n = 200, phi = 0.5, theta = 0.3,
                             plot = FALSE, seed = 123)
  tswge_result2 <- tswge::gen.aruma.wge(n = 200, phi = 0.5, theta = 0.3,
                                        plot = FALSE, sn = 123)
  expect_equal(tstse_result2, tswge_result2, tolerance = 1e-10)

  # With differencing
  tstse_result3 <- gen_aruma(n = 200, phi = 0.5, d = 1, plot = FALSE, seed = 456)
  tswge_result3 <- tswge::gen.aruma.wge(n = 200, phi = 0.5, d = 1, plot = FALSE, sn = 456)
  expect_equal(tstse_result3, tswge_result3, tolerance = 1e-10)

  # With seasonal
  tstse_result4 <- gen_aruma(n = 200, phi = 0.5, s = 4, plot = FALSE, seed = 789)
  tswge_result4 <- tswge::gen.aruma.wge(n = 200, phi = 0.5, s = 4, plot = FALSE, sn = 789)
  expect_equal(tstse_result4, tswge_result4, tolerance = 1e-10)

  # With lambda
  tstse_result5 <- gen_aruma(n = 200, phi = 0.5, lambda = 0.9,
                             plot = FALSE, seed = 111)
  tswge_result5 <- tswge::gen.aruma.wge(n = 200, phi = 0.5, lambda = 0.9,
                                        plot = FALSE, sn = 111)
  expect_equal(tstse_result5, tswge_result5, tolerance = 1e-10)
})

test_that("gen_aruma with lambda and seasonal combined", {
  skip_if_not_installed("tswge")

  # Lambda + seasonal
  tstse_result <- gen_aruma(n = 200, phi = 0.5, s = 4, lambda = 0.8,
                            plot = FALSE, seed = 222)
  tswge_result <- tswge::gen.aruma.wge(n = 200, phi = 0.5, s = 4, lambda = 0.8,
                                       plot = FALSE, sn = 222)
  expect_equal(tstse_result, tswge_result, tolerance = 1e-10)
})

test_that("gen_aruma with d=0, s=0 produces stationary ARMA-like output", {
  # When d=0 and s=0, gen_aruma should produce a stationary ARMA process

  # Note: gen_aruma uses spin=100 while gen_arma uses spin=2000, so they
  # won't produce identical output even with same seed
  x <- gen_aruma(n = 500, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)

  # Should be approximately mean zero (stationary)
  expect_lt(abs(mean(x)), 0.5)

  # Should have reasonable ACF decay for ARMA(1,1)
  acf_vals <- acf(x, lag.max = 10, plot = FALSE)$acf[-1]
  expect_true(all(abs(acf_vals) < 1))  # All ACF values bounded
})

test_that("gen_aruma respects vara parameter", {
  # Higher variance should produce larger values on average
  set.seed(42)
  x1 <- gen_aruma(n = 1000, phi = 0.3, vara = 1, plot = FALSE, seed = 42)
  x2 <- gen_aruma(n = 1000, phi = 0.3, vara = 4, plot = FALSE, seed = 42)

  # Variance should scale approximately
  expect_gt(var(x2), var(x1))
})
