# ==============================================================================
# Parity tests for gen_aruma_flex ARMA filtering pipeline
#
# These tests capture the exact current behavior of gen_aruma_flex as a
# contract. When the C++ ARMA filter replacement is wired in, these same
# tests must pass. If any test fails, the C++ code is wrong — not the test.
# ==============================================================================


# ==============================================================================
# 1A. Pure AR models — golden values
# ==============================================================================

test_that("parity: AR(1) phi=0.7 golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.7, plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-0.792978598562918, -0.221307821560472,
    1.01640965226646, 2.77102599888583, 0.562856600979559,
    -0.756855944941414, -1.23562055621911, -1.91899017143057,
    -1.98903684314389, -1.57770375787723))
})

test_that("parity: AR(1) phi=0.7 seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.7, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.7, plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: AR(1) phi=0.7 ACF matches theoretical", {
  r <- gen_aruma_flex(n = 5000, phi = 0.7, plot = FALSE, seed = 42)
  acf1 <- acf(r$y, lag.max = 1, plot = FALSE)$acf[2]
  expect_true(abs(acf1 - 0.7) < 0.05)
})

test_that("parity: AR(2) phi=c(0.5,0.3) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.3), plot = FALSE, seed = 99)
  r2 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.3), plot = FALSE, seed = 99)
  expect_identical(r1$y, r2$y)
})

test_that("parity: AR(5) seed reproducibility", {
  phi5 <- c(0.3, 0.2, 0.1, 0.05, 0.02)
  r1 <- gen_aruma_flex(n = 200, phi = phi5, plot = FALSE, seed = 123)
  r2 <- gen_aruma_flex(n = 200, phi = phi5, plot = FALSE, seed = 123)
  expect_identical(r1$y, r2$y)
})

test_that("parity: white noise golden values", {
  r <- gen_aruma_flex(n = 200, plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-2.00092923773151, 0.33377719743357,
    1.17132512735879, 2.0595392422993, -1.37686159824052,
    -1.15085556562711, -0.705821394760121, -1.05405578207719,
    -0.645743723142491, -0.185377967676503))
})


# ==============================================================================
# 1B. Pure MA models — golden values
# ==============================================================================

test_that("parity: MA(1) theta=0.4 golden values", {
  r <- gen_aruma_flex(n = 200, theta = 0.4, plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-2.05245780917246, 1.13414889252617,
    1.03781424838537, 1.59100919135579, -2.20067729516024,
    -0.600110926330897, -0.245479168509279, -0.771727224173143,
    -0.224121410311615, 0.0729195215804932))
})

test_that("parity: MA(1) theta=0.4 seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, theta = 0.4, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, theta = 0.4, plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: MA(2) theta=c(0.4,0.2) golden values", {
  r <- gen_aruma_flex(n = 200, theta = c(0.4, 0.2), plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-2.41550349840326, 1.1083846068057,
    1.43800009593167, 1.52425375186907, -2.434942320632,
    -1.01201877479076, 0.0298931511388255, -0.541556111047722,
    -0.0829571313595906, 0.283730677995931))
})

test_that("parity: MA(2) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, theta = c(0.4, 0.2), plot = FALSE, seed = 77)
  r2 <- gen_aruma_flex(n = 200, theta = c(0.4, 0.2), plot = FALSE, seed = 77)
  expect_identical(r1$y, r2$y)
})

test_that("parity: MA(1) has near-zero ACF at lag 2+", {
  # MA(1) has non-zero ACF only at lag 1
  r <- gen_aruma_flex(n = 5000, theta = 0.4, plot = FALSE, seed = 42)
  acf_vals <- acf(r$y, lag.max = 5, plot = FALSE)$acf[-1]
  # lag-1 should be non-zero
  expect_true(abs(acf_vals[1]) > 0.1)
  # lag-2 through lag-5 should be near zero
  expect_true(all(abs(acf_vals[2:5]) < 0.05))
})


# ==============================================================================
# 1C. ARMA models — golden values and ACF
# ==============================================================================

test_that("parity: ARMA(1,1) phi=0.7 theta=0.3 golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3,
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-1.31067172963517, 0.0165857580084031,
    1.08280199873461, 2.46610310320589, -0.26845119868619,
    -0.925712925235282, -1.00856377273669, -1.54830400456484,
    -1.41333979171472, -0.980992704934059))
})

test_that("parity: ARMA(1,1) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: ARMA(1,1) ACF matches theoretical", {
  r <- gen_aruma_flex(n = 5000, phi = 0.7, theta = 0.3,
                      plot = FALSE, seed = 42)
  acf1 <- acf(r$y, lag.max = 1, plot = FALSE)$acf[2]
  # Theoretical ARMA(1,1) lag-1 ACF using arima.sim convention:
  # gen_aruma_flex negates theta (line 150: model$ma <- -theta),

  # so arima.sim uses ma = -0.3. The formula uses theta_arima = -theta_pkg.
  # rho(1) = (1 + phi*ma)(phi + ma) / (1 + ma^2 + 2*phi*ma)
  ma <- -0.3  # arima.sim convention
  theo_rho1 <- (1 + 0.7 * ma) * (0.7 + ma) /
               (1 + ma^2 + 2 * 0.7 * ma)
  expect_true(abs(acf1 - theo_rho1) < 0.05)
})

test_that("parity: ARMA(2,1) phi=c(0.5,0.3) theta=0.2 golden values", {
  r <- gen_aruma_flex(n = 200, phi = c(0.5, 0.3), theta = 0.2,
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-1.18277747038927, 0.345861203502622,
    0.922667048506611, 2.39036610213164, -0.31678628108258,
    -0.3167665558808, -0.729069443899874, -1.37245619183934,
    -1.33988149581669, -1.13790682850815))
})

test_that("parity: ARMA(1,2) phi=0.6 theta=c(0.3,0.1) golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.6, theta = c(0.3, 0.1),
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-1.77762137815162, -0.145399000998187,
    1.18404549130296, 2.38519127913009, -0.68074111618814,
    -1.35219568009776, -1.0341959733066, -1.3477413910704,
    -1.06758968368556, -0.526803082737375))
})

test_that("parity: ARMA(2,2) phi=c(0.5,0.2) theta=c(0.3,0.1) golden values", {
  r <- gen_aruma_flex(n = 200, phi = c(0.5, 0.2), theta = c(0.3, 0.1),
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-1.54811329651523, 0.264133999713241,
    1.09372923245545, 2.27445540051868, -0.755882336915761,
    -0.866801098739024, -0.807455582000601, -1.30431181783455,
    -1.07259187436072, -0.683407573273305))
})

test_that("parity: ARMA(2,1) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.3), theta = 0.2,
                       plot = FALSE, seed = 55)
  r2 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.3), theta = 0.2,
                       plot = FALSE, seed = 55)
  expect_identical(r1$y, r2$y)
})

test_that("parity: ARMA(2,2) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.2), theta = c(0.3, 0.1),
                       plot = FALSE, seed = 88)
  r2 <- gen_aruma_flex(n = 200, phi = c(0.5, 0.2), theta = c(0.3, 0.1),
                       plot = FALSE, seed = 88)
  expect_identical(r1$y, r2$y)
})


# ==============================================================================
# 1D. With differencing (d > 0) — golden values and stationarity
# ==============================================================================

test_that("parity: ARIMA(1,1,0) golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-19.9955892930865, -20.3496781884807,
    -19.355397508819, -16.7987179266888, -16.8972397338642,
    -18.097356203079, -19.4032358324466, -21.1102314292075,
    -22.6094729507305, -23.5444716791685))
})

test_that("parity: ARIMA(1,1,0) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: ARIMA(1,1,0) differenced series is stationary", {
  r <- gen_aruma_flex(n = 500, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  diffs <- diff(r$y)
  # ACF of differenced series should decay (stationary AR(1))
  acf_vals <- acf(diffs, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(all(abs(acf_vals) < 1))
  # Lag-1 ACF should be close to phi=0.5
  expect_true(abs(acf_vals[1] - 0.5) < 0.1)
})

test_that("parity: ARIMA(1,1,1) golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.5, theta = 0.3, d = 1,
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-14.4478156902084, -14.3891849299059,
    -13.2886775816259, -11.0302822033943, -11.8958078852088,
    -13.066367812271, -14.0122125008741, -15.3274442088248,
    -16.3145870513194, -16.7998133233005))
})

test_that("parity: ARIMA(1,1,1) seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, theta = 0.3, d = 1,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, theta = 0.3, d = 1,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})


# ==============================================================================
# 1E. With seasonal + lambda factors — golden values
# ==============================================================================

test_that("parity: AR(1) + s=12 golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.5, s = 12, plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-5.58410374188071, -0.341457586063415,
    -0.515349057098586, 0.31916825098466, -2.77559292088492,
    -2.07885270872757, -2.23234685234362, -2.4760298407872,
    -3.13620333049155, -3.02441111182746))
})

test_that("parity: AR(1) + s=12 seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, s = 12, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, s = 12, plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: AR(1) + lambda=0.9 golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.5, lambda = 0.9,
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-0.399170384500014, 0.635027333611713,
    3.12820418238071, 2.7168619569672, 1.24505929205566,
    -0.185326266517444, -1.87378923662666, -3.18565183448696,
    -3.80208537947625, -5.09059825648761))
})

test_that("parity: AR(1) + lambda=0.9 seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, lambda = 0.9,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, lambda = 0.9,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: full ARUMA (phi + d + s + lambda) golden values", {
  r <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, s = 4, lambda = 0.8,
                      plot = FALSE, seed = 42)
  expect_equal(r$y[1:10], c(-1617.14011231317, -1639.76922669278,
    -1663.3947967426, -1694.37337815755, -1720.84100654388,
    -1747.94313517148, -1773.94450516018, -1806.11464605161,
    -1833.2650898823, -1860.28263515537))
})

test_that("parity: full ARUMA seed reproducibility", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, s = 4, lambda = 0.8,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, s = 4, lambda = 0.8,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})


# ==============================================================================
# 1F. Innovation generator independence
# ==============================================================================

test_that("parity: different innovation generators produce different output", {
  r_norm <- gen_aruma_flex(n = 200, phi = 0.5, plot = FALSE, seed = 42)
  t_gen <- make_gen_t(df = 5, scale = TRUE)
  r_t <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = t_gen,
                        plot = FALSE, seed = 42)
  # Different distributions → different values

  expect_false(identical(r_norm$y, r_t$y))
  # But same structural properties (length, no NA)
  expect_length(r_t$y, 200)
  expect_false(any(is.na(r_t$y)))
})

test_that("parity: t-innovations seed reproducibility", {
  t_gen <- make_gen_t(df = 5, scale = TRUE)
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = t_gen,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = t_gen,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: laplace-innovations seed reproducibility", {
  lap_gen <- make_gen_laplace()
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = lap_gen,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = lap_gen,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: mixture-innovations seed reproducibility", {
  mix_gen <- make_gen_mixnorm(sd1 = 1, sd2 = 3, prob1 = 0.9)
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = mix_gen,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = mix_gen,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: ARMA(1,1) with t-innovations seed reproducibility", {
  t_gen <- make_gen_t(df = 3, scale = TRUE)
  r1 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3, innov_gen = t_gen,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3, innov_gen = t_gen,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})

test_that("parity: hetero-innovations seed reproducibility", {
  h_gen <- make_gen_hetero(shape = "linear", from = 1, to = 5)
  r1 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = h_gen,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = h_gen,
                       plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
})


# ==============================================================================
# 1G. Structural property tests
# ==============================================================================

test_that("parity: all configs produce no NA", {
  configs <- list(
    list(phi = 0.7),
    list(theta = 0.4),
    list(phi = 0.7, theta = 0.3),
    list(phi = 0.5, d = 1),
    list(phi = c(0.5, 0.3), theta = 0.2),
    list(phi = 0.6, theta = c(0.3, 0.1)),
    list(phi = c(0.5, 0.2), theta = c(0.3, 0.1)),
    list(phi = 0.5, s = 12),
    list(phi = 0.5, lambda = 0.9),
    list(phi = 0.5, d = 1, s = 4, lambda = 0.8)
  )
  for (cfg in configs) {
    args <- c(list(n = 100, plot = FALSE, seed = 42), cfg)
    r <- do.call(gen_aruma_flex, args)
    expect_length(r$y, 100)
    expect_false(any(is.na(r$y)), info = paste(names(cfg), collapse = "+"))
  }
})

test_that("parity: output length correct for all configs", {
  for (n in c(10, 50, 200, 1000)) {
    r <- gen_aruma_flex(n = n, phi = 0.7, theta = 0.3,
                        plot = FALSE, seed = 42)
    expect_length(r$y, n)
  }
})

test_that("parity: different seeds produce different output", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3,
                       plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3,
                       plot = FALSE, seed = 43)
  expect_false(identical(r1$y, r2$y))
})

test_that("parity: AR(1) near unit root (phi=0.99) works", {
  r1 <- gen_aruma_flex(n = 200, phi = 0.99, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 200, phi = 0.99, plot = FALSE, seed = 42)
  expect_identical(r1$y, r2$y)
  expect_length(r1$y, 200)
  expect_false(any(is.na(r1$y)))
})

test_that("parity: small n works for ARMA", {
  r <- gen_aruma_flex(n = 5, phi = 0.7, theta = 0.3,
                      plot = FALSE, seed = 42)
  expect_length(r$y, 5)
  expect_false(any(is.na(r$y)))
})


# ==============================================================================
# 1H. Innovation count verification
# ==============================================================================

test_that("parity: innovation generator called with correct n", {
  # Track how many innovations are requested
  call_log <- integer(0)
  counting_gen <- function(n) {
    call_log[length(call_log) + 1L] <<- n
    rnorm(n)
  }

  # AR(1): p=1, q=0, n_start = max(100, 10*1) = 100
  # n_sim = 200 + 0 + 100 (spin) = 300 (for dlams=0)
  # Wait, spin = max(100, 2*dlams) = max(100, 0) = 100
  # total_innov = n_sim + n_start = 300 + 100 = 400
  r <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = counting_gen,
                      plot = FALSE, seed = 42)
  expect_length(call_log, 1)  # Called exactly once
  expect_true(call_log[1] > 200)  # More than n (includes burn-in)
})
