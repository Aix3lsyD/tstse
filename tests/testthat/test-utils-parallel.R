# tests/testthat/test-utils-parallel.R

test_that("get_cores returns valid core count", {

  # Explicit values
  expect_equal(get_cores(1), 1L)
  expect_equal(get_cores(2), min(2L, parallel::detectCores(logical = FALSE)))

  # Zero means all
  all_cores <- parallel::detectCores(logical = FALSE)
  if (!is.na(all_cores)) {
    expect_equal(get_cores(0), all_cores)
  }

  # NULL uses option
  old_opt <- getOption("tstse.cores")
  options(tstse.cores = 2)
  expect_equal(get_cores(NULL), min(2L, parallel::detectCores(logical = FALSE)))
  options(tstse.cores = old_opt)

  # Default when no option set
  options(tstse.cores = NULL)
  expect_equal(get_cores(NULL), 1L)
})


test_that("pmap works sequentially", {

  result <- pmap(1:5, function(x) x^2, cores = 1L)
  expect_equal(result, list(1, 4, 9, 16, 25))
})


test_that("pmap works in parallel", {

  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  result <- pmap(1:5, function(x) x^2, cores = 2L)
  expect_equal(result, list(1, 4, 9, 16, 25))
})


test_that("aic_ar produces same results with parallel", {

  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- rnorm(100)

  capture.output({
    result_seq <- aic_ar(x, p = 1:5, cores = 1)
    result_par <- aic_ar(x, p = 1:5, cores = 2)
  })

  expect_equal(result_seq$p, result_par$p)
  expect_equal(result_seq$phi, result_par$phi)
  expect_equal(result_seq$table, result_par$table)
})


test_that("aic produces same results with parallel", {

  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- rnorm(100)

  result_seq <- aic(x, p = 0:2, q = 0:1, cores = 1)
  result_par <- aic(x, p = 0:2, q = 0:1, cores = 2)

  expect_equal(result_seq$p, result_par$p)
  expect_equal(result_seq$q, result_par$q)
  expect_equal(result_seq$table, result_par$table)
})
