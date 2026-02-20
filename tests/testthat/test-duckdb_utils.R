# Helper: create a fake wbg_boot_fast result
mock_boot_result <- function(obs_stat = 1.5, pvalue = 0.04, nb = 99) {
  list(
    p = 2L,
    phi = c(0.6, -0.1),
    vara = 1.2,
    pvalue = pvalue,
    pvalue_upper = pvalue / 2,
    pvalue_lower = 1 - pvalue / 2,
    pvalue_asymp = 0.03,
    pvalue_adj = NULL,
    tco_obs = obs_stat,
    boot_tstats = rnorm(nb),
    n = 100L,
    nb = nb
  )
}


test_that("mc_db_connect works with in-memory database", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  expect_true(DBI::dbIsValid(con))
  DBI::dbDisconnect(con)
})


test_that("mc_db_init creates tables and views (idempotent)", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  mc_db_init(con)

  tables <- DBI::dbListTables(con)
  expect_true("batches" %in% tables)
  expect_true("simulations" %in% tables)
  expect_true("v_rejection_rates" %in% tables)
  expect_true("v_rejection_rates_by_batch" %in% tables)

  # Idempotent: calling again should not error
  expect_no_error(mc_db_init(con))

  # Tables still exist
  tables2 <- DBI::dbListTables(con)
  expect_true("batches" %in% tables2)
  expect_true("simulations" %in% tables2)

  DBI::dbDisconnect(con)
})


test_that("mc_db_write_batch writes correct rows with field mapping", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  mc_db_init(con)

  batch_id <- .mc_create_batch(con, "test batch")
  expect_true(is.integer(batch_id))
  expect_true(batch_id >= 1)

  set.seed(42)
  results <- lapply(1:10, function(i) mock_boot_result(obs_stat = i * 0.5))

  mc_db_write_batch(con, results, batch_id, n = 100, phi = 0.7, innov_dist = "norm")

  # Verify row count
  count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM simulations")$n
  expect_equal(count, 10)

  # Verify field mapping
  row1 <- DBI::dbGetQuery(con, "SELECT * FROM simulations WHERE iteration = 1")
  expect_equal(row1$n, 100L)
  expect_equal(row1$phi, 0.7)
  expect_equal(row1$innov_dist, "norm")
  expect_equal(row1$obs_stat, 0.5)  # first result: 1 * 0.5
  expect_equal(row1$null_ar_order, 2L)
  expect_equal(row1$batch_id, batch_id)

  DBI::dbDisconnect(con)
})


test_that("batch isolation: pooled vs per-batch views", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  mc_db_init(con)

  # Batch 1: 10 sims, all reject (pvalue < 0.05)
  batch1 <- .mc_create_batch(con, "run 1")
  set.seed(1)
  results1 <- lapply(1:10, function(i) mock_boot_result(pvalue = 0.01))
  mc_db_write_batch(con, results1, batch1, n = 50, phi = 0.8, innov_dist = "t(3)")

  # Batch 2: 10 sims, none reject (pvalue > 0.05)
  batch2 <- .mc_create_batch(con, "run 2")
  set.seed(2)
  results2 <- lapply(1:10, function(i) mock_boot_result(pvalue = 0.50))
  mc_db_write_batch(con, results2, batch2, n = 50, phi = 0.8, innov_dist = "t(3)")

  # Pooled: one row, 20 sims, 2 batches, 50% rejection rate
  pooled <- DBI::dbGetQuery(con,
    "SELECT * FROM v_rejection_rates WHERE n = 50 AND phi = 0.8 AND innov_dist = 't(3)'")
  expect_equal(nrow(pooled), 1)
  expect_equal(pooled$n_sims, 20)
  expect_equal(pooled$n_batches, 2)
  expect_equal(pooled$reject_05, 0.5)

  # Per-batch: two rows
  by_batch <- DBI::dbGetQuery(con,
    "SELECT * FROM v_rejection_rates_by_batch WHERE n = 50 AND phi = 0.8 AND innov_dist = 't(3)'
     ORDER BY batch_id")
  expect_equal(nrow(by_batch), 2)
  expect_equal(by_batch$reject_05[1], 1.0)   # batch 1: all reject
  expect_equal(by_batch$reject_05[2], 0.0)   # batch 2: none reject
  expect_equal(by_batch$batch_label[1], "run 1")
  expect_equal(by_batch$batch_label[2], "run 2")

  DBI::dbDisconnect(con)
})


test_that("mc_study closure round-trip", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  study <- mc_study(":memory:", label = "test session")
  expect_s3_class(study, "mc_study")
  expect_equal(study$info$label, "test session")
  expect_true(study$info$batch_id >= 1)

  set.seed(100)
  results <- lapply(1:5, function(i) mock_boot_result(pvalue = 0.04))
  study$save_batch(results, n = 200, phi = 0.5, innov_dist = "norm")

  # Query
  rows <- study$query(n = 200, phi = 0.5)
  expect_equal(nrow(rows), 5)

  # Rejection rates
  rates <- study$rejection_rates()
  expect_equal(nrow(rates), 1)
  expect_equal(rates$n_sims, 5)
  expect_equal(rates$reject_05, 1.0)  # all have pvalue = 0.04 < 0.05

  study$end()
})


test_that("mc_db_recalc_pvalue recomputes from stored boot_dist", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  mc_db_init(con)

  batch_id <- .mc_create_batch(con, "recalc test")

  # Create a result with known boot_dist
  result <- mock_boot_result(obs_stat = 2.0, nb = 99)
  result$boot_tstats <- c(rep(0.5, 90), rep(3.0, 9))  # 9 of 99 have |t| >= |2.0|

  mc_db_write_batch(con, list(result), batch_id, n = 100, phi = 0.7, innov_dist = "norm")

  # Get sim_id
  sim_id <- DBI::dbGetQuery(con, "SELECT sim_id FROM simulations LIMIT 1")$sim_id

  recalc <- mc_db_recalc_pvalue(con, sim_id)
  expect_equal(recalc$nb, 99)
  # (9 + 1) / (99 + 1) = 0.10
  expect_equal(recalc$pvalue, 0.10)

  DBI::dbDisconnect(con)
})


test_that("transactional rollback on bad data", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- mc_db_connect(":memory:")
  mc_db_init(con)

  batch_id <- .mc_create_batch(con, "rollback test")

  # Create results where one is invalid (missing tco_obs)
  set.seed(1)
  good_result <- mock_boot_result()
  bad_result <- list()  # completely empty -- will fail on r$tco_obs

  expect_error(
    mc_db_write_batch(con, list(good_result, bad_result), batch_id,
                      n = 100, phi = 0.7, innov_dist = "norm"),
    "Batch write failed"
  )

  # No partial rows committed
  count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM simulations")$n
  expect_equal(count, 0)

  DBI::dbDisconnect(con)
})
