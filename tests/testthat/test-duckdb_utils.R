test_that("boot_db_connect works with in-memory database", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  expect_true(DBI::dbIsValid(con))
  DBI::dbDisconnect(con)
})


test_that("boot_db_init creates tables and views", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Check tables exist
  tables <- DBI::dbListTables(con)
  expect_true("studies" %in% tables)
  expect_true("dgp_configs" %in% tables)
  expect_true("method_configs" %in% tables)
  expect_true("trials" %in% tables)
  expect_true("runs" %in% tables)

  # Check views exist
  expect_true("v_run_summary" %in% tables)
  expect_true("v_rejection_rates" %in% tables)
  expect_true("v_rejection_rates_by_trial" %in% tables)

  DBI::dbDisconnect(con)
})


test_that("boot_db_study creates and retrieves studies", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Create study
  study_id <- boot_db_study(con, "Test Study", "size", "A test study")
  expect_type(study_id, "character")
  expect_true(nchar(study_id) > 0)

  # Retrieve same study (should return same ID)
  study_id2 <- boot_db_study(con, "Test Study", "size")
  expect_equal(study_id, study_id2)

  DBI::dbDisconnect(con)
})


test_that("boot_db_dgp creates DGP configurations and deduplicates", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Create simple AR(1) DGP
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  expect_type(dgp_id, "character")

  # Verify it's in the database
  result <- DBI::dbGetQuery(con,
    "SELECT *, LENGTH(ar_phi) as ar_order FROM dgp_configs WHERE dgp_id = ?",
    params = list(dgp_id))
  expect_equal(nrow(result), 1)
  expect_equal(result$n, 100)
  expect_equal(result$ar_order, 1)

  # Creating same config should return same ID (deduplication)
  dgp_id2 <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  expect_equal(dgp_id, dgp_id2)

  # Different config should get different ID
  dgp_id3 <- boot_db_dgp(con, n = 100, ar_phi = 0.8)
  expect_false(dgp_id == dgp_id3)

  DBI::dbDisconnect(con)
})


test_that("boot_db_method creates method configurations and deduplicates", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  method_id <- boot_db_method(con, "wbg_boot", nb = 399, ar_method = "burg")
  expect_type(method_id, "character")

  result <- DBI::dbGetQuery(con,
    "SELECT * FROM method_configs WHERE method_id = ?",
    params = list(method_id))
  expect_equal(nrow(result), 1)
  expect_equal(result$method_name, "wbg_boot")
  expect_equal(result$nb, 399)
  expect_equal(result$ar_method, "burg")

  # Creating same config should return same ID (deduplication)
  method_id2 <- boot_db_method(con, "wbg_boot", nb = 399, ar_method = "burg")
  expect_equal(method_id, method_id2)

  DBI::dbDisconnect(con)
})


test_that("boot_db_trial creates trials", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 99)

  # Create first trial
  trial1_id <- boot_db_trial(con, study_id, dgp_id, method_id,
                             n_planned = 1000, trial_name = "Trial 1")
  expect_type(trial1_id, "character")

  # Create second trial - should be different (trials are never deduplicated)
  trial2_id <- boot_db_trial(con, study_id, dgp_id, method_id,
                             n_planned = 1000, trial_name = "Trial 2")
  expect_false(trial1_id == trial2_id)

  # Verify both trials exist
  trials <- DBI::dbGetQuery(con, "SELECT * FROM trials")
  expect_equal(nrow(trials), 2)

  DBI::dbDisconnect(con)
})


test_that("boot_db_write stores wbg_boot result", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Create test data and run bootstrap
  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- wbg_boot(x, nb = 49, seed = 123)

  # Create trial
  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)

  # Store result
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  expect_type(run_id, "character")

  # Verify run record (results are now in runs table)
  run <- DBI::dbGetQuery(con, "SELECT * FROM runs WHERE run_id = ?",
                         params = list(run_id))
  expect_equal(nrow(run), 1)
  expect_equal(run$iteration_num, 1)
  expect_equal(run$obs_stat, result$tco_obs)
  expect_equal(run$pvalue, result$pvalue)

  DBI::dbDisconnect(con)
})


test_that("boot_db_query retrieves runs", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Store a result
  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- wbg_boot(x, nb = 49, seed = 123)

  study_id <- boot_db_study(con, "Query Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  boot_db_write(con, result, trial_id, iteration_num = 1)

  # Query all runs
  runs <- boot_db_query(con)
  expect_equal(nrow(runs), 1)
  expect_equal(runs$study_name, "Query Test")

  # Query by study
  runs_study <- boot_db_query(con, study_id = study_id)
  expect_equal(nrow(runs_study), 1)

  # Query by trial
  runs_trial <- boot_db_query(con, trial_id = trial_id)
  expect_equal(nrow(runs_trial), 1)

  DBI::dbDisconnect(con)
})


test_that("boot_db_get_replicates retrieves bootstrap distribution", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- wbg_boot(x, nb = 49, seed = 123)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  # Retrieve replicates
  replicates <- boot_db_get_replicates(con, run_id)

  expect_equal(replicates$obs_stat, result$tco_obs)
  expect_equal(length(replicates$boot_dist), 49)
  expect_equal(replicates$boot_dist, result$boot_tstats)

  DBI::dbDisconnect(con)
})


test_that("boot_db_recalc_pvalue matches original", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- wbg_boot(x, nb = 99, seed = 123)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 99)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  # Recalculate p-value
  recalc <- boot_db_recalc_pvalue(con, run_id)

  expect_equal(recalc$pvalue, result$pvalue)
  expect_equal(recalc$pvalue_upper, result$pvalue_upper)
  expect_equal(recalc$pvalue_lower, result$pvalue_lower)

  DBI::dbDisconnect(con)
})


test_that("boot_db_write handles wbg_boot_flex result", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  stat_fn <- make_stat_co()
  result <- wbg_boot_flex(x, stat_fn, nb = 49, seed = 123, verbose = FALSE,
                          bootadj = FALSE)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_flex", nb = 49,
                               stat_fn_name = "make_stat_co")
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  expect_type(run_id, "character")

  # Verify result stored correctly
  res <- DBI::dbGetQuery(con, "SELECT obs_stat, pvalue FROM runs WHERE run_id = ?",
                         params = list(run_id))
  expect_equal(res$obs_stat, result$obs_stat)
  expect_equal(res$pvalue, result$pvalue)

  DBI::dbDisconnect(con)
})


test_that("boot_db_write handles co_tas_boot result", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- co_tas_boot(x, nb = 49, btest = TRUE, seed = 123)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "co_tas_boot", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  expect_type(run_id, "character")

  res <- DBI::dbGetQuery(con, "SELECT obs_stat, pvalue FROM runs WHERE run_id = ?",
                         params = list(run_id))
  expect_equal(res$obs_stat, result$tco)
  expect_equal(res$pvalue, result$pvalue)

  DBI::dbDisconnect(con)
})


test_that("master_seed is stored and retrievable", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  x <- gen_arma(100, phi = 0.7)
  result <- wbg_boot(x, nb = 49, seed = 12345)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)
  run_id <- boot_db_write(con, result, trial_id, iteration_num = 1)

  # Verify master_seed is stored
  run <- DBI::dbGetQuery(con, "SELECT master_seed FROM runs WHERE run_id = ?",
                         params = list(run_id))
  expect_equal(run$master_seed, 12345)

  DBI::dbDisconnect(con)
})


# ==============================================================================
# Batch write tests
# ==============================================================================

test_that("boot_db_write_batch writes multiple results", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Generate multiple results
  set.seed(42)
  results_list <- lapply(1:10, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i)
  })

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id, n_planned = 10)

  run_ids <- boot_db_write_batch(con, results_list, trial_id)

  # Check return value
  expect_length(run_ids, 10)
  expect_type(run_ids, "character")

  # Check all records were inserted
  n_runs <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM runs")$n
  expect_equal(n_runs, 10)

  # Check iteration numbers are correct
  iters <- DBI::dbGetQuery(con, "SELECT iteration_num FROM runs ORDER BY iteration_num")
  expect_equal(iters$iteration_num, 1:10)

  DBI::dbDisconnect(con)
})


test_that("boot_db_write_batch with store_boot_dist=FALSE skips boot_dist", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  results_list <- lapply(1:5, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i)
  })

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)

  # Write without boot distribution
  run_ids <- boot_db_write_batch(con, results_list, trial_id,
                                  store_boot_dist = FALSE)

  # Check boot_dist is NULL
  res <- DBI::dbGetQuery(con, "SELECT boot_dist FROM runs LIMIT 1")
  expect_true(is.na(res$boot_dist) || is.null(res$boot_dist) ||
              length(res$boot_dist[[1]]) == 0)

  # But pvalues should be stored
  pvals <- DBI::dbGetQuery(con, "SELECT pvalue FROM runs")
  expect_equal(nrow(pvals), 5)
  expect_true(all(!is.na(pvals$pvalue)))

  DBI::dbDisconnect(con)
})


test_that("boot_db_write_batch_lite is equivalent to store_boot_dist=FALSE", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  set.seed(42)
  results_list <- lapply(1:3, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i)
  })

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)

  # Use lite version
  run_ids <- boot_db_write_batch_lite(con, results_list, trial_id)

  expect_length(run_ids, 3)

  # Verify pvalues stored
  pvals <- DBI::dbGetQuery(con, "SELECT pvalue, pvalue_asymp FROM runs")
  expect_equal(nrow(pvals), 3)
  expect_true(all(!is.na(pvals$pvalue)))
  expect_true(all(!is.na(pvals$pvalue_asymp)))

  DBI::dbDisconnect(con)
})


test_that("boot_db_write_batch handles empty list", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  study_id <- boot_db_study(con, "Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id)

  # Empty list should return empty vector
  run_ids <- boot_db_write_batch(con, list(), trial_id)
  expect_length(run_ids, 0)

  DBI::dbDisconnect(con)
})


test_that("boot_db_write_batch stores rejection rates correctly", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  # Generate results with bootadj=TRUE to get COBA
  set.seed(42)
  results_list <- lapply(1:20, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, bootadj = TRUE, seed = i)
  })

  study_id <- boot_db_study(con, "Batch Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49, bootadj = TRUE)
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id, n_planned = 20)

  boot_db_write_batch(con, results_list, trial_id, store_boot_dist = FALSE)

  # Complete trial so it shows up in rejection rates view
  boot_db_trial_complete(con, trial_id)

  # Check rejection rates view works
  rates <- boot_db_rejection_rates(con)
  expect_equal(nrow(rates), 1)
  expect_equal(rates$n_runs, 20)
  expect_true(!is.na(rates$reject_co_05))
  expect_true(!is.na(rates$reject_cob_05))
  expect_true(!is.na(rates$reject_coba_05))

  DBI::dbDisconnect(con)
})


# ==============================================================================
# Trial workflow tests
# ==============================================================================

test_that("multiple trials aggregate correctly", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  con <- boot_db_connect(":memory:")
  boot_db_init(con)

  study_id <- boot_db_study(con, "Multi-Trial Test", "size")
  dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
  method_id <- boot_db_method(con, "wbg_boot_fast", nb = 49)

  # Trial 1: 10 simulations
  trial1_id <- boot_db_trial(con, study_id, dgp_id, method_id,
                             n_planned = 10, trial_name = "Trial 1")
  set.seed(100)
  results1 <- lapply(1:10, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i)
  })
  boot_db_write_batch(con, results1, trial1_id, store_boot_dist = FALSE)
  boot_db_trial_complete(con, trial1_id)

  # Trial 2: 10 simulations
  trial2_id <- boot_db_trial(con, study_id, dgp_id, method_id,
                             n_planned = 10, trial_name = "Trial 2")
  set.seed(200)
  results2 <- lapply(1:10, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i + 100)
  })
  boot_db_write_batch(con, results2, trial2_id, store_boot_dist = FALSE)
  boot_db_trial_complete(con, trial2_id)

  # Aggregated view should show 20 runs across 2 trials
  rates <- boot_db_rejection_rates(con)
  expect_equal(nrow(rates), 1)
  expect_equal(rates$n_runs, 20)
  expect_equal(rates$n_trials, 2)

  # By-trial view should show 2 rows
  rates_by_trial <- boot_db_rejection_rates(con, by_trial = TRUE)
  expect_equal(nrow(rates_by_trial), 2)
  expect_equal(sum(rates_by_trial$n_runs), 20)

  DBI::dbDisconnect(con)
})


test_that("boot_study wrapper creates trial and tracks iterations", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  study <- boot_study(":memory:", "Wrapper Test",
                      dgp = list(n = 100, ar_phi = 0.7),
                      method = list(name = "wbg_boot_fast", nb = 49),
                      type = "size",
                      n_planned = 5,
                      trial_name = "Auto Trial")

  # Save some results
  set.seed(42)
  for (i in 1:5) {
    x <- gen_arma(100, phi = 0.7)
    result <- wbg_boot_fast(x, nb = 49, seed = i)
    study$save(result)
  }

  # Complete and query
  study$complete()
  runs <- study$query()

  expect_equal(nrow(runs), 5)
  expect_equal(sort(runs$iteration_num), 1:5)
  expect_equal(unique(runs$trial_name), "Auto Trial")

  study$end()
})


test_that("boot_study save_batch works", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  study <- boot_study(":memory:", "Batch Test",
                      dgp = list(n = 100, ar_phi = 0.7),
                      method = list(name = "wbg_boot_fast", nb = 49),
                      type = "size")

  # Save batch of results
  set.seed(42)
  results <- lapply(1:10, function(i) {
    x <- gen_arma(100, phi = 0.7)
    wbg_boot_fast(x, nb = 49, seed = i)
  })
  study$save_batch(results)

  # Query
  runs <- study$query()
  expect_equal(nrow(runs), 10)
  expect_equal(sort(runs$iteration_num), 1:10)

  study$end()
})
