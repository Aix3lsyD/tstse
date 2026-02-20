#' DuckDB Utilities for Monte Carlo Simulation Studies
#'
#' Functions for storing and querying bootstrap simulation results in a
#' DuckDB database. Flat two-table schema optimized for a single capstone
#' study analyzing WBG bootstrap rejection rates.
#'
#' @name duckdb_utils
#' @keywords internal
NULL


#' Connect to Simulation Database
#'
#' Opens a connection to a DuckDB database file. Creates the file if it
#' doesn't exist, appends if it does.
#'
#' @param path Character. Path to the DuckDB file. Default is "simulations.duckdb".
#' @param read_only Logical. If TRUE, open in read-only mode. Default is FALSE.
#'
#' @return A DBI connection object.
#'
#' @examples
#' \dontrun{
#' con <- mc_db_connect("capstone.duckdb")
#' mc_db_init(con)
#' DBI::dbDisconnect(con)
#' }
#'
#' @seealso [mc_db_init()]
#' @export
mc_db_connect <- function(path = "simulations.duckdb", read_only = FALSE) {
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("Package 'duckdb' is required. Install with: install.packages('duckdb')",
         call. = FALSE)
  }
  duckdb::dbConnect(duckdb::duckdb(), dbdir = path, read_only = read_only)
}


#' Initialize Simulation Database Schema
#'
#' Creates the tables, indexes, and views needed for storing simulation
#' results. Safe to call on an existing database (idempotent).
#'
#' @param con A DBI connection from [mc_db_connect()].
#'
#' @return Invisible NULL.
#'
#' @details
#' Creates two tables:
#' \itemize{
#'   \item \code{batches}: Execution metadata (auto-increment ID, label, timestamp)
#'   \item \code{simulations}: Flat results with DGP parameters and bootstrap output
#' }
#'
#' And two views:
#' \itemize{
#'   \item \code{v_rejection_rates}: Pooled rejection rates across all batches
#'   \item \code{v_rejection_rates_by_batch}: Per-batch rejection rates for comparison
#' }
#'
#' @examples
#' \dontrun{
#' con <- mc_db_connect(":memory:")
#' mc_db_init(con)
#' }
#'
#' @seealso [mc_db_connect()]
#' @export
mc_db_init <- function(con) {
  # ============================================================================
  # Batches table -- execution metadata
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE SEQUENCE IF NOT EXISTS seq_batch_id START 1;
  ")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS batches (
      batch_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_batch_id'),
      label       VARCHAR,
      created_at  TIMESTAMP DEFAULT current_timestamp
    )
  ")

  # ============================================================================
  # Simulations table -- flat results with DGP parameters
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE SEQUENCE IF NOT EXISTS seq_sim_id START 1;
  ")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS simulations (
      sim_id          INTEGER PRIMARY KEY DEFAULT nextval('seq_sim_id'),
      batch_id        INTEGER NOT NULL REFERENCES batches(batch_id),
      iteration       INTEGER NOT NULL,

      -- DGP parameters
      n               INTEGER NOT NULL,
      phi             DOUBLE NOT NULL,
      innov_dist      VARCHAR NOT NULL,

      -- Test results
      obs_stat        DOUBLE NOT NULL,
      boot_dist       DOUBLE[],
      pvalue          DOUBLE,
      pvalue_upper    DOUBLE,
      pvalue_lower    DOUBLE,
      pvalue_asymp    DOUBLE,
      pvalue_adj      DOUBLE,

      -- Null model fitted
      null_ar_order   INTEGER,
      null_ar_phi     DOUBLE[],
      null_vara       DOUBLE,

      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  # Indexes
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_sim_config ON simulations(n, phi, innov_dist)
  ")
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_sim_batch ON simulations(batch_id)
  ")

  # ============================================================================
  # View: Pooled rejection rates across all batches
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_rejection_rates AS
    SELECT
      s.n,
      s.phi,
      s.innov_dist,
      COUNT(*) AS n_sims,
      COUNT(DISTINCT s.batch_id) AS n_batches,

      -- Bootstrap rejection rate
      AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_05,
      SQRT(AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_05_se,

      -- Asymptotic rejection rate
      AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_asymp_05,
      SQRT(AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_asymp_05_se,

      -- COBA adjusted rejection rate
      AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_adj_05,
      SQRT(AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_adj_05_se

    FROM simulations s
    GROUP BY s.n, s.phi, s.innov_dist
  ")

  # ============================================================================
  # View: Per-batch rejection rates for comparison
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_rejection_rates_by_batch AS
    SELECT
      s.batch_id,
      b.label AS batch_label,
      b.created_at AS batch_created_at,
      s.n,
      s.phi,
      s.innov_dist,
      COUNT(*) AS n_sims,

      -- Bootstrap rejection rate
      AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_05,
      SQRT(AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_05_se,

      -- Asymptotic rejection rate
      AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_asymp_05,
      SQRT(AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_asymp_05_se,

      -- COBA adjusted rejection rate
      AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_adj_05,
      SQRT(AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN s.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_adj_05_se

    FROM simulations s
    JOIN batches b ON s.batch_id = b.batch_id
    GROUP BY s.batch_id, b.label, b.created_at, s.n, s.phi, s.innov_dist
  ")

  invisible(NULL)
}


#' Create a New Batch
#'
#' @param con A DBI connection.
#' @param label Optional batch label.
#' @return The auto-incremented batch_id (integer).
#' @noRd
.mc_create_batch <- function(con, label = NULL) {
  label_val <- if (is.null(label)) NA_character_ else label
  DBI::dbExecute(con,
    "INSERT INTO batches (label) VALUES (?)",
    params = list(label_val))
  # Retrieve the last inserted batch_id
  result <- DBI::dbGetQuery(con,
    "SELECT currval('seq_batch_id') AS batch_id")
  as.integer(result$batch_id[1])
}


#' Write Batch of Bootstrap Results to Database
#'
#' Stores multiple bootstrap results in a single transaction.
#'
#' @param con A DBI connection from [mc_db_connect()].
#' @param results A list of `wbg_boot_fast` result objects.
#' @param batch_id Integer. Batch ID from \code{.mc_create_batch()}.
#' @param n Integer. Observation length.
#' @param phi Numeric. AR(1) coefficient.
#' @param innov_dist Character. Innovation distribution (e.g., "norm", "t(3)").
#'
#' @return Invisible NULL.
#'
#' @details
#' Extracts fields directly from `wbg_boot_fast` output:
#' \itemize{
#'   \item \code{tco_obs} -> obs_stat
#'   \item \code{boot_tstats} -> boot_dist
#'   \item \code{phi} -> null_ar_phi
#'   \item \code{p} -> null_ar_order
#'   \item \code{vara} -> null_vara
#'   \item p-value fields map 1:1
#' }
#'
#' @examples
#' \dontrun{
#' con <- mc_db_connect(":memory:")
#' mc_db_init(con)
#' batch_id <- .mc_create_batch(con, "test batch")
#' results <- lapply(1:10, function(i) {
#'   x <- gen_arma(100, phi = 0.7, plot = FALSE)
#'   wbg_boot_fast(x, nb = 99, seed = i)
#' })
#' mc_db_write_batch(con, results, batch_id, n = 100, phi = 0.7, innov_dist = "norm")
#' }
#'
#' @export
mc_db_write_batch <- function(con, results, batch_id, n, phi, innov_dist) {
  if (length(results) == 0) return(invisible(NULL))

  na_real <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)

  DBI::dbBegin(con)
  tryCatch({
    for (i in seq_along(results)) {
      r <- results[[i]]

      # Extract fields from wbg_boot_fast result
      obs_stat <- r$tco_obs
      boot_dist_str <- .vec_to_array(r$boot_tstats)
      null_ar_order <- as.integer(r$p)
      null_ar_phi_str <- if (length(r$phi) > 0) .vec_to_array(r$phi) else "NULL"
      null_vara <- na_real(r$vara)

      DBI::dbExecute(con, sprintf("
        INSERT INTO simulations (batch_id, iteration, n, phi, innov_dist,
                                  obs_stat, boot_dist, pvalue, pvalue_upper,
                                  pvalue_lower, pvalue_asymp, pvalue_adj,
                                  null_ar_order, null_ar_phi, null_vara)
        VALUES (?, ?, ?, ?, ?, ?, %s, ?, ?, ?, ?, ?, ?, %s, ?)",
        boot_dist_str, null_ar_phi_str),
        params = list(
          as.integer(batch_id), as.integer(i), as.integer(n), phi, innov_dist,
          obs_stat,
          na_real(r$pvalue), na_real(r$pvalue_upper), na_real(r$pvalue_lower),
          na_real(r$pvalue_asymp), na_real(r$pvalue_adj),
          null_ar_order, null_vara))
    }
    DBI::dbCommit(con)
  }, error = function(err) {
    DBI::dbRollback(con)
    stop("Batch write failed: ", conditionMessage(err), call. = FALSE)
  })

  invisible(NULL)
}


#' Query Simulation Runs
#'
#' Retrieves simulation runs with optional filtering.
#'
#' @param con A DBI connection from [mc_db_connect()].
#' @param n Integer. Filter by observation length.
#' @param phi Numeric. Filter by AR coefficient.
#' @param innov_dist Character. Filter by innovation distribution.
#' @param batch_id Integer. Filter by batch.
#' @param limit Integer. Maximum rows to return.
#'
#' @return A data frame of simulation results.
#'
#' @export
mc_db_query <- function(con, n = NULL, phi = NULL, innov_dist = NULL,
                        batch_id = NULL, limit = NULL) {
  sql <- "SELECT * FROM simulations WHERE 1=1"
  params <- list()

  if (!is.null(n)) {
    sql <- paste(sql, "AND n = ?")
    params <- c(params, as.integer(n))
  }
  if (!is.null(phi)) {
    sql <- paste(sql, "AND phi = ?")
    params <- c(params, as.numeric(phi))
  }
  if (!is.null(innov_dist)) {
    sql <- paste(sql, "AND innov_dist = ?")
    params <- c(params, innov_dist)
  }
  if (!is.null(batch_id)) {
    sql <- paste(sql, "AND batch_id = ?")
    params <- c(params, as.integer(batch_id))
  }
  if (!is.null(limit)) {
    sql <- paste(sql, "LIMIT", as.integer(limit))
  }

  DBI::dbGetQuery(con, sql, params = params)
}


#' Get Rejection Rates
#'
#' Retrieves aggregated rejection rates, either pooled across all batches
#' or per-batch for comparison.
#'
#' @param con A DBI connection from [mc_db_connect()].
#' @param by_batch Logical. If TRUE, return per-batch rates. Default FALSE (pooled).
#' @param n Integer. Optional filter by observation length.
#' @param phi Numeric. Optional filter by AR coefficient.
#' @param innov_dist Character. Optional filter by innovation distribution.
#'
#' @return A data frame with rejection rates.
#'
#' @export
mc_db_rejection_rates <- function(con, by_batch = FALSE, n = NULL,
                                  phi = NULL, innov_dist = NULL) {
  view_name <- if (by_batch) "v_rejection_rates_by_batch" else "v_rejection_rates"
  sql <- sprintf("SELECT * FROM %s WHERE 1=1", view_name)
  params <- list()

  if (!is.null(n)) {
    sql <- paste(sql, "AND n = ?")
    params <- c(params, as.integer(n))
  }
  if (!is.null(phi)) {
    sql <- paste(sql, "AND phi = ?")
    params <- c(params, as.numeric(phi))
  }
  if (!is.null(innov_dist)) {
    sql <- paste(sql, "AND innov_dist = ?")
    params <- c(params, innov_dist)
  }

  DBI::dbGetQuery(con, sql, params = params)
}


#' Recalculate P-value from Stored Bootstrap Distribution
#'
#' Recomputes p-values using stored bootstrap statistics for a given simulation.
#'
#' @param con A DBI connection from [mc_db_connect()].
#' @param sim_id Integer. The simulation ID.
#'
#' @return A list with `pvalue`, `pvalue_upper`, `pvalue_lower`, and `nb`.
#'
#' @export
mc_db_recalc_pvalue <- function(con, sim_id) {
  result <- DBI::dbGetQuery(con,
    "SELECT obs_stat, boot_dist FROM simulations WHERE sim_id = ?",
    params = list(as.integer(sim_id)))

  if (nrow(result) == 0) {
    stop("Simulation not found: ", sim_id, call. = FALSE)
  }

  obs_stat <- result$obs_stat[1]
  boot_dist <- result$boot_dist[[1]]
  nb <- length(boot_dist)

  pvalue <- (sum(abs(boot_dist) >= abs(obs_stat)) + 1) / (nb + 1)
  pvalue_upper <- (sum(boot_dist >= obs_stat) + 1) / (nb + 1)
  pvalue_lower <- (sum(boot_dist <= obs_stat) + 1) / (nb + 1)

  list(
    pvalue = pvalue,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    nb = nb
  )
}


#' Create a Monte Carlo Study Session
#'
#' Opens a database connection, creates a new batch, and returns a closure
#' with methods for saving results and querying data.
#'
#' @param path Character. Path to the DuckDB database file.
#' @param label Character. Optional label for this batch (e.g., "t(3) n=50 run 2").
#'
#' @return An `mc_study` object (list) with methods:
#' \describe{
#'   \item{`$save_batch(results, n, phi, innov_dist)`}{Save a list of wbg_boot_fast results}
#'   \item{`$query(...)`}{Query simulations with optional filters}
#'   \item{`$rejection_rates(...)`}{Get rejection rates (pooled or by_batch)}
#'   \item{`$end()`}{Close the database connection}
#'   \item{`$info`}{List with path, batch_id, label}
#' }
#'
#' @examples
#' \dontrun{
#' study <- mc_study("capstone.duckdb", label = "t(3) n=50 phi=0.8")
#' on.exit(study$end())
#'
#' results <- lapply(1:1000, function(i) {
#'   x <- gen_arma(50, phi = 0.8, plot = FALSE)
#'   wbg_boot_fast(x, nb = 399, seed = i)
#' })
#' study$save_batch(results, n = 50, phi = 0.8, innov_dist = "t(3)")
#' study$rejection_rates()
#' study$end()
#' }
#'
#' @seealso [mc_db_connect()], [mc_db_write_batch()]
#' @export
mc_study <- function(path, label = NULL) {
  con <- mc_db_connect(path)
  mc_db_init(con)
  batch_id <- .mc_create_batch(con, label)

  structure(
    list(
      save_batch = function(results, n, phi, innov_dist) {
        mc_db_write_batch(con, results, batch_id, n, phi, innov_dist)
        invisible(NULL)
      },

      query = function(...) {
        mc_db_query(con, ...)
      },

      rejection_rates = function(...) {
        mc_db_rejection_rates(con, ...)
      },

      end = function() {
        DBI::dbDisconnect(con)
        invisible(NULL)
      },

      info = list(
        path = path,
        batch_id = batch_id,
        label = label
      )
    ),
    class = "mc_study"
  )
}


#' @export
print.mc_study <- function(x, ...) {
  cat("<mc_study>\n")
  cat("  Database:", x$info$path, "\n")
  cat("  Batch ID:", x$info$batch_id, "\n")
  if (!is.null(x$info$label)) {
    cat("  Label:", x$info$label, "\n")
  }
  cat("Usage: $save_batch(results, n, phi, innov_dist), $query(), $rejection_rates(), $end()\n")
  invisible(x)
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Convert R vector to DuckDB array literal
#' @noRd
.vec_to_array <- function(x) {
  if (is.null(x) || length(x) == 0) return("NULL")
  paste0("[", paste(x, collapse = ", "), "]")
}
