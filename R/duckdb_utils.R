#' DuckDB Utilities for Bootstrap Simulation Studies
#'
#' Functions for storing and querying bootstrap simulation results in a
#' DuckDB database. Enables reproducible research, post-hoc analysis,
#' and systematic comparison of bootstrap methods.
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
#' @details
#' Requires the `duckdb` package to be installed. The database file is created
#' if it doesn't exist. Use [boot_db_init()] after connecting to create the
#' schema if this is a new database.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect("my_simulations.duckdb")
#' boot_db_init(con)
#' # ... run simulations ...
#' DBI::dbDisconnect(con)
#' }
#'
#' @seealso [boot_db_init()], [boot_db_write()]
#' @export
boot_db_connect <- function(path = "simulations.duckdb", read_only = FALSE) {
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("Package 'duckdb' is required. Install with: install.packages('duckdb')",
         call. = FALSE)
  }

  duckdb::dbConnect(duckdb::duckdb(), dbdir = path, read_only = read_only)
}


#' Initialize Simulation Database Schema
#'
#' Creates the tables and views needed for storing simulation results.
#' Safe to call on an existing database (idempotent).
#'
#' @param con A DBI connection from [boot_db_connect()].
#'
#' @return Invisible NULL.
#'
#' @details
#' Creates the following tables if they don't exist:
#' \itemize{
#'   \item \code{studies}: Top-level grouping of simulation studies
#'   \item \code{dgp_configs}: Data generating process configurations (deduplicated)
#'   \item \code{method_configs}: Bootstrap method configurations (deduplicated)
#'   \item \code{trials}: Execution batches within a study (e.g., "run 1 of 1000 sims")
#'   \item \code{runs}: Individual simulation runs with results
#' }
#'
#' Also creates views for easy querying:
#' \itemize{
#'   \item \code{v_run_summary}: Denormalized view of all run data
#'   \item \code{v_rejection_rates}: Aggregated rejection rates (across all trials)
#'   \item \code{v_rejection_rates_by_trial}: Rejection rates per trial
#' }
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#' }
#'
#' @seealso [boot_db_connect()]
#' @export
boot_db_init <- function(con) {
  # Install and load JSON extension for native JSON type support
  DBI::dbExecute(con, "INSTALL json;")
  DBI::dbExecute(con, "LOAD json;")

  # ============================================================================
  # Studies table - top-level grouping
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS studies (
      study_id      VARCHAR PRIMARY KEY,
      study_name    VARCHAR NOT NULL UNIQUE,
      study_type    VARCHAR NOT NULL,
      error_type    VARCHAR,
      description   TEXT,
      created_at    TIMESTAMP DEFAULT current_timestamp,
      metadata      JSON
    )
  ")

  # ============================================================================
  # DGP configurations - deduplicated by config_hash
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS dgp_configs (
      dgp_id          VARCHAR PRIMARY KEY,
      config_hash     VARCHAR NOT NULL UNIQUE,
      dgp_name        VARCHAR,
      n               INTEGER NOT NULL,
      ar_phi          DOUBLE[],
      ma_theta        DOUBLE[],
      vara            DOUBLE DEFAULT 1.0,
      garch_omega     DOUBLE,
      garch_alpha     DOUBLE[],
      garch_beta      DOUBLE[],
      innov_dist      VARCHAR DEFAULT 'norm',
      innov_params    JSON,
      has_trend       BOOLEAN DEFAULT FALSE,
      trend_slope     DOUBLE DEFAULT 0.0,
      trend_intercept DOUBLE DEFAULT 0.0,
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  # ============================================================================
  # Method configurations - deduplicated by config_hash
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS method_configs (
      method_id       VARCHAR PRIMARY KEY,
      config_hash     VARCHAR NOT NULL UNIQUE,
      method_name     VARCHAR NOT NULL,
      nb              INTEGER NOT NULL,
      maxp            INTEGER DEFAULT 5,
      ar_method       VARCHAR DEFAULT 'burg',
      criterion       VARCHAR DEFAULT 'aic',
      bootadj         BOOLEAN DEFAULT FALSE,
      stat_fn_name    VARCHAR,
      stat_fn_params  JSON,
      garch_dist      VARCHAR,
      method_params   JSON,
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  # ============================================================================
  # Trials table - NEW: tracks execution batches
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS trials (
      trial_id        VARCHAR PRIMARY KEY,
      study_id        VARCHAR NOT NULL REFERENCES studies(study_id),
      dgp_id          VARCHAR NOT NULL REFERENCES dgp_configs(dgp_id),
      method_id       VARCHAR NOT NULL REFERENCES method_configs(method_id),
      trial_name      VARCHAR,
      n_planned       INTEGER,
      n_completed     INTEGER DEFAULT 0,
      started_at      TIMESTAMP DEFAULT current_timestamp,
      completed_at    TIMESTAMP,
      elapsed_seconds DOUBLE,
      status          VARCHAR DEFAULT 'running',
      error_message   TEXT,
      trial_metadata  JSON
    )
  ")

  # ============================================================================
  # Runs table - merged with run_results (1:1 relationship)
  # ============================================================================
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS runs (
      run_id          VARCHAR PRIMARY KEY,
      trial_id        VARCHAR NOT NULL REFERENCES trials(trial_id),
      iteration_num   INTEGER NOT NULL,
      master_seed     BIGINT,
      boot_seeds      BIGINT[],

      -- Observed statistic and bootstrap distribution
      obs_stat        DOUBLE NOT NULL,
      boot_dist       DOUBLE[],

      -- P-values
      pvalue          DOUBLE,
      pvalue_upper    DOUBLE,
      pvalue_lower    DOUBLE,
      pvalue_asymp    DOUBLE,

      -- Null model estimates (fitted to observed data)
      null_ar_phi     DOUBLE[],
      null_vara       DOUBLE,

      -- COBA adjustment fields
      obs_stat_adj    DOUBLE,
      pvalue_adj      DOUBLE,
      adj_factor      DOUBLE,

      -- GARCH-specific
      garch_coef      JSON,

      -- Metadata
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  # Create index for fast iteration lookup
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_runs_trial_iter
    ON runs(trial_id, iteration_num)
  ")

  # ============================================================================
  # Views for easy querying
  # ============================================================================

  # Denormalized view of all run data

  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_run_summary AS
    SELECT
      r.run_id, r.iteration_num,
      t.trial_id, t.trial_name, t.n_planned, t.status as trial_status,
      s.study_id, s.study_name, s.study_type, s.error_type,
      d.dgp_id, d.dgp_name, d.n,
      d.ar_phi, d.ar_phi[1] as phi_value, LENGTH(d.ar_phi) as ar_order,
      d.ma_theta, LENGTH(d.ma_theta) as ma_order,
      d.innov_dist, d.has_trend, d.trend_slope,
      m.method_id, m.method_name, m.nb, m.ar_method, m.bootadj,
      r.obs_stat, r.pvalue, r.pvalue_upper, r.pvalue_lower,
      r.pvalue_asymp, r.pvalue_adj,
      r.null_ar_phi, LENGTH(r.null_ar_phi) as null_ar_order,
      r.master_seed
    FROM runs r
    JOIN trials t ON r.trial_id = t.trial_id
    JOIN studies s ON t.study_id = s.study_id
    JOIN dgp_configs d ON t.dgp_id = d.dgp_id
    JOIN method_configs m ON t.method_id = m.method_id
  ")

  # Aggregated rejection rates (across ALL trials)
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_rejection_rates AS
    SELECT
      s.study_name, s.error_type, d.dgp_name, d.n,
      d.ar_phi, d.ar_phi[1] AS phi_value, m.method_name,
      COUNT(*) AS n_runs,
      COUNT(DISTINCT t.trial_id) AS n_trials,

      -- CO (asymptotic) rejection rates
      AVG(CASE WHEN r.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_co_05,
      SQRT(AVG(CASE WHEN r.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN r.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_co_05_se,

      -- COB (bootstrap) rejection rates
      AVG(CASE WHEN r.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_cob_05,
      SQRT(AVG(CASE WHEN r.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN r.pvalue < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_cob_05_se,

      -- COBA (adjusted) rejection rates
      AVG(CASE WHEN r.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_coba_05,
      SQRT(AVG(CASE WHEN r.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) *
           (1 - AVG(CASE WHEN r.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END)) /
           NULLIF(COUNT(*), 0)) AS reject_coba_05_se

    FROM runs r
    JOIN trials t ON r.trial_id = t.trial_id
    JOIN studies s ON t.study_id = s.study_id
    JOIN dgp_configs d ON t.dgp_id = d.dgp_id
    JOIN method_configs m ON t.method_id = m.method_id
    WHERE t.status = 'completed'
    GROUP BY s.study_name, s.error_type, d.dgp_name, d.n, d.ar_phi, m.method_name
  ")

  # Rejection rates BY TRIAL (for comparing trials)
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_rejection_rates_by_trial AS
    SELECT
      t.trial_id, t.trial_name,
      s.study_name, s.error_type, d.dgp_name, d.n,
      d.ar_phi, d.ar_phi[1] AS phi_value, m.method_name,
      COUNT(*) AS n_runs,

      -- CO (asymptotic) rejection rates
      AVG(CASE WHEN r.pvalue_asymp < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_co_05,

      -- COB (bootstrap) rejection rates
      AVG(CASE WHEN r.pvalue < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_cob_05,

      -- COBA (adjusted) rejection rates
      AVG(CASE WHEN r.pvalue_adj < 0.05 THEN 1.0 ELSE 0.0 END) AS reject_coba_05,

      t.started_at, t.completed_at

    FROM runs r
    JOIN trials t ON r.trial_id = t.trial_id
    JOIN studies s ON t.study_id = s.study_id
    JOIN dgp_configs d ON t.dgp_id = d.dgp_id
    JOIN method_configs m ON t.method_id = m.method_id
    WHERE t.status = 'completed'
    GROUP BY t.trial_id, t.trial_name, s.study_name, s.error_type, d.dgp_name, d.n,
             d.ar_phi, m.method_name, t.started_at, t.completed_at
  ")

  invisible(NULL)
}


#' Migrate Database: Add error_type Column
#'
#' One-time migration to add error_type column to studies table and backfill
#' based on study_name patterns. Safe to run multiple times (idempotent).
#'
#' @param con A DBI connection from [boot_db_connect()].
#'
#' @return Invisible NULL. Prints migration status.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Adds error_type column if it doesn't exist
#'   \item Backfills values based on study_name patterns:
#'     \itemize{
#'       \item "Under Null" -> "iid"
#'       \item "ARCH" (but not "Hetero") -> "garch"
#'       \item "Heteroscedastic" -> "hetero"
#'       \item "t(3)" or "t-dist" -> "t_dist"
#'     }
#'   \item Recreates views to include error_type
#' }
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect("simulations.duckdb")
#' boot_db_migrate_error_type(con)
#' DBI::dbGetQuery(con, "SELECT study_name, error_type FROM studies")
#' DBI::dbDisconnect(con)
#' }
#'
#' @export
boot_db_migrate_error_type <- function(con) {
  # Add column if not exists (DuckDB will error if exists, so we catch it)
  tryCatch({
    DBI::dbExecute(con, "ALTER TABLE studies ADD COLUMN error_type VARCHAR")
    message("Added error_type column to studies table")
  }, error = function(e) {
    if (grepl("already exists", conditionMessage(e), ignore.case = TRUE)) {
      message("error_type column already exists")
    } else {
      stop(e)
    }
  })

  # Backfill based on study_name patterns
  n_updated <- DBI::dbExecute(con, "
    UPDATE studies SET error_type = CASE
      WHEN study_name LIKE '%Under Null%' THEN 'iid'
      WHEN study_name LIKE '%ARCH%' AND study_name NOT LIKE '%Hetero%' THEN 'garch'
      WHEN study_name LIKE '%Heteroscedastic%' THEN 'hetero'
      WHEN study_name LIKE '%t(3)%' OR study_name LIKE '%t-dist%' THEN 't_dist'
      ELSE 'unknown'
    END
    WHERE error_type IS NULL
  ")

  if (n_updated > 0) {
    message(sprintf("Backfilled error_type for %d studies", n_updated))
  } else {
    message("No studies needed backfill (all already have error_type)")
  }

  # Recreate views to include error_type (idempotent)
  boot_db_init(con)
  message("Views recreated with error_type and phi_value columns")

  invisible(NULL)
}


#' Create or Get Study
#'
#' Creates a new study or returns an existing one by name.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param name Character. Name of the study.
#' @param type Character. Type of study: "size", "power", "comparison", or "custom".
#' @param error_type Character. Optional innovation/error type: "iid", "garch",
#'   "hetero", "t_dist", etc. Used for cross-DGP comparisons.
#' @param description Character. Optional description of the study.
#'
#' @return The study_id (character UUID).
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#' study_id <- boot_db_study(con, "AR Size Study", "size",
#'                           error_type = "iid",
#'                           description = "Type I error rates for AR(1) DGPs")
#' }
#'
#' @export
boot_db_study <- function(con, name, type = c("size", "power", "comparison", "custom"),
                          error_type = NULL, description = NULL) {
  type <- match.arg(type)

  # Check if study already exists
  existing <- DBI::dbGetQuery(con,
    "SELECT study_id FROM studies WHERE study_name = ?",
    params = list(name))

  if (nrow(existing) > 0) {
    return(existing$study_id[1])
  }

  # Create new study
  study_id <- .generate_uuid()
  description_val <- if (is.null(description)) NA_character_ else description
  error_type_val <- if (is.null(error_type)) NA_character_ else error_type

  DBI::dbExecute(con,
    "INSERT INTO studies (study_id, study_name, study_type, error_type, description)
     VALUES (?, ?, ?, ?, ?)",
    params = list(study_id, name, type, error_type_val, description_val))

  study_id
}


#' Create or Get DGP Configuration
#'
#' Creates a new DGP configuration or returns an existing one with matching
#' parameters. Uses content-based hashing for deduplication.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param n Integer. Series length.
#' @param ar_phi Numeric vector. AR coefficients. Default NULL (no AR component).
#' @param ma_theta Numeric vector. MA coefficients. Default NULL.
#' @param vara Numeric. Innovation variance. Default 1.0.
#' @param innov_dist Character. Innovation distribution. Default "norm".
#' @param innov_params List. Distribution-specific parameters.
#' @param has_trend Logical. Whether series has trend. Default FALSE.
#' @param trend_slope Numeric. Trend slope. Default 0.
#' @param trend_intercept Numeric. Trend intercept. Default 0.
#' @param dgp_name Character. Optional human-readable name.
#' @param garch_omega,garch_alpha,garch_beta GARCH parameters if applicable.
#'
#' @return The dgp_id (character UUID).
#'
#' @details
#' DGP configurations are deduplicated: if you call this function with the
#' same parameters twice, you get the same dgp_id back. This prevents
#' database clutter when running multiple trials with identical DGP settings.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#' dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7, has_trend = FALSE)
#' }
#'
#' @export
boot_db_dgp <- function(con, n, ar_phi = NULL, ma_theta = NULL, vara = 1.0,
                        innov_dist = "norm", innov_params = NULL,
                        has_trend = FALSE, trend_slope = 0, trend_intercept = 0,
                        dgp_name = NULL, garch_omega = NULL,
                        garch_alpha = NULL, garch_beta = NULL) {

  # Compute config hash for deduplication
  config_hash <- .hash_dgp_config(
    n = n, ar_phi = ar_phi, ma_theta = ma_theta, vara = vara,
    innov_dist = innov_dist, innov_params = innov_params,
    has_trend = has_trend, trend_slope = trend_slope,
    trend_intercept = trend_intercept, garch_omega = garch_omega,
    garch_alpha = garch_alpha, garch_beta = garch_beta
  )

 # Check if config already exists
  existing <- DBI::dbGetQuery(con,
    "SELECT dgp_id FROM dgp_configs WHERE config_hash = ?",
    params = list(config_hash))

  if (nrow(existing) > 0) {
    return(existing$dgp_id[1])
  }

  # Generate name if not provided
  if (is.null(dgp_name)) {
    dgp_name <- .generate_dgp_name(ar_phi, ma_theta, has_trend, trend_slope)
  }

  # Create new DGP config
  dgp_id <- .generate_uuid()

  # Prepare array strings for SQL
 ar_phi_str <- if (!is.null(ar_phi) && length(ar_phi) > 0) {
    .vec_to_array(ar_phi)
  } else {
    "NULL"
  }

  ma_theta_str <- if (!is.null(ma_theta) && length(ma_theta) > 0) {
    .vec_to_array(ma_theta)
  } else {
    "NULL"
  }

  garch_alpha_str <- if (!is.null(garch_alpha)) .vec_to_array(garch_alpha) else "NULL"
  garch_beta_str <- if (!is.null(garch_beta)) .vec_to_array(garch_beta) else "NULL"

  # Convert NULL to NA for scalar params
  innov_params_json <- if (!is.null(innov_params)) {
    as.character(jsonlite::toJSON(innov_params, auto_unbox = TRUE))
  } else {
    NA_character_
  }
  garch_omega_val <- if (is.null(garch_omega)) NA_real_ else garch_omega

  DBI::dbExecute(con, sprintf("
    INSERT INTO dgp_configs (dgp_id, config_hash, dgp_name, n, ar_phi, ma_theta,
                             vara, innov_dist, innov_params, has_trend, trend_slope,
                             trend_intercept, garch_omega, garch_alpha, garch_beta)
    VALUES (?, ?, ?, ?, %s, %s, ?, ?, ?, ?, ?, ?, ?, %s, %s)",
    ar_phi_str, ma_theta_str, garch_alpha_str, garch_beta_str),
    params = list(dgp_id, config_hash, dgp_name, as.integer(n),
                  vara, innov_dist, innov_params_json, has_trend, trend_slope,
                  trend_intercept, garch_omega_val))

  dgp_id
}


#' Create or Get Method Configuration
#'
#' Creates a new method configuration or returns an existing one with matching
#' parameters. Uses content-based hashing for deduplication.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param method_name Character. Name of the bootstrap method.
#' @param nb Integer. Number of bootstrap replicates.
#' @param maxp Integer. Maximum AR order. Default 5.
#' @param ar_method Character. AR estimation method. Default "burg".
#' @param criterion Character. Information criterion. Default "aic".
#' @param bootadj Logical. Whether COBA adjustment is used. Default FALSE.
#' @param stat_fn_name Character. Statistic function name (for wbg_boot_flex).
#' @param stat_fn_params List. Parameters for statistic function.
#' @param garch_dist Character. GARCH distribution (for wbg_boot_garch).
#'
#' @return The method_id (character UUID).
#'
#' @details
#' Method configurations are deduplicated: if you call this function with the
#' same parameters twice, you get the same method_id back.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#' method_id <- boot_db_method(con, "wbg_boot", nb = 399, ar_method = "burg")
#' }
#'
#' @export
boot_db_method <- function(con, method_name, nb, maxp = 5L, ar_method = "burg",
                           criterion = "aic", bootadj = FALSE,
                           stat_fn_name = NULL, stat_fn_params = NULL,
                           garch_dist = NULL) {

  # Compute config hash for deduplication
  config_hash <- .hash_method_config(
    method_name = method_name, nb = nb, maxp = maxp, ar_method = ar_method,
    criterion = criterion, bootadj = bootadj, stat_fn_name = stat_fn_name,
    stat_fn_params = stat_fn_params, garch_dist = garch_dist
  )

  # Check if config already exists
  existing <- DBI::dbGetQuery(con,
    "SELECT method_id FROM method_configs WHERE config_hash = ?",
    params = list(config_hash))

  if (nrow(existing) > 0) {
    return(existing$method_id[1])
  }

  # Create new method config
  method_id <- .generate_uuid()

  # Convert NULL to NA for scalar params
  stat_fn_params_json <- if (!is.null(stat_fn_params)) {
    as.character(jsonlite::toJSON(stat_fn_params, auto_unbox = TRUE))
  } else {
    NA_character_
  }
  stat_fn_name_val <- if (is.null(stat_fn_name)) NA_character_ else stat_fn_name
  garch_dist_val <- if (is.null(garch_dist)) NA_character_ else garch_dist

  DBI::dbExecute(con,
    "INSERT INTO method_configs (method_id, config_hash, method_name, nb, maxp, ar_method,
                                 criterion, bootadj, stat_fn_name, stat_fn_params, garch_dist)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(method_id, config_hash, method_name, as.integer(nb), as.integer(maxp),
                  ar_method, criterion, bootadj, stat_fn_name_val,
                  stat_fn_params_json, garch_dist_val))

  method_id
}


#' Create a New Trial
#'
#' Creates a new trial (execution batch) within a study.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param study_id Character. Study ID from [boot_db_study()].
#' @param dgp_id Character. DGP configuration ID from [boot_db_dgp()].
#' @param method_id Character. Method configuration ID from [boot_db_method()].
#' @param n_planned Integer. Planned number of simulations. Optional.
#' @param trial_name Character. Optional human-readable name for the trial.
#' @param trial_metadata List. Optional metadata to attach to the trial.
#'
#' @return The trial_id (character UUID).
#'
#' @details
#' A trial represents a single execution batch of simulations. For example,
#' if you run 1000 simulations today and 1000 more tomorrow with the same
#' study/DGP/method, those are two separate trials.
#'
#' Unlike studies, DGP configs, and method configs, trials are NEVER
#' deduplicated - each call creates a new trial.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#'
#' study_id <- boot_db_study(con, "AR Size Study", "size")
#' dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
#' method_id <- boot_db_method(con, "wbg_boot", nb = 399)
#'
#' # First trial
#' trial1_id <- boot_db_trial(con, study_id, dgp_id, method_id,
#'                            n_planned = 1000, trial_name = "Trial 1")
#'
#' # Second trial (same config, different execution)
#' trial2_id <- boot_db_trial(con, study_id, dgp_id, method_id,
#'                            n_planned = 1000, trial_name = "Trial 2")
#' }
#'
#' @export
boot_db_trial <- function(con, study_id, dgp_id, method_id,
                          n_planned = NULL, trial_name = NULL,
                          trial_metadata = NULL) {

  trial_id <- .generate_uuid()

  # Convert NULL to NA for scalar params
  n_planned_val <- if (is.null(n_planned)) NA_integer_ else as.integer(n_planned)
  trial_name_val <- if (is.null(trial_name)) NA_character_ else trial_name
  trial_metadata_json <- if (!is.null(trial_metadata)) {
    as.character(jsonlite::toJSON(trial_metadata, auto_unbox = TRUE))
  } else {
    NA_character_
  }

  DBI::dbExecute(con,
    "INSERT INTO trials (trial_id, study_id, dgp_id, method_id, trial_name,
                         n_planned, trial_metadata)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
    params = list(trial_id, study_id, dgp_id, method_id, trial_name_val,
                  n_planned_val, trial_metadata_json))

  trial_id
}


#' Complete a Trial
#'
#' Marks a trial as completed and updates summary statistics.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param trial_id Character. The trial ID to complete.
#' @param status Character. Final status. Default "completed".
#' @param error_message Character. Error message if status is "failed".
#'
#' @return Invisible NULL.
#'
#' @export
boot_db_trial_complete <- function(con, trial_id, status = "completed",
                                   error_message = NULL) {
  # Count completed runs
  n_completed <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) as n FROM runs WHERE trial_id = ?",
    params = list(trial_id))$n[1]

  # Get start time to compute elapsed
  start_time <- DBI::dbGetQuery(con,
    "SELECT started_at FROM trials WHERE trial_id = ?",
    params = list(trial_id))$started_at[1]

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  error_msg_val <- if (is.null(error_message)) NA_character_ else error_message

  DBI::dbExecute(con,
    "UPDATE trials SET n_completed = ?, completed_at = ?, elapsed_seconds = ?,
                       status = ?, error_message = ?
     WHERE trial_id = ?",
    params = list(n_completed, Sys.time(), elapsed, status, error_msg_val, trial_id))

  invisible(NULL)
}


#' Write Bootstrap Result to Database
#'
#' Stores a bootstrap result object in the database.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param result A bootstrap result object from one of the boot functions.
#' @param trial_id Character. Trial ID from [boot_db_trial()].
#' @param iteration_num Integer. Iteration number within the trial (1, 2, 3, ...).
#'
#' @return The run_id (character UUID).
#'
#' @details
#' Automatically detects the type of boot result and extracts appropriate fields.
#' Supports results from: [wbg_boot()], [wbg_boot_fast()], [wbg_boot_flex()],
#' [wbg_boot_garch()], [co_tas_boot()].
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#'
#' study_id <- boot_db_study(con, "Test", "size")
#' dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
#' method_id <- boot_db_method(con, "wbg_boot", nb = 99)
#' trial_id <- boot_db_trial(con, study_id, dgp_id, method_id, n_planned = 10)
#'
#' for (i in 1:10) {
#'   x <- gen_arma(100, phi = 0.7, plot = FALSE)
#'   result <- wbg_boot(x, nb = 99, seed = i * 1000)
#'   boot_db_write(con, result, trial_id, iteration_num = i)
#' }
#'
#' boot_db_trial_complete(con, trial_id)
#' }
#'
#' @export
boot_db_write <- function(con, result, trial_id, iteration_num) {

  run_id <- .generate_uuid()

  # Extract fields based on result type
  extracted <- .extract_boot_fields(result)

  # Prepare array strings for SQL
  boot_seeds_str <- if (!is.null(extracted$boot_seeds)) {
    .vec_to_array(as.numeric(extracted$boot_seeds))
  } else {
    "NULL"
  }

  boot_dist_str <- if (!is.null(extracted$boot_tstats)) {
    .vec_to_array(extracted$boot_tstats)
  } else {
    "NULL"
  }

  null_ar_phi_str <- if (!is.null(extracted$null_ar_phi) && length(extracted$null_ar_phi) > 0) {
    .vec_to_array(extracted$null_ar_phi)
  } else {
    "NULL"
  }

  garch_coef_json <- if (!is.null(extracted$garch_coef)) {
    as.character(jsonlite::toJSON(as.list(extracted$garch_coef), auto_unbox = TRUE))
  } else {
    NA_character_
  }

  # Helper to convert NULL to appropriate NA type
  na_real <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  master_seed_val <- if (is.null(extracted$master_seed)) NA_integer_ else as.integer(extracted$master_seed)

  DBI::dbExecute(con, sprintf("
    INSERT INTO runs (run_id, trial_id, iteration_num, master_seed, boot_seeds,
                      obs_stat, boot_dist, pvalue, pvalue_upper, pvalue_lower,
                      pvalue_asymp, null_ar_phi, null_vara, obs_stat_adj,
                      pvalue_adj, adj_factor, garch_coef)
    VALUES (?, ?, ?, ?, %s, ?, %s, ?, ?, ?, ?, %s, ?, ?, ?, ?, ?)",
    boot_seeds_str, boot_dist_str, null_ar_phi_str),
    params = list(run_id, trial_id, as.integer(iteration_num), master_seed_val,
                  extracted$obs_stat, na_real(extracted$pvalue),
                  na_real(extracted$pvalue_upper), na_real(extracted$pvalue_lower),
                  na_real(extracted$pvalue_asymp), na_real(extracted$null_vara),
                  na_real(extracted$obs_stat_adj), na_real(extracted$pvalue_adj),
                  na_real(extracted$adj_factor), garch_coef_json))

  run_id
}


#' Write Batch of Bootstrap Results to Database
#'
#' Efficiently stores multiple bootstrap results in a single transaction.
#' Much faster than calling [boot_db_write()] in a loop.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param results A list of bootstrap result objects.
#' @param trial_id Character. Trial ID from [boot_db_trial()].
#' @param start_iteration Integer. Starting iteration number. Default 1.
#' @param store_boot_dist Logical. If FALSE, skip storing boot_dist array
#'   for faster writes when only p-values are needed. Default TRUE.
#'
#' @return A character vector of run_ids (one per result).
#'
#' @details
#' This function is optimized for batch writes, which is critical for
#' simulation studies. DuckDB is optimized for columnar batch inserts,
#' making this significantly faster than per-iteration writes.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect(":memory:")
#' boot_db_init(con)
#'
#' study_id <- boot_db_study(con, "Test", "size")
#' dgp_id <- boot_db_dgp(con, n = 100, ar_phi = 0.7)
#' method_id <- boot_db_method(con, "wbg_boot_fast", nb = 99)
#' trial_id <- boot_db_trial(con, study_id, dgp_id, method_id, n_planned = 100)
#'
#' # Run simulations and collect results
#' results <- lapply(1:100, function(i) {
#'   x <- gen_arma(100, phi = 0.7, plot = FALSE)
#'   wbg_boot_fast(x, nb = 99, seed = i)
#' })
#'
#' run_ids <- boot_db_write_batch(con, results, trial_id)
#' boot_db_trial_complete(con, trial_id)
#' }
#'
#' @seealso [boot_db_write()]
#' @export
boot_db_write_batch <- function(con, results, trial_id,
                                start_iteration = 1L, store_boot_dist = TRUE) {
  if (length(results) == 0) {
    return(character(0))
  }

  n_results <- length(results)

  # Generate UUIDs for all runs
  run_ids <- vapply(seq_len(n_results), function(i) .generate_uuid(), character(1))

  # Extract fields from all results
  extracted_list <- lapply(results, .extract_boot_fields)

  DBI::dbBegin(con)

  tryCatch({
    for (i in seq_len(n_results)) {
      e <- extracted_list[[i]]
      iteration_num <- start_iteration + i - 1L

      boot_seeds_str <- if (!is.null(e$boot_seeds)) {
        .vec_to_array(as.numeric(e$boot_seeds))
      } else {
        "NULL"
      }

      boot_dist_str <- if (store_boot_dist && !is.null(e$boot_tstats)) {
        .vec_to_array(e$boot_tstats)
      } else {
        "NULL"
      }

      null_ar_phi_str <- if (!is.null(e$null_ar_phi) && length(e$null_ar_phi) > 0) {
        .vec_to_array(e$null_ar_phi)
      } else {
        "NULL"
      }

      garch_coef_json <- if (!is.null(e$garch_coef)) {
        as.character(jsonlite::toJSON(as.list(e$garch_coef), auto_unbox = TRUE))
      } else {
        NA_character_
      }

      na_real <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
      master_seed_val <- if (is.null(e$master_seed)) NA_integer_ else as.integer(e$master_seed)

      DBI::dbExecute(con, sprintf("
        INSERT INTO runs (run_id, trial_id, iteration_num, master_seed, boot_seeds,
                          obs_stat, boot_dist, pvalue, pvalue_upper, pvalue_lower,
                          pvalue_asymp, null_ar_phi, null_vara, obs_stat_adj,
                          pvalue_adj, adj_factor, garch_coef)
        VALUES (?, ?, ?, ?, %s, ?, %s, ?, ?, ?, ?, %s, ?, ?, ?, ?, ?)",
        boot_seeds_str, boot_dist_str, null_ar_phi_str),
        params = list(run_ids[i], trial_id, as.integer(iteration_num), master_seed_val,
                      e$obs_stat, na_real(e$pvalue),
                      na_real(e$pvalue_upper), na_real(e$pvalue_lower),
                      na_real(e$pvalue_asymp), na_real(e$null_vara),
                      na_real(e$obs_stat_adj), na_real(e$pvalue_adj),
                      na_real(e$adj_factor), garch_coef_json))
    }

    DBI::dbCommit(con)
  }, error = function(err) {
    DBI::dbRollback(con)
    stop("Batch write failed: ", conditionMessage(err), call. = FALSE)
  })

  run_ids
}


#' Write Batch of Bootstrap Results (Lite Version)
#'
#' Ultra-fast batch write that only stores p-values, skipping bootstrap
#' distributions. Use when you only need rejection rates.
#'
#' @inheritParams boot_db_write_batch
#'
#' @return A character vector of run_ids.
#'
#' @details
#' This is a convenience wrapper around [boot_db_write_batch()] with
#' `store_boot_dist = FALSE`. Skipping the boot_dist array (typically
#' 399 doubles per run) dramatically reduces write time and database size.
#'
#' @export
boot_db_write_batch_lite <- function(con, results, trial_id,
                                     start_iteration = 1L) {
  boot_db_write_batch(con, results, trial_id, start_iteration,
                      store_boot_dist = FALSE)
}


#' Query Simulation Runs
#'
#' Retrieves simulation runs with optional filtering.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param study_id Character. Filter by study ID.
#' @param trial_id Character. Filter by trial ID.
#' @param dgp_id Character. Filter by DGP configuration ID.
#' @param method_id Character. Filter by method configuration ID.
#' @param limit Integer. Maximum number of rows to return.
#'
#' @return A data frame with run summary information.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect("simulations.duckdb", read_only = TRUE)
#' runs <- boot_db_query(con)
#' }
#'
#' @export
boot_db_query <- function(con, study_id = NULL, trial_id = NULL,
                          dgp_id = NULL, method_id = NULL, limit = NULL) {

  sql <- "SELECT * FROM v_run_summary WHERE 1=1"
  params <- list()

  if (!is.null(study_id)) {
    sql <- paste(sql, "AND study_id = ?")
    params <- c(params, study_id)
  }
  if (!is.null(trial_id)) {
    sql <- paste(sql, "AND trial_id = ?")
    params <- c(params, trial_id)
  }
  if (!is.null(dgp_id)) {
    sql <- paste(sql, "AND dgp_id = ?")
    params <- c(params, dgp_id)
  }
  if (!is.null(method_id)) {
    sql <- paste(sql, "AND method_id = ?")
    params <- c(params, method_id)
  }
  if (!is.null(limit)) {
    sql <- paste(sql, "LIMIT", as.integer(limit))
  }

  DBI::dbGetQuery(con, sql, params = params)
}


#' Get Bootstrap Replicates for a Run
#'
#' Retrieves the full bootstrap distribution for a specific run.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param run_id Character. The run ID.
#'
#' @return A list with `obs_stat` and `boot_dist` (numeric vector).
#'
#' @export
boot_db_get_replicates <- function(con, run_id) {
  result <- DBI::dbGetQuery(con,
    "SELECT obs_stat, boot_dist FROM runs WHERE run_id = ?",
    params = list(run_id))

  if (nrow(result) == 0) {
    stop("Run not found: ", run_id, call. = FALSE)
  }

  list(
    obs_stat = result$obs_stat[1],
    boot_dist = result$boot_dist[[1]]
  )
}


#' Recalculate P-value from Stored Bootstrap Distribution
#'
#' Recomputes the p-value using the stored bootstrap statistics.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param run_id Character. The run ID.
#' @param two_sided Logical. If TRUE (default), compute two-sided p-value.
#'
#' @return A list with `pvalue`, `pvalue_upper`, `pvalue_lower`, and `nb`.
#'
#' @examples
#' \dontrun{
#' con <- boot_db_connect("simulations.duckdb", read_only = TRUE)
#' recalc <- boot_db_recalc_pvalue(con, "some-run-id")
#' }
#'
#' @export
boot_db_recalc_pvalue <- function(con, run_id, two_sided = TRUE) {
  data <- boot_db_get_replicates(con, run_id)

  obs_stat <- data$obs_stat
  boot_dist <- data$boot_dist
  nb <- length(boot_dist)

  # Compute p-values with plus-one correction
  pvalue_two <- (sum(abs(boot_dist) >= abs(obs_stat)) + 1) / (nb + 1)
  pvalue_upper <- (sum(boot_dist >= obs_stat) + 1) / (nb + 1)
  pvalue_lower <- (sum(boot_dist <= obs_stat) + 1) / (nb + 1)

  list(
    pvalue = pvalue_two,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    nb = nb
  )
}


#' Get Rejection Rates
#'
#' Retrieves aggregated rejection rates from the database.
#'
#' @param con A DBI connection from [boot_db_connect()].
#' @param study_name Character. Optional filter by study name.
#' @param by_trial Logical. If TRUE, return rates per trial. Default FALSE.
#'
#' @return A data frame with rejection rates by DGP and method.
#'
#' @export
boot_db_rejection_rates <- function(con, study_name = NULL, by_trial = FALSE) {
  view_name <- if (by_trial) "v_rejection_rates_by_trial" else "v_rejection_rates"

  if (is.null(study_name)) {
    DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", view_name))
  } else {
    DBI::dbGetQuery(con,
      sprintf("SELECT * FROM %s WHERE study_name = ?", view_name),
      params = list(study_name))
  }
}


# ============================================================================
# Internal Helper Functions
# ============================================================================

#' Generate UUID (RNG-independent)
#'
#' Uses system time, process info, and counter to generate unique IDs
#' without affecting or being affected by the global RNG state.
#' @noRd
.generate_uuid <- function() {
  # High-resolution timestamp (microseconds)
  time_us <- as.numeric(Sys.time()) * 1e6
  time_hex <- sprintf("%08x", as.integer(time_us %% 2^31))

  # Process ID
  pid_hex <- sprintf("%04x", Sys.getpid() %% 65536)

  # Monotonic counter - guarantees uniqueness even within same microsecond
  counter <- get0(".uuid_counter", envir = .tstse_env, ifnotfound = 0L) + 1L
  assign(".uuid_counter", counter, envir = .tstse_env)
  counter_hex <- sprintf("%04x", counter %% 65536)

  # Additional entropy from proc.time elapsed (changes every call)
  pt <- proc.time()
  pt_hex <- sprintf("%04x", as.integer((pt[3] * 1e6) %% 65536))

  # Low bits of timestamp for more variation
  low_hex <- sprintf("%08x", as.integer((time_us * 1000) %% 2^31))

  paste(time_hex, pid_hex, counter_hex, pt_hex, low_hex, sep = "-")
}

# Environment for UUID counter
.tstse_env <- new.env(parent = emptyenv())


#' Convert R vector to DuckDB array literal
#' @noRd
.vec_to_array <- function(x) {
  if (is.null(x) || length(x) == 0) return("NULL")
  paste0("[", paste(x, collapse = ", "), "]")
}


#' Generate DGP name from parameters
#' @noRd
.generate_dgp_name <- function(ar_phi, ma_theta, has_trend, trend_slope) {
  parts <- character(0)

  if (!is.null(ar_phi) && length(ar_phi) > 0) {
    parts <- c(parts, sprintf("AR(%d)", length(ar_phi)))
  }
  if (!is.null(ma_theta) && length(ma_theta) > 0) {
    parts <- c(parts, sprintf("MA(%d)", length(ma_theta)))
  }
  if (length(parts) == 0) {
    parts <- "WN"
  }

  if (has_trend && trend_slope != 0) {
    parts <- c(parts, sprintf("trend=%.3f", trend_slope))
  }

  paste(parts, collapse = " + ")
}


#' Compute hash for DGP config deduplication
#' @noRd
.hash_dgp_config <- function(n, ar_phi, ma_theta, vara, innov_dist, innov_params,
                             has_trend, trend_slope, trend_intercept,
                             garch_omega, garch_alpha, garch_beta) {
  # Create canonical representation
  config_list <- list(
    n = as.integer(n),
    ar_phi = if (is.null(ar_phi)) NULL else round(ar_phi, 10),
    ma_theta = if (is.null(ma_theta)) NULL else round(ma_theta, 10),
    vara = round(vara, 10),
    innov_dist = innov_dist,
    innov_params = innov_params,
    has_trend = has_trend,
    trend_slope = round(trend_slope, 10),
    trend_intercept = round(trend_intercept, 10),
    garch_omega = if (is.null(garch_omega)) NULL else round(garch_omega, 10),
    garch_alpha = if (is.null(garch_alpha)) NULL else round(garch_alpha, 10),
    garch_beta = if (is.null(garch_beta)) NULL else round(garch_beta, 10)
  )

  # Convert to JSON for consistent serialization
  json_str <- as.character(jsonlite::toJSON(config_list, auto_unbox = TRUE, digits = 10))

  # Use digest for hashing (MD5 is fast and sufficient for deduplication)
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(json_str, algo = "md5")
  } else {
    # Fallback: use R's built-in hashing
    sprintf("%x", abs(.Internal(nchar(json_str, "bytes", FALSE)) * 31 +
                      sum(utf8ToInt(substr(json_str, 1, min(100, nchar(json_str)))))))
  }
}


#' Compute hash for method config deduplication
#' @noRd
.hash_method_config <- function(method_name, nb, maxp, ar_method, criterion,
                                bootadj, stat_fn_name, stat_fn_params, garch_dist) {
  # Create canonical representation
  config_list <- list(
    method_name = method_name,
    nb = as.integer(nb),
    maxp = as.integer(maxp),
    ar_method = ar_method,
    criterion = criterion,
    bootadj = bootadj,
    stat_fn_name = stat_fn_name,
    stat_fn_params = stat_fn_params,
    garch_dist = garch_dist
  )

  # Convert to JSON for consistent serialization
  json_str <- as.character(jsonlite::toJSON(config_list, auto_unbox = TRUE))

  # Use digest for hashing
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(json_str, algo = "md5")
  } else {
    # Fallback
    sprintf("%x", abs(.Internal(nchar(json_str, "bytes", FALSE)) * 31 +
                      sum(utf8ToInt(substr(json_str, 1, min(100, nchar(json_str)))))))
  }
}


#' Extract fields from boot result object
#' @noRd
.extract_boot_fields <- function(result) {
  # Handle different result types
  # All should have these after our Phase 0 updates

  # Observed statistic (different names in different functions)
  obs_stat <- result$obs_stat %||% result$tco_obs %||% result$tco

  # Bootstrap distribution
  boot_tstats <- result$boot_dist %||% result$boot_tstats %||% result$boot_tco

  # P-values
  pvalue <- result$pvalue
  pvalue_upper <- result$pvalue_upper
  pvalue_lower <- result$pvalue_lower
  pvalue_asymp <- result$pvalue_asymp

  # Null model parameters
  null_ar_phi <- result$ar_phi %||% result$phi %||% result$ar_coef %||% result$phi_null

  null_vara <- result$ar_vara %||% result$vara

  # Reproducibility
  master_seed <- result$master_seed
  boot_seeds <- result$boot_seeds

  # COBA adjustment
  obs_stat_adj <- result$obs_stat_adj %||% result$tco_obs_adj
  pvalue_adj <- result$pvalue_adj
  adj_factor <- result$adj_factor

  # GARCH (only for wbg_boot_garch)
  garch_coef <- result$garch_coef

  list(
    obs_stat = obs_stat,
    boot_tstats = boot_tstats,
    pvalue = pvalue,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    pvalue_asymp = pvalue_asymp,
    null_ar_phi = null_ar_phi,
    null_vara = null_vara,
    master_seed = master_seed,
    boot_seeds = boot_seeds,
    obs_stat_adj = obs_stat_adj,
    pvalue_adj = pvalue_adj,
    adj_factor = adj_factor,
    garch_coef = garch_coef
  )
}


#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


# ============================================================================
# High-Level Wrapper: Closure/Factory Pattern
# ============================================================================

#' Create a Boot Study Session
#'
#' Creates a study session that returns bound functions for saving results
#' and querying the database. Uses a closure/factory pattern for clean,
#' self-contained state management.
#'
#' @param path Character. Path to the DuckDB database file.
#'   Use ":memory:" for an in-memory database (for testing).
#' @param study_name Character. Name of the study.
#' @param dgp List. DGP parameters. Must include `n` (series length).
#'   Other common parameters: `ar_phi`, `has_trend`, `trend_slope`.
#' @param method List. Method parameters. Must include `name` (method name)
#'   and `nb` (number of bootstrap replicates). Other parameters passed
#'   to [boot_db_method()].
#' @param type Character. Study type: "size", "power", "comparison", or "custom".
#' @param description Character. Optional study description.
#' @param n_planned Integer. Planned number of simulations for this trial.
#' @param trial_name Character. Optional name for this trial.
#'
#' @return A `boot_study` object (list) with methods:
#' \describe{
#'   \item{`$save(result)`}{Save a bootstrap result (auto-increments iteration)}
#'   \item{`$save_batch(results)`}{Save multiple results efficiently}
#'   \item{`$query(...)`}{Query runs for this trial}
#'   \item{`$query_study(...)`}{Query all runs in the study (across trials)}
#'   \item{`$complete()`}{Mark the trial as completed}
#'   \item{`$end()`}{Close the database connection}
#'   \item{`$info`}{List with path, study_name, study_id, dgp_id, method_id, trial_id}
#' }
#'
#' @details
#' Each call to `boot_study()` creates a NEW TRIAL. If you want to add
#' results to an existing trial, use the lower-level functions directly.
#'
#' The returned object uses closures to capture the database connection
#' and configuration IDs. This avoids global state issues and makes each
#' study session self-contained.
#'
#' Always call `$end()` when done to close the database connection.
#' Use `on.exit(study$end())` in scripts for automatic cleanup on error.
#'
#' @examples
#' \dontrun{
#' # Create a study session (creates a new trial)
#' study <- boot_study("my_study.duckdb", "AR1 Size Study",
#'                     dgp = list(n = 100, ar_phi = 0.7),
#'                     method = list(name = "wbg_boot", nb = 399),
#'                     type = "size",
#'                     n_planned = 1000,
#'                     trial_name = "Trial 1")
#' on.exit(study$end())
#'
#' # Run simulations
#' for (sim in 1:1000) {
#'   x <- gen_arma(100, phi = 0.7, plot = FALSE)
#'   result <- wbg_boot(x, nb = 399, seed = sim * 1000)
#'   study$save(result)
#' }
#'
#' # Mark trial complete and get results
#' study$complete()
#' runs <- study$query()
#'
#' # Run a second trial with same config
#' study2 <- boot_study("my_study.duckdb", "AR1 Size Study",
#'                      dgp = list(n = 100, ar_phi = 0.7),
#'                      method = list(name = "wbg_boot", nb = 399),
#'                      type = "size",
#'                      n_planned = 1000,
#'                      trial_name = "Trial 2")
#' # ... run more simulations ...
#' study2$complete()
#' study2$end()
#'
#' # Query aggregated results across all trials
#' con <- boot_db_connect("my_study.duckdb", read_only = TRUE)
#' boot_db_rejection_rates(con, "AR1 Size Study")
#' boot_db_rejection_rates(con, "AR1 Size Study", by_trial = TRUE)
#' DBI::dbDisconnect(con)
#' }
#'
#' @seealso [boot_db_connect()], [boot_db_write()], [boot_db_query()], [boot_db_trial()]
#' @export
boot_study <- function(path, study_name, dgp, method,
                       type = c("size", "power", "comparison", "custom"),
                       description = NULL, n_planned = NULL, trial_name = NULL) {
  type <- match.arg(type)

  # Validate inputs
  if (!"n" %in% names(dgp)) {
    stop("dgp must include 'n' (series length)", call. = FALSE)
  }
  if (!"name" %in% names(method)) {
    stop("method must include 'name' (method name)", call. = FALSE)
  }
  if (!"nb" %in% names(method)) {
    stop("method must include 'nb' (number of bootstrap replicates)", call. = FALSE)
  }

  # Create internal context (captured by closures)
  con <- boot_db_connect(path)
  boot_db_init(con)

  study_id <- boot_db_study(con, study_name, type, description)

  dgp_args <- c(list(con = con), dgp)
  dgp_id <- do.call(boot_db_dgp, dgp_args)

  method_name <- method$name
  method$name <- NULL
  method_args <- c(list(con = con, method_name = method_name), method)
  method_id <- do.call(boot_db_method, method_args)

  # Create a NEW trial for this session
  trial_id <- boot_db_trial(con, study_id, dgp_id, method_id,
                            n_planned = n_planned, trial_name = trial_name)

  # Track iteration number
  current_iteration <- 0L

  # Return list of bound functions
  structure(
    list(
      save = function(result) {
        current_iteration <<- current_iteration + 1L
        run_id <- boot_db_write(con, result, trial_id, iteration_num = current_iteration)
        invisible(run_id)
      },

      save_batch = function(results) {
        start_iter <- current_iteration + 1L
        run_ids <- boot_db_write_batch(con, results, trial_id, start_iteration = start_iter)
        current_iteration <<- current_iteration + length(results)
        invisible(run_ids)
      },

      query = function(...) {
        boot_db_query(con, trial_id = trial_id, ...)
      },

      query_study = function(...) {
        boot_db_query(con, study_id = study_id, ...)
      },

      complete = function(status = "completed") {
        boot_db_trial_complete(con, trial_id, status = status)
        invisible(NULL)
      },

      end = function() {
        DBI::dbDisconnect(con)
        invisible(NULL)
      },

      info = list(
        path = path,
        study_name = study_name,
        study_id = study_id,
        dgp_id = dgp_id,
        method_id = method_id,
        trial_id = trial_id
      )
    ),
    class = "boot_study"
  )
}


#' @export
print.boot_study <- function(x, ...) {
  cat("<boot_study>\n")
  cat("  Database:", x$info$path, "\n")
  cat("  Study:", x$info$study_name, "\n")
  cat("  Study ID:", substr(x$info$study_id, 1, 8), "...\n")
  cat("  DGP ID:", substr(x$info$dgp_id, 1, 8), "...\n")
  cat("  Method ID:", substr(x$info$method_id, 1, 8), "...\n")
  cat("  Trial ID:", substr(x$info$trial_id, 1, 8), "...\n")
  cat("Usage: $save(result), $save_batch(results), $query(), $complete(), $end()\n")
  invisible(x)
}
