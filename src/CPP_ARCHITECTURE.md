# C++ Architecture Documentation

This document maps all C++ files, functions, and their usage in the tstse package.

## Overview

The C++ code is organized into three layers:
1. **API Layer** (`api_*.cpp`) - R-facing Rcpp interfaces
2. **Kernel Layer** (`kernel_*.cpp`) - Thread-safe pure C++ implementations (hot path)
3. **Utility Layer** (`util_*.cpp`) - Helper functions for various operations

---

## File Summary Table

| File | Category | Thread-Safe | Purpose |
|------|----------|-------------|---------|
| `kernel_types.h` | Types | Yes | Shared struct definitions |
| `kernel_burg_aic.cpp` | Hot Path | Yes | Burg algorithm with AIC selection |
| `kernel_co_tstat.cpp` | Hot Path | Yes | Cochrane-Orcutt t-statistic |
| `kernel_pw_tstat.cpp` | Hot Path | Yes | Prais-Winsten t-statistic |
| `kernel_gen_ar.cpp` | Hot Path | Yes | Thread-safe AR generation (dqrng) |
| `kernel_wbg_boot.cpp` | Orchestration | Yes | TBB parallel bootstrap |
| `api_burg_aic.cpp` | API | No | Burg AIC R interface |
| `api_co_tstat.cpp` | API | No | Cochrane-Orcutt R interface |
| `api_pw_tstat.cpp` | API | No | Prais-Winsten R interface |
| `api_ols_detrend.cpp` | API | No | OLS detrending R interface |
| `api_ar_transform.cpp` | API | Yes | AR filter R interface |
| `api_burg_fit.cpp` | API | No | Fixed-order Burg (deprecated) |
| `api_backcast.cpp` | API | No | ARMA residual initialization |
| `util_gen_arch.cpp` | Utility | No | ARCH generation |
| `util_gen_garch.cpp` | Utility | No | GARCH generation |
| `util_gegenb.cpp` | Utility | Yes | Gegenbauer coefficients |
| `util_convolve.cpp` | Utility | Yes | Truncated convolution |
| `util_sen_slope.cpp` | Utility | Yes | Sen's slope estimator |
| `RcppExports.cpp` | Auto-gen | N/A | DO NOT EDIT |

---

## Detailed Function Reference

### kernel_types.h
Shared data structures used across modules:
- `BurgResult` - AR fitting output (phi, vara, p, ic)
- `COResult` - Cochrane-Orcutt result (tco, p, phi, vara)
- `PWResult` - Prais-Winsten result (tpw, rho, vara)
- `CoBootstrapWorkspace` - Pre-allocated buffers for CO bootstrap
- `PwBootstrapWorkspace` - Pre-allocated buffers for PW bootstrap
- `CriterionType` enum - AIC, AICc, BIC

---

### kernel_burg_aic.cpp
Pure C++ Burg's algorithm with AIC/BIC selection.

**Functions:**
| Function | Description | Called By |
|----------|-------------|-----------|
| `burg_aic_select_pure()` | Main Burg with IC selection | api_burg_aic.cpp |
| `burg_fit_pure()` | Fixed-order Burg fitting | Internal |
| `burg_aic_select_ws()` | Workspace-aware (zero alloc) | kernel_wbg_boot.cpp |

**R Interface:** `burg_aic_select_cpp()` in api_burg_aic.cpp
**Used by R:** `co_fast()` (R/co_fast.R:81)

---

### kernel_co_tstat.cpp
Cochrane-Orcutt t-statistic with O(1) time transform optimization.

**Functions:**
| Function | Description | Called By |
|----------|-------------|-----------|
| `co_tstat_fused()` | Single-pass CO (O(1) time transform) | co_tstat_ws |
| `co_tstat_pure()` | Entry point | api_co_tstat.cpp |
| `co_tstat_ws()` | Workspace-aware | kernel_wbg_boot.cpp |
| `co_full_pure()` | Full results with AR coefficients | api_co_tstat.cpp |
| `ols_detrend_ws()` | Workspace detrending | Internal |

**Key Optimization:** `co_tstat_fused()` computes time transform in O(1) using: t* = phi(1)*t + sum(j*phi_j)

**R Interface:** `co_tstat_cpp()`, `co_full_cpp()`
**Used by R:** `co_fast()` (R/co_fast.R)

---

### kernel_pw_tstat.cpp
Prais-Winsten t-statistic for AR(1) correlated errors.

**Functions:**
| Function | Description | Called By |
|----------|-------------|-----------|
| `estimate_rho()` | AR(1) from lag-1 autocorrelation | Internal |
| `pw_tstat_fused()` | Single-pass PW transform | pw_tstat_pure |
| `pw_tstat_pure()` | Two-step PW | api_pw_tstat.cpp |
| `pw_tstat_ws()` | Workspace-aware | Future PW bootstrap |
| `pw_full_pure()` | Full results | api_pw_tstat.cpp |
| `pw_full_iterative_pure()` | Iterative rho estimation | api_pw_tstat.cpp |

**R Interface:** `pw_tstat_cpp()`, `pw_full_cpp()`, `pw_tstat_iterative_cpp()`, `pw_full_iterative_cpp()`
**Used by R:** `pw_fast()` (R/pw_fast.R:95-97)

---

### kernel_gen_ar.cpp
Thread-safe AR generation using dqrng xoshiro256+.

**Functions:**
| Function | Description | Called By |
|----------|-------------|-----------|
| `calc_ar_burnin_cpp()` | Adaptive burn-in calculation | kernel_wbg_boot.cpp |
| `gen_ar_dqrng_boxmuller()` | Box-Muller AR generation | Internal |
| `gen_ar_into_precomputed()` | Zero-copy with precomputed constants | kernel_wbg_boot.cpp |

**R Interface:** `gen_ar_cpp()`, `gen_ar_seeded_cpp()`
**Used by C++:** `WBGBootstrapWorker::operator()` in kernel_wbg_boot.cpp

---

### kernel_wbg_boot.cpp
TBB parallel bootstrap orchestration.

**Worker Classes:**
| Class | Description |
|-------|-------------|
| `WBGBootstrapWorker` | Standard CO bootstrap |
| `WBGBootstrapCOBAWorker` | Bootstrap adjustment variant (with AR fitting for variance correction) |

**Exported Functions:**
| Function | Description |
|----------|-------------|
| `wbg_bootstrap_kernel_cpp()` | Main parallel bootstrap |
| `wbg_bootstrap_kernel_grain_cpp()` | With grain size control |
| `wbg_bootstrap_coba_kernel_cpp()` | Bootstrap adjustment (two-stage) |
| `wbg_bootstrap_copw_kernel_cpp()` | CO + PW combined |

**Dependencies:** kernel_gen_ar.cpp, kernel_co_tstat.cpp, kernel_burg_aic.cpp
**Used by R:** `wbg_boot_fast()` (R/wbg_boot_fast.R)

---

### api_burg_aic.cpp
**Exported:** `burg_aic_select_cpp(x, maxp, criterion, min_p)` -> List{p, phi, vara, aic}
**Used by R:** `co_fast()` (R/co_fast.R:81)

---

### api_co_tstat.cpp
**Exported:**
- `co_tstat_cpp(x, maxp, criterion)` -> double
- `co_full_cpp(x, maxp, criterion)` -> List

**Dependencies:** ols_detrend_cpp, burg_aic_select_cpp, ar_transform_cpp

---

### api_pw_tstat.cpp
**Exported:**
- `pw_tstat_cpp(x)` -> double
- `pw_full_cpp(x)` -> List{tpw, rho, vara}
- `pw_tstat_iterative_cpp(x, max_iter, tol)` -> double
- `pw_full_iterative_cpp(x, max_iter, tol)` -> List

**Used by R:** `pw_fast()` (R/pw_fast.R)

---

### api_ols_detrend.cpp
**Exported:**
- `ols_detrend_cpp(x)` -> arma::vec (residuals)
- `ols_detrend_full_cpp(x)` -> List{residuals, intercept, slope}

**Used by R:** `co_fast()` (R/co_fast.R:77)

---

### api_ar_transform.cpp
**Exported:**
- `ar_transform_cpp(x, phi)` -> arma::vec (filtered)
- `co_time_transform_cpp(n, phi)` -> arma::vec (DEPRECATED - hot path uses fused O(1))

**Used by R:** `co_fast()` (R/co_fast.R:88)

---

### api_backcast.cpp
**Exported:** `backcast_cpp(x, phi, theta, n_back)` -> NumericVector
**Used by R:** `backcast()` (R/backcast.R)

---

### util_gen_arch.cpp
**Exported:** `gen_arch_cpp(n, alpha0, alpha, spin, seed)` -> NumericVector
**Note:** Uses R's RNG (not thread-safe)
**Used by R:** `gen_arch()` (R/gen_arch.R)

---

### util_gen_garch.cpp
**Exported:** `gen_garch_cpp(n, alpha0, alpha, beta, spin, seed)` -> NumericVector
**Note:** Uses R's RNG (not thread-safe)
**Used by R:** `gen_garch()` (R/gen_garch.R)

---

### util_gegenb.cpp
**Exported:** `gegenb_cpp(u, d, n_coef)` -> NumericVector
**Purpose:** Gegenbauer polynomial coefficients for long-memory ARUMA
**Used by R:** `gegenb()` (R/gegenb.R)

---

### util_convolve.cpp
**Exported:** `convolve_truncated_cpp(C1, C2, n)` -> NumericVector
**Purpose:** Polynomial multiplication (first n coefficients)
**Used by R:** `convolve_truncated()` (R/macoef_geg.R)

---

### util_sen_slope.cpp
**Exported:** `sen_slope_cpp(x)` -> double
**Purpose:** Median of pairwise slopes (robust trend)
**Used by R:** `sen_slope()` (R/stat_generators.R)

---

## Call Hierarchy

### WBG Bootstrap (Parallel Hot Path)
```
R: wbg_boot_fast()
  |
  v
C++: wbg_bootstrap_kernel_cpp()           [kernel_wbg_boot.cpp]
  |
  +-> WBGBootstrapWorker::operator()      [TBB thread pool]
       |
       +-> gen_ar_into_precomputed()      [kernel_gen_ar.cpp]
       |     Uses: dqrng xoshiro256+ (thread-safe)
       |
       +-> co_tstat_ws()                  [kernel_co_tstat.cpp]
             |
             +-> ols_detrend_ws()
             +-> burg_aic_select_ws()     [kernel_burg_aic.cpp]
             +-> co_tstat_fused()         [O(1) time transform]
```

### Cochrane-Orcutt (Sequential)
```
R: co_fast()
  |
  +-> ols_detrend_cpp()        [api_ols_detrend.cpp]
  +-> burg_aic_select_cpp()    [api_burg_aic.cpp]
  +-> ar_transform_cpp()       [api_ar_transform.cpp]
```

### Prais-Winsten (Sequential)
```
R: pw_fast()
  |
  v
C++: pw_full_cpp() or pw_full_iterative_cpp()   [api_pw_tstat.cpp]
  |
  +-> pw_tstat_pure() or pw_full_iterative_pure()  [kernel_pw_tstat.cpp]
       |
       +-> pw_ols_detrend_ws()
       +-> estimate_rho()
       +-> pw_tstat_fused()
```

---

## Key Optimizations

| Component | Optimization | Impact |
|-----------|--------------|--------|
| AR Generation | dqrng xoshiro256+ | ~2-3x faster than MT |
| Burg Algorithm | Single-pass Levinson | O(n*p) vs O(n*p^2) |
| CO Time Transform | O(1) algebraic identity | Eliminated O(p) per obs |
| OLS Regression | SS_res = SS_y - b^2*SS_t | No residual allocation |
| Workspace Reuse | Thread-local buffers | Zero allocs in hot loop |
| Criterion Enum | String->enum once | No string comparison |

---

## Auto-Generated Files (DO NOT EDIT)

- `RcppExports.cpp` - C++ side Rcpp wrappers
- `R/RcppExports.R` - R side .Call() stubs

Regenerate with: `Rcpp::compileAttributes()`
