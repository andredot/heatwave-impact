# =============================================================================
# Configuration of the validation pipeline.
# Every analytical choice lives here, so the report can print it and targets
# re-runs whatever depends on a changed value.
# =============================================================================

cfg <- list(

  # ---- scope ---------------------------------------------------------------
  country = "IT",

  # ---- inputs --------------------------------------------------------------
  istat_csv  = "data/raw/comuni_giornaliero_30giugno26.csv",
  data_dir   = "data",
  zenodo_dir = "data/zenodo",
  cache_dir  = "data/cache",
  # 491 MB file with 1000 coefficient draws that keep the correlation between
  # age groups. FALSE = draw from vcov.csv (age groups independent).
  use_coef_simu = TRUE,

  # ---- city definition (Urban Audit city -> Istat comuni) -------------------
  urau_year     = 2021,   # release used to label the crosswalk
  urau_years    = c(2021, 2020),  # candidate releases: 2021 merges greater cities into
                                  # CITIES, 2020 still has the smaller core cities
  urau_levels   = c("CITIES", "GREATER_CITIES"),  # candidates; levels absent from a release are skipped
                                                  # (from 2021 greater cities are inside CITIES)
  lau_year      = 2024,   # GISCO LAU release (closest to current Istat codes)
  lau_share_min = 0.5,    # a comune belongs to a city if >= 50% of its area is inside
  point_max_dist = 5000,  # metres: nearest polygon accepted when a point falls just outside
  missing_pop_tol = 0.02, # share of a city's population allowed to have no Istat code
  crosswalk_overrides = "data/raw/crosswalk_overrides.csv",  # optional manual fixes
  pop_ratio_ok    = c(0.85, 1.15),  # LAU population / Masselot population
  deaths_ratio_ok = c(0.80, 1.25),  # Istat deaths 20+ / deaths implied by metadata
  keep_flagged    = TRUE,           # analyse cities flagged only on population/name (never on deaths)

  # ---- periods -------------------------------------------------------------
  test_start = as.Date("2023-01-01"),
  data_end   = NULL,                 # NULL = last day with Istat data
  run_reproduction = TRUE,           # also test 2011-2019 (independent of MCC Italy 2001-2010)
  repro_start = as.Date("2011-01-01"),
  repro_end   = as.Date("2019-12-31"),
  overlap_start = as.Date("2015-01-01"),  # years where our ERA5-Land series is compared
  overlap_end   = as.Date("2019-12-31"),  # with Masselot's own era5series.csv

  # ---- exposure ------------------------------------------------------------
  era5_timezone = "GMT",
  era5_mode     = "polygon",  # "polygon" averages the ERA5-Land cells inside the city
                              # (as Masselot did); "point" uses the city point only
  era5_max_points = 8,        # cells per city in polygon mode (Open-Meteo quota)
  bias_correct  = TRUE,   # remove city x month mean difference vs Masselot's series

  # ---- first-stage model, as in Masselot et al. (05_FirstStage.R) -----------
  lag_max = 21, lag_nknots = 3,
  var_pct = c(10, 75, 90), var_degree = 2,
  df_time_per_year = 7,
  min_deaths = 5000,        # cumulative age binning threshold
  mmp_range  = c(25, 99),   # percentiles where the MMT is searched
  rr_pct     = 99,          # percentile for the heat relative risk
  nsim       = 1000,

  # ---- heatwave episodes ---------------------------------------------------
  hw_pct      = 95,     # city-specific percentile of the 1990-2019 distribution
  hw_min_days = 3,
  hw_max_gap  = 1,      # hot spells separated by <= 1 day are merged
  warm_months = 5:9,
  tail_days   = 7,      # days after the episode counted in the episode window
  post_days   = 14,     # further days used to look at mortality displacement

  # ---- successive time windows (test of a trend in vulnerability) ----------
  period_windows = list(c("2011-01-01", "2015-12-31"), c("2016-01-01", "2019-12-31"),
                        c("2023-01-01", NA)),   # NA = last day with data

  # ---- recalibrated Italian curves -----------------------------------------
  recal_windows = list(c("2011-01-01", "2019-12-31"), c("2016-01-01", "2019-12-31")),
  example_city  = "Milan",
  example_years = 2024:2026,

  # ---- counterfactual baseline: deaths expected without heat --------------
  baseline_exclude_months = c(12, 1, 2, 3),  # influenza season, not used to fit the baseline
  baseline_hot_lag   = 3,     # days after a day above the MMT also excluded
  baseline_harmonics = 1,     # annual sine/cosine pairs (1, as EuroMOMO: summer has few reference days)
  baseline_trend_df_per_year = 0.5,
  baseline_min_daily = 2,     # age groups with fewer deaths/day use the 20+ model x their share
  baseline_cold_offset = TRUE, # remove the published cold effect from reference days (B = deaths at MMT)

  # ---- temperature forecasts (Open-Meteo Previous Runs API) ------------------
  forecast_model    = "ecmwf_ifs025",
  forecast_start    = as.Date("2024-01-01"),
  forecast_leads    = 1:7,
  forecast_timezone = "GMT",
  fc_debias         = TRUE,  # remove city x month x lead mean error (non-episode days)

  seed = 20260922
)
