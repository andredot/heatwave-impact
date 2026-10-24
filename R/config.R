# =============================================================================
# Configuration of the validation pipeline.
# Every analytical choice lives here, so the report can print it and targets
# re-runs whatever depends on a changed value.
# =============================================================================

cfg <- list(

  # ---- scope ---------------------------------------------------------------
  country = "IT",

  # ---- inputs --------------------------------------------------------------
  istat_csv  = "data/raw/comuni_giornaliero_30giugno26.csv",  # used if already present
  istat_url  = paste0("https://www.istat.it/storage/dati_mortalita/giugno-2026/",
                      "decessi-comunali-mese-provvisori-3.zip"),
  istat_dir  = "data/raw",
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
  # city composition:
  #   "istat"    official Istat composition of the Cities 2021 (SITUAS report
  #              "Composizione delle Citta (Cities) 2021"), shipped with the
  #              repository: Urban Audit code -> Istat comune codes, no
  #              inference. See R/istat_cities.R to refresh it.
  #   "official" Eurostat CITY-LAU correspondence table
  #              (https://ec.europa.eu/eurostat/web/nuts/local-administrative-units,
  #              files "CITY - LAU 2021"); paste its URL in city_table_url.
  #   "geometry" the original matching (city point in GISCO polygons).
  # Whatever is chosen, the geometric matching is still computed and the
  # checks report comune by comune how far the two agree.
  city_source = "istat",
  istat_cities_file = "reference/istat_cities_2021.csv",
  city_table_url = "",
  lau_share_min = 0.5,    # a comune belongs to a city if >= 50% of its area is inside
  point_max_dist = 5000,  # metres: nearest polygon accepted when a point falls just outside
  missing_pop_tol = 0.02, # share of a city's population allowed to have no Istat code
  crosswalk_overrides = "data/raw/crosswalk_overrides.csv",  # optional manual fixes
  pop_ratio_ok    = c(0.85, 1.15),  # LAU population / Masselot population
  deaths_ratio_ok = c(0.80, 1.25),  # Istat deaths 20+ / deaths implied by metadata
  keep_flagged    = TRUE,           # analyse cities flagged only on population/name (never on deaths)
  # region shown city by city, including cities below Masselot's death
  # threshold (their curves are shrunk towards the pooled one by the BLUPs).
  # Name as the Istat composition writes it; "" disables the section.
  region_focus = "Lombardia",

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
  hw_merge_windows = TRUE, # episodes whose windows touch the next one are one heatwave
  post_days   = 14,     # further days used to look at mortality displacement

  # ---- Report 1: adaptation across Italy -----------------------------------
  # Every city refitted in consecutive, non-overlapping windows (see R/adaptation.R).
  adaptation = list(
    # 5-year windows where the data allow: the years before COVID-19 give one full
    # window and four years, the years after it three and a half (2026 has data to
    # June, so it is folded into the last window). NA = last day with data.
    windows = list(c("2011-01-01", "2015-12-31"), c("2016-01-01", "2019-12-31"),
                   c("2023-01-01", NA)),
    summer_months = 6:9,       # reference summers of the standardised attributable fraction
    exclude_years = 2020:2022, # COVID-19: no reference summer, no window
    dose_pct      = 90         # heat dose = degrees above this city percentile, summed
  ),
  example_city  = "Milan",     # city of the daily example figures

  # ---- Report 2 and the app: rolling recalibration --------------------------
  # Each target month gets a curve refitted (Masselot's first stage + BLUPs) on
  # the deaths available before it. See R/rolling.R.
  rolling = list(
    train_years   = 3,            # years of deaths behind each curve
    anchor_years  = 0,            # > train_years adds the "shift" variant (city shape on
                                  # anchor_years + common recent shift); 0 = not fitted
    skip_years    = 2020:2022,    # COVID-19: never used to train a curve
    skip_targets  = TRUE,         # ... and not predicted either
    update        = "monthly",    # "monthly": refitted every target month on the deaths
                                  # up to the month before; "yearly": once a year, on the
                                  # deaths up to 31 December of the previous year
    target_months = 4:9,          # April-September: the summers the app shows
    first_target_year = NULL,     # NULL = as early as the Istat history allows (a target
                                  # without train_years of data before it is skipped)
    pool       = "italy",         # cities in the meta-analysis: "italy" (every city with
                                  # min_deaths in the window, plus the scope) or "region"
    min_deaths = 5000,            # Masselot's threshold, for the pooled cities only
    scope      = "all",           # cities given curves (the app needs every city):
                                  # "all" or "region"; the validation scores region_focus
    large_quantile = 2 / 3,       # episodes above this quantile of the published
                                  # prediction are "large"
    method = "window",            # "window" or "shift" (needs anchor_years)
    label = "Rolling 3-year Italian curves",
    label_shift = "10-year city curves + 3-year Italian shift"
  ),

  # ---- Report 2: implementation readiness ----------------------------------
  # Fixed before looking at the results; the verdict is computed from them.
  # core = TRUE criteria must all pass for "ready"; the others turn a pass into
  # "ready with caveats" when they fail.
  readiness = list(
    list(id = "ratio",    label = "Observed / predicted heat deaths, all episodes",
         lo = 0.80, hi = 1.25, core = TRUE),
    list(id = "ratio_lg", label = "Observed / predicted heat deaths, large episodes",
         lo = 0.80, hi = 1.25, core = TRUE),
    list(id = "slope",    label = "Calibration slope, all episodes",
         lo = 0.80, hi = 1.25, core = FALSE),
    list(id = "coverage", label = "Episodes inside their 95% prediction interval",
         lo = 0.85, hi = 1, core = FALSE),
    list(id = "pod",      label = "Probability of detection of significant excess (POD)",
         lo = 0.80, hi = 1, core = TRUE),
    list(id = "far",      label = "False alarm ratio (FAR)",
         lo = 0, hi = 0.30, core = TRUE),
    list(id = "hss",      label = "Heidke skill score",
         lo = 0.40, hi = 1, core = FALSE),
    list(id = "timing",   label = "Deaths follow temperature by (days, peak of the correlation)",
         lo = 0, hi = 2, core = TRUE)
  ),

  # ---- counterfactual baseline: deaths expected without heat --------------
  baseline_exclude_months = c(12, 1, 2, 3),  # influenza season, not used to fit the baseline
  baseline_hot_lag   = 3,     # days after a day above the MMT also excluded
  baseline_harmonics = 1,     # annual sine/cosine pairs (1, as EuroMOMO: summer has few reference days)
  baseline_trend_df_per_year = 0.5,
  baseline_min_daily = 2,     # age groups with fewer deaths/day use the 20+ model x their share
  baseline_cold_offset = TRUE, # remove the published cold effect from reference days (B = deaths at MMT)
  # alternative baselines for the audit of the system validation (overrides of
  # the settings above); the evaluation is repeated with each
  baseline_variants = list(
    `Two seasonal harmonics` = list(baseline_harmonics = 2),
    `No cold offset` = list(baseline_cold_offset = FALSE),
    `Reference days April-October only` = list(baseline_exclude_months = c(11, 12, 1, 2, 3))
  ),

  # ---- temperature forecasts (Open-Meteo Previous Runs API) ------------------
  forecast_model    = "ecmwf_ifs025",
  forecast_start    = as.Date("2024-01-01"),
  forecast_leads    = 1:7,
  forecast_timezone = "GMT",
  fc_debias         = FALSE,  # TRUE: subtract each city x month x lead mean forecast
                              # error (estimated outside the episodes being scored)

  # ---- influenza activity (input to FluMOMO) --------------------------------
  influenza = list(
    country = "Italy", country_code = "ITA",
    file = "data/raw/influenza_activity_IT.csv",   # built by the pipeline
    flunet_file = "data/raw/flunet_italy.csv",     # optional manual FluNet export
    use_existing = FALSE,  # TRUE: keep the file at `file` as it is, download nothing
    allow_zero = FALSE     # TRUE: run FluMOMO with IA = 0 if nothing can be obtained
  ),

  # ---- regional excess mortality with the official FluMOMO code -------------
  flumomo = list(
    region = "Lombardia", country_code = "IT-LOM",
    code_dir = "R/flumomo",
    code_url = "https://euromomo.eu/uploads/data/FluMOMO_version_4_2_R.zip",
    work_dir = "flumomo_run",
    years = 2015:2026, agegrp = 4L, # 4 = total, 3 = 65+
    ia_lags = 2, et_lags = 2, ia_restricted = TRUE,
    chart_year = 2026,
    provinces = data.frame(
      NUTS3 = c("ITC4C", "ITC47", "ITC46", "ITC4D", "ITC42", "ITC41",
                "ITC48", "ITC4B", "ITC4A", "ITC43", "ITC49", "ITC44"),
      name = c("Milano", "Brescia", "Bergamo", "Monza e Brianza", "Como", "Varese",
               "Pavia", "Mantova", "Cremona", "Lecco", "Lodi", "Sondrio"),
      prov = c("015", "017", "016", "108", "013", "012",
               "018", "020", "019", "097", "098", "014"),
      lat = c(45.4642, 45.5416, 45.6983, 45.5845, 45.8081, 45.8206,
              45.1847, 45.1564, 45.1332, 45.8566, 45.3140, 46.1700),
      lon = c(9.1900, 10.2118, 9.6773, 9.2744, 9.0852, 8.8251,
              9.1582, 10.7914, 10.0227, 9.3977, 9.5030, 9.8700),
      pop3 = c(3220000, 1260000, 1110000, 875000, 600000, 880000,
               540000, 410000, 355000, 335000, 230000, 180000),
      stringsAsFactors = FALSE)
  ),

  seed = 20260922
)
