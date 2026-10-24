# =============================================================================
# Heat and mortality in Italian cities, from the Masselot et al. (2023) curves.
# Two reports:
#   reports/adaptation_report.qmd   adaptation to heat across Italy since 2011
#   reports/system_validation.qmd   validation of the impact-based system
#                                   (rolling curves, recorded temperatures)
# and the data bundle of the Shiny app.
#
#   targets::tar_make()        # run (resumable: downloads are cached)
#   targets::tar_visnetwork()  # dependency graph
# =============================================================================

library(targets)
library(tarchetypes)

source("R//config.R")
# Only this project's own files: a vendored copy of the official FluMOMO code
# under R/ must not be sourced (it expects variables set by its own launcher).
tar_source(files = list.files("R", pattern = "[.][Rr]$", full.names = TRUE))

tar_option_set(
  packages = c("data.table", "dlnm", "splines", "mixmeta", "sf", "jsonlite"),
  format = "qs",
  seed = cfg$seed
)
if (!requireNamespace("qs2", quietly = TRUE)) tar_option_set(format = "rds")

list(
  # ---- configuration (tracked: editing config.R re-runs what depends on it) --
  tar_target(config, cfg),

  # ---- published ERFs -------------------------------------------------------
  tar_target(zenodo_files, download_zenodo(config$zenodo_dir, config$use_coef_simu),
             format = "file"),
  tar_target(erf, read_erf_bundle(zenodo_files, config$country)),
  tar_target(published, read_published_results(zenodo_files, erf$metadata$URAU_CODE)),
  tar_target(era5_masselot, read_era5_masselot(zenodo_files, erf$metadata$URAU_CODE)),
  # do we reconstruct the published curves correctly?
  tar_target(recon_check, check_reconstruction(erf, era5_masselot, published, config)),

  # ---- city definition ------------------------------------------------------
  tar_target(gisco_files, download_gisco(config), format = "file"),
  tar_target(overlay, overlay_cities(gisco_files, erf$metadata, config)),
  tar_target(lau_pop, read_lau_population(gisco_files)),
  # official Istat composition of the Cities 2021 (shipped with the repository)
  tar_target(istat_cities_file, locate_istat_cities(config), format = "file"),
  tar_target(istat_cities, read_istat_cities(istat_cities_file)),
  tar_target(istat_map, istat_city_map(istat_cities, erf$metadata, lau_pop)),
  # Eurostat CITY-LAU table, used when config$city_source is "official"
  tar_target(city_table, read_city_lau_table(config)),
  tar_target(city_map, switch(
    config$city_source,
    istat    = istat_map,
    official = official_city_map(city_table, erf$metadata,
                                 unique(lau_pop[, .(PRO_COM, lau_pop)])),
    NULL)),

  # ---- mortality ------------------------------------------------------------
  tar_target(istat_file, download_istat(config), format = "file"),
  tar_target(istat, read_istat(istat_file,
                               unique(c(overlay$members$PRO_COM, city_map$PRO_COM)))),
  tar_target(data_end, detect_data_end(istat, config$data_end)),
  tar_target(periods, make_periods(config, data_end)),
  tar_target(crosswalk, finalize_crosswalk(overlay, istat, erf, config, city_map)),
  tar_target(mort, aggregate_mortality(istat, crosswalk$map,
                                       min(config$repro_start, config$test_start), data_end)),
  tar_target(cities, erf$metadata[URAU_CODE %in% crosswalk$map$URAU_CODE,
                                  .(URAU_CODE, name, lat = as.numeric(lat),
                                    lon = as.numeric(lon), pop = as.numeric(pop), inmcc)]),
  tar_target(region_codes, region_city_codes(istat_cities, crosswalk$map,
                                             config$region_focus)),

  # ---- exposure -------------------------------------------------------------
  tar_target(era5_points, if (config$era5_mode == "polygon")
    city_grid_points(gisco_files, crosswalk, cities, config) else
      cities[, .(URAU_CODE, lat, lon, n_cells = 1L)]),
  tar_target(era5_overlap, fetch_era5_all(era5_points, config$overlap_start,
                                          config$overlap_end, config)),
  tar_target(era5_test, fetch_era5_all(era5_points, periods$test$start - config$lag_max - 1,
                                       min(data_end, Sys.Date() - 7), config)),
  tar_target(exposure_check, compare_exposure(era5_overlap, era5_masselot)),
  tar_target(exposure_test, build_exposure(era5_test, exposure_check, config)),
  # Masselot's own temperatures to 2019, our aligned series from 2020
  tar_target(expo_all, rolling_exposure(era5_masselot, exposure_test)),

  # ---- deaths expected without heat (counterfactual baseline) ---------------
  # fitted separately before and after COVID-19, every city
  tar_target(episodes, detect_episodes(exposure_test, erf, periods$test$start,
                                       periods$test$end, config)),
  tar_target(counterfactual, fit_counterfactual(mort, exposure_test, episodes, erf,
                                                periods$test$start, periods$test$end, config)),
  tar_target(cf_eval, evaluate_counterfactual(counterfactual, config)),
  tar_target(episodes_hist, detect_episodes(era5_masselot, erf, config$repro_start,
                                            config$repro_end, config)),
  tar_target(counterfactual_hist, fit_counterfactual(
    mort, era5_masselot, episodes_hist, erf, config$repro_start, config$repro_end, config)),
  tar_target(cf_eval_hist, evaluate_counterfactual(counterfactual_hist, config)),

  # ==== REPORT 1: adaptation across Italy ====================================
  # every city refitted in consecutive 3-year windows, BLUPs, standardised AF
  tar_target(adaptation, run_adaptation(mort, era5_masselot, exposure_test, erf,
                                        data_end, config)),
  # curve-free check: summer excess deaths against the summer's heat dose
  tar_target(dose_response, summer_dose_response(list(counterfactual_hist, counterfactual),
                                                 erf, config)),
  tar_target(timing_italy, list(
    `2011-2019` = lag_check(counterfactual_hist, episodes_hist),
    `2023 onwards` = lag_check(counterfactual, episodes))),

  # ==== REPORT 2: the impact-based system (rolling curves) ==================
  # curves for every city (the app needs them); scores for the focus region
  tar_target(roll_codes, if (identical(config$rolling$scope, "region")) region_codes
                         else cities$URAU_CODE),
  tar_target(eval_codes, if (length(region_codes)) region_codes else cities$URAU_CODE),
  # cities entering each window's meta-analysis
  tar_target(roll_pool, if (identical(config$rolling$pool, "region")) roll_codes
                        else cities$URAU_CODE),
  tar_target(roll_plan, rolling_windows(min(mort$date), data_end, config)),
  tar_target(roll_windows, roll_plan$windows, iteration = "list"),
  # one branch per training window: cached, so a new Istat release only fits
  # the windows it changes
  tar_target(roll_fits, fit_rolling_window(roll_windows, mort, expo_all, erf, config,
                                           roll_pool, roll_codes),
             pattern = map(roll_windows), iteration = "list"),
  tar_target(rolling_sets, rolling_collect(roll_plan, roll_fits, config, erf)),
  # the set the system uses (rolling$method)
  tar_target(rolling, if (!is.null(rolling_sets[[config$rolling$method]]))
    rolling_sets[[config$rolling$method]] else rolling_sets$window),
  tar_target(roll_eval, rolling_evaluation(
    list(list(label = sprintf("%s-%s", format(config$repro_start, "%Y"),
                              format(config$repro_end, "%Y")),
              cf = counterfactual_hist, ev = cf_eval_hist),
         list(label = sprintf("%s-%s", format(periods$test$start, "%Y"),
                              format(periods$test$end, "%Y")),
              cf = counterfactual, ev = cf_eval)),
    rolling, list(), erf, config, eval_codes)),
  tar_target(roll_example, rolling_example(
    list(list(label = "pre-COVID", cf = counterfactual_hist),
         list(label = "post-COVID", cf = counterfactual)),
    rolling, erf, cities, config, eval_codes)),
  tar_target(timing_region, list(
    `2011-2019` = lag_check(list(daily = counterfactual_hist$daily[URAU_CODE %in% eval_codes]),
                            episodes_hist[URAU_CODE %in% eval_codes]),
    `2023 onwards` = lag_check(list(daily = counterfactual$daily[URAU_CODE %in% eval_codes]),
                               episodes[URAU_CODE %in% eval_codes]))),
  tar_target(readiness, readiness_check(roll_eval, timing_region, rolling$label, config)),
  # is the yardstick sound? cool-day check and alternative baselines
  tar_target(baseline_audit, baseline_sensitivity(
    config$baseline_variants, mort, era5_masselot, exposure_test, episodes_hist, episodes,
    erf, config, eval_codes, rolling, periods, roll_eval, counterfactual_hist, counterfactual)),

  # ---- regional excess mortality (official FluMOMO code, app EuroMOMO tab) --
  tar_target(influenza_file, build_influenza_activity(config), format = "file"),
  tar_target(flumomo_code, download_flumomo(config), format = "file"),
  tar_target(flumomo_weather, build_flumomo_weather(
    config, as.Date(sprintf("%d-01-01", min(config$flumomo$years))),
    min(data_end, Sys.Date() - 7)), format = "file"),
  tar_target(flumomo_results, run_region_flumomo(config, istat_file, flumomo_weather,
                                                 influenza_file)),
  tar_target(region_cities, crosswalk$map[substr(PRO_COM, 1, 3) %in%
                                            config$flumomo$provinces$prov]),
  tar_target(flumomo_cities, run_region_flumomo(config, istat_file, flumomo_weather,
                                                influenza_file,
                                                procom = unique(region_cities$PRO_COM),
                                                label = "cities", code = flumomo_code)),

  # ---- data bundle for the Shiny app ----------------------------------------
  # forecasts only feed the app's time machine and scenarios
  tar_target(forecasts, fetch_forecasts_all(cities, min(data_end, Sys.Date() - 7), config)),
  tar_target(app_data, export_app_bundle(counterfactual, list(), erf, cities,
                                         crosswalk, gisco_files, exposure_check,
                                         forecasts, era5_test, data_end, config,
                                         "app/app_data.rds",
                                         flumomo_results, flumomo_cities, istat_file,
                                         mort, istat_cities, rolling,
                                         cf_hist = counterfactual_hist),
             format = "file"),

  # ---- reports --------------------------------------------------------------
  tar_quarto(report_adaptation, "reports//adaptation_report.qmd"),
  tar_quarto(report_validation, "reports//system_validation.qmd")
)
