# =============================================================================
# Validation of the Masselot et al. (2023) temperature-mortality curves for
# Italian cities with Istat daily deaths, and of Open-Meteo temperature
# forecasts as their input.
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
  tar_target(stage1, read_stage1(zenodo_files, erf$metadata)),
  tar_target(published, read_published_results(zenodo_files, erf$metadata$URAU_CODE)),
  tar_target(era5_masselot, read_era5_masselot(zenodo_files, erf$metadata$URAU_CODE)),

  # ---- do we reconstruct the published curves correctly? --------------------
  tar_target(metamodel, inspect_metamodel(zenodo_files)),
  tar_target(recon_check, check_reconstruction(erf, era5_masselot, published, config)),

  # ---- city definition ------------------------------------------------------
  tar_target(gisco_files, download_gisco(config), format = "file"),
  tar_target(overlay, overlay_cities(gisco_files, erf$metadata, config)),

  # ---- mortality ------------------------------------------------------------
  tar_target(istat_file, config$istat_csv, format = "file"),
  tar_target(istat, read_istat(istat_file, unique(overlay$members$PRO_COM))),
  tar_target(data_end, detect_data_end(istat, config$data_end)),
  tar_target(periods, make_periods(config, data_end)),
  tar_target(crosswalk, finalize_crosswalk(overlay, istat, erf, config)),
  tar_target(mort, aggregate_mortality(istat, crosswalk$map,
                                       min(config$repro_start, config$test_start), data_end)),
  tar_target(cities, erf$metadata[URAU_CODE %in% crosswalk$map$URAU_CODE,
                                  .(URAU_CODE, name, lat = as.numeric(lat),
                                    lon = as.numeric(lon), pop = as.numeric(pop), inmcc)]),
  tar_target(baseline_audit, audit_baseline(mort, erf, periods$test$start, periods$test$end)),

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

  # ---- 1. temperature forecast skill ----------------------------------------
  tar_target(forecasts, fetch_forecasts_all(cities, min(data_end, Sys.Date() - 7), config)),
  tar_target(fc_skill, forecast_skill(forecasts, era5_test, erf, config)),

  # ---- 2A. curve comparison (Masselot first stage refitted) -----------------
  tar_target(curves_test, run_curve_tests(mort, exposure_test, erf,
                                          periods$test$start, periods$test$end, config)),
  tar_target(pooled_test, pool_curve_differences(curves_test)),
  tar_target(blup_test, blup_curves(curves_test, erf, config)),
  # 2011-2019 uses Masselot's own exposure series (exact same temperatures)
  tar_target(curves_repro, if (config$run_reproduction)
    run_curve_tests(mort, era5_masselot, erf, periods$repro$start, periods$repro$end, config)),
  tar_target(pooled_repro, if (config$run_reproduction) pool_curve_differences(curves_repro)),
  tar_target(blup_repro, if (config$run_reproduction) blup_curves(curves_repro, erf, config)),

  tar_target(curves_periods, run_period_tests(mort, era5_masselot, exposure_test, erf,
                                              data_end, config)),

  # ---- 2B. heat deaths vs a counterfactual "no heat" baseline ---------------
  tar_target(episodes, detect_episodes(exposure_test, erf, periods$test$start,
                                       periods$test$end, config)),
  tar_target(counterfactual, fit_counterfactual(mort, exposure_test, episodes, erf,
                                                periods$test$start, periods$test$end, config)),
  tar_target(cf_eval, evaluate_counterfactual(counterfactual, config)),
  tar_target(recal_curves, lapply(config$recal_windows, function(w)
    recalibrate_curves(mort, era5_masselot, erf, w, config))),
  tar_target(curve_sets, compare_curve_sets(counterfactual, cf_eval, recal_curves, erf, config)),
  tar_target(example, example_series(counterfactual, recal_curves, erf, cities, config)),
  tar_target(fc_impact, forecast_episode_impact(counterfactual, cf_eval, forecasts,
                                                era5_test, exposure_check, erf, config)),

  # ---- regional excess mortality (official FluMOMO code) --------------------
  tar_target(influenza_file, build_influenza_activity(config), format = "file"),
  tar_target(flumomo_weather, build_flumomo_weather(
    config, as.Date(sprintf("%d-01-01", min(config$flumomo$years))),
    min(data_end, Sys.Date() - 7)), format = "file"),
  tar_target(flumomo_results, run_region_flumomo(config, istat_file, flumomo_weather,
                                                 influenza_file)),
  tar_target(flumomo_chart, plot_region_excess(flumomo_results, istat_file, config),
             format = "file"),
  # same code, restricted to the comuni of the validated cities of the region,
  # so its baseline can replace ours in test B
  tar_target(region_cities, crosswalk$map[substr(PRO_COM, 1, 3) %in%
                                            config$flumomo$provinces$prov]),
  tar_target(flumomo_cities, run_region_flumomo(config, istat_file, flumomo_weather,
                                                influenza_file,
                                                procom = unique(region_cities$PRO_COM),
                                                label = "cities")),
  tar_target(flumomo_sensitivity, flumomo_test_b(
    flumomo_cities, counterfactual, recal_curves, erf, istat_file,
    unique(region_cities$PRO_COM), unique(region_cities$URAU_CODE), config)),

  # ---- data bundle for the Shiny app ----------------------------------------
  tar_target(app_data, export_app_bundle(counterfactual, recal_curves, erf, cities,
                                         crosswalk, gisco_files, exposure_check,
                                         forecasts, era5_test, data_end, config,
                                         "app/app_data.rds",
                                         flumomo_results, flumomo_cities, istat_file),
             format = "file"),

  # ---- report ---------------------------------------------------------------
  tar_quarto(report, "reports//validation_report.qmd")
)
