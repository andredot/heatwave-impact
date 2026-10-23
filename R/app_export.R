#' Weekly FluMOMO series for selected age groups
#'
#' Keeps the columns the app needs: the baseline, the baseline plus the
#' influenza contribution (the expected series without extreme temperature),
#' and the deaths FluMOMO attributes to extreme temperature.
#'
#' @param res Output of [run_region_flumomo()], or `NULL`.
#' @param agegrps Age groups to sum.
#' @return data.table, or `NULL`.
flumomo_weekly <- function(res, agegrps) {
  if (is.null(res) || !nrow(res)) return(NULL)
  r <- as.data.frame(res)
  r <- r[r$agegrp %in% agegrps & !is.na(r$year) & !is.na(r$week), ]
  d <- data.table::as.data.table(r)[, .(deaths = sum(deaths), EB = sum(EB),
                                        expected = sum(EB + EdIA), EdET = sum(EdET)),
                                    by = .(year, week)]
  d[, week_start := ISOweek::ISOweek2date(sprintf("%d-W%02d-1", year, week))]
  d[order(week_start)]
}

#' Curves of every city in the form the app reads
#'
#' @param codes Urban Audit codes.
#' @param erf Output of [read_erf_bundle()].
#' @param recals List of [recalibrate_curves()] outputs.
#' @param rolling Output of [rolling_collect()], or `NULL`.
#' @param cfg Configuration list.
#' @return Named list (one element per city): `spec`, `published` (by age
#'   group), `recalibrated` (by set), `rolling` (`index` of windows with the
#'   day each ends, `curves` by window), each curve with `beta`, `vcov`, `mmt`.
bundle_curves <- function(codes, erf, recals, rolling, cfg) {
  lapply(stats::setNames(codes, codes), function(cd) {
    sp <- erf_spec(erf, cd, cfg)
    pub <- lapply(stats::setNames(AGE_GROUPS, AGE_GROUPS), function(g) {
      b <- erf_coef(erf, cd, g)
      list(beta = b, vcov = erf_vcov(erf, cd, g), mmt = erf_mmt(sp, b))
    })
    rec <- list()
    for (rc in Filter(Negate(is.null), recals)) if (!is.null(rc$curves[[cd]])) {
      b <- rc$curves[[cd]]$beta
      rec[[rc$label]] <- list(beta = b, vcov = rc$curves[[cd]]$vcov, mmt = erf_mmt(sp, b))
    }
    # rolling curves: one per training window, with the date each window ends
    rol <- NULL
    if (!is.null(rolling) && nrow(rolling$windows)) {
      ids <- rolling$windows$window_id
      has <- vapply(ids, function(w) !is.null(rolling$curves[[w]][[cd]]), TRUE)
      if (any(has)) {
        idx <- rolling$windows[has, .(window_id, train_end, label)]
        rol <- list(index = idx[order(train_end)],
                    curves = stats::setNames(lapply(idx$window_id, function(w) {
                      cv <- rolling$curves[[w]][[cd]]
                      list(beta = cv$beta, vcov = cv$vcov, mmt = erf_mmt(sp, cv$beta),
                           source = cv$source)
                    }), idx$window_id))
      }
    }
    list(spec = list(knots = sp$knots, bound = sp$bound, degree = sp$degree,
                     p99 = sp$p99, pred_grid = sp$pred_grid, pred_pct = sp$pred_pct,
                     hw_threshold = pct_value(erf$tdist[URAU_CODE == cd], cfg$hw_pct)),
         published = pub, recalibrated = rec, rolling = rol)
  })
}

#' Build the data bundle used by the Shiny app
#'
#' Everything the app needs, computed once by the pipeline so that the app
#' only has to download temperatures: city list with region (NUTS-2, which in
#' Italy are the administrative regions), simplified boundaries for the map,
#' the published and the recalibrated curves, the daily observed deaths and
#' the baseline of deaths expected without heat, and the city x month
#' alignment of our ERA5-Land series with Masselot's.
#'
#' @param cf Output of [fit_counterfactual()].
#' @param recals List of [recalibrate_curves()] outputs.
#' @param erf Output of [read_erf_bundle()].
#' @param cities City table.
#' @param crosswalk Output of [finalize_crosswalk()].
#' @param gisco Output of [download_gisco()].
#' @param fc Output of [fetch_forecasts_all()] (archived forecasts).
#' @param era5 Raw ERA5-Land point series ([fetch_era5_all()]).
#' @param check Output of [compare_exposure()].
#' @param data_end Last day with mortality data.
#' @param cfg Configuration list.
#' @param path Destination `.rds`.
#' @param flumomo_region,flumomo_cities Outputs of the FluMOMO runs.
#' @param istat_file,mort Mortality inputs.
#' @param istat_cities Output of [read_istat_cities()]; when given, it decides
#'   each city's region, so the app's region list matches the report's.
#' @param rolling Output of [rolling_collect()]; its curves are exported city
#'   by city with the window index, so the app can use, on every date, the
#'   curve a forecaster would have had then.
#' @return `path`.
export_app_bundle <- function(cf, recals, erf, cities, crosswalk, gisco, check,
                              fc, era5, data_end, cfg, path = "app/app_data.rds",
                              flumomo_region = NULL, flumomo_cities = NULL,
                              istat_file = NULL, mort = NULL, istat_cities = NULL,
                              rolling = NULL) {
  ensure_dir(dirname(path))
  codes <- sort(unique(cf$daily$URAU_CODE))
  ct <- cities[URAU_CODE %in% codes]

  # ---- region (NUTS-2) and boundaries for the map ---------------------------
  nuts_file <- file.path(cfg$data_dir, "gisco", sprintf("nuts2_%s.gpkg", cfg$country))
  if (!file.exists(nuts_file)) {
    n <- if (requireNamespace("giscoR", quietly = TRUE)) {
      giscoR::gisco_get_nuts(year = "2021", nuts_level = 2, epsg = 4326,
                             resolution = "20", country = cfg$country)
    } else {
      x <- sf::st_read(paste0("https://gisco-services.ec.europa.eu/distribution/v2/nuts/",
                              "geojson/NUTS_RG_20M_2021_4326_LEVL_2.geojson"), quiet = TRUE)
      x[x$CNTR_CODE == cfg$country, ]
    }
    sf::st_write(n, nuts_file, quiet = TRUE, delete_dsn = TRUE)
  }
  nuts <- sf::st_make_valid(sf::st_read(nuts_file, quiet = TRUE))
  names(nuts)[names(nuts) == "NAME_LATN"] <- "region"
  pts <- sf::st_as_sf(ct[, .(URAU_CODE, lon, lat)], coords = c("lon", "lat"), crs = 4326)
  # one region per city even if boundaries overlap slightly
  idx <- suppressMessages(sf::st_within(pts, nuts))
  ct[, region := vapply(idx, function(i)
    if (length(i)) as.character(nuts$region[i[1]]) else NA_character_, character(1))]
  # the Istat composition names the region of every comune, so it decides where
  # a city belongs; the app's region list then matches the report's
  if (!is.null(istat_cities)) {
    cr <- city_regions(istat_cities, crosswalk$map)
    if (!is.null(cr) && nrow(cr)) {
      ct <- merge(ct, cr[, .(URAU_CODE, region_istat = region)], by = "URAU_CODE",
                  all.x = TRUE)
      ct[!is.na(region_istat), region := region_istat][, region_istat := NULL]
    }
  }
  ct[is.na(region), region := "(non assegnata)"]

  city_geom <- NULL
  sel <- unique(crosswalk$checks[URAU_CODE %in% codes,
                                 .(URAU_CODE, level = poly_level, poly_code)])
  for (f in gisco[grepl("^urau_", basename(gisco))]) {
    lvl <- sub("^urau_(.*)_([0-9]{4})_[A-Z]{2}\\.gpkg$", "\\1 \\2", basename(f))
    want <- sel[level == lvl]
    if (!nrow(want)) next
    u <- sf::st_read(f, quiet = TRUE)
    names(u) <- sub("^URAU_ID$", "URAU_CODE", names(u))
    u <- u[u$URAU_CODE %in% want$poly_code, "URAU_CODE"]
    if (!nrow(u)) next
    u$URAU_CODE <- want$URAU_CODE[match(u$URAU_CODE, want$poly_code)]
    city_geom <- rbind(city_geom, u)
  }
  simplify <- function(x, d) sf::st_make_valid(sf::st_simplify(x, dTolerance = d,
                                                               preserveTopology = TRUE))

  # ---- curves ---------------------------------------------------------------
  curves <- bundle_curves(codes, erf, recals, rolling, cfg)

  # ---- daily deaths and baseline --------------------------------------------
  # every observed day of every city, so the app can show mortality before the
  # period the counterfactual is fitted on
  deaths_long <- if (!is.null(mort))
    mort[, .(deaths = sum(deaths)), by = .(URAU_CODE, date)][order(URAU_CODE, date)]
  else NULL
  daily <- cf$daily[, .(URAU_CODE, agegroup, date, deaths, B,
                        B_logsd = if ("B_logsd" %in% names(cf$daily)) B_logsd else 0,
                        tmean, logrr, mmt)]
  data.table::setkey(daily, URAU_CODE, date, agegroup)
  # overdispersion of daily deaths around the baseline, by city: the Poisson
  # noise of the grey prediction band
  phi <- if ("phi" %in% names(cf$daily))
    cf$daily[, .(phi = sum(B * phi) / max(sum(B), 1e-9)), by = URAU_CODE] else NULL

  # ---- systematic forecast error, by city, month and lead -------------------
  # The report shows the forecasts running warm on hot days; since the curve is
  # convex, that bias inflates predicted deaths. These offsets let the app
  # subtract it before converting temperatures into deaths.
  fc_bias <- merge(fc, era5, by = c("URAU_CODE", "date"))[
    , .(bias = mean(tmean_fc - tmean_raw), sd = stats::sd(tmean_fc - tmean_raw), n = .N),
    by = .(URAU_CODE, month = data.table::month(date), lead)][n >= 10]

  out <- list(
    cities = ct[, .(URAU_CODE, name, region, lat, lon, pop)],
    fc_bias = fc_bias, forecast_model = cfg$forecast_model,
    regions = simplify(nuts[, "region"], 1000),
    city_geom = if (!is.null(city_geom)) simplify(city_geom, 500) else NULL,
    curves = curves, daily = daily, deaths_long = deaths_long, delta = check$delta,
    phi = phi,
    rolling = if (!is.null(rolling)) list(
      label = rolling$label, method = rolling$method, update = rolling$update,
      train_years = rolling$train_years, anchor_years = rolling$anchor_years,
      skip_years = rolling$skip_years, windows = rolling$windows,
      cities = names(Filter(function(x) !is.null(x$rolling), curves))) else NULL,
    region = list(
      name = cfg$flumomo$region,
      cities = unique(crosswalk$map[substr(PRO_COM, 1, 3) %in% cfg$flumomo$provinces$prov,
                                    URAU_CODE]),
      daily = if (!is.null(istat_file))
        istat_daily_area(istat_file, cfg$flumomo$provinces$prov, 0, cfg$flumomo$years)
      else NULL,
      # observed deaths of the validated cities over the same years: the daily
      # table only covers the test period, so the EuroMOMO tab would otherwise
      # show nothing before it
      cities_daily = if (!is.null(istat_file))
        istat_daily_area(istat_file, cfg$flumomo$provinces$prov, 0, cfg$flumomo$years,
                         procom = unique(crosswalk$map[
                           substr(PRO_COM, 1, 3) %in% cfg$flumomo$provinces$prov, PRO_COM]))
      else NULL),
    flumomo = list(
      region = flumomo_weekly(flumomo_region, c(2L, 3L)),
      cities = flumomo_weekly(flumomo_cities, c(2L, 3L)),
      region_all = flumomo_weekly(flumomo_region, 4L)),
    data_end = data_end, era5_timezone = cfg$era5_timezone,
    bias_correct = isTRUE(cfg$bias_correct),
    built = Sys.Date())
  saveRDS(out, path, compress = "xz")
  path
}
