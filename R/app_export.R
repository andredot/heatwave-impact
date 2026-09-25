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
#' @return `path`.
export_app_bundle <- function(cf, recals, erf, cities, crosswalk, gisco, check,
                              fc, era5, data_end, cfg, path = "app/app_data.rds",
                              flumomo_region = NULL, flumomo_cities = NULL,
                              istat_file = NULL) {
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
  ct[is.na(region), region := "(non assegnata)"]

  city_geom <- NULL
  sel <- unique(crosswalk$checks[URAU_CODE %in% codes, .(URAU_CODE, level, poly_code)])
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
  curves <- lapply(stats::setNames(codes, codes), function(cd) {
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
    list(spec = list(knots = sp$knots, bound = sp$bound, degree = sp$degree,
                     p99 = sp$p99, pred_grid = sp$pred_grid, pred_pct = sp$pred_pct),
         published = pub, recalibrated = rec)
  })

  # ---- daily deaths and baseline --------------------------------------------
  daily <- cf$daily[, .(URAU_CODE, agegroup, date, deaths, B,
                        B_logsd = if ("B_logsd" %in% names(cf$daily)) B_logsd else 0,
                        tmean, logrr, mmt)]
  data.table::setkey(daily, URAU_CODE, date, agegroup)

  # ---- systematic forecast error, by city, month and lead -------------------
  # The report shows the forecasts running warm on hot days; since the curve is
  # convex, that bias inflates predicted deaths. These offsets let the app
  # subtract it before converting temperatures into deaths.
  fc_bias <- merge(fc, era5, by = c("URAU_CODE", "date"))[
    , .(bias = mean(tmean_fc - tmean_raw), n = .N),
    by = .(URAU_CODE, month = data.table::month(date), lead)][n >= 10]

  out <- list(
    cities = ct[, .(URAU_CODE, name, region, lat, lon, pop)],
    fc_bias = fc_bias, forecast_model = cfg$forecast_model,
    regions = simplify(nuts[, "region"], 1000),
    city_geom = if (!is.null(city_geom)) simplify(city_geom, 500) else NULL,
    curves = curves, daily = daily, delta = check$delta,
    region = list(
      name = cfg$flumomo$region,
      cities = unique(crosswalk$map[substr(PRO_COM, 1, 3) %in% cfg$flumomo$provinces$prov,
                                    URAU_CODE]),
      daily = if (!is.null(istat_file))
        istat_daily_area(istat_file, cfg$flumomo$provinces$prov, 0, cfg$flumomo$years)
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
