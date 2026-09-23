#' Archived temperature forecasts at fixed lead times (one point, one year)
#'
#' Uses the Open-Meteo Previous Runs API: `temperature_2m_previous_dayN`
#' is the value predicted N days (24N hours) before the valid time. Hourly
#' values are averaged into daily means; a day is kept only with 24 hours.
#' Lead 0 is the most recent run (close to an analysis).
#'
#' @param lat,lon Coordinates.
#' @param start,end Dates.
#' @param model Open-Meteo model name.
#' @param leads Integer lead days (1-7).
#' @param tz Time zone used to define the day.
#' @return data.table `date`, `lead`, `tmean_fc`.
fetch_forecast_point <- function(lat, lon, start, end, model, leads, tz = "GMT") {
  vars <- c("temperature_2m", sprintf("temperature_2m_previous_day%d", leads))
  url <- sprintf(paste0("https://previous-runs-api.open-meteo.com/v1/forecast",
                        "?latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s",
                        "&hourly=%s&models=%s&timezone=%s"),
                 lat, lon, format(start), format(end),
                 paste(vars, collapse = ","), model, tz)
  js <- get_json(url)
  h <- data.table::as.data.table(js$hourly)
  h[, date := as.Date(substr(time, 1, 10))]
  lv <- setdiff(names(h), c("time", "date"))
  long <- data.table::melt(h, id.vars = "date", measure.vars = lv,
                           variable.name = "var", value.name = "t")
  long[, lead := data.table::fifelse(var == "temperature_2m", 0L,
                                     as.integer(sub(".*previous_day", "", var)))]
  long[, .(tmean_fc = if (sum(!is.na(t)) == 24) mean(t) else NA_real_),
       by = .(date, lead)][!is.na(tmean_fc)]
}

#' Archived forecasts for all cities, with a file cache
#'
#' Requested one calendar year at a time and cached per city and year.
#'
#' @param cities data.table with `URAU_CODE`, `lat`, `lon`.
#' @param end Last date.
#' @param cfg Configuration list.
#' @return data.table `URAU_CODE`, `date`, `lead`, `tmean_fc`.
fetch_forecasts_all <- function(cities, end, cfg) {
  dir <- ensure_dir(file.path(cfg$cache_dir, "forecast"))
  yrs <- seq(as.integer(format(cfg$forecast_start, "%Y")), as.integer(format(end, "%Y")))
  out <- list()
  for (i in seq_len(nrow(cities))) for (y in yrs) {
    s <- max(cfg$forecast_start, as.Date(sprintf("%d-01-01", y)))
    e <- min(end, as.Date(sprintf("%d-12-31", y)))
    cd <- cities$URAU_CODE[i]
    f <- file.path(dir, sprintf("%s_%s_%s_%s.csv", cd, cfg$forecast_model, s, e))
    if (!file.exists(f)) {
      d <- fetch_forecast_point(cities$lat[i], cities$lon[i], s, e, cfg$forecast_model,
                                cfg$forecast_leads, cfg$forecast_timezone)
      data.table::fwrite(d, f)
      Sys.sleep(1)
    }
    d <- data.table::fread(f)
    d[, `:=`(URAU_CODE = cd, date = as.Date(date))]
    out[[length(out) + 1]] <- d
  }
  data.table::rbindlist(out, use.names = TRUE)
}

#' Temperature forecast skill by lead time
#'
#' Forecasts are scored against ERA5-Land at the same point (the temperature
#' the mortality model uses). A persistence forecast (the temperature observed
#' on the issue day) is the reference: skill = 1 - MAE / MAE_persistence.
#' Scores are given for all days, the warm season, and hot days (ERA5-Land
#' above the city's published 90th percentile). Heatwave-day detection is
#' scored with the probability of detection (POD) and false alarm ratio (FAR),
#' using the heatwave threshold of the episode definition.
#'
#' @param fc Output of [fetch_forecasts_all()].
#' @param era5 Output of [fetch_era5_all()] (raw point series).
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `scores` (by lead and subset), `detection` (by lead) and
#'   `by_city` (MAE by city and lead, warm season).
forecast_skill <- function(fc, era5, erf, cfg) {
  thr <- data.table::rbindlist(lapply(unique(fc$URAU_CODE), function(cd) {
    row <- erf$tdist[URAU_CODE == cd]
    data.table::data.table(URAU_CODE = cd, p90 = pct_value(row, 90),
                           phw = pct_value(row, cfg$hw_pct))
  }))
  x <- merge(fc[lead %in% cfg$forecast_leads], era5, by = c("URAU_CODE", "date"))
  x <- merge(x, thr, by = "URAU_CODE")
  pers <- data.table::copy(era5)[, .(URAU_CODE, issue = date, tpers = tmean_raw)]
  x[, issue := date - lead]
  x <- merge(x, pers, by = c("URAU_CODE", "issue"), all.x = TRUE)
  x[, `:=`(err = tmean_fc - tmean_raw, err_pers = tpers - tmean_raw,
           warm = data.table::month(date) %in% cfg$warm_months)]
  score <- function(d) d[, .(n = .N, bias = mean(err), mae = mean(abs(err)),
                             rmse = sqrt(mean(err^2)),
                             mae_persistence = mean(abs(err_pers), na.rm = TRUE)), by = lead]
  sc <- data.table::rbindlist(list(
    cbind(subset = "All days", score(x)),
    cbind(subset = "Warm season", score(x[warm == TRUE])),
    cbind(subset = "Hot days (> P90)", score(x[tmean_raw > p90]))))
  sc[, skill := 1 - mae / mae_persistence]
  det <- x[warm == TRUE, .(hits = sum(tmean_fc >= phw & tmean_raw >= phw),
                           misses = sum(tmean_fc < phw & tmean_raw >= phw),
                           false_alarms = sum(tmean_fc >= phw & tmean_raw < phw)), by = lead]
  det[, `:=`(pod = hits / (hits + misses), far = false_alarms / (hits + false_alarms))]
  byc <- x[warm == TRUE, .(mae = mean(abs(err)), bias = mean(err)), by = .(URAU_CODE, lead)]
  list(scores = sc[order(subset, lead)], detection = det[order(lead)], by_city = byc)
}
