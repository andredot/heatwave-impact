#' ERA5-Land grid cells inside each city
#'
#' Masselot averaged the ERA5-Land cells whose centre falls inside the city
#' boundary; taking a single point instead is a measurement error that
#' flattens a refitted exposure-response curve. This returns the centres of
#' the 0.1-degree ERA5-Land cells inside each city's polygon (at most
#' `cfg$era5_max_points`, spread evenly), falling back to the city point when
#' the polygon holds no cell centre.
#'
#' @param gisco Output of [download_gisco()].
#' @param crosswalk Output of [finalize_crosswalk()].
#' @param cities data.table with `URAU_CODE`, `lat`, `lon`.
#' @param cfg Configuration list.
#' @return data.table `URAU_CODE`, `lat`, `lon`, `n_cells`.
city_grid_points <- function(gisco, crosswalk, cities, cfg) {
  sel <- unique(crosswalk$checks[analysed == TRUE, .(URAU_CODE, level, poly_code)])
  out <- list()
  for (f in gisco[grepl("^urau_", basename(gisco))]) {
    lvl <- sub("^urau_(.*)_([0-9]{4})_[A-Z]{2}\\.gpkg$", "\\1 \\2", basename(f))
    if (!any(sel$level == lvl)) next
    u <- sf::st_read(f, quiet = TRUE)
    names(u) <- sub("^URAU_ID$", "URAU_CODE", names(u))
    want <- sel[level == lvl]
    u <- u[u$URAU_CODE %in% want$poly_code, ]
    if (!nrow(u)) next
    for (i in seq_len(nrow(u))) {
      bb <- sf::st_bbox(u[i, ])
      grid <- expand.grid(
        lon = seq(floor(bb["xmin"] * 10) / 10, ceiling(bb["xmax"] * 10) / 10, by = 0.1),
        lat = seq(floor(bb["ymin"] * 10) / 10, ceiling(bb["ymax"] * 10) / 10, by = 0.1))
      pts <- sf::st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
      inside <- grid[lengths(sf::st_within(pts, u[i, ])) > 0, , drop = FALSE]
      cd <- want$URAU_CODE[match(u$URAU_CODE[i], want$poly_code)]
      if (!nrow(inside)) next
      if (nrow(inside) > cfg$era5_max_points)
        inside <- inside[round(seq(1, nrow(inside), length.out = cfg$era5_max_points)), ]
      out[[length(out) + 1]] <- data.table::data.table(
        URAU_CODE = cd, lat = inside$lat, lon = inside$lon, n_cells = nrow(inside))
    }
  }
  pts <- data.table::rbindlist(out)
  miss <- setdiff(cities$URAU_CODE, pts$URAU_CODE)
  if (length(miss))
    pts <- rbind(pts, cities[URAU_CODE %in% miss, .(URAU_CODE, lat, lon, n_cells = 1L)])
  pts[]
}

#' Daily ERA5-Land mean temperature, averaged over several points
#'
#' One request with all the coordinates of a city (Open-Meteo accepts
#' comma-separated lists) and a simple average of the cells, as in Masselot's
#' city-polygon average.
#'
#' @param lat,lon Coordinates (vectors).
#' @param start,end Dates.
#' @param tz Time zone used to define the day.
#' @return data.table `date`, `tmean`.
fetch_era5_point <- function(lat, lon, start, end, tz = "GMT") {
  one <- function(la, lo) {
    url <- sprintf(paste0("https://archive-api.open-meteo.com/v1/archive",
                          "?latitude=%s&longitude=%s&start_date=%s&end_date=%s",
                          "&daily=temperature_2m_mean&models=era5_land&timezone=%s"),
                   paste(sprintf("%.4f", la), collapse = ","),
                   paste(sprintf("%.4f", lo), collapse = ","),
                   format(start), format(end), tz)
    era5_daily_list(get_json(url))
  }
  d <- one(lat, lon)
  # A multi-coordinate answer that did not parse into one series per point is
  # never accepted silently: fall back to one request per cell.
  bad <- length(d) != length(lat) || any(vapply(d, function(x) length(x$tmean) == 0, TRUE))
  if (bad && length(lat) > 1) {
    d <- unlist(lapply(seq_along(lat), function(i) one(lat[i], lon[i])), recursive = FALSE)
  }
  d <- Filter(function(x) length(x$tmean) > 0, d)
  if (!length(d)) stop("Open-Meteo returned no daily series for ", format(start), "-", format(end))
  dates <- d[[1]]$date
  m <- do.call(cbind, lapply(d, function(x) x$tmean[match(dates, x$date)]))
  data.table::data.table(date = dates, tmean = rowMeans(m, na.rm = TRUE))
}

#' Normalise an Open-Meteo archive answer into one series per location
#'
#' With one coordinate the answer is a list; with several it is an array,
#' which `jsonlite` may simplify to a data frame whose `daily` column holds
#' either lists or matrices. All three shapes are handled here, because a
#' silent mis-parse would drop whole cities from the analysis.
#'
#' @param js Parsed JSON.
#' @return List of `list(date, tmean)`.
era5_daily_list <- function(js) {
  pick <- function(x, i) if (is.list(x)) x[[i]] else if (is.matrix(x)) x[i, ] else x
  as_series <- function(tm, tv) list(date = as.Date(unlist(tm)), tmean = as.numeric(unlist(tv)))
  if (is.data.frame(js) && !is.null(js$daily)) {
    d <- js$daily
    return(lapply(seq_len(nrow(js)), function(i)
      as_series(pick(d$time, i), pick(d$temperature_2m_mean, i))))
  }
  if (!is.null(js$daily) && !is.null(js$daily$time))
    return(list(as_series(js$daily$time, js$daily$temperature_2m_mean)))
  if (is.list(js))
    return(lapply(js, function(x) as_series(x$daily$time, x$daily$temperature_2m_mean)))
  list()
}

#' ERA5-Land series for all cities, with a file cache
#'
#' One file per city and period in `cache_dir/era5`; existing files are
#' reused, so an interrupted run (e.g. API quota) resumes where it stopped.
#' Coordinates are the city points in Masselot's metadata; the difference with
#' his city-polygon averages is measured by [compare_exposure()].
#'
#' @param cities data.table with `URAU_CODE`, `lat`, `lon`: one row per city
#'   (point exposure) or several rows per city ([city_grid_points()]).
#' @param start,end Period.
#' @param cfg Configuration list.
#' @return data.table `URAU_CODE`, `date`, `tmean_raw`.
fetch_era5_all <- function(cities, start, end, cfg) {
  dir <- ensure_dir(file.path(cfg$cache_dir, "era5"))
  codes <- unique(cities$URAU_CODE)
  out <- lapply(codes, function(cd) {
    ci <- cities[URAU_CODE == cd]
    # the number of cells is part of the cache key: switching between point and
    # polygon exposure must trigger a new download, not reuse the old file
    f <- file.path(dir, sprintf("%s_%s_%s_n%d.csv", cd, start, end, nrow(ci)))
    if (!file.exists(f)) {
      d <- fetch_era5_point(ci$lat, ci$lon, start, end, cfg$era5_timezone)
      if (!nrow(d)) stop("Empty ERA5-Land series for ", cd)   # never cache an empty file
      data.table::fwrite(d, f)
      Sys.sleep(1)
    }
    d <- data.table::fread(f)
    d[, `:=`(URAU_CODE = cd, date = as.Date(date))]
    d[, .(URAU_CODE, date, tmean_raw = tmean)]
  })
  out <- data.table::rbindlist(out)
  out <- out[!is.na(tmean_raw)]
  lost <- setdiff(codes, unique(out$URAU_CODE))
  if (length(lost))
    stop("No ERA5-Land series for ", length(lost), " cities (",
         paste(utils::head(lost, 5), collapse = ", "), "). Delete their files in ",
         dir, " and run again.")
  out
}

#' Compare our ERA5-Land series with Masselot's
#'
#' On the overlap years, compares the point series downloaded from
#' Open-Meteo with the city-polygon series in `era5series.csv`. Returns
#' per-city agreement statistics and city x month mean differences, which
#' [build_exposure()] can remove.
#'
#' @param ours Output of [fetch_era5_all()] for the overlap period.
#' @param masselot Output of [read_era5_masselot()].
#' @return List with `summary` (by city) and `delta` (by city and month).
compare_exposure <- function(ours, masselot) {
  x <- merge(ours, masselot, by = c("URAU_CODE", "date"))
  x[, diff := tmean_raw - tmean]
  summ <- x[, .(n = .N, bias = mean(diff), mae = mean(abs(diff)),
                rmse = sqrt(mean(diff^2)), r = stats::cor(tmean_raw, tmean),
                slope = stats::coef(stats::lm(tmean_raw ~ tmean))[2]), by = URAU_CODE]
  delta <- x[, .(delta = mean(diff)), by = .(URAU_CODE, month = data.table::month(date))]
  list(summary = summ, delta = delta, daily = x)
}

#' Exposure series used in the analysis
#'
#' Optionally removes the city x month mean difference with Masselot's
#' series, so the temperatures sit on the same scale as the distribution the
#' ERF knots were computed on.
#'
#' @param ours Output of [fetch_era5_all()] for the analysis period.
#' @param check Output of [compare_exposure()].
#' @param cfg Configuration list.
#' @return data.table `URAU_CODE`, `date`, `tmean_raw`, `tmean`.
build_exposure <- function(ours, check, cfg) {
  x <- data.table::copy(ours)
  x[, month := data.table::month(date)]
  x <- merge(x, check$delta, by = c("URAU_CODE", "month"), all.x = TRUE)
  x[is.na(delta), delta := 0]
  x[, tmean := if (isTRUE(cfg$bias_correct)) tmean_raw - delta else tmean_raw]
  data.table::setkey(x, URAU_CODE, date)
  x[, .(URAU_CODE, date, tmean_raw, tmean)]
}
