# =============================================================================
#  Rolling recalibration: a curve for every forecast month.
#
#  The fixed recalibrated curves (R/recalibrate.R) showed that the level of
#  heat risk is what has gone out of date, and that it keeps falling: a window
#  ending in 2019 fits the test period better than one starting in 2011. A
#  forecasting system should therefore not freeze a curve at all. Here, for
#  every target month, the first stage of Masselot et al. is refitted on the
#  deaths that were available *before* that month - the last `train_years`
#  years of Istat data, skipping the COVID-19 years - and the city fits are
#  pooled and shrunk exactly as in the second stage. The published curves
#  (MCC Italy 2001-2010) are kept untouched as the reference the rolling curves
#  are evaluated against.
#
#  Example, monthly update, train_years = 3, data from January 2011:
#    July 2015   trained on 1 Jul 2012 - 30 Jun 2015
#    July 2024   trained on 1 Jan 2023 - 30 Jun 2024 + 1 Jul 2018 - 31 Dec 2019
#    July 2013   not predicted: three years before it are not in the Istat file
#
#  Every curve is out of sample for the month it is used in, so scoring it on
#  that month's episodes is an honest test of what a forecaster would have
#  had at the time.
# =============================================================================

#' Training windows of the rolling recalibration
#'
#' One target per month of `cfg$rolling$target_months`, from the first year
#' the Istat history can support to the season after the last data. Each
#' target is given the `train_years` years of deaths that end the day before
#' it (`update = "monthly"`) or on 31 December of the previous year
#' (`"yearly"`), never later than the last day with data, counting only days
#' outside `skip_years`. A target whose training days would start before the
#' first day of the mortality data is skipped, and so is a target in a skipped
#' year when `skip_targets` is set. Targets that end up with the same training
#' days share one window, which is fitted once.
#'
#' @param data_start,data_end First and last day with mortality data.
#' @param cfg Configuration list (uses `cfg$rolling`).
#' @return List with `targets` (one row per target month, with `window_id`
#'   and `status`) and `windows` (list, one element per distinct window:
#'   `id`, `train_start`, `train_end`, `n_days`, `blocks`, `label`).
rolling_windows <- function(data_start, data_end, cfg) {
  rc <- cfg$rolling
  skip <- as.integer(rc$skip_years)
  n_days <- round(365.25 * rc$train_years)
  y0 <- as.integer(format(data_start, "%Y"))
  if (!is.null(rc$first_target_year)) y0 <- max(y0, as.integer(rc$first_target_year))
  y1 <- as.integer(format(data_end, "%Y"))
  # once the data pass the last target month, the next season is the one ahead
  if (as.integer(format(data_end, "%m")) >= max(rc$target_months)) y1 <- y1 + 1L
  tg <- data.table::CJ(year = y0:y1, month = as.integer(rc$target_months))
  tg[, target := as.Date(sprintf("%d-%02d-01", year, month))]
  tg[, key := format(target, "%Y-%m")]
  tg[, train_end := if (identical(rc$update, "yearly"))
    as.Date(sprintf("%d-12-31", year - 1L)) else target - 1]
  tg[train_end > data_end, train_end := data_end]
  tg[, `:=`(window_id = NA_character_, status = "ok")]
  if (isTRUE(rc$skip_targets)) tg[year %in% skip, status := "COVID-19 year, not predicted"]

  wins <- list()
  for (te in sort(unique(tg[status == "ok", train_end]))) {
    te <- as.Date(te, origin = "1970-01-01")
    all_days <- seq(te, data_start, by = "-1 day")
    days <- all_days[!as.integer(format(all_days, "%Y")) %in% skip]
    if (length(days) < n_days) {
      tg[status == "ok" & train_end == te,
         status := sprintf("Istat history too short (needs %d years before %s)",
                           rc$train_years, format(te + 1))]
      next
    }
    days <- sort(days[seq_len(n_days)])
    br <- c(TRUE, diff(as.numeric(days)) > 1)
    blocks <- data.table::data.table(start = days[br],
                                     end = days[c(which(br)[-1] - 1, length(days))])
    id <- format(te, "w%Y%m%d")
    wins[[id]] <- list(
      id = id, train_start = min(days), train_end = te, n_days = length(days),
      blocks = blocks,
      label = paste(sprintf("%s to %s", format(blocks$start, "%b %Y"),
                            format(blocks$end, "%b %Y")), collapse = " + "))
    tg[status == "ok" & train_end == te, window_id := id]
  }
  list(targets = tg[], windows = unname(wins))
}

#' Is each date inside one of the blocks?
#'
#' @param dates Dates.
#' @param blocks data.table with `start`, `end`.
#' @return Logical vector.
in_blocks <- function(dates, blocks) {
  out <- logical(length(dates))
  for (k in seq_len(nrow(blocks)))
    out <- out | (dates >= blocks$start[k] & dates <= blocks$end[k])
  out
}

#' Masselot's first stage on a training window with gaps
#'
#' Same model as [fit_stage1()] - quasi-Poisson, cross-basis of temperature
#' (published knots) and lag, day of week, 7 df per year for time - except
#' that the time trend is a separate natural spline within each contiguous
#' block of training days, each with its own intercept. A single spline over
#' a window that jumps from 2019 to 2023 would spend degrees of freedom on
#' three years with no data and tie the two sides of the gap together. The
#' cross-basis is built on the continuous daily series, so the first days of
#' a block see their real lagged temperatures.
#'
#' @param dat data.table `date`, `y`, `tmean`, one row per day, starting
#'   `lag_max` days before the first training day.
#' @param inp Logical: training days (with deaths and temperature).
#' @param spec Output of [erf_spec()].
#' @param cfg Configuration list.
#' @return As [fit_stage1()], or `NULL` if the model fails.
fit_stage1_blocks <- function(dat, inp, spec, cfg) {
  cb <- suppressWarnings(dlnm::crossbasis(dat$tmean, lag = cfg$lag_max,
                         argvar = list(fun = "bs", degree = spec$degree, knots = spec$knots,
                                       Boundary.knots = spec$bound),
                         arglag = list(fun = "ns",
                                       knots = dlnm::logknots(cfg$lag_max, cfg$lag_nknots))))
  inp <- inp & stats::complete.cases(cb)
  if (sum(inp) < 365) return(NULL)
  idx <- which(inp)
  blk <- integer(nrow(dat))
  blk[idx] <- cumsum(c(TRUE, diff(idx) > 1))
  tm <- as.numeric(dat$date)
  trend <- do.call(cbind, lapply(sort(unique(blk[idx])), function(b) {
    rows <- which(blk == b)
    df <- max(1L, round(cfg$df_time_per_year * length(rows) / 365.25))
    m <- matrix(0, nrow(dat), df)
    m[rows, ] <- splines::ns(tm[rows], df = df)
    m
  }))
  block <- factor(blk)
  dow <- factor(weekdays(dat$date))
  y <- dat$y
  # one intercept per block (a factor with a single level would not fit)
  f <- if (length(unique(blk[idx])) > 1) y ~ cb + dow + block + trend else y ~ cb + dow + trend
  fit <- tryCatch(suppressWarnings(stats::glm(
    f, family = stats::quasipoisson, subset = inp)),
    error = function(e) NULL)
  if (is.null(fit) || !fit$converged) return(NULL)
  red <- tryCatch(suppressWarnings(
    dlnm::crossreduce(cb, fit, cen = stats::median(dat$tmean[inp]))),
    error = function(e) NULL)
  if (is.null(red) || anyNA(stats::coef(red))) return(NULL)
  list(coef = stats::coef(red), vcov = stats::vcov(red), converged = TRUE,
       dispersion = summary(fit)$dispersion, totdeath = sum(y[inp]))
}

#' Random-effects pooling that degrades gracefully
#'
#' REML first, as in Masselot's second stage; with few cities (a region on its
#' own) the 15 variance parameters may not be estimable, so the method of
#' moments and then a fixed-effect model are tried.
#'
#' @param Y Coefficients, one row per city.
#' @param S List of covariance matrices.
#' @return List with `model`, `method`, or `NULL`.
pool_mixmeta <- function(Y, S) {
  for (m in c("reml", "mm", "fixed")) {
    fit <- tryCatch(suppressWarnings(mixmeta::mixmeta(Y, S, method = m)),
                    error = function(e) NULL)
    if (!is.null(fit) && all(is.finite(stats::coef(fit)))) return(list(model = fit, method = m))
  }
  NULL
}

#' Heat relative risk at P99 of an all-ages curve, with its 95% interval
#'
#' @param spec Output of [erf_spec()].
#' @param beta Coefficients.
#' @param V Covariance matrix.
#' @return List with `mmt`, `rr99`, `lo`, `hi`.
rr99_of <- function(spec, beta, V) {
  m <- erf_mmt(spec, beta)
  d <- erf_basis(spec$p99, spec) - erf_basis(m, spec)
  lr <- as.numeric(d %*% beta)
  se <- sqrt(max(as.numeric(d %*% V %*% t(d)), 0))
  list(mmt = m, rr99 = exp(lr), lo = exp(lr - 1.96 * se), hi = exp(lr + 1.96 * se))
}

#' Published heat RR at P99 of a city, death-weighted over age groups
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param w Deaths by age group (order of `AGE_GROUPS`).
#' @param spec Output of [erf_spec()].
#' @return Numeric.
published_rr99 <- function(erf, code, w, spec) {
  w <- pmax(w, 1) / sum(pmax(w, 1))
  f <- function(t) {
    acc <- 0
    for (k in seq_along(AGE_GROUPS)) {
      b <- erf_coef(erf, code, AGE_GROUPS[k])
      m <- erf_mmt(spec, b)
      acc <- acc + w[k] * exp(as.numeric(erf_basis(t, spec) %*% b) -
                                as.numeric(erf_basis(m, spec) %*% b))
    }
    log(acc)
  }
  g <- f(spec$mmt_grid)
  exp(f(spec$p99) - min(g))
}

#' Fit one rolling window: first stage in every pooled city, then BLUPs
#'
#' @param win One element of `rolling_windows()$windows`.
#' @param mort Output of [aggregate_mortality()].
#' @param expo Exposure covering the window and `lag_max` days before it (see
#'   [rolling_exposure()]).
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @param pool Urban Audit codes entering the meta-analysis (cities below
#'   `cfg$rolling$min_deaths` in the window are left out unless in `keep`).
#' @param keep Codes whose curves are returned. They are always fitted,
#'   whatever their size, and shrunk by the BLUPs; a city whose model fails
#'   gets the pooled curve.
#' @return List with `id`, `train_end`, `label`, `n_pooled`, `method`,
#'   `curves` (code -> `beta`, `vcov`, `source`) and `summary`.
fit_rolling_window <- function(win, mort, expo, erf, cfg, pool, keep) {
  set.seed(cfg$seed)
  bl <- win$blocks
  first <- min(bl$start); last <- win$train_end
  days <- seq(first - cfg$lag_max, last, by = "day")
  intrain <- in_blocks(days, bl)
  fits <- list(); notes <- list()
  codes <- intersect(union(keep, pool), intersect(unique(mort$URAU_CODE),
                                                  unique(expo$URAU_CODE)))
  for (cd in codes) {
    ex <- expo[URAU_CODE == cd & date >= min(days) & date <= last, .(date, tmean)]
    y <- mort[URAU_CODE == cd & date >= first & date <= last, .(y = sum(deaths)), by = date]
    dat <- merge(data.table::data.table(date = days), ex, by = "date", all.x = TRUE)
    dat <- merge(dat, y, by = "date", all.x = TRUE)
    inp <- intrain & !is.na(dat$y) & !is.na(dat$tmean)
    td <- sum(dat$y[inp])
    if (td < cfg$rolling$min_deaths && !cd %in% keep) next
    f <- fit_stage1_blocks(dat, inp, erf_spec(erf, cd, cfg), cfg)
    if (is.null(f)) {
      notes[[cd]] <- "not converged"
      next
    }
    ages <- mort[URAU_CODE == cd & date >= first & date <= last & in_blocks(date, bl),
                 .(d = sum(deaths)), by = agegroup][match(AGE_GROUPS, agegroup), d]
    ages[is.na(ages)] <- 0
    fits[[cd]] <- c(f, list(URAU_CODE = cd, deaths = td, w = ages))
  }
  out <- list(id = win$id, train_start = win$train_start, train_end = win$train_end,
              label = win$label, blocks = bl, n_pooled = length(fits),
              method = NA_character_, curves = list(), summary = NULL)
  if (length(fits) < 3) return(out)
  pm <- pool_mixmeta(do.call(rbind, lapply(fits, `[[`, "coef")),
                     lapply(fits, `[[`, "vcov"))
  if (is.null(pm)) return(out)
  out$method <- pm$method
  bl_ <- if (pm$method == "fixed") NULL else mixmeta::blup(pm$model, vcov = TRUE)
  psi <- if (is.null(pm$model$Psi)) 0 else pm$model$Psi
  pooled <- list(beta = as.numeric(stats::coef(pm$model)),
                 vcov = stats::vcov(pm$model) + psi)
  rows <- list()
  for (cd in keep) {
    i <- match(cd, names(fits))
    cur <- if (!is.na(i) && !is.null(bl_)) {
      list(beta = as.numeric(bl_[[i]]$blup), vcov = bl_[[i]]$vcov, source = "BLUP")
    } else if (!is.na(i)) {
      list(beta = as.numeric(fits[[i]]$coef), vcov = fits[[i]]$vcov, source = "own fit")
    } else {
      c(pooled, list(source = if (is.null(notes[[cd]])) "pooled (no data)" else
        paste0("pooled (", notes[[cd]], ")")))
    }
    if (!cd %in% erf$tdist$URAU_CODE) next
    out$curves[[cd]] <- cur
    sp <- erf_spec(erf, cd, cfg)
    r <- rr99_of(sp, cur$beta, cur$vcov)
    w <- if (!is.na(i)) fits[[i]]$w else rep(1, length(AGE_GROUPS))
    rows[[length(rows) + 1]] <- data.table::data.table(
      window_id = win$id, train_end = win$train_end, URAU_CODE = cd,
      deaths = if (!is.na(i)) fits[[i]]$deaths else NA_real_, source = cur$source,
      mmt = r$mmt, rr99 = r$rr99, rr99_lo = r$lo, rr99_hi = r$hi,
      rr99_pub = published_rr99(erf, cd, w, sp))
  }
  out$summary <- data.table::rbindlist(rows)
  out
}

#' Exposure series for the rolling windows
#'
#' Masselot's own ERA5-Land series up to the end of 2019 (the exposure the
#' published curves were built on) and our aligned series afterwards.
#'
#' @param hist Output of [read_era5_masselot()].
#' @param recent Output of [build_exposure()].
#' @param switch Date from which `recent` is used.
#' @return data.table `URAU_CODE`, `date`, `tmean`.
rolling_exposure <- function(hist, recent, switch = as.Date("2020-01-01")) {
  x <- data.table::rbindlist(list(hist[date < switch, .(URAU_CODE, date, tmean)],
                                  recent[date >= switch, .(URAU_CODE, date, tmean)]))
  data.table::setkey(x, URAU_CODE, date)
  x
}

#' Collect the fitted windows into one object
#'
#' @param plan Output of [rolling_windows()].
#' @param fits List of [fit_rolling_window()] outputs.
#' @param cfg Configuration list.
#' @return List with `label`, `targets`, `windows` (only those with curves),
#'   `curves` (window id -> code -> curve) and `summary`.
rolling_collect <- function(plan, fits, cfg) {
  fits <- Filter(function(f) !is.null(f) && length(f$curves), fits)
  w <- data.table::rbindlist(lapply(fits, function(f) data.table::data.table(
    window_id = f$id, train_start = f$train_start, train_end = f$train_end,
    label = f$label, n_pooled = f$n_pooled, method = f$method)))
  tg <- data.table::copy(plan$targets)
  if (nrow(w)) tg[status == "ok" & !window_id %in% w$window_id,
                  status := "fewer than 3 cities could be fitted"]
  list(label = cfg$rolling$label, update = cfg$rolling$update,
       train_years = cfg$rolling$train_years, skip_years = cfg$rolling$skip_years,
       targets = tg,
       windows = if (nrow(w)) w[order(train_end)] else w,
       curves = stats::setNames(lapply(fits, `[[`, "curves"), vapply(fits, `[[`, "", "id")),
       summary = data.table::rbindlist(lapply(fits, `[[`, "summary"), fill = TRUE))
}

#' Rolling window in force on each date
#'
#' The curve a forecaster has on day d is the one fitted on the latest window
#' that ended before d. In monthly mode that is the curve of d's month; on a
#' day outside the target months it is the last one fitted. Dates with no
#' earlier window, or whose latest window is older than `max_age` days (e.g.
#' inside the COVID-19 gap), get `NA`.
#'
#' @param roll Output of [rolling_collect()].
#' @param dates Dates.
#' @param max_age Days after which a window is considered stale.
#' @return Character vector of window ids.
rolling_lookup <- function(roll, dates, max_age = 400) {
  w <- roll$windows
  if (is.null(w) || !nrow(w)) return(rep(NA_character_, length(dates)))
  te <- as.numeric(w$train_end)
  i <- findInterval(as.numeric(dates) - 1, te)
  id <- rep(NA_character_, length(dates))
  ok <- i > 0
  id[ok] <- w$window_id[i[ok]]
  id[ok & (as.numeric(dates) - te[pmax(i, 1)]) > max_age] <- NA_character_
  id
}
