# =============================================================================
#  Adaptation to heat in Italian cities, 2011 onwards (Report 1).
#
#  The published curves describe Italian mortality of 2001-2010. Here the same
#  first-stage model is refitted on consecutive, non-overlapping windows of
#  Istat deaths (cfg$adaptation$windows, five years where the data allow), every
#  city in every window, and the
#  city curves are pooled and shrunk as in Masselot's second stage. Three
#  measures of heat vulnerability are compared across windows:
#    1. the heat relative risk at the 99th percentile, against the published
#       curve (city differences pooled by random-effects meta-analysis, trend
#       by meta-regression on the window midpoint);
#    2. the minimum mortality temperature;
#    3. the standardised heat-attributable fraction: each window's curves
#       applied to the same summers and the same deaths, so that differences
#       between windows come from vulnerability alone, not from hotter or
#       cooler summers or from population change.
#  A curve-free check completes them: observed excess deaths of every summer
#  against that summer's heat dose, window by window.
# =============================================================================

#' Fit one adaptation window: every city, then BLUPs
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo Exposure covering the window (and `lag_max` days before).
#' @param erf Output of [read_erf_bundle()].
#' @param start,end Window.
#' @param cfg Configuration list.
#' @return List with `label`, `start`, `end`, `tests` (output of
#'   [run_curve_tests()] on deaths 20+), `pooled` (difference with the
#'   published curves, [pool_curve_differences()]), `curves` (code -> BLUP
#'   `beta`, `vcov`, `source`), `i2`, `n_fitted`.
adaptation_window <- function(mort, expo, erf, start, end, cfg) {
  tst <- run_curve_tests(mort, expo, erf, start, end, cfg, age_series = FALSE,
                         min_fit_deaths = 0)
  lab <- sprintf("%s-%s", format(start, "%Y"), format(end, "%Y"))
  f <- Filter(function(x) x$level == "all20", tst$fits)
  out <- list(label = lab, start = start, end = end, tests = tst,
              pooled = pool_curve_differences(tst), curves = list(), i2 = NA_real_,
              n_fitted = length(f))
  if (length(f) < 3) return(out)
  pm <- pool_mixmeta(do.call(rbind, lapply(f, `[[`, "coef")), lapply(f, `[[`, "vcov"))
  if (is.null(pm)) return(out)
  bl <- if (pm$method == "fixed") NULL else mixmeta::blup(pm$model, vcov = TRUE)
  psi <- if (is.null(pm$model$Psi)) 0 else pm$model$Psi
  pooled <- list(beta = as.numeric(stats::coef(pm$model)),
                 vcov = stats::vcov(pm$model) + psi, source = "pooled (not fitted)")
  codes <- intersect(unique(mort$URAU_CODE), erf$tdist$URAU_CODE)
  have <- vapply(f, `[[`, "", "URAU_CODE")
  for (cd in codes) {
    i <- match(cd, have)
    out$curves[[cd]] <- if (is.na(i)) pooled else if (is.null(bl))
      list(beta = as.numeric(f[[i]]$coef), vcov = f[[i]]$vcov, source = "own fit") else
      list(beta = as.numeric(bl[[i]]$blup), vcov = bl[[i]]$vcov, source = "BLUP")
  }
  out$i2 <- tryCatch(as.numeric(summary(pm$model)$i2stat[1]), error = function(e) NA_real_)
  out
}

#' Reference summers for the standardised attributable fraction
#'
#' Every city's daily temperatures over the summer months of all the years
#' analysed, and its mean daily deaths (20+) over the same days: the "same
#' summers, same deaths" every window's curves are applied to.
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo Exposure over the whole span ([rolling_exposure()]).
#' @param years Years used.
#' @param months Summer months.
#' @return data.table `URAU_CODE`, `date`, `tmean`, `D` (reference daily deaths).
reference_summers <- function(mort, expo, years, months) {
  x <- expo[data.table::year(date) %in% years & data.table::month(date) %in% months,
            .(URAU_CODE, date, tmean)]
  d <- mort[data.table::year(date) %in% years & data.table::month(date) %in% months,
            .(deaths = sum(deaths)), by = .(URAU_CODE, date)]
  D <- d[, .(D = mean(deaths)), by = URAU_CODE]
  merge(x, D, by = "URAU_CODE")
}

#' Standardised heat-attributable fraction of a set of curves
#'
#' For each city, the share of its reference summer deaths attributable to
#' heat (Masselot's formula: on days at or above the MMT, 1 - exp(-(f(T) -
#' f(MMT))), with the curve's own MMT), over the reference summers. The
#' national value is the deaths-weighted mean over cities. Draws of each
#' city's coefficients give its interval (cities drawn independently).
#'
#' @param curves Named list code -> `beta`, `vcov` (or a function `code ->
#'   curve object` for the published curves, see [published_curve_20()]).
#' @param ref Output of [reference_summers()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @param nsim Draws.
#' @return List with `cities` (code, AF, lo, hi, weight) and `draws`
#'   (national AF, one per draw).
standardised_af <- function(curves, ref, erf, cfg, nsim = cfg$nsim) {
  set.seed(cfg$seed)
  rows <- list(); nat <- 0; wsum <- 0
  for (cd in intersect(names(curves), unique(ref$URAU_CODE))) {
    r <- ref[URAU_CODE == cd]
    sp <- erf_spec(erf, cd, cfg)
    cv <- curves[[cd]]
    crv <- if (inherits(cv, "erf_curve")) cv else
      curve_single(sp, cv$beta, rmvn(nsim, cv$beta, cv$vcov))
    mmt <- sp$mmt_grid[which.min(curve_eval(crv, sp$mmt_grid))]
    ft <- curve_eval(crv, r$tmean) - curve_eval(crv, mmt)
    fd <- curve_eval(crv, r$tmean, draws = TRUE) -
      as.numeric(curve_eval(crv, mmt, draws = TRUE))
    hot <- r$tmean >= mmt
    af <- mean(ifelse(hot, 1 - exp(-ft), 0))
    afd <- rowMeans(sweep(1 - exp(-fd), 2, hot, "*"))
    w <- r$D[1]
    rows[[cd]] <- data.table::data.table(
      URAU_CODE = cd, af = af, lo = stats::quantile(afd, .025, names = FALSE),
      hi = stats::quantile(afd, .975, names = FALSE), weight = w, mmt = mmt)
    nat <- nat + w * afd; wsum <- wsum + w
  }
  list(cities = data.table::rbindlist(rows), draws = nat / wsum)
}

#' The published curve of a city for deaths 20+
#'
#' Death-weighted combination of the five published age-specific curves,
#' weights = the city's deaths by age over the reference summers.
#'
#' @param erf Output of [read_erf_bundle()].
#' @param mort Output of [aggregate_mortality()].
#' @param codes Cities.
#' @param cfg Configuration list.
#' @return Named list code -> `erf_curve`.
published_curve_20 <- function(erf, mort, codes, cfg) {
  w <- mort[URAU_CODE %in% codes, .(d = sum(deaths)), by = .(URAU_CODE, agegroup)]
  stats::setNames(lapply(codes, function(cd) {
    ww <- w[URAU_CODE == cd][match(AGE_GROUPS, agegroup), d]
    ww[is.na(ww)] <- 0
    curve_published(erf, cd, AGE_GROUPS, pmax(ww, 1), erf_spec(erf, cd, cfg), cfg$nsim)
  }), codes)
}

#' Adaptation analysis over consecutive windows
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo_hist Masselot's ERA5-Land series (windows ending by 2019).
#' @param expo_recent Our aligned series (later windows).
#' @param erf Output of [read_erf_bundle()].
#' @param data_end Last day with mortality data.
#' @param cfg Configuration list (uses `cfg$adaptation`).
#' @return List with `windows` (one [adaptation_window()] per window),
#'   `pooled` (difference with the published RR99, by window), `cities` (city
#'   x window results), `trend` (meta-regression per decade, all cities and
#'   the balanced panel), `mmt` (MMT and its percentile by window), `af`
#'   (standardised attributable fraction by window and for the published
#'   curves, with the trend per decade) and `reference` (what the
#'   standardisation used).
run_adaptation <- function(mort, expo_hist, expo_recent, erf, data_end, cfg) {
  ac <- cfg$adaptation
  expo_all <- rolling_exposure(expo_hist, expo_recent)
  wins <- list()
  for (w in ac$windows) {
    s0 <- as.Date(w[1]); e0 <- if (is.na(w[2])) data_end else min(as.Date(w[2]), data_end)
    if (s0 >= e0) next
    ex <- if (e0 <= as.Date("2019-12-31")) expo_hist else expo_recent
    wins[[length(wins) + 1]] <- adaptation_window(mort, ex, erf, s0, e0, cfg)
  }
  mid <- function(x) (as.numeric(format(x$start, "%Y")) + as.numeric(format(x$end, "%Y")) + 1) / 2
  pooled <- data.table::rbindlist(lapply(wins, function(x)
    cbind(window = x$label, midyear = mid(x), n_fitted = x$n_fitted, i2_blup = x$i2,
          x$pooled[group == "All cities"])), fill = TRUE)
  cities <- data.table::rbindlist(lapply(wins, function(x)
    cbind(window = x$label, midyear = mid(x),
          x$tests$results[level == "all20",
                          .(URAU_CODE, name, inmcc, deaths, status, below_threshold,
                            d_lr99, se_d_lr99, rr99_pub, rr99_new, rr99_new_lo, rr99_new_hi,
                            mmt_pub, mmt_new, p99)])), fill = TRUE)

  # trend of the city differences with the published curves, per decade
  trend_one <- function(d, lab) {
    d <- d[status == "fitted" & is.finite(se_d_lr99) & se_d_lr99 > 0]
    if (data.table::uniqueN(d$window) < 2 || nrow(d) < 5) return(NULL)
    m <- tryCatch(mixmeta::mixmeta(d_lr99 ~ I((midyear - 2010) / 10), S = d$se_d_lr99^2,
                                   data = d, random = ~ 1 | URAU_CODE, method = "reml"),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    b <- stats::coef(m)[2]; se <- sqrt(stats::vcov(m)[2, 2])
    data.table::data.table(sample = lab, cities = data.table::uniqueN(d$URAU_CODE),
                           fits = nrow(d), change_per_decade = exp(b),
                           lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
                           p = 2 * stats::pnorm(-abs(b / se)))
  }
  nwin <- data.table::uniqueN(cities$window)
  balanced <- cities[status == "fitted" & below_threshold == FALSE,
                     if (.N == nwin) .SD, by = URAU_CODE]
  trend <- data.table::rbindlist(list(
    trend_one(cities, "All cities, every window"),
    trend_one(balanced, sprintf("Cities above %s deaths in every window",
                                format(cfg$min_deaths, big.mark = ",")))))

  # MMT: BLUP curves, temperature and percentile of the city's 1990-2019 distribution
  mmt <- data.table::rbindlist(lapply(wins, function(x)
    data.table::rbindlist(lapply(names(x$curves), function(cd) {
      sp <- erf_spec(erf, cd, cfg); m <- erf_mmt(sp, x$curves[[cd]]$beta)
      data.table::data.table(window = x$label, midyear = mid(x), URAU_CODE = cd, mmt = m,
                             mmp = stats::approx(sp$pred_grid, sp$pred_pct, m, rule = 2)$y)
    }))))

  # standardised attributable fraction: same summers, same deaths, each window's curves
  yrs <- sort(unique(data.table::year(expo_all$date)))
  yrs <- yrs[yrs >= as.integer(format(min(mort$date), "%Y")) &
               yrs <= as.integer(format(data_end, "%Y")) & !yrs %in% ac$exclude_years]
  ref <- reference_summers(mort, expo_all, yrs, ac$summer_months)
  codes <- Reduce(intersect, lapply(wins, function(x) names(x$curves)))
  codes <- intersect(codes, unique(ref$URAU_CODE))
  pubs <- standardised_af(published_curve_20(erf, mort, codes, cfg), ref, erf, cfg)
  afw <- lapply(wins, function(x) standardised_af(x$curves[codes], ref, erf, cfg))
  q <- function(v) stats::quantile(v, c(.025, .975), names = FALSE)
  af <- data.table::rbindlist(c(
    list(data.table::data.table(window = "Published (2001-2010 data)", midyear = 2005.5,
                                af = mean(pubs$draws), lo = q(pubs$draws)[1],
                                hi = q(pubs$draws)[2])),
    lapply(seq_along(wins), function(i) data.table::data.table(
      window = wins[[i]]$label, midyear = mid(wins[[i]]), af = mean(afw[[i]]$draws),
      lo = q(afw[[i]]$draws)[1], hi = q(afw[[i]]$draws)[2]))))
  # trend of the national standardised AF per decade, draw by draw
  mids <- vapply(wins, mid, 0)
  D <- do.call(cbind, lapply(afw, `[[`, "draws"))
  sl <- apply(D, 1, function(y) stats::coef(stats::lm(log(y) ~ I((mids - 2010) / 10)))[2])
  af_trend <- data.table::data.table(change_per_decade = exp(stats::median(sl)),
                                     lo = exp(q(sl)[1]), hi = exp(q(sl)[2]))
  af_cities <- data.table::rbindlist(c(
    list(cbind(window = "Published (2001-2010 data)", pubs$cities)),
    lapply(seq_along(wins), function(i) cbind(window = wins[[i]]$label, afw[[i]]$cities))))

  list(windows = lapply(wins, function(x) x[c("label", "start", "end", "pooled", "i2",
                                              "n_fitted", "curves")]),
       tests = lapply(wins, function(x) x$tests$results),
       pooled = pooled, cities = cities, trend = trend, mmt = mmt,
       af = af, af_trend = af_trend, af_cities = af_cities,
       reference = list(years = yrs, months = ac$summer_months, cities = codes,
                        days = nrow(ref)))
}

#' Curve-free check: summer excess deaths against the summer's heat dose
#'
#' For every city and summer with a counterfactual baseline, the observed
#' excess X = O - B over June-September and the heat dose, the sum over the
#' summer days of the degrees above the city's `dose_pct` percentile of its
#' 1990-2019 distribution. Within each window, excess as a share of the
#' baseline is regressed on the dose through the origin (weighted by the
#' baseline): the slope is the excess mortality per 10 degree-days, which
#' falls over time if populations become more resilient. No curve enters this
#' measure, only deaths, temperatures and the baseline. Only complete summers
#' (all 122 days of June-September) are used.
#'
#' @param cfs List of [fit_counterfactual()] outputs (periods).
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `summers` (city x year) and `slopes` (by window).
summer_dose_response <- function(cfs, erf, cfg) {
  ac <- cfg$adaptation
  d <- data.table::rbindlist(lapply(Filter(Negate(is.null), cfs), function(cf)
    cf$daily[data.table::month(date) %in% 6:9,
             .(O = sum(deaths), B = sum(B), tmean = tmean[1]), by = .(URAU_CODE, date)]))
  if (!nrow(d)) return(NULL)
  thr <- data.table::data.table(URAU_CODE = unique(d$URAU_CODE))
  thr[, t_ref := vapply(URAU_CODE, function(cd)
    pct_value(erf$tdist[URAU_CODE == cd], ac$dose_pct), 0)]
  d <- merge(d, thr, by = "URAU_CODE")
  s <- d[, .(O = sum(O), B = sum(B), dose = sum(pmax(tmean - t_ref, 0)), days = .N),
         by = .(URAU_CODE, year = data.table::year(date))]
  # whole summers only: a summer cut by the end of the data has a partial dose
  s <- s[days == 122L][, days := NULL]
  s[, X := O - B]
  win <- data.table::rbindlist(lapply(ac$windows, function(w) data.table::data.table(
    window = sprintf("%s-%s", substr(w[1], 1, 4),
                     if (is.na(w[2])) "" else substr(w[2], 1, 4)),
    y0 = as.integer(substr(w[1], 1, 4)),
    y1 = if (is.na(w[2])) 9999L else as.integer(substr(w[2], 1, 4)))))
  s[, window := vapply(year, function(y) {
    k <- which(win$y0 <= y & win$y1 >= y)
    if (length(k)) win$window[k[1]] else NA_character_ }, "")]
  s <- s[!is.na(window)]
  s[endsWith(window, "-"), window := paste0(window, max(year)), by = window]
  slopes <- s[, {
    m <- stats::lm(I(X / B) ~ 0 + I(dose / 10), weights = B)
    ci <- suppressWarnings(stats::confint(m))
    .(summers = .N, cities = data.table::uniqueN(URAU_CODE), dose_mean = mean(dose),
      excess_pct = 100 * sum(X) / sum(B),
      slope_pct = 100 * stats::coef(m)[1], lo = 100 * ci[1, 1], hi = 100 * ci[1, 2])
  }, by = window]
  list(summers = s[], slopes = slopes[order(window)], dose_pct = ac$dose_pct)
}
