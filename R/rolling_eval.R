# =============================================================================
#  Scoring the rolling curves, and the prediction intervals behind the grey
#  bands of the report and of the app.
#
#  Every set of curves is scored on the same windows, against the same
#  baseline, with the same uncertainty model, so differences between sets come
#  from the curves alone:
#    - curve uncertainty: draws of the coefficients from their covariance
#      (BLUP covariance for the recalibrated sets);
#    - baseline uncertainty: one log-normal multiplier per draw, shared by the
#      days of a window (an error in the level of B is common to its days);
#    - Poisson noise with the overdispersion of the baseline model (negative
#      binomial), for the prediction interval of the observed deaths.
#  Observed deaths outside the 95% prediction interval are "abnormal": more
#  (or fewer) deaths than the curves can account for.
# =============================================================================

PUBLISHED_LABEL <- "Published (Masselot et al.)"

#' Negative binomial draws with a quasi-Poisson dispersion
#'
#' @param mu Means.
#' @param phi Dispersion (variance / mean), recycled.
#' @return Integer vector.
rnb_draw <- function(mu, phi) {
  mu <- pmax(mu, 1e-9); phi <- pmax(phi, 1 + 1e-6)
  stats::rnbinom(length(mu), mu = mu, size = mu / (phi - 1))
}

#' Daily series of the counterfactual summed over age groups
#'
#' @param cf Output of [fit_counterfactual()].
#' @param codes Cities to keep (`NULL` = all).
#' @return data.table `URAU_CODE`, `date`, `deaths`, `tmean`, `B`, `B_logsd`,
#'   `phi`, `E_pub` (deaths expected with the published curves).
city_daily <- function(cf, codes = NULL) {
  d <- if (is.null(codes)) cf$daily else cf$daily[URAU_CODE %in% codes]
  d <- data.table::copy(d)
  if (!"phi" %in% names(d)) d[, phi := 1]
  if (!"B_logsd" %in% names(d)) d[, B_logsd := 0]
  out <- d[, .(deaths = sum(deaths), tmean = tmean[1], B = sum(B),
               B_logsd = sum(B * B_logsd) / max(sum(B), 1e-9),
               phi = sum(B * phi) / max(sum(B), 1e-9),
               E_pub = sum(B * exp(logrr))), by = .(URAU_CODE, date)]
  data.table::setkey(out, URAU_CODE, date)
  out
}

#' Curve provider for a fixed set of curves (one curve per city)
#'
#' @param rc Output of [recalibrate_curves()].
#' @return Function `(code, dates)` returning `id` (per date) and `curves`.
curves_fixed <- function(rc) function(cd, dates) {
  cv <- rc$curves[[cd]]
  if (is.null(cv)) return(NULL)
  list(id = rep(cd, length(dates)), curves = stats::setNames(list(cv), cd))
}

#' Curve provider for the rolling curves (the window in force on each date)
#'
#' @param roll Output of [rolling_collect()].
#' @return Function `(code, dates)` returning `id` (per date, `NA` where no
#'   curve exists) and `curves`.
curves_rolling <- function(roll) function(cd, dates) {
  wid <- rolling_lookup(roll, dates)
  ids <- unique(wid[!is.na(wid)])
  cvs <- stats::setNames(lapply(ids, function(w) roll$curves[[w]][[cd]]), ids)
  wid[wid %in% ids[vapply(cvs, is.null, TRUE)]] <- NA_character_
  list(id = wid, curves = Filter(Negate(is.null), cvs))
}

#' Log relative risks of a run of days, each day with its own curve
#'
#' @param spec Output of [erf_spec()].
#' @param tmean Temperatures.
#' @param ids Curve id of each day (no `NA`).
#' @param curves Named list of curves (`beta`, `vcov`).
#' @param nsim Draws.
#' @param cache Environment holding coefficient draws, so a curve used by
#'   many windows is drawn once (and its windows share the same draws).
#' @param prefix Key prefix for the cache (the city code).
#' @return List with `lr` (point, centred on each curve's MMT), `lrd`
#'   (`nsim x n`) and `hot` (day at or above the curve's MMT).
curve_days <- function(spec, tmean, ids, curves, nsim, cache, prefix) {
  n <- length(tmean)
  lr <- numeric(n); lrd <- matrix(0, nsim, n); hot <- logical(n)
  for (u in unique(ids)) {
    j <- which(ids == u); cv <- curves[[u]]
    key <- paste(prefix, u)
    if (is.null(cache[[key]])) {
      assign(key, list(draws = rmvn(nsim, cv$beta, cv$vcov), mmt = erf_mmt(spec, cv$beta)),
             envir = cache)
    }
    cc <- cache[[key]]
    X <- erf_basis(tmean[j], spec) - rep(erf_basis(cc$mmt, spec), each = length(j))
    lr[j] <- as.numeric(X %*% cv$beta)
    lrd[, j] <- cc$draws %*% t(X)
    hot[j] <- tmean[j] >= cc$mmt
  }
  list(lr = lr, lrd = lrd, hot = hot)
}

#' Score a set of curves on the counterfactual windows
#'
#' @param daily Output of [city_daily()].
#' @param wins Windows (`ev$windows` rows) with a `phi` column.
#' @param provider Curve provider ([curves_fixed()], [curves_rolling()]).
#' @param label Name of the curve set.
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return data.table, one row per window the set could score.
score_windows <- function(daily, wins, provider, label, erf, cfg) {
  set.seed(cfg$seed)
  cache <- new.env(parent = emptyenv()); specs <- list()
  out <- vector("list", nrow(wins))
  for (k in seq_len(nrow(wins))) {
    w <- wins[k]; cd <- w$URAU_CODE
    d <- daily[URAU_CODE == cd & date >= w$start & date <= w$end]
    if (!nrow(d)) next
    pv <- provider(cd, d$date)
    if (is.null(pv) || anyNA(pv$id)) next
    if (is.null(specs[[cd]])) specs[[cd]] <- erf_spec(erf, cd, cfg)
    cv <- curve_days(specs[[cd]], d$tmean, pv$id, pv$curves, cfg$nsim, cache, cd)
    bm <- exp(stats::rnorm(cfg$nsim, 0, mean(d$B_logsd)))
    rrd <- exp(cv$lrd); rr <- exp(cv$lr)
    P <- sum((d$B * (rr - 1))[cv$hot])
    Pd <- bm * as.numeric((rrd[, cv$hot, drop = FALSE] - 1) %*% d$B[cv$hot])
    E <- sum(d$B * rr)
    Ed <- bm * as.numeric(rrd %*% d$B)
    o <- rnb_draw(Ed, w$phi)
    q <- function(x) stats::quantile(x, c(.025, .975), names = FALSE)
    out[[k]] <- data.table::data.table(
      type = w$type, id = w$id, URAU_CODE = cd, year = w$year, start = w$start,
      end = w$end, O = w$O, B = w$B, X = w$X, v_pub = w$v, phi = w$phi, curves = label,
      P = P, P_lo = q(Pd)[1], P_hi = q(Pd)[2], E = E, var_E = stats::var(Ed),
      pi_lo = q(o)[1], pi_hi = q(o)[2],
      windows_used = paste(unique(pv$id), collapse = " "))
  }
  data.table::rbindlist(out)
}

#' Out-of-sample comparison of published, fixed and rolling curves
#'
#' Scores every set on the same windows of every source (a counterfactual
#' over a period with its evaluation), keeping only the windows that every set
#' could score - a rolling curve exists only from the first target the Istat
#' history allows, and not in skipped years. Episodes are split by size: an
#' episode is "large" when its heat deaths predicted with the published
#' curves are in the top `1 - cfg$rolling$large_quantile` share. Size is set
#' by the published prediction, a function of temperatures and baseline only,
#' so it does not select on the observed deaths and is the same subset for
#' every set of curves.
#'
#' The fixed recalibrated sets are fitted on 2011-2019 windows, so they are
#' scored only where they are out of sample: a source with `recals = FALSE`
#' (the pre-COVID period) compares the published and the rolling curves only.
#' Calibration is reported per period and, for the sets present in every
#' period, over all periods together.
#'
#' @param sources List of `list(label, cf, ev, recals)`; `recals = FALSE`
#'   leaves the fixed recalibrated sets out of that period.
#' @param roll One rolling curve set (`rolling_collect()$window` or `$shift`), or a
#'   list of them (e.g. the whole output of [rolling_collect()]): each is scored.
#' @param recals List of [recalibrate_curves()] outputs (fixed windows).
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @param codes Cities evaluated.
#' A window has a *significant excess* when its observed deaths exceed the
#' upper end of the 95% range of deaths expected without heat
#' ([noheat_band()]); a set of curves *forecasts* one when its expected
#' deaths E exceed that same bound. The two give the warning contingency
#' table of `warnings`.
#'
#' @return List with `windows` (long, with `nh_lo`, `nh_hi`, `obs_sig`,
#'   `fc_sig`), `calibration` (by period and over all periods: all windows,
#'   large episodes, episodes with significant excess), `warnings`,
#'   `by_year`, `displacement`, `large_threshold`, `curve_sets`.
rolling_evaluation <- function(sources, roll, recals, erf, cfg, codes) {
  rows <- list(); nhs <- list()
  rolls <- if (!is.null(roll$windows)) list(roll) else Filter(Negate(is.null), roll)
  recals <- Filter(Negate(is.null), recals)
  for (s in sources) {
    if (is.null(s$ev) || is.null(s$cf)) next
    w <- s$ev$windows[type %in% c("episode", "summer", "post") & URAU_CODE %in% codes]
    if (!nrow(w)) next
    w[, phi := pmax(v / pmax(E, 1e-9), 1)]
    daily <- city_daily(s$cf, codes)
    nhs[[length(nhs) + 1]] <- noheat_band(w, daily, cfg)[, period := s$label]
    rows[[length(rows) + 1]] <- w[, .(period = s$label, type, id, URAU_CODE, year, start,
                                      end, O, B, X, v_pub = v, phi,
                                      curves = PUBLISHED_LABEL, P, P_lo, P_hi, E, var_E,
                                      pi_lo, pi_hi, windows_used = NA_character_)]
    if (!isFALSE(s$recals)) for (rc in recals) {
      r <- score_windows(daily, w, curves_fixed(rc), rc$label, erf, cfg)
      if (nrow(r)) rows[[length(rows) + 1]] <- r[, period := s$label]
    }
    for (rl in rolls) {
      r <- score_windows(daily, w, curves_rolling(rl), rl$label, erf, cfg)
      if (nrow(r)) rows[[length(rows) + 1]] <- r[, period := s$label]
    }
  }
  long <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(long)) return(NULL)
  # within a period, only the windows every set of that period could score
  long[, nset := data.table::uniqueN(curves), by = period]
  both <- long[!is.na(P), .(N = .N, nset = nset[1]), by = .(period, type, id)][
    N == nset, .(period, type, id)]
  long <- long[both, on = .(period, type, id)][, nset := NULL]
  long[, `:=`(v = phi * E)]
  long[curves == PUBLISHED_LABEL, v := v_pub]
  long[, `:=`(in_pi = O >= pi_lo & O <= pi_hi, z = (O - E) / sqrt(v + var_E),
              abnormal = data.table::fifelse(O > pi_hi, "above",
                                             data.table::fifelse(O < pi_lo, "below", "within")))]
  # significant excess: deaths above the 95% range of the no-heat expectation,
  # observed (obs_sig) and predicted by each set of curves (fc_sig)
  long <- merge(long, data.table::rbindlist(nhs), by = c("period", "type", "id"), all.x = TRUE)
  long[, `:=`(obs_sig = O > nh_hi, fc_sig = E > nh_hi)]

  sz <- long[curves == PUBLISHED_LABEL & type == "episode", .(period, id, P_pub = P)]
  thr <- if (nrow(sz)) stats::quantile(sz$P_pub, cfg$rolling$large_quantile, names = FALSE) else Inf
  long <- merge(long, sz[, .(period, id, P_pub, type = "episode")],
                by = c("period", "type", "id"), all.x = TRUE)
  long[, size := data.table::fifelse(type == "episode" & !is.na(P_pub) & P_pub >= thr,
                                     "Large", "Other")]

  cal_one <- function(x) {
    if (nrow(x) < 3) return(NULL)
    m <- stats::lm(X ~ 0 + P, data = x, weights = 1 / pmax(x$v_pub, 1e-6))
    ci <- suppressWarnings(stats::confint(m))
    data.table::data.table(n = nrow(x), observed = sum(x$X), predicted = sum(x$P),
                           ratio = sum(x$X) / sum(x$P), slope = stats::coef(m)[1],
                           lo = ci[1, 1], hi = ci[1, 2],
                           r = suppressWarnings(stats::cor(x$X, x$P)),
                           mae = mean(abs(x$X - x$P)), coverage_95 = mean(x$in_pi),
                           above_pi = mean(x$abnormal == "above"),
                           mean_z = mean(x$z))
  }
  ord <- unique(c(PUBLISHED_LABEL, vapply(recals, `[[`, "", "label"),
                  vapply(rolls, `[[`, "", "label")))
  ord <- intersect(ord, unique(long$curves))
  pers <- unique(long$period)
  groups <- stats::setNames(lapply(pers, function(p) long[period == p]), pers)
  if (length(pers) > 1) {
    in_all <- long[, data.table::uniqueN(period), by = curves][V1 == length(pers), curves]
    groups[["All periods"]] <- long[curves %in% in_all]
  }
  cal <- data.table::rbindlist(lapply(names(groups), function(p) {
    x <- groups[[p]]
    data.table::rbindlist(list(
      x[type %in% c("episode", "summer"), cal_one(.SD), by = .(type, curves)][
        , subset := "All windows"],
      x[type == "episode" & size == "Large", cal_one(.SD), by = .(type, curves)][
        , subset := "Large episodes"],
      x[type == "episode" & obs_sig == TRUE, cal_one(.SD), by = .(type, curves)][
        , subset := "Episodes with significant excess"]), fill = TRUE)[, period := p]
  }), fill = TRUE)
  cal[, `:=`(curves = factor(curves, ord), period = factor(period, names(groups)))]
  by_year <- long[type %in% c("episode", "summer"),
                  .(windows = .N, X = sum(X), P = sum(P), coverage_95 = mean(in_pi)),
                  by = .(type, year, curves)]
  by_year[, curves := factor(curves, ord)]
  disp <- long[type %in% c("episode", "post"), .(X = sum(X), P = sum(P)),
               by = .(period, type, curves)]
  disp[, curves := factor(curves, ord)]
  warn <- data.table::rbindlist(lapply(names(groups), function(p)
    groups[[p]][type == "episode", warning_scores(obs_sig, fc_sig), by = curves][
      , period := p]))
  warn[, `:=`(curves = factor(curves, ord), period = factor(period, names(groups)))]
  list(windows = long[], calibration = cal[order(period, subset, type, curves)],
       warnings = warn[order(period, curves)],
       by_year = by_year[order(type, year, curves)],
       displacement = disp[order(period, curves, type)], large_threshold = thr,
       curve_sets = ord, periods = names(groups))
}

#' 95% range of a window's deaths if there were no heat
#'
#' The deaths expected without heat (B) with the uncertainty of its level and
#' the negative binomial noise of the counts: observed deaths above its upper
#' end are a significant excess, whatever the curves say.
#'
#' @param wins Windows (`ev$windows` rows) with `B` and `phi`.
#' @param daily Output of [city_daily()] (for the baseline uncertainty).
#' @param cfg Configuration list.
#' @return data.table `type`, `id`, `nh_lo`, `nh_hi`.
noheat_band <- function(wins, daily, cfg) {
  set.seed(cfg$seed)
  data.table::rbindlist(lapply(seq_len(nrow(wins)), function(k) {
    w <- wins[k]
    d <- daily[URAU_CODE == w$URAU_CODE & date >= w$start & date <= w$end]
    if (!nrow(d)) return(NULL)
    o <- rnb_draw(w$B * exp(stats::rnorm(cfg$nsim, 0, mean(d$B_logsd))), w$phi)
    q <- stats::quantile(o, c(.025, .975), names = FALSE)
    data.table::data.table(type = w$type, id = w$id, nh_lo = q[1], nh_hi = q[2])
  }))
}

#' Contingency table and skill scores of a yes/no warning
#'
#' @param obs,fc Logical: event observed / event forecast.
#' @return One-row data.table: hits, false alarms, misses, correct negatives,
#'   probability of detection (POD), false alarm ratio (FAR), critical success
#'   index (CSI), frequency bias, Peirce (PSS) and Heidke (HSS) skill scores.
warning_scores <- function(obs, fc) {
  ok <- !is.na(obs) & !is.na(fc); obs <- obs[ok]; fc <- fc[ok]
  a <- sum(obs & fc); b <- sum(!obs & fc); c <- sum(obs & !fc); d <- sum(!obs & !fc)
  dv <- function(x, y) if (y > 0) x / y else NA_real_
  data.table::data.table(
    n = a + b + c + d, hits = a, false_alarms = b, misses = c, correct_negatives = d,
    pod = dv(a, a + c), far = dv(b, a + b), csi = dv(a, a + b + c), bias = dv(a + b, a + c),
    pss = dv(a, a + c) - dv(b, b + d),
    hss = dv(2 * (a * d - b * c), (a + c) * (c + d) + (a + b) * (b + d)))
}

#' Daily prediction band of one city under the rolling curves
#'
#' For each day: deaths expected without heat (B), deaths predicted with the
#' rolling curve in force that day and with the published curves, the 95%
#' interval of the prediction (curve and baseline uncertainty) and the 95%
#' prediction interval of the observed count (adding Poisson noise with the
#' baseline's overdispersion). `flag` says whether the observed deaths fall
#' outside that interval - the days the report and the app mark as abnormal.
#' With `k > 1` the same is computed for k-day moving sums, which is what a
#' 7-day moving-average chart needs (the noise of a sum is smaller, relative
#' to its mean, than that of a day).
#'
#' @param cf Output of [fit_counterfactual()].
#' @param roll Output of [rolling_collect()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @param code Urban Audit code.
#' @param from,to Period.
#' @param k Days in the moving window (1 = daily values).
#' @return data.table, or `NULL` if no day has a rolling curve.
rolling_daily_band <- function(cf, roll, erf, cfg, code, from, to, k = 1L) {
  set.seed(cfg$seed)
  d <- city_daily(cf, code)[date >= from & date <= to]
  if (!nrow(d)) return(NULL)
  pv <- curves_rolling(roll)(code, d$date)
  ok <- !is.na(pv$id)
  if (!any(ok)) return(NULL)
  d <- d[ok]; ids <- pv$id[ok]
  spec <- erf_spec(erf, code, cfg)
  cv <- curve_days(spec, d$tmean, ids, pv$curves, cfg$nsim,
                   new.env(parent = emptyenv()), code)
  bm <- exp(stats::rnorm(cfg$nsim, 0, mean(d$B_logsd)))
  Ed <- sweep(exp(cv$lrd), 2, d$B, "*") * bm
  O <- matrix(rnb_draw(as.vector(Ed), rep(d$phi, each = cfg$nsim)), cfg$nsim)
  roll_k <- function(m) if (k <= 1) m else
    t(apply(m, 1, function(r) data.table::frollsum(r, k, align = "center")))
  qs <- function(m, p) apply(roll_k(m), 2, stats::quantile, p, na.rm = TRUE, names = FALSE)
  sk <- function(x) if (k <= 1) x else data.table::frollsum(x, k, align = "center")
  out <- data.table::data.table(
    URAU_CODE = code, date = d$date, window_id = ids, tmean = d$tmean,
    observed = sk(d$deaths), B = sk(d$B), pred = sk(d$B * exp(cv$lr)),
    pred_pub = sk(d$E_pub), pred_lo = qs(Ed, .025), pred_hi = qs(Ed, .975),
    pi_lo = qs(O, .025), pi_hi = qs(O, .975))
  if (k > 1) {
    cols <- c("observed", "B", "pred", "pred_pub", "pred_lo", "pred_hi", "pi_lo", "pi_hi")
    out[, (cols) := lapply(.SD, function(x) x / k), .SDcols = cols]
  }
  out[, flag := data.table::fifelse(observed > pi_hi, "above",
                                    data.table::fifelse(observed < pi_lo, "below", "within"))]
  out[!is.na(pi_lo)]
}

#' Summers of one city under the rolling curves, for the report figure
#'
#' @param sources List of `list(label, cf, ev)` (counterfactuals by period).
#' @param roll Output of [rolling_collect()].
#' @param erf Output of [read_erf_bundle()].
#' @param cities City table (for the name).
#' @param cfg Configuration list.
#' @param codes Cities evaluated (the example city must be one of them).
#' @return List with `data` and `city`, or `NULL`.
rolling_example <- function(sources, roll, erf, cities, cfg, codes) {
  cd <- cities[URAU_CODE %in% codes & norm_name(name) == norm_name(cfg$example_city), URAU_CODE]
  if (!length(cd)) cd <- codes[1]
  cd <- cd[1]
  out <- list()
  for (s in sources) {
    if (is.null(s$cf)) next
    yrs <- unique(data.table::year(s$cf$daily[URAU_CODE == cd, date]))
    for (y in yrs) {
      b <- rolling_daily_band(s$cf, roll, erf, cfg, cd, as.Date(sprintf("%d-06-01", y)),
                              as.Date(sprintf("%d-09-30", y)), k = 1L)
      if (!is.null(b)) out[[length(out) + 1]] <- b[, `:=`(year = y, period = s$label)]
    }
  }
  if (!length(out)) return(NULL)
  list(data = data.table::rbindlist(out), city = cities[URAU_CODE == cd, name][1])
}
