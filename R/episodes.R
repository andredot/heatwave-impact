#' Detect heatwave episodes
#'
#' An episode is a spell of at least `hw_min_days` days in the warm season
#' with daily mean temperature at or above the city's `hw_pct` percentile of
#' its 1990-2019 distribution (the same distribution behind the ERF);
#' spells separated by at most `hw_max_gap` days are merged. Only episodes
#' whose evaluation window (episode + `tail_days`) ends within the period
#' are kept.
#'
#' @param expo Exposure data.table (`URAU_CODE`, `date`, `tmean`).
#' @param erf Output of [read_erf_bundle()].
#' @param start,end Period.
#' @param cfg Configuration list.
#' @return data.table with one row per episode.
detect_episodes <- function(expo, erf, start, end, cfg) {
  out <- lapply(unique(expo$URAU_CODE), function(cd) {
    thr <- pct_value(erf$tdist[URAU_CODE == cd], cfg$hw_pct)
    x <- merge(data.table::data.table(date = seq(start, end, by = "day")),
               expo[URAU_CODE == cd, .(date, tmean)], by = "date", all.x = TRUE)
    hot <- !is.na(x$tmean) & x$tmean >= thr & data.table::month(x$date) %in% cfg$warm_months
    r <- rle(hot); e_i <- cumsum(r$lengths); s_i <- e_i - r$lengths + 1
    sp <- data.table::data.table(s = x$date[s_i[r$values]], e = x$date[e_i[r$values]])
    if (!nrow(sp)) return(NULL)
    merged <- sp[1]
    if (nrow(sp) > 1) for (k in 2:nrow(sp)) {
      last <- nrow(merged)
      if (as.integer(sp$s[k] - merged$e[last]) - 1 <= cfg$hw_max_gap) {
        merged$e[last] <- sp$e[k]
      } else merged <- rbind(merged, sp[k])
    }
    merged[, `:=`(URAU_CODE = cd, thr = thr)]
    inwin <- function(k) x$date >= merged$s[k] & x$date <= merged$e[k]
    merged[, hot_days := vapply(seq_len(.N), function(k) sum(hot[inwin(k)]), 0L)]
    merged[, t_mean := vapply(seq_len(.N), function(k) mean(x$tmean[inwin(k)]), 0)]
    merged[, t_peak := vapply(seq_len(.N), function(k) max(x$tmean[inwin(k)]), 0)]
    merged[hot_days >= cfg$hw_min_days & e + cfg$tail_days <= end]
  })
  ep <- data.table::rbindlist(out)
  if (!nrow(ep)) return(ep)
  data.table::setnames(ep, c("s", "e"), c("start", "end"))
  ep[, `:=`(ndays = as.integer(end - start) + 1L,
            win_end = end + cfg$tail_days, year = data.table::year(start))]
  ep[, episode := sprintf("%s_%s", URAU_CODE, format(start, "%Y%m%d"))]
  data.table::setcolorder(ep, c("episode", "URAU_CODE", "year", "start", "end", "win_end"))
  ep[]
}

#' Add annual harmonic terms
#'
#' @param d data.table with a `date` column (modified in place).
#' @param K Number of sine/cosine pairs.
#' @return `d`, invisibly.
add_harmonics <- function(d, K) {
  doy <- as.numeric(format(d$date, "%j"))
  for (k in seq_len(K)) {
    data.table::set(d, j = paste0("s", k), value = sin(2 * pi * k * doy / 365.25))
    data.table::set(d, j = paste0("c", k), value = cos(2 * pi * k * doy / 365.25))
  }
  invisible(d)
}

#' Reference days for the counterfactual baseline
#'
#' Days used to estimate "deaths expected without heat": outside the
#' influenza season and not within `baseline_hot_lag` days after a day at or
#' above `mmt`. This follows EuroMOMO-type baselines (fitted away from the
#' winter and summer peaks), refined with the city's own temperatures so that
#' cool summer days also inform the summer level.
#'
#' @param date Dates.
#' @param temp Temperatures.
#' @param mmt Temperature above which heat is considered.
#' @param cfg Configuration list.
#' @return Logical vector.
reference_days <- function(date, temp, mmt, cfg) {
  hot <- temp >= mmt
  recent <- do.call(pmax, c(lapply(0:cfg$baseline_hot_lag, function(k)
    as.integer(data.table::shift(hot, k, fill = FALSE))), na.rm = TRUE)) > 0
  !(data.table::month(date) %in% cfg$baseline_exclude_months) & !recent
}

#' Fit a counterfactual baseline model
#'
#' Quasi-Poisson regression of daily deaths on day of week, a slow trend
#' (natural spline of time) and annual harmonics, fitted on reference days
#' only (see [reference_days()]). Reference days are cool, so their deaths
#' include cold-related deaths; with `baseline_cold_offset` the published
#' log-RR (cold side on those days) enters as an offset, so that the model
#' describes deaths at the MMT. Predicted on every day without the offset,
#' it gives the deaths expected with neither heat nor cold (the reference of
#' the attributable-deaths formula), including the lower mortality of summer.
#'
#' @param d data.table with `date`, `y`, `ref`, `dow`, `time` and harmonics.
#' @param cfg Configuration list.
#' @return The fitted glm, or a character error message.
fit_baseline_model <- function(d, cfg) {
  hk <- paste0(c("s", "c"), rep(seq_len(cfg$baseline_harmonics), each = 2))
  df <- max(1, round(length(unique(data.table::year(d$date))) * cfg$baseline_trend_df_per_year))
  f <- stats::as.formula(paste("y ~ dow + splines::ns(time, df =", df, ") +",
                               paste(hk, collapse = " + "),
                               if (isTRUE(cfg$baseline_cold_offset)) "+ offset(off)" else ""))
  m <- tryCatch(suppressWarnings(stats::glm(f, stats::quasipoisson, data = d[ref == TRUE])),
                error = function(e) conditionMessage(e))
  if (!is.character(m) && !m$converged) m <- "baseline model did not converge"
  m
}

#' Draws of daily baseline deaths
#'
#' @param m Fitted baseline glm.
#' @param newdata data.table of the days to predict.
#' @param nsim Number of draws.
#' @return List with `point` (vector) and `draws` (`nsim x nrow(newdata)`).
baseline_draws <- function(m, newdata, nsim) {
  tt <- stats::delete.response(stats::terms(m))
  X <- stats::model.matrix(tt, stats::model.frame(tt, newdata, xlev = m$xlevels))
  b <- stats::coef(m); ok <- !is.na(b)
  bd <- rmvn(nsim, b[ok], stats::vcov(m)[ok, ok, drop = FALSE])
  list(point = as.numeric(exp(X[, ok, drop = FALSE] %*% b[ok])),
       draws = exp(bd %*% t(X[, ok, drop = FALSE])))
}

#' Counterfactual comparison of heat deaths
#'
#' For every city and age group, estimates B_t, the deaths expected on day t
#' without heat ([fit_baseline_model()]), and compares over each window:
#' \describe{
#'   \item{observed heat deaths}{X = O - B, observed minus expected without heat;}
#'   \item{predicted heat deaths}{P = sum B_t (RR_t - 1) over the days at or
#'     above the MMT, with the published RR centred on the MMT: the quantity
#'     the app reports, computed with a fair, seasonally varying baseline.
#'     Days below the MMT inside a window carry the cold side of the curve
#'     and are reported separately as `P_cold`.}
#' }
#' Windows are heatwave episodes (plus `tail_days`), the `post_days` after
#' them (displacement) and whole summers (June-September). Baseline and ERF
#' uncertainty are propagated with draws. Age groups with fewer than
#' `baseline_min_daily` deaths per day take the 20+ baseline times their
#' share of deaths on reference days. Cities whose 20+ model fails are
#' recorded in `failures`.
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo Exposure data.table.
#' @param episodes Output of [detect_episodes()].
#' @param erf Output of [read_erf_bundle()].
#' @param start,end Period.
#' @param cfg Configuration list.
#' @return List with `windows`, `draws_P`, `draws_E` (windows x nsim),
#'   `daily` and `failures`.
fit_counterfactual <- function(mort, expo, episodes, erf, start, end, cfg) {
  set.seed(cfg$seed)
  dates <- seq(start, end, by = "day")
  yrs <- unique(data.table::year(dates))
  win_rows <- list(); dP <- list(); dE <- list(); daily <- list(); failures <- list()
  for (cd in intersect(unique(mort$URAU_CODE), unique(expo$URAU_CODE))) {
    tt <- merge(data.table::data.table(date = dates), expo[URAU_CODE == cd, .(date, tmean)],
                by = "date", all.x = TRUE)
    if (anyNA(tt$tmean)) {
      failures[[length(failures) + 1]] <- data.table::data.table(
        URAU_CODE = cd, agegroup = "all", reason = "missing temperatures")
      next
    }
    ep <- episodes[URAU_CODE == cd]
    wins <- data.table::rbindlist(list(
      if (nrow(ep)) ep[, .(type = "episode", id = episode, year, start, end = win_end)],
      if (nrow(ep)) ep[, .(type = "post", id = episode, year, start = win_end + 1,
                           end = pmin(win_end + cfg$post_days, max(dates)))],
      data.table::data.table(type = "summer", id = paste0(cd, "_", yrs), year = yrs,
                             start = as.Date(sprintf("%d-06-01", yrs)),
                             end = pmin(as.Date(sprintf("%d-09-30", yrs)), max(dates)))))
    pend <- max(dates)
    wins <- wins[start >= min(dates) & start <= pend]
    wins[, URAU_CODE := cd]
    pls <- lapply(setNames(AGE_GROUPS, AGE_GROUPS), function(g)
      published_logrr(erf, cd, g, tt$tmean, cfg))
    # (named `ya`, not `y`: inside d[...] a column `y` would mask it)
    ya <- data.table::dcast(mort[URAU_CODE == cd & date >= min(dates) & date <= pend],
                            date ~ agegroup, value.var = "deaths")
    ya <- merge(tt[, .(date)], ya, by = "date", all.x = TRUE)
    for (g in AGE_GROUPS) if (!g %in% names(ya)) ya[, (g) := 0L]
    ymat <- as.matrix(ya[, AGE_GROUPS, with = FALSE]); ymat[is.na(ymat)] <- 0
    d_all <- data.table::copy(tt)
    d_all[, `:=`(y = rowSums(ymat), dow = factor(weekdays(date)), time = as.numeric(date))]
    add_harmonics(d_all, cfg$baseline_harmonics)
    d_all[, ref := reference_days(date, tmean, min(vapply(pls, `[[`, 0, "mmt")), cfg)]
    sh <- colSums(ymat) / max(sum(ymat), 1)
    d_all[, off := log(Reduce(`+`, Map(function(p, w) w * exp(p$logrr), pls, sh)))]
    m_all <- fit_baseline_model(d_all, cfg)
    if (is.character(m_all)) {
      failures[[length(failures) + 1]] <- data.table::data.table(
        URAU_CODE = cd, agegroup = "20+", reason = m_all)
      next
    }
    b_all <- baseline_draws(m_all, d_all, cfg$nsim)
    phi_all <- max(summary(m_all)$dispersion, 1)
    z0 <- numeric(nrow(wins))
    acc <- data.table::data.table(O = z0, B = z0, P = z0, P_cold = z0, E = z0, v = z0,
                                  attr_app = z0)
    P_d <- matrix(0, nrow(wins), cfg$nsim); E_d <- P_d
    for (g in AGE_GROUPS) {
      pl <- pls[[g]]
      d <- data.table::copy(d_all)
      d[, y := as.numeric(ymat[, g])]
      d[, `:=`(ref = reference_days(date, tmean, pl$mmt, cfg), off = pl$logrr)]
      own <- mean(d$y) >= cfg$baseline_min_daily
      m <- if (own) fit_baseline_model(d, cfg) else "few deaths"
      if (is.character(m)) {
        share <- sum(d$y[d_all$ref]) / max(sum(d_all$y[d_all$ref]), 1)
        bp <- list(point = share * b_all$point, draws = share * b_all$draws)
        phi <- phi_all
        if (own) failures[[length(failures) + 1]] <- data.table::data.table(
          URAU_CODE = cd, agegroup = g, reason = paste(m, "(20+ model x share used)"))
      } else {
        bp <- baseline_draws(m, d, cfg$nsim); phi <- max(summary(m)$dispersion, 1)
      }
      D <- erf_draws(erf, cd, g, cfg$nsim)
      Bm <- erf_basis(pl$mmt, pl$spec)
      bd_app <- app_baseline(erf, cd, g)
      for (k in seq_len(nrow(wins))) {
        w <- d$date >= wins$start[k] & d$date <= wins$end[k]
        rr <- exp(pl$logrr[w])
        hot <- d$tmean[w] >= pl$mmt   # days below the MMT carry the cold side, not heat
        rrd <- exp(D %*% t(erf_basis(d$tmean[w], pl$spec) - rep(Bm, each = sum(w))))
        bw <- bp$draws[, w, drop = FALSE]
        P_d[k, ] <- P_d[k, ] + rowSums((bw * (rrd - 1))[, hot, drop = FALSE])
        E_d[k, ] <- E_d[k, ] + rowSums(bw * rrd)
        acc[k, `:=`(O = O + sum(d$y[w]), B = B + sum(bp$point[w]),
                    P = P + sum((bp$point[w] * (rr - 1))[hot]),
                    P_cold = P_cold + sum((bp$point[w] * (rr - 1))[!hot]),
                    E = E + sum(bp$point[w] * rr),
                    v = v + phi * sum(bp$point[w] * rr),
                    attr_app = attr_app + sum(((1 - 1 / rr) * bd_app)[pl$logrr[w] > 0]))]
      }
      daily[[length(daily) + 1]] <- d[, .(URAU_CODE = cd, agegroup = g, date, deaths = y, tmean,
                                          B = bp$point, logrr = pl$logrr, mmt = pl$mmt, ref)]
    }
    win_rows[[length(win_rows) + 1]] <- cbind(wins, acc)
    dP[[length(dP) + 1]] <- P_d; dE[[length(dE) + 1]] <- E_d
  }
  list(windows = data.table::rbindlist(win_rows), draws_P = do.call(rbind, dP),
       draws_E = do.call(rbind, dE), daily = data.table::rbindlist(daily),
       failures = data.table::rbindlist(failures, fill = TRUE))
}

#' Evaluate predicted against observed heat deaths
#'
#' For each window type (episodes, whole summers): observed heat deaths
#' X = O - B against predicted P = sum B (RR - 1), with the 95% interval of
#' P and a 95% prediction interval of O (baseline and ERF draws + quasi-
#' Poisson noise, simulated as negative binomial). Across windows:
#' calibration slope of X on P (weighted, through the origin; 1 = calibrated,
#' < 1 = the published curves over-predict), coverage, deviance skill of the
#' ERF over the no-heat baseline, and totals. Post-episode windows show
#' displacement.
#'
#' @param cf Output of [fit_counterfactual()].
#' @param cfg Configuration list.
#' @return List with `windows`, `calibration`, `totals`, `by_year`.
evaluate_counterfactual <- function(cf, cfg) {
  set.seed(cfg$seed)
  w <- data.table::copy(cf$windows)
  if (!nrow(w)) return(NULL)
  rnb <- function(mu, phi) {
    phi <- pmax(phi, 1 + 1e-6); mu <- pmax(mu, 1e-9)
    stats::rnbinom(length(mu), mu = mu, size = mu / (phi - 1))
  }
  sim <- t(vapply(seq_len(nrow(w)), function(k) {
    e <- cf$draws_E[k, ]
    o <- rnb(e, w$v[k] / max(w$E[k], 1e-9))
    c(stats::quantile(cf$draws_P[k, ], c(.025, .975)), stats::quantile(o, c(.025, .975)),
      stats::var(e))
  }, numeric(5)))
  w[, `:=`(X = O - B, P_lo = sim[, 1], P_hi = sim[, 2], pi_lo = sim[, 3], pi_hi = sim[, 4],
           var_E = sim[, 5])]
  w[, `:=`(in_pi = O >= pi_lo & O <= pi_hi, z = (O - E) / sqrt(v + var_E))]
  dev <- function(o, e) 2 * (ifelse(o > 0, o * log(o / e), 0) - (o - e))
  one <- function(x, lab) {
    if (nrow(x) < 3) return(NULL)
    c0 <- stats::lm(X ~ 0 + P, data = x, weights = 1 / (x$v + x$var_E))
    ci <- stats::confint(c0)
    data.table::data.table(windows = lab, n = nrow(x), slope = stats::coef(c0)[1],
                           lo = ci[1, 1], hi = ci[1, 2], r = stats::cor(x$X, x$P),
                           observed = sum(x$X), predicted = sum(x$P),
                           ratio_total = sum(x$X) / sum(x$P),
                           deviance_skill = 1 - sum(dev(x$O, x$E)) / sum(dev(x$O, x$B)),
                           coverage_95 = mean(x$in_pi), mean_z = mean(x$z), sd_z = stats::sd(x$z))
  }
  cal <- data.table::rbindlist(list(one(w[type == "episode"], "Heatwave episodes"),
                                    one(w[type == "summer"], "Whole summers (June-September)")))
  ie <- which(w$type == "episode")
  totals <- data.table::data.table(
    observed_excess = sum(w$X[ie]), predicted = sum(w$P[ie]),
    predicted_lo = stats::quantile(colSums(cf$draws_P[ie, , drop = FALSE]), .025),
    predicted_hi = stats::quantile(colSums(cf$draws_P[ie, , drop = FALSE]), .975),
    attributable_app = sum(w$attr_app[ie]),
    post_observed = sum(w[type == "post", X]), post_predicted = sum(w[type == "post", P]))
  by_year <- w[type %in% c("episode", "summer"),
               .(windows = .N, cities = data.table::uniqueN(URAU_CODE), O = sum(O), B = sum(B),
                 X = sum(X), P = sum(P), coverage_95 = mean(in_pi)), by = .(type, year)]
  list(windows = w, calibration = cal, totals = totals, by_year = by_year[order(type, year)])
}

#' Propagate temperature forecast errors to episode predictions
#'
#' For each episode within the forecast archive and each lead time L, the
#' daily temperatures of the episode window are replaced by the forecasts
#' issued L days earlier (with the same exposure alignment as ERA5-Land) and
#' the predicted heat deaths P_L = sum B (RR(T_L) - 1) are recomputed with the
#' same baselines. Comparing P_L with P isolates the forecast error.
#'
#' @param cf Output of [fit_counterfactual()].
#' @param ev Output of [evaluate_counterfactual()].
#' @param fc Output of [fetch_forecasts_all()].
#' @param era5 Raw ERA5-Land series, for the debiasing option.
#' @param check Output of [compare_exposure()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `episodes` (episode x lead) and `summary` (by lead).
forecast_episode_impact <- function(cf, ev, fc, era5, check, erf, cfg) {
  ep <- ev$windows[type == "episode" & start >= cfg$forecast_start + max(cfg$forecast_leads)]
  if (!nrow(ep)) return(NULL)
  f <- data.table::copy(fc[lead %in% cfg$forecast_leads])
  f[, month := data.table::month(date)]
  if (isTRUE(cfg$bias_correct)) {
    f <- merge(f, check$delta, by = c("URAU_CODE", "month"), all.x = TRUE)
    f[is.na(delta), delta := 0][, tmean_fc := tmean_fc - delta][, delta := NULL]
  }
  if (isTRUE(cfg$fc_debias)) {
    b <- merge(fc, era5, by = c("URAU_CODE", "date"))
    # estimated away from the episodes being scored, so the correction cannot
    # borrow information from the days it is later applied to
    inep <- data.table::rbindlist(lapply(seq_len(nrow(ep)), function(k)
      data.table::data.table(URAU_CODE = ep$URAU_CODE[k],
                             date = seq(ep$start[k], ep$end[k], by = "day"))))
    b <- b[!inep, on = .(URAU_CODE, date)]
    b <- b[, .(fbias = mean(tmean_fc - tmean_raw)),
           by = .(URAU_CODE, month = data.table::month(date), lead)]
    f <- merge(f, b, by = c("URAU_CODE", "month", "lead"), all.x = TRUE)
    f[is.na(fbias), fbias := 0][, tmean_fc := tmean_fc - fbias]
  }
  res <- list()
  for (k in seq_len(nrow(ep))) {
    e <- ep[k]
    dd <- cf$daily[URAU_CODE == e$URAU_CODE & date >= e$start & date <= e$end]
    for (L in cfg$forecast_leads) {
      tf <- f[URAU_CODE == e$URAU_CODE & lead == L & date >= e$start & date <= e$end,
              .(date, tmean_fc)]
      if (nrow(tf) < as.integer(e$end - e$start) + 1) next
      x <- merge(dd, tf, by = "date")
      p <- 0
      for (g in unique(x$agegroup)) {
        xg <- x[agegroup == g]
        lr <- published_logrr(erf, e$URAU_CODE, g, xg$tmean_fc, cfg)$logrr
        p <- p + sum(xg$B * (exp(lr) - 1))
      }
      res[[length(res) + 1]] <- data.table::data.table(
        episode = e$id, URAU_CODE = e$URAU_CODE, lead = L, P_era5 = e$P, P_fc = p,
        X = e$X, t_err = mean(x[agegroup == x$agegroup[1], tmean_fc - tmean]))
    }
  }
  r <- data.table::rbindlist(res)
  if (!nrow(r)) return(NULL)
  r[, rel_err := (P_fc - P_era5) / pmax(abs(P_era5), 1)]
  summ <- r[, .(episodes = .N, ratio_total = sum(P_fc) / sum(P_era5),
                mae_deaths = mean(abs(P_fc - P_era5)),
                median_abs_rel_err = stats::median(abs(rel_err)),
                mean_temp_error = mean(t_err)), by = lead]
  list(episodes = r, summary = summ[order(lead)])
}
