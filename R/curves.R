#' Cumulative age binning (Masselot's rule)
#'
#' Reproduces `MESS::cumsumbinning(x, threshold, cutwhenpassed = TRUE)` as
#' used in Masselot's first stage: age groups are accumulated from the
#' youngest until their total deaths pass the threshold, then a new group
#' starts. The last group may stay below the threshold.
#'
#' @param x Deaths per age group, ordered by age.
#' @param threshold Minimum deaths per group.
#' @return Integer group index for each element of `x`.
cumsum_binning <- function(x, threshold) {
  g <- integer(length(x)); k <- 1L; acc <- 0
  for (i in seq_along(x)) {
    g[i] <- k; acc <- acc + x[i]
    if (acc >= threshold && i < length(x)) { k <- k + 1L; acc <- 0 }
  }
  g
}

#' Label of a merged age group
#'
#' @param ages Consecutive Masselot age groups.
#' @return Label such as `"20-64"` or `"65+"`.
merged_label <- function(ages) {
  lo <- sub("-.*|\\+", "", ages[1])
  last <- ages[length(ages)]
  if (grepl("\\+$", last)) paste0(lo, "+") else paste0(lo, "-", sub(".*-", "", last))
}

#' Fit Masselot's first-stage model on one series
#'
#' Quasi-Poisson regression with a cross-basis of temperature (quadratic
#' B-spline, knots at the 10/75/90th percentiles and boundaries at the range
#' of the city's 1990-2019 distribution, i.e. the knots of the published
#' ERF) and lag (natural spline, 3 knots equally spaced on the log scale up
#' to lag 21), day-of-week indicators and a natural spline of time with 7 df
#' per year. The fit is reduced to the overall cumulative curve, whose 5
#' coefficients live in the same basis as the published ones.
#'
#' @param dat data.table with `date`, `y`, `tmean`, starting `lag_max` days
#'   before the period (for the lags).
#' @param start,end Period over which the model is fitted.
#' @param spec Output of [erf_spec()].
#' @param cfg Configuration list.
#' @return List with `coef`, `vcov`, `converged`, `dispersion`, `totdeath`,
#'   or `NULL` if the model fails.
fit_stage1 <- function(dat, start, end, spec, cfg) {
  cb <- suppressWarnings(dlnm::crossbasis(dat$tmean, lag = cfg$lag_max,
                         argvar = list(fun = "bs", degree = spec$degree, knots = spec$knots,
                                       Boundary.knots = spec$bound),
                         arglag = list(fun = "ns",
                                       knots = dlnm::logknots(cfg$lag_max, cfg$lag_nknots))))
  dat <- data.table::copy(dat)
  dat[, `:=`(dow = factor(weekdays(date)), time = as.numeric(date))]
  inp <- dat$date >= start & dat$date <= end
  nyr <- length(unique(data.table::year(dat$date[inp])))
  fit <- tryCatch(suppressWarnings(stats::glm(
    y ~ cb + dow + splines::ns(time, df = cfg$df_time_per_year * nyr),
    data = dat, family = stats::quasipoisson, subset = inp)), error = function(e) NULL)
  if (is.null(fit) || !fit$converged) return(NULL)
  red <- suppressWarnings(dlnm::crossreduce(cb, fit, cen = stats::median(dat$tmean[inp])))
  list(coef = stats::coef(red), vcov = stats::vcov(red), converged = TRUE,
       dispersion = summary(fit)$dispersion, totdeath = sum(dat$y[inp]))
}

#' Compare refitted and published curves in every city
#'
#' For each city, fits the first-stage model (a) on all deaths aged 20+ and
#' (b) on age groups merged with Masselot's cumulative 5000-death rule, over
#' the period, and compares each fit with the published curve of the same
#' ages (death-weighted when groups are merged). Returned metrics:
#' MMT, log-RR at the 99th percentile vs the published MMT (common
#' reference), heat-attributable fraction of the period's deaths, their
#' differences with standard errors from independent draws, and, for single
#' Masselot age groups, a Wald test of equality of the 5 coefficients.
#'
#' @param mort Output of [aggregate_mortality()] (must cover `start`).
#' @param expo data.table `URAU_CODE`, `date`, `tmean` (must start
#'   `lag_max` days before `start`).
#' @param erf Output of [read_erf_bundle()].
#' @param start,end Period.
#' @param cfg Configuration list.
#' @return List with `results` (one row per city x series) and `bands`
#'   (curves on a temperature grid, for plots).
run_curve_tests <- function(mort, expo, erf, start, end, cfg) {
  set.seed(cfg$seed)
  res <- list(); bands <- list(); fits <- list()
  codes <- intersect(unique(mort$URAU_CODE), unique(expo$URAU_CODE))
  for (cd in codes) {
    spec <- erf_spec(erf, cd, cfg)
    ex <- expo[URAU_CODE == cd & date >= start - cfg$lag_max & date <= end, .(date, tmean)]
    m  <- mort[URAU_CODE == cd & date >= start & date <= end]
    tot <- m[, .(d = sum(deaths)), by = agegroup][match(AGE_GROUPS, agegroup), d]
    tot[is.na(tot)] <- 0
    series <- list(list(level = "all20", ages = AGE_GROUPS))
    if (sum(tot) >= cfg$min_deaths) {
      bins <- cumsum_binning(tot, cfg$min_deaths)
      if (length(unique(bins)) > 1)
        series <- c(series, lapply(unique(bins), function(b)
          list(level = "age", ages = AGE_GROUPS[bins == b])))
    }
    for (s in series) {
      lab <- merged_label(s$ages)
      y <- m[agegroup %in% s$ages, .(y = sum(deaths)), by = date]
      dat <- merge(ex, y, by = "date", all.x = TRUE)
      td <- sum(y$y)
      row <- data.table::data.table(URAU_CODE = cd, level = s$level, series = lab,
                                    n_ages = length(s$ages), deaths = td,
                                    below_threshold = td < cfg$min_deaths)
      if (td < cfg$min_deaths && s$level == "all20") {
        res[[length(res) + 1]] <- cbind(row, status = "not fitted (deaths below threshold)")
        next
      }
      fit <- fit_stage1(dat, start, end, spec, cfg)
      if (is.null(fit)) {
        res[[length(res) + 1]] <- cbind(row, status = "not converged")
        next
      }
      w <- tot[match(s$ages, AGE_GROUPS)]
      pub <- curve_published(erf, cd, s$ages, pmax(w, 1), spec, cfg$nsim)
      new <- curve_single(spec, fit$coef, rmvn(cfg$nsim, fit$coef, fit$vcov))
      per <- dat[date >= start & date <= end & !is.na(y)]
      mp <- curve_metrics(pub, per$tmean, per$y, ref = NA)
      mp <- curve_metrics(pub, per$tmean, per$y, ref = mp$mmt)
      mn <- curve_metrics(new, per$tmean, per$y, ref = mp$mmt)
      d_lr <- mn$lr99_ref - mp$lr99_ref
      se_lr <- sqrt(stats::var(mn$lr99_ref_d) + stats::var(mp$lr99_ref_d))
      d_af <- mn$af_heat - mp$af_heat
      se_af <- sqrt(stats::var(mn$af_heat_d) + stats::var(mp$af_heat_d))
      wald <- wald_p <- NA_real_
      if (length(s$ages) == 1) {
        db <- fit$coef - erf_coef(erf, cd, s$ages)
        V  <- fit$vcov + erf_vcov(erf, cd, s$ages)
        wald <- tryCatch(as.numeric(t(db) %*% solve(V, db)), error = function(e) NA_real_)
        wald_p <- stats::pchisq(wald, df = length(db), lower.tail = FALSE)
      }
      q <- function(v) stats::quantile(v, c(.025, .975))
      res[[length(res) + 1]] <- cbind(row, status = "fitted",
        dispersion = fit$dispersion,
        mmt_pub = mp$mmt, mmt_new = mn$mmt, p99 = spec$p99,
        rr99_pub = exp(mp$lr99_ref), rr99_pub_lo = exp(q(mp$lr99_ref_d))[1],
        rr99_pub_hi = exp(q(mp$lr99_ref_d))[2],
        rr99_new = exp(mn$lr99_ref), rr99_new_lo = exp(q(mn$lr99_ref_d))[1],
        rr99_new_hi = exp(q(mn$lr99_ref_d))[2],
        d_lr99 = d_lr, se_d_lr99 = se_lr,
        af_pub = mp$af_heat, af_new = mn$af_heat, d_af = d_af, se_d_af = se_af,
        wald = wald, wald_p = wald_p)
      fits[[length(fits) + 1]] <- list(URAU_CODE = cd, level = s$level, series = lab,
                                       ages = s$ages, w = pmax(w, 1), coef = fit$coef,
                                       vcov = fit$vcov,
                                       trange = stats::quantile(per$tmean, c(.005, .995)))
      grid <- curve_grid(spec, per$tmean)
      bands[[length(bands) + 1]] <- rbind(
        cbind(URAU_CODE = cd, level = s$level, series = lab, curve = "Published",
              curve_band(pub, grid)),
        cbind(URAU_CODE = cd, level = s$level, series = lab, curve = "Refitted",
              curve_band(new, grid)))
    }
  }
  results <- data.table::rbindlist(res, fill = TRUE)
  results <- merge(results, erf$metadata[, .(URAU_CODE, name, inmcc)], by = "URAU_CODE")
  list(results = results, bands = data.table::rbindlist(bands), fits = fits)
}

#' Temperature grid for curve plots
#'
#' Curves are drawn only where the period has data (0.5th-99.5th percentile
#' of its daily temperatures) and within the 1990-2019 range the basis is
#' defined on. Outside it a refitted curve is pure extrapolation: with three
#' mild winters, the cold end of the basis (between the 1990-2019 minimum and
#' its 10th percentile) is almost unsupported, which produced the spurious
#' cold spikes.
#'
#' @param spec Output of [erf_spec()].
#' @param temps Daily temperatures of the period.
#' @return Numeric grid.
curve_grid <- function(spec, temps) {
  lo <- max(spec$pct(0), stats::quantile(temps, .005))
  hi <- min(spec$pct(100), stats::quantile(temps, .995))
  seq(lo, hi, length.out = 80)
}

#' Second-stage pooling of refitted curves (BLUPs)
#'
#' Short series give noisy city curves. As in Masselot's second stage, the
#' five coefficients of the refitted overall curves (deaths 20+) are pooled
#' in a multivariate random-effects meta-analysis (`mixmeta`, REML) and each
#' city's best linear unbiased prediction (BLUP) is used instead: a weighted
#' compromise between the city's own curve and the average curve, driven by
#' the city's precision. Coefficients are comparable across cities because
#' knots sit at the same percentiles. The pooled model is also used to test
#' whether the cities' curves differ from the published ones on average
#' (meta-regression on the published coefficients is not needed: the
#' comparison is done on the relative-risk scale).
#'
#' @param tests Output of [run_curve_tests()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `results` (BLUP RR99 vs published, by city) and `bands`.
blup_curves <- function(tests, erf, cfg) {
  set.seed(cfg$seed)
  f <- Filter(function(x) x$level == "all20", tests$fits)
  if (length(f) < 3) return(NULL)
  Y <- do.call(rbind, lapply(f, `[[`, "coef"))
  S <- lapply(f, `[[`, "vcov")
  mm <- mixmeta::mixmeta(Y, S, method = "reml")
  bl <- mixmeta::blup(mm, vcov = TRUE)
  res <- list(); bands <- list()
  for (i in seq_along(f)) {
    cd <- f[[i]]$URAU_CODE; spec <- erf_spec(erf, cd, cfg)
    b <- bl[[i]]$blup; V <- bl[[i]]$vcov
    new <- curve_single(spec, b, rmvn(cfg$nsim, b, V))
    pub <- curve_published(erf, cd, f[[i]]$ages, f[[i]]$w, spec, cfg$nsim)
    ref <- spec$mmt_grid[which.min(curve_eval(pub, spec$mmt_grid))]
    pts <- c(spec$p99, ref)
    fn <- curve_eval(new, pts); fnd <- curve_eval(new, pts, draws = TRUE)
    fp <- curve_eval(pub, pts); fpd <- curve_eval(pub, pts, draws = TRUE)
    d <- (fn[1] - fn[2]) - (fp[1] - fp[2])
    se <- sqrt(stats::var(fnd[, 1] - fnd[, 2]) + stats::var(fpd[, 1] - fpd[, 2]))
    res[[i]] <- data.table::data.table(URAU_CODE = cd, rr99_pub = exp(fp[1] - fp[2]),
                                       rr99_blup = exp(fn[1] - fn[2]), d_lr99 = d, se_d_lr99 = se)
    grid <- seq(f[[i]]$trange[1], f[[i]]$trange[2], length.out = 80)
    grid <- grid[grid >= spec$pct(0) & grid <= spec$pct(100)]
    bands[[i]] <- rbind(cbind(URAU_CODE = cd, curve = "Published", curve_band(pub, grid)),
                        cbind(URAU_CODE = cd, curve = "Refitted (BLUP)", curve_band(new, grid)))
  }
  out <- merge(data.table::rbindlist(res), erf$metadata[, .(URAU_CODE, name, inmcc)],
               by = "URAU_CODE")
  list(results = out, bands = data.table::rbindlist(bands),
       heterogeneity = summary(mm)$i2stat)
}

#' Pool the refitted-vs-published differences across cities
#'
#' Random-effects meta-analysis (REML, `mixmeta`) of the city differences in
#' log-RR at the 99th percentile, for all deaths 20+: overall, by whether the
#' city's published curve was informed by its own mortality data (MCC city)
#' or extrapolated, and a meta-regression testing the difference between the
#' two groups. A pooled difference of 0 means the published curves are, on
#' average, neither too steep nor too flat in the new data.
#'
#' @param tests Output of [run_curve_tests()].
#' @return data.table with pooled estimates (log scale and ratio of RRs),
#'   95% CI, p-value, I-squared and number of cities.
pool_curve_differences <- function(tests) {
  r <- tests$results[level == "all20" & status == "fitted" & is.finite(se_d_lr99) & se_d_lr99 > 0]
  pool <- function(d, lab) {
    if (nrow(d) < 2) return(data.table::data.table(group = lab, k = nrow(d)))
    fit <- mixmeta::mixmeta(d_lr99, S = se_d_lr99^2, data = d, method = "reml")
    sm <- summary(fit)
    est <- stats::coef(fit)[1]; se <- sqrt(stats::vcov(fit)[1, 1])
    data.table::data.table(group = lab, k = nrow(d), est = est, se = se,
                           ratio = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                           p = 2 * stats::pnorm(-abs(est / se)),
                           i2 = as.numeric(sm$i2stat[1]))
  }
  out <- list(pool(r, "All cities"))
  if (any(!is.na(r$inmcc))) {
    out <- c(out, list(pool(r[inmcc == TRUE], "Curve fitted on local data (MCC)"),
                       pool(r[inmcc == FALSE], "Curve extrapolated")))
    if (length(unique(r$inmcc)) == 2) {
      mr <- mixmeta::mixmeta(d_lr99 ~ inmcc, S = se_d_lr99^2, data = r, method = "reml")
      z <- stats::coef(mr)[2] / sqrt(stats::vcov(mr)[2, 2])
      out <- c(out, list(data.table::data.table(
        group = "Difference MCC vs extrapolated (meta-regression)", k = nrow(r),
        est = stats::coef(mr)[2], se = sqrt(stats::vcov(mr)[2, 2]),
        ratio = exp(stats::coef(mr)[2]), p = 2 * stats::pnorm(-abs(z)))))
    }
  }
  data.table::rbindlist(out, fill = TRUE)
}


#' Refit the curves in successive time windows
#'
#' Runs [run_curve_tests()] separately in each window of `cfg$period_windows`
#' and pools the city differences with the published curves. Windows ending
#' before 2020 use Masselot's own ERA5-Land series (identical exposure),
#' later ones our aligned series. If heat vulnerability is falling, the
#' pooled ratio should decline monotonically; a meta-regression on the window
#' midpoint tests that trend.
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo_hist Masselot's ERA5-Land series.
#' @param expo_recent Our aligned exposure series.
#' @param erf Output of [read_erf_bundle()].
#' @param data_end Last day with mortality data.
#' @param cfg Configuration list.
#' @return List with `pooled` (one row per window), `cities` (city-level
#'   differences with the window label) and `trend` (meta-regression of the
#'   difference on the window midpoint, per decade).
run_period_tests <- function(mort, expo_hist, expo_recent, erf, data_end, cfg) {
  res <- list(); cit <- list()
  for (w in cfg$period_windows) {
    s0 <- as.Date(w[1]); e0 <- if (is.na(w[2])) data_end else as.Date(w[2])
    if (e0 > data_end) e0 <- data_end
    if (s0 >= e0) next
    expo <- if (e0 <= as.Date("2019-12-31")) expo_hist else expo_recent
    tst <- run_curve_tests(mort, expo, erf, s0, e0, cfg)
    pl <- pool_curve_differences(tst)
    lab <- sprintf("%s-%s", format(s0, "%Y"), format(e0, "%Y"))
    mid <- as.numeric(format(s0, "%Y")) + (as.numeric(format(e0, "%Y")) -
                                             as.numeric(format(s0, "%Y"))) / 2
    res[[length(res) + 1]] <- cbind(window = lab, midyear = mid,
                                    pl[group == "All cities"])
    cit[[length(cit) + 1]] <- cbind(window = lab, midyear = mid,
      tst$results[level == "all20" & status == "fitted",
                  .(URAU_CODE, name, inmcc, deaths, d_lr99, se_d_lr99, mmt_new, mmt_pub)])
  }
  pooled <- data.table::rbindlist(res, fill = TRUE)
  cities <- data.table::rbindlist(cit, fill = TRUE)
  trend <- NULL
  if (nrow(cities) && data.table::uniqueN(cities$window) > 1) {
    d <- cities[is.finite(se_d_lr99) & se_d_lr99 > 0]
    m <- mixmeta::mixmeta(d_lr99 ~ I((midyear - 2010) / 10), S = d$se_d_lr99^2,
                          data = d, random = ~ 1 | URAU_CODE, method = "reml")
    b <- stats::coef(m)[2]; se <- sqrt(stats::vcov(m)[2, 2])
    trend <- data.table::data.table(
      change_per_decade = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
      p = 2 * stats::pnorm(-abs(b / se)))
  }
  list(pooled = pooled, cities = cities, trend = trend)
}
