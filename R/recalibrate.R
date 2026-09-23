#' Recalibrate the curves on a past window
#'
#' Refits Masselot's first-stage model on all deaths aged 20+ in a past
#' window, using his own ERA5-Land series, and pools the cities in a
#' multivariate random-effects meta-analysis, keeping each city's best linear
#' unbiased prediction (BLUP). These are "Italian curves, period X": the same
#' shape family as the published ones, with a level estimated from Italian
#' data of that period. Fitting them on a window that ends before the test
#' period keeps the later comparison out of sample.
#'
#' @param mort Output of [aggregate_mortality()].
#' @param expo Exposure series covering the window.
#' @param erf Output of [read_erf_bundle()].
#' @param window Character vector of two dates.
#' @param cfg Configuration list.
#' @return List with `label`, `curves` (per city: `beta`, `vcov`) and
#'   `results` (the underlying city fits), or `NULL` if too few cities.
recalibrate_curves <- function(mort, expo, erf, window, cfg) {
  s0 <- as.Date(window[1]); e0 <- as.Date(window[2])
  tst <- run_curve_tests(mort, expo, erf, s0, e0, cfg)
  f <- Filter(function(x) x$level == "all20", tst$fits)
  if (length(f) < 3) return(NULL)
  Y <- do.call(rbind, lapply(f, `[[`, "coef"))
  S <- lapply(f, `[[`, "vcov")
  bl <- mixmeta::blup(mixmeta::mixmeta(Y, S, method = "reml"), vcov = TRUE)
  curves <- stats::setNames(lapply(seq_along(f), function(i)
    list(beta = bl[[i]]$blup, vcov = bl[[i]]$vcov)), vapply(f, `[[`, "", "URAU_CODE"))
  list(label = sprintf("Italian curves %s-%s", format(s0, "%Y"), format(e0, "%Y")),
       window = window, curves = curves, results = tst$results)
}

#' Daily deaths expected under a set of curves
#'
#' Applies a city's curve to its daily temperatures and the baseline of
#' [fit_counterfactual()] (summed over age groups), giving the deaths
#' expected with heat, `E = B x RR`, and the heat deaths `P = B (RR - 1)` on
#' days at or above the curve's own MMT.
#'
#' @param daily `cf$daily` from [fit_counterfactual()].
#' @param recal Output of [recalibrate_curves()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return data.table `URAU_CODE`, `date`, `B`, `rr`, `heat`.
daily_under_curves <- function(daily, recal, erf, cfg) {
  b <- daily[, .(B = sum(B), tmean = tmean[1]), by = .(URAU_CODE, date)]
  out <- lapply(intersect(unique(b$URAU_CODE), names(recal$curves)), function(cd) {
    sp <- erf_spec(erf, cd, cfg)
    beta <- recal$curves[[cd]]$beta
    m <- erf_mmt(sp, beta)
    x <- b[URAU_CODE == cd]
    lr <- as.numeric(erf_basis(x$tmean, sp) %*% beta) -
      as.numeric(erf_basis(m, sp) %*% beta)
    x[, `:=`(rr = exp(lr), heat = tmean >= m)]
    x
  })
  data.table::rbindlist(out)
}

#' Calibration of one set of predictions
#'
#' @param x data.table with `X` (observed heat deaths), `P` (predicted) and
#'   `v` (variance of the observed count).
#' @param lab Label.
#' @return One-row data.table.
calibration_row <- function(x, lab) {
  if (nrow(x) < 3) return(NULL)
  m <- stats::lm(X ~ 0 + P, data = x, weights = 1 / pmax(x$v, 1e-6))
  ci <- stats::confint(m)
  data.table::data.table(curves = lab, n = nrow(x), observed = sum(x$X),
                         predicted = sum(x$P), ratio = sum(x$X) / sum(x$P),
                         slope = stats::coef(m)[1], lo = ci[1, 1], hi = ci[1, 2],
                         r = stats::cor(x$X, x$P))
}

#' Compare published and recalibrated curves on the test windows
#'
#' Recomputes the predicted heat deaths of every window with the recalibrated
#' curves, keeping the same baseline, and reports the calibration of each set
#' of curves side by side. The recalibrated curves are estimated on a window
#' that ends before the test period, so this is an out-of-sample check.
#'
#' @param cf Output of [fit_counterfactual()].
#' @param ev Output of [evaluate_counterfactual()].
#' @param recals List of [recalibrate_curves()] outputs.
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `windows` (predictions by curve set) and `calibration`.
compare_curve_sets <- function(cf, ev, recals, erf, cfg) {
  w <- ev$windows
  rows <- list(); wins <- list()
  for (ty in c("episode", "summer")) {
    ws <- w[type == ty]
    rows[[length(rows) + 1]] <- cbind(
      windows = ty, calibration_row(ws[, .(X, P, v)], "Published (Masselot et al.)"))
    for (rc in Filter(Negate(is.null), recals)) {
      d <- daily_under_curves(cf$daily, rc, erf, cfg)
      pr <- data.table::rbindlist(lapply(seq_len(nrow(ws)), function(k) {
        x <- d[URAU_CODE == ws$URAU_CODE[k] & date >= ws$start[k] & date <= ws$end[k]]
        data.table::data.table(id = ws$id[k], type = ty, URAU_CODE = ws$URAU_CODE[k],
                               year = ws$year[k], X = ws$X[k], v = ws$v[k],
                               P = sum(x$B[x$heat] * (x$rr[x$heat] - 1)))
      }))
      if (!nrow(pr)) next
      wins[[length(wins) + 1]] <- cbind(curves = rc$label, pr)
      rows[[length(rows) + 1]] <- cbind(windows = ty, calibration_row(pr[, .(X, P, v)], rc$label))
    }
  }
  list(calibration = data.table::rbindlist(rows, fill = TRUE),
       windows = data.table::rbindlist(wins, fill = TRUE))
}

#' Daily series of one city for the example figure
#'
#' Observed deaths, deaths expected without heat, and deaths expected under
#' the published and the recalibrated curves, for chosen summers.
#'
#' @param cf Output of [fit_counterfactual()].
#' @param recals List of [recalibrate_curves()] outputs.
#' @param erf Output of [read_erf_bundle()].
#' @param cities City table (for the name).
#' @param cfg Configuration list.
#' @return List with `data` (long data.table) and `city` (name).
example_series <- function(cf, recals, erf, cities, cfg) {
  cd <- cities[norm_name(name) == norm_name(cfg$example_city), URAU_CODE]
  if (!length(cd)) {
    tot <- cf$daily[, .(d = sum(deaths)), by = URAU_CODE][order(-d)]
    cd <- tot$URAU_CODE[1]
  }
  cd <- cd[1]
  d <- cf$daily[URAU_CODE == cd & data.table::month(date) %in% 6:9 &
                  data.table::year(date) %in% cfg$example_years]
  if (!nrow(d)) return(NULL)
  base <- d[, .(Observed = sum(deaths), `Expected without heat` = sum(B),
                `Published curves` = sum(B * exp(logrr))), by = date]
  for (rc in Filter(Negate(is.null), recals)) {
    dd <- daily_under_curves(cf$daily[URAU_CODE == cd], rc, erf, cfg)
    base <- merge(base, dd[, .(date, v = B * rr)], by = "date", all.x = TRUE)
    data.table::setnames(base, "v", rc$label)
  }
  long <- data.table::melt(base, id.vars = "date", variable.name = "series",
                           value.name = "deaths")
  long[, year := data.table::year(date)]
  long[, doy := as.numeric(format(date, "%j"))]
  list(data = long, city = cities[URAU_CODE == cd, name][1])
}
