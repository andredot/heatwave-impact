# =============================================================================
#  Test B repeated with the FluMOMO baseline.
#
#  Test B compares heat deaths predicted by the curves with deaths observed
#  above a baseline of deaths expected without heat. That baseline is ours.
#  Here the same comparison is made with the baseline of the official FluMOMO
#  code, estimated on the same population (the comuni of the validated
#  cities of one region) by a different method and a different group.
#
#  FluMOMO's expected series without temperature but with influenza is
#  EB + EdIA: the baseline plus the deaths it attributes to influenza. Using
#  it leaves influenza out of the residual, so what is left above it on hot
#  weeks is the heat signal to compare with the curves.
#
#  Resolution is weekly and the population is all ages; deaths from heat below
#  age 20 are negligible, which is why Masselot's curves start there.
# =============================================================================

#' Weekly heat deaths predicted by a set of curves
#'
#' @param daily `cf$daily` from [fit_counterfactual()], restricted to the cities.
#' @param rr data.table `URAU_CODE`, `date`, `B`, `rr`, `heat` from
#'   [daily_under_curves()], or `NULL` to use the published curves in `daily`.
#' @return data.table `year`, `week`, `P`.
weekly_predicted_heat <- function(daily, rr = NULL) {
  d <- if (is.null(rr)) {
    daily[tmean >= mmt, .(P = sum(B * (exp(logrr) - 1))), by = date]
  } else {
    rr[heat == TRUE, .(P = sum(B * (rr - 1))), by = date]
  }
  d[, week_start := iso_week_start(date)]
  d[, .(P = sum(P)), by = week_start]
}

#' Test B against the FluMOMO baseline
#'
#' @param res FluMOMO results for the cities' comuni ([run_region_flumomo()]).
#' @param cf Output of [fit_counterfactual()].
#' @param recals List of [recalibrate_curves()] outputs.
#' @param erf Output of [read_erf_bundle()].
#' @param istat_file Istat daily CSV.
#' @param procom Comune codes of the cities used in the FluMOMO run.
#' @param codes Urban Audit codes of those cities.
#' @param cfg Configuration list.
#' @return List with `weekly` (one row per summer week and curve set),
#'   `calibration` (slope and totals per curve set) and `n_cities`.
flumomo_test_b <- function(res, cf, recals, erf, istat_file, procom, codes, cfg) {
  fm <- cfg$flumomo
  obs <- istat_weekly_euromomo(istat_file, fm$provinces$prov, fm$years, procom)
  obs <- obs[agegrp == 4L, .(year, week, observed = deaths)]
  r <- as.data.frame(res)
  r <- r[r$agegrp == fm$agegrp & !is.na(r$year) & !is.na(r$week), ]
  fm_w <- data.table::data.table(year = r$year, week = r$week,
                                 expected = r$EB + r$EdIA, EB = r$EB, EdET = r$EdET)
  base <- merge(obs, fm_w, by = c("year", "week"))
  base[, week_start := ISOweek::ISOweek2date(sprintf("%d-W%02d-1", year, week))]
  base[, `:=`(X = observed - expected,
              month = data.table::month(week_start + 3))]

  daily <- cf$daily[URAU_CODE %in% codes]
  # only the weeks the curve predictions cover: FluMOMO is fitted on a longer
  # history, and weeks before it would enter with a prediction of zero
  span <- range(daily$date)
  base <- base[week_start >= span[1] & week_start + 6 <= span[2]]
  sets <- list(`Published (Masselot et al.)` = weekly_predicted_heat(daily))
  for (rc in Filter(Negate(is.null), recals))
    sets[[rc$label]] <- weekly_predicted_heat(
      daily, daily_under_curves(daily, rc, erf, cfg))

  out <- data.table::rbindlist(lapply(names(sets), function(nm) {
    x <- merge(base, sets[[nm]], by = "week_start", all.x = TRUE)
    x[is.na(P), P := 0][, curves := nm][]
  }))
  warm <- out[month %in% cfg$warm_months]
  cal <- data.table::rbindlist(lapply(split(warm, warm$curves), function(x) {
    m <- stats::lm(X ~ 0 + P, data = x)
    ci <- stats::confint(m)
    data.table::data.table(curves = x$curves[1], weeks = nrow(x),
                           observed = sum(x$X), predicted = sum(x$P),
                           ratio = sum(x$X) / sum(x$P), slope = stats::coef(m)[1],
                           lo = ci[1, 1], hi = ci[1, 2], r = stats::cor(x$X, x$P),
                           flumomo_ET = sum(x$EdET))
  }))
  list(weekly = out, calibration = cal, n_cities = length(codes),
       period = format(span))
}
