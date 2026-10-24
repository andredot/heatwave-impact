# =============================================================================
#  Audit of the baseline: is the yardstick of "observed heat deaths" sound?
#
#  Observed heat deaths are X = O - B, so any error in the level of B, the
#  deaths expected without heat, moves X directly: a baseline 3% too high
#  removes 3% of all summer deaths from X, which can be most of the heat
#  signal. Two checks:
#    1. cool days: on the reference days (not hot, not just after heat) the
#       model expects B x the cold relative risk; observed / expected there
#       should be 1. Summer cool days are few and sit between hot spells, so a
#       baseline whose seasonal shape is wrong in summer shows up here, while
#       spring and autumn days, which dominate the fit, stay close to 1;
#    2. sensitivity: the baseline is refitted with alternative, equally
#       defensible choices and the whole evaluation is repeated, to see how
#       much of the calibration ratio depends on the yardstick.
# =============================================================================

#' Observed against expected deaths on the baseline's reference days
#'
#' @param cf Output of [fit_counterfactual()].
#' @param label Period label.
#' @return data.table by city and season: days, `O`, `E` (= B x cold RR, what
#'   the model expects on those days), `V` (variance of O), `ratio`, `lo`, `hi`.
cool_day_audit <- function(cf, label) {
  if (is.null(cf) || !nrow(cf$daily)) return(NULL)
  d <- data.table::copy(cf$daily[ref == TRUE])
  if (!"phi" %in% names(d)) d[, phi := 1]
  m <- data.table::month(d$date)
  d[, season := data.table::fcase(m %in% 6:9, "Summer (Jun-Sep)",
                                  m %in% c(4, 5, 10, 11), "Spring and autumn",
                                  default = NA_character_)]
  d <- d[!is.na(season)]
  d[, E := B * exp(logrr)]
  out <- d[, .(days = data.table::uniqueN(date), O = sum(deaths), E = sum(E),
               V = sum(phi * E)), by = .(URAU_CODE, season)]
  out[, `:=`(period = label, ratio = O / E, lo = (O - 1.96 * sqrt(V)) / E,
             hi = (O + 1.96 * sqrt(V)) / E)]
  out[]
}

#' Pool the cool-day audit over cities
#'
#' @param a Output of [cool_day_audit()] (several periods / baselines bound).
#' @param by Grouping columns.
#' @return data.table with pooled ratio and 95% interval.
pool_audit <- function(a, by) {
  a[, .(cities = .N, days = sum(days), O = sum(O), E = sum(E), V = sum(V)), by = by][
    , `:=`(ratio = O / E, lo = (O - 1.96 * sqrt(V)) / E, hi = (O + 1.96 * sqrt(V)) / E)][]
}

#' Summary of one evaluation, for the sensitivity table
#'
#' @param re Output of [rolling_evaluation()].
#' @param label Baseline label.
#' @return data.table, one row per curve set.
summarise_evaluation <- function(re, label) {
  if (is.null(re)) return(NULL)
  ap <- tail(re$periods, 1)
  cal <- re$calibration[period == ap & type == "episode" & subset == "All windows",
                        .(curves, episodes = n, observed, predicted, ratio, slope, coverage_95)]
  w <- re$warnings[period == ap, .(curves, events = hits + misses, pod, far, hss)]
  sm <- re$windows[type == "summer" & curves == re$curve_sets[1],
                   .(summer_O = sum(O), summer_B = sum(B))]
  out <- merge(cal, w, by = "curves")
  out[, `:=`(baseline = label, summer_excess_pct = 100 * (sm$summer_O / sm$summer_B - 1))]
  data.table::setcolorder(out, "baseline")
  out[]
}

#' Baseline audit: cool days and sensitivity of the evaluation
#'
#' The configured baseline is summarised from the pipeline's own objects; each
#' variant of `cfg$baseline_variants` (a list of configuration overrides)
#' refits the counterfactual of both periods for the evaluated cities and
#' repeats [rolling_evaluation()] with the same curves.
#'
#' @param variants Named list of configuration overrides.
#' @param mort Output of [aggregate_mortality()].
#' @param expo_hist,expo_recent Exposure series of the two periods.
#' @param ep_hist,ep_test Episodes of the two periods.
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @param codes Cities evaluated.
#' @param rolling The system's curve set.
#' @param periods Output of [make_periods()].
#' @param re_default Output of [rolling_evaluation()] with the configured baseline.
#' @param cf_hist,cf_test Counterfactuals of the configured baseline.
#' @return List with `cool_days` (city x season x period x baseline),
#'   `cool_pooled` and `sensitivity` (one row per baseline x curve set).
baseline_sensitivity <- function(variants, mort, expo_hist, expo_recent, ep_hist, ep_test,
                                 erf, cfg, codes, rolling, periods, re_default,
                                 cf_hist, cf_test) {
  lab_h <- sprintf("%s-%s", format(cfg$repro_start, "%Y"), format(cfg$repro_end, "%Y"))
  lab_t <- sprintf("%s-%s", format(periods$test$start, "%Y"), format(periods$test$end, "%Y"))
  def <- "As configured"
  aud <- list(cbind(baseline = def, data.table::rbindlist(list(
    cool_day_audit(list(daily = cf_hist$daily[URAU_CODE %in% codes]), lab_h),
    cool_day_audit(list(daily = cf_test$daily[URAU_CODE %in% codes]), lab_t)))))
  sens <- list(summarise_evaluation(re_default, def))
  m <- mort[URAU_CODE %in% codes]
  for (nm in names(variants)) {
    c2 <- utils::modifyList(cfg, variants[[nm]])
    cfh <- fit_counterfactual(m, expo_hist[URAU_CODE %in% codes], ep_hist[URAU_CODE %in% codes],
                              erf, cfg$repro_start, cfg$repro_end, c2)
    cft <- fit_counterfactual(m, expo_recent[URAU_CODE %in% codes], ep_test[URAU_CODE %in% codes],
                              erf, periods$test$start, periods$test$end, c2)
    re <- rolling_evaluation(
      list(list(label = lab_h, cf = cfh, ev = evaluate_counterfactual(cfh, c2)),
           list(label = lab_t, cf = cft, ev = evaluate_counterfactual(cft, c2))),
      rolling, list(), erf, c2, codes)
    sens[[length(sens) + 1]] <- summarise_evaluation(re, nm)
    aud[[length(aud) + 1]] <- cbind(baseline = nm, data.table::rbindlist(list(
      cool_day_audit(cfh, lab_h), cool_day_audit(cft, lab_t))))
  }
  a <- data.table::rbindlist(aud, fill = TRUE)
  lev <- c(def, names(variants))
  a[, baseline := factor(baseline, lev)]
  s <- data.table::rbindlist(sens, fill = TRUE)
  s[, baseline := factor(baseline, lev)]
  list(cool_days = a, cool_pooled = pool_audit(a, c("baseline", "period", "season")),
       sensitivity = s[order(baseline, curves)])
}
