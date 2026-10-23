# =============================================================================
#  Synthetic end-to-end check of the rolling recalibration.
#
#  Six fake cities, 2011-2026, with a heat risk that falls over time
#  (adaptation) and a published curve that is too steep. Runs the same
#  functions as the pipeline - windows, fits, BLUPs, counterfactuals, scoring,
#  prediction bands - without any download, and checks that:
#    - targets before May 2014 and in 2020-2022 are skipped;
#    - no training window touches 2020-2022;
#    - the rolling curves track the falling risk and beat the published ones;
#    - the 95% prediction band covers roughly 95% of days.
#
#    Rscript tests/test_rolling_synthetic.R      (from the repository root)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dlnm); library(mixmeta); library(splines)
})
for (f in c("R/config.R", "R/utils.R", "R/erf.R", "R/curves.R", "R/episodes.R",
            "R/recalibrate.R", "R/rolling.R", "R/rolling_eval.R")) source(f)
set.seed(1)
cfg$nsim <- 200
cfg$rolling$pool <- "region"
cfg$rolling$min_deaths <- 0
check <- function(ok, msg) {
  cat(if (isTRUE(ok)) "  ok   " else "  FAIL ", msg, "\n", sep = "")
  if (!isTRUE(ok)) failed <<- c(failed, msg)
}
failed <- character()

# ---- fake world ---------------------------------------------------------------
codes <- sprintf("IT9%02dC", 1:6)
dates <- seq(as.Date("2009-01-01"), as.Date("2026-06-30"), by = "day")
doy <- as.numeric(format(dates, "%j"))
temp_of <- function(k) {
  t <- 14 + 10 * sin(2 * pi * (doy - 110) / 365.25) + as.numeric(arima.sim(list(ar = .8), length(dates), sd = 1.3))
  # a few heatwaves every summer
  for (y in unique(year(dates))) for (h in 1:3) {
    s <- as.Date(sprintf("%d-06-15", y)) + sample(0:70, 1)
    i <- which(dates >= s & dates < s + sample(4:12, 1))
    t[i] <- t[i] + runif(1, 3, 6)
  }
  t + k * .3
}
tmp <- rbindlist(lapply(seq_along(codes), function(k)
  data.table(URAU_CODE = codes[k], date = dates, tmean = temp_of(k))))
pcts <- c(0, .1, .5, 1, 2.5, seq(5, 95, by = 5), 97.5, 99, 99.5, 99.9, 100)
pcts <- sort(unique(c(pcts, 1:99)))
tdist <- tmp[year(date) <= 2019, as.list(setNames(quantile(tmean, pcts / 100),
                                                   sprintf("%.1f%%", pcts))), by = URAU_CODE]
tdist[, tmax := `100.0%`]

# true heat slope per degree above the MMT: falls from .06 (2011) to .03 (2026)
heat_slope <- function(d) .06 - .03 * (as.numeric(d - as.Date("2011-01-01")) / (16 * 365.25))
true_logrr <- function(t, mmt, d) ifelse(t > mmt, heat_slope(d) * (t - mmt), .012 * pmax(mmt - t, 0))
# "published" curve: the 2011 slope, steeper than anything after it
fit_beta <- function(code, f) {
  sp <- list(knots = pct_value(tdist[URAU_CODE == code], c(10, 75, 90)),
             bound = pct_value(tdist[URAU_CODE == code], c(0, 100)), degree = 2)
  g <- seq(sp$bound[1], sp$bound[2], length.out = 200)
  X <- unclass(bs(g, knots = sp$knots, degree = 2, Boundary.knots = sp$bound))
  as.numeric(coef(lm(f(g) ~ 0 + X)))
}
ages <- AGE_GROUPS; share <- c(.03, .12, .17, .3, .38)
coefs <- rbindlist(lapply(codes, function(cd) {
  mmt <- pct_value(tdist[URAU_CODE == cd], 80)
  b <- fit_beta(cd, function(t) ifelse(t > mmt, .07 * (t - mmt), .012 * pmax(mmt - t, 0)))
  rbindlist(lapply(ages, function(a) data.table(URAU_CODE = cd, agegroup = a,
                                                b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5])))
}))
vc <- coefs[, .(URAU_CODE, agegroup)]
for (i in 1:5) for (j in 1:5) vc[, (sprintf("v%d%d", i, j)) := if (i == j) 4e-4 else 1e-4]
erf <- list(tdist = tdist, coefs = coefs, vcov = vc, simu = NULL,
            metadata = data.table(URAU_CODE = codes, name = paste("City", 1:6), inmcc = FALSE,
                                  pop = 1, lat = 45, lon = 9))

mort <- rbindlist(lapply(seq_along(codes), function(k) {
  cd <- codes[k]; x <- tmp[URAU_CODE == cd]
  mmt <- pct_value(tdist[URAU_CODE == cd], 80)
  # 2-day mean temperature drives the risk (a short lag)
  t2 <- frollmean(x$tmean, 2, fill = x$tmean[1])
  # two small cities (4 and 8 deaths a day), as Lecco or Cremona, and larger ones
  base <- c(4, 8, 40, 70, 100, 150)[k] * (1 + .15 * cos(2 * pi * (doy - 15) / 365.25))
  covid <- ifelse(year(dates) %in% 2020:2022, 1 + .4 * exp(-((doy - 90) / 25)^2), 1)
  mu <- base * covid * exp(true_logrr(t2, mmt, dates))
  rbindlist(lapply(seq_along(ages), function(a)
    data.table(URAU_CODE = cd, agegroup = ages[a], date = dates, deaths = rpois(length(mu), mu * share[a]))))
}))
mort <- mort[date >= as.Date("2011-01-01")]
setkey(mort, URAU_CODE, agegroup, date)
data_end <- max(mort$date)
hist_expo <- tmp[date <= as.Date("2019-12-31")]
test_expo <- tmp[date >= as.Date("2022-12-01")]

# ---- windows ----------------------------------------------------------------
cat("Windows\n")
plan <- rolling_windows(min(mort$date), data_end, cfg)
tg <- plan$targets
check(all(tg[year == 2013, status != "ok"]) && all(tg[year < 2014, status != "ok"]),
      "targets before 2014 skipped (no three years of history)")
check(tg[key == "2014-05", status] == "ok", "May 2014 is the first target")
check(all(tg[year %in% 2020:2022, status != "ok"]), "COVID years not predicted")
touch <- vapply(plan$windows, function(w) any(year(seq(min(w$blocks$start), w$train_end, by = "day")[
  in_blocks(seq(min(w$blocks$start), w$train_end, by = "day"), w$blocks)]) %in% 2020:2022), TRUE)
check(!any(touch), "no training window uses 2020-2022")
rec_w <- Filter(function(w) w$kind == "recent", plan$windows)
check(all(vapply(rec_w, `[[`, 0, "n_days") == round(365.25 * 3)), "every recent window has 3 years of days")
w24 <- rec_w[[which(vapply(rec_w, function(w) w$train_end == as.Date("2024-06-30"), TRUE))]]
cat("  July 2024 trained on:", w24$label, "\n")
check(nrow(w24$blocks) == 2, "a post-COVID window is two blocks")

# ---- fits ---------------------------------------------------------------------
cat("Fits\n")
expo <- rolling_exposure(hist_expo, test_expo)
t0 <- Sys.time()
fits <- lapply(plan$windows, fit_rolling_window, mort = mort, expo = expo, erf = erf, cfg = cfg,
               pool = codes, keep = codes)
cat(sprintf("  %d windows in %.0f s\n", length(fits), as.numeric(Sys.time() - t0, units = "secs")))
sets <- rolling_collect(plan, fits, cfg, erf)
roll <- sets$window; shift <- sets$shift
n_recent <- sum(vapply(plan$windows, `[[`, "", "kind") == "recent")
check(nrow(roll$windows) == n_recent, "every recent window produced curves")
check(!is.null(shift) && any(!is.na(shift$shifts$anchor_id)), "shift set built on anchor windows")
an <- Filter(function(w) w$kind == "anchor", plan$windows)
cat("  anchors:", paste(vapply(an, `[[`, "", "label"), collapse = " | "), "\n")
check(all(vapply(an, `[[`, 0, "n_days") > round(365.25 * 3)), "anchors longer than the recent window")
s <- roll$summary
tr <- s[, .(rr = mean(rr99), pub = mean(rr99_pub)), by = .(y = year(train_end))]
print(tr[order(y)])
check(tr[y == min(y), rr] > tr[y == max(y), rr], "rolling RR99 falls over time (adaptation)")
check(all(tr$rr < tr$pub), "rolling RR99 below the (too steep) published one")
# stability and accuracy of the two sets, against the true RR at P99
true_rr99 <- function(cd, d) {
  r <- tdist[URAU_CODE == cd]
  exp(heat_slope(d) * (pct_value(r, 99) - pct_value(r, 80)))
}
cmp <- rbindlist(lapply(list(window = roll, shift = shift), function(x)
  x$summary[, .(URAU_CODE, train_end, rr99)]), idcol = "set")
cmp[, truth := mapply(true_rr99, URAU_CODE, train_end + 1)]
cmp[, small := URAU_CODE %in% codes[1:2]]
setorder(cmp, set, URAU_CODE, train_end)
st <- cmp[, .(month_to_month = mean(abs(diff(log(rr99)))),
              rmse_vs_truth = sqrt(mean((log(rr99) - log(truth))^2))), by = .(set, small)]
print(st)
check(st[set == "shift" & small == TRUE, month_to_month] < st[set == "window" & small == TRUE, month_to_month],
      "small cities: shift curves move less from month to month than window curves")
# accuracy over all cities. (In this fake world every city has the same true
# risk, so the heaviest shrinkage - the window BLUP of a small city, which is
# almost the pooled curve - is the best possible answer for small cities; the
# shift set gives them their own 10-year shape instead, which only pays off
# where cities really differ. Accuracy is therefore checked over all cities.)
acc <- cmp[, .(rmse = sqrt(mean((log(rr99) - log(truth))^2))), by = set]
check(acc[set == "shift", rmse] < acc[set == "window", rmse],
      "all cities: shift curves closer to the true risk than window curves")
check(identical(rolling_lookup(roll, as.Date("2024-07-15")), "w20240630"),
      "lookup: 15 July 2024 uses the window ending 30 June 2024")
check(is.na(rolling_lookup(roll, as.Date("2021-07-15"))), "lookup: no curve inside the COVID gap")

# ---- scoring --------------------------------------------------------------------
cat("Scoring\n")
cfg2 <- cfg; cfg2$nsim <- 200
ep_h <- detect_episodes(hist_expo, erf, as.Date("2011-01-01"), as.Date("2019-12-31"), cfg2)
cf_h <- fit_counterfactual(mort, hist_expo, ep_h, erf, as.Date("2011-01-01"), as.Date("2019-12-31"), cfg2)
ev_h <- evaluate_counterfactual(cf_h, cfg2)
ep_t <- detect_episodes(test_expo, erf, as.Date("2023-01-01"), data_end, cfg2)
cf_t <- fit_counterfactual(mort, test_expo, ep_t, erf, as.Date("2023-01-01"), data_end, cfg2)
ev_t <- evaluate_counterfactual(cf_t, cfg2)
check("phi" %in% names(cf_t$daily), "counterfactual daily carries the dispersion")
recal <- list(recalibrate_curves(mort, hist_expo, erf, c("2016-01-01", "2019-12-31"), cfg2))
re <- rolling_evaluation(
  list(list(label = "2011-2019", cf = cf_h, ev = ev_h, recals = FALSE),
       list(label = "2023-2026", cf = cf_t, ev = ev_t, recals = TRUE)),
  sets, recal, erf, cfg2, codes)
cal <- re$calibration
print(cal[type == "episode", .(period, subset, curves, n, ratio = round(ratio, 2),
                               slope = round(slope, 2), cover = round(coverage_95, 2))])
a <- cal[period == "All periods" & subset == "All windows" & type == "episode"]
check(abs(a[curves == roll$label, ratio] - 1) < abs(a[curves == PUBLISHED_LABEL, ratio] - 1),
      "rolling curves closer to calibration than the published ones")
lg <- cal[period == "All periods" & subset == "Large episodes" & type == "episode"]
check(nrow(lg) >= 2 && abs(lg[curves == roll$label, slope] - 1) < abs(lg[curves == PUBLISHED_LABEL, slope] - 1),
      "large episodes: rolling slope closer to 1")
check(setequal(unique(re$windows[period == "2011-2019", curves]),
               c(PUBLISHED_LABEL, roll$label, shift$label)),
      "fixed recalibrated curves not scored where they are in sample")

# ---- bands --------------------------------------------------------------------
cat("Bands\n")
b <- rolling_daily_band(cf_t, roll, erf, cfg2, codes[3], as.Date("2023-06-01"), as.Date("2025-09-30"))
cov <- mean(b$flag == "within")
cat(sprintf("  daily band coverage %.3f (above %.3f, below %.3f)\n", cov,
            mean(b$flag == "above"), mean(b$flag == "below")))
check(cov > .90 && cov < .99, "daily 95% band covers about 95% of days")
b7 <- rolling_daily_band(cf_t, roll, erf, cfg2, codes[3], as.Date("2023-06-01"), as.Date("2025-09-30"), k = 7)
check(mean(b7$pi_hi - b7$pi_lo) < mean(b$pi_hi - b$pi_lo), "7-day band narrower than the daily one")
ex <- rolling_example(list(list(label = "pre", cf = cf_h), list(label = "post", cf = cf_t)),
                      roll, erf, data.table(URAU_CODE = codes, name = paste("City", 1:6)), cfg2, codes)
check(!is.null(ex) && all(c(2014:2019, 2023:2025) %in% ex$data$year), "example covers every scored summer")

out <- Sys.getenv("ROLLING_TEST_OUT")   # keep the objects, e.g. to try the app on them
if (nzchar(out)) saveRDS(list(roll = roll, sets = sets, re = re, ex = ex, cf_t = cf_t, erf = erf,
                              codes = codes, cfg = cfg2, tmp = tmp, recal = recal), out)
cat(if (length(failed)) sprintf("\n%d check(s) failed\n", length(failed)) else "\nall checks passed\n")
if (length(failed)) quit(status = 1)
