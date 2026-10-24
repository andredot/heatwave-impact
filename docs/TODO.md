# Open items

## 1. Replace the geometric city crosswalk with an official correspondence table
Status: DONE (branch istat-city-composition). The source is the Istat SITUAS report
*Composizione delle Città (Cities) 2021* (EU Regulation 2019/1130): 89 Cities, 11 of
them Greater Cities, 298 comuni. It carries both keys we need - the Urban Audit code
of Masselot's curves and the six-digit Istat comune code of the mortality file - so
the crosswalk needs no name matching and no population heuristic. A copy lives in
`reference/istat_cities_2021.csv` (28 kB, fixed 2021 definition), so a fresh clone
reproduces the geography with no download; `R/istat_cities.R` documents how to
re-export it from https://situas.istat.it/web/ for a future release.

`cfg$city_source` selects `"istat"` (default), `"official"` (Eurostat CITY-LAU table
via `city_table_url`) or `"geometry"`. Whichever is chosen, the geometric matching is
still computed and `finalize_crosswalk()` reports, city by city, how many comuni the
two compositions share (`n_comuni`, `n_comuni_geo`, `shared_comuni`, `agreement`), so
the inference this project started from is now measured rather than trusted. All the
earlier checks (codes present in Istat, name, population ratio, deaths ratio against
the metadata rates) are unchanged and the override file still applies.

Not used: the FUA report from the same portal (pfun 446). A Functional Urban Area is
the city plus its commuting belt and is not interchangeable with a City - it would
put a far larger population under the same curve. Masselot's unit is the City.

Still to watch: comune code changes between the 2021 composition and the Istat
mortality release (handled by `in_istat` + `missing_pop_tol`), and the fact that
exposure still averages ERA5-Land over the GISCO Urban Audit polygon chosen
geometrically, not over the union of the official comuni.

## 2. Download the Istat mortality file instead of receiving it by hand
Status: DONE (branch inputs-and-scenarios): `download_istat()` + `istat_url` in config.

Wanted: a `download_istat()` function plus an `istat_file` target that fetches and
unzips the release, so a fresh clone runs end to end. Example of the URL actually
used: https://www.istat.it/storage/dati_mortalita/giugno-2026/decessi-comunali-mese-provvisori-3.zip

Notes for the implementation:
- the path carries the release month ("giugno-2026") and a version suffix ("-3"),
  so put the full URL in config.R (`istat_url`) rather than guessing it, and let a
  manually placed file win if it is already there;
- unzip into data/raw/, find the CSV inside rather than assuming its name (the name
  has changed between releases), and keep `format = "file"` on the target so a new
  release re-runs the pipeline;
- check the last date with data after download (`detect_data_end`) and report it,
  since the release month and the data cutoff are not the same thing.

## 3. Input audit: what a fresh run downloads and what it does not
Status: DONE - the Istat release and the FluMOMO code are now downloaded too;
only the optional FluNet export and the optional crosswalk overrides stay manual.

Downloaded automatically (cached under data/):
- Zenodo files: coefs, vcov, tmean_distribution, metadata, additional_data.zip,
  results.zip (`download_zenodo`)
- GISCO boundaries: Urban Audit city polygons for the configured releases, LAU
  polygons, NUTS-2 for the app (`download_gisco`, `export_app_bundle`)
- ERA5-Land daily temperatures, Open-Meteo archive (`fetch_era5_all`)
- Archived forecasts, Open-Meteo previous runs (`fetch_forecasts_all`)
- Influenza activity from ECDC ERVISS (`build_influenza_activity`)

NOT downloaded, must be put in place by hand (a fresh run fails or degrades):
- Istat daily deaths by comune -> item 2 above (fails)
- Official FluMOMO 4.2 code in the folder set by `flumomo$code_dir` (the FluMOMO
  targets fail; the rest of the pipeline is unaffected). Check whether the EuroMOMO
  terms allow automating the download; if not, keep it manual and documented.
- WHO FluNet export at `influenza$flunet_file`, only needed for influenza activity
  before 2022: the public WHO API returns 403 (degrades: ERVISS alone covers
  2022 onwards)
- `data/raw/crosswalk_overrides.csv` (optional, only to force a city definition)

## 4. Scenario output for policymakers ("impact-based forecast")
Status: FIRST VERSION DONE (branch inputs-and-scenarios): "Scenari" tab in the app.
Still to do: onset probability from a real ensemble rather than from the archive error
spread, and forecast error as a third variance component in the intervals.

Today the app answers "how many heat deaths over the days I selected". A
decision maker needs the question the other way round: how likely is an episode,
over which territory, and what does each extra day of it cost.

Target format, one block per territory:

  "Probability of a heatwave starting in the next N days: P%, over <territory>.
   Current forecast: D deaths (credible interval A-B) over the days forecast.
   If it lasts L days: total D1 (A1-B1).
   If it continues 5 days beyond that: additional D2, total D3 (A3-B3).
   If it continues 5 days further: additional D4, total D5 (A5-B5)."

Design notes:
- probability of onset: from the forecast ensemble if available, otherwise from
  the observed detection skill by lead (POD/FAR are already in `fc_skill`), which
  turns "forecast says hot" into a calibrated probability;
- scenario lengths rather than a single horizon: hold the peak temperature at the
  forecast level and extend the episode, so each extra block of days is a
  conditional increment, not a new forecast;
- intervals must include the forecast error by lead (`fc_impact` has it) on top of
  curve and baseline uncertainty, otherwise the ranges are too narrow;
- state deaths brought forward vs added: the displacement finding means a long
  episode's total is not simply additive, and policymakers should see that;
- use the recalibrated Italian curves for the numbers, the published ones only as
  a comparison;
- the same block should work for a region (all its cities) and for one city.

Later: medium-range forecasts (sub-seasonal) would extend the onset probability
beyond 7 days; not needed for the first version.

## 5. Ratio plot: does over-prediction get worse with episode size?
Status: DONE (branch inputs-and-scenarios): ratio panel in the report, after the
calibration scatter. The asinh-scale variant of the scatter is still open.

The calibration scatter shows observed vs predicted with a fitted slope. It does
not answer directly whether the over-prediction is proportional or worsens with
size, which matters because the recalibrated curves reach a ratio of totals near
1 while their slope stays around 0.79 - meaning the largest episodes stay
over-predicted relative to the small ones, and those are the episodes a warning
system exists for.

Wanted: observed / predicted on the vertical axis against predicted on the
horizontal, one panel per window type, a reference line at 1, one colour per curve
set (published and recalibrated). Points sitting on a downward trend mean the
over-prediction grows with episode size.

Notes:
- keep only windows with a predicted value above a small threshold: the ratio is
  unstable when the denominator is near zero, which is most of the cloud near the
  origin in the current scatter;
- weight or size the points by precision so the eye is not drawn to noisy windows;
- keep the existing scatter as well: it shows the negative observed values, which a
  ratio plot cannot, and those are what demonstrate the noise is symmetric;
- consider a signed square-root (or asinh) scale on the current scatter as a
  separate improvement: it decompresses the crowd near zero while keeping negatives.

## 6. Curves for every city of a region (BLUPs), and then for the app
Status: PARTLY DONE (branch istat-city-composition). The report has a section
"One region in full" driven by `cfg$region_focus` (default Lombardia): every city
of the region is refitted with Masselot's first stage with the death threshold
lifted (`run_curve_tests(..., min_fit_deaths = 0)`), and the shrinkage is left to
the second stage - `blup_curves(..., extra_fits = , keep = )` pools the city with
every Italian fit, so a noisy city gets a large variance and its BLUP sits close
to the pooled curve. The table reports the published, own-fit and BLUP relative
risks with intervals, and the figure shows the three curves per city; Lombardy
resolves to nine cities (Bergamo, Brescia, Como, Cremona, Lecco, Milano, Pavia,
Varese, Vigevano).

Still to do:
- carry the BLUP curves into the app bundle, so the small cities are usable
  operationally and not only visible in the report. Today `app_export.R` exports
  the published and recalibrated curve sets; the regional BLUP set would be a
  third option in the curve selector.
- validate the transfer rather than assume it: leave-one-city-out, refit the
  meta-analysis without a city, predict its curve from the others and compare with
  its own fit. That is the honest test of whether a curve can be lent to a city
  that has too few deaths of its own.
- decide the period the operational BLUP curves are fitted on. The report uses the
  test period, which is the right choice for validating; 2016-2019 is what restored
  calibration out of sample, so it is probably the right choice for forecasting.
- Massa (`IT045C`) has a published curve and no composition: it is a City in the
  earlier Urban Audit releases and not one in 2021. It could be added by hand
  (one comune) through `data/raw/crosswalk_overrides.csv` if the city matters, but
  the override path currently needs the city to have a candidate polygon too.

## 7. Rolling recalibration and prediction bands (Lombardy first)
Status: FIRST VERSION (branch rolling-recalibration-lombardy). Not yet run on the
real data.

What it does:
- `R/rolling.R`: for every target month (May-September by default) Masselot's first
  stage is refitted on the last `rolling$train_years` (3) years of Istat deaths
  available before that month, skipping 2020-2022, and pooled/shrunk with the BLUPs.
  A target without three years of Istat data before it (anything before May 2014
  with the current file) is skipped. The time trend is one natural spline per
  contiguous block of training days, so a window spanning 2019 -> 2023 does not
  spend degrees of freedom on the gap. One `targets` branch per training window.
- `R/rolling_eval.R`: published, fixed recalibrated (post-COVID only, where they are
  out of sample) and rolling curves scored on the same windows with the same
  uncertainty model (curve draws + baseline level + negative binomial noise).
  Calibration overall and for "large" episodes (top third by the published
  prediction, which does not select on observed deaths). Pre-COVID years are
  scored with a second counterfactual on 2011-2019 (Masselot's temperatures).
- Prediction bands: `rolling_daily_band()` for the report, and a grey band in the
  app (around the prediction or around the expectation without heat), with the
  days whose observed 7-day mean leaves it marked as abnormal.
- App: rolling curves are the first curve set when the bundle has them; each day
  uses the curve of the latest window ending before it. Cities without rolling
  curves (outside the scope) fall back to the last fixed recalibrated set.

To check on the first real run:
- run time of `roll_fits` with `pool = "italy"` (about 50 windows x 80 cities);
  `pool = "region"` is the quick option but pools only nine cities;
- whether the large-episode slope actually moves towards 1: if it does not, the
  remaining gap is probably displacement within long episodes (Milan 2024: 29 days)
  rather than the level of the curve, and the next step is a depletion term;
- whether the grey band covers ~95% of summer days in the pre-COVID years; if it is
  too narrow, the baseline level uncertainty (B_logsd) is the first suspect.
- the yardstick itself. In the synthetic check (`tests/test_rolling_synthetic.R`,
  where the truth is known) the counterfactual baseline came out 3.5% above the true
  no-heat level in every month, because the published cold arm used as offset did
  not match the true one; that alone cut the "observed" heat deaths of episodes by
  ~20% (8,800 vs 10,900 true), so even the right curve looked like over-prediction.
  On the real data, re-run the evaluation with `baseline_cold_offset = FALSE` and see
  how much of the remaining gap on large episodes moves.

### 7b. Stability of the rolling curves in small cities
Three years are enough for Milan, not for a city with a few deaths a day: its own
fit is noise and its window BLUP is mostly the pooled curve, moving month to month.
A longer window with a time trend does not help: the curve read at the end of a
linear trend over T years has about the variance of a T/4-year window, so ten
years with a trend is no steadier than three without. What helps is borrowing
the trend itself, which the validation found to be common to all cities.

Now `rolling_collect()` builds two sets from the same fits, both scored in the report:
- `window`: city BLUP on the last 3 years (as before);
- `shift` (default, `rolling$method`): city BLUP on up to `anchor_years` (10) years
  ending on the last 31 December (one anchor per year, so ~10 extra fits per city)
  + the change common to all pooled Italian cities between that anchor and the last
  3 years (meta-analysis of each city's recent-minus-anchor coefficients).
Synthetic check: month-to-month movement of log RR99 down by about a third, and
closer to the true risk over all cities. In a world where every city has the same
risk the heavy shrinkage of the window set is already ideal for small cities, so
the real-data comparison (report: "Stability of the two rolling sets", calibration
table) is what should decide `rolling$method`.

## 8. Two reports and an all-summers app (branch rolling-recalibration-lombardy)
Status: FIRST VERSION, checked on simulated data only.

- `reports/adaptation_report.qmd` (paper): every Italian city refitted in consecutive
  3-year windows (2011-13, 2014-16, 2017-19, 2023-26; COVID years out), BLUPs;
  RR99 vs published with meta-regression trend (all cities and balanced panel),
  MMT percentile, heat-attributable fraction standardised to the same summers and
  deaths, and a curve-free excess-per-degree-day check. Code: `R/adaptation.R`,
  targets `adaptation`, `dose_response`, `timing_italy`.
- `reports/system_validation.qmd`: the rolling 3-year curves (published as reference)
  on Lombardy, 2014-19 + 2023-26, recorded temperatures; MMT, uncertainty,
  calibration, warnings confusion matrix, readiness verdict against criteria fixed
  in `cfg$readiness`. Targets `roll_eval`, `roll_example`, `timing_region`, `readiness`.
- Dropped from the pipeline: fixed 2011-2019 / 2016-2019 curves, the 10-year
  "shift" variant (still available: `rolling$anchor_years > train_years`), forecast
  skill and fc_impact (forecast validation comes back in a later version), the old
  single report. The functions stay in R/ (recalibrate.R, forecast.R,
  run_period_tests, forecast_episode_impact, flumomo_validation.R) for reuse.
- Rolling curves now cover every city (`rolling$scope = "all"`) and April
  (`target_months = 4:9`), and the 2011-2019 counterfactual is fitted for every city:
  the app shows every summer (1 April - 30 September) from the bundle, without
  downloads. Expect a longer first `tar_make()`.

## 9. Validation fixes after the first real run (branch rolling-recalibration-lombardy)
- First real run: 245 city-episodes, only 8 with a significant excess (6 in Milan);
  calibration chart shown on those 8 looked good (ratio 1.08) because they were
  selected on the outcome; all episodes 0.48. Whole summers 2023-2026: observed
  excess ~0 while the curves predict ~2,500 -> yardstick suspected.
- Done: touching/overlapping episode windows merged into whole heatwaves
  (`hw_merge_windows`); scoring stays city by city, with a per-city table
  (`roll_eval$cities`); baseline audit (`R/baseline_audit.R`, target
  `baseline_audit`): cool-day observed/expected and the whole evaluation repeated
  with `cfg$baseline_variants`. Report 1 windows now 5 years where possible
  (2011-2015, 2016-2019, 2023-2026).
- Synthetic check of the audit: cool days pass (99%) while two harmonics move the
  ratio 0.71 -> 0.88: the sensitivity table, not the cool-day check, is the test.
- Still open: probabilistic warnings (Brier, reliability, ROC) and confidence
  intervals on the scores; readiness criteria amendments must be dated.
