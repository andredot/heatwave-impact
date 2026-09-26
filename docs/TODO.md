# Open items

## 1. Replace the geometric city crosswalk with an official correspondence table
Status: DONE (branch inputs-and-scenarios). Set `city_table_url` in R/config.R to the
Eurostat "CITY - LAU 2021" file; with an empty URL the pipeline falls back to the
geometric matching and says so.

Today `overlay_cities()` + `finalize_crosswalk()` infer which comuni make up each
Urban Audit city: the city point is located in the GISCO polygons of two releases
and two levels, and the candidate whose population is closest to Masselot's is
chosen. It works (85 cities, population ratios 0.9-1.1, Milan as a 96-comune
greater city), but it is inference where an authoritative list exists.

Wanted instead: an official file listing, for each city / greater city / metro
area, the member municipalities with their Istat codes. Candidate sources to
check:
- Eurostat Cities (Urban Audit) database: https://ec.europa.eu/eurostat/web/cities/database
  (the LAU-city correspondence table Masselot himself used: "CITY-LAU-2021")
- Eurostat metropolitan regions: https://ec.europa.eu/eurostat/web/metropolitan-regions/database
- Eurostat LAU-NUTS correspondence (annual xlsx, has CITY_ID / GREATER_CITY_ID / FUA_ID columns)
- Istat Urban Audit quality note: https://www.istat.it/scheda-qualita/urban-audit/
- Istat codes of comuni, province and regions:
  https://www.istat.it/classificazione/codici-dei-comuni-delle-province-e-delle-regioni/

Keep: every validation check in `finalize_crosswalk()` (codes present in Istat,
name match, population ratio, deaths ratio against the metadata rates) and the
report section that shows them. They become the check on the official table
rather than on our inference, and the override file stays.

Watch for: the code renumbering between Urban Audit releases (Masselot's codes are
the 2020/2021 numbering, where greater cities are merged into the city level), and
comune code changes between the correspondence table's year and the Istat file.

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
