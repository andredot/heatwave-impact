# =============================================================================
#  Wrapper around the official FluMOMO code (version 4.2, EuroMOMO).
#
#  The published scripts are used unchanged: this file only prepares their
#  input files, sets the parameters that `FluMOMO_v42.R` would set by hand,
#  and runs `Estimation_v42.R`. Nothing here re-implements the model.
#
#  Unzip the official archive into vendor/flumomo/, so that it holds
#  Estimation_v42.R and the Output_*.R scripts. Keep it outside R/, which
#  targets sources automatically.
# =============================================================================

#' Monday of the ISO week of a date
#'
#' @param date Dates.
#' @return Dates of the Mondays.
iso_week_start <- function(date) date - (as.integer(format(date, "%u")) - 1L)

#' EuroMOMO age group of an Istat age class
#'
#' FluMOMO expects `agegrp` 0 to 4 for 0-4, 5-14, 15-64, 65+ and Total.
#' Istat `CL_ETA`: 0 = under 1, 1 = 1-4, 2 = 5-9, ..., 21 = 100+.
#'
#' @param cl Integer vector of `CL_ETA` codes.
#' @return Integer vector 0-3 (`NA` if unknown); group 4 (Total) is added by
#'   [istat_weekly_euromomo()].
euromomo_agegrp <- function(cl) {
  cl <- as.integer(cl)
  data.table::fcase(cl %in% 0:1, 0L, cl %in% 2:3, 1L, cl %in% 4:13, 2L,
                    cl >= 14, 3L, default = NA_integer_)
}

#' Weekly deaths by EuroMOMO age group for a set of provinces
#'
#' Reads the Istat daily file, keeps every comune of the given provinces and
#' aggregates to ISO weeks and to the four EuroMOMO age groups plus the total.
#' Incomplete weeks at the edges of the data are dropped.
#'
#' @param path Istat CSV.
#' @param prov Province codes (three digits), or `NULL` when `procom` is given.
#' @param years Years to keep.
#' @param procom Comune codes (six digits) to keep instead of whole provinces.
#' @return data.table `agegrp`, `year`, `week`, `deaths`.
istat_weekly_euromomo <- function(path, prov, years, procom = NULL) {
  hdr <- names(data.table::fread(path, nrows = 0))
  tcols <- grep("^T_[0-9]{2}$", hdr, value = TRUE)
  tcols <- tcols[(as.integer(sub("T_", "", tcols)) + 2000L) %in% years]
  d <- data.table::fread(path, select = c("COD_PROVCOM", "CL_ETA", "GE", tcols),
                         colClasses = list(character = c("COD_PROVCOM", "GE")),
                         na.strings = c("", "NA", "n.d."))
  d[, PRO_COM := pad_procom(COD_PROVCOM)]
  d <- if (!is.null(procom)) d[PRO_COM %in% procom] else d[substr(PRO_COM, 1, 3) %in% prov]
  d[, agegrp := euromomo_agegrp(CL_ETA)]
  long <- data.table::melt(d[!is.na(agegrp)], id.vars = c("agegrp", "GE"),
                           measure.vars = tcols, variable.name = "yy",
                           value.name = "deaths")
  long[, GE := sprintf("%04d", as.integer(GE))]
  long[, date := as.Date(sprintf("%d-%s-%s", 2000L + as.integer(sub("T_", "", yy)),
                                 substr(GE, 1, 2), substr(GE, 3, 4)))]
  long <- long[!is.na(date) & !is.na(deaths)]
  w <- long[, .(deaths = sum(deaths), days = data.table::uniqueN(date)),
            by = .(agegrp, week_start = iso_week_start(date))]
  tot <- w[, .(agegrp = 4L, deaths = sum(deaths), days = max(days)), by = week_start]
  w <- rbind(w, tot)[days == 7]
  # the official code loops over age groups 0 to 4: all must be present
  w <- merge(data.table::CJ(agegrp = 0:4, week_start = unique(w$week_start)),
             w, by = c("agegrp", "week_start"), all.x = TRUE)
  w[is.na(deaths), deaths := 0L]
  w[, `:=`(year = as.integer(format(week_start, "%G")),
           week = as.integer(format(week_start, "%V")))]
  w[order(agegrp, year, week), .(agegrp, year, week, deaths)]
}

#' Daily deaths of a set of provinces (all ages or from an age)
#'
#' Used for the observed series in the figures, where a daily resolution is
#' wanted rather than the weekly one FluMOMO works on.
#'
#' @param path Istat CSV.
#' @param prov Province codes.
#' @param age_min Minimum age (0 = all ages).
#' @param years Years to keep.
#' @return data.table `date`, `deaths`.
istat_daily_area <- function(path, prov, age_min = 0, years = 2015:2026) {
  hdr <- names(data.table::fread(path, nrows = 0))
  tcols <- grep("^T_[0-9]{2}$", hdr, value = TRUE)
  tcols <- tcols[(as.integer(sub("T_", "", tcols)) + 2000L) %in% years]
  d <- data.table::fread(path, select = c("COD_PROVCOM", "CL_ETA", "GE", tcols),
                         colClasses = list(character = c("COD_PROVCOM", "GE")),
                         na.strings = c("", "NA", "n.d."))
  d[, PRO_COM := pad_procom(COD_PROVCOM)]
  d <- d[substr(PRO_COM, 1, 3) %in% prov]
  if (age_min > 0) d <- d[CL_ETA >= (if (age_min == 20) 5L else floor(age_min / 5) + 1L)]
  long <- data.table::melt(d, id.vars = "GE", measure.vars = tcols,
                           variable.name = "yy", value.name = "deaths")
  long[, GE := sprintf("%04d", as.integer(GE))]
  long[, date := as.Date(sprintf("%d-%s-%s", 2000L + as.integer(sub("T_", "", yy)),
                                 substr(GE, 1, 2), substr(GE, 3, 4)))]
  long[!is.na(date) & !is.na(deaths), .(deaths = sum(deaths)), by = date][order(date)]
}

#' Write the three FluMOMO input files
#'
#' Formats expected by `Estimation_v42.R`: `deaths.txt` (`agegrp`, `year`,
#' `week`, `deaths`), `wdata_<code>.txt` (`date`, `pop3`, `NUTS3`, `temp`:
#' daily temperature by NUTS-3 area with its population, which the official
#' code turns into a population-weighted weekly mean and into the extreme
#' temperature variable) and `IA.txt` (`agegrp`, `year`, `week`, `IA`), all
#' semicolon-separated in `<work_dir>/data`.
#'
#' @param work_dir Working directory for the run.
#' @param deaths Output of [istat_weekly_euromomo()].
#' @param weather data.table `date`, `pop3`, `NUTS3`, `temp` (daily).
#' @param ia data.table `year`, `week`, `IA`, or `NULL` (then IA = 0 and the
#'   influenza term carries no information: the baseline will absorb winter
#'   influenza, which is not FluMOMO as intended).
#' @param country_code Two-letter code used in the weather file name.
#' @return The data directory.
flumomo_write_inputs <- function(work_dir, deaths, weather, ia, country_code) {
  dir <- ensure_dir(file.path(work_dir, "data"))
  utils::write.table(deaths, file.path(dir, "deaths.txt"), sep = ";", dec = ".",
                     row.names = FALSE, quote = FALSE)
  w <- data.table::copy(weather)
  w[, date := format(as.Date(date), "%Y-%m-%d")]
  utils::write.table(w[, .(date, pop3, NUTS3, temp)],
                     file.path(dir, sprintf("wdata_%s.txt", country_code)),
                     sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
  grid <- unique(deaths[, .(agegrp, year, week)])
  ia_tab <- if (is.null(ia)) {
    warning("No influenza activity supplied: IA set to 0, so the baseline will ",
            "contain influenza-related mortality.", call. = FALSE)
    grid[, IA := 0][]
  } else {
    merge(grid, unique(ia[, .(year, week, IA)]), by = c("year", "week"), all.x = TRUE)[
      is.na(IA), IA := 0][]
  }
  utils::write.table(ia_tab[, .(agegrp, year, week, IA)], file.path(dir, "IA.txt"),
                     sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
  dir
}

#' Run the official FluMOMO estimation
#'
#' Sets the parameters that `FluMOMO_v42.R` sets interactively and sources
#' `Estimation_v42.R` unchanged, in its own environment. Deaths, weather and
#' influenza activity are taken from the files written by
#' [flumomo_write_inputs()] (`A_MOMO = 0`, `WeatherData = 0`).
#'
#' @param code_dir Folder with the official scripts.
#' @param work_dir Working directory holding `data/`.
#' @param country,country_code Labels used in file names.
#' @param start,end Dates bounding the period (converted to ISO year/week).
#' @param ia_lags,et_lags Number of weekly lags for IA and ET (official default 2).
#' @param ia_restricted `IArest` in the official code.
#' @param verbose Print the estimation output.
#' @return data.table of the official results (one row per age group and week),
#'   with `EB` the baseline in counts, `EdIA` and `EdET` the deaths attributed
#'   to influenza and to extreme temperature, and their confidence limits.
flumomo_run <- function(code_dir, work_dir, country, country_code, start, end,
                        ia_lags = 2, et_lags = 2, ia_restricted = TRUE,
                        verbose = FALSE) {
  est <- file.path(code_dir, "Estimation_v42.R")
  if (!file.exists(est))
    stop("Estimation_v42.R not found in ", code_dir,
         ": unzip the official FluMOMO 4.2 archive there.")
  if (!requireNamespace("ISOweek", quietly = TRUE))
    stop("The official code needs the ISOweek package: install.packages('ISOweek')")
  start <- as.Date(start); end <- as.Date(end)
  e <- new.env(parent = globalenv())
  e$country <- country; e$country.code <- country_code
  e$wdir <- normalizePath(work_dir, mustWork = TRUE)
  e$indir <- file.path(e$wdir, "data")
  e$outdir <- ensure_dir(file.path(e$wdir, "output"))
  e$start_year <- as.integer(format(start, "%G"))
  e$start_week <- as.integer(format(start, "%V"))
  e$end_year <- as.integer(format(end, "%G"))
  e$end_week <- as.integer(format(end, "%V"))
  e$A_MOMO <- 0; e$WeatherData <- 0; e$population <- FALSE
  e$IArest <- ia_restricted; e$IAlags <- ia_lags; e$ETlags <- et_lags
  run <- function() sys.source(est, envir = e)
  if (verbose) run() else utils::capture.output(suppressWarnings(run()))
  f <- file.path(e$outdir, sprintf("%s_output_v4%s.txt", country,
                                   if (ia_restricted) "_IArestricted" else ""))
  if (!file.exists(f)) stop("FluMOMO produced no output file: ", f)
  data.table::fread(f)
}

#' Daily baseline from the weekly FluMOMO output
#'
#' Spreads the weekly baseline (and, optionally, the weekly attributed deaths)
#' over days with a smooth interpolation, for figures drawn on a daily scale.
#'
#' @param res Output of [flumomo_run()].
#' @param age_group Age group to use (4 = total).
#' @param dates Dates to predict.
#' @param col Column to spread (`"EB"` by default).
#' @return Numeric vector of daily values.
flumomo_daily <- function(res, age_group, dates, col = "EB") {
  # the mask is built outside `[`: inside it, `agegrp` would be the column
  keep <- which(res$agegrp == age_group & !is.na(res$year) & !is.na(res$week) &
                  !is.na(res[[col]]))
  r <- as.data.frame(res)[keep, ]
  # ISO week 53 and week 1 of the next year must not collide
  ws <- ISOweek::ISOweek2date(sprintf("%d-W%02d-1", r$year, r$week))
  o <- order(ws)
  stats::splinefun(as.numeric(ws[o]) + 3, r[[col]][o] / 7)(as.numeric(dates))
}
