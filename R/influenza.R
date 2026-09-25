# =============================================================================
#  Weekly influenza activity (IA) for FluMOMO.
#
#  FluMOMO needs one consistent indicator over the whole fitted period. The
#  authors recommend the Goldstein index, the product of the consultation rate
#  for influenza-like illness and the share of sentinel samples positive for
#  influenza, because it is more conservative and represents circulation in
#  the community better than either part alone.
#
#  Sources, in the order they are used:
#   1. ECDC ERVISS weekly files (machine readable, but only recent seasons);
#   2. WHO FluNet/FluID, or a CSV exported from them, for the earlier years.
# =============================================================================

ERVISS_BASE <- paste0("https://raw.githubusercontent.com/EU-ECDC/",
                      "Respiratory_viruses_weekly_data/refs/heads/main/data/")

#' Year and week from an ISO "2026-W09" label
#'
#' @param x Character vector.
#' @return data.table with `year` and `week`.
split_yearweek <- function(x) {
  data.table::data.table(year = as.integer(substr(x, 1, 4)),
                         week = as.integer(sub(".*W", "", x)))
}

#' Influenza activity from ECDC ERVISS
#'
#' Downloads the consultation rates and the sentinel positivity published by
#' ECDC and returns the Goldstein index for one country. From season 2025-26
#' Italy reports acute respiratory infections instead of influenza-like
#' illness; multiplying by influenza positivity removes most of the
#' non-influenza traffic that the change brought in, but the break is real and
#' is flagged in the returned table.
#'
#' @param country Country name as used by ERVISS (e.g. `"Italy"`).
#' @param cache_dir Where the raw files are cached.
#' @return data.table `year`, `week`, `IA`, `rate_indicator`, `source`.
fetch_erviss_ia <- function(country, cache_dir) {
  dir <- ensure_dir(file.path(cache_dir, "erviss"))
  get <- function(f) {
    p <- file.path(dir, f)
    stale <- !file.exists(p) ||
      difftime(Sys.time(), file.mtime(p), units = "days") > 3
    if (stale) {
      ok <- tryCatch({
        utils::download.file(paste0(ERVISS_BASE, f), p, quiet = TRUE, mode = "wb")
        file.exists(p) && file.size(p) > 1000
      }, error = function(e) {
        message("Could not download ", f, ": ", conditionMessage(e)); FALSE
      })
      if (!ok && !file.exists(p))
        stop("Download of ", paste0(ERVISS_BASE, f), " failed. Save the file by hand ",
             "in ", dir, " (behind a proxy, set options(download.file.method), or ",
             "use a browser), then run again.", call. = FALSE)
    }
    d <- tryCatch(data.table::fread(p), error = function(e)
      stop("Cannot read ", p, " (", conditionMessage(e),
           "). Delete it and run again.", call. = FALSE))
    if (!nrow(d)) stop("The file ", p, " is empty; delete it and run again.", call. = FALSE)
    d
  }
  rates <- get("ILIARIRates.csv")
  pos <- get("sentinelTestsDetectionsPositivity.csv")
  need <- function(d, cols, f) {
    miss <- setdiff(cols, names(d))
    if (length(miss))
      stop("Columns ", paste(miss, collapse = ", "), " missing from ", f,
           ": the ERVISS format has changed, R/influenza.R needs updating.",
           call. = FALSE)
  }
  need(rates, c("countryname", "yearweek", "indicator", "age", "value"), "ILIARIRates.csv")
  need(pos, c("countryname", "yearweek", "pathogen", "indicator", "age", "value",
              "survtype"), "sentinelTestsDetectionsPositivity.csv")
  if (!country %in% rates$countryname)
    stop("Country '", country, "' not in the ERVISS files. Available: ",
         paste(utils::head(sort(unique(rates$countryname)), 10), collapse = ", "),
         ", ...", call. = FALSE)
  r <- rates[countryname == country & age == "total" &
               indicator %in% c("ILIconsultationrate", "ARIconsultationrate"),
             .(yearweek, indicator, rate = as.numeric(value))]
  p <- pos[countryname == country & age == "total" & pathogen == "Influenza" &
             indicator == "positivity" & survtype == "primary care sentinel",
           .(yearweek, positivity = as.numeric(value))]
  if (!nrow(r) || !nrow(p)) return(data.table::data.table())
  x <- merge(r, p, by = "yearweek")
  x <- cbind(split_yearweek(x$yearweek), x[, .(rate_indicator = indicator,
                                               IA = rate * positivity / 100)])
  x[, source := "ERVISS"][order(year, week)]
}

#' Influenza activity from WHO FluNet (or an export of it)
#'
#' Used for the years ERVISS does not cover. Tries the public WHO API and
#' falls back to a local CSV exported from FluNet, whose columns are matched
#' case-insensitively: an ISO year and week, the number of influenza-positive
#' specimens and the number of specimens processed. Without consultation rates
#' the indicator is the positivity share alone, which FluMOMO accepts as an
#' influenza activity indicator, less conservative than the Goldstein index.
#'
#' @param country_code ISO3 code (e.g. `"ITA"`).
#' @param file Optional local FluNet export.
#' @param cache_dir Where the download is cached.
#' @return data.table `year`, `week`, `IA`, `source` (empty if unavailable).
fetch_flunet_ia <- function(country_code, file = NULL, cache_dir = "data/cache") {
  d <- NULL
  if (!is.null(file) && file.exists(file)) {
    d <- data.table::fread(file)
  } else {
    p <- file.path(ensure_dir(cache_dir), sprintf("flunet_%s.csv", country_code))
    url <- sprintf(paste0("https://xmart-api-public.who.int/FLUMART/VIW_FNT?",
                          "$filter=COUNTRY_CODE%%20eq%%20'%s'&$format=csv"), country_code)
    ok <- !inherits(try(utils::download.file(url, p, quiet = TRUE, mode = "wb"),
                        silent = TRUE), "try-error")
    if (ok && file.exists(p) && file.size(p) > 1000) d <- data.table::fread(p)
  }
  if (is.null(d) || !nrow(d)) return(data.table::data.table())
  nm <- function(...) {
    hit <- intersect(tolower(c(...)), tolower(names(d)))
    if (length(hit)) names(d)[match(hit[1], tolower(names(d)))] else NA_character_
  }
  c_y <- nm("iso_year", "mmwr_year", "year"); c_w <- nm("iso_week", "mmwr_week", "week")
  c_pos <- nm("inf_all", "influenza_all", "inf_a", "all_inf")
  c_tot <- nm("spec_processed_nb", "spec_received_nb", "specimens_processed")
  if (any(is.na(c(c_y, c_w, c_pos, c_tot)))) return(data.table::data.table())
  x <- d[, .(year = as.integer(get(c_y)), week = as.integer(get(c_w)),
             pos = as.numeric(get(c_pos)), tot = as.numeric(get(c_tot)))]
  x <- x[!is.na(year) & !is.na(week) & tot > 0,
         .(IA = 100 * sum(pos, na.rm = TRUE) / sum(tot, na.rm = TRUE)), by = .(year, week)]
  x[, source := "FluNet"][order(year, week)]
}

#' Influenza settings, with defaults
#'
#' Keeps the pipeline working when `config.R` predates the `influenza` block,
#' and fails with a readable message when a value is unusable.
#'
#' @param cfg Configuration list.
#' @return The settings list, completed with defaults.
influenza_settings <- function(cfg) {
  def <- list(country = "Italy", country_code = "ITA",
              file = "data/raw/influenza_activity_IT.csv",
              flunet_file = "data/raw/flunet_italy.csv",
              use_existing = FALSE, allow_zero = FALSE)
  ic <- cfg$influenza
  if (is.null(ic)) {
    warning("config.R has no `influenza` settings: using the defaults (country ",
            def$country, "). Update config.R to choose another country.",
            call. = FALSE)
    ic <- list()
  }
  miss <- setdiff(names(def), names(ic))
  ic[miss] <- def[miss]
  if (!is.character(ic$country) || length(ic$country) != 1 || !nzchar(ic$country))
    stop("influenza$country in config.R must be a single country name, ",
         "as spelled by ERVISS (e.g. \"Italy\").", call. = FALSE)
  ic
}

#' Build the influenza activity file used by FluMOMO
#'
#' ERVISS is used wherever it reaches; earlier weeks come from FluNet, rescaled
#' so that the two series have the same mean over their overlap (they measure
#' activity on different scales, and FluMOMO only needs a consistent
#' indicator). Weeks with no data are left out and treated as zero by the
#' official code.
#'
#' @param cfg Configuration list (uses `influenza`).
#' @param path Output CSV.
#' @return `path`.
build_influenza_activity <- function(cfg, path = NULL) {
  ic <- influenza_settings(cfg)
  path <- path %||% ic$file
  # 1. a file prepared by hand always wins
  if (isTRUE(ic$use_existing) && file.exists(path)) {
    message("Using the existing influenza file ", path)
    return(path)
  }
  a <- tryCatch(fetch_erviss_ia(ic$country, cfg$cache_dir), error = function(e) {
    message("ERVISS not available: ", conditionMessage(e)); data.table::data.table()
  })
  b <- tryCatch(fetch_flunet_ia(ic$country_code, ic$flunet_file, cfg$cache_dir),
                error = function(e) data.table::data.table())
  if (!nrow(a) && !nrow(b)) {
    if (!isTRUE(ic$allow_zero))
      stop("No influenza activity could be obtained. Either fix the download, or put ",
           "a file with columns year, week, IA at ", path, " and set ",
           "influenza$use_existing = TRUE in config.R, or set ",
           "influenza$allow_zero = TRUE to run FluMOMO without the influenza term ",
           "(the baseline will then contain influenza-related mortality).",
           call. = FALSE)
    warning("Running with IA = 0: influenza stays inside the baseline.", call. = FALSE)
    out <- data.table::data.table(year = integer(), week = integer(),
                                  IA = numeric(), source = character())
    ensure_dir(dirname(path)); data.table::fwrite(out, path)
    return(path)
  }
  out <- if (!nrow(b)) a else if (!nrow(a)) b else {
    ov <- merge(a[, .(year, week, ia_e = IA)], b[, .(year, week, ia_f = IA)],
                by = c("year", "week"))
    k <- if (nrow(ov) >= 26 && mean(ov$ia_f) > 0) mean(ov$ia_e) / mean(ov$ia_f) else 1
    rbind(a[, .(year, week, IA, source)],
          b[!a, on = .(year, week)][, .(year, week, IA = IA * k,
                                        source = "FluNet (rescaled)")])
  }
  out <- out[order(year, week)]
  ensure_dir(dirname(path))
  data.table::fwrite(out, path)
  message("Influenza activity: ", nrow(out), " weeks, ",
          paste(sort(unique(out$source)), collapse = " + "), " -> ", path)
  path
}

`%||%` <- function(x, y) if (is.null(x)) y else x
