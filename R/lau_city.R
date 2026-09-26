# =============================================================================
#  Official city composition.
#
#  Eurostat publishes, on the Local administrative units page, correspondence
#  tables "CITY - LAU" and "FUA - LAU": which municipalities make up each city,
#  greater city and functional urban area. That list is authoritative, so it
#  replaces the geometric matching (city point located in GISCO polygons, best
#  population match) that this project used first.
#
#  The checks stay where they were, in finalize_crosswalk(): they now test the
#  official table instead of testing our inference.
# =============================================================================

#' Read the official CITY - LAU correspondence table
#'
#' Accepts the Eurostat file as .xlsx (needs readxl) or as .csv, downloads it
#' once and caches it. Columns are matched case-insensitively, since the header
#' wording changes between releases: a LAU code, and a city code in either a
#' city or a greater-city column.
#'
#' @param cfg Configuration list (uses `city_table_url`, `data_dir`, `country`).
#' @return data.table `city_code`, `PRO_COM`, `level`, or `NULL` when no table
#'   is configured or it cannot be read.
read_city_lau_table <- function(cfg) {
  url <- cfg$city_table_url
  if (is.null(url) || !nzchar(url)) return(NULL)
  dir <- ensure_dir(file.path(cfg$data_dir, "eurostat"))
  f <- file.path(dir, basename(sub("\\?.*$", "", url)))
  if (!file.exists(f))
    try(utils::download.file(url, f, mode = "wb", quiet = TRUE), silent = TRUE)
  if (!file.exists(f)) {
    message("City-LAU table not downloaded from ", url,
            ": falling back to the geometric matching.")
    return(NULL)
  }
  d <- if (grepl("\\.xlsx?$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      message("Package readxl is needed to read ", f,
              ": falling back to the geometric matching.")
      return(NULL)
    }
    sheets <- readxl::excel_sheets(f)
    # the sheet holding the correspondence is the one with a city column
    tabs <- lapply(sheets, function(s)
      data.table::as.data.table(readxl::read_excel(f, sheet = s, guess_max = 10000)))
    hit <- vapply(tabs, function(x) any(grepl("city", names(x), ignore.case = TRUE)), TRUE)
    if (!any(hit)) return(NULL)
    tabs[[which(hit)[1]]]
  } else data.table::fread(f)
  nm <- function(...) {
    pats <- c(...)
    for (p in pats) {
      hit <- grep(p, names(d), ignore.case = TRUE, value = TRUE)
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }
  c_lau <- nm("^lau.?code$", "^lau.?id$", "^lau$", "gisco.?id")
  c_city <- nm("^city.?code$", "^city.?id$", "^city$")
  c_great <- nm("greater.?city")
  if (is.na(c_lau) || (is.na(c_city) && is.na(c_great))) {
    message("Unexpected columns in ", f, " (", paste(names(d), collapse = ", "),
            "): falling back to the geometric matching.")
    return(NULL)
  }
  out <- list()
  add <- function(col, level) {
    if (is.na(col)) return(NULL)
    x <- data.table::data.table(city_code = toupper(trimws(as.character(d[[col]]))),
                                PRO_COM = pad_procom(d[[c_lau]]), level = level)
    x[nzchar(city_code) & !city_code %in% c("NA", "-") & !is.na(PRO_COM)]
  }
  out[[1]] <- add(c_city, "city"); out[[2]] <- add(c_great, "greater city")
  res <- data.table::rbindlist(Filter(Negate(is.null), out))
  res <- res[substr(city_code, 1, 2) == toupper(cfg$country)]
  if (!nrow(res)) return(NULL)
  unique(res)
}

#' City-comune map for the analysed cities, from the official table
#'
#' Masselot's codes and the table's codes are matched on the six-character
#' stem, so a renumbered trailing digit does not matter. Where a city appears
#' both as a city and as a greater city, the level whose population is closest
#' to Masselot's is kept, which is the same criterion the geometric version
#' used and reproduces his use of greater cities for the large ones.
#'
#' @param tab Output of [read_city_lau_table()].
#' @param metadata Country metadata from [read_erf_bundle()].
#' @param pop data.table `PRO_COM`, `lau_pop` (from the GISCO overlay).
#' @return data.table `URAU_CODE`, `level`, `PRO_COM`, `lau_pop`, or `NULL`.
official_city_map <- function(tab, metadata, pop) {
  if (is.null(tab) || !nrow(tab)) return(NULL)
  tab <- data.table::copy(tab)[, stem := urau_stem(city_code)]
  meta <- data.table::data.table(URAU_CODE = metadata$URAU_CODE,
                                 stem = urau_stem(metadata$URAU_CODE),
                                 pop = as.numeric(metadata$pop))
  x <- merge(tab, meta, by = "stem", allow.cartesian = TRUE)
  if (!nrow(x)) return(NULL)
  x <- merge(x, pop, by = "PRO_COM", all.x = TRUE)
  cand <- x[, .(n_comuni = .N, lau_pop = sum(lau_pop, na.rm = TRUE)),
            by = .(URAU_CODE, level, pop)]
  cand[, score := abs(log(pmax(lau_pop / pop, 1e-6)))]
  best <- cand[, .SD[which.min(score)], by = URAU_CODE][, .(URAU_CODE, level)]
  out <- merge(x, best, by = c("URAU_CODE", "level"))
  out[, .(URAU_CODE, level = paste("official:", level), PRO_COM, LAU_NAME = NA_character_,
          lau_pop)]
}
