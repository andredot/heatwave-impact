# =============================================================================
#  Official composition of the Italian Cities (Istat, Cities 2021).
#
#  Istat publishes, in the SITUAS portal (Sistema Informativo sulle Unita
#  Territoriali e sugli Indirizzi), the report
#
#      "Composizione delle Citta (Cities) 2021"
#
#  which gives, for every City and Greater City defined by EU Regulation
#  2019/1130, its Urban Audit code and the comuni that compose it with their
#  Istat codes: 89 Cities, 11 of them Greater Cities, 298 comuni in all.
#
#  That is exactly the crosswalk this project had to infer geometrically
#  (city point located in GISCO polygons, comuni by area share, candidate
#  chosen on population). Two columns of the report replace all of it:
#
#      Codice Citta (City/Greater City) 2021  ->  Masselot's URAU_CODE
#      Codice Comune (alfanumerico)           ->  Istat PRO_COM in the
#                                                 mortality file
#
#  so no name matching and no population heuristic are involved any more.
#  The geometric matching is still computed, and [finalize_crosswalk()] now
#  reports comune by comune how far the two agree.
#
#  The definition is fixed (2021 reference date) and the file is 28 kB, so a
#  copy is kept in the repository at reference/istat_cities_2021.csv and a
#  fresh clone needs no download. To refresh it (e.g. for a future Cities
#  release), open
#
#      https://situas.istat.it/web/            Report > Unita territoriali >
#      "Composizione delle Citta (Cities) 2021"
#
#  set the reference date, export as CSV and save it over that file, or drop
#  the export anywhere under data/ and the pipeline will find it.
#
#  Note on FUAs: the same portal publishes the Functional Urban Areas, which
#  add the commuting belt around the city. A FUA is NOT interchangeable with
#  a City: it would put a much larger population under the same curve, and
#  Masselot's unit is the City. The FUA report is therefore not used.
# =============================================================================

#' Locate the Istat Cities composition file
#'
#' Looks for `cfg$istat_cities_file` first, then for any CSV export of the
#' SITUAS report left under `cfg$data_dir` (its default name starts with
#' "Composizione"). Errors rather than falling back silently, since this file
#' now defines the geography of the whole analysis.
#'
#' @param cfg Configuration list (uses `istat_cities_file`, `data_dir`).
#' @return Path to the file.
locate_istat_cities <- function(cfg) {
  f <- cfg$istat_cities_file
  if (!is.null(f) && nzchar(f) && file.exists(f)) return(normalizePath(f, winslash = "/"))
  hits <- list.files(cfg$data_dir, pattern = "^Composizione.*\\.csv$",
                     full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(hits)) return(normalizePath(hits[1], winslash = "/"))
  stop("The Istat Cities composition is missing.\n",
       "  expected at: ", f, "\n",
       "  it is shipped with the repository: restore it with\n",
       "      git checkout -- ", f, "\n",
       "  or export it again from https://situas.istat.it/web/ (report\n",
       "  \"Composizione delle Citta (Cities) 2021\") and save it there.\n",
       "  Setting cfg$city_source to \"geometry\" restores the geometric\n",
       "  matching, but this file is still read for the agreement check.",
       call. = FALSE)
}

#' Read the composition of the Cities: one row per comune
#'
#' Columns are identified by their *content*, not by their header, because the
#' SITUAS export labels them in Italian and renames them between releases:
#' the Urban Audit code is the column whose values look like `IT001C`, the
#' comune code the digit-only column that identifies the row, the city name
#' the column with one value per city code, the comune name the remaining
#' near-unique text column. What was detected is printed, so a change of
#' format is visible instead of silent.
#'
#' @param file Path from [locate_istat_cities()].
#' @return data.table with one row per comune:
#'   `city_code` (Urban Audit code), `city` (name, without the
#'   "(greater city)" suffix), `level` (`"city"` / `"greater city"`),
#'   `PRO_COM` (six-digit Istat comune code), `comune`, and `region`,
#'   `province`, `nuts3` when present.
read_istat_cities <- function(file) {
  d <- data.table::fread(file, colClasses = "character", encoding = "UTF-8",
                         na.strings = c("", "NA"), showProgress = FALSE)
  if (!nrow(d) || ncol(d) < 3)
    stop("Nothing to read in ", file, ": ", nrow(d), " rows, ", ncol(d), " columns.",
         call. = FALSE)
  cols <- names(d)
  val <- lapply(cols, function(cn) {
    x <- trimws(as.character(d[[cn]]))
    x[is.na(x)] <- ""
    x
  })
  names(val) <- cols
  n <- nrow(d)
  nd <- vapply(val, function(x) length(unique(x[nzchar(x)])), 0L)
  share <- function(x, p) {
    x <- x[nzchar(x)]
    if (!length(x)) 0 else mean(grepl(p, x))
  }
  byname <- function(...) {
    for (p in c(...)) {
      hit <- grep(p, cols, ignore.case = TRUE, value = TRUE)
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }

  # Urban Audit code: two letters, three digits, a letter, an optional version
  is_urau <- vapply(val, function(x) share(x, "^[A-Za-z]{2}[0-9]{3}[A-Za-z][0-9]?$") > 0.9, TRUE)
  c_code <- names(which(is_urau))[1]
  # Comune code: digits only, four to six of them, one value per row. Both the
  # alphanumeric (zero-padded) and the numeric column qualify; prefer the
  # padded one, though pad_procom() makes them equivalent.
  is_num <- vapply(val, function(x) share(x, "^[0-9]{4,6}$") > 0.99, TRUE)
  cand <- cols[is_num & nd == n]
  c_procom <- if (!length(cand)) NA_character_ else
    cand[order(-vapply(val[cand], function(x) mean(nchar(x) == 6), 0))][1]
  if (is.na(c_code) || is.na(c_procom))
    stop("Unexpected format in ", file, ": no ",
         if (is.na(c_code)) "Urban Audit code" else "comune code",
         " column found among: ", paste(cols, collapse = ", "), call. = FALSE)

  # City name: one value per city code, as many values as there are cities.
  ncity <- nd[[c_code]]
  one_per_city <- vapply(cols, function(cn) {
    if (nd[[cn]] != ncity) return(FALSE)
    all(tapply(val[[cn]], val[[c_code]], function(z) length(unique(z))) == 1L)
  }, TRUE)
  c_city <- setdiff(cols[one_per_city], c_code)[1]
  if (is.na(c_city)) c_city <- byname("citt", "city")
  # Comune name: text, near-unique, not one of the code columns.
  c_comune <- setdiff(cols[nd >= 0.9 * n & !is_num & !is_urau], c(c_code, c_city))[1]
  if (is.na(c_comune)) c_comune <- byname("^comune", "den.*com")

  squish <- function(x) trimws(gsub("[[:space:]]+", " ", x))
  gc_pat <- "\\(greater[[:space:]]*city\\)"
  nmc <- val[[c_city]]
  out <- data.table::data.table(
    city_code = toupper(val[[c_code]]),
    city      = squish(gsub(gc_pat, "", nmc, ignore.case = TRUE)),
    level     = data.table::fifelse(grepl(gc_pat, nmc, ignore.case = TRUE),
                                    "greater city", "city"),
    PRO_COM   = pad_procom(val[[c_procom]]),
    comune    = if (is.na(c_comune)) NA_character_ else squish(val[[c_comune]]))
  # Region, province and NUTS3 names, when the export carries them: the first
  # matching column that holds names rather than codes.
  for (extra in list(c(nm = "region", pat = "regione"),
                     c(nm = "province", pat = "provincia"),
                     c(nm = "nuts3", pat = "nuts3"))) {
    hits <- setdiff(grep(extra[["pat"]], cols, ignore.case = TRUE, value = TRUE),
                    c(c_code, c_city, c_procom, c_comune))
    hits <- hits[vapply(val[hits], function(x) share(x, "^[0-9]+$") < 0.5, TRUE)]
    named <- grep("^cod", hits, ignore.case = TRUE, invert = TRUE, value = TRUE)
    if (length(named)) hits <- named
    if (length(hits)) out[[extra[["nm"]]]] <- squish(val[[hits[1]]])
  }
  out <- unique(out[nzchar(city_code) & !is.na(PRO_COM)])

  message("Istat Cities composition: ", basename(file), "\n  columns used: city code = ",
          c_code, ", city = ", c_city, ", comune code = ", c_procom, ", comune = ",
          c_comune, "\n  ", nrow(out), " comuni in ",
          data.table::uniqueN(out$city_code), " cities (",
          data.table::uniqueN(out[level == "greater city", city_code]), " greater cities)")
  dup <- out[, .N, by = PRO_COM][N > 1L]
  if (nrow(dup))
    warning("A comune belongs to more than one city in ", basename(file), ": ",
            paste(dup$PRO_COM, collapse = ", "), call. = FALSE)
  lvl <- out[, data.table::uniqueN(level), by = city_code][V1 > 1L]
  if (nrow(lvl))
    warning("A city appears at both levels in ", basename(file), ": ",
            paste(lvl$city_code, collapse = ", "), call. = FALSE)
  attr(out, "file") <- file
  attr(out, "source") <-
    "Istat SITUAS, \"Composizione delle Citta (Cities) 2021\" (EU Regulation 2019/1130)"
  out[]
}

#' City-comune map for the analysed cities, from the Istat composition
#'
#' The report carries Urban Audit codes, so Masselot's cities are matched on
#' the six-character code stem (`IT001C` vs `IT001C1`), with no name matching
#' and no choice of level: the report states whether a city is a City or a
#' Greater City. GISCO is used only for the comune population, which the
#' alignment checks need.
#'
#' @param tab Output of [read_istat_cities()].
#' @param metadata Country metadata from [read_erf_bundle()].
#' @param pop data.table `PRO_COM`, `lau_pop` (from [read_lau_population()]),
#'   or `NULL`.
#' @return data.table `URAU_CODE`, `level`, `PRO_COM`, `LAU_NAME`, `lau_pop`,
#'   with the source in `attr(, "source")`.
istat_city_map <- function(tab, metadata, pop = NULL) {
  if (is.null(tab) || !nrow(tab)) return(NULL)
  x <- data.table::copy(tab)[, stem := urau_stem(city_code)]
  meta <- data.table::data.table(URAU_CODE = metadata$URAU_CODE,
                                 stem = urau_stem(metadata$URAU_CODE))
  x <- merge(x, unique(meta), by = "stem", allow.cartesian = TRUE)
  if (!nrow(x))
    stop("No Urban Audit code of the Istat composition matches the metadata; ",
         "first codes in the file: ",
         paste(utils::head(unique(tab$city_code), 5), collapse = ", "), call. = FALSE)
  x[, LAU_NAME := comune]
  if (!is.null(pop) && nrow(pop)) {
    p <- unique(data.table::as.data.table(pop)[, .(PRO_COM, lau_pop)])
    x <- merge(x, p, by = "PRO_COM", all.x = TRUE)
  } else x[, lau_pop := NA_real_]
  out <- unique(x[, .(URAU_CODE, level = paste("Istat 2021:", level), PRO_COM,
                      LAU_NAME, lau_pop)])
  attr(out, "source") <- attr(tab, "source")
  out[order(URAU_CODE, PRO_COM)]
}

#' Coverage of Masselot's cities by the Istat composition
#'
#' @param tab Output of [read_istat_cities()].
#' @param metadata Country metadata from [read_erf_bundle()].
#' @return List with `cities` (one row per city of the metadata: whether the
#'   Istat report covers it, at which level, with how many comuni) and
#'   `unmatched` (cities of the report absent from the metadata).
istat_city_coverage <- function(tab, metadata) {
  if (is.null(tab) || !nrow(tab)) return(NULL)
  nm <- if ("name" %in% names(metadata)) metadata$name else metadata$URAU_NAME
  meta <- data.table::data.table(URAU_CODE = metadata$URAU_CODE, name = nm,
                                 stem = urau_stem(metadata$URAU_CODE),
                                 pop = as.numeric(metadata$pop))
  s <- data.table::copy(tab)[, stem := urau_stem(city_code)][
    , .(istat_code = city_code[1], istat_city = city[1], level = level[1],
        n_comuni = .N), by = stem]
  m <- merge(meta, s, by = "stem", all.x = TRUE)
  m[, matched := !is.na(n_comuni)]
  list(cities = m[order(-pop)][, !"stem"],
       unmatched = s[!stem %in% meta$stem][order(istat_code)][, !"stem"])
}
