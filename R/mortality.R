#' Read the Istat daily deaths file for selected comuni
#'
#' Reads `comuni_giornaliero_*.csv` (one row per comune x age class x
#' day-of-year, one column `T_yy` per year), keeps the requested comuni,
#' merges Istat 5-year age classes into Masselot's age groups and drops ages
#' below 20. Values marked `n.d.` (not yet available) are kept as `NA` so the
#' end of the data can be detected; days with no deaths are absent from the
#' file and are filled with zeros later, in [aggregate_mortality()].
#'
#' @param path Path to the Istat CSV.
#' @param procom Comune codes (6 digits) to keep.
#' @return List with `long` (data.table `PRO_COM`, `agegroup`, `date`,
#'   `deaths`) and `age_check` (deaths by `CL_ETA`, to verify the coding).
read_istat <- function(path, procom) {
  hdr <- names(data.table::fread(path, nrows = 0))
  tcols <- grep("^T_[0-9]{2}$", hdr, value = TRUE)
  d <- data.table::fread(path, select = c("COD_PROVCOM", "CL_ETA", "GE", tcols),
                         colClasses = list(character = c("COD_PROVCOM", "GE")),
                         na.strings = c("", "NA", "n.d."))
  d[, PRO_COM := pad_procom(COD_PROVCOM)]
  d <- d[PRO_COM %in% procom]
  for (tc in tcols) data.table::set(d, j = tc, value = suppressWarnings(as.integer(d[[tc]])))
  age_check <- d[, lapply(.SD, sum, na.rm = TRUE), by = CL_ETA, .SDcols = tcols]
  age_check <- age_check[, .(CL_ETA, deaths = rowSums(.SD)), .SDcols = tcols][order(CL_ETA)]
  age_check[, agegroup := istat_age_group(CL_ETA)]
  d[, agegroup := istat_age_group(CL_ETA)]
  d <- d[!is.na(agegroup)]
  long <- data.table::melt(d, id.vars = c("PRO_COM", "agegroup", "GE"),
                           measure.vars = tcols, variable.name = "yy",
                           value.name = "deaths")
  long[, GE := sprintf("%04d", as.integer(GE))]
  long[, date := as.Date(sprintf("%d-%s-%s", 2000L + as.integer(sub("T_", "", yy)),
                                 substr(GE, 1, 2), substr(GE, 3, 4)),
                         optional = TRUE)]
  long <- long[!is.na(date)]   # 29 February in non-leap years
  long <- long[, .(deaths = if (all(is.na(deaths))) NA_integer_ else sum(deaths, na.rm = TRUE)),
               by = .(PRO_COM, agegroup, date)]
  list(long = long, age_check = age_check)
}

#' Last day with Istat data
#'
#' @param istat Output of [read_istat()].
#' @param override Optional date that takes precedence.
#' @return A Date.
detect_data_end <- function(istat, override = NULL) {
  if (!is.null(override)) return(as.Date(override))
  max(istat$long[!is.na(deaths), date])
}

#' Daily deaths by city and age group
#'
#' Sums comuni into Urban Audit cities and fills missing days with zeros
#' (the Istat file omits days without deaths).
#'
#' @param istat Output of [read_istat()].
#' @param map City-comune table (`finalize_crosswalk()$map`).
#' @param start First date.
#' @param end Last date.
#' @return data.table `URAU_CODE`, `agegroup`, `date`, `deaths`.
aggregate_mortality <- function(istat, map, start, end) {
  x <- merge(istat$long[date >= start & date <= end],
             map[, .(URAU_CODE, PRO_COM)], by = "PRO_COM", allow.cartesian = TRUE)
  x <- x[, .(deaths = sum(deaths, na.rm = TRUE)), by = .(URAU_CODE, agegroup, date)]
  grid <- data.table::CJ(URAU_CODE = unique(map$URAU_CODE), agegroup = AGE_GROUPS,
                         date = seq(start, end, by = "day"))
  out <- merge(grid, x, by = c("URAU_CODE", "agegroup", "date"), all.x = TRUE)
  out[is.na(deaths), deaths := 0L]
  data.table::setkey(out, URAU_CODE, agegroup, date)
  out
}

#' Compare the app baseline with observed Istat deaths
#'
#' The app converts attributable fractions into deaths with a baseline built
#' from Masselot's metadata. Since attributable deaths scale linearly with the
#' baseline, their ratio to observed mean daily deaths is a direct bias factor.
#'
#' @param mort Output of [aggregate_mortality()].
#' @param erf Output of [read_erf_bundle()].
#' @param start,end Period.
#' @return data.table by city and age group.
audit_baseline <- function(mort, erf, start, end) {
  obs <- mort[date >= start & date <= end, .(istat = mean(deaths)), by = .(URAU_CODE, agegroup)]
  obs[, app := mapply(function(cd, a) app_baseline(erf, cd, a), URAU_CODE, agegroup)]
  obs[, ratio := app / istat]
  obs
}
