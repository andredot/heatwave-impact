# =============================================================================
#  One region in full.
#
#  The validation tests every city that carries enough deaths for Masselot's
#  first stage (5,000 over the period). For operational use a region needs a
#  curve for *every* one of its cities, including the small ones, and the
#  second stage is what makes that possible: a city whose own fit is noisy
#  enters the multivariate random-effects meta-analysis with a large variance,
#  so its best linear unbiased prediction sits close to the pooled Italian
#  curve. The city contributes what its data support and borrows the rest.
#
#  This file only selects the region's cities; the fitting is done by
#  run_curve_tests(..., min_fit_deaths = 0) and the shrinkage by
#  blup_curves(..., extra_fits = , keep = ), so the region's curves come out
#  of exactly the same code as the national ones.
# =============================================================================

#' Region of each analysed city
#'
#' Taken from the Istat composition, which gives the region of every comune;
#' a city is assigned the region holding most of its comuni.
#'
#' @param istat_cities Output of [read_istat_cities()].
#' @param map City-comune map, e.g. `crosswalk$map`.
#' @return data.table `URAU_CODE`, `region`, or `NULL` when the composition
#'   carries no region column.
city_regions <- function(istat_cities, map) {
  if (is.null(istat_cities) || !"region" %in% names(istat_cities)) return(NULL)
  r <- unique(data.table::as.data.table(istat_cities)[, .(PRO_COM, region)])
  x <- merge(data.table::as.data.table(map)[, .(URAU_CODE, PRO_COM)], r, by = "PRO_COM")
  if (!nrow(x)) return(NULL)
  x <- x[, .N, by = .(URAU_CODE, region)]
  x[order(URAU_CODE, -N), .SD[1L], by = URAU_CODE][, .(URAU_CODE, region)]
}

#' Urban Audit codes of one region's cities
#'
#' @param istat_cities Output of [read_istat_cities()].
#' @param map City-comune map, e.g. `crosswalk$map`.
#' @param region Region name as the Istat composition writes it (e.g.
#'   `"Lombardia"`); matched ignoring case, accents and punctuation.
#' @return Character vector of codes (possibly empty).
region_city_codes <- function(istat_cities, map, region) {
  cr <- city_regions(istat_cities, map)
  if (is.null(cr) || is.null(region) || !nzchar(region)) return(character())
  # `region` is both the argument and a column name, so the comparison is
  # computed outside the brackets, where it means the argument.
  want <- norm_name(region)
  hit <- norm_name(cr$region) == want
  if (!any(hit)) hit <- grepl(want, norm_name(cr$region), fixed = TRUE)
  out <- cr$URAU_CODE[which(hit)]
  if (!length(out))
    message("No city of region \"", region, "\" among the analysed cities. Regions present: ",
            paste(sort(unique(cr$region)), collapse = ", "))
  out
}

#' Table of a region's cities: published, own and BLUP curves side by side
#'
#' @param blup Output of [blup_curves()] restricted to the region.
#' @param tests Output of [run_curve_tests()] for the region (for the fitting
#'   status of cities whose model did not converge).
#' @param metadata Country metadata from [read_erf_bundle()].
#' @return data.table, one row per city, ordered by deaths.
region_curve_table <- function(blup, tests, metadata) {
  if (is.null(tests) || !nrow(tests$results)) return(NULL)
  meta <- data.table::as.data.table(metadata)[, .(URAU_CODE, name, inmcc,
                                                  pop = as.numeric(pop))]
  st <- tests$results[level == "all20", .(URAU_CODE, status)]
  out <- merge(meta, st, by = "URAU_CODE")
  res <- if (is.null(blup)) NULL else data.table::copy(blup$results)
  if (!is.null(res)) res[, intersect(c("name", "inmcc"), names(res)) := NULL]
  if (!is.null(res) && nrow(res)) out <- merge(out, res, by = "URAU_CODE", all.x = TRUE)
  # the columns the report expects exist even when nothing could be fitted
  for (cl in c("deaths", "rr99_pub", "rr99_own", "rr99_own_lo", "rr99_own_hi",
               "rr99_blup", "rr99_blup_lo", "rr99_blup_hi", "d_lr99", "se_d_lr99"))
    if (!cl %in% names(out)) data.table::set(out, j = cl, value = NA_real_)
  if (!"below_threshold" %in% names(out))
    data.table::set(out, j = "below_threshold", value = NA)
  out[order(-pop)]
}
