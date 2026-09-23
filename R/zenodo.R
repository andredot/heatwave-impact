ZENODO_URL <- "https://zenodo.org/records/10288665/files/%s?download=1"

#' Download the Masselot et al. (2023) Zenodo files
#'
#' Files are cached in `dir` and downloaded only once.
#'
#' @param dir Destination directory.
#' @param coef_simu Also download `coef_simu.csv` (491 MB)?
#' @return Named character vector of file paths.
download_zenodo <- function(dir, coef_simu = FALSE) {
  files <- c("coefs.csv", "vcov.csv", "tmean_distribution.csv", "metadata.csv",
             "additional_data.zip", "results.zip")
  if (coef_simu) files <- c(files, "coef_simu.csv")
  paths <- vapply(files, function(f)
    download_if_missing(sprintf(ZENODO_URL, f), file.path(dir, f)), character(1))
  setNames(paths, sub("\\.(csv|zip)$", "", files))
}

#' Pick a downloaded file by name
#'
#' File targets lose vector names, so files are matched on their base name.
#'
#' @param paths Character vector of paths.
#' @param name Base name without extension (e.g. `"metadata"`).
#' @return The matching path.
zpath <- function(paths, name) {
  hit <- paths[sub("\\.(csv|zip|gpkg)$", "", basename(paths)) == name]
  if (!length(hit)) stop("File '", name, "' not among: ", paste(basename(paths), collapse = ", "))
  hit[1]
}

#' Read one CSV stored inside a zip archive
#'
#' @param zip Path to the archive.
#' @param name File name inside the archive (matched on the base name).
#' @return A data.table.
read_from_zip <- function(zip, name) {
  inside <- utils::unzip(zip, list = TRUE)$Name
  hit <- inside[basename(inside) == name]
  if (!length(hit)) stop(name, " not found in ", zip)
  ex <- tempfile()
  utils::unzip(zip, files = hit[1], exdir = ex)
  on.exit(unlink(ex, recursive = TRUE))
  data.table::fread(file.path(ex, hit[1]))
}

#' Read the published ERFs for one country
#'
#' Returns coefficients, covariance matrices, temperature distributions and
#' metadata of all cities of `country`, plus the coefficient draws if
#' `coef_simu.csv` was downloaded.
#'
#' @param paths Output of [download_zenodo()].
#' @param country Two-letter country code (`CNTR_CODE` in the metadata).
#' @return A list with elements `metadata`, `coefs`, `vcov`, `tdist`, `simu`
#'   (`NULL` when not available).
read_erf_bundle <- function(paths, country) {
  meta <- data.table::fread(zpath(paths, "metadata"), encoding = "UTF-8")
  # The row filter is computed OUTSIDE `[`: metadata.csv has a column called
  # `country`, which inside `meta[...]` would mask the function argument.
  ctry <- toupper(country)
  is_ctry <- toupper(trimws(meta$CNTR_CODE)) == ctry |
    substr(toupper(trimws(meta$URAU_CODE)), 1, 2) == ctry
  meta <- meta[which(is_ctry)]
  if (!nrow(meta))
    stop("No city of country '", ctry, "' in metadata.csv.")
  if (!"inmcc" %in% names(meta)) meta[, inmcc := NA]
  if (!"LABEL" %in% names(meta)) meta[, LABEL := NA_character_]
  meta[, name := data.table::fifelse(!is.na(LABEL) & nzchar(LABEL), LABEL, URAU_NAME)]
  codes <- meta$URAU_CODE
  keep  <- function(d) d[URAU_CODE %in% codes]
  simu <- NULL
  if (any(basename(paths) == "coef_simu.csv")) {
    simu <- keep(data.table::fread(zpath(paths, "coef_simu")))
    data.table::setkey(simu, URAU_CODE, agegroup, sim)
  }
  list(metadata = meta,
       coefs = keep(data.table::fread(zpath(paths, "coefs"))),
       vcov  = keep(data.table::fread(zpath(paths, "vcov"))),
       tdist = keep(data.table::fread(zpath(paths, "tmean_distribution"))),
       simu  = simu)
}

#' Read Masselot's ERA5-Land city series (1990-2019)
#'
#' `era5series.csv` (in `additional_data.zip`) holds the daily mean
#' temperature averaged over the ERA5-Land cells inside each city boundary:
#' the exact exposure the ERFs were built on.
#'
#' @param paths Output of [download_zenodo()].
#' @param codes Urban Audit codes to keep.
#' @return data.table with `URAU_CODE`, `date`, `tmean`.
read_era5_masselot <- function(paths, codes) {
  d <- read_from_zip(zpath(paths, "additional_data"), "era5series.csv")
  need <- c("URAU_CODE", "date", "era5landtmean")
  if (!all(need %in% names(d)))
    stop("Unexpected columns in era5series.csv: ", paste(names(d), collapse = ", "))
  out <- d[URAU_CODE %in% codes,
           .(URAU_CODE = as.character(URAU_CODE), date = as.Date(date),
             tmean = as.numeric(era5landtmean))]
  if (!nrow(out))
    stop("era5series.csv has no rows for the selected cities (",
         paste(utils::head(codes), collapse = ", "), ")")
  data.table::setkey(out, URAU_CODE, date)
  out
}

#' Read the published city x age results
#'
#' `results.zip` holds the numbers reported in the paper: minimum mortality
#' temperature and percentile, relative risks at the 1st and 99th
#' percentiles, annual deaths and annual excess deaths attributed to heat and
#' cold. They are the reference for [check_reconstruction()].
#'
#' @param paths Output of [download_zenodo()].
#' @param codes Urban Audit codes to keep.
#' @return data.table of published results by city and age group.
read_published_results <- function(paths, codes) {
  d <- read_from_zip(zpath(paths, "results"), "cityage.csv")
  d[URAU_CODE %in% codes]
}

#' Read Masselot's first-stage estimates for the country's MCC cities
#'
#' `stage1res.csv` holds the city-specific curves fitted on observed
#' mortality (MCC data; Italy 2001-2010). The city is identified by the MCC
#' code, mapped to Urban Audit through the `mcc_code` metadata column.
#'
#' @param paths Output of [download_zenodo()].
#' @param metadata Country metadata from [read_erf_bundle()].
#' @return data.table of first-stage results with `URAU_CODE` added
#'   (empty if no city of the country was in MCC).
read_stage1 <- function(paths, metadata) {
  s1 <- read_from_zip(zpath(paths, "additional_data"), "stage1res.csv")
  if (!"mcc_code" %in% names(metadata)) return(s1[0])
  s1[, URAU_CODE := metadata$URAU_CODE[match(city, metadata$mcc_code)]]
  s1[!is.na(URAU_CODE)]
}


#' Inspect the published second-stage model
#'
#' `additional_data.zip` ships `meta-model.RData`, the fitted second-stage
#' model Masselot used to predict the curve of every city from its
#' characteristics. Loading it makes three further checks possible: predicting
#' the Italian curves again from updated city characteristics, leaving Italy
#' out of the meta-model to measure the extrapolation error, and adding our
#' new Italian first-stage estimates to the pool. This function only reports
#' what the file contains, so the objects can be used knowingly.
#'
#' @param paths Output of [download_zenodo()].
#' @return data.table with one row per object (name, class, size), or `NULL`
#'   when the file cannot be loaded.
inspect_metamodel <- function(paths) {
  f <- tryCatch({
    zip <- zpath(paths, "additional_data")
    inside <- utils::unzip(zip, list = TRUE)$Name
    hit <- inside[basename(inside) == "meta-model.RData"]
    if (!length(hit)) return(NULL)
    ex <- tempfile(); utils::unzip(zip, files = hit[1], exdir = ex)
    file.path(ex, hit[1])
  }, error = function(e) NULL)
  if (is.null(f)) return(NULL)
  e <- new.env()
  ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
  if (!ok) return(NULL)
  data.table::rbindlist(lapply(ls(e), function(n) {
    x <- get(n, envir = e)
    data.table::data.table(object = n, class = paste(class(x), collapse = "/"),
                           length = length(x),
                           elements = paste(utils::head(names(x), 8), collapse = ", "))
  }))
}
