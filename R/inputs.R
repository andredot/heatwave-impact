# =============================================================================
#  Inputs that used to be placed by hand: the Istat mortality release and the
#  official FluMOMO code. Both are downloaded and cached, so a fresh clone
#  runs end to end.
# =============================================================================

#' Download and unzip the Istat daily deaths release
#'
#' The URL carries the release month and a version suffix
#' (".../giugno-2026/decessi-comunali-mese-provvisori-3.zip"), so it is taken
#' from `cfg$istat_url` rather than guessed. A file already present at
#' `cfg$istat_csv` wins: put a release there by hand and nothing is downloaded.
#' The CSV inside the archive is found by extension, since its name changes
#' between releases.
#'
#' @param cfg Configuration list (uses `istat_url`, `istat_csv`, `istat_dir`).
#' @return Path to the CSV.
download_istat <- function(cfg) {
  if (!is.null(cfg$istat_csv) && file.exists(cfg$istat_csv)) return(cfg$istat_csv)
  if (is.null(cfg$istat_url) || !nzchar(cfg$istat_url))
    stop("No Istat file at ", cfg$istat_csv, " and no istat_url in config.R.")
  dir <- ensure_dir(cfg$istat_dir %||% "data/raw")
  zip <- file.path(dir, basename(sub("\\?.*$", "", cfg$istat_url)))
  if (!file.exists(zip)) {
    old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
    utils::download.file(cfg$istat_url, zip, mode = "wb", quiet = TRUE)
  }
  inside <- utils::unzip(zip, list = TRUE)$Name
  csv <- inside[grepl("\\.csv$", inside, ignore.case = TRUE)]
  if (!length(csv)) stop("No CSV inside ", zip)
  if (length(csv) > 1) {
    sizes <- utils::unzip(zip, list = TRUE)$Length[match(csv, inside)]
    csv <- csv[which.max(sizes)]          # the comune-level file is the big one
  }
  utils::unzip(zip, files = csv, exdir = dir, junkpaths = TRUE)
  out <- file.path(dir, basename(csv))
  message("Istat release: ", basename(zip), " -> ", out)
  out
}

#' Download the official FluMOMO code
#'
#' Unzips the EuroMOMO archive into `cfg$flumomo$code_dir` when the estimation
#' script is not already there.
#'
#' @param cfg Configuration list (uses `flumomo$code_dir`, `flumomo$code_url`).
#' @return Path to `Estimation_v42.R`.
download_flumomo <- function(cfg) {
  dir <- cfg$flumomo$code_dir
  est <- file.path(dir, "Estimation_v42.R")
  if (file.exists(est)) return(est)
  ensure_dir(dir)
  url <- cfg$flumomo$code_url %||% "https://euromomo.eu/uploads/data/FluMOMO_version_4_2_R.zip"
  zip <- file.path(dir, basename(url))
  utils::download.file(url, zip, mode = "wb", quiet = TRUE)
  utils::unzip(zip, exdir = dir, junkpaths = TRUE)
  if (!file.exists(est))
    stop("Estimation_v42.R not found after unzipping ", url,
         ": check the archive layout in ", dir)
  message("FluMOMO code downloaded into ", dir)
  est
}
