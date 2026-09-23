#' Create a directory if it does not exist
#'
#' @param path Directory path.
#' @return `path`, invisibly.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' Download a file unless it is already on disk
#'
#' @param url Source URL.
#' @param dest Destination path.
#' @return `dest`.
download_if_missing <- function(url, dest) {
  ensure_dir(dirname(dest))
  if (!file.exists(dest)) {
    old <- options(timeout = max(3600, getOption("timeout")))
    on.exit(options(old))
    tmp <- paste0(dest, ".part")
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, dest)
  }
  dest
}

#' Read JSON from a URL, retrying on errors and rate limits
#'
#' Open-Meteo answers HTTP 429 when the daily quota is exhausted. Waits grow
#' linearly; after `tries` failures the error is raised, and since all
#' downloads are cached, re-running the pipeline resumes where it stopped.
#'
#' @param url Request URL.
#' @param tries Number of attempts.
#' @param wait Base waiting time in seconds.
#' @return The parsed JSON (list).
get_json <- function(url, tries = 5, wait = 15) {
  last <- NULL
  for (i in seq_len(tries)) {
    res <- tryCatch(jsonlite::fromJSON(url), error = function(e) e)
    if (!inherits(res, "error")) {
      if (isTRUE(res$error)) stop("API error: ", res$reason, "\n", url)
      return(res)
    }
    last <- conditionMessage(res)
    Sys.sleep(wait * i)
  }
  stop("Request failed after ", tries, " attempts: ", last, "\n", url)
}

#' Stem of an Urban Audit city code
#'
#' Masselot's metadata and GISCO releases may differ in the trailing version
#' digit (e.g. `IT001C` vs `IT001C1`); the first six characters identify the
#' city.
#'
#' @param code Character vector of Urban Audit codes.
#' @return Character vector of 6-character stems.
urau_stem <- function(code) substr(toupper(trimws(code)), 1, 6)

#' Pad Istat comune codes to six digits
#'
#' @param x Codes as character or numeric (e.g. `"15146"`, `15146`, `"IT_015146"`).
#' @return Character vector like `"015146"`.
pad_procom <- function(x) {
  x <- gsub("[^0-9]", "", as.character(x))
  sprintf("%06d", as.integer(x))
}

#' Normalise place names for matching
#'
#' @param x Character vector.
#' @return Lower-case ASCII names without punctuation.
norm_name <- function(x) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^A-Za-z]", "", x))
  x
}

#' Value of a percentile from a Masselot temperature-distribution row
#'
#' `tmean_distribution.csv` stores percentiles as columns named like `"95.0%"`.
#' Values between stored percentiles are linearly interpolated.
#'
#' @param row One-row data frame/data.table from `tmean_distribution.csv`.
#' @param p Percentile(s), 0-100.
#' @return Numeric temperature(s).
pct_value <- function(row, p) {
  pc <- grep("%$", names(row), value = TRUE)
  pn <- as.numeric(sub("%", "", pc))
  v  <- as.numeric(unlist(row[1, pc, with = FALSE]))
  o  <- order(pn)
  stats::approx(pn[o], v[o], xout = p, rule = 2)$y
}

#' Multivariate normal draws
#'
#' @param n Number of draws.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @return `n x length(mu)` matrix.
rmvn <- function(n, mu, Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  L <- tryCatch(chol(Sigma), error = function(e) {
    ev <- eigen(Sigma, symmetric = TRUE)
    ev$values[ev$values < 1e-12] <- 1e-12
    chol(ev$vectors %*% diag(ev$values) %*% t(ev$vectors))
  })
  z <- matrix(stats::rnorm(n * length(mu)), n)
  sweep(z %*% L, 2, mu, "+")
}

#' Periods used in the analysis
#'
#' @param cfg Configuration list.
#' @param data_end Last date with mortality data.
#' @return Named list of periods, each a list with `start`, `end`, `label`.
make_periods <- function(cfg, data_end) {
  p <- list(test = list(start = cfg$test_start, end = data_end,
                        label = sprintf("%s to %s", format(cfg$test_start),
                                        format(data_end))))
  if (isTRUE(cfg$run_reproduction))
    p$repro <- list(start = cfg$repro_start, end = cfg$repro_end,
                    label = sprintf("%s to %s", format(cfg$repro_start),
                                    format(cfg$repro_end)))
  p
}
