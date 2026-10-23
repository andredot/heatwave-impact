#' Timing check: do deaths follow temperature, or precede it?
#'
#' A misaligned date anywhere between the Istat file and the temperature
#' series would shift observed deaths against temperatures by whole days. This
#' measures the shift directly from the data: within each city and summer, the
#' daily excess deaths (observed minus the baseline B) and the daily
#' temperature are centred on their summer means, and their correlation is
#' computed with temperature moved by k days, pooled over cities and summers
#' (weighted by days). With correct dates the peak sits at k = 0 to 2 (deaths
#' on the day of the heat and the following days); a peak at negative k would
#' mean deaths come *before* the temperature that is supposed to cause them,
#' i.e. a date error. The same is reported per heatwave episode: day of the
#' highest excess minus day of the highest temperature.
#'
#' @param cf Output of [fit_counterfactual()] (any period).
#' @param episodes Output of [detect_episodes()] for the same period, or `NULL`.
#' @param max_lag Largest shift, in days, either way.
#' @param months Months kept (the warm season).
#' @return List with `ccf` (`lag`, `r`, `n`; positive lag = deaths after
#'   temperature), `peak_lag`, and `episodes` (one row per episode with
#'   `offset` = day of peak excess - day of peak temperature), or `NULL`.
lag_check <- function(cf, episodes = NULL, max_lag = 10, months = 5:9) {
  if (is.null(cf) || !nrow(cf$daily)) return(NULL)
  d <- cf$daily[, .(e = sum(deaths) - sum(B), t = tmean[1]), by = .(URAU_CODE, date)]
  data.table::setorder(d, URAU_CODE, date)
  d[, yr := data.table::year(date)]
  # temperature k days earlier, within the same city's continuous series
  out <- data.table::rbindlist(lapply(-max_lag:max_lag, function(k) {
    x <- data.table::copy(d)
    x[, tk := data.table::shift(t, abs(k), type = if (k >= 0) "lag" else "lead"),
      by = URAU_CODE]
    x <- x[data.table::month(date) %in% months & !is.na(tk)]
    x[, `:=`(ec = e - mean(e), tc = tk - mean(tk)), by = .(URAU_CODE, yr)]
    # pooled correlation: covariance and variances summed over city-summers
    data.table::data.table(lag = k, r = sum(x$ec * x$tc) /
                             sqrt(sum(x$ec^2) * sum(x$tc^2)), n = nrow(x))
  }))
  ep <- NULL
  if (!is.null(episodes) && nrow(episodes)) {
    ep <- data.table::rbindlist(lapply(seq_len(nrow(episodes)), function(i) {
      w <- d[URAU_CODE == episodes$URAU_CODE[i] & date >= episodes$start[i] - 3 &
               date <= episodes$win_end[i]]
      if (nrow(w) < 5) return(NULL)
      # 3-day means, so a single noisy day does not set the peak
      w[, `:=`(e3 = data.table::frollmean(e, 3, align = "center"),
               t3 = data.table::frollmean(t, 3, align = "center"))]
      w <- w[!is.na(e3)]
      data.table::data.table(URAU_CODE = episodes$URAU_CODE[i], start = episodes$start[i],
                             end = episodes$end[i], ndays = episodes$ndays[i],
                             peak_t = w$date[which.max(w$t3)], peak_e = w$date[which.max(w$e3)])
    }))
    if (nrow(ep)) ep[, offset := as.integer(peak_e - peak_t)]
  }
  list(ccf = out, peak_lag = out$lag[which.max(out$r)], episodes = ep)
}
