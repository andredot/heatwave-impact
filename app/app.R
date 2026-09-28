# =============================================================================
#  Decessi attribuibili al caldo - citta' italiane
#  ---------------------------------------------------------------------------
#  Fonti della temperatura:
#    - Previsione: Open-Meteo, dal giorno scelto in avanti. Con la "macchina
#      del tempo" si puo' scegliere una data passata: in quel caso vengono
#      usate le previsioni realmente emesse allora (archivio delle emissioni
#      precedenti, fino a 7 giorni di anticipo).
#    - Ricostruzione ERA5-Land: temperature osservate fra due date.
#
#  Curve esposizione-risposta: quelle pubblicate (Masselot et al. 2023) e
#  quelle ricalibrate sui dati italiani recenti. Le stime hanno un intervallo
#  al 95% che combina l'incertezza della curva e quella del livello atteso.
#
#  Scheda EuroMOMO: decessi osservati con due attese indipendenti, la nostra e
#  quella del codice ufficiale FluMOMO, per la regione o per le sue citta'.
#
#  Dati precalcolati in app_data.rds (targets::tar_make()).
# =============================================================================
library(shiny)
library(data.table)
library(ggplot2)
library(sf)
# jsonlite is not attached: its validate() would mask shiny's
fromJSON <- jsonlite::fromJSON

`%||%` <- function(x, y) if (is.null(x)) y else x

BUNDLE <- readRDS("app_data.rds")
CITIES <- as.data.table(BUNDLE$cities)
DAILY  <- as.data.table(BUNDLE$daily)
# observed deaths over the whole Istat history, when the bundle carries them
DEATHS <- if (!is.null(BUNDLE$deaths_long)) as.data.table(BUNDLE$deaths_long) else
  DAILY[, .(deaths = sum(deaths)), by = .(URAU_CODE, date)]
BASELINE_FROM <- min(DAILY$date)
# older bundles lack the fields added for the intervals and the EuroMOMO tab
if (!"B_logsd" %in% names(DAILY)) {
  DAILY[, B_logsd := 0]
  message("app_data.rds predates the baseline uncertainty: intervals will show ",
          "only the uncertainty of the curves. Run targets::tar_make() to refresh it.")
}
for (col in c("tmean", "logrr", "mmt")) if (!col %in% names(DAILY)) DAILY[, (col) := NA_real_]
HAS_VCOV <- !is.null(BUNDLE$curves[[1]]$published[[1]]$vcov)
HAS_BSD  <- any(DAILY$B_logsd > 0, na.rm = TRUE)
BUNDLE_INFO <- sprintf(
  "app_data.rds: %s, creato il %s | covarianze curve: %s | incertezza baseline: %s | FluMOMO: %s",
  normalizePath("app_data.rds"), format(BUNDLE$built %||% NA),
  if (HAS_VCOV) "s\u00ec" else "NO", if (HAS_BSD) "s\u00ec" else "NO",
  if (!is.null(BUNDLE$flumomo)) "s\u00ec" else "NO")
message(BUNDLE_INFO)
if (!HAS_VCOV || !HAS_BSD)
  message("Senza questi campi gli intervalli non possono essere calcolati: ",
          "rieseguire targets::tar_make() con la versione aggiornata di R/app_export.R ",
          "e copiare il file prodotto accanto ad app.R.")
REGIONI <- sort(unique(CITIES$region))
NSIM <- 400
FC_MAX_LEAD <- 7      # lead massimo nell'archivio delle emissioni precedenti
PUB <- "Pubblicate (Masselot 2023)"
COL <- c("#B2182B", "#1B7837", "#762A83", "#B35806")

it_label <- function(x) sub("^Italian curves", "Curve italiane", x)
pct_fmt <- function(p) if (is.na(p)) "n.d." else sprintf("%.0f%%", 100 * p)

# ---- utilita' numeriche -----------------------------------------------------
basis <- function(t, spec) {
  suppressWarnings(unclass(splines::bs(t, knots = spec$knots, degree = spec$degree,
                                       Boundary.knots = spec$bound)))
}

#' Log relative risk of a curve at given temperatures, centred on its MMT
#'
#' @param t Temperatures.
#' @param spec Basis specification of the city.
#' @param curve List with `beta` and `mmt`.
#' @return Numeric vector.
rr_of <- function(t, spec, curve) {
  as.numeric(basis(t, spec) %*% curve$beta) -
    as.numeric(basis(curve$mmt, spec) %*% curve$beta)
}

#' Draws of the log relative risk of a curve, centred on its MMT
rr_draws <- function(t, spec, curve, nsim) {
  B <- basis(t, spec); Bm <- basis(curve$mmt, spec)
  V <- curve$vcov
  bd <- if (is.null(V)) matrix(curve$beta, nsim, length(curve$beta), byrow = TRUE) else {
    V <- (V + t(V)) / 2
    L <- tryCatch(chol(V), error = function(e) {
      ev <- eigen(V, symmetric = TRUE); ev$values[ev$values < 1e-12] <- 1e-12
      chol(ev$vectors %*% diag(ev$values) %*% t(ev$vectors))
    })
    sweep(matrix(rnorm(nsim * length(curve$beta)), nsim) %*% L, 2, curve$beta, "+")
  }
  bd %*% t(B) - as.numeric(bd %*% t(Bm))
}

get_json <- function(url, tries = 3, wait = 4) {
  for (i in seq_len(tries)) {
    r <- tryCatch(jsonlite::fromJSON(url), error = function(e) e)
    if (!inherits(r, "error")) {
      if (isTRUE(r$error)) stop(r$reason)
      return(r)
    }
    Sys.sleep(wait)
  }
  stop("Richiesta non riuscita: ", url)
}

# ---- temperature ------------------------------------------------------------
#' Observed ERA5-Land daily means
fetch_era5 <- function(lat, lon, from, to, tz) {
  js <- get_json(sprintf(paste0("https://archive-api.open-meteo.com/v1/archive?",
                                "latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s",
                                "&daily=temperature_2m_mean&models=era5_land&timezone=%s"),
                         lat, lon, format(from), format(to), tz))
  data.table(date = as.Date(js$daily$time),
             tmean = as.numeric(js$daily$temperature_2m_mean), lead = NA_integer_)
}

#' Forecast issued on `as_of`
#'
#' For today, the live forecast; for a past date, the runs actually issued then,
#' taken from the previous-runs archive, where lead = date - as_of.
fetch_forecast <- function(lat, lon, as_of, to, tz, model) {
  if (as_of >= Sys.Date()) {
    js <- get_json(sprintf(paste0("https://api.open-meteo.com/v1/forecast?latitude=%.4f",
                                  "&longitude=%.4f&daily=temperature_2m_mean",
                                  "&forecast_days=16&past_days=7&timezone=%s%s"),
                           lat, lon, tz,
                           if (nzchar(model %||% "")) paste0("&models=", model) else ""))
    d <- data.table(date = as.Date(js$daily$time),
                    tmean = as.numeric(js$daily$temperature_2m_mean))
    d[, lead := as.integer(date - as_of)]
    return(d[date <= to])
  }
  leads <- 1:FC_MAX_LEAD
  vars <- c(sprintf("temperature_2m_previous_day%d", leads))
  js <- get_json(sprintf(paste0("https://previous-runs-api.open-meteo.com/v1/forecast?",
                                "latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s",
                                "&hourly=%s&models=%s&timezone=%s"),
                         lat, lon, format(as_of + 1), format(to),
                         paste(vars, collapse = ","), model, tz))
  h <- as.data.table(js$hourly)
  h[, date := as.Date(substr(time, 1, 10))]
  long <- melt(h, id.vars = "date", measure.vars = setdiff(names(h), c("time", "date")),
               variable.name = "var", value.name = "t")
  long[, lead := as.integer(sub(".*previous_day", "", var))]
  d <- long[, .(tmean = if (sum(!is.na(t)) == 24) mean(t) else NA_real_), by = .(date, lead)]
  d <- d[!is.na(tmean)]
  d[, want := as.integer(date - as_of)]
  d[lead == want, .(date, tmean, lead)]
}

#' Temperatures for one city, aligned to the scale the curves were built on
city_temperature <- function(code, source, as_of, from, to, debias) {
  ci <- CITIES[URAU_CODE == code]
  tt <- if (source == "era5") fetch_era5(ci$lat, ci$lon, from, to, BUNDLE$era5_timezone)
        else fetch_forecast(ci$lat, ci$lon, as_of, to, BUNDLE$era5_timezone,
                            BUNDLE$forecast_model)
  if (!nrow(tt)) return(tt)
  tt[, month := month(date)]
  if (source == "fc" && isTRUE(debias) && !is.null(BUNDLE$fc_bias)) {
    fb <- as.data.table(BUNDLE$fc_bias)[URAU_CODE == code]
    if (nrow(fb)) {
      tt[, lead_c := pmin(pmax(lead, min(fb$lead)), max(fb$lead))]
      tt <- merge(tt, fb[, .(month, lead_c = lead, bias)], by = c("month", "lead_c"),
                  all.x = TRUE)
      tt[is.na(bias), bias := 0][, tmean := tmean - bias]
    }
  }
  if (isTRUE(BUNDLE$bias_correct)) {
    dl <- as.data.table(BUNDLE$delta)[URAU_CODE == code]
    tt <- merge(tt, dl[, .(month, delta)], by = "month", all.x = TRUE)
    tt[is.na(delta), delta := 0][, tmean := tmean - delta]
  }
  setorder(tt, date)
  tt[, .(date, tmean, lead)]
}

#' Expected deaths without heat, by day and age group
#'
#' Istat-based where the data reach, then the seasonal profile of the last
#' twelve months, which is what a forecast has to rely on.
city_baseline <- function(code, dates) {
  bs <- DAILY[URAU_CODE == code]
  if (!nrow(bs)) return(NULL)
  bs[, doy := as.integer(format(date, "%j"))]
  clim <- bs[date > max(date) - 365, .(Bc = mean(B), sdc = mean(B_logsd)), by = .(agegroup, doy)]
  grid <- CJ(agegroup = unique(bs$agegroup), date = dates)
  grid[, doy := as.integer(format(date, "%j"))]
  x <- merge(grid, bs[, .(agegroup, date, B, B_logsd, obs = deaths)],
             by = c("agegroup", "date"), all.x = TRUE)
  x <- merge(x, clim, by = c("agegroup", "doy"), all.x = TRUE)
  x[is.na(B), `:=`(B = Bc, B_logsd = sdc)]
  x[is.na(B_logsd), B_logsd := 0]
  x[!is.na(B), .(agegroup, date, B, B_logsd, obs)]
}

#' Heat deaths per day for one city and one set of curves, with draws
city_attributable <- function(code, temps, curve_name, nsim) {
  cu <- BUNDLE$curves[[code]]
  base <- city_baseline(code, temps$date)
  if (is.null(base) || !nrow(base)) return(NULL)
  ages <- unique(base$agegroup)
  point <- rep(0, nrow(temps)); draws <- matrix(0, nsim, nrow(temps))
  for (g in ages) {
    b <- base[agegroup == g][match(temps$date, date)]
    cv <- if (identical(curve_name, PUB)) cu$published[[g]] else cu$recalibrated[[curve_name]]
    if (is.null(cv)) next
    lr <- rr_draws(temps$tmean, cu$spec, cv, nsim)
    hot <- temps$tmean >= cv$mmt
    bmult <- exp(matrix(rnorm(nsim, 0, mean(b$B_logsd, na.rm = TRUE)), nsim, nrow(temps)))
    pt <- b$B * pmax(exp(as.numeric(basis(temps$tmean, cu$spec) %*% cv$beta) -
                           as.numeric(basis(cv$mmt, cu$spec) %*% cv$beta)) - 1, 0) * hot
    point <- point + ifelse(is.na(pt), 0, pt)
    dr <- sweep(pmax(exp(lr) - 1, 0), 2, b$B * hot, "*") * bmult
    dr[is.na(dr)] <- 0
    draws <- draws + dr
  }
  list(point = point, draws = draws,
       baseline = base[, .(B = sum(B)), by = date][match(temps$date, date), B],
       observed = base[, .(obs = sum(obs)), by = date][match(temps$date, date), obs])
}

# ---- scenari ----------------------------------------------------------------
#' Probability that a heatwave starts within the horizon
#'
#' A day counts as hot when the temperature reaches the city's heatwave
#' threshold (the percentile behind the episode definition). Forecast error is
#' taken from the archive: for each lead, the spread of past errors turns the
#' forecast into a probability. An episode is `hw_min` consecutive hot days;
#' the probability of at least one such run starting in the window is computed
#' assuming days are independent given the forecast, which understates
#' persistence and is therefore conservative.
#'
#' @param temps data.table `date`, `tmean`, `lead` for one city.
#' @param code Urban Audit code.
#' @param hw_min Consecutive days required.
#' @return List with `p_day` (per day) and `p_onset`.
onset_probability <- function(temps, code, hw_min = 3) {
  thr <- BUNDLE$curves[[code]]$spec$hw_threshold
  if (is.null(thr) || !nrow(temps)) return(list(p_day = numeric(0), p_onset = NA_real_))
  sdv <- rep(1.2, nrow(temps))
  fb <- if (!is.null(BUNDLE$fc_bias)) as.data.table(BUNDLE$fc_bias)[URAU_CODE == code] else NULL
  if (!is.null(fb) && nrow(fb) && "sd" %in% names(fb)) {
    m <- month(temps$date)
    lead_c <- pmin(pmax(temps$lead, min(fb$lead)), max(fb$lead))
    key <- data.table(month = m, lead = lead_c)
    got <- merge(key, fb[, .(month, lead, sd)], by = c("month", "lead"), all.x = TRUE,
                 sort = FALSE)$sd
    sdv <- ifelse(is.na(got), 1.2, got)
  }
  p <- stats::pnorm(temps$tmean, mean = thr, sd = pmax(sdv, .3))   # P(T >= threshold)
  runs <- if (length(p) >= hw_min)
    vapply(seq_len(length(p) - hw_min + 1), function(i) prod(p[i:(i + hw_min - 1)]), 0)
  else 0
  list(p_day = p, p_onset = 1 - prod(1 - runs))
}

#' Scenarios: the forecast window, then the episode continuing
#'
#' The forecast window is costed as it stands. Each further block of days
#' assumes the heat continues at the level of the hottest forecast days, so the
#' extra deaths are conditional on the episode lasting that long, not a
#' prediction that it will.
#'
#' @param codes Cities.
#' @param temps_by_city Named list of temperature tables.
#' @param curve_name Curve set to use.
#' @param blocks Lengths of the additional blocks, in days.
#' @param nsim Draws.
#' @return data.table with one row per scenario.
scenario_table <- function(codes, temps_by_city, curve_name, blocks = c(5, 5), nsim = NSIM) {
  acc <- NULL; add <- function(x, y) if (is.null(x)) y else x + y
  point <- 0; draws <- 0
  peak_point <- 0; peak_draws <- 0
  for (cd in codes) {
    tt <- temps_by_city[[cd]]
    if (is.null(tt) || !nrow(tt)) next
    a <- city_attributable(cd, tt, curve_name, nsim)
    if (is.null(a)) next
    point <- point + sum(a$point); draws <- add(draws, rowSums(a$draws))
    # a further day at the level of the three hottest forecast days
    hot <- tt[order(-tmean)][seq_len(min(3, .N))]
    one <- data.table(date = max(tt$date) + 1, tmean = mean(hot$tmean), lead = max(tt$lead))
    b <- city_attributable(cd, one, curve_name, nsim)
    if (!is.null(b)) {
      peak_point <- peak_point + sum(b$point); peak_draws <- add(peak_draws, rowSums(b$draws))
    }
  }
  q <- function(x) stats::quantile(x, c(.025, .975))
  rows <- list(data.table(scenario = "Previsione attuale", giorni = nrow(temps_by_city[[codes[1]]]),
                          stima = point, lo = q(draws)[1], hi = q(draws)[2]))
  cum_p <- point; cum_d <- draws; extra <- 0
  for (b in blocks) {
    extra <- extra + b
    cum_p <- cum_p + b * peak_point; cum_d <- cum_d + b * peak_draws
    rows[[length(rows) + 1]] <- data.table(
      scenario = sprintf("Se prosegue altri %d giorni", extra),
      giorni = nrow(temps_by_city[[codes[1]]]) + extra,
      stima = cum_p, lo = q(cum_d)[1], hi = q(cum_d)[2])
  }
  rbindlist(rows)
}

#' Cumulative paths: reconstruction, forecast, and duration scenarios
#'
#' Builds one trajectory of cumulative heat deaths per scenario, in the style
#' of the IPCC pathway figures: a single observed history up to the issue
#' date, then branches. Each branch keeps its own 95% band, and the bands
#' widen with time because the draws accumulate.
#'
#' @param codes Cities.
#' @param as_of Issue date (the vertical line in the chart).
#' @param horizon Days of forecast.
#' @param back Days of reconstruction shown before the issue date.
#' @param curve_name Curve set.
#' @param debias Correct the forecast bias.
#' @param blocks Additional blocks of days, each a scenario.
#' @param nsim Draws.
#' @return List with `paths` (date, scenario, cum, lo, hi), `as_of`, `temps`.
scenario_paths <- function(codes, as_of, horizon, back, curve_name, debias,
                           blocks = c(5, 5), nsim = NSIM) {
  hist_to <- min(as_of, Sys.Date() - 6)
  hist_from <- hist_to - back
  acc_hist <- NULL; acc_fc <- NULL; peak <- NULL; temps <- list()
  add <- function(a, b) if (is.null(a)) b else
    list(dates = a$dates, draws = a$draws + b$draws)
  for (cd in codes) {
    th <- try(city_temperature(cd, "era5", as_of, hist_from, hist_to, FALSE), silent = TRUE)
    # the forecast call also returns the days between the end of ERA5-Land and
    # the issue date, so the reconstruction runs up to the vertical line
    tf <- try(city_temperature(cd, "fc", as_of, hist_to + 1, as_of + horizon, debias),
              silent = TRUE)
    if (inherits(tf, "try-error") || !nrow(tf)) next
    gap <- tf[date <= as_of]; tf <- tf[date > as_of]
    if (!inherits(th, "try-error") && nrow(th))
      th <- rbind(th, gap, fill = TRUE)[order(date)] else th <- gap
    if (nrow(th)) {
      a <- city_attributable(cd, th, curve_name, nsim)
      if (!is.null(a)) acc_hist <- add(acc_hist, list(dates = th$date, draws = a$draws))
    }
    if (!nrow(tf)) next
    temps[[cd]] <- tf
    b <- city_attributable(cd, tf, curve_name, nsim)
    if (is.null(b)) next
    acc_fc <- add(acc_fc, list(dates = tf$date, draws = b$draws))
    hot <- tf[order(-tmean)][seq_len(min(3, .N))]
    one <- data.table(date = max(tf$date) + 1, tmean = mean(hot$tmean), lead = max(tf$lead))
    c1 <- city_attributable(cd, one, curve_name, nsim)
    if (!is.null(c1)) peak <- if (is.null(peak)) c1$draws[, 1] else peak + c1$draws[, 1]
  }
  if (is.null(acc_fc)) return(NULL)
  q <- function(m) data.table(cum = colMeans(m),
                              lo = apply(m, 2, stats::quantile, .025),
                              hi = apply(m, 2, stats::quantile, .975))
  cum_hist <- if (!is.null(acc_hist)) t(apply(acc_hist$draws, 1, cumsum)) else NULL
  base <- if (is.null(cum_hist)) 0 else cum_hist[, ncol(cum_hist)]
  paths <- list()
  if (!is.null(cum_hist))
    paths[[1]] <- cbind(data.table(date = acc_hist$dates, scenario = "Ricostruzione"),
                        q(cum_hist))
  cum_fc <- sweep(t(apply(acc_fc$draws, 1, cumsum)), 1, base, "+")
  join <- function(dt, date0, m0) if (is.null(m0)) dt else
    rbind(cbind(data.table(date = date0, scenario = dt$scenario[1]), q(m0)), dt)
  last_hist <- if (!is.null(cum_hist)) max(acc_hist$dates) else NULL
  paths[[length(paths) + 1]] <- join(
    cbind(data.table(date = acc_fc$dates, scenario = "Previsione"), q(cum_fc)),
    last_hist, if (!is.null(cum_hist)) cum_hist[, ncol(cum_hist), drop = FALSE] else NULL)
  end <- cum_fc[, ncol(cum_fc)]; last <- max(acc_fc$dates); extra <- 0
  if (!is.null(peak)) for (b in blocks) {
    days <- seq_len(b)
    m <- sapply(days, function(k) end + k * peak)
    if (is.null(dim(m))) m <- matrix(m, nrow = length(end))
    extra <- extra + b
    paths[[length(paths) + 1]] <- join(
      cbind(data.table(date = last + days,
                       scenario = sprintf("Se prosegue (+%d giorni)", extra)), q(m)),
      last, matrix(end, ncol = 1))
    end <- m[, ncol(m)]; last <- last + b
  }
  out <- rbindlist(paths)
  out[, date := as.Date(date, origin = "1970-01-01")]   # rbind can drop the class
  list(paths = out, as_of = as_of, temps = temps,
       hist_from = hist_from, hist_to = hist_to)
}

# ---- the figure, in English, rebuilt on the fly -----------------------------

#' Daily deaths: observed, expected, and what the curves predict
#'
#' Observed deaths and the counterfactual are drawn as 7-day moving averages
#' with the difference shaded (red above the expectation, blue below), and one
#' line per exposure-response function shows the deaths those curves predict
#' (expected plus heat). Temperatures come from ERA5-Land up to the issue date
#' and from the Open-Meteo forecast after it, and a vertical line marks the
#' boundary.
#'
#' @param d data.table: `date`, `observed`, `expected`, `obs7`, `exp7`.
#' @param pred data.table: `date`, `curve`, `pred7` (may be empty).
#' @param as_of Issue date, or `NULL` for a purely historical chart.
#' @param title,subtitle_extra,caption,ylab Text.
#' @return A ggplot.
plot_deaths <- function(d, pred, as_of, title, caption, exp_label,
                        ylab = "Deaths per day", subtitle_extra = NULL) {
  total <- round(sum(d$observed - d$expected, na.rm = TRUE))
  span <- sprintf("%s to %s", format(min(d$date), "%d %b"), format(max(d$date), "%d %b %Y"))
  sub <- sprintf(paste("7-day moving averages. Red: days above the expectation,",
                       "blue: days below (balance %s, %s)"),
                 format(total, big.mark = ","), span)
  if (!is.null(subtitle_extra)) sub <- paste0(sub, ". ", subtitle_extra)
  cols <- c("Observed" = "#1F2933"); cols[exp_label] <- "#2C7FB8"
  if (nrow(pred)) {
    cu <- unique(pred$curve)
    cols[cu] <- c("#B2182B", "#1B7837", "#762A83", "#B35806")[seq_along(cu)]
  }
  p <- ggplot(d, aes(date)) +
    geom_ribbon(aes(ymin = exp7, ymax = pmax(obs7, exp7)), fill = "#C0392B", alpha = .28) +
    geom_ribbon(aes(ymin = pmin(obs7, exp7), ymax = exp7), fill = "#3A8DAE", alpha = .18)
  if (!is.null(as_of) && as_of >= min(d$date) && as_of <= max(d$date))
    p <- p + geom_vline(xintercept = as.numeric(as_of), colour = "grey45", linetype = "22")
  if (nrow(pred))
    p <- p + geom_line(data = pred, aes(date, pred7, colour = curve), linewidth = .8)
  p + geom_line(aes(y = exp7, colour = exp_label), linewidth = .8, linetype = "22") +
    geom_line(aes(y = obs7, colour = "Observed"), linewidth = 0) +
    scale_colour_manual(values = cols, breaks = names(cols),
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_date(date_labels = "%d %b", expand = expansion(mult = c(.01, .01))) +
    scale_y_continuous(expand = expansion(mult = c(.05, .12))) +
    labs(title = title, subtitle = wrap_txt(sub, 110), x = NULL, y = ylab,
         colour = NULL, caption = wrap_txt(caption, 120)) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = rel(1.2)),
          plot.subtitle = element_text(colour = "grey30", margin = margin(b = 12)),
          plot.caption = element_text(colour = "grey45", size = rel(.72), hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          legend.position = "top", legend.justification = "left",
          legend.key.width = unit(28, "pt"), panel.grid.minor = element_blank(),
          axis.title.y = element_text(colour = "grey30", margin = margin(r = 8)))
}

wrap_txt <- function(x, n = 120) paste(strwrap(x, width = n), collapse = "\n")

#' Temperatures of one city over a period spanning past and future
#'
#' ERA5-Land up to the issue date, the forecast issued that day afterwards,
#' and the days between the end of ERA5-Land and the issue date taken from the
#' same forecast run.
#'
#' @param code City.
#' @param from,to Period.
#' @param as_of Issue date.
#' @param debias Correct the forecast bias.
#' @return data.table `date`, `tmean`, `lead`, `source`.
city_temperature_span <- function(code, from, to, as_of, debias) {
  era_to <- min(to, as_of, Sys.Date() - 6)
  out <- list()
  if (from <= era_to) {
    h <- try(city_temperature(code, "era5", as_of, from, era_to, FALSE), silent = TRUE)
    if (!inherits(h, "try-error") && nrow(h)) out[[1]] <- h[, source := "ERA5-Land"]
  }
  if (to > era_to) {
    f <- try(city_temperature(code, "fc", as_of, era_to + 1, to, debias), silent = TRUE)
    if (!inherits(f, "try-error") && nrow(f))
      out[[length(out) + 1]] <- f[, source := fifelse(date <= as_of, "ERA5-Land", "Forecast")]
  }
  if (!length(out)) return(NULL)
  rbindlist(out, fill = TRUE)[order(date)][!duplicated(date)]
}

#' Observed deaths, counterfactual and predicted deaths over a period
#'
#' @param codes Cities (for the FluMOMO region baseline, the region's cities).
#' @param from,to Period.
#' @param as_of Issue date.
#' @param curve_names Curve sets to predict with.
#' @param debias Correct the forecast bias.
#' @param baseline `"project"` (our counterfactual) or `"flumomo"`.
#' @param scope `"cities"` or `"region"` (only with the FluMOMO baseline).
#' @param nsim Draws.
#' @return List with `daily`, `pred`, `temps`.
combined_series <- function(codes, from, to, as_of, curve_names, debias,
                            baseline = "project", scope = "cities", nsim = NSIM) {
  dates <- seq(from, to, by = "day")
  temps <- list(); base <- NULL; pred <- list(); attr_tot <- list(); att_city <- list()
  for (cd in codes) {
    tt <- city_temperature_span(cd, from, to, as_of, debias)
    if (is.null(tt) || !nrow(tt)) next
    temps[[cd]] <- tt[, code := cd][]
    b <- city_baseline(cd, tt$date)
    if (is.null(b) || !nrow(b)) next
    cu <- BUNDLE$curves[[cd]]
    # expected without heat, cold effect kept in, summed over age groups
    exp_city <- rep(0, nrow(tt))
    for (g in unique(b$agegroup)) {
      bg <- b[agegroup == g][match(tt$date, date)]
      cv <- cu$published[[g]]
      lr <- rr_of(tt$tmean, cu$spec, cv)
      cold <- tt$tmean < cv$mmt
      exp_city <- exp_city + ifelse(is.na(bg$B), 0, bg$B * ifelse(cold, exp(lr), 1))
    }
    base <- if (is.null(base)) data.table(date = tt$date, expected = exp_city)
            else merge(base, data.table(date = tt$date, e = exp_city),
                       by = "date", all = TRUE)[
                         , .(date, expected = rowSums(cbind(expected, e), na.rm = TRUE))]
    for (nm in curve_names) {
      a <- city_attributable(cd, tt, nm, nsim)
      if (is.null(a)) next
      attr_tot[[nm]] <- if (is.null(attr_tot[[nm]]))
        data.table(date = tt$date, att = a$point)
      else merge(attr_tot[[nm]], data.table(date = tt$date, a2 = a$point), by = "date",
                 all = TRUE)[, .(date, att = rowSums(cbind(att, a2), na.rm = TRUE))]
      if (nm == curve_names[1])
        att_city[[cd]] <- data.table(code = cd, att = sum(a$point, na.rm = TRUE))
    }
  }
  if (is.null(base)) return(NULL)
  d <- base[order(date)]
  # observed deaths: the daily series, NA where the Istat data do not reach
  obs <- DEATHS[URAU_CODE %in% codes & date >= from & date <= to,
                .(observed = sum(deaths)), by = date]
  d <- merge(d, obs, by = "date", all.x = TRUE)
  # FluMOMO counterfactual: weekly, spread over days and rescaled to this population
  if (identical(baseline, "flumomo")) {
    fm <- if (scope == "region") BUNDLE$flumomo$region_all else BUNDLE$flumomo$cities
    if (is.null(fm)) return(NULL)
    w <- as.data.table(fm)[order(week_start)]
    if (scope == "region" && !is.null(BUNDLE$region$daily)) {
      rd <- as.data.table(BUNDLE$region$daily)[date >= from & date <= to,
                                               .(date, observed = deaths)]
      d <- merge(d[, .(date)], rd, by = "date", all.x = TRUE)
    } else if (!is.null(BUNDLE$region$cities_daily)) {
      rd <- as.data.table(BUNDLE$region$cities_daily)[date >= from & date <= to,
                                                      .(date, observed = deaths)]
      d <- merge(d[, .(date)], rd, by = "date", all.x = TRUE)
    }
    ref <- w[week_start >= from & week_start <= to]
    sc <- if (nrow(ref) && sum(ref$deaths) > 0 && sum(d$observed, na.rm = TRUE) > 0)
      sum(d$observed, na.rm = TRUE) / sum(ref$deaths) else 1
    d[, expected := sc * stats::splinefun(as.numeric(w$week_start) + 3,
                                          w$expected / 7)(as.numeric(date))]
  }
  d[, `:=`(obs7 = frollmean(observed, 7, align = "center"),
           exp7 = frollmean(expected, 7, align = "center"))]
  for (nm in names(attr_tot)) {
    x <- merge(d[, .(date, expected)], attr_tot[[nm]], by = "date", all.x = TRUE)
    x[is.na(att), att := 0]
    pred[[nm]] <- data.table(date = x$date, curve = it_label(nm),
                             pred7 = frollmean(x$expected + x$att, 7, align = "center"))
  }
  list(daily = d[!is.na(exp7)], pred = rbindlist(pred)[!is.na(pred7)],
       temps = rbindlist(temps, fill = TRUE),
       by_city = rbindlist(att_city))
}

# ================================ UI =========================================# ================================ UI =========================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; }
    .titolo { font-weight: 700; font-size: 1.5rem; margin-bottom: .1rem; }
    .sottotitolo { color: #666; margin-bottom: 1rem; }
    .errore { color: #B00020; font-weight: 600; }
    .nota { color: #666; font-size: .85rem; }"))),
  div(class = "titolo", "Decessi attribuibili al caldo \u2014 citt\u00e0 italiane"),
  div(class = "sottotitolo",
      "Curve di Masselot et al. (2023) e loro ricalibrazione sui dati italiani recenti."),
  tabsetPanel(
    id = "scheda",
    tabPanel(
      "Caldo",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          dateRangeInput("periodo", "Periodo mostrato",
                         start = Sys.Date() - 60, end = Sys.Date() + 7,
                         max = Sys.Date() + 16, format = "dd/mm/yyyy", language = "it",
                         separator = " \u2013 "),
          dateInput("asof", "Data di emissione (macchina del tempo)",
                    value = Sys.Date(), max = Sys.Date(), format = "dd/mm/yyyy",
                    language = "it"),
          div(class = "nota",
              "Prima di questa data le temperature vengono da ERA5-Land, dopo dalle",
              "previsioni emesse allora (al massimo 7 giorni di anticipo se la data",
              "\u00e8 nel passato, 16 se \u00e8 oggi)."),
          checkboxInput("debias", "Correggi la distorsione delle previsioni", TRUE),
          sliderInput("sc_back", "Giorni di ricostruzione negli scenari", 0, 60, 30, 5),
          hr(),
          radioButtons("ambito", "Ambito", c("Regione" = "reg", "Singola citt\u00e0" = "city")),
          conditionalPanel("input.ambito == 'reg'",
                           selectInput("regione", "Regione", choices = REGIONI,
                                       selected = if ("Lombardia" %in% REGIONI) "Lombardia"
                                       else REGIONI[1])),
          conditionalPanel("input.ambito == 'city'",
                           selectInput("citta", "Citt\u00e0", choices = sort(CITIES$name))),
          checkboxGroupInput("curve", "Curve", choices = NULL),
          selectInput("sc_curva", "Curva usata negli scenari", choices = NULL),
          actionButton("vai", "Calcola", class = "btn-primary"),
          br(), br(), uiOutput("msg"), div(class = "nota", textOutput("copertura")),
          br(), div(class = "nota", style = "word-break:break-all;", textOutput("bundle"))
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "uscite",
            tabPanel("Giornalieri", br(), plotOutput("g_giorno", height = "470px")),
            tabPanel("Cumulati", br(), plotOutput("g_cum", height = "430px")),
            tabPanel("Scenari", br(), htmlOutput("sc_testo"), br(),
                     plotOutput("sc_plot", height = "420px"), br(), tableOutput("sc_tab")),
            tabPanel("Temperatura", br(), plotOutput("g_temp", height = "430px")),
            tabPanel("Mappa", br(), plotOutput("g_mappa", height = "520px")),
            tabPanel("Tabella", br(), tableOutput("tab"))
          )
        )
      )
    ),
    tabPanel(
      "EuroMOMO",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          radioButtons("mm_ambito", "Popolazione",
                       c("Intera regione (tutti i comuni)" = "reg",
                         "Citt\u00e0 validate della regione" = "cit")),
          dateRangeInput("mm_periodo", "Periodo",
                         start = as.Date(format(BUNDLE$data_end, "%Y-01-01")),
                         end = BUNDLE$data_end, format = "dd/mm/yyyy", language = "it",
                         separator = " \u2013 "),
          radioButtons("mm_riferimento", "Attesa di riferimento (area colorata)",
                       c("FluMOMO" = "flu", "Modello del progetto" = "our")),
          div(class = "nota",
              "Due attese indipendenti: il modello del progetto (giornaliero, et\u00e0 20+,",
              "senza caldo) e il codice ufficiale FluMOMO (settimanale, senza temperature",
              "estreme, influenza inclusa).")
        ),
        mainPanel(width = 9, plotOutput("g_momo", height = "600px"),
                  br(), tableOutput("tab_momo"))
      )
    )
  )
)

# ============================== SERVER =======================================
server <- function(input, output, session) {

  observe({
    nc <- c(PUB, unique(unlist(lapply(BUNDLE$curves, function(x) names(x$recalibrated)))))
    updateCheckboxGroupInput(session, "curve", choices = setNames(nc, it_label(nc)),
                             selected = nc)
  })

  codici <- reactive({
    if (input$ambito == "reg") CITIES[region == input$regione, URAU_CODE]
    else CITIES[name == input$citta, URAU_CODE]
  })

  risultato <- eventReactive(input$vai, {
    cds <- codici()
    if (!length(cds)) stop("Nessuna citt\u00e0 per questa selezione.")
    if (!length(input$curve)) stop("Selezionare almeno una curva.")
    as_of <- input$asof
    from <- input$periodo[1]
    lead_max <- if (as_of >= Sys.Date()) 16 else FC_MAX_LEAD
    to <- min(input$periodo[2], as_of + lead_max)
    if (from >= to) stop("Periodo non valido.")
    set.seed(1)
    out <- withProgress(message = "Scarico le temperature", value = .3, {
      combined_series(cds, from, to, as_of, input$curve, isTRUE(input$debias),
                      baseline = "project", scope = "cities")
    })
    if (is.null(out)) stop("Nessun dato disponibile per il periodo scelto.")
    nomi <- sort(CITIES[URAU_CODE %in% cds, name])
    c(out, list(codes = cds, from = from, to = to, as_of = as_of,
                titolo = if (input$ambito == "reg") input$regione else input$citta,
                citta = nomi))
  })

  pal <- function(n) setNames(COL[seq_len(length(n))], n)
  tema <- function() theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  output$g_giorno <- renderPlot({
    r <- risultato()
    plot_deaths(
      r$daily, r$pred, r$as_of,
      title = sprintf("Daily deaths and heat-attributable prediction \u2014 %s", r$titolo),
      exp_label = "Expected without heat",
      subtitle_extra = if (r$from < BASELINE_FROM)
        sprintf(paste("Before %s the expectation is the seasonal profile of the last",
                      "twelve months of data, not a fitted baseline"),
                format(BASELINE_FROM, "%d %b %Y")) else NULL,
      caption = paste0(
        "Cities: ", paste(r$citta, collapse = ", "), ". ",
        "Observed: Istat daily deaths by municipality of residence, ages 20+. ",
        "Expected: seasonal quasi-Poisson model fitted on days without heat, with the ",
        "cold effect added back. Predicted: expected plus heat deaths from each ",
        "exposure-response function. Temperatures: ERA5-Land up to the issue date ",
        "(vertical line), Open-Meteo forecast afterwards."),
      ylab = "Deaths per day (ages 20+)")
  })

  output$g_cum <- renderPlot({
    r <- risultato()
    d <- copy(r$daily)[!is.na(observed)][order(date)]
    shiny::validate(shiny::need(nrow(d), "Nessun decesso osservato nel periodo."))
    d[, cum_obs := cumsum(observed - expected)]
    p <- ggplot(d, aes(date, cum_obs)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_line(aes(colour = "Observed excess"), linewidth = 0)
    if (nrow(r$pred)) {
      pp <- merge(r$pred, r$daily[, .(date, expected)], by = "date")
      pp <- pp[date %in% d$date][order(curve, date)]
      pp[, cum := cumsum(pred7 - frollmean(expected, 7, align = "center")), by = curve]
      p <- p + geom_line(data = pp[!is.na(cum)], aes(date, cum, colour = curve),
                         linewidth = .9)
    }
    p + scale_colour_manual(values = c("Observed excess" = "#1F2933",
                                       setNames(c("#B2182B", "#1B7837", "#762A83", "#B35806")[
                                         seq_len(uniqueN(r$pred$curve))],
                                         unique(r$pred$curve)))) +
      labs(title = sprintf("Cumulative heat deaths \u2014 %s", r$titolo),
           subtitle = "Observed excess over the expectation, and what each curve predicts",
           x = NULL, y = "Cumulative deaths", colour = NULL) + tema()
  })

  output$g_temp <- renderPlot({
    r <- risultato()
    t <- merge(r$temps, CITIES[, .(code = URAU_CODE, name)], by = "code")
    ggplot(t, aes(date, tmean, colour = name)) + geom_line(linewidth = .7) +
      labs(title = "Temperatura media giornaliera", subtitle = r$fonte,
           x = NULL, y = "\u00b0C") + tema() +
      theme(legend.position = if (uniqueN(t$name) > 12) "none" else "bottom")
  })

  output$g_mappa <- renderPlot({
    r <- risultato()
    per_city <- merge(CITIES[URAU_CODE %in% r$codes], r$by_city,
                      by.x = "URAU_CODE", by.y = "code", all.x = TRUE)
    per_city[is.na(att), att := 0]
    bb <- sf::st_bbox(sf::st_as_sf(per_city, coords = c("lon", "lat"), crs = 4326))
    m <- ggplot() + geom_sf(data = BUNDLE$regions, fill = "grey97", colour = "grey80",
                            linewidth = .2)
    if (!is.null(BUNDLE$city_geom))
      m <- m + geom_sf(data = BUNDLE$city_geom[BUNDLE$city_geom$URAU_CODE %in% r$codes, ],
                       fill = "#B2182B", colour = NA, alpha = .35)
    m + geom_point(data = per_city, aes(lon, lat, size = att), colour = "#B2182B",
                   alpha = .85) + scale_size_area(max_size = 12) +
      geom_text(data = per_city, aes(lon, lat, label = name), size = 3, vjust = -1.1) +
      coord_sf(xlim = c(bb["xmin"] - .6, bb["xmax"] + .6),
               ylim = c(bb["ymin"] - .4, bb["ymax"] + .4)) +
      labs(title = sprintf("Heat deaths in the period \u2014 %s", r$titolo),
           size = "Deaths", x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_line(colour = "grey95"))
  })

  output$tab <- renderTable({
    r <- risultato()
    d <- r$daily
    out <- data.table(Quantity = "Observed deaths", Value = sum(d$observed, na.rm = TRUE))
    out <- rbind(out, data.table(Quantity = "Expected without heat",
                                 Value = round(sum(d$expected, na.rm = TRUE))))
    out <- rbind(out, data.table(Quantity = "Observed heat deaths (observed - expected)",
                                 Value = round(sum(d$observed - d$expected, na.rm = TRUE))))
    if (nrow(r$pred)) {
      pp <- merge(r$pred, d[, .(date, expected)], by = "date")
      pp[, e7 := frollmean(expected, 7, align = "center")]
      s2 <- pp[!is.na(e7), .(Value = round(sum(pred7 - e7))), by = .(Quantity = curve)]
      s2[, Quantity := paste("Predicted heat deaths:", Quantity)]
      out <- rbind(out, s2)
    }
    out
  })

  output$copertura <- renderText({
    r <- tryCatch(risultato(), error = function(e) NULL)
    if (is.null(r)) return("")
    sprintf("Dati Istat fino al %s. Periodo mostrato: %s \u2013 %s.",
            format(BUNDLE$data_end, "%d/%m/%Y"), format(r$from, "%d/%m"),
            format(r$to, "%d/%m/%Y"))
  })

  output$bundle <- renderText(BUNDLE_INFO)

  output$msg <- renderUI({
    r <- tryCatch(risultato(), error = function(e) e)
    if (inherits(r, "error")) div(class = "errore", conditionMessage(r)) else NULL
  })

  # ---- scheda Scenari -------------------------------------------------------
  observe({
    nc <- c(PUB, unique(unlist(lapply(BUNDLE$curves, function(x) names(x$recalibrated)))))
    sel <- if (length(nc) > 1) nc[length(nc)] else nc[1]   # default: curve pi\u00f9 recenti
    updateSelectInput(session, "sc_curva", choices = setNames(nc, it_label(nc)),
                      selected = sel)
  })

  scenari <- eventReactive(input$vai, {
    cds <- codici()
    if (!length(cds)) stop("Nessuna citt\u00e0 per questa selezione.")
    as_of <- input$asof
    lead_max <- if (as_of >= Sys.Date()) 16 else FC_MAX_LEAD
    horizon <- max(1L, min(as.integer(input$periodo[2] - as_of), lead_max))
    to <- as_of + horizon
    set.seed(1)
    temps <- withProgress(message = "Scarico le previsioni", value = 0, {
      setNames(lapply(cds, function(cd) {
        incProgress(1 / length(cds), detail = CITIES[URAU_CODE == cd, name])
        city_temperature(cd, "fc", as_of, as_of + 1, to, isTRUE(input$debias))
      }), cds)
    })
    temps <- Filter(function(x) !is.null(x) && nrow(x), temps)
    if (!length(temps))
      stop("Nessuna previsione disponibile per il ", format(as_of, "%d/%m/%Y"),
           if (as_of < Sys.Date())
             ": l'archivio delle emissioni precedenti copre gli ultimi anni e al massimo 7 giorni di anticipo."
           else ".")
    cds <- names(temps)
    ons <- vapply(cds, function(cd) onset_probability(temps[[cd]], cd)$p_onset, 0)
    curva <- input$sc_curva %||% PUB
    tab <- scenario_table(cds, temps, curva)
    paths <- scenario_paths(cds, as_of, horizon, input$sc_back,
                            curva, isTRUE(input$debias))
    list(tab = tab, paths = paths, onset = ons, cds = cds, as_of = as_of,
         titolo = if (input$ambito == "reg") input$regione else input$citta,
         giorni = horizon,
         curva = it_label(curva),
         tmax = max(vapply(temps, function(x) max(x$tmean), 0)))
  })

  output$sc_testo <- renderUI({
    r <- tryCatch(scenari(), error = function(e) e)
    if (inherits(r, "error")) return(div(class = "errore", conditionMessage(r)))
    p_any <- 1 - prod(1 - r$onset, na.rm = TRUE)
    n <- nrow(r$tab)
    f <- function(i) sprintf("%s (%s \u2013 %s)", round(r$tab$stima[i]),
                             round(r$tab$lo[i]), round(r$tab$hi[i]))
    HTML(sprintf(paste0(
      "<p><b>%s.</b> Probabilit\u00e0 che un'ondata di calore inizi entro %d giorni: ",
      "<b>%s</b> (almeno una citt\u00e0; massimo previsto %.1f \u00b0C).</p>",
      "<p>Con le previsioni attuali i decessi attribuibili al caldo sono <b>%s</b>. ",
      "Se il caldo proseguisse, si aggiungerebbero fino a <b>%s</b> complessivi.</p>",
      "<p class='nota'>Curve: %s. Intervalli al 95%%: incertezza della curva e del ",
      "livello atteso. Parte di questi decessi \u00e8 anticipata di pochi giorni ",
      "piuttosto che aggiunta: la validazione mostra un eccesso molto minore nelle ",
      "due settimane successive rispetto a quanto previsto dalle curve.</p>"),
      r$titolo, format(r$as_of, "%d/%m/%Y"), r$giorni, pct_fmt(p_any), r$tmax,
      f(1), f(n), r$curva))
  })

  output$sc_tab <- renderTable({
    r <- scenari()
    data.table(Scenario = r$tab$scenario, Giorni = r$tab$giorni,
               `Decessi attribuibili` = sprintf("%s (%s \u2013 %s)", round(r$tab$stima),
                                                round(r$tab$lo), round(r$tab$hi)))
  })

  output$sc_plot <- renderPlot({
    r <- scenari()
    shiny::validate(shiny::need(!is.null(r$paths), "Scenari non disponibili."))
    d <- r$paths$paths
    extra <- setdiff(unique(d$scenario), c("Ricostruzione", "Previsione"))
    extra <- extra[order(as.numeric(gsub("\\D", "", extra)))]   # +5 before +10
    ord <- c("Ricostruzione", "Previsione", extra)
    d[, scenario := factor(scenario, ord)]
    cols <- setNames(c("#1F2933", "#B2182B", "#E08214", "#8073AC", "#4393C3")[seq_along(ord)],
                     ord)
    ggplot(d, aes(date, cum, colour = scenario, fill = scenario)) +
      geom_vline(xintercept = as.numeric(r$as_of), colour = "grey45", linetype = "22") +
      annotate("text", x = r$as_of, y = Inf, vjust = 1.5, hjust = -0.05, size = 3.4,
               colour = "grey35", label = "inizio previsione") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
      labs(title = sprintf("Decessi cumulati attribuibili al caldo \u2014 %s", r$titolo),
           subtitle = paste("Ricostruzione ERA5-Land fino alla data di emissione,",
                            "poi previsione e scenari di durata; bande al 95%"),
           x = NULL, y = "Decessi cumulati", colour = NULL, fill = NULL) + tema()
  })

  # ---- scheda EuroMOMO ------------------------------------------------------
  momo <- reactive({
    reg <- BUNDLE$region
    if (is.null(reg) || is.null(BUNDLE$flumomo)) return(NULL)
    from <- input$mm_periodo[1]; to <- input$mm_periodo[2]
    if (is.na(from) || is.na(to) || from >= to) return(NULL)
    fm <- if (input$mm_ambito == "reg") BUNDLE$flumomo$region_all else BUNDLE$flumomo$cities
    if (is.null(fm)) return(NULL)
    if (input$mm_ambito == "reg") {
      if (is.null(reg$daily)) return(NULL)
      d <- as.data.table(reg$daily)[date >= from & date <= to, .(date, observed = deaths)]
    } else {
      # observed over the whole period from the exported city series; our own
      # baseline exists only for the test period and is merged where available
      obs <- if (!is.null(reg$cities_daily))
        as.data.table(reg$cities_daily)[date >= from & date <= to, .(date, observed = deaths)]
      else DAILY[URAU_CODE %in% reg$cities & date >= from & date <= to,
                 .(observed = sum(deaths)), by = date]
      ours <- DAILY[URAU_CODE %in% reg$cities & date >= from & date <= to,
                    .(ours = sum(B * fifelse(tmean < mmt, exp(logrr), 1))), by = date]
      d <- merge(obs, ours, by = "date", all.x = TRUE)
      if (all(is.na(d$ours))) d[, ours := NULL]
    }
    if (!nrow(d)) return(NULL)
    w <- as.data.table(fm)[order(week_start)]
    # the FluMOMO age bands differ slightly from ours: rescale over the same period
    ref <- w[week_start >= from & week_start <= to]
    sc <- if (nrow(ref) && sum(ref$deaths) > 0) sum(d$observed) / sum(ref$deaths) else 1
    d[, flumomo := sc * stats::splinefun(as.numeric(w$week_start) + 3,
                                         w$expected / 7)(as.numeric(date))]
    setorder(d, date)
    d[, `:=`(obs7 = frollmean(observed, 7, align = "center"),
             flu7 = frollmean(flumomo, 7, align = "center"))]
    if ("ours" %in% names(d) && any(!is.na(d$ours)))
      d[, ours7 := frollmean(ours, 7, align = "center")]
    d[!is.na(obs7)]
  })

  momo_series <- reactive({
    reg <- BUNDLE$region
    shiny::validate(
      shiny::need(!is.null(BUNDLE$flumomo),
                  "app_data.rds has no FluMOMO results: run targets::tar_make()."),
      shiny::need(input$mm_ambito != "cit" || !is.null(reg$cities_daily),
                  paste("app_data.rds predates the export of the validated cities'",
                        "deaths: run targets::tar_make() and copy the new file next",
                        "to app.R.")))
    from <- input$mm_periodo[1]; to <- input$mm_periodo[2]
    shiny::validate(shiny::need(!is.na(from) && !is.na(to) && from < to, "Periodo non valido."))
    cds <- reg$cities
    combined_series(cds, from, to, min(to, Sys.Date()), input$curve %||% PUB, FALSE,
                    baseline = "flumomo",
                    scope = if (input$mm_ambito == "reg") "region" else "cities")
  })

  output$g_momo <- renderPlot({
    r <- momo_series()
    shiny::validate(shiny::need(!is.null(r) && nrow(r$daily), "Nessun dato per il periodo."))
    reg <- BUNDLE$region
    scope_txt <- if (input$mm_ambito == "reg") "all municipalities"
                 else paste(sort(CITIES[URAU_CODE %in% reg$cities, name]), collapse = ", ")
    plot_deaths(
      r$daily, r$pred, NULL,
      title = sprintf("Daily deaths against the FluMOMO baseline \u2014 %s", reg$name),
      exp_label = "Expected (FluMOMO baseline)",
      caption = paste0(
        "Population: ", scope_txt, ". Observed: Istat daily deaths by municipality of ",
        "residence. Expected: FluMOMO baseline (EuroMOMO, version 4.2, official code) ",
        "in absolute counts, that is the deaths expected with no influenza activity and ",
        "no extreme temperature, spread over days. Predicted: the FluMOMO expectation ",
        "plus the heat deaths of each exposure-response function."),
      ylab = if (input$mm_ambito == "reg") "Deaths per day" else "Deaths per day (ages 20+)")
  })

  output$tab_momo <- renderTable({
    r <- momo_series()
    if (is.null(r) || !nrow(r$daily)) return(NULL)
    d <- r$daily
    data.table(Quantity = c("Observed deaths", "Expected (FluMOMO)",
                            "Difference (observed - expected)"),
               Value = c(sum(d$observed, na.rm = TRUE), round(sum(d$expected, na.rm = TRUE)),
                         round(sum(d$observed - d$expected, na.rm = TRUE))))
  })
}

shinyApp(ui, server)
