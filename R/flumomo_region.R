# =============================================================================
#  Region-level excess mortality with the official FluMOMO code, as pipeline
#  steps: weather input, estimation, and the Italian-language chart.
# =============================================================================

#' Daily temperature by NUTS-3 for a region
#'
#' Downloads ERA5-Land at each province capital and returns the file the
#' official FluMOMO code expects (`date`, `pop3`, `NUTS3`, `temp`); it does the
#' population weighting and the weekly aggregation itself.
#'
#' @param cfg Configuration list (uses `flumomo$provinces`).
#' @param from,to Period.
#' @param path Output CSV.
#' @return `path`.
build_flumomo_weather <- function(cfg, from, to, path = NULL) {
  path <- path %||% file.path(cfg$cache_dir,
                              sprintf("flumomo_weather_%s.csv", cfg$flumomo$country_code))
  if (file.exists(path)) {
    d <- data.table::fread(path)
    if (min(as.Date(d$date)) <= from && max(as.Date(d$date)) >= to) return(path)
  }
  p <- data.table::as.data.table(cfg$flumomo$provinces)
  wx <- data.table::rbindlist(lapply(seq_len(nrow(p)), function(i) {
    d <- fetch_era5_point(p$lat[i], p$lon[i], from, to, cfg$era5_timezone)
    d[, .(date, pop3 = p$pop3[i], NUTS3 = p$NUTS3[i], temp = tmean)]
  }))
  ensure_dir(dirname(path))
  data.table::fwrite(wx, path)
  path
}

#' Run the official FluMOMO for the configured region
#'
#' Prepares the inputs from the Istat file, the weather file and the influenza
#' activity series, then runs `Estimation_v42.R` unchanged.
#'
#' @param cfg Configuration list.
#' @param istat_file Istat daily CSV.
#' @param weather_file Output of [build_flumomo_weather()].
#' @param influenza_file Output of [build_influenza_activity()] (or `NULL`).
#' @param procom Comune codes to use instead of the whole provinces (for the
#'   run restricted to the validated cities).
#' @param label Name used for the output files when `procom` is given.
#' @return data.table of official results.
run_region_flumomo <- function(cfg, istat_file, weather_file, influenza_file = NULL,
                               procom = NULL, label = NULL) {
  fm <- cfg$flumomo
  deaths <- istat_weekly_euromomo(istat_file, fm$provinces$prov, fm$years, procom)
  weather <- data.table::fread(weather_file)[, date := as.Date(date)]
  ia <- if (!is.null(influenza_file) && file.exists(influenza_file)) {
    x <- data.table::fread(influenza_file)
    data.table::setnames(x, names(x), sub("^ia$", "IA", tolower(names(x))))
    x[, .(year, week, IA)]
  } else NULL
  region <- label %||% fm$region
  work <- ensure_dir(if (is.null(label)) fm$work_dir else paste0(fm$work_dir, "_", label))
  flumomo_write_inputs(work, deaths, weather, ia, fm$country_code)
  flumomo_run(fm$code_dir, work, country = region, country_code = fm$country_code,
              start = min(weather$date), end = max(weather$date),
              ia_lags = fm$ia_lags, et_lags = fm$et_lags,
              ia_restricted = fm$ia_restricted)
}

#' Chart of observed and expected mortality for a region (Italian labels)
#'
#' Observed daily deaths against the FluMOMO baseline in counts, both as
#' 7-day moving averages, with the difference shaded.
#'
#' @param res Output of [run_region_flumomo()].
#' @param istat_file Istat daily CSV.
#' @param cfg Configuration list.
#' @param year Calendar year to draw.
#' @param path Output PNG.
#' @return `path`.
plot_region_excess <- function(res, istat_file, cfg, year = NULL, path = NULL) {
  fm <- cfg$flumomo
  year <- year %||% fm$chart_year
  path <- path %||% sprintf("mortalita_%s_%d.png", tolower(fm$region), year)
  daily <- istat_daily_area(istat_file, fm$provinces$prov, age_min = 0, years = year)
  d <- daily[data.table::year(date) == year][order(date)]
  d[, expected := flumomo_daily(res, fm$agegrp, date, "EB")]
  d[, `:=`(obs7 = data.table::frollmean(deaths, 7, align = "center"),
           exp7 = data.table::frollmean(expected, 7, align = "center"))]
  d <- d[!is.na(obs7)]
  total <- round(sum(d$deaths - d$expected))
  period <- sprintf("%s - %s", format(min(d$date), "%d/%m"), format(max(d$date), "%d/%m/%Y"))
  months_it <- c("gen", "feb", "mar", "apr", "mag", "giu",
                 "lug", "ago", "set", "ott", "nov", "dic")
  wrap <- function(x, n = 135) paste(strwrap(x, width = n), collapse = "\n")
  has_ia <- any(res$IA != 0, na.rm = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = exp7, ymax = pmax(obs7, exp7)),
                         fill = "#C0392B", alpha = .28) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmin(obs7, exp7), ymax = exp7),
                         fill = "#3A8DAE", alpha = .18) +
    ggplot2::geom_line(ggplot2::aes(y = exp7, colour = "Attesi (baseline FluMOMO)"),
                       linewidth = .8, linetype = "22") +
    ggplot2::geom_line(ggplot2::aes(y = obs7, colour = "Osservati"), linewidth = 1) +
    ggplot2::scale_colour_manual(
      values = c("Osservati" = "#1F2933", "Attesi (baseline FluMOMO)" = "#2C7FB8"),
      breaks = c("Osservati", "Attesi (baseline FluMOMO)")) +
    ggplot2::scale_x_date(date_breaks = "1 month",
                          labels = function(x) months_it[as.integer(format(x, "%m"))],
                          expand = ggplot2::expansion(mult = c(.01, .01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.05, .12))) +
    ggplot2::labs(
      title = sprintf("Mortalit\u00e0 giornaliera in %s, %d", fm$region, year),
      subtitle = wrap(sprintf(paste("Tutti i comuni. Medie mobili a 7 giorni.",
                                    "In rosso i giorni sopra l'attesa, in azzurro quelli",
                                    "sotto (saldo %s, %s)"),
                              format(total, big.mark = ".", decimal.mark = ","), period), 105),
      x = NULL, y = "Decessi al giorno", colour = NULL,
      caption = paste0(
        wrap(sprintf("Decessi giornalieri Istat per comune di residenza, province: %s.",
                     paste(fm$provinces$prov, collapse = ", "))), "\n",
        wrap(paste0("Attesi: baseline FluMOMO (EuroMOMO, versione 4.2, codice ufficiale) ",
                    "in conteggi assoluti, cio\u00e8 i decessi attesi in assenza di ",
                    "attivit\u00e0 influenzale e di temperature estreme, ripartiti sui ",
                    "giorni. Temperatura: ERA5-Land nei capoluoghi, pesata per ",
                    "popolazione provinciale (NUTS-3).",
                    if (!has_ia) " Attenzione: nessun dato di attivit\u00e0 influenzale (IA = 0)."
                    else "")))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25)),
      plot.subtitle = ggplot2::element_text(colour = "grey30",
                                            margin = ggplot2::margin(b = 12)),
      plot.caption = ggplot2::element_text(colour = "grey45",
                                           size = ggplot2::rel(.72), hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      legend.position = "top", legend.justification = "left",
      legend.key.width = ggplot2::unit(28, "pt"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(colour = "grey30",
                                           margin = ggplot2::margin(r = 8)))
  ggplot2::ggsave(path, p, width = 10, height = 5.2, dpi = 200, bg = "white")
  path
}
