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

# ---- utilita' numeriche -----------------------------------------------------
basis <- function(t, spec) {
  suppressWarnings(unclass(splines::bs(t, knots = spec$knots, degree = spec$degree,
                                       Boundary.knots = spec$bound)))
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

# ================================ UI =========================================
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
          radioButtons("fonte", "Temperature",
                       c("Previsione (Open-Meteo)" = "fc",
                         "Ricostruzione ERA5-Land" = "era5")),
          dateInput("asof", "Data di riferimento (macchina del tempo)",
                    value = Sys.Date(), max = Sys.Date(), format = "dd/mm/yyyy",
                    language = "it"),
          conditionalPanel("input.fonte == 'fc'",
                           sliderInput("orizzonte", "Giorni di previsione", 1, 16, 7, 1),
                           checkboxInput("debias", "Correggi la distorsione delle previsioni",
                                         TRUE),
                           div(class = "nota",
                               "Con una data passata si usano le previsioni emesse allora",
                               "(massimo 7 giorni di anticipo).")),
          conditionalPanel("input.fonte == 'era5'",
                           dateRangeInput("periodo", "Periodo", start = Sys.Date() - 60,
                                          end = Sys.Date() - 6, max = Sys.Date() - 6,
                                          format = "dd/mm/yyyy", language = "it",
                                          separator = " \u2013 ")),
          hr(),
          radioButtons("ambito", "Ambito", c("Regione" = "reg", "Singola citt\u00e0" = "city")),
          conditionalPanel("input.ambito == 'reg'",
                           selectInput("regione", "Regione", choices = REGIONI,
                                       selected = if ("Lombardia" %in% REGIONI) "Lombardia"
                                       else REGIONI[1])),
          conditionalPanel("input.ambito == 'city'",
                           selectInput("citta", "Citt\u00e0", choices = sort(CITIES$name))),
          checkboxGroupInput("curve", "Curve", choices = NULL),
          actionButton("vai", "Calcola", class = "btn-primary"),
          br(), br(), uiOutput("msg"), div(class = "nota", textOutput("copertura")),
          br(), div(class = "nota", style = "word-break:break-all;", textOutput("bundle"))
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel("Giornalieri", br(), plotOutput("g_giorno", height = "430px")),
            tabPanel("Cumulati", br(), plotOutput("g_cum", height = "430px")),
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
    if (input$fonte == "fc") {
      lead_max <- if (as_of >= Sys.Date()) 16 else FC_MAX_LEAD
      from <- as_of + 1; to <- as_of + min(input$orizzonte, lead_max)
    } else {
      from <- input$periodo[1]; to <- min(input$periodo[2], Sys.Date() - 6)
    }
    if (from > to) stop("Periodo non valido.")
    set.seed(1)
    out <- withProgress(message = "Scarico le temperature", value = 0, {
      lapply(seq_along(cds), function(i) {
        incProgress(1 / length(cds), detail = CITIES[URAU_CODE == cds[i], name])
        tt <- city_temperature(cds[i], input$fonte, as_of, from, to,
                               isTRUE(input$debias))
        if (!nrow(tt)) return(NULL)
        att <- lapply(setNames(input$curve, input$curve), function(nm)
          city_attributable(cds[i], tt, nm, NSIM))
        att <- Filter(Negate(is.null), att)
        if (!length(att)) return(NULL)
        list(code = cds[i], temp = tt, att = att,
             obs = data.table(date = tt$date, baseline = att[[1]]$baseline,
                              observed = att[[1]]$observed))
      })
    })
    out <- Filter(Negate(is.null), out)
    if (!length(out)) stop("Nessun dato disponibile per il periodo scelto.")
    curve_names <- unique(unlist(lapply(out, function(o) names(o$att))))
    dates <- sort(unique(unlist(lapply(out, function(o) o$temp$date))))
    dates <- as.Date(dates, origin = "1970-01-01")
    # somma fra citta', per curva: punto e distribuzione
    daily <- rbindlist(lapply(curve_names, function(nm) {
      pt <- rep(0, length(dates)); dr <- matrix(0, NSIM, length(dates))
      for (o in out) if (!is.null(o$att[[nm]])) {
        j <- match(o$temp$date, dates)
        pt[j] <- pt[j] + o$att[[nm]]$point
        dr[, j] <- dr[, j] + o$att[[nm]]$draws
      }
      data.table(date = dates, curva = it_label(nm), stima = pt,
                 lo = apply(dr, 2, quantile, .025), hi = apply(dr, 2, quantile, .975),
                 cum = cumsum(pt),
                 cum_lo = apply(t(apply(dr, 1, cumsum)), 2, quantile, .025),
                 cum_hi = apply(t(apply(dr, 1, cumsum)), 2, quantile, .975))
    }))
    obs <- rbindlist(lapply(out, `[[`, "obs"))[, .(baseline = sum(baseline),
                                                   observed = sum(observed)), by = date]
    list(daily = daily, obs = obs,
         temp = rbindlist(lapply(out, function(o) cbind(code = o$code, o$temp))),
         codes = vapply(out, `[[`, "", "code"), from = from, to = to, as_of = as_of,
         titolo = if (input$ambito == "reg") input$regione else input$citta,
         fonte = if (input$fonte == "fc")
           sprintf("previsioni emesse il %s", format(as_of, "%d/%m/%Y"))
         else "ricostruzione ERA5-Land")
  })

  pal <- function(n) setNames(COL[seq_len(length(n))], n)
  tema <- function() theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  output$g_giorno <- renderPlot({
    r <- risultato()
    o <- r$obs[!is.na(observed), .(date, v = observed - baseline)]
    p <- ggplot(r$daily, aes(date, stima, colour = curva, fill = curva)) +
      geom_hline(yintercept = 0, colour = "grey70")
    if (nrow(o)) p <- p + geom_col(data = o, aes(date, v), inherit.aes = FALSE,
                                   fill = "grey78", width = 1)
    p + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) +
      geom_line(linewidth = .9) +
      scale_colour_manual(values = pal(unique(r$daily$curva))) +
      scale_fill_manual(values = pal(unique(r$daily$curva))) +
      labs(title = sprintf("Decessi giornalieri attribuibili al caldo \u2014 %s", r$titolo),
           subtitle = sprintf("%s; %s%s", r$fonte,
                              if (HAS_VCOV && HAS_BSD) "intervalli al 95%"
                              else if (HAS_VCOV || HAS_BSD)
                                "intervalli al 95% (parziali: manca una componente)"
                              else "intervalli non disponibili (aggiornare app_data.rds)",
                              if (nrow(o)) ". Barre grigie: eccesso osservato (Istat)" else ""),
           x = NULL, y = "Decessi al giorno") + tema()
  })

  output$g_cum <- renderPlot({
    r <- risultato()
    o <- r$obs[!is.na(observed)][order(date)]
    p <- ggplot(r$daily, aes(date, cum, colour = curva, fill = curva)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_ribbon(aes(ymin = cum_lo, ymax = cum_hi), alpha = .18, colour = NA) +
      geom_line(linewidth = 1)
    if (nrow(o)) {
      o[, cum := cumsum(observed - baseline)]
      p <- p + geom_line(data = o, aes(date, cum), inherit.aes = FALSE,
                         colour = "#1F2933", linetype = "22", linewidth = .9)
    }
    p + scale_colour_manual(values = pal(unique(r$daily$curva))) +
      scale_fill_manual(values = pal(unique(r$daily$curva))) +
      labs(title = sprintf("Decessi cumulati attribuibili al caldo \u2014 %s", r$titolo),
           subtitle = if (nrow(o)) "Linea nera: eccesso osservato cumulato" else r$fonte,
           x = NULL, y = "Decessi cumulati") + tema()
  })

  output$g_temp <- renderPlot({
    r <- risultato()
    t <- merge(r$temp, CITIES[, .(code = URAU_CODE, name)], by = "code")
    ggplot(t, aes(date, tmean, colour = name)) + geom_line(linewidth = .7) +
      labs(title = "Temperatura media giornaliera", subtitle = r$fonte,
           x = NULL, y = "\u00b0C") + tema() +
      theme(legend.position = if (uniqueN(t$name) > 12) "none" else "bottom")
  })

  output$g_mappa <- renderPlot({
    r <- risultato()
    tot <- r$daily[curva == curva[1]]
    per_city <- CITIES[URAU_CODE %in% r$codes]
    per_city[, v := sum(tot$stima) / .N]
    bb <- sf::st_bbox(sf::st_as_sf(per_city, coords = c("lon", "lat"), crs = 4326))
    m <- ggplot() + geom_sf(data = BUNDLE$regions, fill = "grey97", colour = "grey80",
                            linewidth = .2)
    if (!is.null(BUNDLE$city_geom))
      m <- m + geom_sf(data = BUNDLE$city_geom[BUNDLE$city_geom$URAU_CODE %in% r$codes, ],
                       fill = "#B2182B", colour = NA, alpha = .35)
    m + geom_point(data = per_city, aes(lon, lat), colour = "#B2182B", size = 3, alpha = .85) +
      geom_text(data = per_city, aes(lon, lat, label = name), size = 3, vjust = -1.1) +
      coord_sf(xlim = c(bb["xmin"] - .6, bb["xmax"] + .6),
               ylim = c(bb["ymin"] - .4, bb["ymax"] + .4)) +
      labs(title = sprintf("Citt\u00e0 incluse \u2014 %s", r$titolo), x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_line(colour = "grey95"))
  })

  output$tab <- renderTable({
    r <- risultato()
    a <- r$daily[, .(`Decessi attribuibili` = sprintf("%.1f (%.1f \u2013 %.1f)",
                                                      sum(stima), sum(lo), sum(hi))),
                 by = .(Curva = curva)]
    o <- r$obs[!is.na(observed)]
    if (nrow(o)) a <- rbind(a, data.table(Curva = "Eccesso osservato (Istat)",
                                          `Decessi attribuibili` =
                                            sprintf("%.1f", sum(o$observed - o$baseline))))
    a
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
      d <- DAILY[URAU_CODE %in% reg$cities & date >= from & date <= to,
                 .(observed = sum(deaths),
                   ours = sum(B * fifelse(tmean < mmt, exp(logrr), 1))), by = date]
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
    if ("ours" %in% names(d)) d[, ours7 := frollmean(ours, 7, align = "center")]
    d[!is.na(obs7)]
  })

  output$g_momo <- renderPlot({
    d <- momo()
    # shiny::validate explicitly: jsonlite masks validate() when attached later
    shiny::validate(shiny::need(
      !is.null(d) && nrow(d),
      paste("Dati non disponibili per il periodo scelto:",
            "verificare le date e che i target FluMOMO siano stati eseguiti.")))
    lab_flu <- "Attesi FluMOMO"; lab_our <- "Attesi senza caldo (progetto)"
    cols <- c("Osservati" = "#1F2933", "#D95F02", "#2C7FB8")
    names(cols)[2:3] <- c(lab_flu, lab_our)
    p1 <- "Decessi al giorno (media mobile 7 giorni)"
    p2 <- "Eccesso cumulato"
    line <- function(v, serie, pannello) {
      if (is.null(v) || !length(v)) return(NULL)
      data.table(date = d$date, serie = serie, v = as.numeric(v), pannello = pannello)
    }
    parts <- list(line(d$obs7, "Osservati", p1), line(d$flu7, lab_flu, p1),
                  line(cumsum(d$observed - d$flumomo), lab_flu, p2))
    if ("ours" %in% names(d)) parts <- c(parts, list(
      line(d$ours7, lab_our, p1), line(cumsum(d$observed - d$ours), lab_our, p2)))
    pd <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE)[!is.na(v)]
    ggplot(pd, aes(date, v, colour = serie)) +
      geom_hline(yintercept = 0, colour = "grey85") +
      geom_line(linewidth = .9) +
      facet_wrap(~ pannello, ncol = 1, scales = "free_y") +
      scale_colour_manual(values = cols) +
      labs(title = sprintf("%s, %s \u2013 %s \u2014 %s", BUNDLE$region$name,
                           format(min(d$date), "%d/%m/%Y"), format(max(d$date), "%d/%m/%Y"),
                           if (input$mm_ambito == "reg") "tutti i comuni"
                           else "citt\u00e0 validate"),
           x = NULL, y = NULL) + tema() +
      theme(strip.text = element_text(hjust = 0, face = "bold", colour = "grey25"))
  })

  output$tab_momo <- renderTable({
    d <- momo()
    if (is.null(d) || !nrow(d)) return(NULL)
    out <- data.table(Attesa = "FluMOMO (senza temperature estreme)",
                      Osservati = sum(d$observed), Attesi = round(sum(d$flumomo)),
                      Eccesso = round(sum(d$observed - d$flumomo)))
    if ("ours" %in% names(d))
      out <- rbind(out, data.table(Attesa = "Modello del progetto (senza caldo)",
                                   Osservati = sum(d$observed), Attesi = round(sum(d$ours)),
                                   Eccesso = round(sum(d$observed - d$ours))))
    out
  })
}

shinyApp(ui, server)
