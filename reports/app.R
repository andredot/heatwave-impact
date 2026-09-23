# =============================================================================
#  Mortalita' attribuibile al CALDO - citta' italiane
#  ---------------------------------------------------------------------------
#  Due modalita':
#    - Nowcasting: temperature osservate ERA5-Land (Open-Meteo Archive), per
#      guardare indietro o alla situazione attuale (ERA5-Land ha ~6 giorni di
#      ritardo).
#    - Previsione: temperature previste (Open-Meteo Forecast), fino a 16 giorni.
#
#  Curve esposizione-risposta:
#    - pubblicate: Masselot et al. (2023), Lancet Planetary Health;
#    - ricalibrate: stesso modello ristimato sui dati italiani recenti dalla
#      pipeline di validazione (le pubblicate sovrastimano il caldo di oggi).
#
#  Dove esistono i dati Istat, la mortalita' osservata viene mostrata accanto
#  alle stime, in valore giornaliero e cumulato.
#
#  Dati precalcolati in app_data.rds (prodotti da targets::tar_make()).
# =============================================================================
library(shiny)
library(data.table)
library(ggplot2)
library(jsonlite)
library(sf)

BUNDLE <- readRDS("app_data.rds")
CITIES <- as.data.table(BUNDLE$cities)
REGIONI <- sort(unique(CITIES$region))
COL <- c("Pubblicate (Masselot 2023)" = "#B2182B", "Osservati" = "#1F2933")
#' Etichette in italiano per le curve ricalibrate dalla pipeline
it_label <- function(x) sub("^Italian curves", "Curve italiane", x)
REC_COLS <- c("#1B7837", "#762A83", "#B35806")

# ---- utilita' ---------------------------------------------------------------
basis <- function(t, spec) {
  suppressWarnings(unclass(splines::bs(t, knots = spec$knots, degree = spec$degree,
                                       Boundary.knots = spec$bound)))
}

#' Rischio relativo di una curva rispetto alla sua MMT
rr_of <- function(t, spec, curve) {
  as.numeric(basis(t, spec) %*% curve$beta) -
    as.numeric(basis(curve$mmt, spec) %*% curve$beta)
}

get_json <- function(url, tries = 3, wait = 5) {
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

#' Temperature giornaliere per una citta'
#'
#' `mode = "now"` usa l'archivio ERA5-Land, `mode = "fc"` le previsioni.
fetch_temp <- function(lat, lon, from, to, mode, tz, model = NULL) {
  url <- if (mode == "now") {
    sprintf(paste0("https://archive-api.open-meteo.com/v1/archive?latitude=%.4f",
                   "&longitude=%.4f&start_date=%s&end_date=%s",
                   "&daily=temperature_2m_mean&models=era5_land&timezone=%s"),
            lat, lon, format(from), format(to), tz)
  } else {
    # stesso modello dell'archivio usato per stimare la distorsione
    sprintf(paste0("https://api.open-meteo.com/v1/forecast?latitude=%.4f",
                   "&longitude=%.4f&daily=temperature_2m_mean&forecast_days=16",
                   "&past_days=7&timezone=%s%s"), lat, lon, tz,
            if (!is.null(model) && nzchar(model)) paste0("&models=", model) else "")
  }
  js <- get_json(url)
  d <- data.table(date = as.Date(js$daily$time),
                  tmean = as.numeric(js$daily$temperature_2m_mean))
  d[, lead := as.integer(date - Sys.Date())]
  d[!is.na(tmean) & date >= from & date <= to]
}

#' Serie giornaliera di una citta' con decessi attribuibili al caldo
#'
#' Restituisce, per ogni giorno: i decessi attesi in assenza di caldo, i
#' decessi osservati (se disponibili) e i decessi attribuibili al caldo
#' secondo ogni insieme di curve.
city_series <- function(code, from, to, mode, sets, bundle, debias = TRUE) {
  ci <- CITIES[URAU_CODE == code]
  cu <- bundle$curves[[code]]
  tt <- fetch_temp(ci$lat, ci$lon, from, to, mode, bundle$era5_timezone,
                   bundle$forecast_model)
  if (!nrow(tt)) return(NULL)
  if (mode == "fc" && isTRUE(debias) && !is.null(bundle$fc_bias)) {
    fb <- as.data.table(bundle$fc_bias)[URAU_CODE == code]
    if (nrow(fb)) {
      # gli scarti sono stimati per i lead disponibili nell'archivio: fuori da
      # quell'intervallo si usa il lead piu' vicino
      tt[, `:=`(month = month(date),
                lead_c = pmin(pmax(lead, min(fb$lead)), max(fb$lead)))]
      tt <- merge(tt, fb[, .(month, lead_c = lead, bias)],
                  by = c("month", "lead_c"), all.x = TRUE)
      tt[is.na(bias), bias := 0][, tmean := tmean - bias]
      setorder(tt, date)
    }
  }
  if (isTRUE(bundle$bias_correct)) {
    dl <- as.data.table(bundle$delta)[URAU_CODE == code]
    if (!"month" %in% names(tt)) tt[, month := month(date)]
    tt <- merge(tt, dl[, .(month, delta)], by = "month", all.x = TRUE)
    tt[is.na(delta), delta := 0][, tmean := tmean - delta]
    setorder(tt, date)
  }

  # baseline: giorni con dati Istat; oltre, profilo stagionale dell'ultimo anno
  bs <- bundle$daily[URAU_CODE == code, .(B = sum(B), oss = sum(deaths)), by = date]
  bs[, doy := as.integer(format(date, "%j"))]
  clim <- bs[date > max(date) - 365, .(Bc = mean(B)), by = doy]
  x <- merge(tt[, .(date, tmean)], bs[, .(date, B, oss)], by = "date", all.x = TRUE)
  x[, doy := as.integer(format(date, "%j"))]
  x <- merge(x, clim, by = "doy", all.x = TRUE)[order(date)]
  x[is.na(B), B := Bc][, Bc := NULL]
  x <- x[!is.na(B)]
  if (!nrow(x)) return(NULL)

  # decessi attribuibili al caldo, per ogni insieme di curve
  age_share <- bundle$daily[URAU_CODE == code, .(s = sum(B)), by = agegroup]
  age_share[, s := s / sum(s)]
  res <- list()
  for (nm in names(sets)) {
    cv <- sets[[nm]]
    if (is.null(cv)) next
    if (identical(nm, "Pubblicate (Masselot 2023)")) {
      # curve per fascia d'eta': si pesano con la quota di baseline di ciascuna
      att <- rowSums(vapply(age_share$agegroup, function(g) {
        p <- cu$published[[g]]
        lr <- rr_of(x$tmean, cu$spec, p)
        share <- age_share$s[age_share$agegroup == g]
        share * x$B * pmax(exp(lr) - 1, 0) * (x$tmean >= p$mmt)
      }, numeric(nrow(x))))
      freddo <- rowSums(vapply(age_share$agegroup, function(g) {
        p <- cu$published[[g]]
        share <- age_share$s[age_share$agegroup == g]
        share * x$B * (exp(rr_of(x$tmean, cu$spec, p)) - 1) * (x$tmean < p$mmt)
      }, numeric(nrow(x))))
    } else {
      lr <- rr_of(x$tmean, cu$spec, cv)
      att <- x$B * pmax(exp(lr) - 1, 0) * (x$tmean >= cv$mmt)
      freddo <- x$B * (exp(lr) - 1) * (x$tmean < cv$mmt)
    }
    res[[nm]] <- data.table(date = x$date, curve = nm, attribuibili = att, freddo = freddo)
  }
  att <- rbindlist(res)
  # eccesso osservato: osservati meno attesi (freddo incluso, caldo escluso),
  # calcolato con le curve pubblicate per coerenza con la validazione
  ref <- att[curve == names(sets)[1], .(date, freddo)]
  obs <- merge(x[, .(date, B, oss, tmean)], ref, by = "date", all.x = TRUE)
  obs[is.na(freddo), freddo := 0][, attesi := B + freddo]
  obs[, eccesso := oss - attesi]
  list(citta = ci$name, temp = x[, .(date, tmean)], attrib = att,
       obs = obs[, .(date, tmean, oss, attesi, eccesso)])
}

# ================================ UI =========================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; }
    .titolo { font-weight: 700; font-size: 1.5rem; margin-bottom: .1rem; }
    .sottotitolo { color: #666; margin-bottom: 1rem; }
    .errore { color: #B00020; font-weight: 600; }
    .nota { color: #666; font-size: .85rem; }
  "))),
  div(class = "titolo", "Decessi attribuibili al caldo \u2014 citt\u00e0 italiane"),
  div(class = "sottotitolo",
      "Curve esposizione-risposta di Masselot et al. (2023) e loro ricalibrazione sui dati italiani recenti."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("modo", "Modalit\u00e0",
                   c("Nowcasting (ERA5-Land osservato)" = "now",
                     "Previsione (fino a 16 giorni)" = "fc")),
      radioButtons("ambito", "Ambito", c("Regione" = "reg", "Singola citt\u00e0" = "city")),
      conditionalPanel("input.ambito == 'reg'",
                       selectInput("regione", "Regione", choices = REGIONI,
                                   selected = if ("Lombardia" %in% REGIONI) "Lombardia" else REGIONI[1])),
      conditionalPanel("input.ambito == 'city'",
                       selectInput("citta", "Citt\u00e0", choices = sort(CITIES$name))),
      conditionalPanel("input.modo == 'fc'",
                       checkboxInput("debias", "Correggi la distorsione delle previsioni",
                                     value = TRUE)),
      conditionalPanel("input.modo == 'now'",
                       dateRangeInput("periodo", "Periodo",
                                      start = Sys.Date() - 120, end = Sys.Date() - 6,
                                      max = Sys.Date() - 6, format = "dd/mm/yyyy",
                                      language = "it", separator = " \u2013 ")),
      checkboxGroupInput("curve", "Curve", choices = NULL),
      actionButton("vai", "Calcola", class = "btn-primary"),
      br(), br(), uiOutput("msg"),
      div(class = "nota", textOutput("copertura"))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Giornalieri", br(), plotOutput("g_giorno", height = "420px")),
        tabPanel("Cumulati", br(), plotOutput("g_cum", height = "420px")),
        tabPanel("Temperatura", br(), plotOutput("g_temp", height = "420px")),
        tabPanel("Mappa", br(), plotOutput("g_mappa", height = "520px")),
        tabPanel("Tabella", br(), tableOutput("tab"))
      )
    )
  )
)

# ============================== SERVER =======================================
server <- function(input, output, session) {

  nomi_curve <- reactive({
    rec <- unique(unlist(lapply(BUNDLE$curves, function(x) names(x$recalibrated))))
    c("Pubblicate (Masselot 2023)", rec)
  })
  observe({
    nc <- nomi_curve()
    updateCheckboxGroupInput(session, "curve", choices = stats::setNames(nc, it_label(nc)),
                             selected = nc)
  })

  codici <- reactive({
    if (input$ambito == "reg") CITIES[region == input$regione, URAU_CODE]
    else CITIES[name == input$citta, URAU_CODE]
  })

  risultato <- eventReactive(input$vai, {
    cds <- codici()
    if (!length(cds)) stop("Nessuna citt\u00e0 disponibile per questa selezione.")
    if (input$modo == "now") {
      from <- input$periodo[1]; to <- min(input$periodo[2], Sys.Date() - 6)
    } else {
      from <- Sys.Date() - 7; to <- Sys.Date() + 15
    }
    if (from >= to) stop("Periodo non valido.")
    withProgress(message = "Scarico le temperature", value = 0, {
      out <- lapply(seq_along(cds), function(i) {
        incProgress(1 / length(cds), detail = CITIES[URAU_CODE == cds[i], name])
        sets <- list()
        cu <- BUNDLE$curves[[cds[i]]]
        for (nm in input$curve) {
          sets[[nm]] <- if (identical(nm, "Pubblicate (Masselot 2023)")) TRUE
          else cu$recalibrated[[nm]]
        }
        s <- city_series(cds[i], from, to, input$modo, sets, BUNDLE,
                         debias = isTRUE(input$debias))
        if (!is.null(s)) { s$attrib[, code := cds[i]]; s$obs[, code := cds[i]]
                           s$temp[, code := cds[i]] }
        s
      })
    })
    out <- Filter(Negate(is.null), out)
    if (!length(out)) stop("Nessun dato disponibile per il periodo scelto.")
    list(
      attrib = rbindlist(lapply(out, `[[`, "attrib"))[, curve := it_label(curve)][],
      obs    = rbindlist(lapply(out, `[[`, "obs")),
      temp   = rbindlist(lapply(out, `[[`, "temp")),
      codes  = cds, from = from, to = to,
      titolo = if (input$ambito == "reg") input$regione else input$citta)
  })

  palette_curve <- function(nomi) {
    names(COL) <- it_label(names(COL))
    p <- COL[intersect(names(COL), nomi)]
    rec <- setdiff(nomi, names(COL))
    if (length(rec)) p <- c(p, setNames(REC_COLS[seq_along(rec)], rec))
    p
  }

  tema <- function() theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  output$g_giorno <- renderPlot({
    r <- risultato()
    a <- r$attrib[, .(v = sum(attribuibili)), by = .(date, curve)]
    o <- r$obs[!is.na(eccesso), .(v = sum(eccesso)), by = date]
    p <- ggplot(a, aes(date, v, colour = curve)) +
      geom_hline(yintercept = 0, colour = "grey70")
    if (nrow(o)) p <- p + geom_col(data = o, aes(date, v), inherit.aes = FALSE,
                                   fill = "grey75", width = 1)
    p + geom_line(linewidth = .9) +
      scale_colour_manual(values = palette_curve(unique(a$curve))) +
      labs(title = sprintf("Decessi giornalieri attribuibili al caldo \u2014 %s", r$titolo),
           subtitle = if (nrow(o)) "Barre grigie: eccesso osservato (Istat, et\u00e0 20+)" else NULL,
           x = NULL, y = "Decessi al giorno") + tema()
  })

  output$g_cum <- renderPlot({
    r <- risultato()
    a <- r$attrib[, .(v = sum(attribuibili)), by = .(date, curve)][order(date)]
    a[, cum := cumsum(v), by = curve]
    o <- r$obs[!is.na(eccesso), .(v = sum(eccesso)), by = date][order(date)]
    p <- ggplot(a, aes(date, cum, colour = curve)) + geom_hline(yintercept = 0, colour = "grey70")
    if (nrow(o)) {
      o[, cum := cumsum(v)]
      p <- p + geom_line(data = o, aes(date, cum), inherit.aes = FALSE,
                         colour = COL[["Osservati"]], linewidth = 1, linetype = "22")
    }
    p + geom_line(linewidth = 1) +
      scale_colour_manual(values = palette_curve(unique(a$curve))) +
      labs(title = sprintf("Decessi cumulati attribuibili al caldo \u2014 %s", r$titolo),
           subtitle = if (nrow(o)) "Linea nera tratteggiata: eccesso osservato cumulato" else NULL,
           x = NULL, y = "Decessi cumulati") + tema()
  })

  output$g_temp <- renderPlot({
    r <- risultato()
    t <- merge(r$temp, CITIES[, .(code = URAU_CODE, name)], by = "code")
    ggplot(t, aes(date, tmean, group = name, colour = name)) + geom_line(linewidth = .7) +
      labs(title = "Temperatura media giornaliera (ERA5-Land / previsione)",
           x = NULL, y = "\u00b0C") + tema() +
      theme(legend.position = if (uniqueN(t$name) > 12) "none" else "bottom")
  })

  output$g_mappa <- renderPlot({
    r <- risultato()
    tot <- r$attrib[curve == unique(curve)[1], .(v = sum(attribuibili)), by = code]
    pts <- merge(CITIES, tot, by.x = "URAU_CODE", by.y = "code")
    reg <- BUNDLE$regions
    sel <- CITIES[URAU_CODE %in% r$codes]
    bb <- sf::st_bbox(sf::st_as_sf(sel, coords = c("lon", "lat"), crs = 4326))
    m <- ggplot() +
      geom_sf(data = reg, fill = "grey97", colour = "grey80", linewidth = .2)
    if (!is.null(BUNDLE$city_geom))
      m <- m + geom_sf(data = BUNDLE$city_geom[BUNDLE$city_geom$URAU_CODE %in% r$codes, ],
                       fill = "#B2182B", colour = NA, alpha = .35)
    m + geom_point(data = pts, aes(lon, lat, size = v), colour = "#B2182B", alpha = .85) +
      geom_text(data = pts, aes(lon, lat, label = name), size = 3, vjust = -1.1) +
      coord_sf(xlim = c(bb["xmin"] - .6, bb["xmax"] + .6),
               ylim = c(bb["ymin"] - .4, bb["ymax"] + .4)) +
      scale_size_area(max_size = 12) +
      labs(title = sprintf("Decessi attribuibili al caldo nel periodo \u2014 %s", r$titolo),
           size = "Decessi", x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_line(colour = "grey95"), legend.position = "right")
  })

  output$tab <- renderTable({
    r <- risultato()
    a <- r$attrib[, .(`Decessi attribuibili` = round(sum(attribuibili), 1)), by = .(Curva = curve)]
    o <- r$obs[!is.na(eccesso), .(Curva = "Eccesso osservato (Istat)",
                                  `Decessi attribuibili` = round(sum(eccesso), 1))]
    rbind(a, o)
  })

  output$copertura <- renderText({
    r <- tryCatch(risultato(), error = function(e) NULL)
    if (is.null(r)) return("")
    n <- r$obs[!is.na(eccesso), uniqueN(date)]
    paste0(sprintf("Dati Istat disponibili fino al %s (%d giorni del periodo scelto).",
                   format(BUNDLE$data_end, "%d/%m/%Y"), n),
           if (identical(input$modo, "fc"))
             sprintf(" Previsioni: modello %s%s.", BUNDLE$forecast_model,
                     if (isTRUE(input$debias)) ", con correzione della distorsione" else
                       ", senza correzione") else "")
  })

  output$msg <- renderUI({
    r <- tryCatch(risultato(), error = function(e) e)
    if (inherits(r, "error")) div(class = "errore", conditionMessage(r)) else NULL
  })
}

shinyApp(ui, server)
