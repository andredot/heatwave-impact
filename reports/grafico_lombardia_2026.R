# =============================================================================
#  Grafico: mortalita' osservata e attesa nel 2026, citta' lombarde del campione
#  Medie mobili a 7 giorni; l'area colorata e' l'eccesso (o il difetto).
#  Da eseguire dalla cartella del progetto, dopo targets::tar_make().
# =============================================================================
library(targets); library(data.table); library(ggplot2)

# ---- province lombarde (codici Istat) ---------------------------------------
PROV_LOMBARDIA <- c("012", "013", "014", "015", "016", "017", "018", "019",
                    "020", "097", "098", "108")

xw     <- tar_read(crosswalk)
citta  <- tar_read(cities)
daily  <- tar_read(counterfactual)$daily

codici <- unique(xw$map[substr(PRO_COM, 1, 3) %in% PROV_LOMBARDIA, URAU_CODE])
nomi   <- sort(citta[URAU_CODE %in% codici, name])
stopifnot(length(codici) > 0)

# ---- periodo mostrato -------------------------------------------------------
da <- as.Date("2026-02-01"); a <- as.Date("2026-12-31")   # es. "2026-05-01" per la sola estate

# ---- serie giornaliera aggregata (eta' 20+) ---------------------------------
# Attesi = livello stagionale con l'effetto del freddo, senza quello del caldo:
# cosi' l'area colorata nei mesi caldi e' l'eccesso attribuibile al caldo, e non
# si conta come "eccesso" la normale mortalita' invernale.
d <- daily[URAU_CODE %in% codici & date >= da & date <= a,
           .(osservati = sum(deaths),
             attesi = sum(B * fifelse(tmean < mmt, exp(logrr), 1))),
           by = date][order(date)]
d[, `:=`(oss7 = frollmean(osservati, 7, align = "center"),
         att7 = frollmean(attesi,    7, align = "center"))]
d <- d[!is.na(oss7)]
d[, eccesso := oss7 > att7]

ecc_tot <- round(sum(d$osservati - d$attesi))
wrap <- function(x, w = 135) paste(strwrap(x, width = w), collapse = "\n")
periodo <- sprintf("%s - %s", format(min(d$date), "%d/%m"),
                   format(max(d$date), "%d/%m/%Y"))
mesi <- c("gen", "feb", "mar", "apr", "mag", "giu",
          "lug", "ago", "set", "ott", "nov", "dic")

# ---- grafico ----------------------------------------------------------------
col_oss <- "#1F2933"; col_att <- "#2C7FB8"
col_pos <- "#C0392B"; col_neg <- "#3A8DAE"

p <- ggplot(d, aes(date)) +
  geom_ribbon(aes(ymin = att7, ymax = pmax(oss7, att7)), fill = col_pos, alpha = .28) +
  geom_ribbon(aes(ymin = pmin(oss7, att7), ymax = att7), fill = col_neg, alpha = .18) +
  geom_line(aes(y = att7, colour = "Attesi senza effetto del caldo"),
            linewidth = .8, linetype = "22") +
  geom_line(aes(y = oss7, colour = "Osservati"), linewidth = 1) +
  scale_colour_manual(values = c("Osservati" = col_oss,
                                 "Attesi senza effetto del caldo" = col_att),
                      breaks = c("Osservati", "Attesi senza effetto del caldo")) +
  scale_x_date(date_breaks = "1 month", labels = function(x) mesi[as.integer(format(x, "%m"))],
               expand = expansion(mult = c(.01, .01))) +
  scale_y_continuous(expand = expansion(mult = c(.05, .12))) +
  labs(
    title    = "Mortalit\u00e0 giornaliera nelle citt\u00e0 lombarde, 2026",
    subtitle = sprintf(
      "Medie mobili a 7 giorni. In rosso i giorni sopra l'attesa, in azzurro quelli sotto (saldo %s, %s)",
      format(ecc_tot, big.mark = ".", decimal.mark = ","), periodo),
    x = NULL, y = "Decessi al giorno (et\u00e0 20+)", colour = NULL,
    caption = paste0(
      wrap(paste0("Citt\u00e0 incluse: ", paste(nomi, collapse = ", "), ".")), "\n",
      wrap(paste0("Decessi giornalieri Istat per comune di residenza (et\u00e0 20+). ",
                  "Attesi: modello stagionale quasi-Poisson stimato sui giorni senza caldo, ",
                  "con l'effetto del freddo reintrodotto.")))) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = rel(1.25)),
    plot.subtitle   = element_text(colour = "grey30", margin = margin(b = 12)),
    plot.caption    = element_text(colour = "grey45", size = rel(.72), hjust = 0),
    plot.title.position = "plot", plot.caption.position = "plot",
    legend.position = "top", legend.justification = "left",
    legend.key.width = unit(28, "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92"),
    panel.grid.major.y = element_line(colour = "grey92"),
    axis.title.y    = element_text(colour = "grey30", margin = margin(r = 8)))

print(p)
ggsave("reports//mortalita_lombardia_2026.png", p, width = 10, height = 5.2, dpi = 200, bg = "white")
