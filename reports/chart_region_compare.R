# =============================================================================
#  One chart, one set of cities, two baselines.
#
#  Supersedes the two earlier scripts. It draws, for the validated cities of
#  the configured region (ages 20+):
#    - observed daily deaths (Istat);
#    - deaths expected without heat from this project's counterfactual model
#      (daily, by city and age group, cold effect kept in);
#    - deaths expected without extreme temperature from the official FluMOMO
#      code run on the same comuni (weekly, influenza kept in via EB + EdIA).
#  The lower panel shows the excess accumulated against each baseline.
#
#  Needs: targets::tar_make() for `counterfactual`, `crosswalk`, `cities` and
#  `flumomo_cities`. Run from the project root.
# =============================================================================
library(targets); library(data.table); library(ggplot2)
source("R//config.R"); for (f in list.files("R", full.names = TRUE, pattern = "[.]R$")) source(f)

YEAR <- cfg$flumomo$chart_year
OUT  <- sprintf("confronto_baseline_%s_%d.png", tolower(cfg$flumomo$region), YEAR)

cities <- as.data.table(tar_read(cities))
xw     <- tar_read(crosswalk)$map
cf     <- tar_read(counterfactual)
fmres  <- tryCatch(tar_read(flumomo_cities), error = function(e)
  stop("Run targets::tar_make(flumomo_cities) first (official FluMOMO code needed)."))

# cities of the region that the validation covers
codes <- unique(xw[substr(PRO_COM, 1, 3) %in% cfg$flumomo$provinces$prov, URAU_CODE])
names_it <- sort(cities[URAU_CODE %in% codes, name])

# ---- observed and our baseline: daily, ages 20+ ------------------------------
d <- cf$daily[URAU_CODE %in% codes & year(date) == YEAR,
              .(observed = sum(deaths),
                # our expected: baseline at the MMT with the cold side put back
                ours = sum(B * fifelse(tmean < mmt, exp(logrr), 1))),
              by = date][order(date)]

# ---- FluMOMO baseline: weekly, ages 15+, rescaled to the same population -----
r <- as.data.frame(fmres)
r <- r[r$agegrp %in% 2:3 & !is.na(r$year) & !is.na(r$week), ]
wk <- as.data.table(r)[, .(expected = sum(EB + EdIA), deaths = sum(deaths)),
                       by = .(year, week)]
wk[, week_start := ISOweek::ISOweek2date(sprintf("%d-W%02d-1", year, week))]
setorder(wk, week_start)
# 15-19 year olds are in the FluMOMO groups but not in ours: rescale by the
# ratio of the two observed totals over the plotted year
obs20 <- sum(d$observed)
obs15 <- sum(wk[year == YEAR & week_start <= max(d$date), deaths])
scale15 <- if (obs15 > 0) obs20 / obs15 else 1
d[, flumomo := scale15 * stats::splinefun(as.numeric(wk$week_start) + 3,
                                          wk$expected / 7)(as.numeric(date))]

# ---- 7-day moving averages and cumulative excess -----------------------------
d[, `:=`(obs7 = frollmean(observed, 7, align = "center"),
         ours7 = frollmean(ours, 7, align = "center"),
         flu7 = frollmean(flumomo, 7, align = "center"))]
d <- d[!is.na(obs7)]
d[, `:=`(cum_ours = cumsum(observed - ours), cum_flu = cumsum(observed - flumomo))]

lab_ours <- "Attesi senza caldo (modello del progetto)"
lab_flu  <- "Attesi senza temperature estreme (FluMOMO)"
daily_long <- melt(d[, .(date, Osservati = obs7, a = ours7, b = flu7)],
                   id.vars = "date", variable.name = "serie", value.name = "v")
daily_long[, serie := factor(serie, c("Osservati", "a", "b"),
                             c("Osservati", lab_ours, lab_flu))]
cum_long <- melt(d[, .(date, a = cum_ours, b = cum_flu)], id.vars = "date",
                 variable.name = "serie", value.name = "v")
cum_long[, serie := factor(serie, c("a", "b"), c(lab_ours, lab_flu))]
daily_long[, pannello := "Decessi al giorno (media mobile 7 giorni)"]
cum_long[, pannello := "Eccesso cumulato dall'inizio dell'anno"]
plot_dt <- rbind(daily_long, cum_long)

COLS <- setNames(c("#1F2933", "#2C7FB8", "#D95F02"), c("Osservati", lab_ours, lab_flu))
mesi <- c("gen", "feb", "mar", "apr", "mag", "giu",
          "lug", "ago", "set", "ott", "nov", "dic")
wrap <- function(x, n = 135) paste(strwrap(x, width = n), collapse = "\n")
period <- sprintf("%s - %s", format(min(d$date), "%d/%m"), format(max(d$date), "%d/%m/%Y"))

p <- ggplot(plot_dt, aes(date, v, colour = serie, linetype = serie)) +
  geom_hline(data = data.frame(pannello = "Eccesso cumulato dall'inizio dell'anno", y = 0),
             aes(yintercept = y), colour = "grey70", inherit.aes = FALSE) +
  geom_line(linewidth = .9) +
  facet_wrap(~ pannello, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = COLS) +
  scale_linetype_manual(values = setNames(c("solid", "22", "22"), names(COLS))) +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) mesi[as.integer(format(x, "%m"))],
               expand = expansion(mult = c(.01, .01))) +
  labs(
    title = sprintf("Mortalit\u00e0 e attesa a confronto: citt\u00e0 della %s, %d",
                    cfg$flumomo$region, YEAR),
    subtitle = wrap(sprintf(paste("Stesse citt\u00e0, due attese indipendenti (et\u00e0 20+, %s).",
                                  "Saldo a fine periodo: %s con il modello del progetto,",
                                  "%s con FluMOMO"), period,
                            format(round(tail(d$cum_ours, 1)), big.mark = ".",
                                   decimal.mark = ","),
                            format(round(tail(d$cum_flu, 1)), big.mark = ".",
                                   decimal.mark = ",")), 110),
    x = NULL, y = NULL, colour = NULL, linetype = NULL,
    caption = wrap(paste0(
      "Citt\u00e0: ", paste(names_it, collapse = ", "), ". ",
      "Decessi giornalieri Istat per comune di residenza. Attesi del progetto: modello ",
      "stagionale quasi-Poisson stimato sui giorni senza caldo, con l'effetto del freddo ",
      "reintrodotto. Attesi FluMOMO: codice ufficiale EuroMOMO 4.2 sugli stessi comuni ",
      "(baseline + quota attribuita all'influenza, EB + EdIA), settimanale e per fasce ",
      "15-64 e 65+, riportato alla popolazione 20+ e ripartito sui giorni. Le due attese ",
      "differiscono anche d'inverno perch\u00e9 FluMOMO esclude pure il freddo estremo."))) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = rel(1.15)),
        plot.subtitle = element_text(colour = "grey30", margin = margin(b = 10)),
        plot.caption = element_text(colour = "grey45", size = rel(.7), hjust = 0),
        plot.title.position = "plot", plot.caption.position = "plot",
        strip.text = element_text(hjust = 0, face = "bold", colour = "grey25"),
        legend.position = "top", legend.justification = "left",
        legend.key.width = unit(28, "pt"), panel.grid.minor = element_blank())

print(p)
ggsave(OUT, p, width = 11, height = 7, dpi = 200, bg = "white")
message("Written ", OUT)

print(d[, .(giorni = .N, osservati = sum(observed),
            attesi_progetto = round(sum(ours)), attesi_flumomo = round(sum(flumomo)),
            eccesso_progetto = round(sum(observed - ours)),
            eccesso_flumomo = round(sum(observed - flumomo)))])

