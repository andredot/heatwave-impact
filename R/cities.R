#' Download the GISCO boundary layers
#'
#' Urban Audit polygons for every level listed in `cfg$urau_levels`
#' (cities, greater cities, functional urban areas) and LAU (comuni)
#' polygons, via `giscoR` when installed, otherwise from the GISCO
#' distribution API. Files are cached as GeoPackages.
#'
#' @param cfg Configuration list.
#' @return Character vector of file paths (`urau_<level>_...` and `lau_...`).
download_gisco <- function(cfg) {
  dir <- ensure_dir(file.path(cfg$data_dir, "gisco"))
  base <- "https://gisco-services.ec.europa.eu/distribution/v2"
  out <- character(0)
  for (yr in cfg$urau_years) for (lvl in cfg$urau_levels) {
    f <- file.path(dir, sprintf("urau_%s_%s_%s.gpkg", lvl, yr, cfg$country))
    if (!file.exists(f)) {
      x <- tryCatch({
        if (requireNamespace("giscoR", quietly = TRUE)) {
          giscoR::gisco_get_urban_audit(year = yr, epsg = 4326,
                                        level = lvl, country = cfg$country)
        } else {
          y <- sf::st_read(sprintf("%s/urau/geojson/URAU_RG_100K_%s_4326_%s.geojson",
                                   base, yr, lvl), quiet = TRUE)
          y[y$CNTR_CODE == cfg$country, ]
        }
      }, error = function(e) e)
      # Not every level exists in every release: from Urban Audit 2021 greater
      # cities are merged into the city level, so GREATER_CITIES is absent.
      if (inherits(x, "error") || is.null(x) || !nrow(x)) {
        message("GISCO level '", lvl, "' not available for ", yr, ", skipped.")
        next
      }
      sf::st_write(x, f, quiet = TRUE, delete_dsn = TRUE)
    }
    out <- c(out, f)
  }
  if (!length(out))
    stop("No Urban Audit level could be downloaded for ", paste(cfg$urau_years, collapse = "/"), ".")
  f_lau <- file.path(dir, sprintf("lau_%s_%s.gpkg", cfg$lau_year, cfg$country))
  if (!file.exists(f_lau)) {
    x <- if (requireNamespace("giscoR", quietly = TRUE)) {
      giscoR::gisco_get_lau(year = cfg$lau_year, epsg = 4326, country = cfg$country)
    } else {
      x <- sf::st_read(sprintf("%s/lau/geojson/LAU_RG_01M_%s_4326.geojson",
                               base, cfg$lau_year), quiet = TRUE)
      x[x$CNTR_CODE == cfg$country, ]
    }
    sf::st_write(x, f_lau, quiet = TRUE, delete_dsn = TRUE)
  }
  c(out, f_lau)
}

#' Candidate city definitions for each Masselot city
#'
#' Urban Audit codes were renumbered between releases, so matching on codes
#' can silently pick another city (e.g. a 2020 polygon coded like Masselot's
#' Messina is Battipaglia). Matching is therefore geographic: each
#' Masselot city point (the Urban Audit point, `lat`/`lon` in the metadata)
#' is located in the polygons of every Urban Audit level (city, greater
#' city, FUA). Each polygon containing the point is a candidate definition,
#' described by the comuni with at least `cfg$lau_share_min` of their area
#' inside it and their population. [finalize_crosswalk()] picks the candidate
#' whose population matches Masselot's.
#'
#' @param gisco Output of [download_gisco()].
#' @param metadata Country metadata from [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `candidates` (one row per Masselot city x candidate
#'   polygon) and `members` (comuni of each candidate polygon).
overlay_cities <- function(gisco, metadata, cfg) {
  old <- suppressMessages(sf::sf_use_s2(FALSE)); on.exit(suppressMessages(sf::sf_use_s2(old)))
  lau <- sf::st_read(gisco[grepl("^lau_", basename(gisco))][1], quiet = TRUE)
  idcol  <- intersect(c("LAU_ID", "GISCO_ID"), names(lau))[1]
  popcol <- grep("^POP_", names(lau), value = TRUE)[1]
  namecol <- intersect(c("LAU_NAME", "LAU_LABEL"), names(lau))[1]
  lau$PRO_COM <- pad_procom(sub("^[A-Z]{2}_", "", lau[[idcol]]))
  lau$lau_pop <- if (is.na(popcol)) NA_real_ else as.numeric(lau[[popcol]])
  lau$LAU_NAME <- lau[[namecol]]
  lau <- sf::st_make_valid(sf::st_transform(lau[, c("PRO_COM", "LAU_NAME", "lau_pop")], 3035))
  lau$lau_area <- as.numeric(sf::st_area(lau))

  pts <- sf::st_transform(sf::st_as_sf(
    data.frame(URAU_CODE = metadata$URAU_CODE, lon = as.numeric(metadata$lon),
               lat = as.numeric(metadata$lat)),
    coords = c("lon", "lat"), crs = 4326), 3035)

  cand <- list(); memb <- list()
  for (f in gisco[grepl("^urau_", basename(gisco))]) {
    lvl <- sub("^urau_(.*)_([0-9]{4})_[A-Z]{2}\\.gpkg$", "\\1 \\2", basename(f))
    u <- sf::st_read(f, quiet = TRUE)
    names(u) <- sub("^URAU_ID$", "URAU_CODE", names(u))
    u <- sf::st_make_valid(sf::st_transform(u[, c("URAU_CODE", "URAU_NAME")], 3035))
    names(u)[1:2] <- c("poly_code", "poly_name")
    hit <- sf::st_drop_geometry(sf::st_join(pts, u, join = sf::st_within))
    hit <- data.table::as.data.table(hit)
    # A point can fall just outside its polygon (coastline generalisation):
    # fall back to the nearest polygon within `point_max_dist` metres.
    miss <- which(is.na(hit$poly_code))
    if (length(miss) && !is.null(cfg$point_max_dist)) {
      nr <- sf::st_nearest_feature(pts[miss, ], u)
      dd <- as.numeric(sf::st_distance(pts[miss, ], u[nr, ], by_element = TRUE))
      ok <- dd <= cfg$point_max_dist
      hit$poly_code[miss[ok]] <- u$poly_code[nr[ok]]
      hit$poly_name[miss[ok]] <- u$poly_name[nr[ok]]
    }
    hit <- hit[!is.na(poly_code)]
    if (!nrow(hit)) next
    hit[, level := lvl]
    cand[[length(cand) + 1]] <- hit
    polys <- u[u$poly_code %in% hit$poly_code, ]
    polys$poly_area <- as.numeric(sf::st_area(polys))
    ix <- suppressWarnings(sf::st_intersection(lau, polys))
    ix$area <- as.numeric(sf::st_area(ix))
    ix <- data.table::as.data.table(sf::st_drop_geometry(ix))
    ix[, `:=`(share_lau = area / lau_area, level = lvl)]
    memb[[length(memb) + 1]] <- ix[share_lau > 0.001,
                                   .(level, poly_code, PRO_COM, LAU_NAME, lau_pop, share_lau)]
  }
  cand <- data.table::rbindlist(cand); memb <- data.table::rbindlist(memb)
  if (!nrow(cand)) stop("No Masselot city point falls inside any Urban Audit polygon.")
  inside <- memb[share_lau >= cfg$lau_share_min,
                 .(n_comuni = .N, lau_pop = sum(lau_pop, na.rm = TRUE)), by = .(level, poly_code)]
  cand <- merge(cand, inside, by = c("level", "poly_code"), all.x = TRUE)
  list(candidates = cand, members = memb)
}

#' Choose each city's definition and run the alignment checks
#'
#' Among the candidate polygons of a city (see [overlay_cities()]), the one
#' whose comuni population is closest to Masselot's population (on the log
#' scale) is chosen, unless an override fixes the level. Overrides are read
#' from `cfg$crosswalk_overrides` (CSV, columns `URAU_CODE`, `level`,
#' `PRO_COM`, `action`): `action = "level"` forces a level, `"add"` and
#' `"remove"` edit the comuni list. The chosen definition is checked on:
#' \enumerate{
#'   \item comune codes present in the Istat file;
#'   \item polygon name consistent with Masselot's city name;
#'   \item population ratio;
#'   \item Istat deaths aged 20+ (2011-2019 mean) against the deaths implied
#'         by Masselot's population structure and death rates (decisive).
#' }
#'
#' @param overlay Output of [overlay_cities()].
#' @param istat Output of [read_istat()].
#' @param erf Output of [read_erf_bundle()].
#' @param cfg Configuration list.
#' @return List with `map` (city-comune table used downstream), `checks`
#'   (one row per city) and `candidates` (all candidate definitions).
finalize_crosswalk <- function(overlay, istat, erf, cfg) {
  meta <- erf$metadata
  cand <- merge(overlay$candidates, meta[, .(URAU_CODE, meta_name = URAU_NAME,
                                             pop = as.numeric(pop))], by = "URAU_CODE")
  cand[, pop_ratio := lau_pop / pop]
  ov <- if (!is.null(cfg$crosswalk_overrides) && file.exists(cfg$crosswalk_overrides))
    data.table::fread(cfg$crosswalk_overrides, colClasses = "character") else
      data.table::data.table(URAU_CODE = character(), level = character(),
                             PRO_COM = character(), action = character())
  forced <- ov[action == "level"]
  cand[, forced := paste(URAU_CODE, level) %in% paste(forced$URAU_CODE, forced$level)]
  cand[, score := data.table::fifelse(forced, -Inf, abs(log(pmax(pop_ratio, 1e-6))))]
  cand[, chosen := seq_len(.N) == which.min(score), by = URAU_CODE]
  pick <- cand[chosen == TRUE]

  map <- merge(pick[, .(URAU_CODE, level, poly_code)],
               overlay$members[share_lau >= cfg$lau_share_min],
               by = c("level", "poly_code"), allow.cartesian = TRUE)
  map <- map[, .(URAU_CODE, level, PRO_COM, LAU_NAME, lau_pop)]
  if (nrow(ov[action %in% c("add", "remove")])) {
    ov[, PRO_COM := pad_procom(PRO_COM)]
    drop <- ov[action == "remove"]
    if (nrow(drop)) map <- map[!drop, on = .(URAU_CODE, PRO_COM)]
    ad <- merge(ov[action == "add", .(URAU_CODE, PRO_COM)], pick[, .(URAU_CODE, level)],
                by = "URAU_CODE")[, `:=`(LAU_NAME = NA_character_, lau_pop = NA_real_)]
    map <- unique(rbind(map, ad, use.names = TRUE), by = c("URAU_CODE", "PRO_COM"))
  }
  map[, in_istat := PRO_COM %in% unique(istat$long$PRO_COM)]
  # A few comune codes change between the LAU release and the Istat file
  # (mergers). Drop them when they are a negligible share of the city.
  miss <- map[in_istat == FALSE, .(miss_pop = sum(lau_pop, na.rm = TRUE)), by = URAU_CODE]
  tot <- map[, .(tot_pop = sum(lau_pop, na.rm = TRUE)), by = URAU_CODE]
  miss <- merge(miss, tot, by = "URAU_CODE")[, share := miss_pop / tot_pop]
  minor <- miss[share <= cfg$missing_pop_tol, URAU_CODE]
  map <- map[in_istat == TRUE | !(URAU_CODE %in% minor)]

  yrs <- data.table::year(istat$long$date)
  ref <- istat$long[yrs %in% 2011:2019, .(deaths = sum(deaths, na.rm = TRUE)), by = PRO_COM]
  nyr <- length(unique(yrs[yrs %in% 2011:2019]))
  implied <- vapply(meta$URAU_CODE, function(cd)
    sum(vapply(AGE_GROUPS, function(a) app_baseline(erf, cd, a), 0)) * 365.25, 0)

  chk <- map[, .(n_comuni = .N, all_in_istat = all(in_istat),
                 lau_pop = sum(lau_pop, na.rm = TRUE)), by = .(URAU_CODE, level)]
  chk <- merge(meta[, .(URAU_CODE, name, meta_name = URAU_NAME, pop = as.numeric(pop), inmcc)],
               chk, by = "URAU_CODE", all.x = TRUE)
  chk <- merge(chk, pick[, .(URAU_CODE, poly_code, poly_name, n_candidates = NA_integer_)],
               by = "URAU_CODE", all.x = TRUE)
  chk[, n_candidates := vapply(URAU_CODE, function(cd) sum(cand$URAU_CODE == cd), 0L)]
  chk[, istat_deaths_yr := vapply(URAU_CODE, function(cd)
    sum(ref$deaths[ref$PRO_COM %in% map$PRO_COM[map$URAU_CODE == cd]]) / max(nyr, 1), 0)]
  chk[, implied_deaths_yr := implied[URAU_CODE]]
  chk[, `:=`(pop_ratio = lau_pop / pop, deaths_ratio = istat_deaths_yr / implied_deaths_yr)]
  chk[, name_ok := mapply(function(a, b) {
    a <- norm_name(a); b <- norm_name(b)
    !is.na(a) && !is.na(b) && (grepl(a, b, fixed = TRUE) || grepl(b, a, fixed = TRUE))
  }, poly_name, meta_name)]
  # Population is the direct check of the geography; the deaths ratio also
  # depends on the death rates in Masselot's metadata (regional averages), so
  # a deaths mismatch with a matching population is not a wrong city.
  pop_ok <- function(x) !is.na(x) & x >= cfg$pop_ratio_ok[1] & x <= cfg$pop_ratio_ok[2]
  dth_ok <- function(x) !is.na(x) & x >= cfg$deaths_ratio_ok[1] & x <= cfg$deaths_ratio_ok[2]
  chk[, status := data.table::fcase(
    is.na(n_comuni), "no polygon contains the city point",
    !all_in_istat, "codes missing in Istat",
    !pop_ok(pop_ratio) & !dth_ok(deaths_ratio), "wrong city definition",
    !pop_ok(pop_ratio), "check population",
    !dth_ok(deaths_ratio), "deaths differ from metadata rates",
    !name_ok, "check name",
    default = "ok")]
  chk[, analysed := !status %in% c("no polygon contains the city point",
                                   "codes missing in Istat", "wrong city definition")]
  list(map = map[URAU_CODE %in% chk[analysed == TRUE, URAU_CODE] & in_istat],
       checks = chk[order(-pop)], candidates = cand[order(URAU_CODE, level)])
}
