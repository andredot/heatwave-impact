#' Masselot age groups
AGE_GROUPS <- c("20-44", "45-64", "65-74", "75-84", "85+")

#' Map Istat age classes to Masselot age groups
#'
#' Istat `CL_ETA`: 0 = 0 years, 1 = 1-4, 2 = 5-9, ..., 20 = 95-99, 21 = 100+.
#' Classes below 20 years return `NA` (excluded, as in Masselot et al.).
#'
#' @param cl Integer vector of `CL_ETA` codes.
#' @return Character vector of age groups.
istat_age_group <- function(cl) {
  cl <- as.integer(cl)
  out <- rep(NA_character_, length(cl))
  out[cl %in% 5:9]   <- "20-44"
  out[cl %in% 10:13] <- "45-64"
  out[cl %in% 14:15] <- "65-74"
  out[cl %in% 16:17] <- "75-84"
  out[cl %in% 18:21] <- "85+"
  out
}

#' Basis specification of a city's ERF
#'
#' Quadratic B-spline with knots at the 10th/75th/90th percentiles and boundary
#' knots at the minimum/maximum of the city's 1990-2019 ERA5-Land
#' distribution, as in the first stage of Masselot et al. The MMT is searched
#' over the percentiles in `mmp_range` (25th-99th in the paper).
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param cfg Configuration list.
#' @return List with `knots`, `bound`, `mmt_grid`, `p99`, `pct` (function
#'   returning any percentile).
erf_spec <- function(erf, code, cfg) {
  row <- erf$tdist[URAU_CODE == code]
  if (nrow(row) != 1) stop("No temperature distribution for ", code)
  pc <- grep("%$", names(row), value = TRUE)
  pn <- as.numeric(sub("%", "", pc))
  pv <- as.numeric(unlist(row[1, pc, with = FALSE]))
  inr <- pn >= cfg$mmp_range[1] & pn <= cfg$mmp_range[2]
  # Masselot's prediction grid: the percentiles strictly inside 0-100
  gr <- pv[pn > 0 & pn < 100]; gp <- pn[pn > 0 & pn < 100]
  list(knots = pct_value(row, cfg$var_pct),
       pred_grid = gr[order(gp)], pred_pct = sort(gp),
       bound = pct_value(row, c(0, 100)),
       degree = cfg$var_degree,
       mmt_grid = sort(unique(pv[inr])), mmp_range = cfg$mmp_range,
       p99 = pct_value(row, cfg$rr_pct),
       pct = function(p) pct_value(row, p))
}

#' B-spline basis of temperature for a city's ERF
#'
#' Values outside the boundary knots are extrapolated polynomially (as when
#' the published functions are applied to new extremes).
#'
#' @param t Temperatures.
#' @param spec Output of [erf_spec()].
#' @return Basis matrix (`length(t) x 5`).
erf_basis <- function(t, spec) {
  b <- suppressWarnings(splines::bs(t, knots = spec$knots, degree = spec$degree,
                                    Boundary.knots = spec$bound))
  unclass(b)[, , drop = FALSE]
}

#' Published coefficients of one city and age group
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param age Age group.
#' @return Numeric vector of length 5.
erf_coef <- function(erf, code, age) {
  r <- erf$coefs[URAU_CODE == code & agegroup == age]
  as.numeric(unlist(r[1, paste0("b", 1:5), with = FALSE]))
}

#' Published covariance matrix of one city and age group
#'
#' @inheritParams erf_coef
#' @return 5 x 5 matrix.
erf_vcov <- function(erf, code, age) {
  r <- erf$vcov[URAU_CODE == code & agegroup == age]
  vn <- grep("^v[0-9][0-9]$", names(r), value = TRUE)
  V <- matrix(0, 5, 5)
  for (v in vn) {
    i <- as.integer(substr(v, 2, 2)); j <- as.integer(substr(v, 3, 3))
    V[i, j] <- V[j, i] <- as.numeric(r[[v]][1])
  }
  V
}

#' Draws of published coefficients
#'
#' Uses the 1000 published simulations when available (they preserve the
#' correlation between age groups of the same city, since all ages are
#' projected from the same meta-regression draw), otherwise multivariate
#' normal draws from `vcov.csv`.
#'
#' @inheritParams erf_coef
#' @param nsim Number of draws when simulating from `vcov.csv`.
#' @return `nsim x 5` matrix.
erf_draws <- function(erf, code, age, nsim) {
  if (!is.null(erf$simu)) {
    s <- erf$simu[URAU_CODE == code & agegroup == age]
    if (nrow(s)) return(as.matrix(s[order(sim), paste0("b", 1:5), with = FALSE]))
  }
  rmvn(nsim, erf_coef(erf, code, age), erf_vcov(erf, code, age))
}

#' Minimum mortality temperature
#'
#' Computed exactly as in Masselot's `09_ResultsCityAge.R`: the curve is
#' evaluated on the grid of temperature percentiles (0.1th-99.9th) in a basis
#' whose boundary knots are the ends of that grid, and the minimum is taken
#' between the `mmp_range` percentiles. This differs from the basis used for
#' daily attribution (boundary knots at the minimum and maximum of the daily
#' series), which is also what their attribution code does.
#'
#' @param spec Output of [erf_spec()].
#' @param beta Coefficients.
#' @return The temperature with the lowest log-risk.
erf_mmt <- function(spec, beta) {
  g <- spec$pred_grid
  b <- suppressWarnings(splines::bs(g, knots = spec$knots, degree = spec$degree,
                                    Boundary.knots = range(g)))
  inr <- spec$pred_pct >= spec$mmp_range[1] & spec$pred_pct <= spec$mmp_range[2]
  g[inr][which.min((unclass(b) %*% beta)[inr])]
}

#' Build a curve from single-basis coefficients
#'
#' A "curve" is a small object that can be evaluated at any temperature,
#' returning the (uncentred) log-risk for the point estimate or for each draw.
#'
#' @param spec Output of [erf_spec()].
#' @param beta Point coefficients.
#' @param draws `nsim x 5` matrix of coefficient draws.
#' @return A curve object (list of class `erf_curve`).
curve_single <- function(spec, beta, draws) {
  structure(list(spec = spec, parts = list(list(beta = beta, draws = draws,
                                                 mmt = erf_mmt(spec, beta))),
                 w = 1), class = "erf_curve")
}

#' Build the published curve of a (possibly merged) age group
#'
#' For several age groups the risk of the merged series is the death-weighted
#' mean of the age-specific relative risks, each centred on its own MMT:
#' `f(t) = log sum_g w_g exp(f_g(t) - f_g(mmt_g))`. MMTs are held at their
#' point estimates in the draws, as in Masselot et al.
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param ages Age groups making up the series.
#' @param w Weights (e.g. deaths by age group); normalised internally.
#' @param spec Output of [erf_spec()].
#' @param nsim Number of draws.
#' @return A curve object.
curve_published <- function(erf, code, ages, w, spec, nsim) {
  parts <- lapply(ages, function(a) {
    b <- erf_coef(erf, code, a)
    list(beta = b, draws = erf_draws(erf, code, a, nsim), mmt = erf_mmt(spec, b))
  })
  structure(list(spec = spec, parts = parts, w = w / sum(w)), class = "erf_curve")
}

#' Evaluate a curve
#'
#' @param curve An `erf_curve`.
#' @param t Temperatures.
#' @param draws Return draws (`nsim x length(t)`) instead of the point estimate?
#' @return Numeric vector or matrix of log-risks (uncentred).
curve_eval <- function(curve, t, draws = FALSE) {
  B <- erf_basis(t, curve$spec)
  if (length(curve$parts) == 1) {
    p <- curve$parts[[1]]
    return(if (draws) p$draws %*% t(B) else as.numeric(B %*% p$beta))
  }
  acc <- 0
  for (k in seq_along(curve$parts)) {
    p <- curve$parts[[k]]
    Bm <- erf_basis(p$mmt, curve$spec)
    if (draws) {
      lr <- p$draws %*% t(B) - as.numeric(p$draws %*% t(Bm))
    } else {
      lr <- as.numeric(B %*% p$beta) - as.numeric(Bm %*% p$beta)
    }
    acc <- acc + curve$w[k] * exp(lr)
  }
  log(acc)
}

#' Summary metrics of a curve over a period
#'
#' Computes the MMT, the log relative risk at the `rr_pct` percentile
#' (versus a common reference temperature and versus the curve's own MMT) and
#' the heat-attributable fraction of the observed deaths with Masselot's
#' formula `AF = 1 - exp(-(f(T) - f(MMT)))` on days with `T >= MMT`.
#'
#' @param curve An `erf_curve`.
#' @param temps Daily temperatures of the period.
#' @param deaths Daily deaths of the period.
#' @param ref Common reference temperature (published MMT).
#' @return List with point values and draw vectors.
curve_metrics <- function(curve, temps, deaths, ref) {
  sp <- curve$spec
  mmt <- sp$mmt_grid[which.min(curve_eval(curve, sp$mmt_grid))]
  pts <- c(sp$p99, ref, mmt)
  fp <- curve_eval(curve, pts); fd <- curve_eval(curve, pts, draws = TRUE)
  ft <- curve_eval(curve, temps); ftd <- curve_eval(curve, temps, draws = TRUE)
  heat <- temps >= mmt
  an  <- sum((deaths * (1 - exp(-(ft - fp[3]))))[heat])
  and <- rowSums(sweep(1 - exp(-(ftd[, heat, drop = FALSE] - fd[, 3])), 2,
                       deaths[heat], "*"))
  list(mmt = mmt,
       lr99_ref = fp[1] - fp[2], lr99_ref_d = fd[, 1] - fd[, 2],
       lr99_own = fp[1] - fp[3], lr99_own_d = fd[, 1] - fd[, 3],
       af_heat = an / sum(deaths), af_heat_d = and / sum(deaths))
}

#' Curve values on a grid for plotting
#'
#' @param curve An `erf_curve`.
#' @param grid Temperatures.
#' @return data.table with `t`, `rr`, `lo`, `hi` (relative to the curve's MMT).
curve_band <- function(curve, grid) {
  sp <- curve$spec
  mmt <- sp$mmt_grid[which.min(curve_eval(curve, sp$mmt_grid))]
  f  <- curve_eval(curve, grid) - curve_eval(curve, mmt)
  fd <- curve_eval(curve, grid, draws = TRUE) -
    as.numeric(curve_eval(curve, mmt, draws = TRUE))
  data.table::data.table(t = grid, rr = exp(f),
                         lo = exp(apply(fd, 2, stats::quantile, 0.025)),
                         hi = exp(apply(fd, 2, stats::quantile, 0.975)))
}

#' Daily published log relative risk of one age group
#'
#' Point estimate, centred on the age group's MMT (the quantity used by the
#' app).
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param age Age group.
#' @param temps Temperatures.
#' @param cfg Configuration list.
#' @return List with `logrr` (vector), `mmt`, `spec`, `beta`.
published_logrr <- function(erf, code, age, temps, cfg) {
  sp <- erf_spec(erf, code, cfg)
  b <- erf_coef(erf, code, age)
  m <- erf_mmt(sp, b)
  list(logrr = as.numeric((erf_basis(temps, sp) - rep(erf_basis(m, sp), each = length(temps))) %*% b),
       mmt = m, spec = sp, beta = b)
}

#' Check that the published curves are reconstructed correctly
#'
#' Applies our reconstruction to Masselot's own exposure series (1990-2019)
#' and his own annual deaths, and compares the result with the numbers he
#' published in `results.zip`: minimum mortality temperature and percentile,
#' relative risk at the 99th percentile, and annual heat-attributable deaths
#' (his formula: AF = 1 - exp(-(f(T) - f(MMT))), summed over days above the
#' MMT, times annual deaths, divided by the number of days). Agreement means
#' any mismatch with Istat data comes from the curves, not from how we apply
#' them.
#'
#' @param erf Output of [read_erf_bundle()].
#' @param era5 Masselot's ERA5-Land series ([read_era5_masselot()]).
#' @param published Output of [read_published_results()].
#' @param cfg Configuration list.
#' @return data.table comparing our values with the published ones.
check_reconstruction <- function(erf, era5, published, cfg) {
  num <- function(x) suppressWarnings(as.numeric(x))
  col <- function(d, ...) { n <- intersect(c(...), names(d)); if (length(n)) n[1] else NA_character_ }
  c_mmt <- col(published, "mmt"); c_mmp <- col(published, "mmp")
  c_rr <- col(published, "rr_heat", "rrheat")
  c_an <- col(published, "excess_heat_est", "an_heat_est")
  c_death <- col(published, "death", "deaths")
  out <- list()
  for (i in seq_len(nrow(published))) {
    cd <- published$URAU_CODE[i]; ag <- published$agegroup[i]
    t <- era5[URAU_CODE == cd, tmean]
    if (!length(t)) next
    sp <- erf_spec(erf, cd, cfg)
    b <- erf_coef(erf, cd, ag)
    if (!length(b) || anyNA(b)) next
    mmt <- erf_mmt(sp, b)
    f <- function(x) as.numeric(erf_basis(x, sp) %*% b)
    af <- 1 - exp(-(f(t) - f(mmt)))
    heat <- t >= mmt
    death <- num(published[[c_death]][i])
    out[[length(out) + 1]] <- data.table::data.table(
      URAU_CODE = cd, agegroup = ag,
      mmt_ours = mmt, mmt_pub = num(published[[c_mmt]][i]),
      mmp_ours = stats::approx(sp$pred_grid, sp$pred_pct, mmt, rule = 2)$y,
      mmp_pub = if (!is.na(c_mmp)) num(published[[c_mmp]][i]) else NA_real_,
      rr99_ours = exp(f(sp$p99) - f(mmt)), rr99_pub = num(published[[c_rr]][i]),
      heat_ours = sum(af[heat]) * death / length(t),
      heat_pub = num(published[[c_an]][i]))
  }
  r <- data.table::rbindlist(out)
  r[, `:=`(rr_ratio = (rr99_ours - 1) / (rr99_pub - 1), heat_ratio = heat_ours / heat_pub,
           mmt_diff = mmt_ours - mmt_pub)]
  r[]
}

#' Baseline daily deaths used by the app
#'
#' Derived from the population, age structure and annual death rates in the
#' Masselot metadata (2000-2020 Eurostat figures), as in the Shiny app.
#'
#' @param erf Output of [read_erf_bundle()].
#' @param code Urban Audit code.
#' @param age Age group.
#' @return Expected deaths per day.
app_baseline <- function(erf, code, age) {
  bands <- list("20-44" = c("2024", "2529", "3034", "3539", "4044"),
                "45-64" = c("4549", "5054", "5559", "6064"),
                "65-74" = c("6569", "7074"),
                "75-84" = c("7579", "8084"),
                "85+"   = c("8599"))
  m <- erf$metadata[URAU_CODE == code]
  ann <- 0
  for (b in bands[[age]]) {
    pr <- suppressWarnings(as.numeric(m[[paste0("prop_", b)]]))
    dr <- suppressWarnings(as.numeric(m[[paste0("deathrate_", b)]]))
    if (length(pr) && length(dr) && !is.na(pr) && !is.na(dr))
      ann <- ann + as.numeric(m$pop) * pr / 100 * dr
  }
  ann / 365.25
}
