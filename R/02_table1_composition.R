################################################################################
## Table 1: RDS-II weighted sample composition
## Harmonised comparable sample (Maputo and Nampula, age >= 18, injection
## in the previous 12 months), serial BBS 2014 and 2023.
##
## Why this replaces the original Table 1 and the pooled svychisq version
## --------------------------------------------------------------------------
## Original manuscript: unweighted n (%) and Pearson/Fisher tests. Invalid
## for RDS (unequal inclusion probabilities; dependent recruitment).
##
## Previous revision Table 1: pooled the two cities with population-calibrated
## Gile weights and a single survey::svychisq. That is still incorrect because
##   1. City-weighted percentages for CIDADE are the PSE ratios, not data.
##   2. Gile SS was abandoned for prevalence after chain pruning; Table 1
##      must use the same RDS-II (Volz-Heckathorn) estimator.
##   3. Composition tests are city-stratified Wald tests of RDS-II shares,
##      not chain-clustered models. 2014 recruiter IDs are used in the Gile
##      SS Table 1 (PID_TABLE1_POOLED_SSPSE.R), not in these degree-based
##      RDS-II category shares.
##
## Correct approach (same as the prevalence analysis)
## --------------------------------------------------------------------------
##  * Restrict to the harmonised eligible sample.
##  * Stratify by city (independent RDS samples).
##  * Estimate category shares with RDS-II: sum(y/d) / sum(1/d).
##  * Bootstrap the estimator within city-year for SEs and 95% CIs.
##  * Compare years with a Wald test of independent estimates.
##  * Overall p for a K-category variable: multivariate Wald on K-1
##    categories using the bootstrap covariance.
##  * Two-city pool (supplement): population-weighted combination of the
##    city RDS-II estimates, not a sample-size mix. City itself is not tested.
##
## Data analyst: Rachid Muleia
## Date: 2026-09-19
################################################################################

options(stringsAsFactors = FALSE)
set.seed(20260919)

## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))


suppressPackageStartupMessages({
  library(tidyverse)
  library(xtable)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

N_BOOT <- 400
pop_n <- list(
  "2014" = c(MAPUTO = 1684, NAMPULA = 520),
  "2023" = c(MAPUTO = 990, NAMPULA = 1792)
)

fmt_pct <- function(x, d = 1) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", d, "f"), 100 * x))
}
fmt_pct_ci <- function(est, lo, hi, d = 1) {
  ifelse(
    is.na(est) | is.na(lo) | is.na(hi),
    NA_character_,
    sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), 100 * est, 100 * lo, 100 * hi)
  )
}
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
write_csv_safe <- function(x, file) {
  utils::write.csv(as.data.frame(x, stringsAsFactors = FALSE), file,
                   row.names = FALSE, na = "")
  invisible(x)
}
escape_tex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("<", "$<$", x, fixed = TRUE)
  x <- gsub(">=", "$\\\\ge$", x, fixed = TRUE)
  x
}

## RDS-II / Volz-Heckathorn category shares among non-missing values
vh_shares <- function(x, degree, levels) {
  ok <- !is.na(x) & !is.na(degree) & degree > 0
  x <- as.character(x[ok])
  d <- as.numeric(degree[ok])
  n <- length(x)
  if (n < 5) {
    return(list(
      n = setNames(rep(0L, length(levels)), levels),
      n_tot = n,
      p = setNames(rep(NA_real_, length(levels)), levels)
    ))
  }
  w <- 1 / d
  sw <- sum(w)
  n_lev <- vapply(levels, function(lv) as.integer(sum(x == lv)), integer(1))
  p_lev <- vapply(levels, function(lv) sum(w[x == lv]) / sw, numeric(1))
  list(n = n_lev, n_tot = n, p = p_lev)
}

wald_diff <- function(p1, se1, p0, se0) {
  d <- p1 - p0
  se <- sqrt(se1^2 + se0^2)
  z <- ifelse(se > 0, d / se, NA_real_)
  pval <- 2 * stats::pnorm(-abs(z))
  list(diff = d, se = se,
       lo = d - 1.96 * se, hi = d + 1.96 * se, p = pval)
}

## Joint bootstrap of all category shares for one city-year sample
boot_city_year <- function(dat, var_levels, n_boot = N_BOOT) {
  n <- nrow(dat)
  point <- lapply(names(var_levels), function(v) {
    vh_shares(dat[[v]], dat$degree, var_levels[[v]])
  })
  names(point) <- names(var_levels)

  boot <- lapply(names(var_levels), function(v) {
    matrix(NA_real_, n_boot, length(var_levels[[v]]),
           dimnames = list(NULL, var_levels[[v]]))
  })
  names(boot) <- names(var_levels)

  for (b in seq_len(n_boot)) {
    ii <- sample.int(n, n, replace = TRUE)
    db <- dat[ii, ]
    for (v in names(var_levels)) {
      boot[[v]][b, ] <- vh_shares(db[[v]], db$degree, var_levels[[v]])$p
    }
  }
  list(point = point, boot = boot)
}

se_from_boot <- function(mat) {
  apply(mat, 2, stats::sd, na.rm = TRUE)
}

## Multivariate Wald on K-1 categories; years treated as independent RDS samples
overall_wald <- function(p14, V14, p23, V23) {
  k <- length(p14)
  if (k < 2) return(NA_real_)
  idx <- seq_len(k - 1L)
  diff <- p23[idx] - p14[idx]
  V <- V14[idx, idx, drop = FALSE] + V23[idx, idx, drop = FALSE]
  V[!is.finite(V)] <- 0
  ok <- is.finite(diff)
  if (!any(ok)) return(NA_real_)
  diff <- diff[ok]
  V <- V[ok, ok, drop = FALSE]
  ## Ridge if nearly singular (rare categories)
  ev <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (any(!is.finite(ev)) || min(ev) < 1e-12) {
    V <- V + diag(1e-8, nrow(V))
  }
  W <- tryCatch(as.numeric(t(diff) %*% solve(V, diff)), error = function(e) NA_real_)
  if (!is.finite(W) || W < 0) return(NA_real_)
  stats::pchisq(W, df = length(diff), lower.tail = FALSE)
}

cov_from_boot <- function(mat) {
  stats::cov(mat, use = "pairwise.complete.obs")
}

## ---- Load comparable sample with RDS-II weights and degree ------------------

obj <- readRDS(file.path(path_out, "analysis_objects.rds"))
dat <- obj$dat_comp |>
  mutate(
    STI_SELF = case_when(
      STI_SELF %in% c("Yes", "1_SIM") ~ "Yes",
      STI_SELF %in% c("No", "2_NAO", "2_SIM") ~ "No",
      TRUE ~ NA_character_
    ),
    AGE_CAT = recode(AGE_CAT, "18-24" = "18-24", "25+" = "25+"),
    AGE_FIRST_DRUG = recode(AGE_FIRST_DRUG, "<18" = "<18", "18-24" = "18-24",
                            "25+" = "25+")
  )

stopifnot(all(c("YEAR", "CIDADE", "degree", "w_pop_ii") %in% names(dat)))

## Labels and category order. Main Table 1 is shortened (R1.14); the rest
## go to the supplement. TX_SEX and USED_NEEDLE are incomparable across
## questionnaires and are flagged rather than interpreted as change.
var_meta <- tibble::tribble(
  ~variable,          ~label,                                      ~panel,        ~comparable,
  "AGE_CAT",          "Age (years)",                               "main",        TRUE,
  "SEX_CAT",          "Sex",                                       "main",        TRUE,
  "EDUC_CAT",         "Education",                                 "main",        TRUE,
  "MARITALC",         "Marital status",                            "main",        TRUE,
  "RELIGIAO_CAT",     "Religion",                                  "main",        TRUE,
  "AGE_FIRST_DRUG",   "Age at first injection",                    "main",        TRUE,
  "FREQ_INJEC",       "Injection frequency",                       "main",        TRUE,
  "NEW_SYRINGE",      "New syringe access",                        "main",        TRUE,
  "UNPROTECTED_SEX",  "Unprotected sex",                           "main",        TRUE,
  "HIV_TESTED_12M",   "HIV test, past 12 months",                  "main",        TRUE,
  "SOUGHT_CARE",      "Sought health care",                        "main",        TRUE,
  "SHARED_EQUIP",     "Shared other injection equipment",          "supplement",  TRUE,
  "USED_NEEDLE",      "Injected with a used needle/syringe",       "supplement",  FALSE,
  "TX_SEX",           "Transactional sex",                         "supplement",  FALSE,
  "STI_SELF",         "Self-reported STI, past 12 months",         "supplement",  TRUE
)

var_levels <- list(
  AGE_CAT         = c("18-24", "25+"),
  SEX_CAT         = c("Male", "Female"),
  EDUC_CAT        = c("Primary or none", "Secondary or higher"),
  MARITALC        = c("Never married", "Married or union", "Other"),
  RELIGIAO_CAT    = c("Christian", "Muslim", "None/other"),
  AGE_FIRST_DRUG  = c("<18", "18-24", "25+"),
  FREQ_INJEC      = c("Daily", "Less than daily"),
  NEW_SYRINGE     = c("Yes", "No"),
  UNPROTECTED_SEX = c("Yes", "No"),
  HIV_TESTED_12M  = c("Yes", "No"),
  SOUGHT_CARE     = c("Yes", "No"),
  SHARED_EQUIP    = c("Yes", "No"),
  USED_NEEDLE     = c("Yes", "No"),
  TX_SEX          = c("Yes", "No"),
  STI_SELF        = c("Yes", "No")
)

cities <- c("MAPUTO", "NAMPULA")
years  <- c("2014", "2023")

message("Bootstrapping RDS-II composition by city and year (", N_BOOT, " replicates) ...")
city_boot <- list()
for (ct in cities) {
  city_boot[[ct]] <- list()
  for (yr in years) {
    message("  ", yr, " ", ct)
    sub <- dat |> filter(YEAR == yr, CIDADE == ct)
    city_boot[[ct]][[yr]] <- boot_city_year(sub, var_levels, n_boot = N_BOOT)
  }
}

## ---- City-stratified estimates ----------------------------------------------

rows <- list()
overall_p <- list()

for (ct in cities) {
  for (v in names(var_levels)) {
    levs <- var_levels[[v]]
    p14 <- city_boot[[ct]][["2014"]]$point[[v]]
    p23 <- city_boot[[ct]][["2023"]]$point[[v]]
    b14 <- city_boot[[ct]][["2014"]]$boot[[v]]
    b23 <- city_boot[[ct]][["2023"]]$boot[[v]]
    se14 <- se_from_boot(b14)
    se23 <- se_from_boot(b23)
    V14 <- cov_from_boot(b14)
    V23 <- cov_from_boot(b23)
    p_over <- overall_wald(p14$p, V14, p23$p, V23)
    overall_p[[paste(ct, v, sep = "|")]] <- p_over

    for (lv in levs) {
      wd <- wald_diff(p23$p[[lv]], se23[[lv]], p14$p[[lv]], se14[[lv]])
      rows[[length(rows) + 1]] <- tibble(
        city = ct,
        variable = v,
        category = lv,
        cat_n = match(lv, levs),
        n_2014 = unname(p14$n[[lv]]),
        n_tot_2014 = p14$n_tot,
        n_2023 = unname(p23$n[[lv]]),
        n_tot_2023 = p23$n_tot,
        rds_p_2014 = unname(p14$p[[lv]]),
        rds_se_2014 = unname(se14[[lv]]),
        rds_lo_2014 = pmax(0, unname(p14$p[[lv]]) - 1.96 * se14[[lv]]),
        rds_hi_2014 = pmin(1, unname(p14$p[[lv]]) + 1.96 * se14[[lv]]),
        rds_p_2023 = unname(p23$p[[lv]]),
        rds_se_2023 = unname(se23[[lv]]),
        rds_lo_2023 = pmax(0, unname(p23$p[[lv]]) - 1.96 * se23[[lv]]),
        rds_hi_2023 = pmin(1, unname(p23$p[[lv]]) + 1.96 * se23[[lv]]),
        diff = wd$diff,
        diff_lo = wd$lo,
        diff_hi = wd$hi,
        p_category = wd$p,
        p_overall = p_over
      )
    }
  }
}

tab_city <- bind_rows(rows) |>
  left_join(var_meta, by = "variable") |>
  mutate(
    prev_2014 = fmt_pct_ci(rds_p_2014, rds_lo_2014, rds_hi_2014),
    prev_2023 = fmt_pct_ci(rds_p_2023, rds_lo_2023, rds_hi_2023),
    n_pct_2014 = sprintf("%s; %s", n_2014, prev_2014),
    n_pct_2023 = sprintf("%s; %s", n_2023, prev_2023),
    diff_pp = sprintf("%.1f (%.1f to %.1f)", 100 * diff, 100 * diff_lo, 100 * diff_hi),
    p_cat_lab = fmt_p(p_category),
    p_lab = fmt_p(p_overall)
  )

## ---- Population-weighted two-city pool (secondary; city not tested) ---------

pool_rows <- list()
for (v in names(var_levels)) {
  levs <- var_levels[[v]]
  k <- length(levs)
  ## Independent city estimates pooled with that year's PWID population mix
  get_pool <- function(yr) {
    N <- pop_n[[yr]]
    w <- N / sum(N)
    p <- matrix(NA_real_, 2, k, dimnames = list(cities, levs))
    se <- p
    n <- integer(k)
    n_tot <- 0L
    for (ct in cities) {
      pt <- city_boot[[ct]][[yr]]$point[[v]]
      sb <- se_from_boot(city_boot[[ct]][[yr]]$boot[[v]])
      p[ct, ] <- pt$p
      se[ct, ] <- sb
      n <- n + as.integer(pt$n)
      n_tot <- n_tot + pt$n_tot
    }
    names(n) <- levs
    p_hat <- as.numeric(w["MAPUTO"] * p["MAPUTO", ] + w["NAMPULA"] * p["NAMPULA", ])
    se_hat <- sqrt((w["MAPUTO"] * se["MAPUTO", ])^2 + (w["NAMPULA"] * se["NAMPULA", ])^2)
    names(p_hat) <- levs
    names(se_hat) <- levs
    ## Bootstrap covariance of the pooled vector (delta method from city covs)
    Vmap <- cov_from_boot(city_boot[["MAPUTO"]][[yr]]$boot[[v]])
    Vnam <- cov_from_boot(city_boot[["NAMPULA"]][[yr]]$boot[[v]])
    V <- (w["MAPUTO"]^2) * Vmap + (w["NAMPULA"]^2) * Vnam
    list(p = p_hat, se = se_hat, V = V, n = n, n_tot = n_tot)
  }
  a14 <- get_pool("2014")
  a23 <- get_pool("2023")
  p_over <- overall_wald(a14$p, a14$V, a23$p, a23$V)
  for (lv in levs) {
    wd <- wald_diff(a23$p[[lv]], a23$se[[lv]], a14$p[[lv]], a14$se[[lv]])
    pool_rows[[length(pool_rows) + 1]] <- tibble(
      city = "MAPUTO+NAMPULA",
      variable = v,
      category = lv,
      cat_n = match(lv, levs),
      n_2014 = unname(a14$n[[lv]]),
      n_tot_2014 = a14$n_tot,
      n_2023 = unname(a23$n[[lv]]),
      n_tot_2023 = a23$n_tot,
      rds_p_2014 = unname(a14$p[[lv]]),
      rds_se_2014 = unname(a14$se[[lv]]),
      rds_lo_2014 = pmax(0, unname(a14$p[[lv]]) - 1.96 * a14$se[[lv]]),
      rds_hi_2014 = pmin(1, unname(a14$p[[lv]]) + 1.96 * a14$se[[lv]]),
      rds_p_2023 = unname(a23$p[[lv]]),
      rds_se_2023 = unname(a23$se[[lv]]),
      rds_lo_2023 = pmax(0, unname(a23$p[[lv]]) - 1.96 * a23$se[[lv]]),
      rds_hi_2023 = pmin(1, unname(a23$p[[lv]]) + 1.96 * a23$se[[lv]]),
      diff = wd$diff,
      diff_lo = wd$lo,
      diff_hi = wd$hi,
      p_category = wd$p,
      p_overall = p_over
    )
  }
}

tab_pool <- bind_rows(pool_rows) |>
  left_join(var_meta, by = "variable") |>
  mutate(
    prev_2014 = fmt_pct_ci(rds_p_2014, rds_lo_2014, rds_hi_2014),
    prev_2023 = fmt_pct_ci(rds_p_2023, rds_lo_2023, rds_hi_2023),
    n_pct_2014 = sprintf("%s; %s", n_2014, prev_2014),
    n_pct_2023 = sprintf("%s; %s", n_2023, prev_2023),
    diff_pp = sprintf("%.1f (%.1f to %.1f)", 100 * diff, 100 * diff_lo, 100 * diff_hi),
    p_cat_lab = fmt_p(p_category),
    p_lab = fmt_p(p_overall)
  )

tab_all <- bind_rows(tab_city, tab_pool)
write_csv_safe(tab_all, file.path(path_tab, "table1_rdsii_composition.csv"))

## ---- Publication Table 1: city-stratified, shortened, CIs -------------------

wide_city <- function(tab, panel_keep = "main") {
  levs <- var_meta$variable[var_meta$panel == panel_keep]
  sub <- tab |>
    filter(city %in% cities, panel == panel_keep) |>
    mutate(variable = factor(variable, levels = levs)) |>
    arrange(variable, cat_n, city)
  mapu <- sub |> filter(city == "MAPUTO")
  namp <- sub |> filter(city == "NAMPULA")
  stopifnot(nrow(mapu) == nrow(namp), all(as.character(mapu$category) == as.character(namp$category)))
  tibble(
    Characteristic = ifelse(!duplicated(mapu$variable), mapu$label, ""),
    Category = as.character(mapu$category),
    `Maputo 2014 n; RDS% (95% CI)` = mapu$n_pct_2014,
    `Maputo 2023 n; RDS% (95% CI)` = mapu$n_pct_2023,
    `Maputo p` = ifelse(!duplicated(mapu$variable), mapu$p_lab, ""),
    `Nampula 2014 n; RDS% (95% CI)` = namp$n_pct_2014,
    `Nampula 2023 n; RDS% (95% CI)` = namp$n_pct_2023,
    `Nampula p` = ifelse(!duplicated(namp$variable), namp$p_lab, "")
  )
}

tex_table1 <- wide_city(tab_city, "main")
tex_supp   <- wide_city(tab_city, "supplement")

write_csv_safe(tex_table1, file.path(path_tab, "table1_rdsii_main.csv"))
write_csv_safe(tex_supp,   file.path(path_tab, "table1_rdsii_supplement.csv"))

tex_pool_main <- tab_pool |>
  filter(panel == "main") |>
  mutate(variable = factor(variable, levels = var_meta$variable[var_meta$panel == "main"])) |>
  arrange(variable, cat_n) |>
  transmute(
    Characteristic = ifelse(!duplicated(variable), label, ""),
    Category = category,
    `2014 n; RDS% (95% CI)` = n_pct_2014,
    `2023 n; RDS% (95% CI)` = n_pct_2023,
    `Difference pp (95% CI)` = diff_pp,
    `p` = ifelse(!duplicated(variable), p_lab, "")
  )
write_csv_safe(tex_pool_main, file.path(path_tab, "table1_rdsii_pooled.csv"))

## ---- LaTeX fragments --------------------------------------------------------

to_xtable <- function(df, file, caption, label) {
  df_out <- df
  for (nm in names(df_out)) {
    if (is.character(df_out[[nm]])) df_out[[nm]] <- escape_tex(df_out[[nm]])
  }
  xt <- xtable::xtable(df_out, caption = caption, label = label)
  sink(file)
  print(
    xt,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    caption.placement = "top",
    size = "\\footnotesize",
    floating = TRUE,
    table.placement = "htbp"
  )
  sink()
}

## ---- Publication-quality LaTeX (booktabs + makecell) ------------------------

cell_tex <- function(n, pct_point) {
  sprintf("%s (%.1f)", n, 100 * pct_point)
}

write_city_tex <- function(tab, panel_keep, file, caption, label, notes) {
  levs <- var_meta$variable[var_meta$panel == panel_keep]
  sub <- tab |>
    filter(city %in% cities, panel == panel_keep) |>
    mutate(variable = factor(variable, levels = levs)) |>
    arrange(variable, cat_n)
  mapu <- sub |> filter(city == "MAPUTO")
  namp <- sub |> filter(city == "NAMPULA")
  lines <- c(
    "% Generated by PID_TABLE1_RDS.R",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\begin{threeparttable}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\small",
    "\\setlength{\\tabcolsep}{3.2pt}",
    "\\begin{tabular}{llcccccc}",
    "\\toprule",
    " &  & \\multicolumn{3}{c}{\\textbf{Maputo}} & \\multicolumn{3}{c}{\\textbf{Nampula}} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
    "Characteristic & Category & 2014 & 2023 & $p$ & 2014 & 2023 & $p$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(mapu))) {
    char <- if (!duplicated(mapu$variable)[i]) escape_tex(mapu$label[i]) else ""
    p_m <- if (!duplicated(mapu$variable)[i]) mapu$p_lab[i] else ""
    p_n <- if (!duplicated(namp$variable)[i]) namp$p_lab[i] else ""
    p_m <- gsub("<", "$<$", p_m, fixed = TRUE)
    p_n <- gsub("<", "$<$", p_n, fixed = TRUE)
    cat_lab <- dplyr::recode(
      as.character(mapu$category[i]),
      "Primary or none" = "Primary/none",
      "Secondary or higher" = "Secondary+",
      "Married or union" = "Married/union",
      .default = as.character(mapu$category[i])
    )
    if (!duplicated(mapu$variable)[i] && i > 1) {
      lines <- c(lines, "\\addlinespace")
    }
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      char,
      escape_tex(cat_lab),
      cell_tex(mapu$n_2014[i], mapu$rds_p_2014[i]),
      cell_tex(mapu$n_2023[i], mapu$rds_p_2023[i]),
      p_m,
      cell_tex(namp$n_2014[i], namp$rds_p_2014[i]),
      cell_tex(namp$n_2023[i], namp$rds_p_2023[i]),
      p_n
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    notes,
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
  writeLines(lines, file)
}

notes_main <- paste(
  "\\item Cells are unweighted $n$ (RDS-II weighted \\%). Percentages are Volz--Heckathorn shares, not sample proportions.",
  "95\\% bootstrap CIs are in Supplementary Table~\\ref{tab:comp-supp} and \\texttt{table1\\_rdsii\\_composition.csv}.",
  "\\item $p$-values are multivariate Wald tests of the full categorical distribution (2023 vs 2014),",
  "independent RDS samples. Reported as $<$0.001 when below 0.001.",
  "\\item Harmonised sample: Maputo $n=335$ (2014), $500$ (2023); Nampula $n=132$ (2014), $427$ (2023);",
  "age $\\ge 18$ and injection in the previous 12 months. Item non-response reduces some denominators.",
  collapse = " "
)

notes_supp <- paste(
  "\\item Same estimator as Table~\\ref{tab:comp}.",
  "\\item Used-needle and transactional-sex items have incompatible wording or recall windows across rounds;",
  "the apparent change should not be interpreted as a behavioural trend.",
  collapse = " "
)

write_city_tex(
  tab_city, "main",
  file.path(path_tex, "table1_composition.tex"),
  caption = "RDS-II (Volz--Heckathorn) weighted characteristics of people who inject drugs in the harmonised comparable sample, Maputo and Nampula, Mozambique, 2014 and 2023.",
  label = "tab:comp",
  notes = notes_main
)
## Keep the filename expected by the revision document
file.copy(
  file.path(path_tex, "table1_composition.tex"),
  file.path(path_tex, "table_composition.tex"),
  overwrite = TRUE
)
write_city_tex(
  tab_city, "supplement",
  file.path(path_tex, "table1_composition_supplement.tex"),
  caption = "Additional RDS-II weighted characteristics in the comparable sample. Used-needle and transactional-sex items are not comparable across survey instruments.",
  label = "tab:comp-supp",
  notes = notes_supp
)
to_xtable(
  tex_pool_main,
  file.path(path_tex, "table1_composition_pooled.tex"),
  caption = "Population-weighted two-city RDS-II composition (secondary). City mixing weights are the published PWID population-size estimates for each year, not sample size. This table is not a like-with-like geographic comparison; city-stratified estimates are primary.",
  label = "tab:comp-pool"
)

## ---- Key numeric notes for the manuscript -----------------------------------

key_notes <- tab_city |>
  filter(panel == "main") |>
  select(city, label, category, n_2014, n_2023, prev_2014, prev_2023, diff_pp, p_lab)

write_csv_safe(key_notes, file.path(path_tab, "table1_key_notes.csv"))

message("Table 1 written to ", path_tab)
print(tex_table1, n = 80)
message("\n--- Pooled (secondary) ---")
print(tex_pool_main, n = 40)
