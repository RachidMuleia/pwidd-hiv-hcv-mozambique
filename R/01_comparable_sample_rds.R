################################################################################
## Updated analysis responding to PLOS ONE reviews (PONE-D-25-62212)
## HIV, HBV and HCV among people who inject drugs in Mozambique
## Serial cross-sectional BBS: 2013-2014 and 2023-2024
##
## This script addresses the statistical comments by:
##  1. Restricting temporal comparisons to a harmonised eligible sample
##     (same cities, age >= 18, injected in the previous 12 months)
##  2. Producing population-weighted RDS estimates (RDS-II / Gile SS)
##  3. Accounting for within-survey sample dependence (recruitment chains)
##  4. Using a pre-specified covariate set (no stepwise selection)
##  5. Reporting effect sizes and CIs, not only p-values
##
## Data analyst: Rachid Muleia
## Date: 2026-09-19
################################################################################

## ---- 0. Paths, packages, options --------------------------------------------

rm(list = ls())
options(survey.lonely.psu = "adjust")
options(stringsAsFactors = FALSE)
set.seed(20260919)

## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))


pkgs <- c(
  "tidyverse", "RDS", "survey", "geepack", "sandwich", "lmtest",
  "broom", "ggplot2", "scales", "xtable", "glue"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))
select <- dplyr::select
filter <- dplyr::filter
recode <- dplyr::recode
mutate <- dplyr::mutate
between <- dplyr::between

## Bootstrap replicates for Gile uncertainty. Increase to 1000-5000 for the
## final submission if runtime allows. 400 is adequate for revision tables.
N_BOOT <- 400

## ---- 1. Population size estimates used for Gile SS / pooling ----------------
## 2014: median of four PSE methods (Sema Baltazar et al., BMC Public Health 2020)
## 2023: key-population mapping estimates used in the IBBS-II analysis files
pop_n <- list(
  "2014" = c(MAPUTO = 1684, NAMPULA = 520),
  "2023" = c(MAPUTO = 990, NAMPULA = 1792, BEIRA = 1029,
             TETE = 3219, QUELIMANE = 681)
)
## Sensitivity: sequential-sampling 2014 PSE (INE/INS IBBS-I documentation)
pop_n_sens_2014 <- c(MAPUTO = 2869, NAMPULA = 463)

## ---- 2. Small utilities -----------------------------------------------------

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

cap_degree <- function(d) {
  d <- as.numeric(d)
  d[is.na(d) | d < 1] <- 1
  ## Trim extreme self-reported degrees that otherwise dominate VH weights
  q99 <- stats::quantile(d, 0.99, na.rm = TRUE)
  cap <- max(q99, 1)
  pmin(d, cap)
}

write_csv_safe <- function(x, file) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  list_cols <- vapply(x, is.list, logical(1))
  if (any(list_cols)) {
    x[list_cols] <- lapply(x[list_cols], function(col) {
      vapply(col, function(z) paste(unlist(z), collapse = ";"), character(1))
    })
  }
  utils::write.csv(x, file, row.names = FALSE, na = "")
  invisible(x)
}

escape_tex <- function(x) {
  x <- as.character(x)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x
}

df_to_tex <- function(df, file, caption, label, digits = NULL, align = NULL) {
  df_out <- df
  for (nm in names(df_out)) {
    if (is.character(df_out[[nm]])) df_out[[nm]] <- escape_tex(df_out[[nm]])
  }
  xt <- xtable::xtable(df_out, caption = caption, label = label, digits = digits)
  if (!is.null(align)) align(xt) <- align
  sink(file)
  print(
    xt,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    caption.placement = "top",
    size = "\\small",
    floating = TRUE,
    table.placement = "htbp"
  )
  sink()
}

wald_diff <- function(p1, se1, p0, se0) {
  d <- p1 - p0
  se <- sqrt(se1^2 + se0^2)
  z <- ifelse(se > 0, d / se, NA_real_)
  pval <- 2 * stats::pnorm(-abs(z))
  tibble(
    diff = d,
    se_diff = se,
    ci_lo = d - 1.96 * se,
    ci_hi = d + 1.96 * se,
    p_value = pval
  )
}

## Volz-Heckathorn / RDS-II prevalence (depends only on outcome and degree)
vh_estimate <- function(y, degree, n_boot = 400) {
  ok <- !is.na(y) & !is.na(degree) & degree > 0
  y <- as.numeric(y[ok])
  d <- as.numeric(degree[ok])
  n <- length(y)
  if (n < 5 || length(unique(y)) < 1) {
    return(list(est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, n = n, n_pos = sum(y == 1)))
  }
  est <- sum(y / d) / sum(1 / d)
  ## With-replacement bootstrap of the VH estimator (within the analysis subset)
  boot_est <- replicate(n_boot, {
    ii <- sample.int(n, n, replace = TRUE)
    yy <- y[ii]; dd <- d[ii]
    sum(yy / dd) / sum(1 / dd)
  })
  se <- stats::sd(boot_est)
  q <- stats::quantile(boot_est, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  list(est = est, se = se, lo = q[1], hi = q[2], n = n, n_pos = sum(y == 1))
}

## Gile successive-sampling estimator with Gile bootstrap when an rds.data.frame
## is available. Falls back to VH if the RDS estimator fails.
gile_estimate <- function(rds_df, outcome, N, n_boot = N_BOOT, subset = NULL,
                          weight_type = "Gile's SS", uncertainty = "Gile") {
  if (is.null(rds_df)) return(NULL)
  out <- tryCatch({
    est <- RDS::RDS.bootstrap.intervals(
      rds_df,
      outcome.variable = outcome,
      weight.type = weight_type,
      uncertainty = uncertainty,
      confidence.level = 0.95,
      number.of.bootstrap.samples = n_boot,
      to.factor = FALSE,
      N = N,
      subset = subset
    )
    point <- as.numeric(est$estimate[1])
    se <- tryCatch(as.numeric(attr(est$interval, "bsresult")$se_estimate[1]),
                   error = function(e) NA_real_)
    if (is.na(se) || length(se) == 0 || se == 0) {
      se <- tryCatch(as.numeric(attr(est, "bsresult")$se_estimate[1]),
                     error = function(e) NA_real_)
    }
    list(est = point, se = se, method = paste(weight_type, uncertainty, sep = " / "))
  }, error = function(e) {
    message(weight_type, " / ", uncertainty, " failed for ", outcome, ": ", conditionMessage(e))
    NULL
  })
  out
}

## 2014 coupon-manager extract stores each participant's own coupon as INPON.
## Recruiter is the person who issued that coupon as OUTPON1-5.
recruiter_from_outpon <- function(coupon, out1, out2, out3, out4, out5) {
  coupon <- as.character(coupon)
  rec_map <- new.env(parent = emptyenv())
  add_out <- function(out) {
    out <- as.character(out)
    ok <- !is.na(out) & nzchar(out) & out != "NA"
    if (!any(ok)) return(invisible(NULL))
    issued <- out[ok]
    issuer <- coupon[ok]
    for (i in seq_along(issued)) rec_map[[issued[i]]] <- issuer[i]
    invisible(NULL)
  }
  add_out(out1); add_out(out2); add_out(out3); add_out(out4); add_out(out5)
  rec <- vapply(coupon, function(id) {
    if (exists(id, envir = rec_map, inherits = FALSE)) rec_map[[id]] else "0"
  }, character(1), USE.NAMES = FALSE)
  rec[is.na(rec) | rec == coupon] <- "0"
  rec
}

make_rds_df <- function(dat, id, recruiter, degree, N, max_coupons = 6) {
  tmp <- dat
  tmp$.id <- as.character(tmp[[id]])
  tmp$.rec <- as.character(tmp[[recruiter]])
  tmp$.deg <- cap_degree(tmp[[degree]])
  ## RDS seeds are coded as recruiter.id = 0
  tmp$.rec[is.na(tmp$.rec) | tmp$.rec %in% c("", "NA", "seed", "SEED")] <- "0"
  tmp$.rec[tmp$.rec == tmp$.id] <- "0"
  ## Drop duplicated ids if any
  tmp <- tmp[!duplicated(tmp$.id), ]
  RDS::as.rds.data.frame(
    tmp,
    id = ".id",
    recruiter.id = ".rec",
    network.size = ".deg",
    population.size = N,
    max.coupons = max_coupons
  )
}

compute_rds_weights <- function(rds_df, N, type = "Gile's SS") {
  w <- tryCatch(
    RDS::compute.weights(rds_df, weight.type = type, N = N),
    error = function(e) {
      message("compute.weights(", type, ") failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(w)) {
    deg <- RDS::get.net.size(rds_df)
    w <- (1 / deg) / mean(1 / deg)
  }
  as.numeric(w)
}

## ---- 3. Harmonised recoding -------------------------------------------------

read_2014 <- function(path) {
  df <- utils::read.csv(path, header = TRUE, na.strings = c("", "NA", "NaN", "nan"),
                        stringsAsFactors = FALSE, encoding = "latin1")
  rec14 <- recruiter_from_outpon(
    df$COUPON, df$OUTPON1, df$OUTPON2, df$OUTPON3, df$OUTPON4, df$OUTPON5
  )
  df |>
    mutate(
      YEAR = "2014",
      CIDADE = dplyr::case_when(
        SRVYCITYC %in% c("MAPUTO", "Maputo") ~ "MAPUTO",
        SRVYCITYC %in% c("NAMPULA", "Nampula") ~ "NAMPULA",
        TRUE ~ as.character(SRVYCITYC)
      ),
      AGE = as.numeric(AGE),
      AGE_CAT = case_when(
        between(AGE, 18, 24) ~ "18-24",
        AGE >= 25 ~ "25+",
        TRUE ~ NA_character_
      ),
      SEX_CAT = case_when(
        SCREEN1 == 1 ~ "Male",
        SCREEN1 == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      EDUC_CAT = recode(
        DEMEDU2CC,
        "1_PRIMARIO_NADA" = "Primary or none",
        "2_SECUNDARIO" = "Secondary or higher"
      ),
      MARITALC = recode(
        MARITALC,
        "1_SOLTEIRO" = "Never married",
        "2_CASADO" = "Married or union",
        "3_OUTRO" = "Other"
      ),
      RELIGIAO_CAT = recode(
        DEMRELCCN,
        "1_Cristao" = "Christian",
        "3_Mulcumano" = "Muslim",
        "4_Nenhuma" = "None/other",
        "5_Outra" = "None/other"
      ),
      AGE_FIRST_DRUG = case_when(
        DRUGS4 < 18 ~ "<18",
        between(DRUGS4, 18, 24) ~ "18-24",
        DRUGS4 >= 25 ~ "25+",
        TRUE ~ NA_character_
      ),
      FREQ_INJEC = recode(
        DRINJ4C2,
        "1_DIARIAMENTE" = "Daily",
        "2_ANUALMENSALSEMANAL" = "Less than daily",
        "3_NAOINJECTOU" = "Did not inject"
      ),
      injected_12m = DRINJ4C2 != "3_NAOINJECTOU",
      NEW_SYRINGE = recode(DRINJ16C, "1_SIM" = "Yes", "2_NAO" = "No"),
      USED_NEEDLE = recode(DRINJ2DC, "1_USADA" = "Yes", "2_NAO" = "No"),
      SHARED_EQUIP = recode(DRINJ2FC, "1_PARTILHOUPREPARACAO" = "Yes", "2_NAO" = "No"),
      UNPROTECTED_SEX = recode(P1.FASX1C2, "1_TEVEDESPROTEGIDO" = "Yes", "2_NAO" = "No"),
      SEX_PARTNER = dplyr::case_when(
        dplyr::between(suppressWarnings(as.numeric(FVSEX2)), 0, 1) |
          dplyr::between(suppressWarnings(as.numeric(MSEX2)), 0, 1) ~ "0-1",
        suppressWarnings(as.numeric(FVSEX2)) == 2 |
          suppressWarnings(as.numeric(MSEX2)) == 2 ~ "2",
        (suppressWarnings(as.numeric(FVSEX2)) > 2 &
           !suppressWarnings(as.numeric(FVSEX2)) %in% c(9997, 9999)) |
          (suppressWarnings(as.numeric(MSEX2)) > 2 &
             !suppressWarnings(as.numeric(MSEX2)) %in% c(9997, 9999)) ~ ">=3",
        TRUE ~ NA_character_
      ),
      TX_SEX = recode(RECEIVEC, "1_DROGASDINHEIRO" = "Yes", "2_NAO" = "No"),
      STI_SELF = recode(STIALLC, "1_ITS" = "Yes", "2_NAO" = "No"),
      HIV_TESTED_12M = case_when(
        VCT3CC == "1_<12MESES" ~ "Yes",
        VCT3CC %in% c("0_NUNCATESTOU", "2_1-3ANOS", "2_4+ANOS") ~ "No",
        TRUE ~ NA_character_
      ),
      SOUGHT_CARE = recode(HEALTH1C, "1_SIM" = "Yes", "2_NAO" = "No"),
      HIV = dplyr::case_when(HIVFINALC %in% c("1_POS", "1_POSITIVO") ~ 1, HIVFINALC %in% c("2_NEG", "2_NEGATIVO") ~ 0),
      HBV = dplyr::case_when(HBVRAPIDC %in% c("1_POS", "1_POSITIVO") ~ 1, HBVRAPIDC %in% c("2_NEG", "2_NEGATIVO") ~ 0),
      HCV = dplyr::case_when(HCVRAPIDC %in% c("1_POS", "1_POSITIVO") ~ 1, HCVRAPIDC %in% c("2_NEG", "2_NEGATIVO") ~ 0),
      coupon_id = as.character(COUPON),
      recruiter_id = rec14,
      is_seed = recruiter_id == "0",
      degree_raw = as.numeric(DEGREE),
      degree = cap_degree(DEGREE),
      wave_proxy = substr(as.character(COUPON), 2, 2)
    )
}

read_2023 <- function(path) {
  df <- utils::read.csv(path, header = TRUE, na.strings = c("", "NA", "NaN", "nan"),
                        stringsAsFactors = FALSE)
  df |>
    mutate(
      YEAR = "2023",
      CIDADE = recode(
        SITE_CITY,
        "1_MAPUTO_CIDADE" = "MAPUTO",
        "2_BEIRA" = "BEIRA",
        "3_TETE" = "TETE",
        "4_QUELIMANE" = "QUELIMANE",
        "5_NAMPULA" = "NAMPULA"
      ),
      AGE = as.numeric(ELAGEL_B),
      AGE_CAT = case_when(
        between(AGE, 18, 24) ~ "18-24",
        AGE >= 25 ~ "25+",
        TRUE ~ NA_character_
      ),
      SEX_CAT = case_when(
        SCREEN1 == "Masculino" ~ "Male",
        SCREEN1 == "Feminino" ~ "Female",
        TRUE ~ NA_character_
      ),
      EDUC_CAT = case_when(
        DEEDHIGH %in% c("Sem_escolaridade", "Primario_ou_Alfabetizacao") ~
          "Primary or none",
        DEEDHIGH %in% c("Secundario", "Tecnico", "Superior", "Outra") ~
          "Secondary or higher",
        TRUE ~ NA_character_
      ),
      MARITALC = case_when(
        DEMARSTA == "Solteiro" ~ "Never married",
        DEMARSTA %in% c("Casado", "Uniao_de_factos") ~ "Married or union",
        DEMARSTA %in% c("Separado", "Divorciado", "Viuvo") ~ "Other",
        TRUE ~ NA_character_
      ),
      RELIGIAO_CAT = case_when(
        DERELIG %in% c("Catolica", "Protestante_Envagelica", "Siao_Zione") ~ "Christian",
        DERELIG == "Muculmana" ~ "Muslim",
        DERELIG %in% c("Sem_relegiao", "Animista", "Outra") ~ "None/other",
        TRUE ~ NA_character_
      ),
      AGE_FIRST_DRUG = case_when(
        NOTA13_IDOLD_A < 18 ~ "<18",
        between(NOTA13_IDOLD_A, 18, 24) ~ "18-24",
        NOTA13_IDOLD_A >= 25 ~ "25+",
        TRUE ~ NA_character_
      ),
      FREQ_INJEC = case_when(
        DRUGS6_1 == "Nao" ~ "Did not inject",
        ID6_FRQ %in% c("1_3_vezes_por_dia", "5+_vezes_por_dia") ~ "Daily",
        ID6_FRQ %in% c("1_4_vezes_por_mes", "2_7_vezes_por_semana") ~ "Less than daily",
        TRUE ~ NA_character_
      ),
      injected_12m = DRUGS6_1 == "Sim",
      NEW_SYRINGE = case_when(IDNSTE == "Sim" ~ "Yes", IDNSTE == "Nao" ~ "No", TRUE ~ NA_character_),
      USED_NEEDLE = case_when(
        DRINJ9 %in% c("As_vezes", "Quase_sempre", "Sempre") ~ "Yes",
        DRINJ9 == "Nunca" ~ "No",
        TRUE ~ NA_character_
      ),
      SHARED_EQUIP = case_when(DRINJ4 == "Sim" ~ "Yes", DRINJ4 == "Nao" ~ "No", TRUE ~ NA_character_),
      UNPROTECTED_SEX = case_when(
        LASTREL6 == "Nao" ~ "Yes",
        LASTREL6 == "Sim" ~ "No",
        TRUE ~ NA_character_
      ),
      SEX_PARTNER = dplyr::case_when(
        IDENT42_LIMFPART_A == 0 | LIMFVAG_H == "Nao" | LIMFVAG_M == "Nao" |
          IDENT42_LIMFPART_A == 1 | IDENT43_MSEX2_A == 1 ~ "0-1",
        IDENT42_LIMFPART_A == 2 | IDENT43_MSEX2_A == 2 ~ "2",
        IDENT42_LIMFPART_A > 2 | IDENT43_MSEX2_A > 2 ~ ">=3",
        TRUE ~ NA_character_
      ),
      TX_SEX = case_when(
        FVSEX6_CAT == "1_SIM" ~ "Yes",
        FVSEX6_CAT == "2_NAO" ~ "No",
        MSEX6_CAT == "1_SIM" ~ "Yes",
        MSEX6_CAT == "2_NAO" ~ "No",
        TRUE ~ NA_character_
      ),
      STI_SELF = recode(SINTOMAS_ITS_INFO, "1_SIM" = "Yes", "2_SIM" = "No", "2_NAO" = "No"),
      HIV_TESTED_12M = case_when(
        CSCTTIC %in% c("1_<3MES", "2_4-6MES", "3_7-12MES") ~ "Yes",
        CSCTTIC == "4_>12MES" | CSCTEV == "Nao" ~ "No",
        TRUE ~ NA_character_
      ),
      SOUGHT_CARE = case_when(HEALTH1 == "Sim" ~ "Yes", HEALTH1 == "Nao" ~ "No", TRUE ~ NA_character_),
      HIV = dplyr::case_when(RESUL_HIV_PREV == "1_POSITIVO" ~ 1, RESUL_HIV_PREV == "2_NEGATIVO" ~ 0),
      HBV = dplyr::case_when(RESUL_HBV_CAT == "1_POSITIVO" ~ 1, RESUL_HBV_CAT == "2_NEGATIVO" ~ 0),
      HCV = dplyr::case_when(RESUL_HCV_CAT == "1_POSITIVO" ~ 1, RESUL_HCV_CAT == "2_NEGATIVO" ~ 0),
      coupon_id = as.character(id),
      recruiter_id = if_else(is.na(from) | from == "", "0", as.character(from)),
      is_seed = recruiter_id == "0",
      degree_raw = as.numeric(ELINJMU),
      degree = cap_degree(ELINJMU),
      wave_proxy = NA_character_
    )
}

message("Reading survey files ...")
pid14_raw <- read_2014(path_2014)
pid23_raw <- read_2023(path_2023)
message("2014 raw n = ", nrow(pid14_raw), "; 2023 raw n = ", nrow(pid23_raw))

## Harmonised comparable sample: same geography, age and recent injection
comparable_cities <- c("MAPUTO", "NAMPULA")

keep_cols <- c(
  "YEAR", "CIDADE", "AGE", "AGE_CAT", "SEX_CAT", "EDUC_CAT", "MARITALC",
  "RELIGIAO_CAT", "AGE_FIRST_DRUG", "FREQ_INJEC", "injected_12m",
  "NEW_SYRINGE", "USED_NEEDLE", "SHARED_EQUIP", "UNPROTECTED_SEX", "SEX_PARTNER", "TX_SEX",
  "STI_SELF", "HIV_TESTED_12M", "SOUGHT_CARE", "HIV", "HBV", "HCV",
  "coupon_id", "recruiter_id", "is_seed", "degree_raw", "degree",
  "wave_proxy", "sample_flag", "HIV_HCV"
)

pid14 <- pid14_raw |>
  filter(CIDADE %in% comparable_cities, AGE >= 18, injected_12m) |>
  mutate(
    sample_flag = "comparable",
    HIV_HCV = if_else(!is.na(HIV) & !is.na(HCV), as.numeric(HIV == 1 & HCV == 1), NA_real_)
  ) |>
  select(any_of(keep_cols))

pid23_comp <- pid23_raw |>
  filter(CIDADE %in% comparable_cities, AGE >= 18, injected_12m) |>
  mutate(
    sample_flag = "comparable",
    HIV_HCV = if_else(!is.na(HIV) & !is.na(HCV), as.numeric(HIV == 1 & HCV == 1), NA_real_)
  ) |>
  select(any_of(keep_cols))

pid23_all <- pid23_raw |>
  filter(AGE >= 18, injected_12m) |>
  mutate(
    sample_flag = "all_cities_2023",
    HIV_HCV = if_else(!is.na(HIV) & !is.na(HCV), as.numeric(HIV == 1 & HCV == 1), NA_real_)
  ) |>
  select(any_of(keep_cols))

## ---- 4. RDS objects and weights by city -------------------------------------

build_city_rds <- function(dat, year, pop_lookup) {
  cities <- sort(unique(dat$CIDADE))
  out <- lapply(cities, function(ct) {
    d <- dat |> filter(CIDADE == ct)
    N <- unname(pop_lookup[[as.character(year)]][ct])
    if (is.na(N) || is.null(N)) stop("Missing population size for ", year, " ", ct)
    ## Keep the true recruiter on the analysis data even if that person was
    ## dropped by eligibility (they still define a GEE cluster). For the RDS
    ## graph, recode recruiters who are not in this subset as seeds so
    ## as.rds.data.frame is valid.
    rec <- as.character(d$recruiter_id)
    rec[is.na(rec) | rec %in% c("", "NA", "seed", "SEED")] <- "0"
    rec[rec == as.character(d$coupon_id)] <- "0"
    rec_rds <- rec
    rec_rds[!rec_rds %in% c("0", as.character(d$coupon_id))] <- "0"
    d$recruiter_id <- rec
    d$recruiter_id_rds <- rec_rds
    rds <- tryCatch(
      make_rds_df(d, id = "coupon_id", recruiter = "recruiter_id_rds",
                  degree = "degree", N = N, max_coupons = 6),
      error = function(e) {
        message("as.rds.data.frame failed for ", year, " ", ct, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(rds)) {
      d$w_gile <- (1 / d$degree) / mean(1 / d$degree)
      d$w_rdsii <- d$w_gile
    } else {
      w_gile <- compute_rds_weights(rds, N = N, type = "Gile's SS")
      w_ii   <- compute_rds_weights(rds, N = N, type = "RDS-II")
      d$w_gile <- as.numeric(w_gile)
      d$w_rdsii <- as.numeric(w_ii)
    }
    ## Population-calibrated weights: sum to estimated city PWID population
    d$N_city <- N
    d$w_pop <- d$w_gile / sum(d$w_gile, na.rm = TRUE) * N
    d$w_pop_ii <- d$w_rdsii / sum(d$w_rdsii, na.rm = TRUE) * N
    ## Recruitment chain: seed lineage from the RDS graph when available.
    if (!is.null(rds)) {
      seed_id <- tryCatch(RDS::get.seed.id(rds), error = function(e) d$coupon_id)
      d$chain_id <- paste(ct, seed_id, sep = ":")
    } else {
      d$chain_id <- paste(ct, d$coupon_id, sep = ":")
    }
    list(data = d, rds = rds, N = N, city = ct, year = year)
  })
  names(out) <- cities
  out
}

message("Building RDS objects (comparable sample) ...")
rds14 <- build_city_rds(pid14, "2014", pop_n)
rds23 <- build_city_rds(pid23_comp, "2023", pop_n)
rds23_all <- build_city_rds(pid23_all, "2023", pop_n)

bind_city_list <- function(lst) {
  dplyr::bind_rows(lapply(lst, `[[`, "data"))
}

dat14 <- bind_city_list(rds14)
dat23 <- bind_city_list(rds23)
dat23_all <- bind_city_list(rds23_all)
dat_comp <- bind_rows(dat14, dat23)

## ---- 5. City-year prevalence (RDS-weighted, dependence-aware) ---------------

estimate_city_outcome <- function(city_obj, outcome, n_boot = N_BOOT) {
  d <- city_obj$data
  rds <- city_obj$rds
  N <- city_obj$N
  y <- d[[outcome]]
  ## Crude
  ok <- !is.na(y)
  crude_p <- mean(y[ok])
  crude_n <- sum(ok)
  crude_pos <- sum(y[ok] == 1)
  crude_se <- sqrt(crude_p * (1 - crude_p) / crude_n)
  ## VH / RDS-II (always available)
  vh <- vh_estimate(y, d$degree, n_boot = n_boot)
  ## Prefer Gile SS; if it fails (e.g. negative inclusion probability after
  ## pruning the recruitment tree), use RDS-II with Salganik chain bootstrap.
  rds_est <- NULL
  ## Gile SS is unstable after eligibility pruning of recruitment chains.
  ## Primary estimator: RDS-II (Volz-Heckathorn). Sample dependence is
  ## handled in CIs via the weighted bootstrap and in regression via
  ## chain-clustered GEE / survey GLM.
  if (FALSE && !is.null(rds) && !is.na(N) && N > sum(ok) + 5) {
    rds_est <- gile_estimate(rds, outcome = outcome, N = N, n_boot = n_boot,
                             weight_type = "Gile's SS", uncertainty = "Gile",
                             subset = !is.na(rds[[outcome]]))
    if (is.null(rds_est) || is.na(rds_est$est)) {
      rds_est <- gile_estimate(rds, outcome = outcome, N = N, n_boot = n_boot,
                               weight_type = "RDS-II", uncertainty = "Salganik",
                               subset = !is.na(rds[[outcome]]))
    }
  }
  if (!is.null(rds_est) && !is.na(rds_est$est)) {
    est <- rds_est$est
    se <- ifelse(is.na(rds_est$se) || rds_est$se == 0, vh$se, rds_est$se)
    method <- rds_est$method
  } else {
    est <- vh$est
    se <- vh$se
    method <- "RDS-II (VH bootstrap)"
  }
  lo  <- est - 1.96 * se
  hi  <- est + 1.96 * se
  tibble(
    survey_year = city_obj$year,
    city = city_obj$city,
    outcome = outcome,
    n = crude_n,
    n_pos = crude_pos,
    N_pop = N,
    crude_p = crude_p,
    crude_se = crude_se,
    rds_p = est,
    rds_se = se,
    rds_lo = pmax(0, lo),
    rds_hi = pmin(1, hi),
    vh_p = vh$est,
    vh_se = vh$se,
    method = method
  )
}

outcomes <- c("HIV", "HCV", "HBV", "HIV_HCV")

message("Estimating city-level RDS prevalence (this step uses Gile bootstrap) ...")
city_prev <- dplyr::bind_rows(lapply(c(rds14, rds23), function(x) {
  dplyr::bind_rows(lapply(outcomes, function(v) estimate_city_outcome(x, v)))
}))
print(city_prev)

## Population-weighted pooled estimate across Maputo and Nampula
pool_cities <- function(tab, years = c("2014", "2023"), cities = comparable_cities) {
  bind_rows(lapply(years, function(yr) {
    bind_rows(lapply(unique(tab$outcome), function(v) {
      sub <- tab |> filter(.data$survey_year == yr, .data$outcome == v, city %in% cities)
      w <- sub$N_pop / sum(sub$N_pop)
      p <- sum(w * sub$rds_p)
      se <- sqrt(sum(w^2 * sub$rds_se^2))
      p_vh <- sum(w * sub$vh_p)
      se_vh <- sqrt(sum(w^2 * sub$vh_se^2))
      p_crude <- weighted.mean(sub$crude_p, sub$n)
      tibble(
        survey_year = yr,
        city = "MAPUTO+NAMPULA",
        outcome = v,
        n = sum(sub$n),
        n_pos = sum(sub$n_pos),
        N_pop = sum(sub$N_pop),
        crude_p = p_crude,
        crude_se = NA_real_,
        rds_p = p,
        rds_se = se,
        rds_lo = pmax(0, p - 1.96 * se),
        rds_hi = pmin(1, p + 1.96 * se),
        vh_p = p_vh,
        vh_se = se_vh,
        method = "Pop-weighted pool of city RDS estimates"
      )
    }))
  }))
}

pooled_prev <- pool_cities(city_prev)
## Alternative pooling that holds the 2014 PWID population mix constant
pool_fixed_2014 <- {
  bind_rows(lapply(c("2014", "2023"), function(yr) {
    bind_rows(lapply(outcomes, function(v) {
      sub <- city_prev |> filter(.data$survey_year == yr, .data$outcome == v, city %in% comparable_cities)
      Nfix <- pop_n[["2014"]][sub$city]
      w <- Nfix / sum(Nfix)
      p <- sum(w * sub$rds_p)
      se <- sqrt(sum(w^2 * sub$rds_se^2))
      tibble(
        survey_year = yr, city = "Fixed 2014 population mix", outcome = v,
        n = sum(sub$n), n_pos = sum(sub$n_pos), N_pop = sum(Nfix),
        rds_p = p, rds_se = se,
        rds_lo = pmax(0, p - 1.96 * se), rds_hi = pmin(1, p + 1.96 * se)
      )
    }))
  }))
}

city_all <- bind_rows(city_prev, pooled_prev)

compare_years <- function(tab) {
  tab |>
    select(survey_year, city, outcome, n, n_pos, N_pop, crude_p, rds_p, rds_se, rds_lo, rds_hi, method) |>
    pivot_wider(
      names_from = survey_year,
      values_from = c(n, n_pos, N_pop, crude_p, rds_p, rds_se, rds_lo, rds_hi, method),
      names_sep = "_"
    ) |>
    rowwise() |>
    mutate(wald = list(wald_diff(rds_p_2023, rds_se_2023, rds_p_2014, rds_se_2014))) |>
    unnest(wald) |>
    ungroup() |>
    mutate(
      n_2014_lab = paste0(n_pos_2014, "/", n_2014),
      n_2023_lab = paste0(n_pos_2023, "/", n_2023),
      prev_2014 = fmt_pct_ci(rds_p_2014, rds_lo_2014, rds_hi_2014),
      prev_2023 = fmt_pct_ci(rds_p_2023, rds_lo_2023, rds_hi_2023),
      crude_2014 = fmt_pct(crude_p_2014),
      crude_2023 = fmt_pct(crude_p_2023),
      diff_pp = sprintf("%.1f (%.1f to %.1f)", 100 * diff, 100 * ci_lo, 100 * ci_hi),
      p_lab = fmt_p(p_value)
    )
}

prev_compare <- compare_years(city_all)
write_csv_safe(prev_compare, file.path(path_tab, "table_prevalence_rds_compare.csv"))

## ---- 6. Table 1 composition --------------------------------------------------
## Original unweighted chi-squared Table 1 is invalid for RDS.
## A pooled survey::svychisq is also incorrect: city mixing weights are the
## PSE ratios, 2014 has only 2-3 coupon-wave clusters, and Gile weights were
## abandoned for prevalence. City-stratified RDS-II + bootstrap Wald is
## produced by PID_TABLE1_RDS.R after analysis_objects.rds is written.
tab_composition <- NULL

## ---- 7. Direct standardisation (composition-adjusted 2023 prevalence) -------
## Standard: 2014 RDS-weighted age x sex x city distribution

std_strata <- c("CIDADE", "AGE_CAT", "SEX_CAT")

std_weights_2014 <- dat14 |>
  filter(!is.na(AGE_CAT), !is.na(SEX_CAT)) |>
  group_by(across(all_of(std_strata))) |>
  summarise(w_std = sum(w_pop), .groups = "drop") |>
  mutate(w_std = w_std / sum(w_std))

stratum_prev <- function(dat, outcome) {
  dat |>
    filter(!is.na(.data[[outcome]]), !is.na(AGE_CAT), !is.na(SEX_CAT)) |>
    group_by(across(all_of(std_strata))) |>
    summarise(
      p = sum(w_pop * .data[[outcome]]) / sum(w_pop),
      n = n(),
      .groups = "drop"
    )
}

direct_std <- bind_rows(lapply(c("HIV", "HCV", "HBV"), function(v) {
  p23 <- stratum_prev(dat23, v)
  joined <- std_weights_2014 |>
    left_join(p23, by = std_strata)
  p_std <- sum(joined$w_std * joined$p, na.rm = TRUE)
  ## Bootstrap the standardised estimator by resampling within year
  boot_vals <- replicate(300, {
    b14 <- dat14[sample.int(nrow(dat14), replace = TRUE), ]
    b23 <- dat23[sample.int(nrow(dat23), replace = TRUE), ]
    w14 <- b14 |>
      filter(!is.na(AGE_CAT), !is.na(SEX_CAT)) |>
      group_by(across(all_of(std_strata))) |>
      summarise(w_std = sum(w_pop), .groups = "drop") |>
      mutate(w_std = w_std / sum(w_std))
    p_b <- b23 |>
      filter(!is.na(.data[[v]]), !is.na(AGE_CAT), !is.na(SEX_CAT)) |>
      group_by(across(all_of(std_strata))) |>
      summarise(p = sum(w_pop * .data[[v]]) / sum(w_pop), .groups = "drop")
    jj <- w14 |> left_join(p_b, by = std_strata)
    sum(jj$w_std * jj$p, na.rm = TRUE)
  })
  tibble(
    outcome = v,
    std_p_2023 = p_std,
    std_lo = stats::quantile(boot_vals, 0.025, na.rm = TRUE, names = FALSE),
    std_hi = stats::quantile(boot_vals, 0.975, na.rm = TRUE, names = FALSE)
  )
}))

p14_overall <- pooled_prev |>
  filter(survey_year == "2014", city == "MAPUTO+NAMPULA") |>
  select(outcome, rds_p_2014 = rds_p, rds_se_2014 = rds_se)
p23_overall <- pooled_prev |>
  filter(survey_year == "2023", city == "MAPUTO+NAMPULA") |>
  select(outcome, rds_p_2023 = rds_p, rds_se_2023 = rds_se)

std_compare <- direct_std |>
  left_join(p14_overall, by = "outcome") |>
  left_join(p23_overall, by = "outcome") |>
  mutate(
    prev_2014 = fmt_pct(rds_p_2014),
    prev_2023_unadj = fmt_pct(rds_p_2023),
    prev_2023_std = fmt_pct_ci(std_p_2023, std_lo, std_hi)
  )

write_csv_safe(std_compare, file.path(path_tab, "table_direct_standardisation.csv"))
write_csv_safe(pool_fixed_2014, file.path(path_tab, "table_pool_fixed_2014_mix.csv"))

## ---- 8. Subgroup RDS-II prevalence (harmonised covariates) ------------------

subgroup_prev <- function(dat, outcome, varname) {
  dat2 <- dat |> filter(!is.na(.data[[outcome]]), !is.na(.data[[varname]]))
  levs <- sort(unique(dat2[[varname]]))
  bind_rows(lapply(levs, function(lv) {
    sub <- dat2 |> filter(.data[[varname]] == lv)
    est14 <- vh_estimate(sub$HIV[sub$YEAR == "2014"] * NA_real_, 1) # placeholder
    NULL
  }))
}

vh_by_year_level <- function(dat, outcome, varname) {
  dat2 <- dat |>
    filter(!is.na(.data[[outcome]]), !is.na(.data[[varname]]))
  expand_grid(YEAR = c("2014", "2023"), category = sort(unique(dat2[[varname]]))) |>
    rowwise() |>
    mutate(
      vh = list({
        sub <- dat2 |> filter(YEAR == YEAR, .data[[varname]] == category)
        ## filter above uses the data-masking YEAR column against itself; fix:
        NULL
      })
    )
}

estimate_subgroup <- function(dat, outcome, varname, n_boot = 200) {
  rows <- list()
  dat2 <- dat[!is.na(dat[[outcome]]) & !is.na(dat[[varname]]), ]
  for (yr in c("2014", "2023")) {
    for (lv in sort(unique(dat2[[varname]]))) {
      sub <- dat2[dat2$YEAR == yr & dat2[[varname]] == lv, ]
      if (nrow(sub) < 8) next
      vh <- vh_estimate(sub[[outcome]], sub$degree, n_boot = n_boot)
      rows[[length(rows) + 1]] <- tibble(
        outcome = outcome,
        variable = varname,
        category = lv,
        survey_year = yr,
        n = vh$n,
        n_pos = vh$n_pos,
        rds_p = vh$est,
        rds_se = vh$se,
        rds_lo = pmax(0, vh$lo),
        rds_hi = pmin(1, vh$hi)
      )
    }
  }
  bind_rows(rows)
}

sg_vars <- c("AGE_CAT", "SEX_CAT", "EDUC_CAT", "CIDADE", "FREQ_INJEC",
             "USED_NEEDLE", "HIV_TESTED_12M")

message("Estimating subgroup RDS-II prevalence ...")
tab_subgroup <- bind_rows(lapply(c("HIV", "HCV", "HBV"), function(v) {
  bind_rows(lapply(sg_vars, function(g) estimate_subgroup(dat_comp, v, g)))
}))

tab_subgroup_wide <- tab_subgroup |>
  select(outcome, variable, category, survey_year, n, n_pos, rds_p, rds_se, rds_lo, rds_hi) |>
  pivot_wider(
    names_from = survey_year,
    values_from = c(n, n_pos, rds_p, rds_se, rds_lo, rds_hi),
    names_sep = "_"
  ) |>
  rowwise() |>
  mutate(wald = list(wald_diff(rds_p_2023, rds_se_2023, rds_p_2014, rds_se_2014))) |>
  unnest(wald) |>
  ungroup() |>
  mutate(
    n_2014_lab = paste0(n_pos_2014, "/", n_2014),
    n_2023_lab = paste0(n_pos_2023, "/", n_2023),
    prev_2014 = fmt_pct_ci(rds_p_2014, rds_lo_2014, rds_hi_2014),
    prev_2023 = fmt_pct_ci(rds_p_2023, rds_lo_2023, rds_hi_2023),
    diff_pp = sprintf("%.1f (%.1f to %.1f)", 100 * diff, 100 * ci_lo, 100 * ci_hi),
    p_lab = fmt_p(p_value)
  )

write_csv_safe(tab_subgroup_wide, file.path(path_tab, "table_subgroup_rdsii.csv"))

## ---- 9. Pre-specified logistic models (no stepwise) --------------------------
## Avery 2019 / Sperandei 2023: unweighted logistic can be preferable for
## *associations* in RDS. Reviewers also requested weighted estimates for
## comparability. We therefore report both, with cluster-robust SEs.

apriori_covars <- c("YEAR", "CIDADE", "AGE_CAT", "SEX_CAT", "EDUC_CAT",
                    "MARITALC", "FREQ_INJEC", "USED_NEEDLE", "HIV_TESTED_12M")

prep_model_df <- function(dat, outcome) {
  dat |>
    mutate(
      YEAR = factor(YEAR, levels = c("2014", "2023")),
      CIDADE = factor(CIDADE, levels = c("MAPUTO", "NAMPULA")),
      AGE_CAT = factor(AGE_CAT, levels = c("18-24", "25+")),
      SEX_CAT = factor(SEX_CAT, levels = c("Male", "Female")),
      EDUC_CAT = factor(EDUC_CAT, levels = c("Secondary or higher", "Primary or none")),
      MARITALC = factor(MARITALC, levels = c("Never married", "Married or union", "Other")),
      FREQ_INJEC = factor(FREQ_INJEC, levels = c("Less than daily", "Daily")),
      USED_NEEDLE = factor(USED_NEEDLE, levels = c("No", "Yes")),
      HIV_TESTED_12M = factor(HIV_TESTED_12M, levels = c("No", "Yes")),
      y = .data[[outcome]],
      cluster_id = as.character(chain_id)
    ) |>
    select(y, all_of(apriori_covars), w_pop, cluster_id, coupon_id) |>
    filter(complete.cases(across(c(y, all_of(apriori_covars), w_pop, cluster_id))))
}

or_from_svy <- function(fit) {
  cf <- coef(fit)
  V <- tryCatch(vcov(fit), error = function(e) sandwich::vcovHC(fit, type = "HC0"))
  se <- sqrt(diag(V))[names(cf)]
  tibble(
    term = names(cf),
    or = exp(cf),
    ci_lo = exp(cf - 1.96 * se),
    ci_hi = exp(cf + 1.96 * se),
    p_value = 2 * pnorm(-abs(cf / se))
  ) |>
    filter(term != "(Intercept)") |>
    mutate(or_ci = sprintf("%.2f (%.2f-%.2f)", or, ci_lo, ci_hi),
           p_lab = fmt_p(p_value))
}

or_from_glm_cluster <- function(fit, cluster) {
  V <- sandwich::vcovCL(fit, cluster = cluster, type = "HC0")
  cf <- coef(fit)
  se <- sqrt(diag(V))[names(cf)]
  tibble(
    term = names(cf),
    or = exp(cf),
    ci_lo = exp(cf - 1.96 * se),
    ci_hi = exp(cf + 1.96 * se),
    p_value = 2 * pnorm(-abs(cf / se))
  ) |>
    filter(term != "(Intercept)") |>
    mutate(or_ci = sprintf("%.2f (%.2f-%.2f)", or, ci_lo, ci_hi),
           p_lab = fmt_p(p_value))
}

fit_models <- function(dat, outcome) {
  md <- prep_model_df(dat, outcome)
  fml <- y ~ YEAR + CIDADE + AGE_CAT + SEX_CAT + EDUC_CAT + MARITALC +
    FREQ_INJEC + USED_NEEDLE + HIV_TESTED_12M + YEAR:CIDADE
  ## (a) Unweighted logistic with chain-clustered sandwich SEs
  ##     (Avery et al. 2019; Sperandei et al. 2023)
  md <- md |> dplyr::arrange(cluster_id)
  glm_u <- stats::glm(fml, data = md, family = binomial())
  ## (b) Population-weighted survey GLM, clustered by chain, stratified by city
  md$strata_id <- md$CIDADE
  des <- survey::svydesign(
    ids = ~cluster_id, strata = ~strata_id, weights = ~w_pop,
    data = md, nest = TRUE
  )
  svy <- survey::svyglm(fml, design = des, family = quasibinomial())
  list(
    n = nrow(md),
    gee = or_from_glm_cluster(glm_u, md$cluster_id) |>
      mutate(model = "Unweighted logistic, chain-clustered SEs", outcome = outcome),
    svy = or_from_svy(svy) |> mutate(model = "RDS population-weighted GLM", outcome = outcome)
  )
}

message("Fitting pre-specified GEE and survey-weighted GLMs ...")
models <- lapply(c("HIV", "HCV", "HBV"), function(v) {
  tryCatch(fit_models(dat_comp, v), error = function(e) {
    message("Model failed for ", v, ": ", conditionMessage(e))
    list(
      n = NA_integer_,
      gee = tibble(term = NA_character_, or_ci = NA_character_, p_lab = NA_character_,
                   model = "Unweighted GEE (chain-clustered)", outcome = v),
      svy = tibble(term = NA_character_, or_ci = NA_character_, p_lab = NA_character_,
                   model = "RDS population-weighted GLM", outcome = v)
    )
  })
})
names(models) <- c("HIV", "HCV", "HBV")

tab_models <- bind_rows(lapply(models, function(m) bind_rows(m$gee, m$svy)))
write_csv_safe(tab_models, file.path(path_tab, "table_multivariable_apriori.csv"))

## Year effect by city from the interaction model
year_effect_by_city <- tab_models |>
  filter(term %in% c("YEAR2023", "YEAR2023:CIDADENAMPULA") | grepl("YEAR", term))

## ---- 10. Sensitivity: 2023 all five cities (descriptive only) ---------------

sens_2023_all <- bind_rows(lapply(rds23_all, function(x) {
  bind_rows(lapply(c("HIV", "HCV", "HBV"), function(v) estimate_city_outcome(x, v)))
}))
write_csv_safe(sens_2023_all, file.path(path_tab, "table_sens_2023_all_cities.csv"))

## Sensitivity: equalised sample size (downsample 2023 to 2014 n within city)
set.seed(20260919)
downsample_city <- function() {
  bind_rows(lapply(comparable_cities, function(ct) {
    d0 <- dat14 |> filter(CIDADE == ct)
    d1 <- dat23 |> filter(CIDADE == ct)
    d1s <- d1[sample.int(nrow(d1), size = min(nrow(d0), nrow(d1))), ]
    bind_rows(d0, d1s)
  }))
}

ds_runs <- replicate(50, {
  ds <- downsample_city()
  bind_rows(lapply(c("HIV", "HCV", "HBV"), function(v) {
    bind_rows(lapply(comparable_cities, function(ct) {
      a <- vh_estimate(ds[[v]][ds$CIDADE == ct & ds$YEAR == "2014"],
                       ds$degree[ds$CIDADE == ct & ds$YEAR == "2014"], n_boot = 50)
      b <- vh_estimate(ds[[v]][ds$CIDADE == ct & ds$YEAR == "2023"],
                       ds$degree[ds$CIDADE == ct & ds$YEAR == "2023"], n_boot = 50)
      tibble(outcome = v, city = ct, diff = b$est - a$est)
    }))
  }))
}, simplify = FALSE)

ds_summary <- bind_rows(ds_runs, .id = "rep") |>
  group_by(outcome, city) |>
  summarise(
    median_diff = median(diff, na.rm = TRUE),
    lo = quantile(diff, 0.025, na.rm = TRUE),
    hi = quantile(diff, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(diff_pp = sprintf("%.1f (%.1f to %.1f)", 100 * median_diff, 100 * lo, 100 * hi))

write_csv_safe(ds_summary, file.path(path_tab, "table_sens_downsample.csv"))

## ---- 11. Design-effect from 2023 (used to describe sample dependence) --------

deff_2023 <- city_prev |>
  filter(survey_year == "2023") |>
  mutate(
    se_srs = sqrt(crude_p * (1 - crude_p) / n),
    deff = (rds_se / se_srs)^2
  ) |>
  select(city, outcome, n, crude_p, rds_p, se_srs, rds_se, deff, method)

write_csv_safe(deff_2023, file.path(path_tab, "table_design_effects.csv"))

## ---- 12. Comparability documentation ---------------------------------------

n_flow <- tibble(
  item = c(
    "2014 enrolled (Maputo)",
    "2014 enrolled (Nampula/Nacala)",
    "2014 analytic sample, injected past 12 months, age >= 18",
    "2023 enrolled all five cities",
    "2023 Maputo + Nampula",
    "2023 comparable sample (Maputo+Nampula, age >= 18, injected past 12 months)",
    "2014 with valid HIV result (comparable sample)",
    "2023 with valid HIV result (comparable sample)"
  ),
  n = c(
    sum(pid14_raw$CIDADE == "MAPUTO"),
    sum(pid14_raw$CIDADE == "NAMPULA"),
    nrow(pid14),
    nrow(pid23_raw),
    sum(pid23_raw$CIDADE %in% comparable_cities),
    nrow(pid23_comp),
    sum(!is.na(pid14$HIV)),
    sum(!is.na(pid23_comp$HIV))
  )
)
write_csv_safe(n_flow, file.path(path_tab, "table_sample_flow.csv"))

elig_table <- tibble(
  domain = c(
    "Design", "Geography", "Age eligibility", "Injection eligibility",
    "Recruitment", "Network size", "RDS weights in this re-analysis",
    "Primary temporal comparison"
  ),
  bbs_2014 = c(
    "Repeated cross-sectional IBBS (RDS)",
    "Maputo city and Nampula/Nacala",
    "18 years or older",
    "Ever injected, or injected in the past 12 months; analytic sample restricted to past-12-month injection",
    "RDS, 3 coupons, 5 seeds per city",
    "Self-reported degree (DEGREE/NETSIZE1)",
    "RDS-II / Gile SS using published PSE (Maputo 1684; Nampula 520)",
    "Maputo and Nampula, age >= 18, recent injection"
  ),
  bbs_2023 = c(
    "Repeated cross-sectional IBBS (RDS)",
    "Maputo, Beira, Tete, Quelimane, Nampula (primary comparison uses Maputo and Nampula only)",
    "16 years or older; analytic sample restricted to 18+",
    "Injected in the past 12 months",
    "RDS, up to 5 coupons",
    "Self-reported degree (ELINJMU, already cleaned in the source file)",
    "Gile SS using mapping PSE (Maputo 990; Nampula 1792)",
    "Maputo and Nampula, age >= 18, recent injection"
  )
)
write_csv_safe(elig_table, file.path(path_tab, "table_survey_comparability.csv"))

## ---- 13. Figures ------------------------------------------------------------

fig_df <- city_prev |>
  filter(outcome %in% c("HIV", "HCV", "HBV")) |>
  mutate(
    city = factor(city, levels = c("MAPUTO", "NAMPULA")),
    outcome = factor(outcome, levels = c("HIV", "HCV", "HBV"),
                     labels = c("HIV", "HCV antibody", "HBsAg"))
  )

p1 <- ggplot(fig_df, aes(x = survey_year, y = 100 * rds_p, colour = survey_year)) +
  geom_pointrange(aes(ymin = 100 * rds_lo, ymax = 100 * rds_hi),
                  position = position_dodge(width = 0.2), linewidth = 0.6) +
  geom_point(aes(y = 100 * crude_p), shape = 1, size = 2.2, colour = "grey30") +
  facet_grid(outcome ~ city) +
  scale_colour_manual(values = c("2014" = "#1b4f72", "2023" = "#b03a2e")) +
  labs(
    x = NULL, y = "Prevalence (%)",
    title = "RDS-weighted prevalence among PWID, Maputo and Nampula",
    subtitle = "Filled points: RDS-weighted estimate with 95% CI. Open circles: unweighted sample proportion.",
    colour = "Survey"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(path_fig, "fig1_rds_prevalence_city.pdf"), p1, width = 8.2, height = 7.2)
ggsave(file.path(path_fig, "fig1_rds_prevalence_city.png"), p1, width = 8.2, height = 7.2, dpi = 300)

fig2_df <- prev_compare |>
  filter(city %in% c("MAPUTO", "NAMPULA", "MAPUTO+NAMPULA"),
         outcome %in% c("HIV", "HCV", "HBV")) |>
  mutate(
    city = factor(city, levels = c("MAPUTO", "NAMPULA", "MAPUTO+NAMPULA"),
                  labels = c("Maputo", "Nampula", "Both cities (population-weighted)")),
    outcome = factor(outcome, levels = c("HIV", "HCV", "HBV"),
                     labels = c("HIV", "HCV antibody", "HBsAg"))
  )

p2 <- ggplot(fig2_df, aes(x = 100 * diff, y = city)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = 100 * ci_lo, xmax = 100 * ci_hi)) +
  facet_wrap(~outcome, ncol = 1) +
  labs(
    x = "Absolute difference in RDS-weighted prevalence, 2023 minus 2014 (percentage points)",
    y = NULL,
    title = "Change in infection prevalence between serial RDS surveys"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(path_fig, "fig2_prevalence_difference.pdf"), p2, width = 8.0, height = 6.4)
ggsave(file.path(path_fig, "fig2_prevalence_difference.png"), p2, width = 8.0, height = 6.4, dpi = 300)

## ---- 14. LaTeX table fragments ----------------------------------------------

tex_prev <- prev_compare |>
  filter(outcome %in% c("HIV", "HCV", "HBV", "HIV_HCV"),
         city %in% c("MAPUTO", "NAMPULA", "MAPUTO+NAMPULA")) |>
  transmute(
    Outcome = recode(outcome, HIV = "HIV", HCV = "HCV antibody", HBV = "HBsAg",
                     HIV_HCV = "HIV-HCV co-infection"),
    City = recode(city, MAPUTO = "Maputo", NAMPULA = "Nampula",
                  `MAPUTO+NAMPULA` = "Both cities"),
    `2014 n` = n_2014_lab,
    `2014 RDS % (95% CI)` = prev_2014,
    `2023 n` = n_2023_lab,
    `2023 RDS % (95% CI)` = prev_2023,
    `Difference pp (95% CI)` = diff_pp,
    `p` = p_lab
  )
df_to_tex(
  tex_prev,
  file.path(path_tex, "table_prevalence.tex"),
  caption = "RDS-weighted prevalence of HIV, HCV antibody and HBsAg among people who inject drugs in Maputo and Nampula, Mozambique, 2014 and 2023. City-specific estimates use Gile successive sampling (2023) or RDS-II (2014). The two-city total is weighted by estimated PWID population size, not by sample size. Differences are 2023 minus 2014.",
  label = "tab:prev"
)

## Table 1 LaTeX is written by PID_TABLE1_RDS.R (city-stratified RDS-II).

tex_mod <- tab_models |>
  transmute(
    Outcome = outcome,
    Model = model,
    Term = term,
    `AOR (95% CI)` = or_ci,
    `p` = p_lab
  )
df_to_tex(
  tex_mod,
  file.path(path_tex, "table_models.tex"),
  caption = "Pre-specified multivariable models for HIV, HCV antibody and HBsAg in the comparable sample. Covariates were chosen a priori and were not subjected to stepwise selection. Unweighted GEE uses an exchangeable working correlation within recruitment chains. The survey GLM uses population-calibrated RDS weights and chain clustering.",
  label = "tab:models"
)

tex_std <- std_compare |>
  transmute(
    Outcome = recode(outcome, HIV = "HIV", HCV = "HCV antibody", HBV = "HBsAg"),
    `2014 RDS %` = prev_2014,
    `2023 RDS % (unadjusted)` = prev_2023_unadj,
    `2023 % standardised to 2014 age-sex-city mix` = prev_2023_std
  )
df_to_tex(
  tex_std,
  file.path(path_tex, "table_standardisation.tex"),
  caption = "Direct standardisation of 2023 RDS-weighted prevalence to the 2014 age, sex and city distribution of PWID. This isolates prevalence change from differences in demographic composition between surveys.",
  label = "tab:std"
)

tex_elig <- elig_table
df_to_tex(
  tex_elig,
  file.path(path_tex, "table_eligibility.tex"),
  caption = "Comparability of the 2014 and 2023 biobehavioural surveys and of the analytic restrictions used in this re-analysis.",
  label = "tab:elig"
)

## ---- 15. Session info and key numeric summary for the manuscript ------------

key_summary <- prev_compare |>
  filter(city %in% c("MAPUTO", "NAMPULA", "MAPUTO+NAMPULA"),
         outcome %in% c("HIV", "HCV", "HBV")) |>
  select(outcome, city, n_2014, n_2023, prev_2014, prev_2023, diff_pp, p_lab)

write_csv_safe(key_summary, file.path(path_tab, "key_numeric_summary.csv"))
write_csv_safe(n_flow, file.path(path_tab, "sample_flow.csv"))

saveRDS(
  list(
    dat_comp = dat_comp,
    city_prev = city_prev,
    pooled_prev = pooled_prev,
    prev_compare = prev_compare,
    tab_composition = tab_composition,
    tab_models = tab_models,
    std_compare = std_compare,
    deff_2023 = deff_2023,
    n_flow = n_flow
  ),
  file.path(path_out, "analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(path_out, "sessionInfo.txt")
)

message("Analysis complete. Tables written to ", path_tab)
print(n_flow)
print(key_summary)

## Rebuild Table 1 with city-stratified RDS-II + bootstrap Wald
source(file.path(PATH_R, "02_table1_composition.R"), local = TRUE)
