################################################################################
## Table 2. Population-weighted HIV and HCV prevalence (Gile SS + SS-PSE N)
##
## Sample: Maputo and Nampula, age >= 18, injected in the previous 12 months.
## City PWID size N_c from SS-PSE (sspse_city_fits.rds).
## Weights: RDS::compute.weights(..., weight.type = "Gile's SS", N = N_c),
##          scaled so they sum to N_c. Pooled estimate is
##              p = sum(w * y) / sum(w)
## Difference = 2023 minus 2014. Wald p from independent-year bootstrap SEs
## (within-city resampling; Gile weights kept as sampling weights).
################################################################################

options(stringsAsFactors = FALSE)
set.seed(20260919)

## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))


library(RDS)
library(openxlsx)
library(officer)
library(flextable)

N_BOOT <- 300

obj <- readRDS(file.path(path_out, "analysis_objects.rds"))
dat <- obj$dat_comp
sspse_fits <- readRDS(file.path(path_out, "sspse_city_fits.rds"))

N_hat <- list(
  "2014" = c(MAPUTO = sspse_fits[["2014_MAPUTO"]]$N_use,
             NAMPULA = sspse_fits[["2014_NAMPULA"]]$N_use),
  "2023" = c(MAPUTO = sspse_fits[["2023_MAPUTO"]]$N_use,
             NAMPULA = sspse_fits[["2023_NAMPULA"]]$N_use)
)
print(N_hat)

## ---- extra items from the questionnaires (Table 2) --------------------------

d14 <- read.csv(path_2014,
                stringsAsFactors = FALSE, encoding = "latin1",
                na.strings = c("", "NA", "NaN", "nan"))
d23 <- read.csv(path_2023,
                stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN", "nan"))

extra14 <- data.frame(
  YEAR = "2014",
  coupon_id = as.character(d14$COUPON),
  fv = suppressWarnings(as.numeric(d14$FVSEX2)),
  ms = suppressWarnings(as.numeric(d14$MSEX2)),
  share = d14$DRINJ3C,
  risk = d14$VCTRISK1C,
  stringsAsFactors = FALSE
)
extra14$SEX_PARTNER <- NA_character_
extra14$SEX_PARTNER[extra14$fv %in% 0:1 | extra14$ms %in% 0:1] <- "0-1"
extra14$SEX_PARTNER[extra14$fv == 2 | extra14$ms == 2] <- "2"
extra14$SEX_PARTNER[(extra14$fv > 2 & !extra14$fv %in% c(9997, 9999)) |
                      (extra14$ms > 2 & !extra14$ms %in% c(9997, 9999))] <- ">=3"
extra14$SHARE_SYRINGE <- NA_character_
extra14$SHARE_SYRINGE[extra14$share == 1] <- "Yes"
extra14$SHARE_SYRINGE[extra14$share == 0] <- "No"
extra14$RISK_PERCEPTION <- NA_character_
extra14$RISK_PERCEPTION[extra14$risk == "1_NONELOW"] <- "No/low risk"
extra14$RISK_PERCEPTION[extra14$risk == "2_MODHIGH"] <- "Moderate/high risk"
extra14$RISK_PERCEPTION[extra14$risk == "3_REFUSED"] <- "Refused"

extra23 <- data.frame(
  YEAR = "2023",
  coupon_id = as.character(d23$id),
  ident42 = suppressWarnings(as.numeric(d23$IDENT42_LIMFPART_A)),
  ident43 = suppressWarnings(as.numeric(d23$IDENT43_MSEX2_A)),
  limh = d23$LIMFVAG_H,
  limm = d23$LIMFVAG_M,
  drinj3 = d23$DRINJ3,
  idshare = d23$IDSHARE,
  idrel1 = d23$IDREL1,
  risk = d23$VCTRISK1_CAT,
  stringsAsFactors = FALSE
)
extra23$SEX_PARTNER <- NA_character_
extra23$SEX_PARTNER[extra23$ident42 == 0 | extra23$limh == "Nao" |
                      extra23$limm == "Nao" | extra23$ident42 == 1 |
                      extra23$ident43 == 1] <- "0-1"
extra23$SEX_PARTNER[extra23$ident42 == 2 | extra23$ident43 == 2] <- "2"
extra23$SEX_PARTNER[extra23$ident42 > 2 | extra23$ident43 > 2] <- ">=3"
extra23$SHARE_SYRINGE <- NA_character_
extra23$SHARE_SYRINGE[extra23$drinj3 == "Sim" | extra23$idshare == "Sim" |
                        extra23$idrel1 == "Sim"] <- "Yes"
extra23$SHARE_SYRINGE[is.na(extra23$SHARE_SYRINGE) &
                        (extra23$drinj3 == "Nao" | extra23$idshare == "Nao")] <- "No"
extra23$RISK_PERCEPTION <- NA_character_
extra23$RISK_PERCEPTION[extra23$risk == "1_SEM_RISCO_POUCO"] <- "No/low risk"
extra23$RISK_PERCEPTION[extra23$risk == "2_RISCO_MODERADO_ALTO"] <- "Moderate/high risk"
extra23$RISK_PERCEPTION[extra23$risk == "3_SEM_RESPOSTA"] <- "Refused"

extra <- rbind(
  extra14[, c("YEAR", "coupon_id", "SEX_PARTNER", "SHARE_SYRINGE", "RISK_PERCEPTION")],
  extra23[, c("YEAR", "coupon_id", "SEX_PARTNER", "SHARE_SYRINGE", "RISK_PERCEPTION")]
)
extra <- extra[!duplicated(extra[, c("YEAR", "coupon_id")]), ]
dat$SEX_PARTNER <- NULL
dat$SHARE_SYRINGE <- NULL
dat$RISK_PERCEPTION <- NULL
dat <- merge(dat, extra, by = c("YEAR", "coupon_id"), all.x = TRUE, sort = FALSE)

dat$STI_SELF[dat$STI_SELF %in% c("1_SIM", "Yes")] <- "Yes"
dat$STI_SELF[dat$STI_SELF %in% c("2_NAO", "2_SIM", "No")] <- "No"
dat$FREQ_INJEC[dat$FREQ_INJEC == "Less than daily"] <- "Weekly/monthly"
dat$EDUC_CAT[dat$EDUC_CAT == "Primary or none"] <- "No formal/Primary"
dat$EDUC_CAT[dat$EDUC_CAT == "Secondary or higher"] <- "Secondary/Higher"
dat$MARITALC[dat$MARITALC == "Married or union"] <- "Married/Living in union"
dat$AGE_CAT[dat$AGE_CAT == "25+"] <- ">=25"
dat$AGE_FIRST_DRUG[dat$AGE_FIRST_DRUG == "25+"] <- ">=25"

stopifnot(all(dat$CIDADE %in% c("MAPUTO", "NAMPULA")))
stopifnot(min(dat$AGE, na.rm = TRUE) >= 18)

## ---- Gile SS weights, scaled to SS-PSE N_c ---------------------------------

dat$w <- NA_real_
for (yr in c("2014", "2023")) {
  for (ct in c("MAPUTO", "NAMPULA")) {
    ii <- which(dat$YEAR == yr & dat$CIDADE == ct)
    d <- dat[ii, ]
    N <- unname(N_hat[[yr]][ct])
    tmp <- d
    tmp$.id <- as.character(tmp$coupon_id)
    tmp$.rec <- as.character(tmp$recruiter_id)
    tmp$.rec[is.na(tmp$.rec) | tmp$.rec %in% c("", "NA", "seed", "SEED")] <- "0"
    tmp$.rec[!tmp$.rec %in% c("0", tmp$.id)] <- "0"
    tmp$.rec[tmp$.rec == tmp$.id] <- "0"
    tmp$.deg <- pmax(1, as.numeric(tmp$degree))
    tmp <- tmp[!duplicated(tmp$.id), ]
    rds <- RDS::as.rds.data.frame(
      tmp, id = ".id", recruiter.id = ".rec",
      network.size = ".deg", population.size = N, max.coupons = 5
    )
    w <- as.numeric(RDS::compute.weights(rds, weight.type = "Gile's SS", N = N))
    names(w) <- as.character(rds[[".id"]])
    w_out <- unname(w[as.character(d$coupon_id)])
    w_out <- w_out / sum(w_out) * N
    dat$w[ii] <- w_out
    message(yr, " ", ct, ": n=", length(ii), " N=", N,
            " sum(w)=", round(sum(w_out)))
  }
}

## ---- bootstrap samples (within city, by year) -------------------------------

boot <- list("2014" = vector("list", N_BOOT), "2023" = vector("list", N_BOOT))
for (yr in c("2014", "2023")) {
  d_yr <- dat[dat$YEAR == yr, ]
  for (b in seq_len(N_BOOT)) {
    b_map <- d_yr[d_yr$CIDADE == "MAPUTO", ]
    b_nam <- d_yr[d_yr$CIDADE == "NAMPULA", ]
    boot[[yr]][[b]] <- rbind(
      b_map[sample.int(nrow(b_map), replace = TRUE), ],
      b_nam[sample.int(nrow(b_nam), replace = TRUE), ]
    )
  }
}

## ---- Table 2 rows (same order as the manuscript) ----------------------------

rows_spec <- rbind(
  data.frame(domain = "Sociodemographic characteristics",
             variable = "AGE_CAT", category = "18-24",
             label = "Age group: 18-24"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "AGE_CAT", category = ">=25",
             label = "Age group: >=25"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "SEX_CAT", category = "Male",
             label = "Sex: Male"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "SEX_CAT", category = "Female",
             label = "Sex: Female"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "EDUC_CAT", category = "No formal/Primary",
             label = "Education: No formal/Primary"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "EDUC_CAT", category = "Secondary/Higher",
             label = "Education: Secondary/Higher"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "MARITALC", category = "Never married",
             label = "Marital status: Never married"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "MARITALC", category = "Married/Living in union",
             label = "Marital status: Married/Living in union"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "MARITALC", category = "Other",
             label = "Marital status: Other"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "CIDADE", category = "MAPUTO",
             label = "City of residence: Maputo"),
  data.frame(domain = "Sociodemographic characteristics",
             variable = "CIDADE", category = "NAMPULA",
             label = "City of residence: Nampula/Nacala"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "AGE_FIRST_DRUG", category = "<18",
             label = "Age at first drug use: <18"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "AGE_FIRST_DRUG", category = "18-24",
             label = "Age at first drug use: 18-24"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "AGE_FIRST_DRUG", category = ">=25",
             label = "Age at first drug use: >=25"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "FREQ_INJEC", category = "Daily",
             label = "Injection frequency: Daily"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "FREQ_INJEC", category = "Weekly/monthly",
             label = "Injection frequency: Weekly/monthly"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "NEW_SYRINGE", category = "Yes",
             label = "Access to new syringes: Yes"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "NEW_SYRINGE", category = "No",
             label = "Access to new syringes: No"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "SHARE_SYRINGE", category = "Yes",
             label = "Shared syringes (past 12 months): Yes"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "SHARE_SYRINGE", category = "No",
             label = "Shared syringes (past 12 months): No"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "USED_NEEDLE", category = "Yes",
             label = "Ever injected with previously used syringe: Yes"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "USED_NEEDLE", category = "No",
             label = "Ever injected with previously used syringe: No"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "SHARED_EQUIP", category = "Yes",
             label = "Ever shared any other injection equipment: Yes"),
  data.frame(domain = "Drug use and injection-related characteristics",
             variable = "SHARED_EQUIP", category = "No",
             label = "Ever shared any other injection equipment: No"),
  data.frame(domain = "Sexual behaviours",
             variable = "SEX_PARTNER", category = "0-1",
             label = "Number of sexual partners (past 12 months): 0-1"),
  data.frame(domain = "Sexual behaviours",
             variable = "SEX_PARTNER", category = "2",
             label = "Number of sexual partners (past 12 months): 2"),
  data.frame(domain = "Sexual behaviours",
             variable = "SEX_PARTNER", category = ">=3",
             label = "Number of sexual partners (past 12 months): >=3"),
  data.frame(domain = "Sexual behaviours",
             variable = "UNPROTECTED_SEX", category = "Yes",
             label = "Unprotected sex (past 12 months): Yes"),
  data.frame(domain = "Sexual behaviours",
             variable = "UNPROTECTED_SEX", category = "No",
             label = "Unprotected sex (past 12 months): No"),
  data.frame(domain = "Sexual behaviours",
             variable = "TX_SEX", category = "Yes",
             label = "Received money, goods, or services for sex: Yes"),
  data.frame(domain = "Sexual behaviours",
             variable = "TX_SEX", category = "No",
             label = "Received money, goods, or services for sex: No"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "STI_SELF", category = "Yes",
             label = "Self-reported STI: Yes"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "STI_SELF", category = "No",
             label = "Self-reported STI: No"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "HIV_TESTED_12M", category = "Yes",
             label = "HIV test in the last 12 months: Yes"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "HIV_TESTED_12M", category = "No",
             label = "HIV test in the last 12 months: No"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "RISK_PERCEPTION", category = "No/low risk",
             label = "HIV risk perception: No/low risk"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "RISK_PERCEPTION", category = "Moderate/high risk",
             label = "HIV risk perception: Moderate/high risk"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "RISK_PERCEPTION", category = "Refused",
             label = "HIV risk perception: Refused"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "SOUGHT_CARE", category = "Yes",
             label = "Sought healthcare (past 12 months): Yes"),
  data.frame(domain = "Healthcare use and HIV risk perception",
             variable = "SOUGHT_CARE", category = "No",
             label = "Sought healthcare (past 12 months): No")
)

fmt_pct <- function(p, lo, hi) {
  if (is.na(p)) return(NA_character_)
  sprintf("%.1f (%.1f-%.1f)", 100 * p, 100 * lo, 100 * hi)
}
fmt_diff <- function(d, lo, hi) {
  if (is.na(d)) return(NA_character_)
  sprintf("%.1f (%.1f to %.1f)", 100 * d, 100 * lo, 100 * hi)
}
fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

out_hiv <- vector("list", nrow(rows_spec))
out_hcv <- vector("list", nrow(rows_spec))

for (i in seq_len(nrow(rows_spec))) {
  v <- rows_spec$variable[i]
  lv <- rows_spec$category[i]
  message("Row ", i, "/", nrow(rows_spec), ": ", v, " = ", lv)

  for (outc in c("HIV", "HCV")) {
    rec <- list()
    rec$domain <- rows_spec$domain[i]
    rec$characteristic <- rows_spec$label[i]
    rec$variable <- v
    rec$category <- lv
    rec$outcome <- outc

    for (yr in c("2014", "2023")) {
      d_yr <- dat[dat$YEAR == yr, ]
      ok <- !is.na(d_yr[[v]]) & d_yr[[v]] == lv &
        !is.na(d_yr[[outc]]) & !is.na(d_yr$w) & d_yr$w > 0
      n <- sum(ok)
      n_pos <- sum(d_yr[[outc]][ok] == 1)
      p <- if (n >= 8) sum(d_yr$w[ok] * d_yr[[outc]][ok]) / sum(d_yr$w[ok]) else NA_real_
      rec[[paste0("n_", yr)]] <- n
      rec[[paste0("n_pos_", yr)]] <- n_pos
      rec[[paste0("p_", yr)]] <- p

      if (is.na(p)) {
        rec[[paste0("se_", yr)]] <- NA_real_
        rec[[paste0("lo_", yr)]] <- NA_real_
        rec[[paste0("hi_", yr)]] <- NA_real_
      } else {
        bb <- numeric(N_BOOT)
        for (b in seq_len(N_BOOT)) {
          db <- boot[[yr]][[b]]
          okb <- !is.na(db[[v]]) & db[[v]] == lv &
            !is.na(db[[outc]]) & !is.na(db$w) & db$w > 0
          if (sum(okb) < 8) {
            bb[b] <- NA_real_
          } else {
            bb[b] <- sum(db$w[okb] * db[[outc]][okb]) / sum(db$w[okb])
          }
        }
        rec[[paste0("se_", yr)]] <- sd(bb, na.rm = TRUE)
        q <- quantile(bb, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
        rec[[paste0("lo_", yr)]] <- max(0, q[1])
        rec[[paste0("hi_", yr)]] <- min(1, q[2])
      }
    }

    if (!is.na(rec$p_2014) && !is.na(rec$p_2023)) {
      rec$diff <- rec$p_2023 - rec$p_2014
      rec$se_diff <- sqrt(rec$se_2014^2 + rec$se_2023^2)
      rec$diff_lo <- rec$diff - 1.96 * rec$se_diff
      rec$diff_hi <- rec$diff + 1.96 * rec$se_diff
      rec$p_value <- 2 * pnorm(-abs(rec$diff / rec$se_diff))
    } else {
      rec$diff <- rec$se_diff <- rec$diff_lo <- rec$diff_hi <- rec$p_value <- NA_real_
    }

    rec$n_2014_lab <- paste0(rec$n_pos_2014, "/", rec$n_2014)
    rec$n_2023_lab <- paste0(rec$n_pos_2023, "/", rec$n_2023)
    rec$prev_2014 <- fmt_pct(rec$p_2014, rec$lo_2014, rec$hi_2014)
    rec$prev_2023 <- fmt_pct(rec$p_2023, rec$lo_2023, rec$hi_2023)
    rec$diff_lab <- fmt_diff(rec$diff, rec$diff_lo, rec$diff_hi)
    rec$p_lab <- fmt_p(rec$p_value)

    if (outc == "HIV") out_hiv[[i]] <- rec else out_hcv[[i]] <- rec
  }
}

hiv <- as.data.frame(do.call(rbind, lapply(out_hiv, as.data.frame)),
                     stringsAsFactors = FALSE)
hcv <- as.data.frame(do.call(rbind, lapply(out_hcv, as.data.frame)),
                     stringsAsFactors = FALSE)

make_wide <- function(d) {
  data.frame(
    Domain = d$domain,
    Characteristic = d$characteristic,
    `2014 n` = d$n_2014_lab,
    `2014 % (95% CI)` = d$prev_2014,
    `2023 n` = d$n_2023_lab,
    `2023 % (95% CI)` = d$prev_2023,
    `Difference pp (95% CI)` = d$diff_lab,
    `p` = d$p_lab,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

tab_hiv <- make_wide(hiv)
tab_hcv <- make_wide(hcv)

write.csv(tab_hiv, file.path(path_tab, "table2_hiv_gile_sspse.csv"), row.names = FALSE)
write.csv(tab_hcv, file.path(path_tab, "table2_hcv_gile_sspse.csv"), row.names = FALSE)
write.csv(cbind(hiv, outcome = "HIV"),
          file.path(path_tab, "table2_hiv_gile_sspse_long.csv"), row.names = FALSE)
write.csv(cbind(hcv, outcome = "HCV"),
          file.path(path_tab, "table2_hcv_gile_sspse_long.csv"), row.names = FALSE)

## Formatted Excel/Word with merged domain header rows
source(file.path(PATH_R, "03b_table2_format.R"), local = TRUE)

