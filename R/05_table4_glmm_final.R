################################################################################
## Final GLMM after LRT selection
## - Force YEAR x CIDADE into the HIV model (design / sampling contrast)
## - LRT test of recruiter random intercept vs ordinary logistic (both outcomes)
## - Publication-ready Excel table
################################################################################

options(stringsAsFactors = FALSE)
set.seed(20260919)

## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))


library(lme4)
library(openxlsx)

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_or <- function(or, lo, hi) {
  if (is.na(or)) return(NA_character_)
  sprintf("%.2f (%.2f\u2013%.2f)", or, lo, hi)
}

obj <- readRDS(file.path(path_out, "analysis_objects.rds"))
dat <- obj$dat_comp

d14 <- read.csv(path_2014,
                stringsAsFactors = FALSE, encoding = "latin1",
                na.strings = c("", "NA", "NaN", "nan"))
d23 <- read.csv(path_2023,
                stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN", "nan"))
p14 <- data.frame(YEAR = "2014", coupon_id = as.character(d14$COUPON),
                  fv = suppressWarnings(as.numeric(d14$FVSEX2)),
                  ms = suppressWarnings(as.numeric(d14$MSEX2)),
                  share = d14$DRINJ3C, stringsAsFactors = FALSE)
p14$SEX_PARTNER <- NA_character_
p14$SEX_PARTNER[p14$fv %in% 0:1 | p14$ms %in% 0:1] <- "0-1"
p14$SEX_PARTNER[p14$fv == 2 | p14$ms == 2] <- "2"
p14$SEX_PARTNER[(p14$fv > 2 & !p14$fv %in% c(9997, 9999)) |
                  (p14$ms > 2 & !p14$ms %in% c(9997, 9999))] <- ">=3"
p14$SHARE_SYRINGE <- NA_character_
p14$SHARE_SYRINGE[p14$share == 1] <- "Yes"
p14$SHARE_SYRINGE[p14$share == 0] <- "No"
p23 <- data.frame(YEAR = "2023", coupon_id = as.character(d23$id),
                  ident42 = suppressWarnings(as.numeric(d23$IDENT42_LIMFPART_A)),
                  ident43 = suppressWarnings(as.numeric(d23$IDENT43_MSEX2_A)),
                  limh = d23$LIMFVAG_H, limm = d23$LIMFVAG_M,
                  drinj3 = d23$DRINJ3, idshare = d23$IDSHARE, idrel1 = d23$IDREL1,
                  stringsAsFactors = FALSE)
p23$SEX_PARTNER <- NA_character_
p23$SEX_PARTNER[p23$ident42 == 0 | p23$limh == "Nao" | p23$limm == "Nao" |
                  p23$ident42 == 1 | p23$ident43 == 1] <- "0-1"
p23$SEX_PARTNER[p23$ident42 == 2 | p23$ident43 == 2] <- "2"
p23$SEX_PARTNER[p23$ident42 > 2 | p23$ident43 > 2] <- ">=3"
p23$SHARE_SYRINGE <- NA_character_
p23$SHARE_SYRINGE[p23$drinj3 == "Sim" | p23$idshare == "Sim" |
                    p23$idrel1 == "Sim"] <- "Yes"
p23$SHARE_SYRINGE[is.na(p23$SHARE_SYRINGE) &
                    (p23$drinj3 == "Nao" | p23$idshare == "Nao")] <- "No"
p <- rbind(p14[, c("YEAR", "coupon_id", "SEX_PARTNER", "SHARE_SYRINGE")],
           p23[, c("YEAR", "coupon_id", "SEX_PARTNER", "SHARE_SYRINGE")])
p <- p[!duplicated(p[, c("YEAR", "coupon_id")]), ]
dat$SEX_PARTNER <- dat$SHARE_SYRINGE <- NULL
dat <- merge(dat, p, by = c("YEAR", "coupon_id"), all.x = TRUE, sort = FALSE)

dat$STI_SELF[dat$STI_SELF %in% c("1_SIM", "Yes")] <- "Yes"
dat$STI_SELF[dat$STI_SELF %in% c("2_NAO", "2_SIM", "No")] <- "No"
dat$FREQ_INJEC[dat$FREQ_INJEC == "Less than daily"] <- "Weekly/monthly"
dat$EDUC_CAT[dat$EDUC_CAT == "Primary or none"] <- "No formal/Primary"
dat$EDUC_CAT[dat$EDUC_CAT == "Secondary or higher"] <- "Secondary/Higher"
dat$MARITALC[dat$MARITALC == "Married or union"] <- "Married/Living in union"
dat$AGE_CAT[dat$AGE_CAT == "25+"] <- ">=25"
dat$AGE_FIRST_DRUG[dat$AGE_FIRST_DRUG == "25+"] <- ">=25"

make_cluster <- function(year, city, recruiter, coupon) {
  rec <- as.character(recruiter)
  coup <- as.character(coupon)
  bad <- is.na(rec) | rec %in% c("", "NA", "seed", "SEED", "0")
  rec[bad] <- coup[bad]
  paste(year, city, rec, sep = ":")
}

## LRT-selected terms (from PID_TABLE4_GLMM_LRT.R), with YEAR x CIDADE forced for HIV
hiv_mains <- c("CIDADE", "AGE_CAT", "SEX_CAT", "EDUC_CAT", "MARITALC",
               "FREQ_INJEC", "SHARE_SYRINGE", "UNPROTECTED_SEX", "SEX_PARTNER",
               "HIV_TESTED_12M", "SOUGHT_CARE")
hiv_ints  <- c("CIDADE", "AGE_CAT", "MARITALC", "FREQ_INJEC",
               "SEX_PARTNER", "HIV_TESTED_12M")  ## CIDADE forced

hcv_mains <- c("CIDADE", "AGE_CAT", "SEX_CAT", "AGE_FIRST_DRUG", "FREQ_INJEC",
               "NEW_SYRINGE", "SHARE_SYRINGE", "UNPROTECTED_SEX", "SEX_PARTNER",
               "STI_SELF")
hcv_ints  <- c("CIDADE", "AGE_CAT", "SEX_CAT", "NEW_SYRINGE", "SEX_PARTNER")

all_covars <- unique(c(hiv_mains, hiv_ints, hcv_mains, hcv_ints))

prep_data <- function(dat, outcome, covars) {
  md <- dat
  md$YEAR <- factor(md$YEAR, levels = c("2014", "2023"))
  md$CIDADE <- factor(md$CIDADE, levels = c("MAPUTO", "NAMPULA"))
  md$AGE_CAT <- factor(md$AGE_CAT, levels = c("18-24", ">=25"))
  md$SEX_CAT <- factor(md$SEX_CAT, levels = c("Male", "Female"))
  md$EDUC_CAT <- factor(md$EDUC_CAT, levels = c("No formal/Primary", "Secondary/Higher"))
  md$MARITALC <- factor(md$MARITALC, levels = c("Never married",
                                               "Married/Living in union", "Other"))
  md$AGE_FIRST_DRUG <- factor(md$AGE_FIRST_DRUG, levels = c("<18", "18-24", ">=25"))
  md$SEX_PARTNER <- factor(md$SEX_PARTNER, levels = c("0-1", "2", ">=3"))
  md$FREQ_INJEC <- factor(md$FREQ_INJEC, levels = c("Daily", "Weekly/monthly"))
  md$NEW_SYRINGE <- factor(md$NEW_SYRINGE, levels = c("Yes", "No"))
  md$SHARE_SYRINGE <- factor(md$SHARE_SYRINGE, levels = c("Yes", "No"))
  md$UNPROTECTED_SEX <- factor(md$UNPROTECTED_SEX, levels = c("No", "Yes"))
  md$SOUGHT_CARE <- factor(md$SOUGHT_CARE, levels = c("Yes", "No"))
  md$TX_SEX <- factor(md$TX_SEX, levels = c("Yes", "No"))
  md$HIV_TESTED_12M <- factor(md$HIV_TESTED_12M, levels = c("Yes", "No"))
  md$STI_SELF <- factor(md$STI_SELF, levels = c("No", "Yes"))
  md$y <- as.numeric(md[[outcome]])
  md$cluster_id <- factor(make_cluster(md$YEAR, md$CIDADE, md$recruiter_id, md$coupon_id))
  keep <- complete.cases(md[, c("y", "YEAR", covars, "cluster_id")])
  md[keep, ]
}

ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_glmm <- function(md, mains, ints) {
  env <- new.env(parent = globalenv())
  env$md <- md
  rhs <- unique(c("YEAR", mains))
  if (length(ints)) rhs <- c(rhs, paste0("YEAR:", ints))
  rhs <- c(rhs, "(1 | cluster_id)")
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = env)
  lme4::glmer(fml, data = md, family = binomial(link = "logit"),
              nAGQ = 0, control = ctrl)
}

fit_glm <- function(md, mains, ints) {
  env <- new.env(parent = globalenv())
  env$md <- md
  rhs <- unique(c("YEAR", mains))
  if (length(ints)) rhs <- c(rhs, paste0("YEAR:", ints))
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = env)
  stats::glm(fml, data = md, family = binomial(link = "logit"))
}

## Boundary-corrected LRT for H0: sigma_u = 0 (mixture 0.5*chi^2_0 + 0.5*chi^2_1)
lrt_random <- function(fit_re, fit_fe) {
  ll_re <- as.numeric(logLik(fit_re))
  ll_fe <- as.numeric(logLik(fit_fe))
  lrt <- 2 * (ll_re - ll_fe)
  if (lrt < 0) lrt <- 0
  p_raw <- stats::pchisq(lrt, df = 1, lower.tail = FALSE)
  p_bound <- 0.5 * p_raw
  list(LRT = lrt, df = 1, p_raw = p_raw, p_boundary = p_bound,
       ll_re = ll_re, ll_fe = ll_fe,
       aic_re = AIC(fit_re), aic_fe = AIC(fit_fe))
}

or_table <- function(fit) {
  cf <- lme4::fixef(fit)
  se <- sqrt(diag(as.matrix(vcov(fit))))[names(cf)]
  data.frame(
    term = names(cf),
    or = unname(exp(cf)),
    ci_lo = unname(exp(cf - 1.96 * se)),
    ci_hi = unname(exp(cf + 1.96 * se)),
    p_value = unname(2 * stats::pnorm(-abs(cf / se))),
    stringsAsFactors = FALSE
  )
}

run_one <- function(dat, outcome, mains, ints) {
  md <- prep_data(dat, outcome, unique(c(mains, ints)))
  message(outcome, ": n = ", nrow(md), "; fitting GLMM ...")
  fit_re <- fit_glmm(md, mains, ints)
  fit_fe <- fit_glm(md, mains, ints)
  re_test <- lrt_random(fit_re, fit_fe)
  vc <- as.data.frame(lme4::VarCorr(fit_re))
  message(outcome, ": RE LRT = ", sprintf("%.2f", re_test$LRT),
          "; p (boundary) = ", fmt_p(re_test$p_boundary),
          "; sd(u) = ", sprintf("%.3f", vc$sdcor[1]))
  list(
    outcome = outcome, md = md, fit = fit_re, fit_fe = fit_fe,
    ot = or_table(fit_re), re_test = re_test,
    sigma_u = vc$sdcor[1],
    n = nrow(md),
    n_2014 = sum(md$YEAR == "2014"),
    n_2023 = sum(md$YEAR == "2023"),
    n_pos = sum(md$y == 1),
    n_cl = length(unique(md$cluster_id)),
    mains = mains, ints = ints,
    aic = AIC(fit_re)
  )
}

hiv <- run_one(dat, "HIV", hiv_mains, hiv_ints)
hcv <- run_one(dat, "HCV", hcv_mains, hcv_ints)

## ---- Publication table layout -----------------------------------------------

spec <- rbind(
  data.frame(kind = "section", header = "Sociodemographic characteristics",
             category = "", term = NA_character_, var = NA_character_,
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Age group (ref. 18\u201324)",
             category = "", term = NA_character_, var = "AGE_CAT",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Age group (ref. 18\u201324)",
             category = "\u226525", term = "AGE_CAT>=25", var = "AGE_CAT",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Sex (ref. Male)",
             category = "", term = NA_character_, var = "SEX_CAT",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Sex (ref. Male)",
             category = "Female", term = "SEX_CATFemale", var = "SEX_CAT",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Education (ref. No formal/Primary)",
             category = "", term = NA_character_, var = "EDUC_CAT",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Education (ref. No formal/Primary)",
             category = "Secondary/Higher", term = "EDUC_CATSecondary/Higher",
             var = "EDUC_CAT", interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "City (ref. Maputo)",
             category = "", term = NA_character_, var = "CIDADE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "City (ref. Maputo)",
             category = "Nampula/Nacala", term = "CIDADENAMPULA", var = "CIDADE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Marital status (ref. Never married)",
             category = "", term = NA_character_, var = "MARITALC",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Marital status (ref. Never married)",
             category = "Married/Living in union",
             term = "MARITALCMarried/Living in union", var = "MARITALC",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Marital status (ref. Never married)",
             category = "Other", term = "MARITALCOther", var = "MARITALC",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Age at first drug use (ref. <18)",
             category = "", term = NA_character_, var = "AGE_FIRST_DRUG",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Age at first drug use (ref. <18)",
             category = "18\u201324", term = "AGE_FIRST_DRUG18-24",
             var = "AGE_FIRST_DRUG", interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Age at first drug use (ref. <18)",
             category = "\u226525", term = "AGE_FIRST_DRUG>=25",
             var = "AGE_FIRST_DRUG", interact = FALSE, stringsAsFactors = FALSE),

  data.frame(kind = "section", header = "Injection and sexual behaviours",
             category = "", term = NA_character_, var = NA_character_,
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Injection frequency, past 12 months (ref. Daily)",
             category = "", term = NA_character_, var = "FREQ_INJEC",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Injection frequency, past 12 months (ref. Daily)",
             category = "Weekly/monthly", term = "FREQ_INJECWeekly/monthly",
             var = "FREQ_INJEC", interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Access to new syringes (ref. Yes)",
             category = "", term = NA_character_, var = "NEW_SYRINGE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Access to new syringes (ref. Yes)",
             category = "No", term = "NEW_SYRINGENo", var = "NEW_SYRINGE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Shared syringes, past 12 months (ref. Yes)",
             category = "", term = NA_character_, var = "SHARE_SYRINGE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Shared syringes, past 12 months (ref. Yes)",
             category = "No", term = "SHARE_SYRINGENo", var = "SHARE_SYRINGE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Unprotected sex (ref. No)",
             category = "", term = NA_character_, var = "UNPROTECTED_SEX",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Unprotected sex (ref. No)",
             category = "Yes", term = "UNPROTECTED_SEXYes", var = "UNPROTECTED_SEX",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Number of sexual partners, past 12 months (ref. 0\u20131)",
             category = "", term = NA_character_, var = "SEX_PARTNER",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Number of sexual partners, past 12 months (ref. 0\u20131)",
             category = "2", term = "SEX_PARTNER2", var = "SEX_PARTNER",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Number of sexual partners, past 12 months (ref. 0\u20131)",
             category = "\u22653", term = "SEX_PARTNER>=3", var = "SEX_PARTNER",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Self-reported STI, past 12 months (ref. No)",
             category = "", term = NA_character_, var = "STI_SELF",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Self-reported STI, past 12 months (ref. No)",
             category = "Yes", term = "STI_SELFYes", var = "STI_SELF",
             interact = FALSE, stringsAsFactors = FALSE),

  data.frame(kind = "section", header = "Service use",
             category = "", term = NA_character_, var = NA_character_,
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "HIV tested, past 12 months (ref. Yes)",
             category = "", term = NA_character_, var = "HIV_TESTED_12M",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "HIV tested, past 12 months (ref. Yes)",
             category = "No", term = "HIV_TESTED_12MNo", var = "HIV_TESTED_12M",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Sought healthcare, past 12 months (ref. Yes)",
             category = "", term = NA_character_, var = "SOUGHT_CARE",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Sought healthcare, past 12 months (ref. Yes)",
             category = "No", term = "SOUGHT_CARENo", var = "SOUGHT_CARE",
             interact = FALSE, stringsAsFactors = FALSE),

  data.frame(kind = "section", header = "Survey year",
             category = "", term = NA_character_, var = NA_character_,
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Year (ref. 2014)",
             category = "", term = NA_character_, var = "YEAR",
             interact = FALSE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Year (ref. 2014)",
             category = "2023", term = "YEAR2023", var = "YEAR",
             interact = FALSE, stringsAsFactors = FALSE),

  data.frame(kind = "section", header = "Year interactions",
             category = "", term = NA_character_, var = NA_character_,
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "City \u00d7 year",
             category = "", term = NA_character_, var = "CIDADE",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "City \u00d7 year",
             category = "Nampula/Nacala \u00d7 2023", term = "YEAR2023:CIDADENAMPULA",
             var = "CIDADE", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Age group \u00d7 year",
             category = "", term = NA_character_, var = "AGE_CAT",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Age group \u00d7 year",
             category = "\u226525 \u00d7 2023", term = "YEAR2023:AGE_CAT>=25",
             var = "AGE_CAT", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Sex \u00d7 year",
             category = "", term = NA_character_, var = "SEX_CAT",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Sex \u00d7 year",
             category = "Female \u00d7 2023", term = "YEAR2023:SEX_CATFemale",
             var = "SEX_CAT", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Marital status \u00d7 year",
             category = "", term = NA_character_, var = "MARITALC",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Marital status \u00d7 year",
             category = "Married/Living in union \u00d7 2023",
             term = "YEAR2023:MARITALCMarried/Living in union",
             var = "MARITALC", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Marital status \u00d7 year",
             category = "Other \u00d7 2023", term = "YEAR2023:MARITALCOther",
             var = "MARITALC", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Injection frequency \u00d7 year",
             category = "", term = NA_character_, var = "FREQ_INJEC",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Injection frequency \u00d7 year",
             category = "Weekly/monthly \u00d7 2023",
             term = "YEAR2023:FREQ_INJECWeekly/monthly",
             var = "FREQ_INJEC", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Sexual partners \u00d7 year",
             category = "", term = NA_character_, var = "SEX_PARTNER",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Sexual partners \u00d7 year",
             category = "2 \u00d7 2023", term = "YEAR2023:SEX_PARTNER2",
             var = "SEX_PARTNER", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Sexual partners \u00d7 year",
             category = "\u22653 \u00d7 2023", term = "YEAR2023:SEX_PARTNER>=3",
             var = "SEX_PARTNER", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "HIV tested \u00d7 year",
             category = "", term = NA_character_, var = "HIV_TESTED_12M",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "HIV tested \u00d7 year",
             category = "No \u00d7 2023", term = "YEAR2023:HIV_TESTED_12MNo",
             var = "HIV_TESTED_12M", interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "header", header = "Access to new syringes \u00d7 year",
             category = "", term = NA_character_, var = "NEW_SYRINGE",
             interact = TRUE, stringsAsFactors = FALSE),
  data.frame(kind = "data", header = "Access to new syringes \u00d7 year",
             category = "No \u00d7 2023", term = "YEAR2023:NEW_SYRINGENo",
             var = "NEW_SYRINGE", interact = TRUE, stringsAsFactors = FALSE)
)

fill_side <- function(res) {
  ot <- res$ot
  aor <- p <- character(nrow(spec))
  for (i in seq_len(nrow(spec))) {
    if (spec$kind[i] != "data") { aor[i] <- ""; p[i] <- ""; next }
    if (spec$interact[i] && !(spec$var[i] %in% res$ints)) {
      aor[i] <- "\u2014"; p[i] <- "\u2014"; next
    }
    if (!spec$interact[i] && spec$var[i] != "YEAR" && !(spec$var[i] %in% res$mains)) {
      aor[i] <- "\u2014"; p[i] <- "\u2014"; next
    }
    j <- match(spec$term[i], ot$term)
    if (is.na(j)) { aor[i] <- "\u2014"; p[i] <- "\u2014" } else {
      aor[i] <- fmt_or(ot$or[j], ot$ci_lo[j], ot$ci_hi[j])
      p[i] <- fmt_p(ot$p_value[j])
    }
  }
  list(aor = aor, p = p)
}

hiv_side <- fill_side(hiv)
hcv_side <- fill_side(hcv)

tab <- data.frame(
  kind = spec$kind,
  Characteristic = spec$header,
  Category = spec$category,
  `HIV AOR (95% CI)` = hiv_side$aor,
  `HIV p` = hiv_side$p,
  `HCV AOR (95% CI)` = hcv_side$aor,
  `HCV p` = hcv_side$p,
  check.names = FALSE, stringsAsFactors = FALSE
)

footnote <- paste0(
  "Comparable sample: Maputo and Nampula/Nacala; age \u2265 18 years; injection in the previous 12 months. ",
  "Logistic mixed model (binomial logit) with random recruiter intercept nested in year \u00d7 city. ",
  "Fixed effects selected by backward LRT (\u03b1 = 0.05) from all mains and pairwise year interactions; ",
  "YEAR \u00d7 city was forced into the HIV model. ",
  "Main-effect AORs are associations in 2014; year-interaction AORs describe change in 2023. ",
  "Wald 95% CI for fixed effects. Em dash (\u2014) = term not in the selected model. ",
  "HIV: n = ", hiv$n, " (", hiv$n_2014, " in 2014, ", hiv$n_2023, " in 2023; ",
  hiv$n_pos, " positive; ", hiv$n_cl, " clusters). ",
  "HCV: n = ", hcv$n, " (", hcv$n_2014, " in 2014, ", hcv$n_2023, " in 2023; ",
  hcv$n_pos, " positive; ", hcv$n_cl, " clusters). ",
  "Random-effect LRT (H0: \u03c3u = 0; boundary-corrected p = \u00bd\u00d7\u03c7\u00b2\u2081): ",
  "HIV LRT = ", sprintf("%.2f", hiv$re_test$LRT), ", p = ", fmt_p(hiv$re_test$p_boundary),
  " (\u03c3u = ", sprintf("%.3f", hiv$sigma_u), "); ",
  "HCV LRT = ", sprintf("%.2f", hcv$re_test$LRT), ", p = ", fmt_p(hcv$re_test$p_boundary),
  " (\u03c3u = ", sprintf("%.3f", hcv$sigma_u), ")."
)

title <- paste0(
  "Table. Factors associated with HIV and HCV infection among people who inject drugs: ",
  "logistic mixed models with year interactions (Maputo and Nampula, 2014 and 2023)"
)

## ---- Excel ------------------------------------------------------------------

wb <- createWorkbook()
font_name <- "Times New Roman"

title_style <- createStyle(fontName = font_name, fontSize = 12, textDecoration = "bold",
                           wrapText = TRUE, valign = "center")
note_style <- createStyle(fontName = font_name, fontSize = 8, wrapText = TRUE,
                          valign = "top")
header_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            wrapText = TRUE, valign = "center", halign = "center",
                            border = "TopBottom", borderColour = "black")
section_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                             wrapText = TRUE, valign = "center")
group_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                           wrapText = TRUE, valign = "center")
var_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         wrapText = TRUE)
cat_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         wrapText = TRUE, indent = 1)
num_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         halign = "center")
bottom_style <- createStyle(border = "bottom", borderStyle = "medium",
                            borderColour = "black")

sh <- "Table"
addWorksheet(wb, sh, gridLines = FALSE)
pageSetup(wb, sh, orientation = "landscape", fitToWidth = TRUE, paperSize = 9,
          left = 0.4, right = 0.4, top = 0.5, bottom = 0.5)
setColWidths(wb, sh, cols = 1:6, widths = c(44, 28, 22, 10, 22, 10))

writeData(wb, sh, title, startRow = 1, colNames = FALSE)
mergeCells(wb, sh, cols = 1:6, rows = 1)
addStyle(wb, sh, title_style, rows = 1, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, sh, rows = 1, heights = 32)

hdr <- 3
writeData(wb, sh, "Characteristic", startRow = hdr, startCol = 1, colNames = FALSE)
writeData(wb, sh, "Category", startRow = hdr, startCol = 2, colNames = FALSE)
writeData(wb, sh, "HIV", startRow = hdr, startCol = 3, colNames = FALSE)
writeData(wb, sh, "HCV", startRow = hdr, startCol = 5, colNames = FALSE)
writeData(wb, sh, "AOR (95% CI)", startRow = hdr + 1, startCol = 3, colNames = FALSE)
writeData(wb, sh, "p", startRow = hdr + 1, startCol = 4, colNames = FALSE)
writeData(wb, sh, "AOR (95% CI)", startRow = hdr + 1, startCol = 5, colNames = FALSE)
writeData(wb, sh, "p", startRow = hdr + 1, startCol = 6, colNames = FALSE)
mergeCells(wb, sh, cols = 1, rows = hdr:(hdr + 1))
mergeCells(wb, sh, cols = 2, rows = hdr:(hdr + 1))
mergeCells(wb, sh, cols = 3:4, rows = hdr)
mergeCells(wb, sh, cols = 5:6, rows = hdr)
addStyle(wb, sh, header_style, rows = hdr:(hdr + 1), cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, sh, rows = hdr:(hdr + 1), heights = 20)

r <- hdr + 1
for (i in seq_len(nrow(tab))) {
  r <- r + 1
  if (tab$kind[i] == "section") {
    writeData(wb, sh, tab$Characteristic[i], startRow = r, startCol = 1, colNames = FALSE)
    mergeCells(wb, sh, cols = 1:6, rows = r)
    addStyle(wb, sh, section_style, rows = r, cols = 1:6, gridExpand = TRUE)
  } else if (tab$kind[i] == "header") {
    writeData(wb, sh, tab$Characteristic[i], startRow = r, startCol = 1, colNames = FALSE)
    mergeCells(wb, sh, cols = 1:6, rows = r)
    addStyle(wb, sh, group_style, rows = r, cols = 1:6, gridExpand = TRUE)
  } else {
    vals <- c("", tab$Category[i],
              tab$`HIV AOR (95% CI)`[i], tab$`HIV p`[i],
              tab$`HCV AOR (95% CI)`[i], tab$`HCV p`[i])
    for (j in seq_along(vals))
      writeData(wb, sh, vals[j], startRow = r, startCol = j, colNames = FALSE)
    addStyle(wb, sh, var_style, rows = r, cols = 1)
    addStyle(wb, sh, cat_style, rows = r, cols = 2)
    addStyle(wb, sh, num_style, rows = r, cols = 3:6, gridExpand = TRUE)
  }
}
addStyle(wb, sh, bottom_style, rows = r, cols = 1:6, gridExpand = TRUE, stack = TRUE)

note_row <- r + 2
writeData(wb, sh, paste0("Footnote. ", footnote), startRow = note_row, colNames = FALSE)
mergeCells(wb, sh, cols = 1:6, rows = note_row)
addStyle(wb, sh, note_style, rows = note_row, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, sh, rows = note_row, heights = 72)
freezePane(wb, sh, firstActiveRow = hdr + 2)

## Random-effect LRT sheet
addWorksheet(wb, "Random-effect LRT", gridLines = FALSE)
re_tab <- data.frame(
  Outcome = c("HIV", "HCV"),
  `n` = c(hiv$n, hcv$n),
  Clusters = c(hiv$n_cl, hcv$n_cl),
  `sd(u)` = c(sprintf("%.3f", hiv$sigma_u), sprintf("%.3f", hcv$sigma_u)),
  `logLik GLMM` = c(sprintf("%.2f", hiv$re_test$ll_re), sprintf("%.2f", hcv$re_test$ll_re)),
  `logLik GLM` = c(sprintf("%.2f", hiv$re_test$ll_fe), sprintf("%.2f", hcv$re_test$ll_fe)),
  LRT = c(sprintf("%.2f", hiv$re_test$LRT), sprintf("%.2f", hcv$re_test$LRT)),
  `p (chi2_1)` = c(fmt_p(hiv$re_test$p_raw), fmt_p(hcv$re_test$p_raw)),
  `p (boundary)` = c(fmt_p(hiv$re_test$p_boundary), fmt_p(hcv$re_test$p_boundary)),
  `AIC GLMM` = c(sprintf("%.1f", hiv$re_test$aic_re), sprintf("%.1f", hcv$re_test$aic_re)),
  `AIC GLM` = c(sprintf("%.1f", hiv$re_test$aic_fe), sprintf("%.1f", hcv$re_test$aic_fe)),
  check.names = FALSE, stringsAsFactors = FALSE
)
writeData(wb, "Random-effect LRT",
          "Likelihood-ratio test of recruiter random intercept (H0: sigma_u = 0)",
          startRow = 1, colNames = FALSE)
mergeCells(wb, "Random-effect LRT", cols = 1:11, rows = 1)
addStyle(wb, "Random-effect LRT", title_style, rows = 1, cols = 1:11, gridExpand = TRUE)
writeData(wb, "Random-effect LRT", re_tab, startRow = 3)
addStyle(wb, "Random-effect LRT", header_style, rows = 3, cols = 1:11, gridExpand = TRUE)
addStyle(wb, "Random-effect LRT", num_style, rows = 4:5, cols = 1:11, gridExpand = TRUE)
setColWidths(wb, "Random-effect LRT", cols = 1:11, widths = c(10, 8, 10, 10, 12, 12, 8, 12, 14, 10, 10))
writeData(wb, "Random-effect LRT",
          paste0("Boundary p-value uses the 50:50 mixture of chi-square(0) and chi-square(1) ",
                 "appropriate for testing a variance component on the boundary."),
          startRow = 7, colNames = FALSE)
mergeCells(wb, "Random-effect LRT", cols = 1:11, rows = 7)
addStyle(wb, "Random-effect LRT", note_style, rows = 7, cols = 1:11, gridExpand = TRUE)

## Notes
addWorksheet(wb, "Notes")
notes <- data.frame(
  Item = c("HIV fixed effects", "HIV year interactions",
           "HCV fixed effects", "HCV year interactions",
           "HIV YEAR x city", "Random effect", "RE LRT HIV", "RE LRT HCV", "Script"),
  Detail = c(
    paste(c("YEAR", hiv_mains), collapse = "; "),
    paste(paste0("YEAR x ", hiv_ints), collapse = "; "),
    paste(c("YEAR", hcv_mains), collapse = "; "),
    paste(paste0("YEAR x ", hcv_ints), collapse = "; "),
    "Forced into HIV model (not retained by LRT alone; LRT p was 0.112 in selection).",
    "Recruiter nested in year x city.",
    sprintf("LRT=%.2f; boundary p=%s; sd(u)=%.3f",
            hiv$re_test$LRT, fmt_p(hiv$re_test$p_boundary), hiv$sigma_u),
    sprintf("LRT=%.2f; boundary p=%s; sd(u)=%.3f",
            hcv$re_test$LRT, fmt_p(hcv$re_test$p_boundary), hcv$sigma_u),
    "PID_TABLE4_GLMM_FINAL.R"
  ), stringsAsFactors = FALSE
)
writeData(wb, "Notes", notes)
setColWidths(wb, "Notes", cols = 1:2, widths = c(24, 110))
addStyle(wb, "Notes", note_style, rows = 1:(nrow(notes) + 1), cols = 1:2, gridExpand = TRUE)

xlsx_file <- file.path(path_out, "Table4_GLMM_final_publish.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

write.csv(tab, file.path(path_tab, "table4_glmm_final_publish.csv"), row.names = FALSE)
write.csv(re_tab, file.path(path_tab, "table4_glmm_re_lrt.csv"), row.names = FALSE)
write.csv(hiv$ot, file.path(path_tab, "table4_glmm_final_hiv_coef.csv"), row.names = FALSE)
write.csv(hcv$ot, file.path(path_tab, "table4_glmm_final_hcv_coef.csv"), row.names = FALSE)

message("Wrote ", xlsx_file)
print(re_tab)
print(tab[tab$kind == "data", c("Characteristic", "Category",
                                 "HIV AOR (95% CI)", "HIV p",
                                 "HCV AOR (95% CI)", "HCV p")])
