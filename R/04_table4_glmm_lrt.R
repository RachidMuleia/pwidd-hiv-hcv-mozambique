################################################################################
## Logistic mixed model (random recruiter intercept) + backward LRT selection
##
## 1. Drop USED_NEEDLE, SHARED_EQUIP, RISK_PERCEPTION.
## 2. Full GLMM: y ~ YEAR + all X + all YEAR:X + (1 | cluster)
##    cluster = recruiter nested in year x city.
## 3. Backward LRT (alpha = 0.05): drop YEAR interactions first (joint per
##    factor), then mains not protected by a retained interaction. YEAR fixed.
## 4. Report fixed-effect AORs from the selected GLMM (Wald 95% CI).
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
library(officer)
library(flextable)

alpha_lrt <- 0.05

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_or <- function(or, lo, hi) {
  if (is.na(or)) return(NA_character_)
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
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

covars <- c("CIDADE", "AGE_CAT", "SEX_CAT", "EDUC_CAT", "MARITALC",
            "AGE_FIRST_DRUG", "FREQ_INJEC", "NEW_SYRINGE", "SHARE_SYRINGE",
            "UNPROTECTED_SEX", "SEX_PARTNER", "TX_SEX",
            "STI_SELF", "HIV_TESTED_12M", "SOUGHT_CARE")

int_lab <- c(
  CIDADE = "City", AGE_CAT = "Age group", SEX_CAT = "Sex",
  EDUC_CAT = "Education", MARITALC = "Marital status",
  AGE_FIRST_DRUG = "Age at first drug use", FREQ_INJEC = "Injection frequency",
  NEW_SYRINGE = "Access to new syringes",
  SHARE_SYRINGE = "Shared syringes, past 12 months",
  UNPROTECTED_SEX = "Unprotected sex", SEX_PARTNER = "Number of sexual partners",
  TX_SEX = "Transactional sex", STI_SELF = "Self-reported STI",
  HIV_TESTED_12M = "HIV tested, past 12 months", SOUGHT_CARE = "Sought healthcare"
)

prep_data <- function(dat, outcome) {
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
  md <- md[keep, ]
  md
}

make_fml <- function(mains, ints) {
  rhs <- unique(c("YEAR", mains))
  if (length(ints)) rhs <- c(rhs, paste0("YEAR:", ints))
  rhs <- c(rhs, "(1 | cluster_id)")
  stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_glmm <- function(md, mains, ints) {
  ## Keep data in the formula environment so nested refits / anova work
  env <- new.env(parent = globalenv())
  env$md <- md
  rhs <- unique(c("YEAR", mains))
  if (length(ints)) rhs <- c(rhs, paste0("YEAR:", ints))
  rhs <- c(rhs, "(1 | cluster_id)")
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = env)
  lme4::glmer(fml, data = md, family = binomial(link = "logit"),
              nAGQ = 0, control = ctrl)
}

or_table <- function(fit) {
  cf <- lme4::fixef(fit)
  V <- as.matrix(vcov(fit))
  se <- sqrt(diag(V))[names(cf)]
  data.frame(
    term = names(cf),
    or = unname(exp(cf)),
    ci_lo = unname(exp(cf - 1.96 * se)),
    ci_hi = unname(exp(cf + 1.96 * se)),
    p_value = unname(2 * stats::pnorm(-abs(cf / se))),
    stringsAsFactors = FALSE
  )
}

## Nested LRT: current fit vs fit without one term
lrt_p <- function(fit, fit_red) {
  an <- anova(fit, fit_red)
  list(
    Chisq = as.numeric(an$Chisq[2]),
    df = as.integer(an$Df[2]),
    p = as.numeric(an$`Pr(>Chisq)`[2])
  )
}

## Backward LRT (joint tests for multi-df factors via whole YEAR:X / X terms)
backward_lrt <- function(md, outcome) {
  mains <- covars
  ints <- covars
  path <- list()
  step <- 0L

  message(outcome, ": fitting full GLMM ...")
  fit <- fit_glmm(md, mains, ints)
  message(outcome, ": full AIC = ", sprintf("%.1f", AIC(fit)))

  ## Phase 1: YEAR interactions
  repeat {
    if (!length(ints)) break
    message(outcome, ": LRT interactions (", length(ints), " candidates) ...")
    tests <- lapply(ints, function(v) {
      fit_red <- tryCatch(fit_glmm(md, mains, setdiff(ints, v)), error = function(e) NULL)
      if (is.null(fit_red)) {
        return(data.frame(term = paste0("YEAR:", v), Chisq = NA_real_,
                          df = NA_integer_, p = NA_real_, stringsAsFactors = FALSE))
      }
      r <- tryCatch(lrt_p(fit, fit_red), error = function(e) list(Chisq = NA, df = NA, p = NA))
      data.frame(term = paste0("YEAR:", v), Chisq = r$Chisq, df = r$df, p = r$p,
                 stringsAsFactors = FALSE)
    })
    tt <- do.call(rbind, tests)
    tt <- tt[order(-tt$p, na.last = TRUE), ]
    worst <- tt[1, ]
    message(outcome, ": worst interaction = ", worst$term, " p = ", fmt_p(worst$p))
    if (is.na(worst$p) || worst$p > alpha_lrt) {
      drop_v <- sub("^YEAR:", "", worst$term)
      step <- step + 1L
      path[[length(path) + 1L]] <- data.frame(
        step = step, action = paste0("- ", worst$term),
        Chisq = worst$Chisq, df = worst$df, p = worst$p,
        stringsAsFactors = FALSE
      )
      message(outcome, ": drop ", worst$term, " (LRT p = ", fmt_p(worst$p), ")")
      ints <- setdiff(ints, drop_v)
      fit <- fit_glmm(md, mains, ints)
    } else break
  }

  ## Phase 2: mains not protected by a retained interaction
  repeat {
    cands <- setdiff(mains, ints)
    if (!length(cands)) break
    message(outcome, ": LRT mains (", length(cands), " candidates) ...")
    tests <- lapply(cands, function(v) {
      fit_red <- tryCatch(fit_glmm(md, setdiff(mains, v), ints), error = function(e) NULL)
      if (is.null(fit_red)) {
        return(data.frame(term = v, Chisq = NA_real_, df = NA_integer_,
                          p = NA_real_, stringsAsFactors = FALSE))
      }
      r <- tryCatch(lrt_p(fit, fit_red), error = function(e) list(Chisq = NA, df = NA, p = NA))
      data.frame(term = v, Chisq = r$Chisq, df = r$df, p = r$p, stringsAsFactors = FALSE)
    })
    tt <- do.call(rbind, tests)
    tt <- tt[order(-tt$p, na.last = TRUE), ]
    worst <- tt[1, ]
    message(outcome, ": worst main = ", worst$term, " p = ", fmt_p(worst$p))
    if (is.na(worst$p) || worst$p > alpha_lrt) {
      step <- step + 1L
      path[[length(path) + 1L]] <- data.frame(
        step = step, action = paste0("- ", worst$term),
        Chisq = worst$Chisq, df = worst$df, p = worst$p,
        stringsAsFactors = FALSE
      )
      message(outcome, ": drop ", worst$term, " (LRT p = ", fmt_p(worst$p), ")")
      mains <- setdiff(mains, worst$term)
      fit <- fit_glmm(md, mains, ints)
    } else break
  }

  path_df <- if (length(path)) do.call(rbind, path) else
    data.frame(step = integer(), action = character(), Chisq = numeric(),
               df = integer(), p = numeric(), stringsAsFactors = FALSE)

  vc <- as.data.frame(lme4::VarCorr(fit))
  list(
    fit = fit, mains = mains, ints = ints, path = path_df,
    ot = or_table(fit), aic = AIC(fit),
    fml = paste(deparse(formula(fit), width.cutoff = 500), collapse = " "),
    sigma_u = vc$sdcor[1],
    n_cl = length(unique(md$cluster_id))
  )
}

run_outcome <- function(dat, outcome) {
  md <- prep_data(dat, outcome)
  message(outcome, " complete-case n = ", nrow(md),
          "; clusters = ", length(unique(md$cluster_id)))
  sel <- backward_lrt(md, outcome)
  message(outcome, " selected AIC = ", sprintf("%.1f", sel$aic))
  message(outcome, " kept YEAR x: ",
          if (length(sel$ints)) paste(sel$ints, collapse = ", ") else "none")
  message(outcome, " kept mains: ", paste(c("YEAR", sel$mains), collapse = ", "))
  list(
    outcome = outcome, md = md, fit = sel$fit, final = sel$fit,
    aic = sel$aic, path = sel$path, fml = sel$fml,
    ints = sel$ints, mains = sel$mains, ot = sel$ot,
    sigma_u = sel$sigma_u, n_cl = sel$n_cl,
    n = nrow(md),
    n_2014 = sum(md$YEAR == "2014"),
    n_2023 = sum(md$YEAR == "2023"),
    n_pos = sum(md$y == 1)
  )
}

message("HIV: GLMM + LRT ...")
hiv <- run_outcome(dat, "HIV")
message("HCV: GLMM + LRT ...")
hcv <- run_outcome(dat, "HCV")

## ---- Table shell (same layout as before) ------------------------------------

spec <- rbind(
  data.frame(kind = "header", var = "AGE_CAT", interact = FALSE,
             header = "Age group (Ref. = 18-24)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_CAT", interact = FALSE,
             header = "Age group (Ref. = 18-24)", category = ">=25",
             term = "AGE_CAT>=25", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "YEAR", interact = FALSE,
             header = "Year (Ref. = 2014)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "YEAR", interact = FALSE,
             header = "Year (Ref. = 2014)", category = "2023",
             term = "YEAR2023", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SEX_CAT", interact = FALSE,
             header = "Sex (Ref. = Male)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_CAT", interact = FALSE,
             header = "Sex (Ref. = Male)", category = "Female",
             term = "SEX_CATFemale", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "EDUC_CAT", interact = FALSE,
             header = "Education (Ref. = No formal/Primary)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "EDUC_CAT", interact = FALSE,
             header = "Education (Ref. = No formal/Primary)", category = "Secondary/Higher",
             term = "EDUC_CATSecondary/Higher", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "CIDADE", interact = FALSE,
             header = "City (Ref. = Maputo)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "CIDADE", interact = FALSE,
             header = "City (Ref. = Maputo)", category = "Nampula/Nacala",
             term = "CIDADENAMPULA", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "MARITALC", interact = FALSE,
             header = "Marital status (Ref. = Never married)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "MARITALC", interact = FALSE,
             header = "Marital status (Ref. = Never married)",
             category = "Married/Living in union",
             term = "MARITALCMarried/Living in union", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "MARITALC", interact = FALSE,
             header = "Marital status (Ref. = Never married)", category = "Other",
             term = "MARITALCOther", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "AGE_FIRST_DRUG", interact = FALSE,
             header = "Age at first drug use (Ref. = <18)", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_FIRST_DRUG", interact = FALSE,
             header = "Age at first drug use (Ref. = <18)", category = "18-24",
             term = "AGE_FIRST_DRUG18-24", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_FIRST_DRUG", interact = FALSE,
             header = "Age at first drug use (Ref. = <18)", category = ">=25",
             term = "AGE_FIRST_DRUG>=25", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SEX_PARTNER", interact = FALSE,
             header = "Number of sexual partners, past 12 months (Ref. = 0-1)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_PARTNER", interact = FALSE,
             header = "Number of sexual partners, past 12 months (Ref. = 0-1)",
             category = "2", term = "SEX_PARTNER2", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_PARTNER", interact = FALSE,
             header = "Number of sexual partners, past 12 months (Ref. = 0-1)",
             category = ">=3", term = "SEX_PARTNER>=3", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "FREQ_INJEC", interact = FALSE,
             header = "Injection frequency, past 12 months (Ref. = Daily)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "FREQ_INJEC", interact = FALSE,
             header = "Injection frequency, past 12 months (Ref. = Daily)",
             category = "Weekly/monthly", term = "FREQ_INJECWeekly/monthly",
             stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "NEW_SYRINGE", interact = FALSE,
             header = "Access to new syringes (Ref. = Yes)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "NEW_SYRINGE", interact = FALSE,
             header = "Access to new syringes (Ref. = Yes)",
             category = "No", term = "NEW_SYRINGENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SHARE_SYRINGE", interact = FALSE,
             header = "Shared syringes, past 12 months (Ref. = Yes)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SHARE_SYRINGE", interact = FALSE,
             header = "Shared syringes, past 12 months (Ref. = Yes)",
             category = "No", term = "SHARE_SYRINGENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "UNPROTECTED_SEX", interact = FALSE,
             header = "Unprotected sex (Ref. = No)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "UNPROTECTED_SEX", interact = FALSE,
             header = "Unprotected sex (Ref. = No)",
             category = "Yes", term = "UNPROTECTED_SEXYes", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SOUGHT_CARE", interact = FALSE,
             header = "Sought healthcare, past 12 months (Ref. = Yes)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SOUGHT_CARE", interact = FALSE,
             header = "Sought healthcare, past 12 months (Ref. = Yes)",
             category = "No", term = "SOUGHT_CARENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "TX_SEX", interact = FALSE,
             header = "Received money, goods or services for sex, past 12 months (Ref. = Yes)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "TX_SEX", interact = FALSE,
             header = "Received money, goods or services for sex, past 12 months (Ref. = Yes)",
             category = "No", term = "TX_SEXNo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "HIV_TESTED_12M", interact = FALSE,
             header = "HIV tested, past 12 months (Ref. = Yes)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "HIV_TESTED_12M", interact = FALSE,
             header = "HIV tested, past 12 months (Ref. = Yes)",
             category = "No", term = "HIV_TESTED_12MNo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "STI_SELF", interact = FALSE,
             header = "Self-reported STI, past 12 months (Ref. = No)",
             category = "", term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "STI_SELF", interact = FALSE,
             header = "Self-reported STI, past 12 months (Ref. = No)",
             category = "Yes", term = "STI_SELFYes", stringsAsFactors = FALSE),
  data.frame(kind = "section", var = NA_character_, interact = TRUE,
             header = "Year interactions retained by LRT", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SEX_CAT", interact = TRUE,
             header = "Sex x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_CAT", interact = TRUE,
             header = "Sex x Year", category = "Female x 2023",
             term = "YEAR2023:SEX_CATFemale", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "AGE_CAT", interact = TRUE,
             header = "Age group x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_CAT", interact = TRUE,
             header = "Age group x Year", category = ">=25 x 2023",
             term = "YEAR2023:AGE_CAT>=25", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "EDUC_CAT", interact = TRUE,
             header = "Education x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "EDUC_CAT", interact = TRUE,
             header = "Education x Year", category = "Secondary/Higher x 2023",
             term = "YEAR2023:EDUC_CATSecondary/Higher", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "MARITALC", interact = TRUE,
             header = "Marital status x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "MARITALC", interact = TRUE,
             header = "Marital status x Year", category = "Married/Living in union x 2023",
             term = "YEAR2023:MARITALCMarried/Living in union", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "MARITALC", interact = TRUE,
             header = "Marital status x Year", category = "Other x 2023",
             term = "YEAR2023:MARITALCOther", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "CIDADE", interact = TRUE,
             header = "City x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "CIDADE", interact = TRUE,
             header = "City x Year", category = "Nampula/Nacala x 2023",
             term = "YEAR2023:CIDADENAMPULA", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "FREQ_INJEC", interact = TRUE,
             header = "Injection frequency x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "FREQ_INJEC", interact = TRUE,
             header = "Injection frequency x Year", category = "Weekly/monthly x 2023",
             term = "YEAR2023:FREQ_INJECWeekly/monthly", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SEX_PARTNER", interact = TRUE,
             header = "Number of sexual partners x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_PARTNER", interact = TRUE,
             header = "Number of sexual partners x Year", category = "2 x 2023",
             term = "YEAR2023:SEX_PARTNER2", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SEX_PARTNER", interact = TRUE,
             header = "Number of sexual partners x Year", category = ">=3 x 2023",
             term = "YEAR2023:SEX_PARTNER>=3", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "TX_SEX", interact = TRUE,
             header = "Transactional sex x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "TX_SEX", interact = TRUE,
             header = "Transactional sex x Year", category = "No x 2023",
             term = "YEAR2023:TX_SEXNo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "HIV_TESTED_12M", interact = TRUE,
             header = "HIV tested x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "HIV_TESTED_12M", interact = TRUE,
             header = "HIV tested x Year", category = "No x 2023",
             term = "YEAR2023:HIV_TESTED_12MNo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "AGE_FIRST_DRUG", interact = TRUE,
             header = "Age at first drug use x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_FIRST_DRUG", interact = TRUE,
             header = "Age at first drug use x Year", category = "18-24 x 2023",
             term = "YEAR2023:AGE_FIRST_DRUG18-24", stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "AGE_FIRST_DRUG", interact = TRUE,
             header = "Age at first drug use x Year", category = ">=25 x 2023",
             term = "YEAR2023:AGE_FIRST_DRUG>=25", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SOUGHT_CARE", interact = TRUE,
             header = "Sought healthcare x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SOUGHT_CARE", interact = TRUE,
             header = "Sought healthcare x Year", category = "No x 2023",
             term = "YEAR2023:SOUGHT_CARENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "SHARE_SYRINGE", interact = TRUE,
             header = "Shared syringes x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "SHARE_SYRINGE", interact = TRUE,
             header = "Shared syringes x Year", category = "No x 2023",
             term = "YEAR2023:SHARE_SYRINGENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "NEW_SYRINGE", interact = TRUE,
             header = "Access to new syringes x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "NEW_SYRINGE", interact = TRUE,
             header = "Access to new syringes x Year", category = "No x 2023",
             term = "YEAR2023:NEW_SYRINGENo", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "UNPROTECTED_SEX", interact = TRUE,
             header = "Unprotected sex x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "UNPROTECTED_SEX", interact = TRUE,
             header = "Unprotected sex x Year", category = "Yes x 2023",
             term = "YEAR2023:UNPROTECTED_SEXYes", stringsAsFactors = FALSE),
  data.frame(kind = "header", var = "STI_SELF", interact = TRUE,
             header = "Self-reported STI x Year", category = "",
             term = NA_character_, stringsAsFactors = FALSE),
  data.frame(kind = "data", var = "STI_SELF", interact = TRUE,
             header = "Self-reported STI x Year", category = "Yes x 2023",
             term = "YEAR2023:STI_SELFYes", stringsAsFactors = FALSE)
)

fill_side <- function(res) {
  ot <- res$ot
  aor <- p <- character(nrow(spec))
  for (i in seq_len(nrow(spec))) {
    if (spec$kind[i] != "data") { aor[i] <- ""; p[i] <- ""; next }
    if (spec$interact[i] && !(spec$var[i] %in% res$ints)) {
      aor[i] <- "*"; p[i] <- "*"; next
    }
    if (!spec$interact[i] && spec$var[i] != "YEAR" && !(spec$var[i] %in% res$mains)) {
      aor[i] <- "*"; p[i] <- "*"; next
    }
    j <- match(spec$term[i], ot$term)
    if (is.na(j)) { aor[i] <- "*"; p[i] <- "*" } else {
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
  Variable = gsub(">=", "\u2265", spec$header, fixed = TRUE),
  Category = gsub(">=", "\u2265", spec$category, fixed = TRUE),
  `HIV AOR (95% CI)` = hiv_side$aor,
  `HIV p` = hiv_side$p,
  `HCV AOR (95% CI)` = hcv_side$aor,
  `HCV p` = hcv_side$p,
  check.names = FALSE, stringsAsFactors = FALSE
)

sel_wide <- data.frame(
  Covariate = unname(int_lab[covars]),
  `HIV main` = ifelse(covars %in% hiv$mains, "Yes", "No"),
  `HIV YEAR x` = ifelse(covars %in% hiv$ints, "Yes", "No"),
  `HCV main` = ifelse(covars %in% hcv$mains, "Yes", "No"),
  `HCV YEAR x` = ifelse(covars %in% hcv$ints, "Yes", "No"),
  check.names = FALSE, stringsAsFactors = FALSE
)

write.csv(tab, file.path(path_tab, "table4_glmm_lrt.csv"), row.names = FALSE)
write.csv(hiv$ot, file.path(path_tab, "table4_glmm_lrt_hiv_coef.csv"), row.names = FALSE)
write.csv(hcv$ot, file.path(path_tab, "table4_glmm_lrt_hcv_coef.csv"), row.names = FALSE)
write.csv(sel_wide, file.path(path_tab, "table4_glmm_lrt_selected.csv"), row.names = FALSE)
write.csv(hiv$path, file.path(path_tab, "table4_glmm_lrt_hiv_path.csv"), row.names = FALSE)
write.csv(hcv$path, file.path(path_tab, "table4_glmm_lrt_hcv_path.csv"), row.names = FALSE)

n_note <- paste0(
  "Comparable sample: age \u2265 18 years, Maputo and Nampula/Nacala, injection in the previous 12 months. ",
  "Logistic mixed model (binomial logit) with random recruiter intercept nested in year \u00d7 city; ",
  "all listed covariates plus all pairwise YEAR \u00d7 interactions entered; ",
  "backward LRT selection (\u03b1 = 0.05) dropped non-significant interactions first, then unprotected mains. ",
  "YEAR always retained; hierarchy preserved. Dropped a priori: used syringe, shared other equipment, risk perception. ",
  "HIV n = ", hiv$n, " (", hiv$n_cl, " clusters); HCV n = ", hcv$n, " (", hcv$n_cl, " clusters). ",
  "* Not retained by LRT. AOR from fixed effects; Wald 95% CI."
)
title <- "Table. HIV and HCV: logistic mixed model after LRT selection of year interactions and mains"

wb <- createWorkbook()
font_name <- "Times New Roman"
header_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            wrapText = TRUE, valign = "center", halign = "center",
                            fgFill = "#1F4E79", fontColour = "white")
group_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                           wrapText = TRUE, fgFill = "#D6E3F0", fontColour = "#1F4E79")
section_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                             fgFill = "#1F4E79", fontColour = "white")
num_style <- createStyle(fontName = font_name, fontSize = 10, halign = "center", valign = "center")
var_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center", wrapText = TRUE)
note_style <- createStyle(fontName = font_name, fontSize = 8, wrapText = TRUE)
title_style <- createStyle(fontName = font_name, fontSize = 12, textDecoration = "bold", wrapText = TRUE)
yes_style <- createStyle(fontName = font_name, fontSize = 10, halign = "center",
                         textDecoration = "bold", fontColour = "#1F4E79")

sh <- "Final model"
addWorksheet(wb, sh, gridLines = FALSE)
writeData(wb, sh, title, startRow = 1, colNames = FALSE)
mergeCells(wb, sh, cols = 1:6, rows = 1)
addStyle(wb, sh, title_style, rows = 1, cols = 1:6, gridExpand = TRUE)
writeData(wb, sh, n_note, startRow = 2, colNames = FALSE)
mergeCells(wb, sh, cols = 1:6, rows = 2)
addStyle(wb, sh, note_style, rows = 2, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, sh, rows = 2, heights = 56)
writeData(wb, sh, data.frame(
  Variable = tab$Variable, Category = tab$Category,
  `HIV AOR (95% CI)` = tab$`HIV AOR (95% CI)`, `HIV p` = tab$`HIV p`,
  `HCV AOR (95% CI)` = tab$`HCV AOR (95% CI)`, `HCV p` = tab$`HCV p`,
  check.names = FALSE
), startRow = 4)
addStyle(wb, sh, header_style, rows = 4, cols = 1:6, gridExpand = TRUE)
for (i in seq_len(nrow(tab))) {
  r <- 4 + i
  if (tab$kind[i] == "section") addStyle(wb, sh, section_style, rows = r, cols = 1:6, gridExpand = TRUE)
  else if (tab$kind[i] == "header") addStyle(wb, sh, group_style, rows = r, cols = 1:6, gridExpand = TRUE)
  else {
    addStyle(wb, sh, var_style, rows = r, cols = 1:2, gridExpand = TRUE)
    addStyle(wb, sh, num_style, rows = r, cols = 3:6, gridExpand = TRUE)
  }
}
setColWidths(wb, sh, cols = 1:6, widths = c(40, 28, 22, 10, 22, 10))

addWorksheet(wb, "LRT selection")
writeData(wb, "LRT selection", sel_wide)
addStyle(wb, "LRT selection", header_style, rows = 1, cols = 1:5, gridExpand = TRUE)
for (i in seq_len(nrow(sel_wide))) {
  for (j in 2:5) if (sel_wide[i, j] == "Yes")
    addStyle(wb, "LRT selection", yes_style, rows = 1 + i, cols = j)
}
setColWidths(wb, "LRT selection", cols = 1:5, widths = c(32, 12, 12, 12, 12))

addWorksheet(wb, "HIV LRT path")
writeData(wb, "HIV LRT path", hiv$path)
addWorksheet(wb, "HCV LRT path")
writeData(wb, "HCV LRT path", hcv$path)

addWorksheet(wb, "Notes")
notes <- data.frame(
  Item = c("Model", "Random effect", "Selection", "Alpha",
           "HIV n", "HIV clusters", "HIV sd(u)", "HIV AIC", "HIV YEAR x",
           "HCV n", "HCV clusters", "HCV sd(u)", "HCV AIC", "HCV YEAR x", "Script"),
  Detail = c(
    "glmer binomial logit, nAGQ = 0, bobyqa.",
    "Recruiter nested in year x city (seeds are singleton clusters).",
    "Backward LRT: interactions first, then mains without retained interaction. YEAR fixed.",
    as.character(alpha_lrt),
    paste0(hiv$n, " (", hiv$n_2014, "/", hiv$n_2023, "); ", hiv$n_pos, " positive."),
    as.character(hiv$n_cl), sprintf("%.3f", hiv$sigma_u), sprintf("%.1f", hiv$aic),
    if (length(hiv$ints)) paste(int_lab[hiv$ints], collapse = "; ") else "none",
    paste0(hcv$n, " (", hcv$n_2014, "/", hcv$n_2023, "); ", hcv$n_pos, " positive."),
    as.character(hcv$n_cl), sprintf("%.3f", hcv$sigma_u), sprintf("%.1f", hcv$aic),
    if (length(hcv$ints)) paste(int_lab[hcv$ints], collapse = "; ") else "none",
    "PID_TABLE4_GLMM_LRT.R"
  ), stringsAsFactors = FALSE
)
writeData(wb, "Notes", notes)
setColWidths(wb, "Notes", cols = 1:2, widths = c(18, 100))

xlsx_file <- file.path(path_out, "Table4_GLMM_LRT.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

docx_file <- file.path(path_out, "Table4_GLMM_LRT.docx")
d <- data.frame(
  Variable = tab$Variable,
  Category = ifelse(tab$kind == "data", tab$Category, ""),
  `HIV AOR (95% CI)` = tab$`HIV AOR (95% CI)`, `HIV p` = tab$`HIV p`,
  `HCV AOR (95% CI)` = tab$`HCV AOR (95% CI)`, `HCV p` = tab$`HCV p`,
  check.names = FALSE, stringsAsFactors = FALSE
)
ft <- flextable(d)
ft <- set_header_labels(ft, Variable = "Variable", Category = "Category",
                        `HIV AOR (95% CI)` = "AOR (95% CI)", `HIV p` = "p",
                        `HCV AOR (95% CI)` = "AOR (95% CI)", `HCV p` = "p")
ft <- add_header_row(ft, values = c("Variable", "Category", "HIV", "HCV"),
                     colwidths = c(1, 1, 2, 2))
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 8, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, j = 3:6, align = "center", part = "all")
ft <- bg(ft, part = "header", bg = "#1F4E79")
ft <- color(ft, part = "header", color = "white")
hdr_i <- which(tab$kind == "header")
sec_i <- which(tab$kind == "section")
for (ii in hdr_i) ft <- merge_at(ft, i = ii, j = 1:6, part = "body")
if (length(hdr_i)) {
  ft <- bg(ft, i = hdr_i, bg = "#D6E3F0", part = "body")
  ft <- bold(ft, i = hdr_i, part = "body")
}
for (ii in sec_i) ft <- merge_at(ft, i = ii, j = 1:6, part = "body")
if (length(sec_i)) {
  ft <- bg(ft, i = sec_i, bg = "#1F4E79", part = "body")
  ft <- color(ft, i = sec_i, color = "white", part = "body")
  ft <- bold(ft, i = sec_i, part = "body")
}
data_i <- which(tab$kind == "data")
if (length(data_i))
  ft <- compose(ft, i = data_i, j = "Variable", value = as_paragraph(""), part = "body")
ft <- width(ft, j = 1, width = 2.6)
ft <- width(ft, j = 2, width = 1.7)
ft <- width(ft, j = 3:6, width = 1.2)
ft <- footnote(ft, i = 1, j = 1, part = "header",
               value = as_paragraph(n_note), ref_symbols = "a")
doc <- read_docx()
doc <- body_add_par(doc, title, style = "heading 1")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "")
doc <- body_add_par(doc, "Terms retained by backward LRT", style = "heading 2")
ft2 <- flextable(sel_wide)
ft2 <- font(ft2, fontname = "Times New Roman", part = "all")
ft2 <- fontsize(ft2, size = 9, part = "all")
ft2 <- bold(ft2, part = "header")
ft2 <- bg(ft2, part = "header", bg = "#1F4E79")
ft2 <- color(ft2, part = "header", color = "white")
ft2 <- align(ft2, j = 2:5, align = "center", part = "all")
doc <- body_add_flextable(doc, ft2)
print(doc, target = docx_file)

message("Wrote ", xlsx_file)
message("Wrote ", docx_file)
print(sel_wide)
print(tab[tab$kind == "data", c("Variable", "Category",
                                 "HIV AOR (95% CI)", "HIV p",
                                 "HCV AOR (95% CI)", "HCV p")])
