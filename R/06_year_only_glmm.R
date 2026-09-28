################################################################################
## Year-only association models: Maputo, Nampula, and overall
## y ~ YEAR  (no other covariates)
## Logistic GLMM with recruiter random intercept nested in year x city
## Outcomes: HIV and HCV
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
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f\u2013%.2f)", or, lo, hi)

obj <- readRDS(file.path(path_out, "analysis_objects.rds"))
dat <- obj$dat_comp

make_cluster <- function(year, city, recruiter, coupon) {
  rec <- as.character(recruiter)
  coup <- as.character(coupon)
  bad <- is.na(rec) | rec %in% c("", "NA", "seed", "SEED", "0")
  rec[bad] <- coup[bad]
  paste(year, city, rec, sep = ":")
}

ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

prep <- function(dat, outcome, city = NULL) {
  md <- dat
  if (!is.null(city)) md <- md[md$CIDADE == city, ]
  md$YEAR <- factor(md$YEAR, levels = c("2014", "2023"))
  md$y <- as.numeric(md[[outcome]])
  md$cluster_id <- factor(make_cluster(md$YEAR, md$CIDADE, md$recruiter_id, md$coupon_id))
  md <- md[complete.cases(md[, c("y", "YEAR", "cluster_id")]), ]
  md
}

fit_year <- function(md, label, outcome) {
  ## Prefer GLMM; fall back to glm if RE variance collapses / singular
  fit <- tryCatch(
    lme4::glmer(y ~ YEAR + (1 | cluster_id), data = md,
                family = binomial(link = "logit"), nAGQ = 0, control = ctrl),
    error = function(e) NULL
  )
  method <- "GLMM"
  if (is.null(fit) || isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
    fit <- stats::glm(y ~ YEAR, data = md, family = binomial(link = "logit"))
    method <- "GLM (singular RE)"
    cf <- coef(fit)
    V <- vcov(fit)
    se <- sqrt(diag(V))["YEAR2023"]
    b <- unname(cf["YEAR2023"])
    sigma_u <- 0
  } else {
    cf <- lme4::fixef(fit)
    V <- as.matrix(vcov(fit))
    se <- sqrt(V["YEAR2023", "YEAR2023"])
    b <- unname(cf["YEAR2023"])
    sigma_u <- as.data.frame(lme4::VarCorr(fit))$sdcor[1]
  }
  or <- exp(b)
  lo <- exp(b - 1.96 * se)
  hi <- exp(b + 1.96 * se)
  p <- 2 * pnorm(-abs(b / se))
  data.frame(
    Outcome = outcome,
    Sample = label,
    n = nrow(md),
    n_2014 = sum(md$YEAR == "2014"),
    n_2023 = sum(md$YEAR == "2023"),
    events = sum(md$y == 1),
    events_2014 = sum(md$y[md$YEAR == "2014"] == 1),
    events_2023 = sum(md$y[md$YEAR == "2023"] == 1),
    AOR = or,
    `AOR (95% CI)` = fmt_or(or, lo, hi),
    `p-value` = fmt_p(p),
    `sd(u)` = round(sigma_u, 3),
    Method = method,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

run_outcome <- function(outcome) {
  rbind(
    fit_year(prep(dat, outcome, "MAPUTO"), "Maputo", outcome),
    fit_year(prep(dat, outcome, "NAMPULA"), "Nampula/Nacala", outcome),
    fit_year(prep(dat, outcome, NULL), "Overall", outcome)
  )
}

message("HIV ...")
hiv <- run_outcome("HIV")
message("HCV ...")
hcv <- run_outcome("HCV")
tab <- rbind(hiv, hcv)
print(tab[, c("Outcome", "Sample", "n", "AOR (95% CI)", "p-value", "Method")])

footnote <- paste0(
  "Comparable sample: age \u2265 18 years; injection in the previous 12 months. ",
  "Unadjusted logistic mixed model: outcome ~ year (2023 vs 2014), no other covariates. ",
  "Random intercept: recruiter nested in year \u00d7 city (ordinary logistic used if the random effect is singular). ",
  "AOR, adjusted odds ratio for survey year (only year is in the model). Wald 95% CI."
)

wb <- createWorkbook()
font_name <- "Times New Roman"
title_style <- createStyle(fontName = font_name, fontSize = 12, textDecoration = "bold",
                           wrapText = TRUE)
header_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            halign = "center", valign = "center",
                            border = "TopBottom", borderColour = "black")
body_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center")
num_style <- createStyle(fontName = font_name, fontSize = 10, halign = "center",
                         valign = "center")
note_style <- createStyle(fontName = font_name, fontSize = 8, wrapText = TRUE)
section_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold")

write_outcome_sheet <- function(wb, sh, d, outcome) {
  addWorksheet(wb, sh, gridLines = FALSE)
  writeData(wb, sh,
            paste0("Table. Year effect (2023 vs 2014) on ", outcome,
                   ": Maputo, Nampula/Nacala, and overall"),
            startRow = 1, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:8, rows = 1)
  addStyle(wb, sh, title_style, rows = 1, cols = 1:8, gridExpand = TRUE)

  out <- data.frame(
    Sample = d$Sample,
    n = d$n,
    `n 2014` = d$n_2014,
    `n 2023` = d$n_2023,
    Events = paste0(d$events, " (", d$events_2014, "/", d$events_2023, ")"),
    `AOR (95% CI)` = d$`AOR (95% CI)`,
    `p-value` = d$`p-value`,
    Method = d$Method,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  writeData(wb, sh, out, startRow = 3)
  addStyle(wb, sh, header_style, rows = 3, cols = 1:8, gridExpand = TRUE)
  addStyle(wb, sh, body_style, rows = 4:(3 + nrow(out)), cols = 1:8, gridExpand = TRUE)
  addStyle(wb, sh, num_style, rows = 4:(3 + nrow(out)), cols = 2:7, gridExpand = TRUE)
  setColWidths(wb, sh, cols = 1:8, widths = c(16, 8, 10, 10, 16, 20, 10, 18))
  writeData(wb, sh, paste0("Footnote. ", footnote), startRow = 5 + nrow(out),
            colNames = FALSE)
  mergeCells(wb, sh, cols = 1:8, rows = 5 + nrow(out))
  addStyle(wb, sh, note_style, rows = 5 + nrow(out), cols = 1:8, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = 5 + nrow(out), heights = 48)
}

write_outcome_sheet(wb, "HIV", hiv, "HIV")
write_outcome_sheet(wb, "HCV", hcv, "HCV")

addWorksheet(wb, "Combined")
writeData(wb, "Combined",
          "Year-only models (2023 vs 2014), no covariates",
          startRow = 1, colNames = FALSE)
writeData(wb, "Combined", tab, startRow = 3)
addStyle(wb, "Combined", header_style, rows = 3, cols = 1:ncol(tab), gridExpand = TRUE)
setColWidths(wb, "Combined", cols = 1:ncol(tab), widths = c(8, 14, 8, 10, 10, 8, 10, 10, 18, 10, 8, 16))

xlsx_file <- file.path(path_out, "Table_year_only_GLMM.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
write.csv(tab, file.path(path_tab, "table_year_only_glmm.csv"), row.names = FALSE)
message("Wrote ", xlsx_file)
