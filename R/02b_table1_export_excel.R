
## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))
################################################################################
## Publication Excel for comparable Table 1 (HIV, HCV Gile SS; RDS-II composition)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

hiv <- read.csv(file.path(path_tab, "table1_hiv_pooled_sspse.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
hcv <- read.csv(file.path(path_tab, "table1_hcv_pooled_sspse.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
long <- read.csv(file.path(path_tab, "table1_pooled_sspse.csv"),
                 stringsAsFactors = FALSE)
Ntab <- read.csv(file.path(path_tab, "table_sspse_N.csv"),
                 stringsAsFactors = FALSE)
comp_city <- read.csv(file.path(path_tab, "table1_rdsii_main.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
comp_pool <- read.csv(file.path(path_tab, "table1_rdsii_pooled.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
n_cri <- function(year, city) {
  r <- Ntab[Ntab$Year == year & Ntab$City == city, ]
  sprintf("%s (prior median %s; 95%% CrI %s)",
          r$N_used_for_pooling[1], r$Prior_median[1], r$CrI_95[1])
}

long_out <- long |>
  filter(outcome %in% c("HIV", "HCV")) |>
  transmute(
    Outcome = outcome,
    Domain = domain,
    Characteristic = label,
    Category = category_lab,
    n_positive_2014 = n_pos_2014,
    n_tested_2014 = n_2014,
    n_2014 = n_2014_lab,
    prevalence_2014 = p_2014,
    se_2014 = se_2014,
    ci_lo_2014 = lo_2014,
    ci_hi_2014 = hi_2014,
    prevalence_2014_pct_CI = prev_2014,
    n_positive_2023 = n_pos_2023,
    n_tested_2023 = n_2023,
    n_2023 = n_2023_lab,
    prevalence_2023 = p_2023,
    se_2023 = se_2023,
    ci_lo_2023 = lo_2023,
    ci_hi_2023 = hi_2023,
    prevalence_2023_pct_CI = prev_2023,
    difference = diff,
    se_difference = se_diff,
    difference_lo = diff_lo,
    difference_hi = diff_hi,
    difference_pp_CI = diff_pp,
    p_value = p_value,
    p = p_lab
  )

notes <- data.frame(
  Item = c(
    "Manuscript",
    "Sample",
    "Outcomes",
    "Estimator",
    "Pooling",
    "N_c 2014 Maputo",
    "N_c 2014 Nampula",
    "N_c 2023 Maputo",
    "N_c 2023 Nampula",
    "2014 recruiter graph",
    "Difference",
    "Kept characteristics",
    "Omitted (not comparable)",
    "Relabel",
    "Footnote",
    "Scripts"
  ),
  Detail = c(
    "PONE-D-25-62212; PWID_AB V2.docx Table 1",
    "Maputo and Nampula; age >= 18; injection in the previous 12 months. Beira, Tete, Quelimane and ages 16-17 excluded.",
    "HIV laboratory result; HCV antibody. Missing laboratory result is not coded as positive.",
    "City-specific Gile successive-sampling (SS) estimator with N = SS-PSE posterior median; RDS compute.weights weight.type = Gile SS. 300 within-city bootstrap replicates of the Gile-weighted Hajek mean.",
    "Two cities pooled by SS-PSE posterior median of the PWID population (sspse::posteriorsize).",
    n_cri("2014", "MAPUTO"),
    n_cri("2014", "NAMPULA"),
    n_cri("2023", "MAPUTO"),
    n_cri("2023", "NAMPULA"),
    "Reconstructed by matching each 2014 coupon to outgoing coupons (OUTPON1-5) in PID_ELIGIBLE_ALL.csv. Seeds and recruits whose recruiter is not in the comparable sample are seeds in the RDS graph.",
    "2023 minus 2014, in percentage points. Wald p from independent-year SEs. p < 0.001 printed as <0.001, never as 0.",
    "Overall; city; age; sex; education; marital status; age at first injection; injection frequency; access to new syringes; number of sexual partners past 12 months; unprotected sex; HIV test past 12 months; sought health care.",
    "Shared syringes; ever used a previously used syringe; ever shared other injection equipment; transactional sex; self-reported STI; HIV risk perception.",
    "Manuscript 'HIV test in the last month' is a past-12-month item in both questionnaires.",
    "Access to new syringes: wording close across rounds, not proven identical. Female 2014 cells are small (HIV n=19; HCV n=16).",
    "PID_RDS_COMPARABLE_ANALYSIS.R; PID_TABLE1_POOLED_SSPSE.R; PID_TABLE1_RDS.R; PID_TABLE1_EXPORT_EXCEL.R"
  ),
  stringsAsFactors = FALSE
)

hiv_all <- long[long$outcome == "HIV" & long$variable == "OVERALL", ][1, ]
hcv_all <- long[long$outcome == "HCV" & long$variable == "OVERALL", ][1, ]

font_name <- "Times New Roman"
header_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            wrapText = TRUE, valign = "center", halign = "center",
                            fgFill = "#F2F2F2", border = "TopBottom", borderStyle = "medium")
header_left <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                           wrapText = TRUE, valign = "center", halign = "left",
                           fgFill = "#F2F2F2", border = "TopBottom", borderStyle = "medium")
title_style <- createStyle(fontName = font_name, fontSize = 11, textDecoration = "bold",
                           wrapText = TRUE, valign = "center")
sub_style <- createStyle(fontName = font_name, fontSize = 9, wrapText = TRUE,
                         valign = "center", fontColour = "#404040")
body_left <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         halign = "left")
body_ctr <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                        halign = "center")
note_style <- createStyle(fontName = font_name, fontSize = 9, wrapText = TRUE,
                          valign = "top", halign = "left")
wrap_style <- createStyle(fontName = font_name, fontSize = 10, wrapText = TRUE,
                          valign = "center")
pct_style <- createStyle(fontName = font_name, numFmt = "0.0%")
num_style <- createStyle(fontName = font_name, numFmt = "0.000")
p_style <- createStyle(fontName = font_name, numFmt = "0.000")

style_prev_sheet <- function(wb, sh, df, title, sub, notes_txt, col_widths) {
  ncol <- ncol(df)
  n <- nrow(df)
  hdr <- 3
  first <- hdr + 1
  last <- hdr + n
  showGridLines(wb, sh, showGridLines = FALSE)
  pageSetup(wb, sh, orientation = "landscape", fitToWidth = TRUE, paperSize = 9,
            left = 0.4, right = 0.4, top = 0.5, bottom = 0.5)
  setColWidths(wb, sh, cols = seq_len(ncol), widths = col_widths)
  writeData(wb, sh, title, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:ncol, rows = 1)
  addStyle(wb, sh, title_style, rows = 1, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = 1, heights = 28)
  writeData(wb, sh, sub, startRow = 2, startCol = 1, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:ncol, rows = 2)
  addStyle(wb, sh, sub_style, rows = 2, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = 2, heights = 32)
  hdr_names <- names(df)
  for (j in seq_along(hdr_names)) {
    writeData(wb, sh, hdr_names[j], startRow = hdr, startCol = j, colNames = FALSE)
  }
  addStyle(wb, sh, header_left, rows = hdr, cols = 1:3)
  addStyle(wb, sh, header_style, rows = hdr, cols = 4:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = hdr, heights = 28)
  writeData(wb, sh, df, startRow = first, startCol = 1, colNames = FALSE)
  addStyle(wb, sh, body_left, rows = first:last, cols = 1:3, gridExpand = TRUE)
  addStyle(wb, sh, body_ctr, rows = first:last, cols = 4:ncol, gridExpand = TRUE)
  addStyle(wb, sh, createStyle(fontName = font_name, fontSize = 10, border = "bottom",
                               borderStyle = "medium"),
           rows = last, cols = 1:ncol, gridExpand = TRUE, stack = TRUE)
  note_start <- last + 2
  for (k in seq_along(notes_txt)) {
    rr <- note_start + k - 1
    writeData(wb, sh, paste0(letters[k], ". ", notes_txt[k]),
              startRow = rr, startCol = 1, colNames = FALSE)
    mergeCells(wb, sh, cols = 1:ncol, rows = rr)
    addStyle(wb, sh, note_style, rows = rr, cols = 1:ncol, gridExpand = TRUE)
    setRowHeights(wb, sh, rows = rr, heights = 36)
  }
  freezePane(wb, sh, firstActiveRow = first)
}

hiv_notes <- c(
  paste("n is HIV-positive/tested (unweighted). Percentages are Gile successive-sampling estimates pooled by SS-PSE city size:",
        "2014 Maputo", n_cri("2014", "MAPUTO"), "; Nampula", n_cri("2014", "NAMPULA"),
        "; 2023 Maputo", n_cri("2023", "MAPUTO"), "; Nampula", n_cri("2023", "NAMPULA"), "."),
  "2014 recruiter IDs were reconstructed from outgoing coupons (OUTPON1-5). Difference is 2023 minus 2014, in percentage points.",
  "Female 2014 n = 19; the interval is wide. Access to new syringes: wording is close across rounds, not proven identical. HIV test is a past-12-month item in both questionnaires. Unprotected sex: 2014 is unprotected intercourse in the past 12 months; 2023 is no condom at last sex (LASTREL6 recoded)."
)
hcv_notes <- c(
  "Same sample, estimator and N_c as the HIV table. n is HCV-antibody-positive/tested (unweighted).",
  paste0("The pooled HCV difference (", hcv_all$diff_pp, ", p = ", hcv_all$p_lab,
         ") is a cancellation across cities: HCV fell in Maputo and rose in Nampula."),
  "Female 2014 n = 16. Access to new syringes and HIV testing footnotes as in the HIV table."
)

wb <- createWorkbook()
options("openxlsx.dateFormat" = "yyyy-mm-dd")
modifyBaseFont(wb, fontSize = 10, fontName = font_name)

addWorksheet(wb, "HIV", gridLines = FALSE)
style_prev_sheet(
  wb, "HIV", hiv,
  title = "Table 1a. HIV prevalence among people who inject drugs, Maputo and Nampula, 2014 and 2023.",
  sub = sprintf("Gile successive-sampling estimates pooled by SS-PSE city size. Complete HIV results: %s (2014) and %s (2023).",
                hiv_all$n_2014_lab, hiv_all$n_2023_lab),
  notes_txt = hiv_notes,
  col_widths = c(18, 28, 18, 12, 26, 12, 26, 24, 10)
)

addWorksheet(wb, "HCV", gridLines = FALSE)
style_prev_sheet(
  wb, "HCV", hcv,
  title = "Table 1b. HCV-antibody prevalence among people who inject drugs, Maputo and Nampula, 2014 and 2023.",
  sub = sprintf("Same sample, estimator and N_c as Table 1a. Complete HCV-antibody results: %s (2014) and %s (2023).",
                hcv_all$n_2014_lab, hcv_all$n_2023_lab),
  notes_txt = hcv_notes,
  col_widths = c(18, 28, 18, 12, 26, 12, 26, 24, 10)
)

addWorksheet(wb, "Composition_RDSII")
writeData(wb, "Composition_RDSII",
          "RDS-II (Volz-Heckathorn) weighted sample composition, Maputo and Nampula. Degree-based; not Gile SS.",
          startRow = 1, colNames = FALSE)
mergeCells(wb, "Composition_RDSII", cols = 1:ncol(comp_city), rows = 1)
addStyle(wb, "Composition_RDSII", title_style, rows = 1, cols = 1:ncol(comp_city), gridExpand = TRUE)
writeData(wb, "Composition_RDSII", comp_city, startRow = 3)
addStyle(wb, "Composition_RDSII", header_style, rows = 3, cols = 1:ncol(comp_city), gridExpand = TRUE)
setColWidths(wb, "Composition_RDSII", cols = 1:ncol(comp_city), widths = 22)
freezePane(wb, "Composition_RDSII", firstActiveRow = 4)

addWorksheet(wb, "Composition_pooled")
writeData(wb, "Composition_pooled",
          "RDS-II composition pooled by published/mapping PWID population mix (secondary).",
          startRow = 1, colNames = FALSE)
mergeCells(wb, "Composition_pooled", cols = 1:ncol(comp_pool), rows = 1)
addStyle(wb, "Composition_pooled", title_style, rows = 1, cols = 1:ncol(comp_pool), gridExpand = TRUE)
writeData(wb, "Composition_pooled", comp_pool, startRow = 3)
addStyle(wb, "Composition_pooled", header_style, rows = 3, cols = 1:ncol(comp_pool), gridExpand = TRUE)
setColWidths(wb, "Composition_pooled", cols = 1:ncol(comp_pool), widths = 22)
freezePane(wb, "Composition_pooled", firstActiveRow = 4)

addWorksheet(wb, "Estimates_numeric")
writeData(wb, "Estimates_numeric", long_out)
setColWidths(wb, "Estimates_numeric", cols = 1:ncol(long_out), widths = 16)
addStyle(wb, "Estimates_numeric", header_style, rows = 1, cols = 1:ncol(long_out), gridExpand = TRUE)
addStyle(wb, "Estimates_numeric", pct_style,
         rows = 2:(nrow(long_out) + 1),
         cols = which(names(long_out) %in% c(
           "prevalence_2014", "ci_lo_2014", "ci_hi_2014",
           "prevalence_2023", "ci_lo_2023", "ci_hi_2023",
           "difference", "difference_lo", "difference_hi"
         )),
         gridExpand = TRUE)
addStyle(wb, "Estimates_numeric", num_style,
         rows = 2:(nrow(long_out) + 1),
         cols = which(names(long_out) %in% c("se_2014", "se_2023", "se_difference")),
         gridExpand = TRUE)
addStyle(wb, "Estimates_numeric", p_style,
         rows = 2:(nrow(long_out) + 1),
         cols = which(names(long_out) == "p_value"),
         gridExpand = TRUE)
freezePane(wb, "Estimates_numeric", firstRow = TRUE)

addWorksheet(wb, "Population_size_SSPSE")
writeData(wb, "Population_size_SSPSE", Ntab)
setColWidths(wb, "Population_size_SSPSE", cols = 1:ncol(Ntab), widths = 22)
addStyle(wb, "Population_size_SSPSE", header_style, rows = 1, cols = 1:ncol(Ntab), gridExpand = TRUE)

addWorksheet(wb, "Notes")
writeData(wb, "Notes", notes)
setColWidths(wb, "Notes", cols = 1, widths = 28)
setColWidths(wb, "Notes", cols = 2, widths = 110)
addStyle(wb, "Notes", header_style, rows = 1, cols = 1:2)
addStyle(wb, "Notes", wrap_style, rows = 2:(nrow(notes) + 1), cols = 1:2,
         gridExpand = TRUE)
setRowHeights(wb, "Notes", rows = 2:(nrow(notes) + 1), heights = 36)

out_xlsx <- file.path(path_tab, "Table1_HIV_HCV_comparable.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Wrote ", out_xlsx)
