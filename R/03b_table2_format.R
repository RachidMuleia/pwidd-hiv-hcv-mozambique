################################################################################
## Format Table 2a (HIV) and Table 2b (HCV): merged domain rows, then variables
################################################################################

options(stringsAsFactors = FALSE)

## Project paths (repo-relative outputs; local microdata via 00_paths.R)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
r_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
source(file.path(r_dir, "00_paths.R"))

library(openxlsx)
library(officer)
library(flextable)

hiv <- read.csv(file.path(path_tab, "table2_hiv_gile_sspse.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
hcv <- read.csv(file.path(path_tab, "table2_hcv_gile_sspse.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
Ntab <- read.csv(file.path(path_tab, "table_sspse_N.csv"), stringsAsFactors = FALSE)
n_use <- function(year, city) {
  Ntab$N_used_for_pooling[Ntab$Year == year & Ntab$City == city][1]
}

split_char <- function(x) {
  x <- gsub(">=", "\u2265", x, fixed = TRUE)
  pos <- vapply(gregexpr(": ", x, fixed = TRUE), function(z) {
    if (all(z < 0)) NA_integer_ else max(z)
  }, integer(1))
  variable <- ifelse(is.na(pos), x, trimws(substr(x, 1, pos - 1)))
  category <- ifelse(is.na(pos), "", trimws(substr(x, pos + 2, nchar(x))))
  data.frame(Variable = variable, Category = category, stringsAsFactors = FALSE)
}

prep <- function(d) {
  sc <- split_char(d$Characteristic)
  data.frame(
    Domain = d$Domain,
    Variable = sc$Variable,
    Category = sc$Category,
    `2014 n` = d$`2014 n`,
    `2014 % (95% CI)` = d$`2014 % (95% CI)`,
    `2023 n` = d$`2023 n`,
    `2023 % (95% CI)` = d$`2023 % (95% CI)`,
    `Difference pp (95% CI)` = d$`Difference pp (95% CI)`,
    p = d$p,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

hiv_p <- prep(hiv)
hcv_p <- prep(hcv)

hdr <- c("Variable", "Category", "2014 n", "2014 % (95% CI)",
         "2023 n", "2023 % (95% CI)", "Difference pp (95% CI)", "p")

build_rows <- function(d) {
  out <- list()
  domains <- unique(d$Domain)
  for (dom in domains) {
    out[[length(out) + 1]] <- list(kind = "domain", text = dom)
    sub <- d[d$Domain == dom, ]
    for (i in seq_len(nrow(sub))) {
      show_var <- i == 1 || sub$Variable[i] != sub$Variable[i - 1]
      out[[length(out) + 1]] <- list(
        kind = "data",
        Variable = if (show_var) sub$Variable[i] else "",
        Category = sub$Category[i],
        `2014 n` = sub$`2014 n`[i],
        `2014 % (95% CI)` = sub$`2014 % (95% CI)`[i],
        `2023 n` = sub$`2023 n`[i],
        `2023 % (95% CI)` = sub$`2023 % (95% CI)`[i],
        `Difference pp (95% CI)` = sub$`Difference pp (95% CI)`[i],
        p = sub$p[i],
        Variable_full = sub$Variable[i]
      )
    }
  }
  out
}

n_note <- paste0(
  "Comparable sample: age \u2265 18 years, Maputo and Nampula/Nacala, injection in the previous 12 months. ",
  "Percentages are Gile successive-sampling estimates pooled by SS-PSE city population size. ",
  "N used: 2014 Maputo ", n_use("2014", "MAPUTO"),
  ", Nampula ", n_use("2014", "NAMPULA"),
  "; 2023 Maputo ", n_use("2023", "MAPUTO"),
  ", Nampula ", n_use("2023", "NAMPULA"),
  ". n is unweighted positive/tested. Difference is 2023 minus 2014, in percentage points."
)
title_hiv <- "Table 2a. Population-weighted HIV prevalence among PWID in Maputo and Nampula, 2014 and 2023"
title_hcv <- "Table 2b. Population-weighted HCV prevalence among PWID in Maputo and Nampula, 2014 and 2023"

font_name <- "Times New Roman"
title_style <- createStyle(fontName = font_name, fontSize = 12, textDecoration = "bold",
                           wrapText = TRUE, valign = "center")
sub_style <- createStyle(fontName = font_name, fontSize = 9, wrapText = TRUE,
                         valign = "center", fontColour = "#333333")
header_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            wrapText = TRUE, valign = "center", halign = "center",
                            fgFill = "#1F4E79", fontColour = "white",
                            border = "TopBottom", borderColour = "#1F4E79")
domain_style <- createStyle(fontName = font_name, fontSize = 10, textDecoration = "bold",
                            wrapText = TRUE, valign = "center", halign = "left",
                            fgFill = "#D6E3F0", fontColour = "#1F4E79",
                            border = "TopBottom", borderColour = "#8FAADC")
var_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         wrapText = TRUE, halign = "left")
cat_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         wrapText = TRUE, halign = "left", indent = 1)
num_style <- createStyle(fontName = font_name, fontSize = 10, valign = "center",
                         halign = "center")
note_style <- createStyle(fontName = font_name, fontSize = 8, wrapText = TRUE,
                          valign = "top", fontColour = "#333333")
bottom_style <- createStyle(border = "bottom", borderStyle = "medium",
                            borderColour = "#1F4E79")

write_formatted_sheet <- function(wb, sh, d, title, sub) {
  rows <- build_rows(d)
  ncol <- 8
  addWorksheet(wb, sh, gridLines = FALSE)
  pageSetup(wb, sh, orientation = "landscape", fitToWidth = TRUE, paperSize = 9,
            left = 0.4, right = 0.4, top = 0.5, bottom = 0.5)
  setColWidths(wb, sh, cols = 1:ncol, widths = c(32, 22, 12, 20, 12, 20, 22, 10))

  writeData(wb, sh, title, startRow = 1, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:ncol, rows = 1)
  addStyle(wb, sh, title_style, rows = 1, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = 1, heights = 26)

  writeData(wb, sh, sub, startRow = 2, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:ncol, rows = 2)
  addStyle(wb, sh, sub_style, rows = 2, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = 2, heights = 48)

  hdr_row <- 3
  for (j in seq_along(hdr)) {
    writeData(wb, sh, hdr[j], startRow = hdr_row, startCol = j, colNames = FALSE)
  }
  addStyle(wb, sh, header_style, rows = hdr_row, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = hdr_row, heights = 28)

  r <- hdr_row
  merge_starts <- list()
  for (item in rows) {
    r <- r + 1
    if (item$kind == "domain") {
      writeData(wb, sh, item$text, startRow = r, startCol = 1, colNames = FALSE)
      mergeCells(wb, sh, cols = 1:ncol, rows = r)
      addStyle(wb, sh, domain_style, rows = r, cols = 1:ncol, gridExpand = TRUE)
      setRowHeights(wb, sh, rows = r, heights = 20)
    } else {
      vals <- c(item$Variable, item$Category, item$`2014 n`, item$`2014 % (95% CI)`,
                item$`2023 n`, item$`2023 % (95% CI)`, item$`Difference pp (95% CI)`,
                item$p)
      for (j in seq_along(vals)) {
        writeData(wb, sh, vals[j], startRow = r, startCol = j, colNames = FALSE)
      }
      addStyle(wb, sh, var_style, rows = r, cols = 1)
      addStyle(wb, sh, cat_style, rows = r, cols = 2)
      addStyle(wb, sh, num_style, rows = r, cols = 3:ncol, gridExpand = TRUE)
    }
  }
  last <- r
  addStyle(wb, sh, bottom_style, rows = last, cols = 1:ncol, gridExpand = TRUE, stack = TRUE)

  ## vertical merge of Variable when the same name repeats
  data_rows <- integer()
  data_vars <- character()
  r <- hdr_row
  for (item in rows) {
    r <- r + 1
    if (item$kind == "data") {
      data_rows <- c(data_rows, r)
      data_vars <- c(data_vars, item$Variable_full)
    }
  }
  i <- 1
  while (i <= length(data_vars)) {
    j <- i
    while (j < length(data_vars) && data_vars[j + 1] == data_vars[i] &&
           data_rows[j + 1] == data_rows[j] + 1) {
      j <- j + 1
    }
    if (j > i) mergeCells(wb, sh, cols = 1, rows = data_rows[i]:data_rows[j])
    i <- j + 1
  }

  note_row <- last + 2
  writeData(wb, sh, paste0("a. ", sub), startRow = note_row, colNames = FALSE)
  mergeCells(wb, sh, cols = 1:ncol, rows = note_row)
  addStyle(wb, sh, note_style, rows = note_row, cols = 1:ncol, gridExpand = TRUE)
  setRowHeights(wb, sh, rows = note_row, heights = 42)
  freezePane(wb, sh, firstActiveRow = hdr_row + 1)
}

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 10, fontName = font_name)
write_formatted_sheet(wb, "Table 2a HIV", hiv_p, title_hiv, n_note)
write_formatted_sheet(wb, "Table 2b HCV", hcv_p, title_hcv, n_note)
xlsx_file <- file.path(path_out, "Table2_HIV_HCV_Gile_SSPSE.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

## ---- Word: grouped domain header, then variable / category ------------------

make_ft <- function(d) {
  g <- as_grouped_data(d, groups = "Domain")
  ft <- as_flextable(g, hide_grouplabel = TRUE)
  ft <- set_header_labels(
    ft,
    Variable = "Variable",
    Category = "Category",
    `2014 n` = "2014 n",
    `2014 % (95% CI)` = "2014 % (95% CI)",
    `2023 n` = "2023 n",
    `2023 % (95% CI)` = "2023 % (95% CI)",
    `Difference pp (95% CI)` = "Difference pp (95% CI)",
    p = "p"
  )
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- fontsize(ft, size = 8, part = "footer")
  ft <- bold(ft, part = "header")
  ft <- align(ft, j = 3:8, align = "center", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- merge_v(ft, j = "Variable")
  ft <- width(ft, j = "Variable", width = 2.5)
  ft <- width(ft, j = "Category", width = 1.6)
  ft <- width(ft, j = 3:8, width = 1.15)
  ft <- bg(ft, part = "header", bg = "#1F4E79")
  ft <- color(ft, part = "header", color = "white")
  ft <- bg(ft, i = ~ !is.na(Domain), bg = "#D6E3F0", part = "body")
  ft <- bold(ft, i = ~ !is.na(Domain), part = "body")
  ft <- color(ft, i = ~ !is.na(Domain), color = "#1F4E79", part = "body")
  ft <- italic(ft, j = "Category", italic = TRUE, part = "body")
  ft <- border_inner_h(ft, border = fp_border(color = "#C5D4E8", width = 0.4))
  ft <- border_outer(ft, border = fp_border(color = "#1F4E79", width = 0.9))
  ft <- padding(ft, padding.top = 3, padding.bottom = 3, part = "body")
  ft <- footnote(
    ft, i = 1, j = 1, part = "header",
    value = as_paragraph(n_note),
    ref_symbols = "a"
  )
  ft <- fontsize(ft, size = 8, part = "footer")
  ft
}

docx_file <- file.path(path_out, "Table2_HIV_HCV_Gile_SSPSE.docx")
doc <- read_docx()
doc <- body_add_par(doc, title_hiv, style = "heading 1")
doc <- body_add_flextable(doc, make_ft(hiv_p))
doc <- body_add_par(doc, "")
doc <- body_add_par(doc, title_hcv, style = "heading 1")
doc <- body_add_flextable(doc, make_ft(hcv_p))
print(doc, target = docx_file)

message("Wrote ", xlsx_file)
message("Wrote ", docx_file)
