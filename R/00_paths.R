################################################################################
## Project paths
## Edit PATH_DATA_2014 / PATH_DATA_2023 if your local data live elsewhere.
################################################################################

## Repository root = parent of R/
.path_r <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  ## Interactive fallback: assume working directory is repo root or R/
  wd <- normalizePath(getwd())
  if (basename(wd) == "R") return(wd)
  file.path(wd, "R")
}

PATH_R    <- .path_r()
PATH_ROOT <- dirname(PATH_R)
PATH_OUT  <- file.path(PATH_ROOT, "outputs")
PATH_TAB  <- file.path(PATH_OUT, "tables")
PATH_FIG  <- file.path(PATH_OUT, "figures")
PATH_TEX  <- file.path(PATH_OUT, "tex")

## Local microdata (not shipped with the repository)
## Default: sibling Dropbox layout used for the manuscript revision.
PATH_DATA_2014 <- Sys.getenv(
  "PWID_DATA_2014",
  unset = "/Users/rachidmuleia/Dropbox/INS/PID/ARTIGO_PRINCIPAL_PID/DADOS/PID_ELIGIBLE_ALL.csv"
)
PATH_DATA_2023 <- Sys.getenv(
  "PWID_DATA_2023",
  unset = "/Users/rachidmuleia/Dropbox/INS/PID/ARTIGO_AURIA_SSR/DADOS_PID.csv"
)

## Optional: raw 2014 folder when a script needs the same CSV via path_root/DADOS
PATH_DATA_DIR_2014 <- dirname(PATH_DATA_2014)

for (p in c(PATH_OUT, PATH_TAB, PATH_FIG, PATH_TEX)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

## Compatibility aliases used by the analysis scripts
path_root <- PATH_ROOT
path_out  <- PATH_OUT
path_tab  <- PATH_TAB
path_fig  <- PATH_FIG
path_tex  <- PATH_TEX
path_2014 <- PATH_DATA_2014
path_2023 <- PATH_DATA_2023
