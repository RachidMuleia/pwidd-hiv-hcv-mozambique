# HIV and HCV among people who inject drugs in Mozambique

Analysis code for the serial cross-sectional bio-behavioural surveys (BBS) among people who inject drugs (PWID) in Maputo and Nampula/Nacala, Mozambique (2014 and 2023).

Comparable sample: age ≥ 18 years; injection in the previous 12 months; Maputo and Nampula/Nacala only.

## Repository layout

```
R/
  00_paths.R                 # Paths and local data locations
  01_comparable_sample_rds.R # Harmonised sample, RDS objects, prevalence pipeline
  02_table1_composition.R    # Table 1: Gile's weighted sample composition
  02b_table1_export_excel.R  # Table 1 Excel export
  03_table2_prevalence.R     # Table 2: HIV/HCV prevalence (Gile SS + SS-PSE N)
  03b_table2_format.R        # Table 2 formatting / export
  04_table4_glmm_lrt.R       # Table 4: GLMM + backward LRT selection
  05_table4_glmm_final.R     # Table 4: final publish models (+ RE LRT)
  06_year_only_glmm.R        # Unadjusted year effect (Maputo, Nampula, overall)
data/                        # Place microdata here (not committed)
outputs/                     # Tables, figures, RDS objects written here
```

## Analysis pipeline

Run scripts in order from the repository root (or from `R/`):

```bash
Rscript R/01_comparable_sample_rds.R
Rscript R/02_table1_composition.R
Rscript R/02b_table1_export_excel.R   # optional Excel polish
Rscript R/03_table2_prevalence.R
Rscript R/04_table4_glmm_lrt.R
Rscript R/05_table4_glmm_final.R
Rscript R/06_year_only_glmm.R
```

| Step | What it does |
|------|----------------|
| 01 | Builds the comparable sample, RDS weights/objects (`analysis_objects.rds`), city prevalence |
| 02 | City-stratified RDS-II composition (Table 1) |
| 03 | Population-weighted HIV/HCV prevalence (Table 2) |
| 04 | Unweighted logistic GLMM with year × covariate interactions; backward LRT |
| 05 | Final Table 4 models; forces year × city for HIV; LRT for recruiter random intercept |
| 06 | Year-only GLMM (no covariates) by city and overall |

## Methods summary

- **Descriptives (Tables 1–2):** RDS-II (Volz–Heckathorn) / Gile successive sampling; city-stratified; not raw sample proportions.
- **Associations (Table 4):** Unweighted logistic mixed models (`lme4::glmer`), random intercept for recruiter nested in year × city; year × covariate interactions; backward LRT (α = 0.05); boundary-corrected LRT for \(H_0:\sigma_u=0\).

## Data (not included)

Survey microdata are confidential and are **not** in this repository. See [`data/README.md`](data/README.md).

Set paths via environment variables or edit `R/00_paths.R`:

```bash
export PWID_DATA_2014="/path/to/PID_ELIGIBLE_ALL.csv"
export PWID_DATA_2023="/path/to/DADOS_PID.csv"
```

## R packages

`tidyverse`, `RDS`, `survey`, `lme4`, `openxlsx`, `xtable`, `ggplot2`, `scales`, `broom`, `glue`, `officer`, `flextable` (and dependencies). Optional older drafts also used `geepack`, `sandwich`, `lmtest`.


