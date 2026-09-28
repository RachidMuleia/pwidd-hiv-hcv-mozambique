# Data

This folder is intentionally empty in the public repository.

## Required inputs

| Survey | File | Role |
|--------|------|------|
| 2014 IBBS PWID | `PID_ELIGIBLE_ALL.csv` | Eligible 2014 respondents (Maputo, Nampula) |
| 2023 IBBS PWID | `DADOS_PID.csv` | 2023 respondents |

Point `R/00_paths.R` (or the `PWID_DATA_2014` / `PWID_DATA_2023` environment variables) to the local copies.

## Outputs written by the pipeline

Scripts write intermediate objects under `outputs/`, including:

- `analysis_objects.rds` — comparable sample and RDS structures
- `sspse_city_fits.rds` — city population-size fits used for Gile SS
- `tables/` — CSV and Excel publication tables

Do not commit microdata or identifiable recruitment files.
