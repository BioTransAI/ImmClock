# Imm-PathwayAge

Imm-PathwayAge estimates immune age from plasma proteins organised into biological pathways. It combines pathway-specific age predictions into an integrated clock, while retaining pathway-level age estimates. This repository contains the R analyses for model development, associations with disease and mortality, and validation in independent cohorts.

## Analysis overview

The workflow begins with the definition of immune pathways from Reactome v97 and their annotation into eight functional classes. A two-stage model combines the resulting pathway clocks. Model development uses an 80:20 training–test split and five-fold outer and inner cross-validation; prediction metrics are accompanied by confidence intervals from 2,000 bootstrap resamples.

Pathway contributions describe how the component clocks contribute to the integrated age estimate and how these contributions are distributed across immune functions. UK Biobank associations use age acceleration derived from full-cohort out-of-fold predictions. Logistic and Cox regression relate age acceleration to prevalent disease, incident disease and mortality. Two-year landmark analyses, restricted cubic splines, ten-year absolute risk estimates and disease–mortality pathway profiles provide further characterisation of these associations. A clock trained on all 2,920 proteins serves as a separate comparator.

For external validation, the selected learner combination is refitted in UK Biobank using proteins shared with the Gong and Chia cohorts. Prediction performance and pathway contributions are evaluated in Gong participants and Chia healthy controls. Participants with ALS and neurological controls are each compared with healthy controls, using age acceleration calibrated across all eligible Chia participants.

## Using the code

All analyses are contained in `Imm_PathwayAge_Analysis.R`. The workflow starts from prepared protein measurements, covariates, collection and outcome dates, and external cohort data. Protein KNN imputation and median covariate imputation should follow the study preprocessing. Required input fields are documented alongside the main function.

The example below assumes that these inputs and the Reactome pathway bundle have been loaded into R.

```r
source("Imm_PathwayAge_Analysis.R")

analysis <- run_imm_pathway_analysis(
  all_data = all_data,
  covariates = covariates,
  blood = blood,
  outcome_dates = outcome_dates,
  external_cohorts = external_cohorts,
  bundle = bundle
)
```

To construct the pathway bundle within the workflow, omit `bundle` and supply a `reactome_cache_dir`. A previously fitted clock can be supplied through `result`. Analysis settings are defined in `imm_pathway_config()`.

The returned `analysis` object contains model results, pathway annotations and contributions, clinical associations, external validation results and session information. Results remain available in R for inspection and further use.

## Software requirements

The script requires R 4.1 or later, with `data.table`, `mlr3verse`, `mlr3learners`, `mlr3extralearners`, `callr`, `survival` and `digest`. Model fitting also uses `glmnet`, `LiblineaR`, `lightgbm` and `FNN`; the all-protein comparator requires LightGBM 4.0 or later. Parallel processing uses the R `parallel` package, with optional support from `parallelly` and `RhpcBLASctl`.
