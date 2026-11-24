# Imm-PathwayAge

Imm-PathwayAge estimates immune age from plasma proteins organised into biological pathways. It combines pathway-specific age predictions into an integrated clock, allowing immune ageing to be examined both overall and within individual pathways. The analyses are organised into seven modules within a single R script.

## Analysis modules

### 1. Immune pathway definition

Immune pathways are curated from Reactome v97, filtered for protein coverage and redundancy, and grouped into eight functional classes. The resulting pathway definitions provide the biological basis for model development and interpretation.

### 2. Clock development and evaluation

A two-stage model integrates pathway-specific age predictions into Imm-PathwayAge. Model selection uses five-fold outer and inner cross-validation within an 80:20 training–test split. Held-out prediction performance is summarised by MAE, RMSE and Pearson correlation, with confidence intervals from 2,000 bootstrap resamples.

### 3. Pathway contributions

Additive contributions quantify the role of each pathway in the integrated age estimate and their distribution across immune functions. Further analyses examine the relationship between protein count and pathway contribution, together with shared proteins among the ten leading pathways.

### 4. Disease associations

Logistic and Cox regression assess associations with prevalent and incident disease. UK Biobank age acceleration is defined by standardised residuals of full-cohort out-of-fold predicted age regressed on chronological age. Pathway-level associations, two-year landmark analyses and disease-category summaries characterise the patterns across clinical outcomes.

### 5. Mortality analyses

Cox models evaluate associations with all-cause mortality, including a two-year landmark analysis. Restricted cubic splines describe nonlinear associations, while ten-year absolute risk estimates aid interpretation. Disease–mortality comparisons examine similarities in the corresponding pathway association profiles.

### 6. All-protein comparison

A clock trained on all 2,920 proteins provides a comparator for Imm-PathwayAge. Separate disease models and mortality spline analyses compare their associations under the same outcome definitions and covariate adjustment.

### 7. External validation

The selected learner combination is refitted in UK Biobank using proteins shared with the Gong and Chia cohorts. Prediction performance and pathway contributions are evaluated in Gong participants and Chia healthy controls. Within Chia, participants with ALS and neurological controls are each compared with healthy controls after age acceleration is calibrated across all eligible participants.

## Software requirements

The script requires R 4.1 or later, with `data.table`, `mlr3verse`, `mlr3learners`, `mlr3extralearners`, `callr`, `survival` and `digest`. Model fitting also uses `glmnet`, `LiblineaR`, `lightgbm` and `FNN`; the all-protein comparator requires LightGBM 4.0 or later. Parallel processing uses the R `parallel` package, with optional support from `parallelly` and `RhpcBLASctl`.
