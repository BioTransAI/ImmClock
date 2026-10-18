# Imm-PathwayAge

**Manuscript analysis code** · R · Single integrated script

Analyses correspond to the current manuscript and its supplementary figures
and tables. All results remain in R objects; visualization and publication-file
exports are excluded.

## Analysis map

| Module | Retained analyses | Manuscript location |
|:--|:--|:--|
| **Immune pathways** | Reactome v97 selection, redundancy filtering and eight-class annotation | Figure 1; Table S2 |
| **Clock development** | Two-stage models, nested cross-validation, held-out evaluation and full-cohort cross-fitting | Figure 2; Table S3 |
| **Pathway contributions** | Additive contributions, category allocation, protein-count correlation and shared proteins among the ten leading pathways | Figure 3; Figure S2; Table S4 |
| **Disease associations** | Prevalent and incident disease, pathway follow-up, two-year landmark, category composition and association breadth | Figures 4–5; Tables S5–S10 |
| **Mortality** | Cox models, two-year landmark, restricted cubic splines, ten-year absolute risks and disease–mortality pathway profiles | Figure 6; Table S11 |
| **All-protein comparator** | A 2,920-protein clock; separate disease models and mortality spline comparisons | Figures S3–S5 |
| **External validation** | Fixed-model prediction, pathway contributions and Chia neurological-disease associations | Figure 7; Tables S12–S15 |

## Use

Load `Imm_PathwayAge_Analysis.R` into the main R session, then call
`run_imm_pathway_analysis()` with the original study objects. Argument comments
describe their required fields. Loading the script does not start an analysis.
For individual module calls, first run `load_imm_pathway_packages()`.

Inputs include the preprocessed protein-and-age matrix, prepared covariates,
collection and outcome dates, and Gong/Chia cohort objects. Protein KNN
imputation and median covariate imputation are upstream requirements; their
implementations are not reconstructed here. Previously fitted clock and
Reactome objects may be reused. Otherwise, provide a reference-cache directory
for the public Reactome v97 inputs. No personal file locations are embedded.

## Analysis conventions

The main clock retains the original 80:20 split, five-fold outer and inner
cross-validation, and 2,000-resample bootstrap confidence intervals for prediction
metrics. UK Biobank
associations use full-cohort out-of-fold age acceleration and the manuscript's
covariates, endpoints and multiple-testing families. External refitting uses
the already selected learner combination; prediction is evaluated in Gong and
Chia controls. Chia age correction pools all eligible participants before the
two separate disease comparisons.

Figure S1 compares literature-reported correlations and requires no additional
participant-level model fitting. Study dimensions and numerical integrity
checks are retained throughout the script.

## Software

R syntax requires version 4.1 or later. Packages: `data.table`, `mlr3verse`,
`mlr3learners`, `mlr3extralearners`, `callr`, `survival`, `digest`, and the original
learner backends `glmnet`, `LiblineaR`, `lightgbm` and `FNN`. The comparator uses
LightGBM 4.0 or later. `parallel` is included with R; `parallelly` and
`RhpcBLASctl` are optional.

The returned object contains the analysis results and session information.
Reproduction requires the original study inputs and computational environment.
