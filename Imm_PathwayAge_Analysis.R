# Imm-PathwayAge: integrated manuscript analysis
# English analysis code aligned with the V8 manuscript.
# Load the definitions, then call run_imm_pathway_analysis().
# Analyses return R objects; plotting and publication-file exports are excluded.


# ==============================================================================
# Workflow and configuration
# ==============================================================================

# Unified analysis workflow
# Supply the original preprocessed study objects. Cohort exclusions, protein
# preprocessing and linkage performed upstream are not reimplemented here.
# Protein KNN imputation and median covariate imputation described in the
# manuscript must already be represented in the supplied study objects.
# No participant data, local filenames or personal directories are embedded.
imm_pathway_config <- function() {
  list(seed = 666L, training_fraction = 0.8, outer_folds = 5L, inner_folds = 5L, model_cores = 30L, fit_timeout =
    600, retain_models = FALSE, clinical_cores = 30L, external_cores = 1L, bootstrap = 2000L, ph_risk_bootstrap =
    200L, ph_risk_seed = 20260910L, expected_test_n = 10602L, minimum_events = 200L, alpha = 0.05, disease_censor =
    "2022-05-31", death_censor = "2024-08-31", comparator = all_protein_config())
}
load_imm_pathway_packages <- function() {
  packages <- c("data.table", "mlr3verse", "mlr3learners", "mlr3extralearners", "survival", "digest", "callr",
    "parallel")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required R packages: ", paste(missing, collapse = ", "))
  for (package in packages) suppressPackageStartupMessages(library(package, character.only = TRUE))
  invisible(TRUE)
}
imm_pathway_learners <- function() {
  list(stage1 = list(Ridge = mlr3::lrn("regr.glmnet", alpha = 0), LinearSVR = mlr3::lrn("regr.liblinear"),
    LightGBM = mlr3::lrn("regr.lightgbm"), FNN = mlr3::lrn("regr.fnn"), LASSO = mlr3::lrn("regr.glmnet", alpha =
    1)), stage2 = list(Ridge = mlr3::lrn("regr.glmnet", alpha = 0), LM = mlr3::lrn("regr.lm"), NotNegRidge =
    mlr3::lrn("regr.glmnet", alpha = 0, lower.limits = 0)), fusion = c("EqualWeight", "MAEWeight"))
}
imm_pathway_metadata <- function(result, bundle, performance) {
  ann <- ia_reactome8_annotation(result, bundle)
  d <- performance$contribution
  list(annotation = ann, external_annotation = data.frame(reactome_id = ann$id, pathway_name = trimws(ann$pathway),
    class_id = ann$primary_class, class_name = ann$primary_class_name), external_contribution = data.frame(reactome_id =
    d$id, n = d$N, rank = d$rank, coefficient = d$coefficient, contribution_years = d$contribution_years,
    contribution_pct = d$contribution_pct))
}

# Object contracts
# all_data: data.frame, first column label (chronological age), remaining columns
#           numeric protein measurements, unique participant IDs as row names.
# covariates: participant IDs (eid or row names) and the original model covariates,
#             after the study's upstream median imputation.
# blood: eid and blood_date. outcome_dates: list(hosp_events, death_data), with
#        hospital fields eid, icd3, event_date and recognized death-date fields.
# external_cohorts: list(imm_of_aging = list(data, meta),
#                       chia_als = list(data, meta)); sample IDs in row names.
#        Each data object has label first; metadata retains original diagnosis,
#        sex, tube and participant identifiers.
# bundle / result: optional previously fitted objects from these same functions;
#        pass NULL to construct the pathway sets and fit the main clock.
# reactome_cache_dir: a user-selected cache for public Reactome v97 inputs,
#        required only when bundle is NULL. No analysis results are written.
# disease_catalog: the original frozen ICD-3 endpoint metadata, if supplied.
run_imm_pathway_analysis <- function(all_data, covariates, blood, outcome_dates, external_cohorts, bundle =
  NULL, result = NULL, reactome_cache_dir = NULL, disease_catalog = NULL, baseline_age = NULL, config = imm_pathway_config()) {
  load_imm_pathway_packages()
  learners <- imm_pathway_learners()
  if (is.null(bundle)) {
    if (is.null(reactome_cache_dir)) stop("Supply bundle or a Reactome reference cache directory.")
    bundle <- build_reactome_immune_sets(all_data, cache_dir = reactome_cache_dir, parameters = reactome_immune_parameters(min_depth =
      2L, min_measured = 10L, max_measured = 100L, min_coverage = 0.10, hierarchy_jaccard = 0.80, min_closest_descendants =
      2L, closest_descendant_union_coverage = 0.80, nonhierarchical_jaccard = 0.80))
  }
  if (is.null(result)) result <- run_pathway_clock(all_data, bundle$sets, L1 = learners$stage1, L2 = learners$stage2,
    FUS = learners$fusion, S = config$seed, PTRAIN = config$training_fraction, KOUT = config$outer_folds,
    KIN = config$inner_folds, NC = config$model_cores, TLIM = config$fit_timeout, KEEP = config$retain_models)
  performance <- analyze_clock_performance(result, bundle, bootstrap_B = config$bootstrap, seed = config$seed,
    expected_n = config$expected_test_n)
  metadata <- imm_pathway_metadata(result, bundle, performance)
  overview <- ia_reactome8_overview(result, bundle, metadata$annotation)
  clinical <- run_clinical_analysis(result, covariates, blood, outcome_dates, baseline_age, config = list(cores =
    config$clinical_cores, minimum_events = config$minimum_events, alpha = config$alpha, disease_censor =
    config$disease_censor, death_censor = config$death_censor))
  disease_summaries <- list(prevalent = summarize_prevalent_diseases(clinical, bundle, disease_n = 20L), incident =
    summarize_incident_diseases(clinical, bundle$sets, bundle))
  mortality <- list(spline_risk = run_mortality_ph_extensions(clinical$prepared$BASE, as.Date(config$death_censor),
    clinical$results$mort_int, cores = config$clinical_cores, boot = config$ph_risk_bootstrap, seed = config$ph_risk_seed),
    disease_profiles = run_disease_mortality_profiles(clinical$results$inc_path, clinical$results$mort_path,
    select_mortality_profile_diseases(clinical$results$inc_int)))
  if (mortality$spline_risk$manifest$horizon != 10) stop("The manuscript requires a supported ten-year mortality risk horizon.")
  if (length(mortality$disease_profiles$common_pathways) != 58L) stop("The manuscript profile comparison requires all 58 matched pathways.")
  external <- run_external_validation(all_data, bundle$sets, L1 = learners$stage1, L2 = learners$stage2, annotation =
    metadata$external_annotation, ukb_contribution = metadata$external_contribution, external_cohorts = external_cohorts,
    selected_pipeline = list(stage1 = result$best$stage1, stage2 = result$best$stage2, pipeline = result$best$pipeline),
    S = config$seed, KIN = config$inner_folds, NC = config$external_cores, TLIM = config$fit_timeout, FUS =
    learners$fusion, POST_BOOT = config$bootstrap)

  # Frozen endpoint definitions preserve the original comparison families.
  if (is.null(disease_catalog)) disease_catalog <- sf_frozen_catalog()
  comparator <- run_all_protein_comparator(result, all_data, covariates, blood, outcome_dates, disease_catalog =
    disease_catalog, config = config$comparator)
  disease_summaries$comparator_prevalent <- summarize_disease_comparator(comparator$prevalent$associations,
    cfg = config$comparator, type = "prevalent")
  disease_summaries$comparator_incident <- summarize_disease_comparator(comparator$incident$associations,
    cfg = config$comparator, type = "incident")
  list(config = config, pathways = bundle, annotation = metadata$annotation, pathway_overview = overview,
    clock = result, performance = performance, clinical = clinical, disease_summaries = disease_summaries,
    mortality = mortality, external = external, all_protein_comparator = comparator, session = sessionInfo())
}

# Example: supply the study objects already loaded in the R session.
# analysis <- run_imm_pathway_analysis(
#   all_data = all_data, covariates = covariates, blood = blood,
#   outcome_dates = outcome_dates, external_cohorts = external_cohorts,
#   bundle = bundle
# )
# Results remain in analysis; plotting and publication-file export are excluded.


# ==============================================================================
# 01 | Immune pathway definitions
# ==============================================================================

# Reactome pathway definition and biological annotation
# ----------------------------------------------------
# Eligibility and pruning use numeric protein feature names, not protein values,
# chronological age, or clinical outcomes. Reactome v97 files are cached locally
# and verified against the checksums from the supplied analysis.
# The active study configuration uses coverage >= 0.10 and nonhierarchical
# Jaccard >= 0.80; the original parameter-constructor defaults remain available.
# Results are returned in memory: sets, audit, membership, overlap pairs,
# source checksums, and parameters. No analysis tables are exported here.
.reactome_immune_raw_defaults <- function() {
  list(min_depth = 2L, min_measured = 10L, max_measured = 100L, min_coverage = 0.20, hierarchy_jaccard =
    0.80, min_closest_descendants = 2L, closest_descendant_union_coverage = 0.80, nonhierarchical_jaccard =
    0.70)
} 

.reactome_immune_scalar <- function(x, name, lower, upper, integer = FALSE, lower_open = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x)) stop(name, " must be one finite numeric value")
  lower_ok <- if (lower_open) x > lower else x >= lower
  if (!lower_ok || x > upper) stop(name, " must be ", if (lower_open) "> " else ">= ", lower, " and <= ",
    upper)
  if (integer && x != floor(x)) stop(name, " must be a whole number")
  if (integer) as.integer(x) else as.numeric(x)
} 

.normalize_reactome_immune_parameters <- function(parameters = list()) {
  defaults <- .reactome_immune_raw_defaults()
  if (is.null(parameters)) parameters <- list()
  if (!is.list(parameters)) stop("parameters must be a named list")
  if (length(parameters)) {
    nms <- names(parameters)
    if (is.null(nms) || any(is.na(nms) | !nzchar(nms)) || anyDuplicated(nms)) stop("parameters must have unique, non-empty names")
    unknown <- setdiff(nms, names(defaults))
    if (length(unknown)) stop("Unknown parameter(s): ", paste(unknown, collapse = ", "))
    defaults[nms] <- parameters
  }
  p <- defaults
  p$min_depth <- .reactome_immune_scalar(p$min_depth, "min_depth", 0, .Machine$integer.max, integer = TRUE)
  p$min_measured <- .reactome_immune_scalar(p$min_measured, "min_measured", 1, .Machine$integer.max, integer =
    TRUE)
  p$max_measured <- .reactome_immune_scalar(p$max_measured, "max_measured", 1, .Machine$integer.max, integer =
    TRUE)
  if (p$max_measured < p$min_measured) stop("max_measured must be greater than or equal to min_measured")
  p$min_coverage <- .reactome_immune_scalar(p$min_coverage, "min_coverage", 0, 1)
  p$hierarchy_jaccard <- .reactome_immune_scalar(p$hierarchy_jaccard, "hierarchy_jaccard", 0, 1, lower_open =
    TRUE)
  p$min_closest_descendants <- .reactome_immune_scalar(p$min_closest_descendants, "min_closest_descendants",
    1, .Machine$integer.max, integer = TRUE)
  p$closest_descendant_union_coverage <- .reactome_immune_scalar(p$closest_descendant_union_coverage, "closest_descendant_union_coverage",
    0, 1, lower_open = TRUE)
  p$nonhierarchical_jaccard <- .reactome_immune_scalar(p$nonhierarchical_jaccard, "nonhierarchical_jaccard",
    0, 1, lower_open = TRUE)
  p
}
reactome_immune_parameters <- function(min_depth = 2L, min_measured = 10L, max_measured = 100L, min_coverage =
  0.20, hierarchy_jaccard = 0.80, min_closest_descendants = 2L, closest_descendant_union_coverage = 0.80,
  nonhierarchical_jaccard = 0.70) {
  .normalize_reactome_immune_parameters(list(min_depth = min_depth, min_measured = min_measured, max_measured =
    max_measured, min_coverage = min_coverage, hierarchy_jaccard = hierarchy_jaccard, min_closest_descendants =
    min_closest_descendants, closest_descendant_union_coverage = closest_descendant_union_coverage, nonhierarchical_jaccard =
    nonhierarchical_jaccard))
} 

.reactome_immune_parameter_metadata <- function(parameters) {
  p <- .normalize_reactome_immune_parameters(parameters)
  fmt <- function(x) format(x, trim = TRUE, scientific = FALSE, digits = 15)
  tag <- function(x) {
    z <- fmt(x)
    z <- gsub("-", "m", z, fixed = TRUE)
    gsub(".", "p", z, fixed = TRUE)
  }
  parameter_set_name <- paste0("d", tag(p$min_depth), "__m", tag(p$min_measured), "_", tag(p$max_measured),
    "__cov", tag(p$min_coverage), "__hj", tag(p$hierarchy_jaccard), "__cd", tag(p$min_closest_descendants),
    "__uc", tag(p$closest_descendant_union_coverage), "__nj", tag(p$nonhierarchical_jaccard))
  list(parameter_set_name = parameter_set_name, eligibility = paste0("minimum directed depth >=", fmt(p$min_depth),
    "; non-disease; ", fmt(p$min_measured), "-", fmt(p$max_measured), " measured Olink proteins; measured/Reactome coverage >=",
    fmt(100 * p$min_coverage), "%"), hierarchy_pruning = paste0("any retained descendant has ancestor-descendant Jaccard >=",
    fmt(100 * p$hierarchy_jaccard), "%, or >=", fmt(p$min_closest_descendants), " closest retained descendants whose ",
    "union covers >=", fmt(100 * p$closest_descendant_union_coverage), "% of the parent's measured proteins"),
    nonhierarchical_pruning = paste0("non-hierarchical measured-protein Jaccard >=", fmt(100 * p$nonhierarchical_jaccard),
    "%"))
}
build_reactome_immune_sets <- function(all_data, label = "label", cache_dir = "<REACTOME_V97_CACHE_DIRECTORY>",
  release = 97L, root = "R-HSA-168256", parameters = reactome_immune_parameters(min_coverage = 0.10, nonhierarchical_jaccard =
  0.80)) {
  if (length(release) != 1L || !(is.integer(release) || is.double(release)) || is.na(release) || !is.finite(release) ||
    release != 97) stop("This reproducible script is fixed to Reactome v97")
  release <- 97L
  if (length(root) != 1L || !is.character(root) || is.na(root) || root != "R-HSA-168256") stop("This script is fixed to Reactome Immune System R-HSA-168256")
  root <- "R-HSA-168256"
  standard <- .normalize_reactome_immune_parameters(parameters)
  base_url <- sprintf("https://download.reactome.org/%d/", release)
  src <- c("ReactomePathways.gmt.zip", "ReactomePathwaysRelation.txt", "ReactomePathways.txt", "HumanDiseasePathways.txt")
  expected_md5 <- c("ReactomePathways.gmt.zip" = "f660bae1d7a48a8b5e92df1409cb2eff", "ReactomePathwaysRelation.txt" =
    "4cd76055c4ea7d18308d5993f160533d", "ReactomePathways.txt" = "63953262f948a43232e01ac18c9606d2", "HumanDiseasePathways.txt" =
    "0d0255a7b74789141ff78f3987c82cbe")
  source_dir <- cache_dir
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  fp <- setNames(file.path(source_dir, src), src)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(600, old_timeout))
  for (f in src) {
    needs_download <- !file.exists(fp[[f]]) || is.na(file.info(fp[[f]])$size) || file.info(fp[[f]])$size ==
    0
    if (!needs_download) {
      cached_md5 <- unname(tools::md5sum(fp[[f]]))
      needs_download <- is.na(cached_md5) || cached_md5 != expected_md5[[f]]
    }
    if (needs_download) {
      z <- try(utils::download.file(paste0(base_url, f), fp[[f]], mode = "wb", method = "libcurl"), silent =
        TRUE)
      bad <- inherits(z, "try-error") || length(z) != 1L || as.integer(z) != 0L || !file.exists(fp[[f]]) ||
      is.na(file.info(fp[[f]])$size) || file.info(fp[[f]])$size == 0
      if (bad) {
        if (file.exists(fp[[f]])) unlink(fp[[f]])
        stop("Download failed: ", f)
      }
    }
  }
  actual_md5 <- unname(tools::md5sum(fp))
  if (any(is.na(actual_md5) | actual_md5 != unname(expected_md5[src]))) stop("Reactome v97 source checksum mismatch; delete the affected cache file and rerun")
  md5 <- data.frame(release, file = src, md5 = actual_md5, expected_md5 = unname(expected_md5[src]), row.names =
    NULL)
  rd <- function(file, nms) {
    x <- utils::read.delim(file, header = FALSE, sep = "\t", quote = "", comment.char = "", colClasses =
      "character", check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, fileEncoding = "UTF-8")
    if (ncol(x) != length(nms)) stop("Unexpected format: ", basename(file))
    names(x) <- nms;
    x
  }
  paths <- rd(fp[["ReactomePathways.txt"]], c("id", "pathway", "species"))
  paths <- paths[paths$species == "Homo sapiens", , drop = FALSE]
  paths <- paths[!duplicated(paths$id), , drop = FALSE]
  rel <- rd(fp[["ReactomePathwaysRelation.txt"]], c("parent", "child"))
  rel <- unique(rel[startsWith(rel$parent, "R-HSA-") & startsWith(rel$child, "R-HSA-"), , drop = FALSE])
  disease <- unique(rd(fp[["HumanDiseasePathways.txt"]], c("id", "disease_name"))$id)
  td <- tempfile("reactome_gmt_");
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(fp[["ReactomePathways.gmt.zip"]], exdir = td)
  gf <- list.files(td, "^ReactomePathways\\.gmt$", recursive = TRUE, full.names = TRUE)
  if (length(gf) != 1L) stop("ReactomePathways.gmt not found in zip")
  zz <- strsplit(readLines(gf, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  zz <- zz[lengths(zz) >= 3L]
  gid <- vapply(zz, function(x) x[[2L]], character(1))
  if (!all(startsWith(gid, "R-HSA-"))) stop("Unexpected GMT stable-ID column")
  ix <- split(seq_along(zz), gid)
  parse_genes <- function(x) {
    z <- trimws(x[ - (1:2)])
    if (length(z) && identical(tolower(z[[1L]]), "reactome pathway")) z <- z[ - 1L]
    unique(toupper(z[nzchar(z)]))
  }
  gmt <- lapply(ix, function(ii) unique(unlist(lapply(zz[ii], parse_genes), use.names = FALSE)))
  gmt_title <- vapply(ix, function(ii) zz[[ii[[1L]]]][[1L]], character(1))
  gmt <- gmt[startsWith(names(gmt), "R-HSA-")]

  ## Directed BFS: minimum depth below Reactome Immune System.
  kids <- split(rel$child, rel$parent)
  depth <- setNames(0L, root);
  frontier <- root;
  d <- 0L
  repeat {
    nx <- setdiff(unique(unlist(kids[frontier], use.names = FALSE)), names(depth))
    if (!length(nx)) break
    d <- d + 1L;
    depth[nx] <- d;
    frontier <- nx
  }
  ids <- names(depth)
  if (length(ids) == 1L) stop("Immune System root has no descendants")
  rel <- rel[rel$parent %in% ids & rel$child %in% ids, , drop = FALSE]
  kids <- split(rel$child, rel$parent);
  parents <- split(rel$parent, rel$child)
  descendants <- function(x) {
    out <- character();
    q <- kids[[x]]
    while (length(q)) {
      q <- setdiff(unique(q), out)
      if (!length(q)) break
      out <- c(out, q);
      q <- unlist(kids[q], use.names = FALSE)
    }
    setdiff(intersect(out, ids), x)
  }
  desc <- setNames(lapply(ids, descendants), ids)

  ## Only column names/types define the measured Olink universe here; neither
  ## protein values nor the age label enter pathway eligibility or pruning.
  data_names <- names(all_data)
  if (length(label) != 1L || !is.character(label) || is.na(label) || !nzchar(label)) stop("label must be one non-empty column name")
  if (is.null(data_names) || !label %in% data_names) stop("Missing outcome column: ", label)
  if (anyDuplicated(data_names)) stop("Duplicated column names in all_data")
  feat <- setdiff(data_names, label)
  feat <- feat[vapply(feat, function(x) is.numeric(all_data[[x]]), logical(1))]
  all_symbols <- unique(unlist(gmt, use.names = FALSE))
  hit <- toupper(feat) %in% all_symbols
  if (!any(hit)) stop("No Olink column matched a Reactome gene symbol. ", "Keep official gene symbols in colnames(all_data).")
  symbol <- toupper(feat[hit])
  if (anyDuplicated(symbol)) stop("Duplicated case-insensitive Olink gene symbols")
  col_of <- setNames(feat[hit], symbol)
  gs <- setNames(lapply(ids, function(x) unique(gmt[[x]])), ids)
  gs <- lapply(gs, function(x) x[!is.na(x) & nzchar(x)])
  measured_symbol <- lapply(gs, function(x) intersect(x, names(col_of)))
  measured <- lapply(measured_symbol, function(x) unname(col_of[x]))
  ttl <- paths$pathway[match(ids, paths$id)]
  miss <- is.na(ttl) | !nzchar(ttl);
  ttl[miss] <- gmt_title[ids[miss]]
  miss <- is.na(ttl) | !nzchar(ttl);
  ttl[miss] <- ids[miss]
  top <- intersect(kids[[root]], ids)
  top_name <- setNames(ttl[match(top, ids)], top)
  branch <- vapply(ids, function(x) {
    b <- top[vapply(top, function(y) identical(x, y) || x %in% desc[[y]], logical(1))]
    paste(unname(top_name[b]), collapse = "; ")
  }, character(1))
  audit <- data.frame(release, id = ids, pathway = ttl, min_depth = unname(depth[ids]), branch, parent_ids =
    vapply(ids, function(x) paste(parents[[x]], collapse = ";"), character(1)), disease = ids %in% disease,
    n_reactome = lengths(gs), n_olink = lengths(measured), stringsAsFactors = FALSE)
  audit$coverage <- with(audit, ifelse(n_reactome > 0, n_olink / n_reactome, 0))
  eligible <- with(audit, min_depth >= standard$min_depth & !disease & n_olink >= standard$min_measured &
    n_olink <= standard$max_measured & coverage >= standard$min_coverage)
  rule <- rep("excluded_below_min_coverage", nrow(audit))
  rule[audit$n_olink < standard$min_measured] <- "excluded_below_min_measured"
  rule[audit$n_olink > standard$max_measured] <- "excluded_above_max_measured"
  rule[audit$n_reactome == 0L] <- "excluded_no_gene_symbols"
  rule[audit$min_depth < standard$min_depth] <- "excluded_below_min_depth"
  rule[audit$disease] <- "excluded_disease"
  rule[eligible] <- "eligible_active_parameters"
  audit$eligibility_rule <- rule
  audit$eligible <- eligible
  cand <- audit$id[audit$eligible]
  if (!length(cand)) stop("No eligible immune pathway; inspect feature symbols/thresholds")

  ## Reachability and similarity are computed only on measured Olink proteins.
  nc <- length(cand)
  ancestor <- matrix(FALSE, nc, nc, dimnames = list(cand, cand))
  for (x in cand) ancestor[x, intersect(desc[[x]], cand)] <- TRUE
  J <- matrix(0, nc, nc, dimnames = list(cand, cand))
  I <- matrix(0L, nc, nc, dimnames = list(cand, cand))
  diag(J) <- 1;
  diag(I) <- lengths(measured[cand])
  if (nc > 1L) for (i in seq_len(nc - 1L)) for (j in (i + 1L):nc) {
    a <- measured[[cand[[i]]]];
    b <- measured[[cand[[j]]]]
    ni <- length(intersect(a, b))
    I[i, j] <- I[j, i] <- ni
    J[i, j] <- J[j, i] <- ni / length(union(a, b))
  }
  hierarchy <- ancestor | t(ancestor)
  red <- J >= standard$nonhierarchical_jaccard & !hierarchy
  diag(red) <- FALSE

  ## Among the already retained descendants of p, the closest frontier contains
  ## nodes with no other retained pathway between themselves and p in the DAG.
  closest_retained_descendants <- function(p, retained) {
    ds <- intersect(retained, desc[[p]])
    if (!length(ds)) return(character())
    ds[colSums(ancestor[ds, ds, drop = FALSE]) == 0L]
  }

  ## Hierarchy pruning: remove p if any retained descendant reaches the active
  ## hierarchy-Jaccard threshold, or if the required number of closest retained
  ## descendants jointly reaches the active union-coverage threshold.
  hierarchy_evidence <- function(p, retained) {
    ds <- intersect(retained, desc[[p]])
    near <- ds[J[p, ds] >= standard$hierarchy_jaccard]
    if (length(near)) {
      ii <- match(near, audit$id)
      q <- near[order( - J[p, near], - audit$coverage[ii], near)][[1L]]
      return(list(valid = TRUE, reason = "descendant_high_jaccard", reference = q, coverage = NA_real_))
    }
    nearest <- closest_retained_descendants(p, retained)
    if (!length(nearest)) return(list(valid = FALSE, reason = NA_character_, reference = character(), coverage =
      NA_real_))
    cov <- length(intersect(measured[[p]], unique(unlist(measured[nearest], use.names = FALSE)))) / length(measured[[p]])
    if (length(nearest) >= standard$min_closest_descendants && cov >= standard$closest_descendant_union_coverage) {
      return(list(valid = TRUE, reason = "closest_descendants_cover_parent", reference = nearest, coverage =
        cov))
    }
    list(valid = FALSE, reason = NA_character_, reference = nearest, coverage = cov)
  }

  ## A reverse topological order guarantees that all descendants have a final
  ## keep/remove decision before an ancestor is examined. Ties are deterministic.
  reverse_topological_order <- function(nodes) {
    remaining <- nodes;
    out <- character()
    while (length(remaining)) {
      sinks <- remaining[rowSums(ancestor[remaining, remaining, drop = FALSE]) == 0L]
      if (!length(sinks)) stop("Cycle detected in candidate Reactome hierarchy")
      ii <- match(sinks, audit$id)
      sinks <- sinks[order( - audit$min_depth[ii], - audit$coverage[ii], sinks)]
      out <- c(out, sinks)
      remaining <- setdiff(remaining, sinks)
    }
    out
  }

  ## Single-pass bottom-up greedy selection. Hierarchy redundancy takes priority;
  ## otherwise a node is kept only if it has no non-hierarchical Jaccard conflict
  ## with an already retained node. Retained nodes are never removed later.
  keep <- character()
  h_remove <- h_reason <- h_rep <- setNames(character(), character())
  h_cover <- setNames(numeric(), character())
  o_remove <- o_rep <- setNames(character(), character())
  for (x in reverse_topological_order(cand)) {
    ev <- hierarchy_evidence(x, keep)
    if (ev$valid) {
      h_remove[x] <- x
      h_reason[x] <- ev$reason
      h_rep[x] <- paste(ev$reference, collapse = ";")
      h_cover[x] <- ev$coverage
      next
    }
    conflict <- keep[red[x, keep]]
    if (length(conflict)) {
      ii <- match(conflict, audit$id)
      q <- conflict[order( - J[x, conflict], - audit$coverage[ii], conflict)][[1L]]
      o_remove[x] <- x;
      o_rep[x] <- q
    } else {
      keep <- c(keep, x)
    }
  }
  h_active <- names(h_remove)

  ## Shared-protein metrics are descriptive audit fields only. They never remove
  ## a pathway and are not model penalty factors or observation weights.
  shared_audit <- function(z) {
    if (!length(z)) return(list(k = table(character()), adjusted = setNames(numeric(), character()), unique =
      setNames(numeric(), character())))
    k <- table(unlist(measured[z], use.names = FALSE))
    ad <- vapply(z, function(x) sum(1 / as.numeric(k[measured[[x]]])), numeric(1))
    uq <- vapply(z, function(x) sum(k[measured[[x]]] == 1L), numeric(1))
    list(k = k, adjusted = ad, unique = uq)
  }
  final <- shared_audit(keep)

  ## Internal invariants for the active parameter set.
  if (!all(keep %in% cand) || anyDuplicated(keep) > 0L) stop("Internal error: retained pathways are not a unique candidate subset")
  if (length(keep) > 1L && any(red[keep, keep, drop = FALSE])) stop("Internal error: retained non-hierarchical Jaccard conflict")
  if (length(keep) && any(vapply(keep, function(x) hierarchy_evidence(x, keep)$valid, logical(1)))) stop("Internal error: retained hierarchical redundancy")
  if (length(h_active) && any(!vapply(h_active, function(x) hierarchy_evidence(x, keep)$valid, logical(1)))) stop("Internal error: unsupported hierarchy removal")
  if (length(o_remove) && any(!unname(o_rep) %in% keep)) stop("Internal error: non-hierarchical representative was not retained")
  accounted <- c(keep, h_active, names(o_remove))
  if (!setequal(accounted, cand) || anyDuplicated(accounted) > 0L) stop("Internal error: candidate pathways were not partitioned exactly once")
  status <- ifelse(audit$eligible, "candidate", "excluded")
  reason <- audit$eligibility_rule
  representative_ids <- rep(NA_character_, nrow(audit))
  status[audit$id %in% h_active] <- "removed_hierarchy"
  reason[match(h_active, audit$id)] <- unname(h_reason[h_active])
  representative_ids[match(h_active, audit$id)] <- unname(h_rep[h_active])
  status[audit$id %in% names(o_remove)] <- "removed_overlap"
  reason[match(names(o_remove), audit$id)] <- "nonhierarchical_overlap"
  representative_ids[match(names(o_remove), audit$id)] <- unname(o_rep[names(o_remove)])
  status[audit$id %in% keep] <- "retained";
  reason[audit$id %in% keep] <- "retained"
  audit$status <- status;
  audit$reason <- reason
  audit$representative_ids <- representative_ids
  active_cover <- h_cover[h_active]
  audit$observed_closest_descendant_union_coverage <- active_cover[audit$id]
  audit$final_unique_protein_count <- audit$final_shared_protein_fraction <- audit$shared_membership_adjusted_count <-
  NA_real_
  if (length(keep)) {
    ii <- match(keep, audit$id)
    audit$final_unique_protein_count[ii] <- final$unique[keep]
    audit$final_shared_protein_fraction[ii] <- 1 - final$unique[keep] / audit$n_olink[ii]
    audit$shared_membership_adjusted_count[ii] <- final$adjusted[keep]
  }
  audit$max_jaccard <- NA_real_
  audit$max_jaccard_with <- NA_character_
  if (nc > 1L) for (x in cand) {
    other <- setdiff(cand, x);
    y <- other[which.max(J[x, other])]
    audit$max_jaccard[audit$id == x] <- J[x, y]
    audit$max_jaccard_with[audit$id == x] <- y
  }
  slug <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
  set_name <- paste0(keep, "__", slug(audit$pathway[match(keep, audit$id)]))
  reactome_immune_sets <- setNames(unname(measured[keep]), set_name)
  audit$set_name <- NA_character_;
  audit$set_name[match(keep, audit$id)] <- set_name
  candidate_k <- table(unlist(measured[cand], use.names = FALSE))
  membership <- if (length(keep)) do.call(rbind, lapply(keep, function(x) data.frame(id = x, pathway = audit$pathway[match(x,
    audit$id)], feature = measured[[x]], stringsAsFactors = FALSE))) else data.frame(id = character(), pathway =
    character(), feature = character(), stringsAsFactors = FALSE)
  rownames(membership) <- NULL
  membership$symbol <- names(col_of)[match(membership$feature, col_of)]
  final_k <- table(membership$feature)
  membership$k_candidate <- as.integer(candidate_k[membership$feature])
  membership$k_final <- as.integer(final_k[membership$feature])
  membership$inverse_membership_audit <- 1 / membership$k_final
  if (nc > 1L) {
    ij <- utils::combn(cand, 2L)
    overlap_pairs <- data.frame(id1 = ij[1L, ], id2 = ij[2L, ], intersection = I[cbind(match(ij[1L, ], cand),
      match(ij[2L, ], cand))], jaccard = J[cbind(match(ij[1L, ], cand), match(ij[2L, ], cand))], hierarchical =
      hierarchy[cbind(match(ij[1L, ], cand), match(ij[2L, ], cand))], stringsAsFactors = FALSE)
    overlap_pairs$redundant_nonhierarchical <- with(overlap_pairs, !hierarchical & jaccard >= standard$nonhierarchical_jaccard)
  } else overlap_pairs <- data.frame(id1 = character(), id2 = character(), intersection = integer(), jaccard =
    numeric(), hierarchical = logical(), redundant_nonhierarchical = logical())
  parameter_metadata <- .reactome_immune_parameter_metadata(standard)
  params <- list(release = release, root = root, standard_name = parameter_metadata$parameter_set_name, parameter_set_name =
    parameter_metadata$parameter_set_name, standard = standard, input_audit = list(n_data_columns = length(data_names),
    n_numeric_nonlabel_columns = length(feat), n_reactome_symbol_matched_columns = length(col_of)), eligibility =
    parameter_metadata$eligibility, closest_descendant = paste("retained descendant with no other retained pathway between it and",
    "the ancestor in the Reactome DAG"), hierarchy_pruning = parameter_metadata$hierarchy_pruning, nonhierarchical_pruning =
    parameter_metadata$nonhierarchical_pruning, selection_algorithm = paste("single-pass reverse-topological greedy: descendants before ancestors;",
    "ties by greater minimum depth, greater coverage and stable ID; test", "hierarchy redundancy first, then non-hierarchical conflict"),
    shared_proteins = paste("retained; k, 1/k and shared-membership-adjusted counts are audit-only", "and do not enter pathway selection or model weighting"),
    outcome_blinding = paste("pathway selection uses feature names only; age and downstream", "phenotype, disease and mortality outcomes are not used"))
  message(sprintf("Reactome v%d: %d immune nodes -> %d eligible -> %d retained; %d proteins", release, length(ids),
    length(cand), length(keep), length(final$k)))
  invisible(list(sets = reactome_immune_sets, audit = audit, membership = membership, overlap_pairs = overlap_pairs,
    source_md5 = md5, params = params))
}

# Study-defined eight-class annotation from stable-ID ancestry in the same v97
# hierarchy. A pathway can have multiple memberships. Its primary class uses
# the nearest anchor, with C1-C8 order breaking ties.
reactome8_classes <- function() data.frame(class_id = paste0("C", 1:8), class_name = c("Pattern recognition and innate immune sensing",
  "Complement system", "Innate immune effector functions", "Antigen processing and presentation", "Lymphocyte signaling and regulation",
  "Interleukin signaling", "Interferon signaling", "Other cytokine and related receptor signaling"), short_label =
  c("Innate sensing", "Complement", "Innate effectors", "Antigen presentation", "Lymphocyte regulation", "Interleukins",
  "Interferons", "Other cytokines"), stringsAsFactors = FALSE)
reactome8_anchors <- function() {
  a <- utils::read.delim(text = "class_id\tid\tpathway
C1\tR-HSA-168898\tToll-like Receptor Cascades
C1\tR-HSA-168643\tNucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways
C1\tR-HSA-879415\tAdvanced glycosylation endproduct receptor signaling
C1\tR-HSA-168928\tDDX58/IFIH1-mediated induction of interferon-alpha/beta
C1\tR-HSA-1834949\tCytosolic sensors of pathogen-associated DNA
C1\tR-HSA-5621481\tC-type lectin receptors (CLRs)
C1\tR-HSA-9645460\tAlpha-protein kinase 1 signaling pathway
C1\tR-HSA-9686114\tNon-canonical inflammasome activation
C2\tR-HSA-166658\tComplement cascade
C3\tR-HSA-2029480\tFcgamma receptor (FCGR) dependent phagocytosis
C3\tR-HSA-2172127\tDAP12 interactions
C3\tR-HSA-2454202\tFc epsilon receptor (FCERI) signaling
C3\tR-HSA-6803157\tAntimicrobial peptides
C3\tR-HSA-6798695\tNeutrophil degranulation
C3\tR-HSA-1222556\tROS and RNS production in phagocytes
C3\tR-HSA-9970672\tFXIIa activates plasma kallikrein-kinin system
C4\tR-HSA-983169\tClass I MHC mediated antigen processing & presentation
C4\tR-HSA-2132295\tMHC class II antigen presentation
C5\tR-HSA-202403\tTCR signaling
C5\tR-HSA-388841\tRegulation of T cell activation by CD28 family
C5\tR-HSA-983705\tSignaling by the B Cell Receptor (BCR)
C5\tR-HSA-198933\tImmunoregulatory interactions between a Lymphoid and a non-Lymphoid cell
C5\tR-HSA-392517\tRap1 signalling
C5\tR-HSA-8851680\tButyrophilin (BTN) family interactions
C5\tR-HSA-9938026\tImmunoglobulin maturation
C6\tR-HSA-449147\tSignaling by Interleukins
C7\tR-HSA-913531\tInterferon Signaling
C8\tR-HSA-982772\tGrowth hormone receptor signaling
C8\tR-HSA-1170546\tProlactin receptor signaling
C8\tR-HSA-5668541\tTNFR2 non-canonical NF-kB pathway
C8\tR-HSA-9607240\tFLT3 Signaling
C8\tR-HSA-9674555\tSignaling by CSF3 (G-CSF)
C8\tR-HSA-9680350\tSignaling by CSF1 (M-CSF) in myeloid cells",
    quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  a$url <- paste0("https://reactome.org/content/detail/", a$id);
  a
}
classify_reactome8 <- function(bundle, modules, anchors = reactome8_anchors()) {
  if (!is.list(bundle) || is.null(bundle$audit)) stop("Supply the complete Reactome v97 pathway bundle, including its full audit.",
    call. = FALSE)
  a <- as.data.frame(bundle$audit, stringsAsFactors = FALSE);
  req <- c("release", "id", "pathway", "parent_ids")
  if (!all(req %in% names(a))) stop("bundle$audit lacks: ", paste(setdiff(req, names(a)), collapse = ", "),
    call. = FALSE)
  if (anyNA(a$id) || anyDuplicated(a$id) || any(!grepl("^R-HSA-[0-9]+$", a$id))) stop("Invalid or duplicate stable IDs in audit.",
    call. = FALSE)
  if (any(is.na(a$release) | as.character(a$release) != "97") || (!is.null(bundle$params$release) && !identical(as.character(bundle$params$release),
    "97"))) stop("This annotation requires consistent Reactome v97 release metadata.", call. = FALSE)
  cls <- reactome8_classes();
  cc <- cls$class_id;
  modules <- as.character(modules);
  ids <- sub("__.*$", "", modules)
  if (!length(modules) || anyNA(modules) || anyDuplicated(ids) || any(!grepl("^R-HSA-[0-9]+$", ids))) stop("modules must contain unique stable IDs or original R-HSA-...__... module names.",
    call. = FALSE)
  if (!all(ids %in% a$id)) stop("Modules missing from the supplied audit: ", paste(setdiff(ids, a$id), collapse =
    ", "), call. = FALSE)
  if (!"R-HSA-168256" %in% a$id) stop("The audit is truncated: Immune System root R-HSA-168256 is missing.",
    call. = FALSE)
  if (!all(c("class_id", "id", "pathway") %in% names(anchors)) || anyNA(anchors[, c("class_id", "id")]) ||
    any(!anchors$class_id %in% cc) || any(!grepl("^R-HSA-[0-9]+$", anchors$id)) || anyDuplicated(anchors[,
    c("class_id", "id")])) stop("Invalid anchor table.", call. = FALSE)
  parents <- strsplit(ifelse(is.na(a$parent_ids), "", as.character(a$parent_ids)), ";", fixed = TRUE)
  parents <- lapply(parents, function(x) unique(trimws(x[nzchar(trimws(x))])));
  names(parents) <- a$id
  outside <- setdiff(unique(unlist(parents, use.names = FALSE)), a$id)
  if (length(outside)) stop("The audit is incomplete: parent nodes are missing: ", paste(outside, collapse =
    ", "), call. = FALSE)
  edge <- data.frame(parent = unlist(parents, use.names = FALSE), child = rep(a$id, lengths(parents)), stringsAsFactors =
    FALSE)
  kids <- split(edge$child, edge$parent);
  indeg <- stats::setNames(lengths(parents), a$id);
  q <- names(indeg)[indeg == 0L];
  visited <- character()
  while (length(q)) {
    z <- q[1L];
    q <- q[ - 1L];
    visited <- c(visited, z);
    nx <- kids[[z]];
    if (length(nx)) {
      indeg[nx] <- indeg[nx] - 1L;
      q <- c(q, nx[indeg[nx] == 0L])
    }
  }
  if (length(visited) != nrow(a)) stop("The audit contains a hierarchy cycle; annotation stopped.", call. =
    FALSE)
  bfs <- function(root) {
    d <- stats::setNames(rep(Inf, nrow(a)), a$id);
    if (!root %in% a$id) return(d);
    d[root] <- 0;
    q <- root;
    k <- 0L
    while (length(q)) {
      nx <- setdiff(unique(unlist(kids[q], use.names = FALSE)), names(d)[is.finite(d)]);
      if (!length(nx)) break;
      k <- k + 1L;
      d[nx] <- k;
      q <- nx
    };
    d
  }
  reachable <- bfs("R-HSA-168256")
  if (any(!is.finite(reachable))) stop("The audit contains nodes disconnected from the Immune System root.",
    call. = FALSE)
  da <- vapply(anchors$id, bfs, stats::setNames(numeric(nrow(a)), a$id));
  colnames(da) <- paste(anchors$class_id, anchors$id, sep = ":")
  dc <- vapply(cc, function(cl) {
    j <- which(anchors$class_id == cl);
    if (length(j)) apply(da[, j, drop = FALSE], 1, min) else rep(Inf, nrow(a))
  }, numeric(nrow(a)))
  rownames(dc) <- a$id;
  dc <- dc[ids, , drop = FALSE];
  mem <- is.finite(dc);
  bad <- rowSums(mem) == 0L
  if (any(bad)) {
    diag <- data.frame(module = modules[bad], a[match(ids[bad], a$id), intersect(c("id", "pathway", "parent_ids",
      "branch", "status"), names(a)), drop = FALSE], row.names = NULL, check.names = FALSE)
    msg <- paste0("No verified eight-class ancestry for ", sum(bad), " module(s). Do not assign these to C8 by default. Check the complete v97 audit or add reviewed stable-ID anchors:\n",
      paste(paste(diag$id, diag$pathway, sep = " | "), collapse = "\n"))
    stop(structure(list(message = msg, call = NULL, diagnostics = diag, anchors = anchors), class = c("reactome8_unclassified",
      "error", "condition")))
  }
  primary <- cc[max.col( - dc, ties.method = "first")];
  mind <- apply(dc, 1, min);
  hit <- da[match(ids, a$id), , drop = FALSE]
  out <- data.frame(module = modules, id = ids, pathway = as.character(a$pathway[match(ids, a$id)]), primary_class =
    primary, classes = apply(mem, 1, function(z) paste(cc[z], collapse = ";")), n_classes = rowSums(mem),
    primary_tie = rowSums(dc == mind) > 1L, matched_anchors = apply(is.finite(hit), 1, function(z) paste(unique(anchors$id[z]),
    collapse = ";")), mem, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  anchors$present_in_audit <- anchors$id %in% a$id;
  anchors$audit_pathway <- as.character(a$pathway[match(anchors$id, a$id)])
  attr(out, "classes") <- cls;
  attr(out, "anchors") <- anchors;
  attr(out, "class_distances") <- dc
  attr(out, "classification_note") <- "Study-defined eight-class grouping of Reactome stable-ID anchor descendants, using only the supplied v97 DAG; multiple membership retained. Primary classification uses the nearest anchor, with C1-C8 order breaking ties."
  out
}

# Validate loaded analysis objects and extract Reactome stable identifiers.
ia_read_saved <- function(x, label) {
  if (!is.list(x)) stop(label, " must be a loaded analysis list.", call. = FALSE)
  x
}
ia_reactome_id <- function(x) {
  x <- as.character(x)
  vapply(x, function(z) {
    if (is.na(z) || !nzchar(trimws(z)) || z %in% c("NA", "N/A")) return(NA_character_)
    hits <- regmatches(z, gregexpr("R-HSA-[0-9]+", z))[[1L]]
    if (length(unique(hits)) > 1L) stop("Multiple stable pathway IDs in one key: ", z, call. = FALSE)
    if (!length(hits)) NA_character_ else hits[[1L]]
  }, character(1), USE.NAMES = FALSE)
}
ia_reactome8_annotation <- function(result, bundle) {
  result <- ia_read_saved(result, "result");
  bundle <- ia_read_saved(bundle, "bundle")
  modules <- as.character(result$mapping$module)
  if (length(modules) != 58L || anyNA(modules) || anyDuplicated(ia_reactome_id(modules))) stop("result$mapping$module must contain all 58 unique fitted pathways.",
    call. = FALSE)
  ann <- classify_reactome8(bundle, modules)
  cls <- reactome8_classes()
  ann$primary_class_name <- cls$class_name[match(ann$primary_class, cls$class_id)]
  ann$class_names <- vapply(strsplit(ann$classes, ";", fixed = TRUE), function(z) paste(cls$class_name[match(z,
    cls$class_id)], collapse = ";"), character(1))
  ann
}
ia_reactome8_overview <- function(result, bundle, annotation = NULL) {
  result <- ia_read_saved(result, "result");
  bundle <- ia_read_saved(bundle, "bundle")
  if (is.null(annotation)) annotation <- ia_reactome8_annotation(result, bundle)
  modules <- annotation$module;
  features <- result$module_features
  if (!is.list(features) || !all(modules %in% names(features))) stop("Overview counts require result$module_features for all 58 fitted pathways.",
    call. = FALSE)
  features <- lapply(features[modules], function(x) {
    x <- as.character(x)
    if (!length(x) || anyNA(x) || any(!nzchar(x))) stop("Empty or missing model feature identifiers.", call. =
      FALSE)
    unique(x)
  })
  edges <- unique(do.call(rbind, lapply(modules, function(m) data.frame(module = m, reactome_id = annotation$id[match(m,
    modules)], feature = features[[m]], stringsAsFactors = FALSE))))
  mapping <- as.data.frame(bundle$membership, stringsAsFactors = FALSE)
  if (!all(c("feature", "symbol") %in% names(mapping))) stop("Overview protein counts require bundle$membership feature-to-symbol mapping.",
    call. = FALSE)
  symbols <- split(as.character(mapping$symbol), as.character(mapping$feature))
  edges$protein <- vapply(edges$feature, function(f) {
    z <- unique(symbols[[f]]);
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z) != 1L) stop("Expected exactly one protein symbol for fitted feature: ", f, call. = FALSE)
    z
  }, character(1))
  cls <- reactome8_classes()
  counts <- cls
  counts$n_pathways_primary <- vapply(cls$class_id, function(c) sum(annotation$primary_class == c), integer(1))
  counts$n_pathways_membership <- vapply(cls$class_id, function(c) sum(annotation[[c]]), integer(1))
  counts$n_unique_proteins_primary <- vapply(cls$class_id, function(c) length(unique(edges$protein[edges$module %in%
    modules[annotation$primary_class == c]])), integer(1))
  counts$n_unique_proteins_membership <- vapply(cls$class_id, function(c) length(unique(edges$protein[edges$module %in%
    modules[annotation[[c]]]])), integer(1))
  counts$counting_note <- "Primary pathway counts are exclusive; full-membership pathway and within-class protein counts overlap between classes."
  totals <- data.frame(quantity = c("fitted_pathways", "unique_fitted_features", "unique_fitted_proteins",
    "protein_pathway_memberships"), value = c(nrow(annotation), length(unique(edges$feature)), length(unique(edges$protein)),
    nrow(unique(edges[c("protein", "module")]))), source = c("result$mapping$module", "result$module_features",
    "result$module_features + bundle$membership", "result$module_features + bundle$membership"), stringsAsFactors =
    FALSE)
  list(class_counts = counts, protein_membership = edges, totals = totals)
}


# ==============================================================================
# 02 | Pathway clock development
# ==============================================================================

# ============================================================================
# Pathway age model development and evaluation
# ============================================================================
# Input contract:
#   all_data: data frame with numeric chronological age in its first column,
#             named "label", numeric protein features, and participant IDs as
#             unique row names. Values must be finite and non-missing.
#   immlist: uniquely named list of pathway-to-protein mappings.
#   L1, L2: named stage-1 and stage-2 regression learner specifications.
#   FUS: names of the fixed stage-2 fusion methods.
#   S, PTRAIN, KOUT, KIN, NC, TLIM, KEEP: original modeling settings.
#     Defaults retain the source configuration: seed 666, 80:20 split,
#     five outer and five inner folds, 30 cores, 600-second fit timeout,
#     and fitted-model retention disabled.
#
# One shared stage-1 algorithm is fitted separately to every pathway. The shared
# stage-1 algorithm and stage-2 integrator are selected jointly using pooled
# training outer-fold out-of-fold predictions. The held-out test set is then
# evaluated, and the selected combination is locked for full-cohort nested OOF.
# No feature selection or imputation is performed within this modeling block.
#
# The supplied source includes R2 as a metric and the final tie-break criterion;
# these calculations are preserved here together with MAE, PCC, and RMSE.
#
run_pathway_clock <- function(all_data, immlist, L1 = list(Ridge = lrn("regr.glmnet", alpha = 0), LinearSVR =
  lrn("regr.liblinear"), LightGBM = lrn("regr.lightgbm"), FNN = lrn("regr.fnn"), LASSO = lrn("regr.glmnet",
  alpha = 1)), L2 = list(Ridge = lrn("regr.glmnet", alpha = 0), LM = lrn("regr.lm"), NotNegRidge = lrn("regr.glmnet",
  alpha = 0, lower.limits = 0)), FUS = c("EqualWeight", "MAEWeight"), S = 666L, PTRAIN = 0.8, KOUT = 5L, KIN =
  5L, NC = 30L, TLIM = 600, KEEP = FALSE) {

  ## ---------- audit, packages and mapped modules ----------
  normL <- function(x, s) {
    if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) || any(!nzchar(trimws(names(x)))) ||
      anyDuplicated(names(x))) stop(s, " must be a uniquely named non-empty list")
    z <- lapply(x, \(q) tryCatch({
      if (is.character(q)) {
        if (length(q) != 1L) stop("learner ID must have length 1")
        q <- lrn(if (startsWith(q, "regr.")) q else paste0("regr.", q))
      } else if (is.list(q) && !inherits(q, "Learner")) {
        id <- q$id
        q$id <- NULL
        if (!is.character(id) || length(id) != 1L) stop("learner specification needs one id")
        id <- if (startsWith(id, "regr.")) id else paste0("regr.", id)
        q <- do.call(lrn, c(list(id), q))
      }
      if (!inherits(q, "Learner") || q$task_type != "regr") stop("not a regression learner")
      q
    }, error = identity))
    ok <- !vapply(z, inherits, logical(1), "error")
    st <- data.table(name = names(x), registered = ok, registry_error = vapply(z, \(q) if (inherits(q, "error")) conditionMessage(q) else "",
      ""))
    if (any(!ok)) warning("Skipping unregistered ", s, ": ", paste(st$name[!ok], collapse = ", "))
    z <- z[ok]
    if (!length(z)) stop("No registered learner in ", s)
    list(learners = z, status = st)
  }
  cap <- function(z) {
    fb <- lrn("regr.featureless")
    if ("when" %in% names(formals(z$encapsulate))) z$encapsulate("callr", fallback = fb, when = \(cond, stage) FALSE) else z$encapsulate("callr",
      fallback = fb)
    z$timeout <- c(train = TLIM, predict = TLIM)
    z
  }
  a1 <- normL(L1, "L1")
  a2 <- normL(L2, "L2")
  L1 <- lapply(a1$learners, cap)
  L2 <- lapply(a2$learners, cap)
  l1_registry_status <- a1$status
  l2_registry_status <- a2$status
  setnames(l1_registry_status, "name", "stage1")
  setnames(l2_registry_status, "name", "stage2")
  rm(a1, a2)
  if (any(grepl("__", c(names(L1), names(L2)), fixed = TRUE)) || any(FUS %in% names(L2))) stop("Learner names cannot contain '__' or duplicate mandatory fusion names")
  stopifnot(KOUT >= 2L, KIN >= 2L, NC >= 1L, PTRAIN > 0, PTRAIN < 1, TLIM > 0, names(all_data)[1L] == "label",
    is.numeric(all_data$label), !anyNA(names(all_data)), !anyDuplicated(names(all_data)), is.list(immlist),
    length(immlist) > 0L)
  if (is.null(names(immlist)) || anyNA(names(immlist)) || any(!nzchar(trimws(names(immlist)))) || anyDuplicated(names(immlist))) stop("Module names are missing, blank or duplicated")
  eid <- rownames(all_data)
  if (is.null(eid) || anyNA(eid) || any(!nzchar(eid)) || anyDuplicated(eid)) stop("EID is missing, blank or duplicated")
  ip <- rownames(installed.packages())
  if (!"callr" %in% ip) stop("Missing mandatory package: callr")
  l1miss <- setNames(lapply(L1, \(z) setdiff(z$packages, ip)), names(L1))
  l1_package_status <- data.table(stage1 = names(L1), available = !lengths(l1miss), missing = vapply(l1miss,
    paste, "", collapse = ", "))
  badL1 <- l1_package_status$stage1[!l1_package_status$available]
  if (length(badL1)) warning("Skipping unavailable L1: ", paste(badL1, collapse = ", "))
  L1 <- L1[l1_package_status$stage1[l1_package_status$available]]
  if (!length(L1)) stop("No available L1 learner")
  l2miss <- setNames(lapply(L2, \(z) setdiff(z$packages, ip)), names(L2))
  l2_package_status <- data.table(stage2 = names(L2), available = !lengths(l2miss), missing = vapply(l2miss,
    paste, "", collapse = ", "))
  badL2 <- l2_package_status$stage2[!l2_package_status$available]
  if (length(badL2)) warning("Skipping unavailable L2: ", paste(badL2, collapse = ", "))
  L2 <- L2[l2_package_status$stage2[l2_package_status$available]]
  ann <- lapply(immlist, \(z) unique(as.character(z)))
  mods0 <- lapply(ann, intersect, names(all_data)[ - 1L])
  unmapped <- Map(setdiff, ann, mods0)
  if (any(!lengths(mods0))) stop("No mapped feature: ", paste(names(mods0)[!lengths(mods0)], collapse = ", "))
  used <- unique(unlist(mods0, use.names = FALSE))
  fi <- sprintf("F%05d", seq_along(used))
  mi <- sprintf("M%03d", seq_along(mods0))
  f2i <- setNames(fi, used)
  mods <- lapply(mods0, \(z) unname(f2i[z]))
  names(mods) <- mi
  dat <- all_data[, c("label", used), drop = FALSE]
  names(dat) <- c("label", fi)
  if (!all(vapply(dat, is.numeric, logical(1))) || !all(vapply(dat, \(z) all(is.finite(z)), logical(1)))) stop("label/features must be finite numeric values")
  constant <- used[vapply(dat[ - 1L], \(z) length(unique(z)) < 2L, logical(1))]
  if (length(constant)) warning("Constant features retained because feature selection is disabled: ", paste(constant,
    collapse = ", "))
  feature_name_map <- data.table(original = used, internal = fi)
  module_name_map <- data.table(original = names(mods0), internal = mi)
  module_feature_map <- rbindlist(Map(\(m, o, i) data.table(module = m, original = o, internal = i), names(mods0),
    mods0, mods))
  reserved <- c("eid", "label", "fold", "prediction", "immune_age", "age_gap_raw", "age_acceleration")
  if (any(names(mods0) %in% reserved)) stop("Module names conflict with output columns: ", paste(intersect(names(mods0),
    reserved), collapse = ", "))
  mapping <- data.table(module = names(mods0), module_internal = mi, annotated = lengths(ann), mapped = lengths(mods0),
    excluded_unmapped = lengths(unmapped))
  y <- dat$label
  n <- length(y)
  tasks <- Map(\(z, i) TaskRegr$new(paste0("module", i), as.data.table(dat[, c("label", z), drop = FALSE]),
    target = "label"), mods, seq_along(mods))
  names(tasks) <- names(mods)
  stopifnot(all(vapply(tasks, \(z) identical(z$row_ids, seq_len(n)), logical(1))))

  ## ---------- folds, metrics and helpers ----------
  mkfold <- function(a, k, seed, b = 10L) {
    if (k > length(a)) stop("folds > samples")
    set.seed(seed)
    q <- unique(quantile(a, seq(0, 1, length.out = b + 1L), type = 8))
    g <- if (length(q) > 2L) cut(a, q, include.lowest = TRUE) else factor(rep(1L, length(a)))
    f <- integer(length(a))
    for (ii in split(seq_along(a), g, drop = TRUE)) {
      v <- rep_len(seq_len(k), length(ii))
      f[ii] <- v[sample.int(length(ii))]
    }
    if (any(tabulate(f, k) == 0L)) {
      v <- rep_len(seq_len(k), length(a))
      f <- v[sample.int(length(a))]
    }
    f
  }
  mksplit <- function(a, p, seed, b = 10L) {
    set.seed(seed)
    q <- unique(quantile(a, seq(0, 1, length.out = b + 1L), type = 8))
    g <- if (length(q) > 2L) cut(a, q, include.lowest = TRUE) else factor(rep(1L, length(a)))
    sort(unlist(lapply(split(seq_along(a), g, drop = TRUE), \(ii) ii[sample.int(length(ii), max(1L, floor(p *
      length(ii))))]), use.names = FALSE))
  }
  met <- function(a, p) {
    s <- sum((a - mean(a)) ^ 2)
    data.table(N = length(a), MAE = mean(abs(a - p)), PCC = if (length(a) > 1L && sd(a) > 0 && sd(p) > 0) cor(a,
      p) else NA_real_, R2 = if (s > 0) 1 - sum((a - p) ^ 2) / s else NA_real_, RMSE = sqrt(mean((a - p) ^
      2)))
  }
  restoreM <- function(M) {
    j <- match(colnames(M), module_name_map$internal)
    if (anyNA(j)) stop("Unknown internal module name")
    colnames(M) <- module_name_map$original[j]
    M
  }
  restoreS <- function(D) {
    D <- copy(D)
    D$module <- module_name_map$original[match(D$module_internal, module_name_map$internal)]
    if (anyNA(D$module)) stop("Unknown module in stage-1 results")
    setcolorder(D, c("module", setdiff(names(D), "module")))
    D
  }
  restoreW <- function(W) {
    if (is.null(W) || !nrow(W)) return(data.table())
    W <- copy(W)
    W$module <- module_name_map$original[match(W$module_internal, module_name_map$internal)]
    if (anyNA(W$module)) stop("Unknown module in weights")
    setcolorder(W, c("module", setdiff(names(W), "module")))
    W
  }
  mperf <- function(a, M, set, stage1, f = NULL) {
    one <- function(ii, v) rbindlist(lapply(seq_len(ncol(M)), \(j) cbind(data.table(set = set, module = colnames(M)[j],
      stage1 = stage1, fold = v), met(a[ii], M[ii, j]))))
    if (is.null(f)) one(seq_along(a), NA_integer_) else rbindlist(c(lapply(sort(unique(f)), \(v) one(which(f ==
      v), v)), list(one(seq_along(a), NA_integer_))))
  }
  MAP <- function(x, f, ...) {
    if (NC == 1L || .Platform$OS.type == "windows") lapply(x, f, ...) else mclapply(x, f, ..., mc.cores =
      min(NC, length(x)), mc.set.seed = TRUE, mc.preschedule = FALSE)
  }
  badfit <- function(l) {
    e <- tryCatch(l$errors, error = \(e) NULL)
    if ((is.numeric(e) && any(e > 0)) || (!is.null(e) && !is.numeric(e) && length(e) > 0L)) return(TRUE)
    z <- tryCatch(l$log, error = \(e) NULL)
    if (is.null(z) || !NROW(z)) return(FALSE)
    if ("class" %in% names(z) && any(tolower(as.character(z$class)) == "error")) return(TRUE)
    x <- if ("msg" %in% names(z)) z$msg else if ("condition" %in% names(z)) vapply(z$condition, \(q) if (inherits(q,
      "condition")) conditionMessage(q) else paste(as.character(q), collapse = " "), "") else unlist(z, use.names =
      FALSE)
    any(grepl("using fallback learner", paste(x, collapse = " "), ignore.case = TRUE))
  }
  chk <- function(l) {
    if (is.null(l$model) || badfit(l)) stop("Learner failed or timed out: ", l$id, call. = FALSE)
    invisible(l)
  }
  trn <- function(l, t, id = NULL) {
    l$train(t, row_ids = id)
    chk(l)
  }
  pp <- function(l, t, id) {
    q <- l$predict(t, row_ids = id)
    chk(l)
    j <- match(id, q$row_ids)
    if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Invalid/misaligned stage-1 row_ids")
    p <- q$response[j]
    if (length(p) != length(id) || any(!is.finite(p))) stop("Invalid stage-1 predictions")
    p
  }

  ## ---------- stage 1: one shared learner is fitted separately to every module ----------
  oof1 <- function(task, tr, base, inner, seed) {
    p <- numeric(length(tr))
    for (v in seq_len(KIN)) {
      va <- which(inner == v)
      l <- base$clone(deep = TRUE)
      l$reset()
      set.seed(seed + v)
      trn(l, task, tr[inner != v])
      p[va] <- pp(l, task, tr[va])
    }
    if (any(!is.finite(p))) stop("Invalid inner OOF predictions")
    p
  }
  fitmods <- function(tr, te, base, stage1, inner = mkfold(y[tr], KIN, 1L), seed = 1L, keep = FALSE) {
    z <- MAP(seq_along(tasks), \(m) tryCatch({
      t0 <- proc.time()[3L]
      p0 <- oof1(tasks[[m]], tr, base, inner, seed + 1000L * m)
      s0 <- unname(proc.time()[3L] - t0)
      t0 <- proc.time()[3L]
      l <- base$clone(deep = TRUE)
      l$reset()
      set.seed(seed + 90000L + m)
      trn(l, tasks[[m]], tr)
      p1 <- pp(l, tasks[[m]], te)
      s1 <- unname(proc.time()[3L] - t0)
      list(train = p0, test = p1, row = cbind(data.table(stage1 = stage1, inner_seconds = s0, refit_seconds =
        s1), met(y[tr], p0)), model = if (keep) l else NULL)
    }, error = identity))
    ok <- vapply(z, \(q) is.list(q) && !inherits(q, "condition") && all(c("train", "test", "row") %in% names(q)),
      logical(1))
    if (any(!ok)) {
      m <- which(!ok)[1L]
      msg <- if (inherits(z[[m]], "condition")) conditionMessage(z[[m]]) else paste(as.character(z[[m]]),
        collapse = " ")
      stop("Module ", names(tasks)[m], ": ", msg)
    }
    z0 <- do.call(cbind, lapply(z, `[[`, "train"))
    z1 <- do.call(cbind, lapply(z, `[[`, "test"))
    colnames(z0) <- colnames(z1) <- names(tasks)
    sc <- rbindlist(Map(\(q, m) {
      x <- copy(q$row)
      x$module_internal <- m
      x
    }, z, names(tasks)))
    setcolorder(sc, c("module_internal", "stage1", setdiff(names(sc), c("module_internal", "stage1"))))
    list(train = z0, test = z1, scores = sc, models = if (keep) setNames(lapply(z, `[[`, "model"), names(tasks)) else NULL)
  }

  ## ---------- stage 2: fixed fusions or user-defined learners ----------
  fit2 <- function(z0, z1, a, base, seed, keep = FALSE) {
    set.seed(seed)
    t <- TaskRegr$new("stage2", as.data.table(cbind(label = a, z0)), target = "label")
    l <- base$clone(deep = TRUE)
    l$reset()
    trn(l, t)
    nd <- as.data.table(z1)
    setcolorder(nd, t$feature_names)
    q <- l$predict_newdata(nd, task = t)
    chk(l)
    j <- match(seq_len(nrow(z1)), q$row_ids)
    if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Invalid/misaligned stage-2 row_ids")
    p <- q$response[j]
    if (length(p) != nrow(z1) || any(!is.finite(p))) stop("Invalid stage-2 predictions")
    list(pred = p, model = if (keep) l else NULL)
  }
  fuse2 <- function(z0, z1, a, method) {
    if (!identical(colnames(z0), colnames(z1)) || any(!is.finite(z0)) || any(!is.finite(z1))) stop("Invalid/misaligned module predictions")
    e <- colMeans(abs(sweep(z0, 1L, a, "-")))
    if (method == "EqualWeight") w <- rep(1 / ncol(z0), ncol(z0)) else {
      zero <- e <= sqrt(.Machine$double.eps)
      if (any(zero)) w <- zero / sum(zero) else {
        w <- 1 / (e + sqrt(.Machine$double.eps) * max(1, median(e)))
        w <- w / sum(w)
      }
    }
    names(w) <- colnames(z0)
    list(pred = drop(z1[, names(w), drop = FALSE] %*% w), weights = data.table(module_internal = names(w),
      module_MAE = e[names(w)], weight = as.numeric(w)))
  }
  run2 <- function(z0, z1, a, method, B, seed, keep = FALSE) {
    if (method %in% FUS) {
      q <- fuse2(z0, z1, a, method)
      q$model <- if (keep) list(method = method, weights = q$weights) else NULL
      q
    } else {
      if (is.null(B[[method]])) stop("Unknown stage-2 method: ", method)
      q <- fit2(z0, z1, a, B[[method]], seed, keep)
      q$weights <- NULL
      q
    }
  }
  once <- function(tr, te, stage1, stage2, seed, keep = FALSE, A = L1, B = L2) {
    if (is.null(A[[stage1]])) stop("Unknown stage-1 learner: ", stage1)
    inn <- mkfold(y[tr], KIN, seed)
    z <- fitmods(tr, te, A[[stage1]], stage1, inn, seed, keep)
    q <- run2(z$train, z$test, y[tr], stage2, B, seed + 50000L, keep)
    list(pred = q$pred, module = z$test, weights = q$weights, stage1_scores = z$scores, models = if (keep) list(stage1 =
      z$models, stage2 = q$model) else NULL)
  }

  ## ---------- nested stacking: compare globally shared L1 x stage-2 pipelines ----------
  cvstack <- function(ids, A = L1, B = L2, seed = S, keep = FALSE, methods = c(FUS, names(B))) {
    methods <- unique(methods)
    badm <- setdiff(methods, c(FUS, names(B)))
    if (length(badm)) stop("Unknown stage-2 method: ", paste(badm, collapse = ", "))
    f <- mkfold(y[ids], KOUT, seed)
    g <- rbindlist(lapply(names(A), \(an) data.table(stage1 = an, stage2 = methods, pipeline = paste0(an,
      "__", methods))))
    P <- matrix(NA_real_, length(ids), nrow(g), dimnames = list(NULL, g$pipeline))
    M <- setNames(lapply(names(A), \(z) matrix(NA_real_, length(ids), length(tasks), dimnames = list(NULL,
      names(tasks)))), names(A))
    S1 <- W <- E1 <- E2 <- list()
    fm <- if (keep) vector("list", KOUT) else NULL
    for (v in seq_len(KOUT)) {
      iv <- which(f == v)
      tr <- ids[f != v]
      te <- ids[iv]
      inn <- mkfold(y[tr], KIN, seed + 100L * v)
      if (keep) fm[[v]] <- list()
      for (i in seq_along(A)) {
        an <- names(A)[i]
        z <- tryCatch(fitmods(tr, te, A[[i]], an, inn, seed + 100000L * v + 1000L * i, keep), error = identity)
        if (inherits(z, "error")) {
          E1[[length(E1) + 1L]] <- data.table(fold = v, stage1 = an, message = conditionMessage(z))
          next
        }
        M[[an]][iv, ] <- z$test
        sx <- copy(z$scores)
        sx$fold <- v
        S1[[length(S1) + 1L]] <- sx
        if (keep) fm[[v]][[an]] <- list(stage1 = z$models, stage2 = list())
        for (j in seq_along(methods)) {
          bn <- methods[j]
          pn <- paste0(an, "__", bn)
          q <- tryCatch(run2(z$train, z$test, y[tr], bn, B, seed + 500000L * v + 100L * i + j, keep), error =
            identity)
          if (inherits(q, "error")) {
            E2[[length(E2) + 1L]] <- data.table(fold = v, stage1 = an, stage2 = bn, pipeline = pn, message =
              conditionMessage(q))
            next
          }
          P[iv, pn] <- q$pred
          if (!is.null(q$weights)) W[[length(W) + 1L]] <- cbind(data.table(fold = v, stage1 = an, stage2 =
            bn, pipeline = pn), q$weights)
          if (keep) fm[[v]][[an]]$stage2[[bn]] <- q$model
        }
      }
    }
    ok <- vapply(seq_len(ncol(P)), \(j) all(is.finite(P[, j])), logical(1))
    if (!any(ok)) stop("All shared-stage1 pipelines failed")
    P <- P[, ok, drop = FALSE]
    g <- g[match(colnames(P), g$pipeline)]
    M <- M[unique(g$stage1)]
    if (any(vapply(M, \(z) !all(is.finite(z)), logical(1)))) stop("Selected pipeline has invalid module predictions")
    W <- if (length(W)) rbindlist(W, fill = TRUE) else data.table()
    if (nrow(W)) W <- W[W$pipeline %in% g$pipeline]
    E1 <- if (length(E1)) rbindlist(E1, fill = TRUE) else data.table()
    E2 <- if (length(E2)) rbindlist(E2, fill = TRUE) else data.table()
    sc <- if (length(S1)) rbindlist(S1, fill = TRUE) else data.table()
    list(pred = P, module = M, fold = f, grid = g, weights = W, stage1_errors = E1, stage2_errors = E2, stage1_scores =
      sc, models = fm)
  }
  score <- function(a, P, g) rbindlist(lapply(seq_len(nrow(g)), \(j) cbind(g[as.integer(j)], met(a, P[, j]))))
  scorefold <- function(a, P, g, f) rbindlist(lapply(sort(unique(f)), \(v) cbind(data.table(fold = v), score(a[f ==
    v], P[f == v, , drop = FALSE], g))))

  ## ---------- 80:20 development/test evaluation ----------
  ids <- seq_len(n)
  tr <- mksplit(y, PTRAIN, S)
  te <- setdiff(ids, tr)
  if (!length(te) || length(tr) < KOUT || floor(length(tr) * (KOUT - 1) / KOUT) < KIN) stop("Invalid split/fold settings")
  t0 <- proc.time()[3L]
  sel <- cvstack(tr, L1, L2, S, FALSE)
  selection_seconds <- unname(proc.time()[3L] - t0)
  leaderboard <- score(y[tr], sel$pred, sel$grid)
  setorder(leaderboard, MAE, RMSE, - PCC, - R2, na.last = TRUE)
  leaderboard[, `:=`(role = "selection_secondary", test_refit_attempted = FALSE, test_refit_ok = NA, selected =
    FALSE)]
  test_refit_errors <- list()
  testfit <- NULL
  for (i in seq_len(nrow(leaderboard))) {
    b1 <- leaderboard$stage1[i]
    b2 <- leaderboard$stage2[i]
    q <- tryCatch(once(tr, te, b1, b2, S + 9000L + 100L * i, KEEP, L1, L2), error = identity)
    leaderboard$test_refit_attempted[i] <- TRUE
    if (inherits(q, "error")) {
      leaderboard$test_refit_ok[i] <- FALSE
      test_refit_errors[[length(test_refit_errors) + 1L]] <- data.table(rank = i, stage1 = b1, stage2 = b2,
        pipeline = leaderboard$pipeline[i], message = conditionMessage(q))
    } else {
      leaderboard$test_refit_ok[i] <- TRUE
      BEST1 <- b1
      BEST2 <- b2
      key <- leaderboard$pipeline[i]
      testfit <- q
      break
    }
  }
  if (is.null(testfit)) stop("No ranked pipeline could be refitted on the full training set")
  leaderboard[pipeline == key, selected := TRUE]
  test_refit_errors <- if (length(test_refit_errors)) rbindlist(test_refit_errors) else data.table(rank =
    integer(), stage1 = character(), stage2 = character(), pipeline = character(), message = character())
  fold_performance <- scorefold(y[tr], sel$pred, sel$grid, sel$fold)
  fold_performance$role <- "selection_secondary"
  split_performance <- rbindlist(list(cbind(data.table(set = "train_selection_oof", role = "selection_secondary",
    stage1 = BEST1, stage2 = BEST2, pipeline = key), met(y[tr], sel$pred[, key])), cbind(data.table(set =
    "heldout_test", role = "primary", stage1 = BEST1, stage2 = BEST2, pipeline = key), met(y[te], testfit$pred))))
  train_prediction <- data.frame(eid = eid[tr], label = y[tr], fold = sel$fold, prediction = sel$pred[, key])
  test_prediction <- data.frame(eid = eid[te], label = y[te], prediction = testfit$pred)
  trM <- restoreM(sel$module[[BEST1]])
  teM <- restoreM(testfit$module)
  train_module_prediction <- data.frame(eid = eid[tr], label = y[tr], fold = sel$fold, trM, check.names =
    FALSE)
  test_module_prediction <- data.frame(eid = eid[te], label = y[te], teM, check.names = FALSE)
  train_stage1_scores <- restoreS(sel$stage1_scores)
  train_stage1_scores[, `:=`(selected_global = stage1 == BEST1, scope = "training_outer-fold_inner-OOF")]
  test_stage1_scores <- restoreS(testfit$stage1_scores)
  test_stage1_scores[, `:=`(selected_global = TRUE, scope = "training_inner-OOF_for_test-refit")]
  train_module_performance <- mperf(y[tr], trM, "train_selection_oof", BEST1, sel$fold)
  train_module_performance$role <- "selection_secondary"
  test_module_performance <- mperf(y[te], teM, "heldout_test", BEST1)
  test_module_performance$role <- "heldout_descriptive_secondary"
  selection_fusion_weights_all <- restoreW(sel$weights)
  selection_fusion_weights <- restoreW(if (nrow(sel$weights)) sel$weights[sel$weights$pipeline == key] else data.table())
  test_fusion_weights <- restoreW(if (is.null(testfit$weights)) data.table() else cbind(data.table(stage1 =
    BEST1, stage2 = BEST2, pipeline = key), testfit$weights))

  ## ---------- lock BEST1 + BEST2; full-sample nested OOF ----------
  full <- cvstack(ids, L1[BEST1], L2, S + 20000L, KEEP, methods = BEST2)
  stopifnot(identical(full$grid$pipeline, key))
  fp <- full$pred[, 1L]
  fmout <- restoreM(full$module[[BEST1]])
  full_performance <- rbindlist(list(scorefold(y, full$pred, full$grid, full$fold), cbind(data.table(fold =
    NA_integer_), score(y, full$pred, full$grid))), fill = TRUE)
  full_performance$role <- "post_selection_secondary"
  full_stage1_scores <- restoreS(full$stage1_scores)
  full_stage1_scores[, `:=`(selected_global = TRUE, scope = "full_outer-fold_inner-OOF_locked")]
  module_performance <- mperf(y, fmout, "full_nested_oof", BEST1, full$fold)
  module_performance$role <- "post_selection_secondary"
  full_module_prediction <- data.frame(eid = eid, label = y, fold = full$fold, fmout, check.names = FALSE)
  qmae <- vapply(seq_len(ncol(teM)), \(j) mean(abs(y[te] - teM[, j])), 0)
  stopifnot(all(vapply(list(train_module_performance, test_module_performance, module_performance), \(z) !("stage2" %in%
    names(z)) && all(z$stage1 == BEST1), TRUE)), isTRUE(all.equal(unname(qmae), test_module_performance$MAE[match(colnames(teM),
    test_module_performance$module)], tolerance = 1e-12)))
  rm(qmae)
  full_fusion_weights <- restoreW(full$weights)
  full_oof <- data.frame(eid = eid, label = y, fold = full$fold, immune_age = fp, age_gap_raw = fp - y, age_acceleration =
    resid(lm(fp ~ y)), fmout, check.names = FALSE)

  ## ---------- model outputs in memory ----------
  desc <- function(x) lapply(x, \(z) list(id = z$id, params = z$param_set$values, packages = z$packages, timeout =
    z$timeout, encapsulation = z$encapsulation))
  result <- list(best = list(stage1 = BEST1, stage2 = BEST2, pipeline = key), mapping = mapping, excluded_unmapped =
    unmapped, feature_name_map = feature_name_map, module_name_map = module_name_map, module_feature_map =
    module_feature_map, module_features = mods0, registry_status = list(stage1 = l1_registry_status, stage2 =
    l2_registry_status), package_status = list(stage1 = l1_package_status, stage2 = l2_package_status), split =
    list(train_eid = eid[tr], test_eid = eid[te], selection_seconds = selection_seconds, stage1_errors =
    sel$stage1_errors, stage2_errors = sel$stage2_errors, test_refit_errors = test_refit_errors, leaderboard =
    leaderboard, fold_performance = fold_performance, performance = split_performance, stage1_module_diagnostics =
    train_stage1_scores, test_refit_stage1_module_diagnostics = test_stage1_scores, train_prediction = train_prediction,
    test_prediction = test_prediction, train_module_prediction = train_module_prediction, train_module_performance =
    train_module_performance, test_module_prediction = test_module_prediction, test_module_performance =
    test_module_performance, fusion_weights_all = selection_fusion_weights_all, fusion_weights = selection_fusion_weights,
    test_fusion_weights = test_fusion_weights, test_models = if (KEEP) testfit$models else NULL), full_cv =
    list(oof = full_oof, performance = full_performance, module_prediction = full_module_prediction, module_performance =
    module_performance, stage1_module_diagnostics = full_stage1_scores, fusion_weights = full_fusion_weights,
    stage1_errors = full$stage1_errors, stage2_errors = full$stage2_errors, fold_models = full$models), config =
    list(seed = S, train_fraction = PTRAIN, outer_folds = KOUT, inner_folds = KIN, cores = NC, per_fit_timeout_seconds =
    TLIM, stage1_selection = "one globally shared L1 jointly selected with L2 by pooled training outer-OOF MAE",
    pipeline_tiebreak = c("RMSE", "PCC(desc)", "R2(desc)"), refit_rule = "first ranked pipeline successfully refitted without using test labels",
    stage2_selection = "joint with shared L1 by pooled training outer-OOF MAE", metrics = c("MAE", "PCC",
    "R2", "RMSE"), mandatory_fusions = FUS, mae_weighting = "normalized inverse inner-OOF module MAE", external_feature_selection =
    "none", unmapped_features = "excluded", imputation = "none", test_refit_stage1_scores = "training inner-OOF scores, not test-set candidate comparison",
    test_role = "primary", full_cv_role = "post-selection secondary"), learners = list(stage1 = desc(L1),
    stage2 = desc(L2)), session = sessionInfo())
  print(mapping)
  print(l1_registry_status)
  print(l1_package_status)
  print(l2_registry_status)
  print(l2_package_status)
  print(leaderboard)
  print(result$best)
  print(test_module_performance)
  print(split_performance)
  print(full_performance[is.na(full_performance$fold)])
  return(result)
}


# ==============================================================================
# 03 | Prediction performance and pathway contributions
# ==============================================================================

# ============================================================================
# Held-out prediction performance and additive pathway contributions
# ============================================================================
# Input objects: the pathway model result and its matching training-time
# Reactome v97 bundle. The annotation module must define classify_reactome8().
# This block uses the existing held-out predictions without refitting models.
# It returns MAE, RMSE, and PCC with paired participant-bootstrap intervals,
# pathway contributions, fractional contributions to eight Reactome classes,
# the protein-count correlation test, and top-ten protein membership counts.
#
# Additive contributions are centered on the mean held-out pathway prediction.
# The original affine equation must reproduce every integrated prediction within
# the source tolerance before contributions are accepted.
ia_clock_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
ia_clock_need <- function(x, columns, label) {
  ia_clock_assert(all(columns %in% names(x)), paste0(label, " lacks: ", paste(setdiff(columns, names(x)),
    collapse = ", ")))
}
ia_clock_affine <- function(R, X, f, reconstruction_tol = 1e-7, condition_limit = 1e6) {

  # The reconstruction target is the saved INTEGRATED PREDICTION f, not age.
  # Coefficients are accepted only for known affine integrations, with exact
  # reproduction and unique, numerically stable varying coefficients.
  S <- R$split
  modules <- colnames(X)
  m <- ncol(X)
  n <- nrow(X)
  method <- as.character(R$best$stage2)
  constant <- vapply(seq_len(m), function(j) all(X[, j] == X[1L, j]), logical(1))
  recovery <- list(status = "not_verified", stage2 = method, participants = n, pathways = m, constant_pathways =
    modules[constant], original_coefficient_nonidentifiable = character(), prediction_tolerance_years = reconstruction_tol,
    condition_limit = condition_limit)
  verify <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }
  W <- if (is.null(S$test_fusion_weights)) data.frame() else as.data.frame(S$test_fusion_weights)
  if (method %in% c("EqualWeight", "MAEWeight")) {
    if (nrow(W)) {
      verify(all(c("module", "weight") %in% names(W)), "test_fusion_weights lacks module or weight.")
      if ("stage1" %in% names(W)) verify(all(W$stage1 == R$best$stage1), "Saved fusion weights use a different stage-1 learner.")
      if ("stage2" %in% names(W)) verify(all(W$stage2 == method), "Saved fusion weights use a different integration method.")
      verify(!anyNA(W$module) && !anyDuplicated(W$module) && setequal(W$module, modules), "Saved fusion weights do not match the model pathways.")
      beta <- stats::setNames(as.numeric(W$weight[match(modules, W$module)]), modules)
      origin <- "saved held-out refit fusion weights"
    } else if (method == "EqualWeight") {
      beta <- stats::setNames(rep(1 / m, m), modules)
      origin <- "exact equal-weight integration rule"
    } else {
      verify(FALSE, "Original MAEWeight weights are missing; inverse held-out MAE cannot replace them.")
    }
    verify(all(is.finite(beta)) && all(beta >= 0) && abs(sum(beta) - 1) <= 1e-8, "Invalid saved fusion weights.")
    intercept <- 0
  } else {
    known <- c(Ridge = "regr.glmnet", LM = "regr.lm", NotNegRidge = "regr.glmnet")
    id <- R$learners$stage2[[method]]$id
    verify(method %in% names(known) && length(id) == 1L && identical(as.character(id), unname(known[method])),
      "Unsupported or unverified integration learner: a known affine model is required.")
    center <- colMeans(X)
    spread <- apply(X, 2L, stats::sd)
    varying <- !constant
    center[constant] <- X[1L, constant]
    recovery$original_coefficient_nonidentifiable <- modules[constant]
    recovery$constant_coefficient_convention <- if (any(constant)) {
      paste("Effective zero; the original coefficient is not identifiable.", "Its centered contribution is zero and its uncentered term is absorbed in the intercept.")
    } else "Not applicable"
    recovery$varying_pathways <- sum(varying)
    recovery$columns <- sum(varying) + 1L
    verify(n > sum(varying) + 1L, "Coefficient recovery requires more participants than varying pathway columns plus the intercept.")
    verify(any(varying) && all(is.finite(spread[varying])) && all(spread[varying] > 0), "Varying pathway predictions cannot be stably standardized.")
    U <- cbind(Intercept = 1, sweep(sweep(X[, varying, drop = FALSE], 2L, center[varying], "-"), 2L, spread[varying],
      "/"))
    qrX <- qr(U, tol = 1e-10, LAPACK = FALSE)
    cond <- kappa(U, exact = TRUE)
    recovery$rank <- qrX$rank
    recovery$condition_number <- cond
    verify(qrX$rank == ncol(U) && is.finite(cond) && cond <= condition_limit, paste("Varying pathway predictions do not identify stable unique coefficients.",
      "Provide the original integration coefficients from the saved model;", "no approximate surrogate or pseudo-inverse is used."))
    q <- qr.coef(qrX, f)
    beta <- stats::setNames(rep(0, m), modules)
    beta[varying] <- unname(q[ - 1L] / spread[varying])
    intercept <- unname(q[1L] - sum(beta * center))
    if (method == "NotNegRidge") verify(min(beta) >= - 1e-9 * max(1, max(abs(beta))), "Recovered coefficients violate the non-negative integration constraint.")
    origin <- if (any(constant)) {
      "unique reconstruction of varying terms; constant terms use effective zero and an absorbed intercept"
    } else "unique reconstruction of the known affine equation from saved pathway and integrated predictions"
  }
  recovery$origin <- origin
  verify(all(is.finite(beta)) && is.finite(intercept), "Non-finite integration coefficients.")
  recovery$max_absolute_prediction_error_years <- max(abs(intercept + drop(X %*% beta) - f))
  verify(recovery$max_absolute_prediction_error_years <= reconstruction_tol, "Affine coefficients do not reproduce saved integrated predictions within 1e-7 years; contribution analysis stopped.")
  center <- colMeans(X)
  center[constant] <- X[1L, constant]
  Phi <- sweep(sweep(X, 2L, center, "-"), 2L, beta, "*")
  recovery$max_centered_sum_error_years <- max(abs(rowSums(Phi) - (f - mean(f))))
  verify(recovery$max_centered_sum_error_years <= 2 * reconstruction_tol, "Centered additive contributions fail the participant-level sum identity; contribution analysis stopped.")
  importance <- colMeans(abs(Phi))
  verify(sum(importance) > 0, "All pathway contributions are zero.")
  recovery$status <- "verified"
  list(beta = beta, intercept = intercept, center = center, importance = importance, mean_signed = colMeans(Phi),
    recovery = recovery)
}
ia_clock_bootstrap <- function(y, P, B, seed) {

  # Each multinomial count vector is exactly a size-n bootstrap of participants.
  # The SAME count vector is used for all 59 clocks and all three metrics.
  old_rng_kind <- RNGkind()
  old_rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  on.exit({
    do.call(RNGkind, as.list(old_rng_kind))
    if (is.null(old_rng)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_rng, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  n <- length(y)
  k <- ncol(P)
  E <- sweep(P, 1L, y, "-")
  a <- y - mean(y)
  p <- sweep(P, 2L, colMeans(P), "-")
  G <- cbind(abs(E), E ^ 2, p, p ^ 2, sweep(p, 1L, a, "*"), a, a ^ 2)
  boot <- array(NA_real_, dim = c(B, k, 3L), dimnames = list(NULL, colnames(P), c("MAE", "RMSE", "PCC")))
  for (start in seq.int(1L, B, by = 32L)) {
    ii <- start:min(start + 31L, B)
    W <- stats::rmultinom(length(ii), n, rep(1 / n, n))
    storage.mode(W) <- "double"
    v <- crossprod(W, G) / n
    boot[ii, , "MAE"] <- v[, seq_len(k), drop = FALSE]
    boot[ii, , "RMSE"] <- sqrt(pmax(v[, k + seq_len(k), drop = FALSE], 0))
    mp <- v[, 2L * k + seq_len(k), drop = FALSE]
    var_p <- pmax(v[, 3L * k + seq_len(k), drop = FALSE] - mp ^ 2, 0)
    ma <- v[, 5L * k + 1L]
    var_a <- pmax(0, v[, 5L * k + 2L] - ma ^ 2)
    den <- sqrt(sweep(var_p, 1L, var_a, "*"))
    covariance <- v[, 4L * k + seq_len(k), drop = FALSE] - sweep(mp, 1L, ma, "*")
    rho <- covariance / den
    rho[den <= 0 | !is.finite(rho)] <- NA_real_
    rho[is.finite(rho)] <- pmin(1, pmax( - 1, rho[is.finite(rho)]))
    boot[ii, , "PCC"] <- rho
    if (start == 1L || max(ii) == B || (start - 1L) %% 512L == 0L) {
      message("Held-out paired bootstrap: ", max(ii), "/", B, " resamples")
    }
  }
  ci <- array(NA_real_, dim = c(k, 3L, 2L), dimnames = list(colnames(P), c("MAE", "RMSE", "PCC"), c("low",
    "high")))
  valid <- matrix(0L, k, 3L, dimnames = dimnames(ci)[1:2])
  for (j in seq_len(k)) for (metric in dimnames(ci)[[2L]]) {
    z <- boot[, j, metric]
    valid[j, metric] <- sum(is.finite(z))

    # Do not silently condition the interval on successful resamples.
    if (all(is.finite(z))) ci[j, metric, ] <- unname(stats::quantile(z, c(.025, .975), type = 7))
  }
  if (any(valid[, "PCC"] < B)) warning("At least one clock has undefined bootstrap PCC values: its PCC interval is NA; see PCC_boot_valid.",
    call. = FALSE)
  list(ci = ci, valid = valid)
}
ia_clock_protein_count_test <- function(contribution) {

  # Exploratory continuous protein-count test; no protein-count groups are tested.
  d <- contribution
  ia_clock_need(d, c("module", "n_proteins", "contribution_pct"), "Pathway contribution")
  ia_clock_assert(nrow(d) == 58L && !anyDuplicated(d$module), "The contribution table must contain all 58 unique pathways.")
  ia_clock_assert(all(is.finite(d$n_proteins) & d$n_proteins > 0 & d$n_proteins == round(d$n_proteins)), "Pathway protein counts must be positive integers.")
  ia_clock_assert(all(is.finite(d$contribution_pct) & d$contribution_pct >= 0) && abs(sum(d$contribution_pct) -
    100) < .01, "Pathway contribution percentages must be finite, nonnegative and sum to 100.")
  estimable <- length(unique(d$n_proteins)) > 1L && length(unique(d$contribution_pct)) > 1L
  test <- if (estimable) stats::cor.test(d$n_proteins, d$contribution_pct, method = "spearman", alternative =
    "two.sided", exact = FALSE) else NULL
  pvalue <- if (is.null(test)) NA_real_ else unname(test$p.value)
  data.frame(n_pathways = nrow(d), p_value = pvalue, test = "Two-sided Spearman correlation; asymptotic P (ties allowed)",
    interpretation = "Nominal exploratory P; pathways can share proteins", contribution_source = "Original held-out-test mean absolute centered additive contribution (%)",
    stringsAsFactors = FALSE)
}
ia_clock_top_pathway_memberships <- function(result, bundle, contribution) {

  # Numerical membership summary for the ten largest pathway contributions.
  # Protein identities and membership counts follow the original contribution
  # analysis; no network coordinates or graphical objects are constructed.
  R <- result
  Z <- bundle
  D <- contribution
  feature_edges <- unique(do.call(rbind, lapply(D$module, function(k) {
    data.frame(feature = as.character(R$module_features[[k]]), module = k, stringsAsFactors = FALSE)
  })))
  feature_nodes <- data.frame(feature = sort(unique(feature_edges$feature)), stringsAsFactors = FALSE)
  feature_nodes$protein <- feature_nodes$feature
  zm <- as.data.frame(Z$membership)
  if (all(c("feature", "symbol") %in% names(zm))) {
    symbols <- split(as.character(zm$symbol), as.character(zm$feature))
    feature_nodes$protein <- vapply(feature_nodes$feature, function(k) {
      v <- unique(symbols[[k]])
      v <- v[!is.na(v) & nzchar(v)]
      ia_clock_assert(length(v) <= 1L, paste("Ambiguous protein identity for feature", k))
      if (length(v)) v else k
    }, character(1))
  }
  all_edges <- unique(data.frame(protein = feature_nodes$protein[match(feature_edges$feature, feature_nodes$feature)],
    module = feature_edges$module, stringsAsFactors = FALSE))
  pathways <- D[seq_len(10L), c("module", "id", "pathway", "rank", "contribution_pct"), drop = FALSE]
  pathways <- pathways[order(pathways$rank, pathways$module, method = "radix"), , drop = FALSE]
  edges <- unique(all_edges[all_edges$module %in% pathways$module, c("protein", "module"), drop = FALSE])
  modules <- pathways$module
  proteins <- sort(unique(edges$protein), method = "radix")
  incidence <- matrix(FALSE, length(proteins), length(modules), dimnames = list(proteins, modules))
  incidence[cbind(match(edges$protein, proteins), match(edges$module, modules))] <- TRUE
  degree <- rowSums(incidence)
  sizes <- colSums(incidence)
  private_n <- colSums(incidence[degree == 1L, , drop = FALSE])
  pathways$n_proteins <- unname(sizes)
  pathways$n_private_proteins <- unname(private_n)
  pathways$n_shared_proteins <- unname(sizes - private_n)
  protein_counts <- data.frame(protein = proteins, n_pathways = unname(degree), stringsAsFactors = FALSE)
  list(totals = data.frame(n_pathways = length(modules), n_unique_proteins = length(proteins), n_memberships =
    nrow(edges), n_shared_proteins = sum(degree > 1L)), pathways = pathways, proteins = protein_counts, membership =
    edges)
}
analyze_clock_performance <- function(result, bundle, bootstrap_B = 2000L, seed = 666L, expected_n = 10602L) {
  ia_clock_assert(exists("classify_reactome8", mode = "function", inherits = TRUE), "The Reactome annotation module must define classify_reactome8().")
  ia_clock_assert(length(bootstrap_B) == 1L && is.finite(bootstrap_B) && bootstrap_B == as.integer(bootstrap_B) &&
    bootstrap_B >= 1000L, "bootstrap_B must be an integer >= 1000; the original setting is 2000.")
  ia_clock_assert(length(seed) == 1L && is.finite(seed) && seed >= 0 && seed < .Machine$integer.max && seed ==
    as.integer(seed), "Invalid bootstrap seed.")
  ia_clock_assert(is.null(expected_n) || (length(expected_n) == 1L && is.finite(expected_n) && expected_n >=
    20L && expected_n == as.integer(expected_n)), "expected_n must be NULL or an integer >= 20.")
  message("Checking held-out predictions, matching Reactome v97 definitions and eight-class annotations.")
  R <- result
  Z <- bundle
  ia_clock_need(R, c("best", "split", "module_features", "learners"), "result")
  ia_clock_need(R$best, c("stage1", "stage2", "pipeline"), "result$best")
  ia_clock_assert(all(vapply(R$best[c("stage1", "stage2", "pipeline")], function(z) length(z) == 1L && !is.na(z) &&
    nzchar(as.character(z)), logical(1))), "Invalid selected pipeline metadata.")
  S <- R$split
  ia_clock_need(S, c("test_prediction", "test_module_prediction", "test_eid", "train_eid"), "result$split")
  T <- as.data.frame(S$test_prediction)
  M <- as.data.frame(S$test_module_prediction)
  ia_clock_need(T, c("eid", "label", "prediction"), "test_prediction")
  ia_clock_need(M, c("eid", "label"), "test_module_prediction")
  ia_clock_assert(!anyDuplicated(names(T)) && !anyDuplicated(names(M)), "Duplicate prediction column names.")
  T$eid <- as.character(T$eid)
  M$eid <- as.character(M$eid)
  ids <- list(T$eid, M$eid, as.character(S$test_eid), as.character(S$train_eid))
  ia_clock_assert(all(vapply(ids, function(z) length(z) > 0L && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z),
    logical(1))), "Missing, blank or duplicated split/prediction IDs.")
  ia_clock_assert(setequal(T$eid, M$eid) && setequal(T$eid, as.character(S$test_eid)), "Held-out prediction IDs disagree across result objects.")
  ia_clock_assert(!length(intersect(T$eid, as.character(S$train_eid))), "Train/test IDs overlap.")
  M <- M[match(T$eid, M$eid), , drop = FALSE]
  ia_clock_assert(is.numeric(T$label) && is.numeric(M$label) && isTRUE(all.equal(T$label, M$label, tolerance =
    1e-12)), "Actual ages disagree after ID alignment.")
  modules <- names(R$module_features)
  m <- length(modules)
  n <- nrow(T)
  ia_clock_assert(n >= 20L && (is.null(expected_n) || n == expected_n), paste0("Held-out N is ", n, "; expected ",
    if (is.null(expected_n)) "at least 20" else expected_n, ". Check the matching model result; full-sample OOF cannot substitute for test predictions."))
  ia_clock_assert(m == 58L && !anyNA(modules) && !anyDuplicated(modules), "This analysis requires the original 58 uniquely named pathways.")
  ia_clock_assert(setequal(modules, setdiff(names(M), c("eid", "label", "fold"))), "Model pathways and held-out pathway prediction columns disagree.")
  ia_clock_assert(is.list(Z$sets) && setequal(modules, names(Z$sets)), "The bundle and model have different retained pathways; use the training-time v97 bundle.")
  ia_clock_assert(all(vapply(modules, function(k) length(R$module_features[[k]]) > 0L && !anyDuplicated(R$module_features[[k]]) &&
    all(R$module_features[[k]] %in% Z$sets[[k]]), logical(1))), "Model protein definitions and Reactome sets disagree.")
  ia_clock_assert(is.numeric(T$prediction) && all(vapply(M[, modules, drop = FALSE], is.numeric, logical(1))),
    "Predictions must be numeric.")
  X <- as.matrix(M[, modules, drop = FALSE])
  y <- T$label
  P <- cbind(Integrated = T$prediction, X)
  ia_clock_assert(all(is.finite(y)) && all(is.finite(P)) && stats::sd(y) > 0 && stats::sd(P[, 1L]) > 0, "Non-finite predictions/ages or a constant age/integrated prediction.")
  ann <- classify_reactome8(Z, modules)
  ia_clock_assert(identical(as.character(ann$module), modules), "Classifier changed the pathway order.")
  classes <- attr(ann, "classes")
  ia_clock_assert(is.data.frame(classes) && nrow(classes) == 8L && identical(as.character(classes$class_id),
    paste0("C", 1:8)), "Missing or inconsistent eight-class metadata.")
  cc <- classes$class_id
  E <- sweep(P, 1L, y, "-")
  estimates <- cbind(MAE = colMeans(abs(E)), RMSE = sqrt(colMeans(E ^ 2)), PCC = vapply(seq_len(ncol(P)),
    function(j) {
    if (stats::sd(P[, j]) > 0) stats::cor(y, P[, j]) else NA_real_
  }, numeric(1)))

  # Saved aggregate tables are checks, never substitutes for saved predictions.
  same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-7, check.attributes =
    FALSE))
  if (!is.null(S$test_module_performance) && NROW(S$test_module_performance)) {
    Q <- as.data.frame(S$test_module_performance)
    ia_clock_need(Q, c("module", "set", "stage1", "N", "MAE", "RMSE", "PCC"), "test_module_performance")
    ia_clock_assert(!anyDuplicated(Q$module) && setequal(Q$module, modules), "Stored pathway metric rows disagree.")
    Q <- Q[match(modules, Q$module), , drop = FALSE]
    ia_clock_assert(all(Q$set == "heldout_test") && all(Q$stage1 == R$best$stage1) && all(Q$N == n) && all(vapply(c("MAE",
      "RMSE", "PCC"), function(k) same(Q[[k]], estimates[ - 1L, k]), logical(1))), "Stored pathway held-out metrics disagree with participant predictions.")
  }
  if (!is.null(S$performance) && NROW(S$performance)) {
    Q <- as.data.frame(S$performance)
    ia_clock_need(Q, c("set", "pipeline", "N", "MAE", "RMSE", "PCC"), "split$performance")
    Q <- Q[Q$set == "heldout_test", , drop = FALSE]
    ia_clock_assert(nrow(Q) == 1L && Q$pipeline == R$best$pipeline && Q$N == n && all(vapply(c("MAE", "RMSE",
      "PCC"), function(k) same(Q[[k]], estimates[1L, k]), logical(1))), "Stored integrated held-out metrics disagree with participant predictions.")
  }
  message("Verifying the original affine integration equation and centered additive contributions.")
  A <- ia_clock_affine(R, X, P[, 1L])
  D <- ann
  D$N <- n
  D$n_proteins <- lengths(R$module_features[modules])
  D$coefficient <- unname(A$beta)
  D$coefficient_identified <- !D$module %in% A$recovery$original_coefficient_nonidentifiable
  D$background_age <- unname(A$center)
  D$contribution_years <- unname(A$importance)
  D$contribution_pct <- 100 * D$contribution_years / sum(D$contribution_years)
  D$mean_signed_contribution <- unname(A$mean_signed)
  D <- D[order( - D$contribution_years, D$id), , drop = FALSE]
  rownames(D) <- NULL
  D$rank <- seq_len(m)
  H <- as.matrix(D[, cc, drop = FALSE])
  ia_clock_assert(all(H %in% c(FALSE, TRUE)) && all(rowSums(H) >= 1L), "Invalid class membership matrix.")
  F <- sweep(H, 1L, rowSums(H), "/")
  L <- data.frame(module = rep(D$module, each = 8L), id = rep(D$id, each = 8L), rank = rep(D$rank, each =
    8L), class_id = rep(cc, m), fraction = as.vector(t(F)), stringsAsFactors = FALSE)
  L <- L[L$fraction > 0, , drop = FALSE]
  L$allocated_pct <- D$contribution_pct[match(L$module, D$module)] * L$fraction
  L$allocated_years <- D$contribution_years[match(L$module, D$module)] * L$fraction
  rownames(L) <- NULL
  G <- classes
  G$n_pathways <- colSums(H)
  G$fractional_pathway_count <- colSums(F)
  G$contribution_years <- colSums(sweep(F, 1L, D$contribution_years, "*"))
  G$contribution_pct <- colSums(sweep(F, 1L, D$contribution_pct, "*"))
  ia_clock_assert(abs(sum(G$contribution_pct) - 100) < 1e-8 && abs(sum(G$contribution_years) - sum(D$contribution_years)) <
    1e-8, "Fractional eight-class allocation does not conserve the total contribution.")
  boot <- ia_clock_bootstrap(y, P, as.integer(bootstrap_B), as.integer(seed))
  Q <- data.frame(clock_type = c("Integrated", rep("Pathway", m)), module = c("Integrated", modules), id =
    c(NA_character_, ann$id), pathway = c("Integrated immune age", ann$pathway), primary_class = c(NA_character_,
    ann$primary_class), classes = c(NA_character_, ann$classes), N = n, set = "heldout_test", stage1 = R$best$stage1,
    stage2 = c(as.character(R$best$stage2), rep(NA_character_, m)), stringsAsFactors = FALSE)
  for (metric in c("MAE", "RMSE", "PCC")) {
    Q[[metric]] <- estimates[, metric]
    Q[[paste0(metric, "_low")]] <- boot$ci[, metric, "low"]
    Q[[paste0(metric, "_high")]] <- boot$ci[, metric, "high"]
    Q[[paste0(metric, "_boot_valid")]] <- boot$valid[, metric]
  }
  Q$bootstrap_B <- as.integer(bootstrap_B)
  Q$bootstrap_seed <- as.integer(seed)
  Q$CI_method <- "paired-participant percentile 95%; fixed saved predictions; quantile type 7"

  # Results stay in memory; no graphics, spreadsheets, or supplementary exports
  # are generated by this analysis block.
  list(heldout = Q, contribution = D, contribution_classes = G, class_allocation = L, protein_count_test =
    ia_clock_protein_count_test(D), top_pathway_memberships = ia_clock_top_pathway_memberships(R, Z, D), integration =
    list(beta = A$beta, intercept = A$intercept, center = A$center), recovery = A$recovery, config = list(reactome_release =
    97L, pathway_count = m, unique_input_proteins = length(unique(unlist(R$module_features, use.names = FALSE))),
    protein_pathway_memberships = sum(lengths(R$module_features)), test_n = n, training_n = length(S$train_eid),
    stage1 = R$best$stage1, stage2 = R$best$stage2, pipeline = R$best$pipeline, model_training_seed = R$config$seed,
    bootstrap_seed = seed, bootstrap_B = bootstrap_B, bootstrap_rng = "Mersenne-Twister; Inversion; Rejection",
    resampling_unit = "participant; same multinomial sample across all clocks and three metrics", bootstrap_refitting =
    "none; intervals condition on saved fitted predictions", percentile_interval = "2.5th and 97.5th percentiles; quantile type 7",
    undefined_bootstrap_metric = "CI is NA if any resample is undefined; valid count retained", contribution_definition =
    "phi_ij=beta_j*(X_ij-mean_test(X_j)); I_j=mean_test(abs(phi_ij)); percent=100*I_j/sum(I)", class_allocation =
    "Each pathway contribution is divided equally among its ancestry-supported classes", classification_note =
    attr(ann, "classification_note"), actual_age_used_for_coefficient_recovery = FALSE))
}


# ==============================================================================
# 04 | Clinical associations
# ==============================================================================

# Clinical associations: prevalence, incidence, mortality, and two-year landmarks.
# Age acceleration is defined from pooled full-cohort out-of-fold predictions.
# No figures, supplementary workbooks, input discovery, or result migration.
# All functions below operate on explicitly supplied analysis objects.
# Covariate inputs must be the original tables after upstream median imputation.
# Disease and mortality models, including two-year landmarks, follow the paper.
clinical_need <- function(x, need, label) {
  miss <- setdiff(need, names(x))
  if (length(miss)) stop(label, " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
clinical_dt <- function(x) data.table::as.data.table(data.table::copy(x))
clinical_id <- function(x, label = "EID") {
  y <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(y)) stop("Invalid ", label, ".", call. = FALSE)
  y
}
clinical_ids <- function(x, label) {
  if (anyNA(x$eid) || anyDuplicated(x$eid)) stop(label, " requires one unique non-missing EID per row.", call. =
    FALSE)
}
clinical_days <- function(x) {
  if (inherits(x, "Date")) return(as.integer(x))
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.integer(as.Date(x)))
  if (is.numeric(x)) return(as.integer(as.Date(x, origin = "1970-01-01")))
  as.integer(as.Date(as.character(x)))
}
clinical_covariates <- function(cov, age = NULL) {

  # Prespecified seven-region assessment-centre mapping.
  cov <- as.data.frame(cov, check.names = FALSE)
  ids <- if ("eid" %in% names(cov)) as.character(cov$eid) else rownames(cov)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || identical(ids, as.character(seq_len(nrow(cov))))) stop("Covariates require participant EIDs, not default row numbers.",
    call. = FALSE)
  cov$eid <- NULL
  rownames(cov) <- ids
  if (!"age" %in% names(cov)) {
    if (is.null(age)) stop("Supply baseline_age or covariates containing age.", call. = FALSE)
    age <- as.data.frame(age)
    ai <- if ("eid" %in% names(age)) as.character(age$eid) else rownames(age)
    if (is.null(ai) || anyNA(ai) || anyDuplicated(ai) || !"label" %in% names(age)) stop("Baseline age must contain unique EIDs and label.",
      call. = FALSE)

    # Intersect covariate and baseline-age identifiers before matching OOF predictions.
    ids <- intersect(ids, ai)
    if (!length(ids)) stop("Covariates and baseline age have no shared EIDs.", call. = FALSE)
    cov <- cov[ids, , drop = FALSE]
    cov$age <- age$label[match(ids, ai)]
  }
  COV <- c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
    "bmi")
  clinical_need(cov, COV, "Covariates")
  cov <- cov[, COV, drop = FALSE]
  reg <- list(Scotland = c("11004", "11005"), Wales = c("11003", "11022", "11023"), `NE/Yorkshire` = c("11009",
    "11017", "11027", "11010", "11014"), `North West` = c("11008", "11001", "11016", "10003", "11024", "11025"),
    Midlands = c("11021", "11013", "11006"), South = c("11011", "11028", "11002", "11007", "11026"), London =
    c("11012", "11020", "11018"))
  lut <- setNames(rep(names(reg), lengths(reg)), unlist(reg, use.names = FALSE))
  x <- as.character(cov$assessment_centre)
  x[x == "North East and Yorkshire"] <- "NE/Yorkshire"
  k <- !is.na(x) & x %in% names(lut)
  x[k] <- unname(lut[x[k]])
  bad <- unique(x[!is.na(x) & !x %in% names(reg)])
  if (length(bad)) stop("Unmapped assessment regions: ", paste(bad, collapse = ", "), call. = FALSE)
  cov$assessment_centre <- factor(x, levels = sort(names(reg)))
  for (nm in c("age", "sex", "town_index", "bmi")) {
    z <- suppressWarnings(as.numeric(as.character(cov[[nm]])))
    z[!is.finite(z)] <- NA_real_
    cov[[nm]] <- z
  }
  for (nm in setdiff(COV, c("age", "sex", "town_index", "bmi"))) cov[[nm]] <- factor(cov[[nm]])
  if (any(!stats::na.omit(unique(cov$sex)) %in% c(0, 1))) stop("sex must retain the original 0/1 coding.",
    call. = FALSE)

  # Supply the original median-imputed covariates described in the manuscript.
  # This analysis does not implement imputation; complete-case restrictions below
  # remain safeguards for unresolved missing values or unmatched records.
  cov$eid <- clinical_id(rownames(cov))
  list(cov = clinical_dt(cov), COV = COV)
}
clinical_prepare <- function(result, covariates, blood, outcome_dates, baseline_age, cfg) {
  if (is.null(result$full_cv$oof) || is.null(result$mapping$module)) stop("Clock result requires full_cv$oof and mapping$module; held-out predictions are not valid here.")
  O <- clinical_dt(result$full_cv$oof);
  paths <- as.character(result$mapping$module)
  if (length(paths) != cfg$expected_pathways || anyDuplicated(paths)) stop("Unexpected pathway count/order in clock result.")
  clinical_need(O, c("eid", "label", "fold", "immune_age", paths), "Full OOF predictions")
  P <- c("immune_age", paths)
  if (!is.numeric(O$label) || !all(vapply(O[, ..P], is.numeric, logical(1)))) stop("OOF predictions and age must be numeric.")
  M <- as.matrix(O[, ..P]);
  A <- as.numeric(O$label)
  if (any(!is.finite(A)) || any(!is.finite(M))) stop("OOF predictions contain non-finite values.")
  R <- stats::lm.fit(cbind(1, A), M)$residuals
  Z <- scale(R)
  if (any(!is.finite(Z)) || any(attr(Z, "scaled:scale") <= 0)) stop("Invalid OOF residual acceleration.")
  zn <- sprintf("Z%03d", 0:cfg$expected_pathways)
  E <- data.table::data.table(eid = clinical_id(O$eid), label = A, fold = O$fold)
  ZZ <- data.table::as.data.table(Z);
  data.table::setnames(ZZ, zn);
  E <- cbind(E, ZZ)
  clock_map <- data.table::data.table(exposure = zn, clock = c("Integrated immune age", paths), layer = c("Integrated",
    rep("Pathway", length(paths))))
  age <- baseline_age
  cp <- clinical_covariates(covariates, age)
  CV <- cp$cov;
  COV <- cp$COV
  blood <- clinical_dt(blood)
  clinical_need(blood, c("eid", "blood_date"), "Blood dates")
  B <- blood[, .(eid = clinical_id(eid), blood_day = clinical_days(blood_date))]
  dates <- outcome_dates
  clinical_need(dates, c("hosp_events", "death_data"), "Outcome dates")
  DD <- clinical_dt(dates$death_data);
  clinical_need(DD, "eid", "Death dates")
  if (!any(c("40000-0.0", "40000-1.0", "death_date") %in% names(DD))) stop("Death date columns are absent.")
  d0 <- if ("40000-0.0" %in% names(DD)) clinical_days(DD[["40000-0.0"]]) else rep(NA_integer_, nrow(DD))
  d1 <- if ("40000-1.0" %in% names(DD)) clinical_days(DD[["40000-1.0"]]) else rep(NA_integer_, nrow(DD))
  dm <- ifelse(is.na(d0), d1, ifelse(is.na(d1), d0, pmin(d0, d1)))
  dc <- if ("death_date" %in% names(DD)) clinical_days(DD$death_date) else rep(NA_integer_, nrow(DD))
  death_mismatch <- sum(!(is.na(dm) & is.na(dc)) & (is.na(dm) | is.na(dc) | dm != dc))
  dm[is.na(dm)] <- dc[is.na(dm)]
  D <- data.table::data.table(eid = clinical_id(DD$eid), death_day = as.integer(dm))
  for (nm in c("E", "B", "D", "CV")) clinical_ids(get(nm), nm)

  # Accept either a full death linkage table or death-only records. An absent
  # death record becomes NA after the original left join. Blood/covariate gaps
  # follow the fixed complete-case restriction.
  input_coverage <- data.table::rbindlist(lapply(c("B", "D", "CV"), function(nm) {
    q <- get(nm)
    absent <- length(setdiff(E$eid, q$eid))
    data.table::data.table(input = c(B = "blood", D = "death", CV = "covariates")[[nm]], OOF_participants =
      nrow(E), EIDs_present = nrow(E) - absent, EIDs_absent = absent, handling = if (nm == "D") "Absent death record becomes NA after left join" else "Absent required data follow complete-case exclusion")
  }))
  BASE <- Reduce(function(x, y) merge(x, y, by = "eid", all.x = TRUE, sort = FALSE), list(E, B, D, CV[, c("eid",
    COV), with = FALSE]))
  complete <- stats::complete.cases(BASE[, c("label", "blood_day", zn, COV), with = FALSE])
  baddeath <- !is.na(BASE$death_day) & !is.na(BASE$blood_day) & BASE$death_day <= BASE$blood_day
  audit <- data.table::data.table(stage = c("full_oof", "complete_9_covariates", "invalid_death_on_or_before_blood",
    "analysis_baseline"), N = c(nrow(E), sum(complete), sum(complete & baddeath), sum(complete & !baddeath)))
  BASE <- BASE[complete & !baddeath];
  data.table::setkey(BASE, eid)
  if (!nrow(BASE)) stop("No eligible participants after original preparation.")
  DS <- clinical_dt(dates$hosp_events);
  clinical_need(DS, c("eid", "icd3", "event_date"), "Hospital events")
  DX <- DS[eid %in% BASE$eid, .(eid = clinical_id(eid), icd3 = toupper(trimws(as.character(icd3))), dx_day =
    clinical_days(event_date))]
  DX <- DX[!is.na(dx_day) & grepl("^[A-Z][0-9]{2}$", icd3) & substr(icd3, 1, 1) %in% LETTERS[1:14]]
  FIRST <- DX[, .(dx_day = min(dx_day)), by = .(icd3, eid)]
  data.table::setkey(FIRST, icd3, eid)
  age_label_audit <- BASE[, .(N = .N, mean_difference = mean(age - label), median_absolute_difference = median(abs(age -
    label)), max_absolute_difference = max(abs(age - label)))]
  acceleration_check <- if ("age_acceleration" %in% names(O) && all(is.finite(O$age_acceleration))) max(abs(R[,
    1] - O$age_acceleration)) else NA_real_
  list(BASE = BASE, FIRST = FIRST, clock_map = clock_map, COV = COV, age_label_audit = age_label_audit, acceleration_check =
    acceleration_check, preparation_audit = audit, input_coverage = input_coverage, death_date_mismatches =
    death_mismatch, acceleration = "Pooled full-OOF predicted age residualized on chronological age; sample SD standardization before outcome cohort restriction")
}
clinical_scope <- function(cd) {
  if (grepl("^(C6[0-3]|D29|D40|N(4[0-9]|5[0-3]))$", cd)) "Male" else if (grepl("^(C5[1-8]|D06|D2[5-8]|D39|N(7[0-9]|8[0-9]|9[0-8]))$",
    cd)) "Female" else "All"
}
clinical_warntext <- function(w) paste(unique(w[!is.na(w) & nzchar(w)]), collapse = " | ")
clinical_status <- function(w) {
  z <- clinical_warntext(w)
  if (!nzchar(z)) "ok" else if (grepl("did not converge|ran out of iterations|coefficients? may be infinite|fitted probabilities numerically 0 or 1|separation|boundary",
    z, ignore.case = TRUE)) "warning_review" else "warning"
}
clinical_empty <- function(q = NULL, status = "error", warning = "", error = "") {
  data.table::data.table(N = if (is.null(q)) NA_integer_ else length(q$y), Events = if (is.null(q)) NA_integer_ else sum(q$y),
    PersonYears = if (is.null(q)) NA_real_ else sum(q$time), Competing = if (is.null(q)) NA_integer_ else sum(q$cr ==
    2L), Beta = NA_real_, SE = NA_real_, Effect = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_, P = NA_real_,
    status = status, warning = warning, error = error, aliased = "")
}
clinical_warning_indices <- function(w) {
  hits <- regmatches(w, gregexpr("variables?\\s+[0-9 ,]+", w, ignore.case = TRUE))
  text <- unlist(hits, use.names = FALSE)
  if (!length(text)) return(integer())
  unique(as.integer(unlist(regmatches(text, gregexpr("[0-9]+", text)))))
}
clinical_design_audit <- function(fit, data, warnings) {
  mm <- fit$x
  co <- stats::coef(fit);
  se <- sqrt(diag(fit$var));
  names(se) <- names(co)
  ij <- clinical_warning_indices(warnings)
  if (is.null(mm)) return(list(columns = data.table::data.table(), categories = data.table::data.table()))
  columns <- data.table::rbindlist(lapply(seq_len(ncol(mm)), function(j) {
    v <- mm[, j];
    binary <- all(unique(v) %in% c(0, 1));
    nm <- colnames(mm)[j]
    data.table::data.table(design_index = j, design_column = nm, warning_variable_index = j %in% ij, coefficient =
      unname(co[nm]), standard_error = unname(se[nm]), aliased = is.na(co[nm]), n_unique = length(unique(v)),
      mean = mean(v), sd = stats::sd(v), minimum = min(v), maximum = max(v), binary_indicator = binary, N_value1 =
      if (binary) sum(v == 1) else NA_integer_, Events_value1 = if (binary) sum(data$outcome[v == 1]) else NA_integer_,
      N_value0 = if (binary) sum(v == 0) else NA_integer_, Events_value0 = if (binary) sum(data$outcome[v ==
      0]) else NA_integer_)
  }), fill = TRUE)
  categorical <- names(data)[vapply(data, is.factor, logical(1))]
  if ("sex" %in% names(data)) categorical <- unique(c(categorical, "sex"))
  categories <- data.table::rbindlist(lapply(categorical, function(nm) {
    z <- data.table::data.table(level = as.character(data[[nm]]), event = data$outcome)
    z[, .(covariate = nm, N = .N, Events = sum(event), Non_events = sum(event == 0)), by = level]
  }), fill = TRUE)
  list(columns = columns, categories = categories)
}
clinical_fit_incident <- function(q, xn, COV, diagnostic = FALSE) {
  wf <- character()
  design <- list(columns = data.table::data.table(), categories = data.table::data.table())
  result <- tryCatch({
    z <- data.frame(time = q$time, outcome = as.integer(q$y), exposure = q$d[[xn]], q$d[, ..COV], check.names =
      FALSE)
    f <- withCallingHandlers(survival::coxph(stats::reformulate(c("exposure", COV), "survival::Surv(time,outcome)"),
      data = z, ties = "efron", x = diagnostic, model = FALSE, na.action = stats::na.fail, singular.ok =
      TRUE), warning = function(w) {
      wf <<- c(wf, conditionMessage(w));
      invokeRestart("muffleWarning")
    })
    s <- summary(f)$coefficients
    if (!"exposure" %in% rownames(s) || any(!is.finite(s["exposure", c("coef", "se(coef)")]))) stop("Exposure coefficient is not estimable")
    b <- s["exposure", "coef"];
    se <- s["exposure", "se(coef)"]
    if (!is.finite(se) || se <= 0) stop("Exposure standard error is not positive")
    if (diagnostic) design <- clinical_design_audit(f, z, wf)
    data.table::data.table(N = length(q$y), Events = sum(q$y), PersonYears = sum(q$time), Competing = sum(q$cr ==
      2L), Beta = b, SE = se, Effect = exp(b), Lower95 = exp(b - 1.96 * se), Upper95 = exp(b + 1.96 * se),
      P = s["exposure", "Pr(>|z|)"], status = clinical_status(wf), warning = clinical_warntext(wf), error =
      "", aliased = paste(names(stats::coef(f))[is.na(stats::coef(f))], collapse = ";"))
  }, error = function(e) clinical_empty(q, warning = clinical_warntext(wf), error = conditionMessage(e)))
  list(result = result, design = design)
}
clinical_estimate_usable <- function(z) {
  need <- c("Beta", "SE", "Effect", "Lower95", "Upper95", "P")
  clinical_need(z, need, "Model table")
  is.finite(z$Beta) & is.finite(z$SE) & z$SE > 0 & is.finite(z$Effect) & z$Effect > 0 & is.finite(z$Lower95) &
  z$Lower95 > 0 & is.finite(z$Upper95) & z$Upper95 >= z$Lower95 & is.finite(z$P) & z$P >= 0 & z$P <= 1
}
clinical_parallel <- function(X, FUN, cores) {
  if (!length(X)) return(list())
  if (.Platform$OS.type != "windows" && cores > 1L) {
    ans <- parallel::mclapply(X, FUN, mc.cores = min(cores, length(X)), mc.preschedule = FALSE, mc.set.seed =
      FALSE)
  } else ans <- lapply(X, FUN)
  ans
}
clinical_model_usable <- function(d, ci = FALSE) {
  need <- c("Beta", "Effect", "P_Bonf", "status")
  if (length(setdiff(need, names(d)))) stop("Missing model columns: ", paste(setdiff(need, names(d)), collapse =
    ", "))
  ok <- as.character(d$status) %in% c("ok", "warning", "warning_review") & is.finite(d$Beta) & is.finite(d$Effect) &
  d$Effect > 0 & is.finite(d$P_Bonf) & d$P_Bonf >= 0 & d$P_Bonf <= 1
  if ("error" %in% names(d)) ok <- ok & (is.na(d$error) | trimws(as.character(d$error)) %in% c("", "NA", "<NA>"))
  if ("SE" %in% names(d)) ok <- ok & is.finite(d$SE) & d$SE > 0
  if ("P" %in% names(d)) ok <- ok & is.finite(d$P) & d$P >= 0 & d$P <= 1
  if (ci) ok <- ok & is.finite(d$Lower95) & is.finite(d$Upper95) & d$Lower95 > 0 & d$Lower95 <= d$Effect &
  d$Upper95 >= d$Effect
  ok[is.na(ok)] <- FALSE;
  ok
}

# Reconstruct the complete clinical analysis from participant-level inputs.
# Planned incident multiplicity remains 274 endpoints, including pathway
# follow-up tested against the full 58 x 274 family rather than the gated subset.
# Finite estimates with review warnings are retained for incident pathway
# eligibility, with their warnings preserved. Prevalent pathway follow-up
# retains its original, narrower warning-status gate.
run_clinical_analysis <- function(result, covariates, blood, outcome_dates, baseline_age = NULL, config =
  list()) {
  required <- c("data.table", "survival")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  cfg <- utils::modifyList(list(cores = 30L, minimum_events = 200L, alpha = 0.05, disease_censor = "2022-05-31",
    death_censor = "2024-08-31", expected_pathways = 58L, expected_incident = 274L), config)
  cfg$cores <- max(1L, as.integer(cfg$cores))
  cfg$expected_pathways <- as.integer(cfg$expected_pathways)
  cfg$expected_incident <- as.integer(cfg$expected_incident)
  cfg$minimum_events <- as.integer(cfg$minimum_events)
  if (cfg$alpha <= 0 || cfg$alpha >= 1 || cfg$expected_pathways < 1L || cfg$expected_incident < 1L) stop("Invalid clinical analysis configuration.")
  NCORE <- cfg$cores;
  MIN_CASES <- cfg$minimum_events;
  ALPHA <- cfg$alpha
  DIS_CENSOR <- as.Date(cfg$disease_censor);
  DEATH_CENSOR <- as.Date(cfg$death_censor)
  if (anyNA(c(DIS_CENSOR, DEATH_CENSOR))) stop("Invalid administrative censor date.")
  INC_CENSOR <- min(DIS_CENSOR, DEATH_CENSOR)
  MALE <- 1;
  FEMALE <- 0;
  KEEP_CHAPTERS <- LETTERS[1:14]
  old_options <- options(datatable.aware = TRUE)
  old_threads <- data.table::getDTthreads()
  thread_variables <- c("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  old_environment <- Sys.getenv(thread_variables, unset = NA_character_)
  on.exit({
    options(old_options);
    data.table::setDTthreads(old_threads)
    for (i in seq_along(thread_variables)) {
      if (is.na(old_environment[i])) Sys.unsetenv(thread_variables[i]) else do.call(Sys.setenv, setNames(list(old_environment[i]),
        thread_variables[i]))
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_variables)), thread_variables)))
  data.table::setDTthreads(1L)
  `%chin%` <- data.table::`%chin%`
  prep <- clinical_prepare(result, covariates, blood, outcome_dates, baseline_age, cfg)
  BASE <- prep$BASE;
  FIRST <- prep$FIRST;
  clock_map <- prep$clock_map;
  COV <- prep$COV
  zn <- as.character(clock_map$exposure)
  BDIS <- BASE[blood_day <= as.integer(DIS_CENSOR)]
  MC <- function(X, FUN) clinical_parallel(X, FUN, NCORE)
  scope <- function(x) if (grepl("^(C6[0-3]|D29|D40|N(4[0-9]|5[0-3]))$", x)) "Male" else if (grepl("^(C5[1-8]|D06|D2[5-8]|D39|N(7[0-9]|8[0-9]|9[0-8]))$",
    x)) "Female" else "All"
  riskbase <- function(cd) {
    s <- scope(cd);
    if (s == "Male") BDIS[sex == MALE] else if (s == "Female") BDIS[sex == FEMALE] else BDIS
  }
  dxvec <- function(cd, eid) {
    q <- FIRST[icd3 == cd];
    q$dx_day[match(eid, q$eid)]
  }
  make_prev <- function(cd) {
    d <- riskbase(cd);
    x <- dxvec(cd, d$eid);
    list(d = d, y = !is.na(x) & x <= d$blood_day)
  }
  make_inc <- function(cd, lag = 0) {
    d <- riskbase(cd);
    x <- dxvec(cd, d$eid);
    k <- (is.na(x) | x > d$blood_day) & d$blood_day < as.integer(INC_CENSOR);
    d <- d[k];
    x <- x[k];
    b <- d$blood_day;
    q <- d$death_day;
    cc <- as.integer(INC_CENSOR)
    ev <- !is.na(x) & x > b & x <= cc & (is.na(q) | x <= q);
    exit <- ifelse(ev, x, ifelse(!is.na(q) & q > b & q <= cc, q, cc));
    cr <- ifelse(ev, 1L, ifelse(!is.na(q) & q > b & q <= cc, 2L, 0L));
    t0 <- (exit - b) / 365.25;
    k <- is.finite(t0) & t0 > lag;
    list(d = d[k], time = t0[k] - lag, y = ev[k], cr = cr[k])
  }
  make_mort <- function(lag = 0) {
    d <- BASE[blood_day < as.integer(DEATH_CENSOR)];
    b <- d$blood_day;
    q <- d$death_day;
    cc <- as.integer(DEATH_CENSOR);
    ev <- !is.na(q) & q > b & q <= cc;
    exit <- ifelse(ev, q, cc);
    t0 <- (exit - b) / 365.25;
    k <- is.finite(t0) & t0 > lag;
    list(d = d[k], time = t0[k] - lag, y = ev[k])
  }
  emptyfit <- function(n = NA_integer_, e = NA_integer_, py = NA_real_, comp = NA_integer_, status = "error",
    warning = "", error = "") data.table::data.table(N = n, Events = e, PersonYears = py, Competing = comp,
    Beta = NA_real_, SE = NA_real_, Effect = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_, P = NA_real_,
    status = status, warning = warning, error = error, aliased = "")
  warntext <- function(w) paste(unique(w[!is.na(w) & nzchar(w)]), collapse = " | ")
  fitstatus <- function(w) {
    z <- warntext(w);
    if (!nzchar(z)) "ok" else if (grepl("did not converge|ran out of iterations|coefficients? may be infinite|fitted probabilities numerically 0 or 1|separation|boundary",
      z, ignore.case = TRUE)) "warning_review" else "warning"
  }
  fit_logit <- function(d, y, xn) {
    n <- length(y);
    e <- sum(y);
    w <- character();
    tryCatch({
      q <- data.frame(outcome = as.integer(y), exposure = d[[xn]], d[, ..COV], check.names = FALSE);
      f <- withCallingHandlers(glm(reformulate(c("exposure", COV), "outcome"), data = q, family = binomial(),
        na.action = na.fail, control = glm.control(maxit = 50)), warning = function(a) {
        w <<- c(w, conditionMessage(a));
        invokeRestart("muffleWarning")
      });
      if (!isTRUE(f$converged)) w <- c(w, "glm did not converge");
      s <- coef(summary(f));
      if (!"exposure" %in% rownames(s) || any(!is.finite(s["exposure", 1:2]))) stop("Exposure coefficient is not estimable");
      b <- s["exposure", 1];
      se <- s["exposure", 2];
      data.table::data.table(N = n, Events = e, PersonYears = NA_real_, Competing = NA_integer_, Beta = b,
        SE = se, Effect = exp(b), Lower95 = exp(b - 1.96 * se), Upper95 = exp(b + 1.96 * se), P = s["exposure",
        4], status = fitstatus(w), warning = warntext(w), error = "", aliased = paste(names(coef(f))[is.na(coef(f))],
        collapse = ";"))
    }, error = function(a) emptyfit(n, e, status = "error", warning = warntext(w), error = conditionMessage(a)))
  }
  fit_cox <- function(d, time, y, xn) {
    n <- length(y);
    e <- sum(y);
    py <- sum(time);
    wf <- character()
    tryCatch({
      q <- data.frame(time = time, outcome = as.integer(y), exposure = d[[xn]], d[, ..COV], check.names =
        FALSE)
      f <- withCallingHandlers(survival::coxph(reformulate(c("exposure", COV), "survival::Surv(time,outcome)"),
        data = q, ties = "efron", x = FALSE, model = FALSE, na.action = na.fail, singular.ok = TRUE), warning =
        function(a) {
        wf <<- c(wf, conditionMessage(a));
        invokeRestart("muffleWarning")
      })
      s <- summary(f)$coefficients
      if (!"exposure" %in% rownames(s) || any(!is.finite(s["exposure", c("coef", "se(coef)")]))) stop("Exposure coefficient is not estimable")
      b <- s["exposure", "coef"];
      se <- s["exposure", "se(coef)"]
      data.table::data.table(N = n, Events = e, PersonYears = py, Competing = NA_integer_, Beta = b, SE =
        se, Effect = exp(b), Lower95 = exp(b - 1.96 * se), Upper95 = exp(b + 1.96 * se), P = s["exposure",
        "Pr(>|z|)"], status = fitstatus(wf), warning = warntext(wf), error = "", aliased = paste(names(coef(f))[is.na(coef(f))],
        collapse = ";"))
    }, error = function(a) emptyfit(n, e, py, status = "error", warning = warntext(wf), error = conditionMessage(a)))
  }
  failmany <- function(cd, xs, msg, status = "error") cbind(data.table::data.table(icd3 = cd, exposure =
    xs), emptyfit(status = status, error = msg)[rep(1L, length(xs))])
  annotate <- function(z) if (!nrow(z)) z else merge(z, clock_map, by = "exposure", all.x = TRUE, sort =
    FALSE)
  adjust_bonf <- function(z, fam, mplan = nrow(z)) {
    if (!nrow(z)) return(z)
    if (length(mplan) != 1L || is.na(mplan) || mplan < nrow(z)) stop("Invalid Bonferroni family size: ", fam)
    adj <- function(v) {
      a <- rep(NA_real_, length(v));
      ok <- is.finite(v)
      a[ok] <- p.adjust(v[ok], "bonferroni", n = mplan);
      a
    }
    z[, `:=`(family = fam, m_tests = as.integer(mplan), P_Bonf = adj(P))]
    z[]
  }
  run_prev <- function(cc, xs, fam, mplan) {
    if (!length(cc)) return(data.table::data.table());
    z <- data.table::rbindlist(MC(cc, function(cd) tryCatch({
      q <- make_prev(cd);
      data.table::rbindlist(lapply(xs, function(xn) cbind(data.table::data.table(icd3 = cd, exposure = xn),
        fit_logit(q$d, q$y, xn))))
    }, error = function(e) failmany(cd, xs, conditionMessage(e)))), fill = TRUE);
    z[, `:=`(analysis = "Prevalent", measure = "OR", landmark_years = 0)];
    adjust_bonf(annotate(z), fam, mplan)
  }
  run_inc <- function(cc, xs, fam, mplan, lag = 0) {
    if (!length(cc)) return(data.table::data.table());
    z <- data.table::rbindlist(MC(cc, function(cd) tryCatch({
      q <- make_inc(cd, lag);
      if (sum(q$y) < MIN_CASES) return(failmany(cd, xs, paste0("Events < ", MIN_CASES), "insufficient_events"));
      data.table::rbindlist(lapply(xs, function(xn) {
        a <- clinical_fit_incident(q, xn, COV)$result;
        a[, Competing := sum(q$cr == 2L)];
        cbind(data.table::data.table(icd3 = cd, exposure = xn), a)
      }))
    }, error = function(e) failmany(cd, xs, conditionMessage(e)))), fill = TRUE);
    z[, `:=`(analysis = "Incident", measure = "HR", landmark_years = lag)];
    adjust_bonf(annotate(z), fam, mplan)
  }
  run_mort <- function(xs, fam, mplan, lag = 0) {
    q <- make_mort(lag);
    if (sum(q$y) < MIN_CASES) {
      z <- failmany("All-cause mortality", xs, paste0("Deaths < ", MIN_CASES), "insufficient_events")
    } else z <- data.table::rbindlist(MC(xs, function(xn) cbind(data.table::data.table(icd3 = "All-cause mortality",
      exposure = xn), fit_cox(q$d, q$time, q$y, xn))), fill = TRUE);
    z[, `:=`(analysis = "Mortality", measure = "HR", landmark_years = lag)];
    adjust_bonf(annotate(z), fam, mplan)
  }
  hitcodes <- function(z) if (nrow(z) && all(c("P_Bonf", "status", "icd3") %in% names(z))) unique(z[is.finite(P_Bonf) &
    P_Bonf < ALPHA & status %chin% c("ok", "warning"), icd3]) else character()
  codes <- sort(unique(FIRST$icd3));
  counts <- data.table::rbindlist(MC(codes, function(cd) tryCatch({
    p <- make_prev(cd);
    i <- make_inc(cd);
    data.table::data.table(icd3 = cd, scope = scope(cd), baseline_n = nrow(p$d), prevalent_cases = sum(p$y),
      prevalent_controls = sum(!p$y), incident_n = nrow(i$d), incident_events = sum(i$y), competing_deaths =
      sum(i$cr == 2L))
  }, error = function(e) data.table::data.table(icd3 = cd, scope = scope(cd), count_error = conditionMessage(e)))),
    fill = TRUE)
  if (!length(codes)) counts <- data.table::data.table(icd3 = character(), scope = character(), baseline_n =
    integer(), prevalent_cases = integer(), prevalent_controls = integer(), incident_n = integer(), incident_events =
    integer(), competing_deaths = integer());
  if (!"count_error" %in% names(counts)) counts[, count_error := NA_character_]
  counts[, `:=`(eligible_prevalent = is.na(count_error) & prevalent_cases >= MIN_CASES & prevalent_controls >=
    MIN_CASES, eligible_incident = is.na(count_error) & incident_events >= MIN_CASES)]
  prev_codes <- counts[eligible_prevalent == TRUE, icd3];
  inc_codes <- counts[eligible_incident == TRUE, icd3]

  # Integrated associations and the original two-year landmark sensitivities.
  message("Clinical analysis: ", length(prev_codes), " eligible prevalent endpoints; ", length(inc_codes),
    " eligible incident endpoints.")
  prev_int <- run_prev(prev_codes, zn[1], "prevalent_integrated", length(prev_codes))
  inc_int <- run_inc(inc_codes, zn[1], "incident_integrated", cfg$expected_incident)
  mort_int <- run_mort(zn[1], "mortality_integrated", 1L)
  inc_lag2 <- run_inc(inc_codes, zn[1], "incident_integrated_landmark2", cfg$expected_incident, lag = 2)
  mort_lag2 <- run_mort(zn[1], "mortality_integrated_landmark2", 1L, lag = 2)
  mort_path <- run_mort(zn[ - 1], "mortality_pathway", cfg$expected_pathways)

  # Preserve the original prevalent pathway gate. The repaired incident
  # pathway gate also retains finite estimates requiring warning review.
  prev_hit <- hitcodes(prev_int)
  inc_hit <- if (nrow(inc_int)) unique(inc_int[clinical_estimate_usable(inc_int) & is.finite(P_Bonf) & P_Bonf <
    ALPHA, icd3]) else character()
  selected <- data.table::rbindlist(list(data.table::data.table(analysis = rep("Prevalent", length(prev_hit)),
    icd3 = prev_hit), data.table::data.table(analysis = rep("Incident", length(inc_hit)), icd3 = inc_hit)))
  prev_path <- run_prev(prev_hit, zn[ - 1], "prevalent_pathway_gated", cfg$expected_pathways * length(prev_codes))
  inc_path <- run_inc(inc_hit, zn[ - 1], "incident_pathway_gated", cfg$expected_pathways * cfg$expected_incident)

  # Retain association estimates and model convergence/error records.
  sets <- list(prev_int = prev_int, inc_int = inc_int, mort_int = mort_int, inc_lag2 = inc_lag2, mort_lag2 =
    mort_lag2, mort_path = mort_path, prev_path = prev_path, inc_path = inc_path)
  all_models <- data.table::rbindlist(Map(function(z, n) {
    if (!nrow(z)) return(NULL)
    z <- data.table::copy(z);
    z[, result_set := n];
    z
  }, sets, names(sets)), fill = TRUE)
  warnings_table <- all_models[!is.na(warning) & nzchar(warning)]
  errors_table <- all_models[status %chin% c("error", "insufficient_events", "package_missing") | (!is.na(error) &
    nzchar(error))]
  analysis_config <- c(cfg, list(incident_censor = INC_CENSOR, chapters = KEEP_CHAPTERS, covariates = COV,
    incident_pathway_family_size = cfg$expected_pathways * cfg$expected_incident, incident_pathway_gate =
    "Finite integrated Bonferroni-significant estimates, including retained review warnings", other_gates =
    "Prevalent pathways retain the original ok/warning status gate", death_date_mismatches = prep$death_date_mismatches,
    existing_acceleration_max_abs_difference = prep$acceleration_check))
  list(config = analysis_config, audit = prep$preparation_audit, age_label_audit = prep$age_label_audit, clock_map =
    clock_map, endpoint_counts = counts, selected_endpoints = selected, results = sets, warnings = warnings_table,
    errors = errors_table, prepared = prep, session = sessionInfo())
}


# ==============================================================================
# 05 | Disease association summaries
# ==============================================================================

# Descriptive summaries of prevalent and incident disease associations.
# Models and saved multiplicity corrections are supplied by the clinical modules.
# Public ICD metadata are embedded; no participant results are embedded.
# Dependencies: data.table, classify_reactome8, reactome8_classes, and the
# mortality_profile_eligibility / mortality_disease_labels metadata helpers.
disease_icd_dictionary <- function() {

  # Frozen public UK Biobank Data-Coding 19, restricted to three-character A-N codes.
  # Columns retain the official disease and chapter names from the supplied source.
  d <- utils::read.delim(text = "code	disease_name	chapter	chapter_name
A00	Cholera	Chapter I	Certain infectious and parasitic diseases
A01	Typhoid and paratyphoid fevers	Chapter I	Certain infectious and parasitic diseases
A02	Other Salmonella infections	Chapter I	Certain infectious and parasitic diseases
A03	Shigellosis	Chapter I	Certain infectious and parasitic diseases
A04	Other bacterial intestinal infections	Chapter I	Certain infectious and parasitic diseases
A05	Other bacterial foodborne intoxications	Chapter I	Certain infectious and parasitic diseases
A06	Amoebiasis	Chapter I	Certain infectious and parasitic diseases
A07	Other protozoal intestinal diseases	Chapter I	Certain infectious and parasitic diseases
A08	Viral and other specified intestinal infections	Chapter I	Certain infectious and parasitic diseases
A09	Diarrhoea and gastro-enteritis of presumed infectious origin	Chapter I	Certain infectious and parasitic diseases
A15	Respiratory tuberculosis, bacteriologically and histologically confirmed	Chapter I	Certain infectious and parasitic diseases
A16	Respiratory tuberculosis, not confirmed bacteriologically or histologically	Chapter I	Certain infectious and parasitic diseases
A17	Tuberculosis of nervous system	Chapter I	Certain infectious and parasitic diseases
A18	Tuberculosis of other organs	Chapter I	Certain infectious and parasitic diseases
A19	Miliary tuberculosis	Chapter I	Certain infectious and parasitic diseases
A20	Plague	Chapter I	Certain infectious and parasitic diseases
A21	Tularaemia	Chapter I	Certain infectious and parasitic diseases
A22	Anthrax	Chapter I	Certain infectious and parasitic diseases
A23	Brucellosis	Chapter I	Certain infectious and parasitic diseases
A24	Glanders and melioidosis	Chapter I	Certain infectious and parasitic diseases
A25	Rat-bite fevers	Chapter I	Certain infectious and parasitic diseases
A26	Erysipeloid	Chapter I	Certain infectious and parasitic diseases
A27	Leptospirosis	Chapter I	Certain infectious and parasitic diseases
A28	Other zoonotic bacterial diseases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A30	Leprosy [Hansen's disease]	Chapter I	Certain infectious and parasitic diseases
A31	Infection due to other mycobacteria	Chapter I	Certain infectious and parasitic diseases
A32	Listeriosis	Chapter I	Certain infectious and parasitic diseases
A33	Tetanus neonatorum	Chapter I	Certain infectious and parasitic diseases
A34	Obstetrical tetanus	Chapter I	Certain infectious and parasitic diseases
A35	Other tetanus	Chapter I	Certain infectious and parasitic diseases
A36	Diphtheria	Chapter I	Certain infectious and parasitic diseases
A37	Whooping cough	Chapter I	Certain infectious and parasitic diseases
A38	Scarlet fever	Chapter I	Certain infectious and parasitic diseases
A39	Meningococcal infection	Chapter I	Certain infectious and parasitic diseases
A40	Streptococcal septicaemia	Chapter I	Certain infectious and parasitic diseases
A41	Other septicaemia	Chapter I	Certain infectious and parasitic diseases
A42	Actinomycosis	Chapter I	Certain infectious and parasitic diseases
A43	Nocardiosis	Chapter I	Certain infectious and parasitic diseases
A44	Bartonellosis	Chapter I	Certain infectious and parasitic diseases
A46	Erysipelas	Chapter I	Certain infectious and parasitic diseases
A48	Other bacterial diseases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A49	Bacterial infection of unspecified site	Chapter I	Certain infectious and parasitic diseases
A50	Congenital syphilis	Chapter I	Certain infectious and parasitic diseases
A51	Early syphilis	Chapter I	Certain infectious and parasitic diseases
A52	Late syphilis	Chapter I	Certain infectious and parasitic diseases
A53	Other and unspecified syphilis	Chapter I	Certain infectious and parasitic diseases
A54	Gonococcal infection	Chapter I	Certain infectious and parasitic diseases
A55	Chlamydial lymphogranuloma (venereum)	Chapter I	Certain infectious and parasitic diseases
A56	Other sexually transmitted chlamydial diseases	Chapter I	Certain infectious and parasitic diseases
A57	Chancroid	Chapter I	Certain infectious and parasitic diseases
A58	Granuloma inguinale	Chapter I	Certain infectious and parasitic diseases
A59	Trichomoniasis	Chapter I	Certain infectious and parasitic diseases
A60	Anogenital herpesviral [herpes simplex] infections	Chapter I	Certain infectious and parasitic diseases
A63	Other predominantly sexually transmitted diseases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A64	Unspecified sexually transmitted disease	Chapter I	Certain infectious and parasitic diseases
A65	Non-venereal syphilis	Chapter I	Certain infectious and parasitic diseases
A66	Yaws	Chapter I	Certain infectious and parasitic diseases
A67	Pinta [carate]	Chapter I	Certain infectious and parasitic diseases
A68	Relapsing fevers	Chapter I	Certain infectious and parasitic diseases
A69	Other spirochaetal infections	Chapter I	Certain infectious and parasitic diseases
A70	Chlamydia psittaci infection	Chapter I	Certain infectious and parasitic diseases
A71	Trachoma	Chapter I	Certain infectious and parasitic diseases
A74	Other diseases caused by chlamydiae	Chapter I	Certain infectious and parasitic diseases
A75	Typhus fever	Chapter I	Certain infectious and parasitic diseases
A77	Spotted fever [tick-borne rickettsioses]	Chapter I	Certain infectious and parasitic diseases
A78	Q fever	Chapter I	Certain infectious and parasitic diseases
A79	Other rickettsioses	Chapter I	Certain infectious and parasitic diseases
A80	Acute poliomyelitis	Chapter I	Certain infectious and parasitic diseases
A81	Atypical virus infections of central nervous system	Chapter I	Certain infectious and parasitic diseases
A82	Rabies	Chapter I	Certain infectious and parasitic diseases
A83	Mosquito-borne viral encephalitis	Chapter I	Certain infectious and parasitic diseases
A84	Tick-borne viral encephalitis	Chapter I	Certain infectious and parasitic diseases
A85	Other viral encephalitis, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A86	Unspecified viral encephalitis	Chapter I	Certain infectious and parasitic diseases
A87	Viral meningitis	Chapter I	Certain infectious and parasitic diseases
A88	Other viral infections of central nervous system, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A89	Unspecified viral infection of central nervous system	Chapter I	Certain infectious and parasitic diseases
A90	Dengue fever [classical dengue]	Chapter I	Certain infectious and parasitic diseases
A91	Dengue haemorrhagic fever	Chapter I	Certain infectious and parasitic diseases
A92	Other mosquito-borne viral fevers	Chapter I	Certain infectious and parasitic diseases
A93	Other arthropod-borne viral fevers, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A94	Unspecified arthropod-borne viral fever	Chapter I	Certain infectious and parasitic diseases
A95	Yellow fever	Chapter I	Certain infectious and parasitic diseases
A96	Arenaviral haemorrhagic fever	Chapter I	Certain infectious and parasitic diseases
A97	Dengue	Chapter I	Certain infectious and parasitic diseases
A98	Other viral haemorrhagic fevers, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
A99	Unspecified viral haemorrhagic fever	Chapter I	Certain infectious and parasitic diseases
B00	Herpesviral [herpes simplex] infections	Chapter I	Certain infectious and parasitic diseases
B01	Varicella [chickenpox]	Chapter I	Certain infectious and parasitic diseases
B02	Zoster [herpes zoster]	Chapter I	Certain infectious and parasitic diseases
B03	Smallpox	Chapter I	Certain infectious and parasitic diseases
B04	Monkeypox	Chapter I	Certain infectious and parasitic diseases
B05	Measles	Chapter I	Certain infectious and parasitic diseases
B06	Rubella [German measles]	Chapter I	Certain infectious and parasitic diseases
B07	Viral warts	Chapter I	Certain infectious and parasitic diseases
B08	Other viral infections characterised by skin and mucous membrane lesions, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
B09	Unspecified viral infection characterised by skin and mucous membrane lesions	Chapter I	Certain infectious and parasitic diseases
B15	Acute hepatitis A	Chapter I	Certain infectious and parasitic diseases
B16	Acute hepatitis B	Chapter I	Certain infectious and parasitic diseases
B17	Other acute viral hepatitis	Chapter I	Certain infectious and parasitic diseases
B18	Chronic viral hepatitis	Chapter I	Certain infectious and parasitic diseases
B19	Unspecified viral hepatitis	Chapter I	Certain infectious and parasitic diseases
B20	Human immunodeficiency virus [HIV] disease resulting in infectious and parasitic diseases	Chapter I	Certain infectious and parasitic diseases
B21	Human immunodeficiency virus [HIV] disease resulting in malignant neoplasms	Chapter I	Certain infectious and parasitic diseases
B22	Human immunodeficiency virus [HIV] disease resulting in other specified diseases	Chapter I	Certain infectious and parasitic diseases
B23	Human immunodeficiency virus [HIV] disease resulting in other conditions	Chapter I	Certain infectious and parasitic diseases
B24	Unspecified human immunodeficiency virus [HIV] disease	Chapter I	Certain infectious and parasitic diseases
B25	Cytomegaloviral disease	Chapter I	Certain infectious and parasitic diseases
B26	Mumps	Chapter I	Certain infectious and parasitic diseases
B27	Infectious mononucleosis	Chapter I	Certain infectious and parasitic diseases
B30	Viral conjunctivitis	Chapter I	Certain infectious and parasitic diseases
B33	Other viral diseases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
B34	Viral infection of unspecified site	Chapter I	Certain infectious and parasitic diseases
B35	Dermatophytosis	Chapter I	Certain infectious and parasitic diseases
B36	Other superficial mycoses	Chapter I	Certain infectious and parasitic diseases
B37	Candidiasis	Chapter I	Certain infectious and parasitic diseases
B38	Coccidioidomycosis	Chapter I	Certain infectious and parasitic diseases
B39	Histoplasmosis	Chapter I	Certain infectious and parasitic diseases
B40	Blastomycosis	Chapter I	Certain infectious and parasitic diseases
B41	Paracoccidioidomycosis	Chapter I	Certain infectious and parasitic diseases
B42	Sporotrichosis	Chapter I	Certain infectious and parasitic diseases
B43	Chromomycosis and phaeomycotic abscess	Chapter I	Certain infectious and parasitic diseases
B44	Aspergillosis	Chapter I	Certain infectious and parasitic diseases
B45	Cryptococcosis	Chapter I	Certain infectious and parasitic diseases
B46	Zygomycosis	Chapter I	Certain infectious and parasitic diseases
B47	Mycetoma	Chapter I	Certain infectious and parasitic diseases
B48	Other mycoses, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
B49	Unspecified mycosis	Chapter I	Certain infectious and parasitic diseases
B50	Plasmodium falciparum malaria	Chapter I	Certain infectious and parasitic diseases
B51	Plasmodium vivax malaria	Chapter I	Certain infectious and parasitic diseases
B52	Plasmodium malariae malaria	Chapter I	Certain infectious and parasitic diseases
B53	Other parasitologically confirmed malaria	Chapter I	Certain infectious and parasitic diseases
B54	Unspecified malaria	Chapter I	Certain infectious and parasitic diseases
B55	Leishmaniasis	Chapter I	Certain infectious and parasitic diseases
B56	African trypanosomiasis	Chapter I	Certain infectious and parasitic diseases
B57	Chagas' disease	Chapter I	Certain infectious and parasitic diseases
B58	Toxoplasmosis	Chapter I	Certain infectious and parasitic diseases
B59	Pneumocystosis	Chapter I	Certain infectious and parasitic diseases
B60	Other protozoal diseases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
B64	Unspecified protozoal disease	Chapter I	Certain infectious and parasitic diseases
B65	Schistosomiasis [bilharziasis]	Chapter I	Certain infectious and parasitic diseases
B66	Other fluke infections	Chapter I	Certain infectious and parasitic diseases
B67	Echinococcosis	Chapter I	Certain infectious and parasitic diseases
B68	Taeniasis	Chapter I	Certain infectious and parasitic diseases
B69	Cysticercosis	Chapter I	Certain infectious and parasitic diseases
B70	Diphyllobothriasis and sparganosis	Chapter I	Certain infectious and parasitic diseases
B71	Other cestode infections	Chapter I	Certain infectious and parasitic diseases
B72	Dracunculiasis	Chapter I	Certain infectious and parasitic diseases
B73	Onchocerciasis	Chapter I	Certain infectious and parasitic diseases
B74	Filariasis	Chapter I	Certain infectious and parasitic diseases
B75	Trichinellosis	Chapter I	Certain infectious and parasitic diseases
B76	Hookworm diseases	Chapter I	Certain infectious and parasitic diseases
B77	Ascariasis	Chapter I	Certain infectious and parasitic diseases
B78	Strongyloidiasis	Chapter I	Certain infectious and parasitic diseases
B79	Trichuriasis	Chapter I	Certain infectious and parasitic diseases
B80	Enterobiasis	Chapter I	Certain infectious and parasitic diseases
B81	Other intestinal helminthiases, not elsewhere classified	Chapter I	Certain infectious and parasitic diseases
B82	Unspecified intestinal parasitism	Chapter I	Certain infectious and parasitic diseases
B83	Other helminthiases	Chapter I	Certain infectious and parasitic diseases
B85	Pediculosis and phthiriasis	Chapter I	Certain infectious and parasitic diseases
B86	Scabies	Chapter I	Certain infectious and parasitic diseases
B87	Myiasis	Chapter I	Certain infectious and parasitic diseases
B88	Other infestations	Chapter I	Certain infectious and parasitic diseases
B89	Unspecified parasitic disease	Chapter I	Certain infectious and parasitic diseases
B90	Sequelae of tuberculosis	Chapter I	Certain infectious and parasitic diseases
B91	Sequelae of poliomyelitis	Chapter I	Certain infectious and parasitic diseases
B92	Sequelae of leprosy	Chapter I	Certain infectious and parasitic diseases
B94	Sequelae of other and unspecified infectious and parasitic diseases	Chapter I	Certain infectious and parasitic diseases
B95	Streptococcus and staphylococcus as the cause of diseases classified to other chapters	Chapter I	Certain infectious and parasitic diseases
B96	Other bacterial agents as the cause of diseases classified to other chapters	Chapter I	Certain infectious and parasitic diseases
B97	Viral agents as the cause of diseases classified to other chapters	Chapter I	Certain infectious and parasitic diseases
B98	Other specified infectious agents as the cause of diseases classified to other chapters	Chapter I	Certain infectious and parasitic diseases
B99	Other and unspecified infectious diseases	Chapter I	Certain infectious and parasitic diseases
C00	Malignant neoplasm of lip	Chapter II	Neoplasms
C01	Malignant neoplasm of base of tongue	Chapter II	Neoplasms
C02	Malignant neoplasm of other and unspecified parts of tongue	Chapter II	Neoplasms
C03	Malignant neoplasm of gum	Chapter II	Neoplasms
C04	Malignant neoplasm of floor of mouth	Chapter II	Neoplasms
C05	Malignant neoplasm of palate	Chapter II	Neoplasms
C06	Malignant neoplasm of other and unspecified parts of mouth	Chapter II	Neoplasms
C07	Malignant neoplasm of parotid gland	Chapter II	Neoplasms
C08	Malignant neoplasm of other and unspecified major salivary glands	Chapter II	Neoplasms
C09	Malignant neoplasm of tonsil	Chapter II	Neoplasms
C10	Malignant neoplasm of oropharynx	Chapter II	Neoplasms
C11	Malignant neoplasm of nasopharynx	Chapter II	Neoplasms
C12	Malignant neoplasm of pyriform sinus	Chapter II	Neoplasms
C13	Malignant neoplasm of hypopharynx	Chapter II	Neoplasms
C14	Malignant neoplasm of other and ill-defined sites in the lip, oral cavity and pharynx	Chapter II	Neoplasms
C15	Malignant neoplasm of oesophagus	Chapter II	Neoplasms
C16	Malignant neoplasm of stomach	Chapter II	Neoplasms
C17	Malignant neoplasm of small intestine	Chapter II	Neoplasms
C18	Malignant neoplasm of colon	Chapter II	Neoplasms
C19	Malignant neoplasm of rectosigmoid junction	Chapter II	Neoplasms
C20	Malignant neoplasm of rectum	Chapter II	Neoplasms
C21	Malignant neoplasm of anus and anal canal	Chapter II	Neoplasms
C22	Malignant neoplasm of liver and intrahepatic bile ducts	Chapter II	Neoplasms
C23	Malignant neoplasm of gallbladder	Chapter II	Neoplasms
C24	Malignant neoplasm of other and unspecified parts of biliary tract	Chapter II	Neoplasms
C25	Malignant neoplasm of pancreas	Chapter II	Neoplasms
C26	Malignant neoplasm of other and ill-defined digestive organs	Chapter II	Neoplasms
C30	Malignant neoplasm of nasal cavity and middle ear	Chapter II	Neoplasms
C31	Malignant neoplasm of accessory sinuses	Chapter II	Neoplasms
C32	Malignant neoplasm of larynx	Chapter II	Neoplasms
C33	Malignant neoplasm of trachea	Chapter II	Neoplasms
C34	Malignant neoplasm of bronchus and lung	Chapter II	Neoplasms
C37	Malignant neoplasm of thymus	Chapter II	Neoplasms
C38	Malignant neoplasm of heart, mediastinum and pleura	Chapter II	Neoplasms
C39	Malignant neoplasm of other and ill-defined sites in the respiratory system and intrathoracic organs	Chapter II	Neoplasms
C40	Malignant neoplasm of bone and articular cartilage of limbs	Chapter II	Neoplasms
C41	Malignant neoplasm of bone and articular cartilage of other and unspecified sites	Chapter II	Neoplasms
C42	hematopoietic and reticuloendothelial systems (ICD-O-3 specific)	Chapter II	Neoplasms
C43	Malignant melanoma of skin	Chapter II	Neoplasms
C44	Other malignant neoplasms of skin	Chapter II	Neoplasms
C45	Mesothelioma	Chapter II	Neoplasms
C46	Kaposi's sarcoma	Chapter II	Neoplasms
C47	Malignant neoplasm of peripheral nerves and autonomic nervous system	Chapter II	Neoplasms
C48	Malignant neoplasm of retroperitoneum and peritoneum	Chapter II	Neoplasms
C49	Malignant neoplasm of other connective and soft tissue	Chapter II	Neoplasms
C50	Malignant neoplasm of breast	Chapter II	Neoplasms
C51	Malignant neoplasm of vulva	Chapter II	Neoplasms
C52	Malignant neoplasm of vagina	Chapter II	Neoplasms
C53	Malignant neoplasm of cervix uteri	Chapter II	Neoplasms
C54	Malignant neoplasm of corpus uteri	Chapter II	Neoplasms
C55	Malignant neoplasm of uterus, part unspecified	Chapter II	Neoplasms
C56	Malignant neoplasm of ovary	Chapter II	Neoplasms
C57	Malignant neoplasm of other and unspecified female genital organs	Chapter II	Neoplasms
C58	Malignant neoplasm of placenta	Chapter II	Neoplasms
C60	Malignant neoplasm of penis	Chapter II	Neoplasms
C61	Malignant neoplasm of prostate	Chapter II	Neoplasms
C62	Malignant neoplasm of testis	Chapter II	Neoplasms
C63	Malignant neoplasm of other and unspecified male genital organs	Chapter II	Neoplasms
C64	Malignant neoplasm of kidney, except renal pelvis	Chapter II	Neoplasms
C65	Malignant neoplasm of renal pelvis	Chapter II	Neoplasms
C66	Malignant neoplasm of ureter	Chapter II	Neoplasms
C67	Malignant neoplasm of bladder	Chapter II	Neoplasms
C68	Malignant neoplasm of other and unspecified urinary organs	Chapter II	Neoplasms
C69	Malignant neoplasm of eye and adnexa	Chapter II	Neoplasms
C70	Malignant neoplasm of meninges	Chapter II	Neoplasms
C71	Malignant neoplasm of brain	Chapter II	Neoplasms
C72	Malignant neoplasm of spinal cord, cranial nerves and other parts of central nervous system	Chapter II	Neoplasms
C73	Malignant neoplasm of thyroid gland	Chapter II	Neoplasms
C74	Malignant neoplasm of adrenal gland	Chapter II	Neoplasms
C75	Malignant neoplasm of other endocrine glands and related structures	Chapter II	Neoplasms
C76	Malignant neoplasm of other and ill-defined sites	Chapter II	Neoplasms
C77	Secondary and unspecified malignant neoplasm of lymph nodes	Chapter II	Neoplasms
C78	Secondary malignant neoplasm of respiratory and digestive organs	Chapter II	Neoplasms
C79	Secondary malignant neoplasm of other sites	Chapter II	Neoplasms
C80	Malignant neoplasm without specification of site	Chapter II	Neoplasms
C81	Hodgkin's disease	Chapter II	Neoplasms
C82	Follicular [nodular] non-Hodgkin's lymphoma	Chapter II	Neoplasms
C83	Diffuse non-Hodgkin's lymphoma	Chapter II	Neoplasms
C84	Peripheral and cutaneous T-cell lymphomas	Chapter II	Neoplasms
C85	Other and unspecified types of non-Hodgkin's lymphoma	Chapter II	Neoplasms
C86	Other specified types of T/NK-cell lymphoma	Chapter II	Neoplasms
C88	Malignant immunoproliferative diseases	Chapter II	Neoplasms
C90	Multiple myeloma and malignant plasma cell neoplasms	Chapter II	Neoplasms
C91	Lymphoid leukaemia	Chapter II	Neoplasms
C92	Myeloid leukaemia	Chapter II	Neoplasms
C93	Monocytic leukaemia	Chapter II	Neoplasms
C94	Other leukaemias of specified cell type	Chapter II	Neoplasms
C95	Leukaemia of unspecified cell type	Chapter II	Neoplasms
C96	Other and unspecified malignant neoplasms of lymphoid, haematopoietic and related tissue	Chapter II	Neoplasms
C97	Malignant neoplasms of independent (primary) multiple sites	Chapter II	Neoplasms
D00	Carcinoma in situ of oral cavity, oesophagus and stomach	Chapter II	Neoplasms
D01	Carcinoma in situ of other and unspecified digestive organs	Chapter II	Neoplasms
D02	Carcinoma in situ of middle ear and respiratory system	Chapter II	Neoplasms
D03	Melanoma in situ	Chapter II	Neoplasms
D04	Carcinoma in situ of skin	Chapter II	Neoplasms
D05	Carcinoma in situ of breast	Chapter II	Neoplasms
D06	Carcinoma in situ of cervix uteri	Chapter II	Neoplasms
D07	Carcinoma in situ of other and unspecified genital organs	Chapter II	Neoplasms
D09	Carcinoma in situ of other and unspecified sites	Chapter II	Neoplasms
D10	Benign neoplasm of mouth and pharynx	Chapter II	Neoplasms
D11	Benign neoplasm of major salivary glands	Chapter II	Neoplasms
D12	Benign neoplasm of colon, rectum, anus and anal canal	Chapter II	Neoplasms
D13	Benign neoplasm of other and ill-defined parts of digestive system	Chapter II	Neoplasms
D14	Benign neoplasm of middle ear and respiratory system	Chapter II	Neoplasms
D15	Benign neoplasm of other and unspecified intrathoracic organs	Chapter II	Neoplasms
D16	Benign neoplasm of bone and articular cartilage	Chapter II	Neoplasms
D17	Benign lipomatous neoplasm	Chapter II	Neoplasms
D18	Haemangioma and lymphangioma, any site	Chapter II	Neoplasms
D19	Benign neoplasm of mesothelial tissue	Chapter II	Neoplasms
D20	Benign neoplasm of soft tissue of retroperitoneum and peritoneum	Chapter II	Neoplasms
D21	Other benign neoplasms of connective and other soft tissue	Chapter II	Neoplasms
D22	Melanocytic naevi	Chapter II	Neoplasms
D23	Other benign neoplasms of skin	Chapter II	Neoplasms
D24	Benign neoplasm of breast	Chapter II	Neoplasms
D25	Leiomyoma of uterus	Chapter II	Neoplasms
D26	Other benign neoplasms of uterus	Chapter II	Neoplasms
D27	Benign neoplasm of ovary	Chapter II	Neoplasms
D28	Benign neoplasm of other and unspecified female genital organs	Chapter II	Neoplasms
D29	Benign neoplasm of male genital organs	Chapter II	Neoplasms
D30	Benign neoplasm of urinary organs	Chapter II	Neoplasms
D31	Benign neoplasm of eye and adnexa	Chapter II	Neoplasms
D32	Benign neoplasm of meninges	Chapter II	Neoplasms
D33	Benign neoplasm of brain and other parts of central nervous system	Chapter II	Neoplasms
D34	Benign neoplasm of thyroid gland	Chapter II	Neoplasms
D35	Benign neoplasm of other and unspecified endocrine glands	Chapter II	Neoplasms
D36	Benign neoplasm of other and unspecified sites	Chapter II	Neoplasms
D37	Neoplasm of uncertain or unknown behaviour of oral cavity and digestive organs	Chapter II	Neoplasms
D38	Neoplasm of uncertain or unknown behaviour of middle ear and respiratory and intrathoracic organs	Chapter II	Neoplasms
D39	Neoplasm of uncertain or unknown behaviour of female genital organs	Chapter II	Neoplasms
D40	Neoplasm of uncertain or unknown behaviour of male genital organs	Chapter II	Neoplasms
D41	Neoplasm of uncertain or unknown behaviour of urinary organs	Chapter II	Neoplasms
D42	Neoplasm of uncertain or unknown behaviour of meninges	Chapter II	Neoplasms
D43	Neoplasm of uncertain or unknown behaviour of brain and central nervous system	Chapter II	Neoplasms
D44	Neoplasm of uncertain or unknown behaviour of endocrine glands	Chapter II	Neoplasms
D45	Polycythaemia vera	Chapter II	Neoplasms
D46	Myelodysplastic syndromes	Chapter II	Neoplasms
D47	Other neoplasms of uncertain or unknown behaviour of lymphoid, haematopoietic and related tissue	Chapter II	Neoplasms
D48	Neoplasm of uncertain or unknown behaviour of other and unspecified sites	Chapter II	Neoplasms
D50	Iron deficiency anaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D51	Vitamin B12 deficiency anaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D52	Folate deficiency anaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D53	Other nutritional anaemias	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D55	Anaemia due to enzyme disorders	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D56	Thalassaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D57	Sickle-cell disorders	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D58	Other hereditary haemolytic anaemias	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D59	Acquired haemolytic anaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D60	Acquired pure red cell aplasia [erythroblastopenia]	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D61	Other aplastic anaemias	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D62	Acute posthaemorrhagic anaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D63	Anaemia in chronic diseases classified elsewhere	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D64	Other anaemias	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D65	Disseminated intravascular coagulation [defibrination syndrome]	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D66	Hereditary factor VIII deficiency	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D67	Hereditary factor IX deficiency	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D68	Other coagulation defects	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D69	Purpura and other haemorrhagic conditions	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D70	Agranulocytosis	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D71	Functional disorders of polymorphonuclear neutrophils	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D72	Other disorders of white blood cells	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D73	Diseases of spleen	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D74	Methaemoglobinaemia	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D75	Other diseases of blood and blood-forming organs	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D76	Certain diseases involving lymphoreticular tissue and reticulohistiocytic system	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D77	Other disorders of blood and blood-forming organs in diseases classified elsewhere	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D80	Immunodeficiency with predominantly antibody defects	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D81	Combined immunodeficiencies	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D82	Immunodeficiency associated with other major defects	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D83	Common variable immunodeficiency	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D84	Other immunodeficiencies	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D86	Sarcoidosis	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
D89	Other disorders involving the immune mechanism, not elsewhere classified	Chapter III	Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism
E00	Congenital iodine-deficiency syndrome	Chapter IV	Endocrine, nutritional and metabolic diseases
E01	Iodine-deficiency-related thyroid disorders and allied conditions	Chapter IV	Endocrine, nutritional and metabolic diseases
E02	Subclinical iodine-deficiency hypothyroidism	Chapter IV	Endocrine, nutritional and metabolic diseases
E03	Other hypothyroidism	Chapter IV	Endocrine, nutritional and metabolic diseases
E04	Other non-toxic goitre	Chapter IV	Endocrine, nutritional and metabolic diseases
E05	Thyrotoxicosis [hyperthyroidism]	Chapter IV	Endocrine, nutritional and metabolic diseases
E06	Thyroiditis	Chapter IV	Endocrine, nutritional and metabolic diseases
E07	Other disorders of thyroid	Chapter IV	Endocrine, nutritional and metabolic diseases
E10	Insulin-dependent diabetes mellitus	Chapter IV	Endocrine, nutritional and metabolic diseases
E11	Non-insulin-dependent diabetes mellitus	Chapter IV	Endocrine, nutritional and metabolic diseases
E12	Malnutrition-related diabetes mellitus	Chapter IV	Endocrine, nutritional and metabolic diseases
E13	Other specified diabetes mellitus	Chapter IV	Endocrine, nutritional and metabolic diseases
E14	Unspecified diabetes mellitus	Chapter IV	Endocrine, nutritional and metabolic diseases
E15	Nondiabetic hypoglycaemic coma	Chapter IV	Endocrine, nutritional and metabolic diseases
E16	Other disorders of pancreatic internal secretion	Chapter IV	Endocrine, nutritional and metabolic diseases
E20	Hypoparathyroidism	Chapter IV	Endocrine, nutritional and metabolic diseases
E21	Hyperparathyroidism and other disorders of parathyroid gland	Chapter IV	Endocrine, nutritional and metabolic diseases
E22	Hyperfunction of pituitary gland	Chapter IV	Endocrine, nutritional and metabolic diseases
E23	Hypofunction and other disorders of pituitary gland	Chapter IV	Endocrine, nutritional and metabolic diseases
E24	Cushing's syndrome	Chapter IV	Endocrine, nutritional and metabolic diseases
E25	Adrenogenital disorders	Chapter IV	Endocrine, nutritional and metabolic diseases
E26	Hyperaldosteronism	Chapter IV	Endocrine, nutritional and metabolic diseases
E27	Other disorders of adrenal gland	Chapter IV	Endocrine, nutritional and metabolic diseases
E28	Ovarian dysfunction	Chapter IV	Endocrine, nutritional and metabolic diseases
E29	Testicular dysfunction	Chapter IV	Endocrine, nutritional and metabolic diseases
E30	Disorders of puberty, not elsewhere classified	Chapter IV	Endocrine, nutritional and metabolic diseases
E31	Polyglandular dysfunction	Chapter IV	Endocrine, nutritional and metabolic diseases
E32	Diseases of thymus	Chapter IV	Endocrine, nutritional and metabolic diseases
E34	Other endocrine disorders	Chapter IV	Endocrine, nutritional and metabolic diseases
E35	Disorders of endocrine glands in diseases classified elsewhere	Chapter IV	Endocrine, nutritional and metabolic diseases
E40	Kwashiorkor	Chapter IV	Endocrine, nutritional and metabolic diseases
E41	Nutritional marasmus	Chapter IV	Endocrine, nutritional and metabolic diseases
E42	Marasmic kwashiorkor	Chapter IV	Endocrine, nutritional and metabolic diseases
E43	Unspecified severe protein-energy malnutrition	Chapter IV	Endocrine, nutritional and metabolic diseases
E44	Protein-energy malnutrition of moderate and mild degree	Chapter IV	Endocrine, nutritional and metabolic diseases
E45	Retarded development following protein-energy malnutrition	Chapter IV	Endocrine, nutritional and metabolic diseases
E46	Unspecified protein-energy malnutrition	Chapter IV	Endocrine, nutritional and metabolic diseases
E50	Vitamin A deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E51	Thiamine deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E52	Niacin deficiency [pellagra]	Chapter IV	Endocrine, nutritional and metabolic diseases
E53	Deficiency of other B group vitamins	Chapter IV	Endocrine, nutritional and metabolic diseases
E54	Ascorbic acid deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E55	Vitamin D deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E56	Other vitamin deficiencies	Chapter IV	Endocrine, nutritional and metabolic diseases
E58	Dietary calcium deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E59	Dietary selenium deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E60	Dietary zinc deficiency	Chapter IV	Endocrine, nutritional and metabolic diseases
E61	Deficiency of other nutrient elements	Chapter IV	Endocrine, nutritional and metabolic diseases
E63	Other nutritional deficiencies	Chapter IV	Endocrine, nutritional and metabolic diseases
E64	Sequelae of malnutrition and other nutritional deficiencies	Chapter IV	Endocrine, nutritional and metabolic diseases
E65	Localised adiposity	Chapter IV	Endocrine, nutritional and metabolic diseases
E66	Obesity	Chapter IV	Endocrine, nutritional and metabolic diseases
E67	Other hyperalimentation	Chapter IV	Endocrine, nutritional and metabolic diseases
E68	Sequelae of hyperalimentation	Chapter IV	Endocrine, nutritional and metabolic diseases
E70	Disorders of aromatic amino-acid metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E71	Disorders of branched-chain amino-acid metabolism and fatty-acid metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E72	Other disorders of amino-acid metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E73	Lactose intolerance	Chapter IV	Endocrine, nutritional and metabolic diseases
E74	Other disorders of carbohydrate metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E75	Disorders of sphingolipid metabolism and other lipid storage disorders	Chapter IV	Endocrine, nutritional and metabolic diseases
E76	Disorders of glycosaminoglycan metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E77	Disorders of glycoprotein metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E78	Disorders of lipoprotein metabolism and other lipidaemias	Chapter IV	Endocrine, nutritional and metabolic diseases
E79	Disorders of purine and pyrimidine metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E80	Disorders of porphyrin and bilirubin metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E83	Disorders of mineral metabolism	Chapter IV	Endocrine, nutritional and metabolic diseases
E84	Cystic fibrosis	Chapter IV	Endocrine, nutritional and metabolic diseases
E85	Amyloidosis	Chapter IV	Endocrine, nutritional and metabolic diseases
E86	Volume depletion	Chapter IV	Endocrine, nutritional and metabolic diseases
E87	Other disorders of fluid, electrolyte and acid-base balance	Chapter IV	Endocrine, nutritional and metabolic diseases
E88	Other metabolic disorders	Chapter IV	Endocrine, nutritional and metabolic diseases
E89	Postprocedural endocrine and metabolic disorders, not elsewhere classified	Chapter IV	Endocrine, nutritional and metabolic diseases
E90	Nutritional and metabolic disorders in diseases classified elsewhere	Chapter IV	Endocrine, nutritional and metabolic diseases
F00	Dementia in Alzheimer's disease	Chapter V	Mental and behavioural disorders
F01	Vascular dementia	Chapter V	Mental and behavioural disorders
F02	Dementia in other diseases classified elsewhere	Chapter V	Mental and behavioural disorders
F03	Unspecified dementia	Chapter V	Mental and behavioural disorders
F04	Organic amnesic syndrome, not induced by alcohol and other psychoactive substances	Chapter V	Mental and behavioural disorders
F05	Delirium, not induced by alcohol and other psychoactive substances	Chapter V	Mental and behavioural disorders
F06	Other mental disorders due to brain damage and dysfunction and to physical disease	Chapter V	Mental and behavioural disorders
F07	Personality and behavioural disorders due to brain disease, damage and dysfunction	Chapter V	Mental and behavioural disorders
F09	Unspecified organic or symptomatic mental disorder	Chapter V	Mental and behavioural disorders
F10	Mental and behavioural disorders due to use of alcohol	Chapter V	Mental and behavioural disorders
F11	Mental and behavioural disorders due to use of opioids	Chapter V	Mental and behavioural disorders
F12	Mental and behavioural disorders due to use of cannabinoids	Chapter V	Mental and behavioural disorders
F13	Mental and behavioural disorders due to use of sedatives or hypnotics	Chapter V	Mental and behavioural disorders
F14	Mental and behavioural disorders due to use of cocaine	Chapter V	Mental and behavioural disorders
F15	Mental and behavioural disorders due to use of other stimulants, including caffeine	Chapter V	Mental and behavioural disorders
F16	Mental and behavioural disorders due to use of hallucinogens	Chapter V	Mental and behavioural disorders
F17	Mental and behavioural disorders due to use of tobacco	Chapter V	Mental and behavioural disorders
F18	Mental and behavioural disorders due to use of volatile solvents	Chapter V	Mental and behavioural disorders
F19	Mental and behavioural disorders due to multiple drug use and use of other psychoactive substances	Chapter V	Mental and behavioural disorders
F20	Schizophrenia	Chapter V	Mental and behavioural disorders
F21	Schizotypal disorder	Chapter V	Mental and behavioural disorders
F22	Persistent delusional disorders	Chapter V	Mental and behavioural disorders
F23	Acute and transient psychotic disorders	Chapter V	Mental and behavioural disorders
F24	Induced delusional disorder	Chapter V	Mental and behavioural disorders
F25	Schizoaffective disorders	Chapter V	Mental and behavioural disorders
F28	Other nonorganic psychotic disorders	Chapter V	Mental and behavioural disorders
F29	Unspecified nonorganic psychosis	Chapter V	Mental and behavioural disorders
F30	Manic episode	Chapter V	Mental and behavioural disorders
F31	Bipolar affective disorder	Chapter V	Mental and behavioural disorders
F32	Depressive episode	Chapter V	Mental and behavioural disorders
F33	Recurrent depressive disorder	Chapter V	Mental and behavioural disorders
F34	Persistent mood [affective] disorders	Chapter V	Mental and behavioural disorders
F38	Other mood [affective] disorders	Chapter V	Mental and behavioural disorders
F39	Unspecified mood [affective] disorder	Chapter V	Mental and behavioural disorders
F40	Phobic anxiety disorders	Chapter V	Mental and behavioural disorders
F41	Other anxiety disorders	Chapter V	Mental and behavioural disorders
F42	Obsessive-compulsive disorder	Chapter V	Mental and behavioural disorders
F43	Reaction to severe stress, and adjustment disorders	Chapter V	Mental and behavioural disorders
F44	Dissociative [conversion] disorders	Chapter V	Mental and behavioural disorders
F45	Somatoform disorders	Chapter V	Mental and behavioural disorders
F48	Other neurotic disorders	Chapter V	Mental and behavioural disorders
F50	Eating disorders	Chapter V	Mental and behavioural disorders
F51	Nonorganic sleep disorders	Chapter V	Mental and behavioural disorders
F52	Sexual dysfunction, not caused by organic disorder or disease	Chapter V	Mental and behavioural disorders
F53	Mental and behavioural disorders associated with the puerperium, not elsewhere classified	Chapter V	Mental and behavioural disorders
F54	Psychological and behavioural factors associated with disorders or diseases classified elsewhere	Chapter V	Mental and behavioural disorders
F55	Abuse of non-dependence-producing substances	Chapter V	Mental and behavioural disorders
F59	Unspecified behavioural syndromes associated with physiological disturbances and physical factors	Chapter V	Mental and behavioural disorders
F60	Specific personality disorders	Chapter V	Mental and behavioural disorders
F61	Mixed and other personality disorders	Chapter V	Mental and behavioural disorders
F62	Enduring personality changes, not attributable to brain damage and disease	Chapter V	Mental and behavioural disorders
F63	Habit and impulse disorders	Chapter V	Mental and behavioural disorders
F64	Gender identity disorders	Chapter V	Mental and behavioural disorders
F65	Disorders of sexual preference	Chapter V	Mental and behavioural disorders
F66	Psychological and behavioural disorders associated with sexual development and orientation	Chapter V	Mental and behavioural disorders
F68	Other disorders of adult personality and behaviour	Chapter V	Mental and behavioural disorders
F69	Unspecified disorder of adult personality and behaviour	Chapter V	Mental and behavioural disorders
F70	Mild mental retardation	Chapter V	Mental and behavioural disorders
F71	Moderate mental retardation	Chapter V	Mental and behavioural disorders
F72	Severe mental retardation	Chapter V	Mental and behavioural disorders
F73	Profound mental retardation	Chapter V	Mental and behavioural disorders
F78	Other mental retardation	Chapter V	Mental and behavioural disorders
F79	Unspecified mental retardation	Chapter V	Mental and behavioural disorders
F80	Specific developmental disorders of speech and language	Chapter V	Mental and behavioural disorders
F81	Specific developmental disorders of scholastic skills	Chapter V	Mental and behavioural disorders
F82	Specific developmental disorder of motor function	Chapter V	Mental and behavioural disorders
F83	Mixed specific developmental disorders	Chapter V	Mental and behavioural disorders
F84	Pervasive developmental disorders	Chapter V	Mental and behavioural disorders
F88	Other disorders of psychological development	Chapter V	Mental and behavioural disorders
F89	Unspecified disorder of psychological development	Chapter V	Mental and behavioural disorders
F90	Hyperkinetic disorders	Chapter V	Mental and behavioural disorders
F91	Conduct disorders	Chapter V	Mental and behavioural disorders
F92	Mixed disorders of conduct and emotions	Chapter V	Mental and behavioural disorders
F93	Emotional disorders with onset specific to childhood	Chapter V	Mental and behavioural disorders
F94	Disorders of social functioning with onset specific to childhood and adolescence	Chapter V	Mental and behavioural disorders
F95	Tic disorders	Chapter V	Mental and behavioural disorders
F98	Other behavioural and emotional disorders with onset usually occurring in childhood and adolescence	Chapter V	Mental and behavioural disorders
F99	Mental disorder, not otherwise specified	Chapter V	Mental and behavioural disorders
G00	Bacterial meningitis, not elsewhere classified	Chapter VI	Diseases of the nervous system
G01	Meningitis in bacterial diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G02	Meningitis in other infectious and parasitic diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G03	Meningitis due to other and unspecified causes	Chapter VI	Diseases of the nervous system
G04	Encephalitis, myelitis and encephalomyelitis	Chapter VI	Diseases of the nervous system
G05	Encephalitis, myelitis and encephalomyelitis in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G06	Intracranial and intraspinal abscess and granuloma	Chapter VI	Diseases of the nervous system
G07	Intracranial and intraspinal abscess and granuloma in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G08	Intracranial and intraspinal phlebitis and thrombophlebitis	Chapter VI	Diseases of the nervous system
G09	Sequelae of inflammatory diseases of central nervous system	Chapter VI	Diseases of the nervous system
G10	Huntington's disease	Chapter VI	Diseases of the nervous system
G11	Hereditary ataxia	Chapter VI	Diseases of the nervous system
G12	Spinal muscular atrophy and related syndromes	Chapter VI	Diseases of the nervous system
G13	Systemic atrophies primarily affecting central nervous system in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G14	Postpolio syndrome	Chapter VI	Diseases of the nervous system
G20	Parkinson's disease	Chapter VI	Diseases of the nervous system
G21	Secondary Parkinsonism	Chapter VI	Diseases of the nervous system
G22	Parkinsonism in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G23	Other degenerative diseases of basal ganglia	Chapter VI	Diseases of the nervous system
G24	Dystonia	Chapter VI	Diseases of the nervous system
G25	Other extrapyramidal and movement disorders	Chapter VI	Diseases of the nervous system
G26	Extrapyramidal and movement disorders in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G30	Alzheimer's disease	Chapter VI	Diseases of the nervous system
G31	Other degenerative diseases of nervous system, not elsewhere classified	Chapter VI	Diseases of the nervous system
G32	Other degenerative disorders of nervous system in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G35	Multiple sclerosis	Chapter VI	Diseases of the nervous system
G36	Other acute disseminated demyelination	Chapter VI	Diseases of the nervous system
G37	Other demyelinating diseases of central nervous system	Chapter VI	Diseases of the nervous system
G40	Epilepsy	Chapter VI	Diseases of the nervous system
G41	Status epilepticus	Chapter VI	Diseases of the nervous system
G43	Migraine	Chapter VI	Diseases of the nervous system
G44	Other headache syndromes	Chapter VI	Diseases of the nervous system
G45	Transient cerebral ischaemic attacks and related syndromes	Chapter VI	Diseases of the nervous system
G46	Vascular syndromes of brain in cerebrovascular diseases	Chapter VI	Diseases of the nervous system
G47	Sleep disorders	Chapter VI	Diseases of the nervous system
G50	Disorders of trigeminal nerve	Chapter VI	Diseases of the nervous system
G51	Facial nerve disorders	Chapter VI	Diseases of the nervous system
G52	Disorders of other cranial nerves	Chapter VI	Diseases of the nervous system
G53	Cranial nerve disorders in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G54	Nerve root and plexus disorders	Chapter VI	Diseases of the nervous system
G55	Nerve root and plexus compressions in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G56	Mononeuropathies of upper limb	Chapter VI	Diseases of the nervous system
G57	Mononeuropathies of lower limb	Chapter VI	Diseases of the nervous system
G58	Other mononeuropathies	Chapter VI	Diseases of the nervous system
G59	Mononeuropathy in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G60	Hereditary and idiopathic neuropathy	Chapter VI	Diseases of the nervous system
G61	Inflammatory polyneuropathy	Chapter VI	Diseases of the nervous system
G62	Other polyneuropathies	Chapter VI	Diseases of the nervous system
G63	Polyneuropathy in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G64	Other disorders of peripheral nervous system	Chapter VI	Diseases of the nervous system
G70	Myasthenia gravis and other myoneural disorders	Chapter VI	Diseases of the nervous system
G71	Primary disorders of muscles	Chapter VI	Diseases of the nervous system
G72	Other myopathies	Chapter VI	Diseases of the nervous system
G73	Disorders of myoneural junction and muscle in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G80	Infantile cerebral palsy	Chapter VI	Diseases of the nervous system
G81	Hemiplegia	Chapter VI	Diseases of the nervous system
G82	Paraplegia and tetraplegia	Chapter VI	Diseases of the nervous system
G83	Other paralytic syndromes	Chapter VI	Diseases of the nervous system
G90	Disorders of autonomic nervous system	Chapter VI	Diseases of the nervous system
G91	Hydrocephalus	Chapter VI	Diseases of the nervous system
G92	Toxic encephalopathy	Chapter VI	Diseases of the nervous system
G93	Other disorders of brain	Chapter VI	Diseases of the nervous system
G94	Other disorders of brain in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
G95	Other diseases of spinal cord	Chapter VI	Diseases of the nervous system
G96	Other disorders of central nervous system	Chapter VI	Diseases of the nervous system
G97	Postprocedural disorders of nervous system, not elsewhere classified	Chapter VI	Diseases of the nervous system
G98	Other disorders of nervous system, not elsewhere classified	Chapter VI	Diseases of the nervous system
G99	Other disorders of nervous system in diseases classified elsewhere	Chapter VI	Diseases of the nervous system
H00	Hordeolum and chalazion	Chapter VII	Diseases of the eye and adnexa
H01	Other inflammation of eyelid	Chapter VII	Diseases of the eye and adnexa
H02	Other disorders of eyelid	Chapter VII	Diseases of the eye and adnexa
H03	Disorders of eyelid in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H04	Disorders of lachrymal system	Chapter VII	Diseases of the eye and adnexa
H05	Disorders of orbit	Chapter VII	Diseases of the eye and adnexa
H06	Disorders of lachrymal system and orbit in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H10	Conjunctivitis	Chapter VII	Diseases of the eye and adnexa
H11	Other disorders of conjunctiva	Chapter VII	Diseases of the eye and adnexa
H13	Disorders of conjunctiva in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H15	Disorders of sclera	Chapter VII	Diseases of the eye and adnexa
H16	Keratitis	Chapter VII	Diseases of the eye and adnexa
H17	Corneal scars and opacities	Chapter VII	Diseases of the eye and adnexa
H18	Other disorders of cornea	Chapter VII	Diseases of the eye and adnexa
H19	Disorders of sclera and cornea in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H20	Iridocyclitis	Chapter VII	Diseases of the eye and adnexa
H21	Other disorders of iris and ciliary body	Chapter VII	Diseases of the eye and adnexa
H22	Disorders of iris and ciliary body in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H25	Senile cataract	Chapter VII	Diseases of the eye and adnexa
H26	Other cataract	Chapter VII	Diseases of the eye and adnexa
H27	Other disorders of lens	Chapter VII	Diseases of the eye and adnexa
H28	Cataract and other disorders of lens in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H30	Chorioretinal inflammation	Chapter VII	Diseases of the eye and adnexa
H31	Other disorders of choroid	Chapter VII	Diseases of the eye and adnexa
H32	Chorioretinal disorders in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H33	Retinal detachments and breaks	Chapter VII	Diseases of the eye and adnexa
H34	Retinal vascular occlusions	Chapter VII	Diseases of the eye and adnexa
H35	Other retinal disorders	Chapter VII	Diseases of the eye and adnexa
H36	Retinal disorders in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H40	Glaucoma	Chapter VII	Diseases of the eye and adnexa
H42	Glaucoma in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H43	Disorders of vitreous body	Chapter VII	Diseases of the eye and adnexa
H44	Disorders of globe	Chapter VII	Diseases of the eye and adnexa
H45	Disorders of vitreous body and globe in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H46	Optic neuritis	Chapter VII	Diseases of the eye and adnexa
H47	Other disorders of optic [2nd] nerve and visual pathways	Chapter VII	Diseases of the eye and adnexa
H48	Disorders of optic [2nd] nerve and visual pathways in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H49	Paralytic strabismus	Chapter VII	Diseases of the eye and adnexa
H50	Other strabismus	Chapter VII	Diseases of the eye and adnexa
H51	Other disorders of binocular movement	Chapter VII	Diseases of the eye and adnexa
H52	Disorders of refraction and accommodation	Chapter VII	Diseases of the eye and adnexa
H53	Visual disturbances	Chapter VII	Diseases of the eye and adnexa
H54	Blindness and low vision	Chapter VII	Diseases of the eye and adnexa
H55	Nystagmus and other irregular eye movements	Chapter VII	Diseases of the eye and adnexa
H57	Other disorders of eye and adnexa	Chapter VII	Diseases of the eye and adnexa
H58	Other disorders of eye and adnexa in diseases classified elsewhere	Chapter VII	Diseases of the eye and adnexa
H59	Postprocedural disorders of eye and adnexa, not elsewhere classified	Chapter VII	Diseases of the eye and adnexa
H60	Otitis externa	Chapter VIII	Diseases of the ear and mastoid process
H61	Other disorders of external ear	Chapter VIII	Diseases of the ear and mastoid process
H62	Disorders of external ear in diseases classified elsewhere	Chapter VIII	Diseases of the ear and mastoid process
H65	Nonsuppurative otitis media	Chapter VIII	Diseases of the ear and mastoid process
H66	Suppurative and unspecified otitis media	Chapter VIII	Diseases of the ear and mastoid process
H67	Otitis media in diseases classified elsewhere	Chapter VIII	Diseases of the ear and mastoid process
H68	Eustachian salpingitis and obstruction	Chapter VIII	Diseases of the ear and mastoid process
H69	Other disorders of Eustachian tube	Chapter VIII	Diseases of the ear and mastoid process
H70	Mastoiditis and related conditions	Chapter VIII	Diseases of the ear and mastoid process
H71	Cholesteatoma of middle ear	Chapter VIII	Diseases of the ear and mastoid process
H72	Perforation of tympanic membrane	Chapter VIII	Diseases of the ear and mastoid process
H73	Other disorders of tympanic membrane	Chapter VIII	Diseases of the ear and mastoid process
H74	Other disorders of middle ear and mastoid	Chapter VIII	Diseases of the ear and mastoid process
H75	Other disorders of middle ear and mastoid in diseases classified elsewhere	Chapter VIII	Diseases of the ear and mastoid process
H80	Otosclerosis	Chapter VIII	Diseases of the ear and mastoid process
H81	Disorders of vestibular function	Chapter VIII	Diseases of the ear and mastoid process
H82	Vertiginous syndromes in diseases classified elsewhere	Chapter VIII	Diseases of the ear and mastoid process
H83	Other diseases of inner ear	Chapter VIII	Diseases of the ear and mastoid process
H90	Conductive and sensorineural hearing loss	Chapter VIII	Diseases of the ear and mastoid process
H91	Other hearing loss	Chapter VIII	Diseases of the ear and mastoid process
H92	Otalgia and effusion of ear	Chapter VIII	Diseases of the ear and mastoid process
H93	Other disorders of ear, not elsewhere classified	Chapter VIII	Diseases of the ear and mastoid process
H94	Other disorders of ear in diseases classified elsewhere	Chapter VIII	Diseases of the ear and mastoid process
H95	Postprocedural disorders of ear and mastoid process, not elsewhere classified	Chapter VIII	Diseases of the ear and mastoid process
I00	Rheumatic fever without mention of heart involvement	Chapter IX	Diseases of the circulatory system
I01	Rheumatic fever with heart involvement	Chapter IX	Diseases of the circulatory system
I02	Rheumatic chorea	Chapter IX	Diseases of the circulatory system
I05	Rheumatic mitral valve diseases	Chapter IX	Diseases of the circulatory system
I06	Rheumatic aortic valve diseases	Chapter IX	Diseases of the circulatory system
I07	Rheumatic tricuspid valve diseases	Chapter IX	Diseases of the circulatory system
I08	Multiple valve diseases	Chapter IX	Diseases of the circulatory system
I09	Other rheumatic heart diseases	Chapter IX	Diseases of the circulatory system
I10	Essential (primary) hypertension	Chapter IX	Diseases of the circulatory system
I11	Hypertensive heart disease	Chapter IX	Diseases of the circulatory system
I12	Hypertensive renal disease	Chapter IX	Diseases of the circulatory system
I13	Hypertensive heart and renal disease	Chapter IX	Diseases of the circulatory system
I15	Secondary hypertension	Chapter IX	Diseases of the circulatory system
I20	Angina pectoris	Chapter IX	Diseases of the circulatory system
I21	Acute myocardial infarction	Chapter IX	Diseases of the circulatory system
I22	Subsequent myocardial infarction	Chapter IX	Diseases of the circulatory system
I23	Certain current complications following acute myocardial infarction	Chapter IX	Diseases of the circulatory system
I24	Other acute ischaemic heart diseases	Chapter IX	Diseases of the circulatory system
I25	Chronic ischaemic heart disease	Chapter IX	Diseases of the circulatory system
I26	Pulmonary embolism	Chapter IX	Diseases of the circulatory system
I27	Other pulmonary heart diseases	Chapter IX	Diseases of the circulatory system
I28	Other diseases of pulmonary vessels	Chapter IX	Diseases of the circulatory system
I30	Acute pericarditis	Chapter IX	Diseases of the circulatory system
I31	Other diseases of pericardium	Chapter IX	Diseases of the circulatory system
I32	Pericarditis in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I33	Acute and subacute endocarditis	Chapter IX	Diseases of the circulatory system
I34	Nonrheumatic mitral valve disorders	Chapter IX	Diseases of the circulatory system
I35	Nonrheumatic aortic valve disorders	Chapter IX	Diseases of the circulatory system
I36	Nonrheumatic tricuspid valve disorders	Chapter IX	Diseases of the circulatory system
I37	Pulmonary valve disorders	Chapter IX	Diseases of the circulatory system
I38	Endocarditis, valve unspecified	Chapter IX	Diseases of the circulatory system
I39	Endocarditis and heart valve disorders in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I40	Acute myocarditis	Chapter IX	Diseases of the circulatory system
I41	Myocarditis in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I42	Cardiomyopathy	Chapter IX	Diseases of the circulatory system
I43	Cardiomyopathy in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I44	Atrioventricular and left bundle-branch block	Chapter IX	Diseases of the circulatory system
I45	Other conduction disorders	Chapter IX	Diseases of the circulatory system
I46	Cardiac arrest	Chapter IX	Diseases of the circulatory system
I47	Paroxysmal tachycardia	Chapter IX	Diseases of the circulatory system
I48	Atrial fibrillation and flutter	Chapter IX	Diseases of the circulatory system
I49	Other cardiac arrhythmias	Chapter IX	Diseases of the circulatory system
I50	Heart failure	Chapter IX	Diseases of the circulatory system
I51	Complications and ill-defined descriptions of heart disease	Chapter IX	Diseases of the circulatory system
I52	Other heart disorders in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I60	Subarachnoid haemorrhage	Chapter IX	Diseases of the circulatory system
I61	Intracerebral haemorrhage	Chapter IX	Diseases of the circulatory system
I62	Other nontraumatic intracranial haemorrhage	Chapter IX	Diseases of the circulatory system
I63	Cerebral infarction	Chapter IX	Diseases of the circulatory system
I64	Stroke, not specified as haemorrhage or infarction	Chapter IX	Diseases of the circulatory system
I65	Occlusion and stenosis of precerebral arteries, not resulting in cerebral infarction	Chapter IX	Diseases of the circulatory system
I66	Occlusion and stenosis of cerebral arteries, not resulting in cerebral infarction	Chapter IX	Diseases of the circulatory system
I67	Other cerebrovascular diseases	Chapter IX	Diseases of the circulatory system
I68	Cerebrovascular disorders in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I69	Sequelae of cerebrovascular disease	Chapter IX	Diseases of the circulatory system
I70	Atherosclerosis	Chapter IX	Diseases of the circulatory system
I71	Aortic aneurysm and dissection	Chapter IX	Diseases of the circulatory system
I72	Other aneurysm	Chapter IX	Diseases of the circulatory system
I73	Other peripheral vascular diseases	Chapter IX	Diseases of the circulatory system
I74	Arterial embolism and thrombosis	Chapter IX	Diseases of the circulatory system
I77	Other disorders of arteries and arterioles	Chapter IX	Diseases of the circulatory system
I78	Diseases of capillaries	Chapter IX	Diseases of the circulatory system
I79	Disorders of arteries, arterioles and capillaries in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I80	Phlebitis and thrombophlebitis	Chapter IX	Diseases of the circulatory system
I81	Portal vein thrombosis	Chapter IX	Diseases of the circulatory system
I82	Other venous embolism and thrombosis	Chapter IX	Diseases of the circulatory system
I83	Varicose veins of lower extremities	Chapter IX	Diseases of the circulatory system
I84	Haemorrhoids	Chapter IX	Diseases of the circulatory system
I85	Oesophageal varices	Chapter IX	Diseases of the circulatory system
I86	Varicose veins of other sites	Chapter IX	Diseases of the circulatory system
I87	Other disorders of veins	Chapter IX	Diseases of the circulatory system
I88	Nonspecific lymphadenitis	Chapter IX	Diseases of the circulatory system
I89	Other non-infective disorders of lymphatic vessels and lymph nodes	Chapter IX	Diseases of the circulatory system
I95	Hypotension	Chapter IX	Diseases of the circulatory system
I97	Postprocedural disorders of circulatory system, not elsewhere classified	Chapter IX	Diseases of the circulatory system
I98	Other disorders of circulatory system in diseases classified elsewhere	Chapter IX	Diseases of the circulatory system
I99	Other and unspecified disorders of circulatory system	Chapter IX	Diseases of the circulatory system
J00	Acute nasopharyngitis [common cold]	Chapter X	Diseases of the respiratory system
J01	Acute sinusitis	Chapter X	Diseases of the respiratory system
J02	Acute pharyngitis	Chapter X	Diseases of the respiratory system
J03	Acute tonsillitis	Chapter X	Diseases of the respiratory system
J04	Acute laryngitis and tracheitis	Chapter X	Diseases of the respiratory system
J05	Acute obstructive laryngitis [croup] and epiglottitis	Chapter X	Diseases of the respiratory system
J06	Acute upper respiratory infections of multiple and unspecified sites	Chapter X	Diseases of the respiratory system
J09	Influenza due to certain identified influenza virus	Chapter X	Diseases of the respiratory system
J10	Influenza due to identified influenza virus	Chapter X	Diseases of the respiratory system
J11	Influenza, virus not identified	Chapter X	Diseases of the respiratory system
J12	Viral pneumonia, not elsewhere classified	Chapter X	Diseases of the respiratory system
J13	Pneumonia due to Streptococcus pneumoniae	Chapter X	Diseases of the respiratory system
J14	Pneumonia due to Haemophilus influenzae	Chapter X	Diseases of the respiratory system
J15	Bacterial pneumonia, not elsewhere classified	Chapter X	Diseases of the respiratory system
J16	Pneumonia due to other infectious organisms, not elsewhere classified	Chapter X	Diseases of the respiratory system
J17	Pneumonia in diseases classified elsewhere	Chapter X	Diseases of the respiratory system
J18	Pneumonia, organism unspecified	Chapter X	Diseases of the respiratory system
J20	Acute bronchitis	Chapter X	Diseases of the respiratory system
J21	Acute bronchiolitis	Chapter X	Diseases of the respiratory system
J22	Unspecified acute lower respiratory infection	Chapter X	Diseases of the respiratory system
J30	Vasomotor and allergic rhinitis	Chapter X	Diseases of the respiratory system
J31	Chronic rhinitis, nasopharyngitis and pharyngitis	Chapter X	Diseases of the respiratory system
J32	Chronic sinusitis	Chapter X	Diseases of the respiratory system
J33	Nasal polyp	Chapter X	Diseases of the respiratory system
J34	Other disorders of nose and nasal sinuses	Chapter X	Diseases of the respiratory system
J35	Chronic diseases of tonsils and adenoids	Chapter X	Diseases of the respiratory system
J36	Peritonsillar abscess	Chapter X	Diseases of the respiratory system
J37	Chronic laryngitis and laryngotracheitis	Chapter X	Diseases of the respiratory system
J38	Diseases of vocal cords and larynx, not elsewhere classified	Chapter X	Diseases of the respiratory system
J39	Other diseases of upper respiratory tract	Chapter X	Diseases of the respiratory system
J40	Bronchitis, not specified as acute or chronic	Chapter X	Diseases of the respiratory system
J41	Simple and mucopurulent chronic bronchitis	Chapter X	Diseases of the respiratory system
J42	Unspecified chronic bronchitis	Chapter X	Diseases of the respiratory system
J43	Emphysema	Chapter X	Diseases of the respiratory system
J44	Other chronic obstructive pulmonary disease	Chapter X	Diseases of the respiratory system
J45	Asthma	Chapter X	Diseases of the respiratory system
J46	Status asthmaticus	Chapter X	Diseases of the respiratory system
J47	Bronchiectasis	Chapter X	Diseases of the respiratory system
J60	Coalworker's pneumoconiosis	Chapter X	Diseases of the respiratory system
J61	Pneumoconiosis due to asbestos and other mineral fibres	Chapter X	Diseases of the respiratory system
J62	Pneumoconiosis due to dust containing silica	Chapter X	Diseases of the respiratory system
J63	Pneumoconiosis due to other inorganic dusts	Chapter X	Diseases of the respiratory system
J64	Unspecified pneumoconiosis	Chapter X	Diseases of the respiratory system
J65	Pneumoconiosis associated with tuberculosis	Chapter X	Diseases of the respiratory system
J66	Airway disease due to specific organic dust	Chapter X	Diseases of the respiratory system
J67	Hypersensitivity pneumonitis due to organic dust	Chapter X	Diseases of the respiratory system
J68	Respiratory conditions due to inhalation of chemicals, gases, fumes and vapours	Chapter X	Diseases of the respiratory system
J69	Pneumonitis due to solids and liquids	Chapter X	Diseases of the respiratory system
J70	Respiratory conditions due to other external agents	Chapter X	Diseases of the respiratory system
J80	Adult respiratory distress syndrome	Chapter X	Diseases of the respiratory system
J81	Pulmonary oedema	Chapter X	Diseases of the respiratory system
J82	Pulmonary eosinophilia, not elsewhere classified	Chapter X	Diseases of the respiratory system
J84	Other interstitial pulmonary diseases	Chapter X	Diseases of the respiratory system
J85	Abscess of lung and mediastinum	Chapter X	Diseases of the respiratory system
J86	Pyothorax	Chapter X	Diseases of the respiratory system
J90	Pleural effusion, not elsewhere classified	Chapter X	Diseases of the respiratory system
J91	Pleural effusion in conditions classified elsewhere	Chapter X	Diseases of the respiratory system
J92	Pleural plaque	Chapter X	Diseases of the respiratory system
J93	Pneumothorax	Chapter X	Diseases of the respiratory system
J94	Other pleural conditions	Chapter X	Diseases of the respiratory system
J95	Postprocedural respiratory disorders, not elsewhere classified	Chapter X	Diseases of the respiratory system
J96	Respiratory failure, not elsewhere classified	Chapter X	Diseases of the respiratory system
J98	Other respiratory disorders	Chapter X	Diseases of the respiratory system
J99	Respiratory disorders in diseases classified elsewhere	Chapter X	Diseases of the respiratory system
K00	Disorders of tooth development and eruption	Chapter XI	Diseases of the digestive system
K01	Embedded and impacted teeth	Chapter XI	Diseases of the digestive system
K02	Dental caries	Chapter XI	Diseases of the digestive system
K03	Other diseases of hard tissues of teeth	Chapter XI	Diseases of the digestive system
K04	Diseases of pulp and periapical tissues	Chapter XI	Diseases of the digestive system
K05	Gingivitis and periodontal diseases	Chapter XI	Diseases of the digestive system
K06	Other disorders of gingiva and edentulous alveolar ridge	Chapter XI	Diseases of the digestive system
K07	Dentofacial anomalies [including malocclusion]	Chapter XI	Diseases of the digestive system
K08	Other disorders of teeth and supporting structures	Chapter XI	Diseases of the digestive system
K09	Cysts of oral region, not elsewhere classified	Chapter XI	Diseases of the digestive system
K10	Other diseases of jaws	Chapter XI	Diseases of the digestive system
K11	Diseases of salivary glands	Chapter XI	Diseases of the digestive system
K12	Stomatitis and related lesions	Chapter XI	Diseases of the digestive system
K13	Other diseases of lip and oral mucosa	Chapter XI	Diseases of the digestive system
K14	Diseases of tongue	Chapter XI	Diseases of the digestive system
K20	Oesophagitis	Chapter XI	Diseases of the digestive system
K21	Gastro-oesophageal reflux disease	Chapter XI	Diseases of the digestive system
K22	Other diseases of oesophagus	Chapter XI	Diseases of the digestive system
K23	Disorders of oesophagus in diseases classified elsewhere	Chapter XI	Diseases of the digestive system
K25	Gastric ulcer	Chapter XI	Diseases of the digestive system
K26	Duodenal ulcer	Chapter XI	Diseases of the digestive system
K27	Peptic ulcer, site unspecified	Chapter XI	Diseases of the digestive system
K28	Gastrojejunal ulcer	Chapter XI	Diseases of the digestive system
K29	Gastritis and duodenitis	Chapter XI	Diseases of the digestive system
K30	Dyspepsia	Chapter XI	Diseases of the digestive system
K31	Other diseases of stomach and duodenum	Chapter XI	Diseases of the digestive system
K35	Acute appendicitis	Chapter XI	Diseases of the digestive system
K36	Other appendicitis	Chapter XI	Diseases of the digestive system
K37	Unspecified appendicitis	Chapter XI	Diseases of the digestive system
K38	Other diseases of appendix	Chapter XI	Diseases of the digestive system
K40	Inguinal hernia	Chapter XI	Diseases of the digestive system
K41	Femoral hernia	Chapter XI	Diseases of the digestive system
K42	Umbilical hernia	Chapter XI	Diseases of the digestive system
K43	Ventral hernia	Chapter XI	Diseases of the digestive system
K44	Diaphragmatic hernia	Chapter XI	Diseases of the digestive system
K45	Other abdominal hernia	Chapter XI	Diseases of the digestive system
K46	Unspecified abdominal hernia	Chapter XI	Diseases of the digestive system
K50	Crohn's disease [regional enteritis]	Chapter XI	Diseases of the digestive system
K51	Ulcerative colitis	Chapter XI	Diseases of the digestive system
K52	Other non-infective gastro-enteritis and colitis	Chapter XI	Diseases of the digestive system
K55	Vascular disorders of intestine	Chapter XI	Diseases of the digestive system
K56	Paralytic ileus and intestinal obstruction without hernia	Chapter XI	Diseases of the digestive system
K57	Diverticular disease of intestine	Chapter XI	Diseases of the digestive system
K58	Irritable bowel syndrome	Chapter XI	Diseases of the digestive system
K59	Other functional intestinal disorders	Chapter XI	Diseases of the digestive system
K60	Fissure and fistula of anal and rectal regions	Chapter XI	Diseases of the digestive system
K61	Abscess of anal and rectal regions	Chapter XI	Diseases of the digestive system
K62	Other diseases of anus and rectum	Chapter XI	Diseases of the digestive system
K63	Other diseases of intestine	Chapter XI	Diseases of the digestive system
K64	Haemorrhoids and perianal venous thrombosis	Chapter XI	Diseases of the digestive system
K65	Peritonitis	Chapter XI	Diseases of the digestive system
K66	Other disorders of peritoneum	Chapter XI	Diseases of the digestive system
K67	Disorders of peritoneum in infectious diseases classified elsewhere	Chapter XI	Diseases of the digestive system
K70	Alcoholic liver disease	Chapter XI	Diseases of the digestive system
K71	Toxic liver disease	Chapter XI	Diseases of the digestive system
K72	Hepatic failure, not elsewhere classified	Chapter XI	Diseases of the digestive system
K73	Chronic hepatitis, not elsewhere classified	Chapter XI	Diseases of the digestive system
K74	Fibrosis and cirrhosis of liver	Chapter XI	Diseases of the digestive system
K75	Other inflammatory liver diseases	Chapter XI	Diseases of the digestive system
K76	Other diseases of liver	Chapter XI	Diseases of the digestive system
K77	Liver disorders in diseases classified elsewhere	Chapter XI	Diseases of the digestive system
K80	Cholelithiasis	Chapter XI	Diseases of the digestive system
K81	Cholecystitis	Chapter XI	Diseases of the digestive system
K82	Other diseases of gallbladder	Chapter XI	Diseases of the digestive system
K83	Other diseases of biliary tract	Chapter XI	Diseases of the digestive system
K85	Acute pancreatitis	Chapter XI	Diseases of the digestive system
K86	Other diseases of pancreas	Chapter XI	Diseases of the digestive system
K87	Disorders of gallbladder, biliary tract and pancreas in diseases classified elsewhere	Chapter XI	Diseases of the digestive system
K90	Intestinal malabsorption	Chapter XI	Diseases of the digestive system
K91	Postprocedural disorders of digestive system, not elsewhere classified	Chapter XI	Diseases of the digestive system
K92	Other diseases of digestive system	Chapter XI	Diseases of the digestive system
K93	Disorders of other digestive organs in diseases classified elsewhere	Chapter XI	Diseases of the digestive system
L00	Staphylococcal scalded skin syndrome	Chapter XII	Diseases of the skin and subcutaneous tissue
L01	Impetigo	Chapter XII	Diseases of the skin and subcutaneous tissue
L02	Cutaneous abscess, furuncle and carbuncle	Chapter XII	Diseases of the skin and subcutaneous tissue
L03	Cellulitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L04	Acute lymphadenitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L05	Pilonidal cyst	Chapter XII	Diseases of the skin and subcutaneous tissue
L08	Other local infections of skin and subcutaneous tissue	Chapter XII	Diseases of the skin and subcutaneous tissue
L10	Pemphigus	Chapter XII	Diseases of the skin and subcutaneous tissue
L11	Other acantholytic disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L12	Pemphigoid	Chapter XII	Diseases of the skin and subcutaneous tissue
L13	Other bullous disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L14	Bullous disorders in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
L20	Atopic dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L21	Seborrhoeic dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L22	Diaper [napkin] dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L23	Allergic contact dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L24	Irritant contact dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L25	Unspecified contact dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L26	Exfoliative dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L27	Dermatitis due to substances taken internally	Chapter XII	Diseases of the skin and subcutaneous tissue
L28	Lichen simplex chronicus and prurigo	Chapter XII	Diseases of the skin and subcutaneous tissue
L29	Pruritus	Chapter XII	Diseases of the skin and subcutaneous tissue
L30	Other dermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L40	Psoriasis	Chapter XII	Diseases of the skin and subcutaneous tissue
L41	Parapsoriasis	Chapter XII	Diseases of the skin and subcutaneous tissue
L42	Pityriasis rosea	Chapter XII	Diseases of the skin and subcutaneous tissue
L43	Lichen planus	Chapter XII	Diseases of the skin and subcutaneous tissue
L44	Other papulosquamous disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L45	Papulosquamous disorders in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
L50	Urticaria	Chapter XII	Diseases of the skin and subcutaneous tissue
L51	Erythema multiforme	Chapter XII	Diseases of the skin and subcutaneous tissue
L52	Erythema nodosum	Chapter XII	Diseases of the skin and subcutaneous tissue
L53	Other erythematous conditions	Chapter XII	Diseases of the skin and subcutaneous tissue
L54	Erythema in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
L55	Sunburn	Chapter XII	Diseases of the skin and subcutaneous tissue
L56	Other acute skin changes due to ultraviolet radiation	Chapter XII	Diseases of the skin and subcutaneous tissue
L57	Skin changes due to chronic exposure to nonionising radiation	Chapter XII	Diseases of the skin and subcutaneous tissue
L58	Radiodermatitis	Chapter XII	Diseases of the skin and subcutaneous tissue
L59	Other disorders of skin and subcutaneous tissue related to radiation	Chapter XII	Diseases of the skin and subcutaneous tissue
L60	Nail disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L62	Nail disorders in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
L63	Alopecia areata	Chapter XII	Diseases of the skin and subcutaneous tissue
L64	Androgenic alopecia	Chapter XII	Diseases of the skin and subcutaneous tissue
L65	Other nonscarring hair loss	Chapter XII	Diseases of the skin and subcutaneous tissue
L66	Cicatricial alopecia [scarring hair loss]	Chapter XII	Diseases of the skin and subcutaneous tissue
L67	Hair colour and hair shaft abnormalities	Chapter XII	Diseases of the skin and subcutaneous tissue
L68	Hypertrichosis	Chapter XII	Diseases of the skin and subcutaneous tissue
L70	Acne	Chapter XII	Diseases of the skin and subcutaneous tissue
L71	Rosacea	Chapter XII	Diseases of the skin and subcutaneous tissue
L72	Follicular cysts of skin and subcutaneous tissue	Chapter XII	Diseases of the skin and subcutaneous tissue
L73	Other follicular disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L74	Eccrine sweat disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L75	Apocrine sweat disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L80	Vitiligo	Chapter XII	Diseases of the skin and subcutaneous tissue
L81	Other disorders of pigmentation	Chapter XII	Diseases of the skin and subcutaneous tissue
L82	Seborrhoeic keratosis	Chapter XII	Diseases of the skin and subcutaneous tissue
L83	Acanthosis nigricans	Chapter XII	Diseases of the skin and subcutaneous tissue
L84	Corns and callosities	Chapter XII	Diseases of the skin and subcutaneous tissue
L85	Other epidermal thickening	Chapter XII	Diseases of the skin and subcutaneous tissue
L86	Keratoderma in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
L87	Transepidermal elimination disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L88	Pyoderma gangrenosum	Chapter XII	Diseases of the skin and subcutaneous tissue
L89	Decubitus ulcer	Chapter XII	Diseases of the skin and subcutaneous tissue
L90	Atrophic disorders of skin	Chapter XII	Diseases of the skin and subcutaneous tissue
L91	Hypertrophic disorders of skin	Chapter XII	Diseases of the skin and subcutaneous tissue
L92	Granulomatous disorders of skin and subcutaneous tissue	Chapter XII	Diseases of the skin and subcutaneous tissue
L93	Lupus erythematosus	Chapter XII	Diseases of the skin and subcutaneous tissue
L94	Other localised connective tissue disorders	Chapter XII	Diseases of the skin and subcutaneous tissue
L95	Vasculitis limited to skin, not elsewhere classified	Chapter XII	Diseases of the skin and subcutaneous tissue
L97	Ulcer of lower limb, not elsewhere classified	Chapter XII	Diseases of the skin and subcutaneous tissue
L98	Other disorders of skin and subcutaneous tissue, not elsewhere classified	Chapter XII	Diseases of the skin and subcutaneous tissue
L99	Other disorders of skin and subcutaneous tissue in diseases classified elsewhere	Chapter XII	Diseases of the skin and subcutaneous tissue
M00	Pyogenic arthritis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M01	Direct infections of joint in infectious and parasitic diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M02	Reactive arthropathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M03	Postinfective and reactive arthropathies in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M05	Seropositive rheumatoid arthritis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M06	Other rheumatoid arthritis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M07	Psoriatic and enteropathic arthropathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M08	Juvenile arthritis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M09	Juvenile arthritis in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M10	Gout	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M11	Other crystal arthropathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M12	Other specific arthropathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M13	Other arthritis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M14	Arthropathies in other diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M15	Polyarthrosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M16	Coxarthrosis [arthrosis of hip]	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M17	Gonarthrosis [arthrosis of knee]	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M18	Arthrosis of first carpometacarpal joint	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M19	Other arthrosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M20	Acquired deformities of fingers and toes	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M21	Other acquired deformities of limbs	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M22	Disorders of patella	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M23	Internal derangement of knee	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M24	Other specific joint derangements	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M25	Other joint disorders, not elsewhere classified	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M30	Polyarteritis nodosa and related conditions	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M31	Other necrotising vasculopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M32	Systemic lupus erythematosus	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M33	Dermatopolymyositis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M34	Systemic sclerosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M35	Other systemic involvement of connective tissue	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M36	Systemic disorders of connective tissue in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M40	Kyphosis and lordosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M41	Scoliosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M42	Spinal osteochondrosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M43	Other deforming dorsopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M45	Ankylosing spondylitis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M46	Other inflammatory spondylopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M47	Spondylosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M48	Other spondylopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M49	Spondylopathies in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M50	Cervical disk disorders	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M51	Other intervertebral disk disorders	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M53	Other dorsopathies, not elsewhere classified	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M54	Dorsalgia	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M60	Myositis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M61	Calcification and ossification of muscle	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M62	Other disorders of muscle	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M63	Disorders of muscle in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M65	Synovitis and tenosynovitis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M66	Spontaneous rupture of synovium and tendon	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M67	Other disorders of synovium and tendon	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M68	Disorders of synovium and tendon in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M70	Soft tissue disorders related to use, overuse and pressure	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M71	Other bursopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M72	Fibroblastic disorders	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M73	Soft tissue disorders in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M75	Shoulder lesions	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M76	Enthesopathies of lower limb, excluding foot	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M77	Other enthesopathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M79	Other soft tissue disorders, not elsewhere classified	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M80	Osteoporosis with pathological fracture	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M81	Osteoporosis without pathological fracture	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M82	Osteoporosis in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M83	Adult osteomalacia	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M84	Disorders of continuity of bone	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M85	Other disorders of bone density and structure	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M86	Osteomyelitis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M87	Osteonecrosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M88	Paget's disease of bone [osteitis deformans]	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M89	Other disorders of bone	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M90	Osteopathies in diseases classified elsewhere	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M91	Juvenile osteochondrosis of hip and pelvis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M92	Other juvenile osteochondrosis	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M93	Other osteochondropathies	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M94	Other disorders of cartilage	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M95	Other acquired deformities of musculoskeletal system and connective tissue	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M96	Postprocedural musculoskeletal disorders, not elsewhere classified	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
M99	Biomechanical lesions, not elsewhere classified	Chapter XIII	Diseases of the musculoskeletal system and connective tissue
N00	Acute nephritic syndrome	Chapter XIV	Diseases of the genitourinary system
N01	Rapidly progressive nephritic syndrome	Chapter XIV	Diseases of the genitourinary system
N02	Recurrent and persistent haematuria	Chapter XIV	Diseases of the genitourinary system
N03	Chronic nephritic syndrome	Chapter XIV	Diseases of the genitourinary system
N04	Nephrotic syndrome	Chapter XIV	Diseases of the genitourinary system
N05	Unspecified nephritic syndrome	Chapter XIV	Diseases of the genitourinary system
N06	Isolated proteinuria with specified morphological lesion	Chapter XIV	Diseases of the genitourinary system
N07	Hereditary nephropathy, not elsewhere classified	Chapter XIV	Diseases of the genitourinary system
N08	Glomerular disorders in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N10	Acute tubulo-interstitial nephritis	Chapter XIV	Diseases of the genitourinary system
N11	Chronic tubulo-interstitial nephritis	Chapter XIV	Diseases of the genitourinary system
N12	Tubulo-interstitial nephritis, not specified as acute or chronic	Chapter XIV	Diseases of the genitourinary system
N13	Obstructive and reflux uropathy	Chapter XIV	Diseases of the genitourinary system
N14	Drug- and heavy-metal-induced tubulo-interstitial and tubular conditions	Chapter XIV	Diseases of the genitourinary system
N15	Other renal tubulo-interstitial diseases	Chapter XIV	Diseases of the genitourinary system
N16	Renal tubulo-interstitial disorders in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N17	Acute renal failure	Chapter XIV	Diseases of the genitourinary system
N18	Chronic renal failure	Chapter XIV	Diseases of the genitourinary system
N19	Unspecified renal failure	Chapter XIV	Diseases of the genitourinary system
N20	Calculus of kidney and ureter	Chapter XIV	Diseases of the genitourinary system
N21	Calculus of lower urinary tract	Chapter XIV	Diseases of the genitourinary system
N22	Calculus of urinary tract in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N23	Unspecified renal colic	Chapter XIV	Diseases of the genitourinary system
N25	Disorders resulting from impaired renal tubular function	Chapter XIV	Diseases of the genitourinary system
N26	Unspecified contracted kidney	Chapter XIV	Diseases of the genitourinary system
N27	Small kidney of unknown cause	Chapter XIV	Diseases of the genitourinary system
N28	Other disorders of kidney and ureter, not elsewhere classified	Chapter XIV	Diseases of the genitourinary system
N29	Other disorders of kidney and ureter in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N30	Cystitis	Chapter XIV	Diseases of the genitourinary system
N31	Neuromuscular dysfunction of bladder, not elsewhere classified	Chapter XIV	Diseases of the genitourinary system
N32	Other disorders of bladder	Chapter XIV	Diseases of the genitourinary system
N33	Bladder disorders in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N34	Urethritis and urethral syndrome	Chapter XIV	Diseases of the genitourinary system
N35	Urethral stricture	Chapter XIV	Diseases of the genitourinary system
N36	Other disorders of urethra	Chapter XIV	Diseases of the genitourinary system
N37	Urethral disorders in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N39	Other disorders of urinary system	Chapter XIV	Diseases of the genitourinary system
N40	Hyperplasia of prostate	Chapter XIV	Diseases of the genitourinary system
N41	Inflammatory diseases of prostate	Chapter XIV	Diseases of the genitourinary system
N42	Other disorders of prostate	Chapter XIV	Diseases of the genitourinary system
N43	Hydrocele and spermatocele	Chapter XIV	Diseases of the genitourinary system
N44	Torsion of testis	Chapter XIV	Diseases of the genitourinary system
N45	Orchitis and epididymitis	Chapter XIV	Diseases of the genitourinary system
N46	Male infertility	Chapter XIV	Diseases of the genitourinary system
N47	Redundant prepuce, phimosis and paraphimosis	Chapter XIV	Diseases of the genitourinary system
N48	Other disorders of penis	Chapter XIV	Diseases of the genitourinary system
N49	Inflammatory disorders of male genital organs, not elsewhere classified	Chapter XIV	Diseases of the genitourinary system
N50	Other disorders of male genital organs	Chapter XIV	Diseases of the genitourinary system
N51	Disorders of male genital organs in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N60	Benign mammary dysplasia	Chapter XIV	Diseases of the genitourinary system
N61	Inflammatory disorders of breast	Chapter XIV	Diseases of the genitourinary system
N62	Hypertrophy of breast	Chapter XIV	Diseases of the genitourinary system
N63	Unspecified lump in breast	Chapter XIV	Diseases of the genitourinary system
N64	Other disorders of breast	Chapter XIV	Diseases of the genitourinary system
N70	Salpingitis and oophoritis	Chapter XIV	Diseases of the genitourinary system
N71	Inflammatory disease of uterus, except cervix	Chapter XIV	Diseases of the genitourinary system
N72	Inflammatory disease of cervix uteri	Chapter XIV	Diseases of the genitourinary system
N73	Other female pelvic inflammatory diseases	Chapter XIV	Diseases of the genitourinary system
N74	Female pelvic inflammatory disorders in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N75	Diseases of Bartholin's gland	Chapter XIV	Diseases of the genitourinary system
N76	Other inflammation of vagina and vulva	Chapter XIV	Diseases of the genitourinary system
N77	Vulvovaginal ulceration and inflammation in diseases classified elsewhere	Chapter XIV	Diseases of the genitourinary system
N80	Endometriosis	Chapter XIV	Diseases of the genitourinary system
N81	Female genital prolapse	Chapter XIV	Diseases of the genitourinary system
N82	Fistulae involving female genital tract	Chapter XIV	Diseases of the genitourinary system
N83	Noninflammatory disorders of ovary, Fallopian tube and broad ligament	Chapter XIV	Diseases of the genitourinary system
N84	Polyp of female genital tract	Chapter XIV	Diseases of the genitourinary system
N85	Other noninflammatory disorders of uterus, except cervix	Chapter XIV	Diseases of the genitourinary system
N86	Erosion and ectropion of cervix uteri	Chapter XIV	Diseases of the genitourinary system
N87	Dysplasia of cervix uteri	Chapter XIV	Diseases of the genitourinary system
N88	Other noninflammatory disorders of cervix uteri	Chapter XIV	Diseases of the genitourinary system
N89	Other noninflammatory disorders of vagina	Chapter XIV	Diseases of the genitourinary system
N90	Other noninflammatory disorders of vulva and perineum	Chapter XIV	Diseases of the genitourinary system
N91	Absent, scanty and rare menstruation	Chapter XIV	Diseases of the genitourinary system
N92	Excessive, frequent and irregular menstruation	Chapter XIV	Diseases of the genitourinary system
N93	Other abnormal uterine and vaginal bleeding	Chapter XIV	Diseases of the genitourinary system
N94	Pain and other conditions associated with female genital organs and menstrual cycle	Chapter XIV	Diseases of the genitourinary system
N95	Menopausal and other perimenopausal disorders	Chapter XIV	Diseases of the genitourinary system
N96	Habitual aborter	Chapter XIV	Diseases of the genitourinary system
N97	Female infertility	Chapter XIV	Diseases of the genitourinary system
N98	Complications associated with artificial fertilisation	Chapter XIV	Diseases of the genitourinary system
N99	Postprocedural disorders of genito-urinary system, not elsewhere classified	Chapter XIV	Diseases of the genitourinary system",
    sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) != 1081L || anyDuplicated(d$code)) stop("Embedded ICD dictionary integrity check failed.", call. =
    FALSE)
  attr(d, "source") <- list(url = "https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=19", download_url = "https://biobank.ndph.ox.ac.uk/ukb/codown.cgi?id=19",
    retrieved = "2026-09-09", csv_sha256 = "3e1b7af7e326552c4a83c3ff62310a36e750453e5b802e3a43526e7ba8c02ad8",
    transformation = "Three-character A-N codes; chapter traced from parent nodes; leading code removed from official disease name.")
  d
}
disease_chapter_metadata <- function() {
  roman <- c(as.character(utils::as.roman(1:14)), "?")
  short <- c("Infections", "Neoplasms", "Blood / immune", "Endocrine", "Mental", "Nervous", "Eye", "Ear",
    "Circulatory", "Respiratory", "Digestive", "Skin", "Musculoskeletal", "Genitourinary", "Unmapped")
  data.frame(chapter = c(paste("Chapter", roman[1:14]), "Unmapped"), roman = roman, short_label = short, label =
    paste(roman, short), order = seq_along(roman), stringsAsFactors = FALSE)
}
prevalent_disease_composition <- function(I, D = disease_chapter_metadata()) {
  need_i <- c("icd3", "gated", "chapter");
  need_d <- c("chapter", "label", "order")
  if (!all(need_i %in% names(I)) || !all(need_d %in% names(D))) stop("Prevalent composition requires I: icd3/gated/chapter and D: chapter/label/order.",
    call. = FALSE)
  I <- as.data.frame(I, stringsAsFactors = FALSE);
  D <- as.data.frame(D, stringsAsFactors = FALSE)
  I$chapter <- as.character(I$chapter);
  D$chapter <- as.character(D$chapter)
  if (anyDuplicated(I$icd3) || anyDuplicated(D$chapter) || anyNA(D$chapter)) stop("Prevalent composition requires unique ICD codes and unique non-missing chapter keys.",
    call. = FALSE)
  if (anyNA(I$gated) || !is.logical(I$gated)) stop("Prevalent composition requires a complete logical integrated gate.",
    call. = FALSE)
  if (anyNA(match(I$chapter, D$chapter))) stop("Prevalent composition chapter dictionary does not cover all eligible diseases.",
    call. = FALSE)
  D <- D[order(D$order, D$chapter), , drop = FALSE]
  idx <- match(I$chapter, D$chapter);
  N <- sum(I$gated)
  D$n_eligible <- tabulate(idx, nbins = nrow(D))
  D$n_significant <- tabulate(idx[I$gated], nbins = nrow(D))
  D$n_significant_total <- N
  D$share_percent <- if (N > 0L) 100 * D$n_significant / N else rep(NA_real_, nrow(D))
  D$rank <- NA_integer_
  j <- which(D$n_significant > 0L)
  j <- j[order( - D$n_significant[j], D$order[j], D$chapter[j])]
  D$rank[j] <- seq_along(j)
  D
}
summarize_prevalent_diseases <- function(results, bundle, disease_n = 20L, selected_diseases = NULL) {
  fail <- function(...) stop(..., call. = FALSE)
  tab <- function(x) as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  required <- function(x, nm, where) {
    if (!all(nm %in% names(x))) fail(where, " lacks: ", paste(setdiff(nm, names(x)), collapse = ", "));
    invisible(x)
  }
  if (!is.list(results) || is.null(results$results)) fail("Supply the complete in-memory clinical analysis result.")
  required(results, c("clock_map", "endpoint_counts", "results"), "Complete result")
  config <- results$config;
  if (is.null(config)) config <- list()
  alpha <- config$alpha
  if (is.null(alpha)) {
    alpha <- .05;
    config$alpha_source <- "Default ALPHA=0.05 in the supplied analysis; saved alpha was unavailable."
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) fail("Invalid saved alpha.")
  config$alpha <- alpha
  if (!is.numeric(disease_n) || length(disease_n) != 1L || is.na(disease_n) || disease_n < 1 || (!is.infinite(disease_n) &&
    disease_n != floor(disease_n))) fail("disease_n must be a positive integer or Inf.")
  cm <- tab(results$clock_map);
  counts <- tab(results$endpoint_counts)
  required(cm, c("exposure", "clock", "layer"), "clock_map");
  required(counts, c("icd3", "eligible_prevalent"), "endpoint_counts")
  for (v in c("exposure", "clock", "layer")) cm[[v]] <- as.character(cm[[v]])
  if (anyNA(cm[, c("exposure", "clock", "layer")]) || anyDuplicated(cm$exposure) || anyDuplicated(cm$clock) ||
    nrow(cm) != 59L || sum(cm$layer == "Integrated") != 1L || sum(cm$layer == "Pathway") != 58L) fail("clock_map must contain one integrated clock and 58 unique pathways.")
  counts$icd3 <- as.character(counts$icd3)
  if (anyNA(counts$icd3) || anyDuplicated(counts$icd3) || any(!grepl("^[A-Z][0-9]{2}$", counts$icd3))) fail("Invalid or duplicate ICD-10 three-character codes in endpoint_counts.")
  eligible <- counts$eligible_prevalent
  if (!is.logical(eligible)) {
    zz <- toupper(as.character(eligible));
    if (any(!is.na(zz) & !zz %in% c("TRUE", "FALSE", "1", "0"))) fail("eligible_prevalent must be logical.")
    eligible <- ifelse(is.na(zz), NA, zz %in% c("TRUE", "1"))
  }
  codes <- sort(counts$icd3[eligible %in% TRUE]);
  M <- length(codes)
  if (!M) fail("No eligible prevalent endpoints are recorded in endpoint_counts.")
  int_exposure <- cm$exposure[cm$layer == "Integrated"];
  pm <- cm[cm$layer == "Pathway", , drop = FALSE]
  ann <- classify_reactome8(bundle, pm$clock);
  ann$exposure <- pm$exposure[match(ann$module, pm$clock)]
  classes <- reactome8_classes()
  numbers <- c("N", "Events", "Beta", "SE", "Effect", "Lower95", "Upper95", "P", "P_Bonf", "m_tests")
  fields <- c("icd3", "exposure", "clock", "layer", "analysis", "measure", "family", "status", numbers)
  canonical <- function() {
    z <- as.data.frame(stats::setNames(replicate(length(fields), character(), simplify = FALSE), fields),
      stringsAsFactors = FALSE);
    for (v in numbers) z[[v]] <- numeric();
    z
  }
  validate <- function(z, family, mplan, exposures) {
    z <- tab(z);
    if (!nrow(z)) return(canonical())
    required(z, fields, family)
    for (v in setdiff(fields, numbers)) z[[v]] <- as.character(z[[v]])
    for (v in numbers) {
      old <- z[[v]];
      z[[v]] <- suppressWarnings(as.numeric(as.character(old)))
      if (any(is.na(z[[v]]) & !is.na(old) & nzchar(as.character(old)))) fail(family, " has non-numeric ",
        v, " values.")
    }
    if (anyNA(z[, c("icd3", "exposure", "clock", "layer", "analysis", "measure", "family", "status")]) ||
      anyDuplicated(z[, c("icd3", "exposure")])) fail("Missing identifiers/status or duplicate endpoint-exposure rows in ",
      family, ".")
    if (any(!z$icd3 %in% codes) || any(!z$exposure %in% exposures)) fail("Unplanned endpoint/exposure in ",
      family, ".")
    if (any(z$analysis != "Prevalent") || any(z$measure != "OR") || any(z$family != family) || any(!is.finite(z$m_tests) |
      z$m_tests != mplan)) fail("Wrong analysis, measure, family or fixed m_tests in ", family, ".")
    k <- match(z$exposure, cm$exposure)
    if (any(z$clock != cm$clock[k]) || any(z$layer != cm$layer[k])) fail("Clock identity does not match clock_map in ",
      family, ".")
    if (any(is.finite(z$P) & (z$P < 0 | z$P > 1)) || any(is.finite(z$P_Bonf) & (z$P_Bonf < 0 | z$P_Bonf >
      1))) fail("P values outside [0,1] in ", family, ".")
    if (any(is.finite(z$P) != is.finite(z$P_Bonf))) fail("Raw and Bonferroni P availability differs in ",
      family, ".")
    ok <- is.finite(z$P);
    expected <- pmin(1, z$P[ok] * mplan)
    if (any(abs(z$P_Bonf[ok] - expected) > 1e-300 + 1e-8 * abs(expected))) fail("Saved P_Bonf disagrees with the fixed planned family in ",
      family, "; no P values were replaced.")
    z
  }
  pi <- validate(results$results$prev_int, "prevalent_integrated", M, int_exposure)
  pp <- validate(results$results$prev_path, "prevalent_pathway_gated", 58L * M, pm$exposure)
  attach_rows <- function(keys, z) {
    kk <- function(d) paste(d$icd3, d$exposure, sep = "\r")
    j <- match(kk(keys), kk(z));
    out <- z[j, , drop = FALSE];
    rownames(out) <- NULL
    out$icd3 <- keys$icd3;
    out$exposure <- keys$exposure;
    out$row_present <- !is.na(j)
    q <- match(out$exposure, cm$exposure);
    out$clock <- cm$clock[q];
    out$layer <- cm$layer[q];
    out
  }
  evaluate <- function(z) {
    n <- nrow(z);
    z$estimate_status <- ifelse(!z$row_present, "missing_row", ifelse(z$status == "warning_review", "warning_review",
      ifelse(z$status %in% c("ok", "warning"), "invalid_estimate", "error")))
    z$log2OR <- z$log2Lower <- z$log2Upper <- rep(NA_real_, n);
    z$ci_recovered <- rep(FALSE, n)
    good <- is.finite(z$Beta) & is.finite(z$SE) & z$SE > 0 & is.finite(z$Beta - 1.96 * z$SE) & is.finite(z$Beta +
      1.96 * z$SE) & is.finite(z$P) & is.finite(z$P_Bonf)

    # The supplied GLM writes exp(Beta +/- 1.96*SE); verify saved values before summarizing.
    available <- rep(TRUE, n)
    for (item in c("Effect", "Lower95", "Upper95")) {
      target <- if (item == "Effect") z$Beta else z$Beta + if (item == "Lower95") - 1.96 * z$SE else 1.96 *
      z$SE
      actual <- z[[item]];
      direct <- good & is.finite(actual) & actual > 0
      mismatch <- direct & abs(suppressWarnings(log(actual)) - target) > 1e-7 + 1e-7 * abs(target)
      if (any(mismatch)) fail("Saved ", item, " disagrees with Beta/SE for ", paste(utils::head(paste(z$icd3[mismatch],
        z$exposure[mismatch]), 5), collapse = ", "), ".")
      implied <- suppressWarnings(exp(target));
      ranged <- good & !is.na(actual) & ((actual == 0 & implied == 0) | (actual == Inf & implied == Inf))
      ranged[is.na(ranged)] <- FALSE;
      available <- available & (direct | ranged)
      z$ci_recovered <- z$ci_recovered | ranged
    }
    usable <- good & available & z$row_present & z$status %in% c("ok", "warning")
    z$estimate_status[usable] <- "valid"
    z$log2OR[usable] <- z$Beta[usable] / log(2)
    z$log2Lower[usable] <- (z$Beta[usable] - 1.96 * z$SE[usable]) / log(2)
    z$log2Upper[usable] <- (z$Beta[usable] + 1.96 * z$SE[usable]) / log(2)
    z$model_warning <- z$row_present & z$status %in% "warning"
    z$significant <- usable & is.finite(z$P_Bonf) & z$P_Bonf < alpha
    z
  }
  I <- evaluate(attach_rows(data.frame(icd3 = codes, exposure = int_exposure, stringsAsFactors = FALSE), pi))
  I$gated <- I$row_present & I$status %in% c("ok", "warning") & is.finite(I$P_Bonf) & I$P_Bonf < alpha
  gate <- I[I$gated, , drop = FALSE];
  gate <- gate[order(gate$P_Bonf, gate$P, gate$icd3), , drop = FALSE]
  if (any(!pp$icd3 %in% gate$icd3)) fail("Pathway results contain disease(s) outside the saved integrated significance gate.")
  if (is.null(selected_diseases)) chosen <- utils::head(gate$icd3, if (is.infinite(disease_n)) nrow(gate) else disease_n) else {
    if (!is.character(selected_diseases) || !length(selected_diseases) || anyNA(selected_diseases) || anyDuplicated(selected_diseases) ||
      any(!selected_diseases %in% gate$icd3)) fail("selected_diseases must be nonempty unique ICD-10 codes from the integrated significance gate.")
    chosen <- selected_diseases
  }
  S <- I[match(chosen, I$icd3), , drop = FALSE];
  S$disease_index <- seq_len(nrow(S));
  S$selection_rank <- match(S$icd3, gate$icd3)
  keys <- expand.grid(exposure = pm$exposure, icd3 = chosen, stringsAsFactors = FALSE)
  H <- evaluate(attach_rows(keys, pp));
  k <- match(H$exposure, ann$exposure)
  for (v in setdiff(names(ann), "exposure")) H[[v]] <- ann[[v]][k]
  H$disease_index <- match(H$icd3, S$icd3)
  if (nrow(H) != 58L * length(chosen) || anyDuplicated(H[, c("icd3", "exposure")])) fail("Internal pathway grid integrity check failed.")
  diagnostics <- list(alpha = alpha, n_eligible = M, integrated_family = M, pathway_family = 58L * M, n_gated =
    nrow(gate), n_selected = length(chosen), gated_codes = gate$icd3, selected_codes = chosen, omitted_gated_codes =
    setdiff(gate$icd3, chosen), empty_gate = !nrow(gate), integrated_status = table(I$estimate_status), pathway_status =
    table(H$estimate_status), missing_integrated_rows = sum(!I$row_present), missing_selected_pathway_rows =
    sum(!H$row_present), missing_all_gated_pathway_rows = 58L * nrow(gate) - nrow(pp), numeric_range_recoveries =
    sum(I$ci_recovered) + sum(H$ci_recovered), selection = "Lowest saved integrated P_Bonf, then raw P, then ICD-10; explicit selected_diseases retains its supplied order.",
    gate = "Saved integrated status ok/warning and P_Bonf < saved alpha; no pathway-significance selection.",
    effect = "Adjusted prevalent-disease odds ratio per 1 SD full-OOF age acceleration; saved logistic models only.",
    missingness = "Missing rows, review warnings and failed/non-estimable results never receive OR=1 or a significance mark.",
    classification = attr(ann, "classification_note"))
  dictionary <- disease_icd_dictionary();
  j <- match(I$icd3, dictionary$code)
  I$disease_name <- dictionary$disease_name[j];
  I$chapter <- dictionary$chapter[j]
  missing_names <- is.na(I$disease_name) | !nzchar(I$disease_name)
  I$disease_name[missing_names] <- I$icd3[missing_names];
  I$chapter[is.na(I$chapter)] <- "Unmapped"
  list(composition = prevalent_disease_composition(I), integrated = I, grid = H, pathways = ann, classes =
    classes, selected_df = S, clock_map = cm, endpoint_counts = counts, diagnostics = diagnostics, inputconfig =
    config)
}
disease_summary_require <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(d)
}
disease_summary_names <- function(d, labels = NULL) {
  d <- mortality_disease_labels(d, labels)
  dictionary <- disease_icd_dictionary();
  j <- match(d$icd3, dictionary$code)
  d[, `:=`(disease = dictionary$disease_name[j], chapter_label = dictionary$chapter_name[j], chapter_order =
    match(dictionary$chapter[j], disease_chapter_metadata()$chapter))]
  if (!is.null(labels)) {
    disease_summary_require(labels, c("icd3", "disease"), "Disease labels")
    if (anyDuplicated(labels$icd3)) stop("Duplicate ICD3 in disease labels.")
    j <- match(d$icd3, toupper(trimws(as.character(labels$icd3))))
    ok <- !is.na(j) & !is.na(labels$disease[j]) & nzchar(trimws(labels$disease[j]))
    data.table::set(d, i = which(ok), j = "disease", value = as.character(labels$disease[j[ok]]))
  }
  if (anyNA(d$disease) || any(!nzchar(trimws(d$disease)))) stop("Missing official disease names.")
  d
}
disease_pathway_metadata <- function(cm, bundle) {
  cm <- data.table::copy(data.table::as.data.table(cm))[layer == "Pathway"]
  ann <- classify_reactome8(bundle, cm$clock)
  id <- if ("reactome_id" %in% names(ann)) ann$reactome_id else ann$id
  nm <- if ("pathway_name" %in% names(ann)) ann$pathway_name else ann$pathway
  j <- match(cm$clock, ann$module)
  cm[, `:=`(reactome_id = id[j], full_name = nm[j], group = reactome8_classes()$class_name[match(ann$primary_class[j],
    reactome8_classes()$class_id)])]
  cm
}
incident_disease_composition <- function(a) {
  disease_summary_require(a, c('icd3', 'chapter_short', 'sig', 'eligible'), 'Integrated disease composition');
  a <- data.table::copy(data.table::as.data.table(a));
  if (anyNA(a$sig) || anyNA(a$chapter_short) || anyDuplicated(a$icd3)) stop('Composition requires one named, non-missing significance record per disease');
  if (any(!a$chapter_short %in% mortality_disease_chapters)) stop('Unmapped disease category in composition');
  if (any(a$sig & !a$eligible %in% TRUE)) stop('Significant composition records need available finite log HRs')
  ch <- data.table::data.table(chapter_short = union(mortality_disease_chapters[seq_len(14)], unique(a$chapter_short)));
  cnt <- a[, .(tested = .N, significant = sum(sig)), by = chapter_short];
  ch <- merge(ch, cnt, by = 'chapter_short', all.x = TRUE, sort = FALSE);
  ch[is.na(tested), `:=`(tested = 0L, significant = 0L)];
  total <- sum(ch$significant);
  stopifnot(total == sum(a$sig));
  ch[, `:=`(percent = if (total) 100 * significant / total else 0, within_chapter_percent = ifelse(tested >
    0, 100 * significant / tested, NA_real_))];
  data.table::setorder(ch, - significant, chapter_short);
  ch[, `:=`(rank = seq_len(.N))]
  members <- a[sig == TRUE, intersect(c('icd3', 'disease', 'chapter_short', 'Beta', 'Effect', 'P_Bonf', 'status',
    'warning', 'error', 'has_warning', 'eligible'), names(a)), with = FALSE];
  data.table::set(members, j = 'category_rank', value = match(members$chapter_short, ch$chapter_short));
  data.table::setorderv(members, c('category_rank', 'icd3'));
  stopifnot(nrow(members) == total, identical(sort(members$icd3), sort(a[sig == TRUE, icd3])))
  list(data = ch, total = total, members = members)
}
incident_landmark_concordance <- function(a, l) {
  need <- c('icd3', 'Beta', 'P_Bonf', 'status', 'trusted', 'review', 'eligible', 'has_warning', 'sig', 'disease',
    'chapter_short');
  for (x in list(a, l)) {
    disease_summary_require(x, need, 'Integrated concordance table');
    if (anyDuplicated(x$icd3)) stop('Duplicate ICD3 in integrated concordance table')
  }
  suffix <- function(x, s) {
    x <- data.table::copy(data.table::as.data.table(x));
    data.table::set(x, j = 'has_row', value = rep(TRUE, nrow(x)));
    n <- setdiff(names(x), 'icd3');
    data.table::setnames(x, n, paste0(n, '_', s));
    x
  };
  d <- merge(suffix(a, 'main'), suffix(l, 'landmark'), by = 'icd3', all = TRUE, sort = FALSE)
  take <- function(n) {
    v <- d[[paste0(n, '_main')]];
    z <- d[[paste0(n, '_landmark')]];
    j <- is.na(v);
    if (is.character(v)) j <- j | !nzchar(v);
    v[j] <- z[j];
    v
  };
  for (n in c('disease', 'chapter_short')) data.table::set(d, j = n, value = take(n))
  d[, `:=`(has_main = has_row_main %in% TRUE, has_landmark = has_row_landmark %in% TRUE, paired_finite =
    is.finite(Beta_main) & is.finite(Beta_landmark))];
  d[, `:=`(paired_trusted = paired_finite & trusted_main %in% TRUE & trusted_landmark %in% TRUE, paired_eligible =
    paired_finite & eligible_main %in% TRUE & eligible_landmark %in% TRUE, paired_review = paired_finite &
    (!(trusted_main %in% TRUE) | !(trusted_landmark %in% TRUE)), absolute_logHR_difference = ifelse(paired_finite,
    abs(Beta_main - Beta_landmark), NA_real_))]
  d[, `:=`(eligible_agreement = paired_eligible & sig_main %in% TRUE, omitted_reason = ifelse(!has_main, 'Main analysis missing',
    ifelse(!has_landmark, 'Landmark analysis missing', ifelse(!paired_finite, 'Non-finite paired estimate',
    ifelse(!paired_eligible, 'Failed or unavailable model', '')))))]
  eligible <- data.table::copy(d[eligible_agreement == TRUE]);
  data.table::setorder(eligible, absolute_logHR_difference, icd3)
  eligible[, agreement_rank := seq_len(.N)]
  d[, agreement_rank := eligible$agreement_rank[match(icd3, eligible$icd3)]];
  data.table::setorder(d, icd3)
  list(data = d, eligible = eligible)
}
disease_pathway_input_counts <- function(immlist, cm, bundle) {
  if (!is.list(immlist) || is.data.frame(immlist) || is.null(names(immlist)) || anyNA(names(immlist)) ||
    any(!nzchar(names(immlist))) || anyDuplicated(names(immlist))) stop('Supply the named input-protein list: immlist.')
  id <- vapply(names(immlist), function(s) {
    m <- regexpr('R-HSA-[0-9]+', s);
    if (m[1] < 0) NA_character_ else regmatches(s, m)
  }, character(1));
  pm <- disease_pathway_metadata(cm, bundle);
  needed <- pm$reactome_id
  if (anyDuplicated(id[!is.na(id) & id %in% needed])) stop('immlist contains duplicate Reactome pathway IDs');
  j <- match(needed, id);
  if (anyNA(j)) stop('Input proteins missing for: ', paste(needed[is.na(j)], collapse = ', '))
  z <- data.table::rbindlist(lapply(seq_along(j), function(k) {
    v <- immlist[[j[k]]];
    if (!is.character(v)) stop('Protein sets must contain character IDs');
    trim <- trimws(v);
    valid <- !is.na(trim) & nzchar(trim);
    clean <- unique(trim[valid]);
    data.table::data.table(exposure = pm$exposure[k], reactome_id = needed[k], protein_n = length(clean),
      input_entries = length(v), excluded_blank_or_na = sum(!valid), duplicate_entries = sum(valid) - length(clean),
      input_name = names(immlist)[j[k]], protein_count_source = 'Unique input-protein identifiers in supplied immlist')
  }))
  if (any(z$protein_n < 1L)) stop('At least one pathway has no valid input protein');
  if (any(z$excluded_blank_or_na + z$duplicate_entries > 0)) warning('NA, blank or duplicate protein identifiers excluded from size counts; inspect the returned input-protein counts.',
    call. = FALSE);
  z
}
incident_pathway_breadth <- function(h, cm, a, protein_counts, bundle) {
  pm <- disease_pathway_metadata(cm, bundle);
  if (nrow(pm) != 58L || anyDuplicated(pm$reactome_id)) stop('The clock map must contain all 58 distinct Reactome pathways');
  pm <- merge(pm, protein_counts[, .(exposure, reactome_id, protein_n, protein_count_source)], by = c('exposure',
    'reactome_id'), all.x = TRUE, sort = FALSE);
  if (anyNA(pm$protein_n)) stop('Missing pathway input-protein counts')
  target <- sort(unique(as.character(a[sig == TRUE, icd3])))
  source_scope <- sort(unique(as.character(h$icd3)))
  if (any(!source_scope %in% a$icd3)) stop('Pathway results contain diseases absent from integrated results')
  if (any(!h$exposure %in% pm$exposure)) stop('Unexpected exposure in incident pathway table')
  scope <- intersect(target, source_scope)
  coverage <- data.table::copy(a[sig == TRUE])
  observed <- h[, .(n_pathway_rows = .N, n_eligible_pathways = sum(eligible), n_significant_pathways = sum(sig)),
    by = icd3]
  coverage <- merge(coverage, observed, by = 'icd3', all.x = TRUE, sort = FALSE)
  for (v in c('n_pathway_rows', 'n_eligible_pathways', 'n_significant_pathways')) coverage[is.na(get(v)),
    (v) := 0L]
  coverage[, `:=`(expected_pathways = nrow(pm), has_saved_pathway_results = n_pathway_rows > 0L, complete_eligible_pathways =
    n_eligible_pathways == nrow(pm), missing_reason = ifelse(n_pathway_rows == 0L, 'No saved pathway analysis for this integrated-significant disease',
    ifelse(n_eligible_pathways < nrow(pm), 'One or more pathway estimates unavailable', '')))]
  category_coverage <- coverage[, .(integrated_significant_diseases = .N, diseases_with_saved_pathway_results =
    sum(has_saved_pathway_results), diseases_with_all_58_eligible_pathways = sum(complete_eligible_pathways),
    diseases_without_saved_pathway_results = sum(!has_saved_pathway_results)), by = chapter_short]
  if (!length(target)) return(list(data = data.table::data.table(), scope = scope, grid = data.table::data.table(),
    coverage = coverage, category_coverage = category_coverage))
  grid <- merge(data.table::CJ(exposure = pm$exposure, icd3 = target), h, by = c('exposure', 'icd3'), all.x =
    TRUE, sort = FALSE)
  for (v in c('trusted', 'eligible', 'effect_eligible', 'has_warning', 'sig', 'review')) grid[is.na(get(v)),
    (v) := FALSE]
  grid[, has_row := !is.na(status)]
  j <- match(grid$icd3, a$icd3)
  for (v in c('disease', 'chapter_short', 'chapter_label')) if (v %in% names(a)) data.table::set(grid, j =
    v, value = a[[v]][j])
  grid[, `:=`(available_model = eligible, omitted_reason = ifelse(!has_row, 'No saved pathway analysis', ifelse(!eligible,
    'Failed or non-finite estimate', '')))]
  b <- grid[, .(n_expected = .N, n_rows = sum(has_row), n_trusted = sum(trusted), n_eligible = sum(eligible),
    n_warning_included = sum(eligible & has_warning), n_review = sum(review), n_missing = sum(!has_row), n_failed =
    sum(has_row & !eligible), n_significant = sum(sig), n_categories = data.table::uniqueN(chapter_short[sig ==
    TRUE]), n_positive = sum(sig & Beta > 0, na.rm = TRUE), n_negative = sum(sig & Beta < 0, na.rm = TRUE)),
    by = exposure]
  b <- merge(pm, b, by = 'exposure', all.x = TRUE, sort = FALSE)
  b[, `:=`(n_integrated_significant_diseases = length(target), n_diseases_with_saved_pathway_results = length(scope),
    n_integrated_significant_without_pathway_results = length(setdiff(target, scope)), coverage = ifelse(n_eligible ==
    n_expected, 'Complete', 'Incomplete'), available = n_eligible > 0L)]
  data.table::setorder(b, n_significant, n_categories, exposure)
  b[, omitted_reason := ifelse(available, '', 'No available pathway-disease model')]
  list(data = b, scope = scope, grid = grid, coverage = coverage, category_coverage = category_coverage)
}
summarize_incident_diseases <- function(clinical, immlist, bundle, disease_labels = NULL) {
  a <- disease_summary_names(mortality_profile_eligibility(clinical$results$inc_int), disease_labels)
  l <- disease_summary_names(mortality_profile_eligibility(clinical$results$inc_lag2), disease_labels)
  h <- mortality_profile_eligibility(clinical$results$inc_path)
  protein_counts <- disease_pathway_input_counts(immlist, clinical$clock_map, bundle)
  list(composition = incident_disease_composition(a), landmark_concordance = incident_landmark_concordance(a,
    l), pathway_breadth = incident_pathway_breadth(h, clinical$clock_map, a, protein_counts, bundle), input_protein_counts =
    protein_counts)
}
disease_summary_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
disease_summary_numeric <- function(d, columns, label) {
  disease_summary_assert(all(vapply(d[columns], is.numeric, logical(1))), paste(label, "requires numeric columns:",
    paste(columns, collapse = ", ")))
}
summarize_disease_category_counts <- function(comparison, alpha = .05, type = "disease") {
  cats <- unique(comparison$category)
  counts <- do.call(rbind, lapply(cats, function(cat) {
    z <- comparison[comparison$category == cat, , drop = FALSE]
    incomplete <- sum(!z$valid_I)
    observed <- sum(z$immune_positive_significant)
    data.frame(category = cat, candidate_diseases = nrow(z), evaluable_immune = sum(z$valid_I), missing_immune =
      incomplete, significant_positive_immune = if (incomplete) NA_integer_ else observed, observed_positive_immune =
      observed, stringsAsFactors = FALSE)
  }))
  ok <- counts$missing_immune == 0L
  dd <- counts[ok, , drop = FALSE]
  estimable <- nrow(dd) >= 3L && length(unique(dd$candidate_diseases)) > 1L && length(unique(dd$significant_positive_immune)) >
  1L
  pv <- if (estimable) unname(stats::cor.test(dd$candidate_diseases, dd$significant_positive_immune, method =
    'spearman', alternative = 'two.sided', exact = FALSE)$p.value) else NA_real_
  status <- if (!estimable) 'Not estimable: fewer than three complete categories or constant counts' else if (any(!ok)) 'Complete categories only; incomplete immune categories excluded' else 'All categories complete'
  test <- data.frame(n_categories = nrow(dd), n_categories_total = nrow(counts), p_value = pv, method = 'Spearman correlation; asymptotic P',
    alternative = 'two.sided', status = status, significance = 'Saved immune clock Bonferroni P < alpha and effect > 1',
    alpha = alpha, interpretation = 'Descriptive category-count association; does not establish enrichment or causation')
  if (any(!ok)) message(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)), ': ', sum(!ok), ' categories have incomplete immune inference; shown as NA and omitted from the correlation test.')
  list(data = counts, test = test)
}
summarize_disease_comparator <- function(associations, cfg = list(alpha = .05), type = c("prevalent", "incident")) {
  type <- match.arg(type)
  req <- c("code", "disease", "category", "clock", "effect", "lower", "upper", "p_bonf", "n", "cases")
  disease_summary_require(associations, req, paste(type, "disease associations"))
  d <- as.data.frame(associations, stringsAsFactors = FALSE)
  if ("type" %in% names(d)) d <- d[!is.na(d$type) & d$type == type, , drop = FALSE]
  disease_summary_assert(nrow(d) > 0L, paste("No", type, "disease association rows were supplied."))
  for (nm in c("code", "disease", "category", "clock")) {
    d[[nm]] <- as.character(d[[nm]])
    disease_summary_assert(all(!is.na(d[[nm]]) & nzchar(trimws(d[[nm]]))), paste("Missing", nm, "in", type,
      "disease associations."))
  }
  disease_summary_assert(all(d$clock %in% c("I", "G")), "Disease clocks must be coded I or G.")
  disease_summary_assert(!anyDuplicated(d[c("code", "clock")]), "Duplicate disease code/clock rows: supply one estimate per frozen code and clock.")
  disease_summary_numeric(d, c("effect", "lower", "upper", "p_bonf", "n", "cases"), "Disease associations")
  disease_summary_assert(all(is.finite(d$n) & d$n > 0 & d$n == floor(d$n) & is.finite(d$cases) & d$cases >=
    0 & d$cases <= d$n & d$cases == floor(d$cases)), "Disease n/cases must be valid participant counts.")
  alpha <- if (is.null(cfg$alpha)) .05 else cfg$alpha
  disease_summary_assert(length(alpha) == 1L && is.finite(alpha) && alpha > 0 && alpha < 1, "cfg$alpha must lie strictly between 0 and 1.")
  cfg$alpha <- alpha
  codes <- sort(unique(d$code))
  expected <- cfg[[if (type == "prevalent") "expected_prev" else "expected_inc"]]
  if (!is.null(expected)) {
    disease_summary_assert(length(expected) == 1L && is.finite(expected) && expected > 0 && expected == floor(expected),
      "Expected frozen disease count must be a positive integer or NULL.")
    disease_summary_assert(length(codes) == expected, sprintf("%s frozen family contains %d unique codes; expected %d. Supply the complete saved family, not a significance-filtered table.",
      type, length(codes), expected))
  }
  if ("m_tests" %in% names(d)) {
    disease_summary_numeric(d, "m_tests", "Disease associations")
    disease_summary_assert(all(is.finite(d$m_tests) & d$m_tests == length(codes)), "Saved m_tests does not match the complete frozen disease family; do not recompute P values after filtering.")
  }

  # Join by the stable code, never by a potentially duplicated disease label.
  ii <- match(codes, d$code[d$clock == "I"])
  gg <- match(codes, d$code[d$clock == "G"])
  di <- d[d$clock == "I", , drop = FALSE][ii, , drop = FALSE]
  dg <- d[d$clock == "G", , drop = FALSE][gg, , drop = FALSE]
  paired <- !is.na(ii) & !is.na(gg)
  for (nm in c("disease", "category", "n", "cases")) {
    same <- di[[nm]][paired] == dg[[nm]][paired]
    disease_summary_assert(all(!is.na(same) & same), paste0("Paired clocks disagree on ", nm, " for a disease code; the same endpoint and participants are required."))
  }
  valid <- function(x) {
    is.finite(x$effect) & is.finite(x$lower) & is.finite(x$upper) & x$lower > 0 & x$lower <= x$effect & x$effect <=
    x$upper & is.finite(x$p_bonf) & x$p_bonf >= 0 & x$p_bonf <= 1
  }
  vi <- !is.na(ii) & valid(di)
  vg <- !is.na(gg) & valid(dg)
  meta <- d[match(codes, d$code), c("code", "disease", "category"), drop = FALSE]
  comparison <- meta
  for (nm in c("effect", "lower", "upper", "p_bonf", "n", "cases")) {
    comparison[[paste0(nm, "_I")]] <- di[[nm]]
    comparison[[paste0(nm, "_G")]] <- dg[[nm]]
  }
  comparison$has_I <- !is.na(ii)
  comparison$has_G <- !is.na(gg)
  comparison$valid_I <- vi
  comparison$valid_G <- vg
  comparison$immune_significant <- vi & di$p_bonf < alpha
  comparison$all_protein_nonsignificant <- vg & dg$p_bonf >= alpha
  comparison$immune_positive_significant <- comparison$immune_significant & di$effect > 1
  comparison$selected <- comparison$immune_positive_significant & comparison$all_protein_nonsignificant

  # Missing or invalid G inference is never treated as a negative comparison.
  comparison$reason <- vapply(seq_along(codes), function(k) {
    reason <- character()
    if (!comparison$has_I[k]) reason <- c(reason, "missing_I")
    else if (!vi[k]) reason <- c(reason, "invalid_I_estimate_CI_or_P")
    if (!comparison$has_G[k]) reason <- c(reason, "missing_G")
    else if (!vg[k]) reason <- c(reason, "invalid_G_estimate_CI_or_P")
    if (length(reason)) return(paste(reason, collapse = ";"))
    if (comparison$selected[k]) return("selected")
    if (di$p_bonf[k] >= alpha) reason <- c(reason, "I_not_significant")
    if (di$effect[k] <= 1) reason <- c(reason, "I_not_positive")
    if (dg$p_bonf[k] < alpha) reason <- c(reason, "G_significant")
    paste(reason, collapse = ";")
  }, character(1))
  category_order <- cfg$disease_category_order
  if (is.null(category_order)) category_order <- sort(unique(comparison$category))
  category_order <- unique(c(as.character(category_order), sort(unique(comparison$category))))
  comparison <- comparison[order(match(comparison$category, category_order), comparison$p_bonf_I, comparison$code,
    na.last = TRUE), , drop = FALSE]
  rownames(comparison) <- NULL
  chosen <- comparison[comparison$selected, , drop = FALSE]
  selected <- d[d$code %in% chosen$code, , drop = FALSE]
  selected <- selected[order(match(selected$code, chosen$code), match(selected$clock, c("I", "G"))), , drop =
    FALSE]
  rownames(selected) <- NULL
  list(comparison = comparison, selected = selected, category_counts = summarize_disease_category_counts(comparison,
    alpha, type))
}


# ==============================================================================
# 06 | Mortality splines, absolute risks and pathway profiles
# ==============================================================================

# Mortality associations, standardized risks, and disease association profiles.
# The mortality analysis accepts the prepared clinical baseline and the
# integrated mortality result. Calculations return R objects in memory.
# Participant bootstrap confidence intervals follow the manuscript analysis.
mortality_covariates <- c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol",
  "physical_activity", "bmi")
mortality_capture <- function(fun) {
  w <- character();
  err <- ""
  z <- tryCatch(withCallingHandlers(fun(), warning = function(e) {
    w <<- c(w, conditionMessage(e));
    invokeRestart("muffleWarning")
  }), error = function(e) {
    err <<- conditionMessage(e);
    NULL
  })
  w <- paste(unique(w), collapse = " | ");
  status <- if (nzchar(err)) "error" else if (!nzchar(w)) "ok" else if (grepl("did not converge|ran out of iterations|may be infinite|separation|singular|not positive|non-finite",
    w, ignore.case = TRUE)) "warning_review" else "warning"
  list(value = z, status = status, warning = w, error = err)
}
mortality_spline_basis <- function(x, k) {
  if (length(k) != 4L || any(!is.finite(k)) || any(diff(k) <= 0)) stop("Four distinct finite spline knots are required.")
  h <- matrix(vapply(k[1:2], function(t)(pmax(x - t, 0) ^ 3 - pmax(x - k[3], 0) ^ 3 * (k[4] - t) / (k[4] -
    k[3]) + pmax(x - k[4], 0) ^ 3 * (k[3] - t) / (k[4] - k[3])) / (k[4] - k[1]) ^ 2, numeric(length(x))),
    nrow = length(x), ncol = 2L)
  z <- cbind(rcs1 = x, h);
  colnames(z) <- c("rcs1", "rcs2", "rcs3");
  z
}
fit_mortality_ph <- function(q, knots, kind = "spline") {
  terms <- mortality_covariates
  if (kind == "linear") terms <- c("exposure", terms)
  if (kind == "spline") {
    q[c("rcs1", "rcs2", "rcs3")] <- as.data.frame(mortality_spline_basis(q$exposure, knots));
    terms <- c("rcs1", "rcs2", "rcs3", terms)
  }
  form <- stats::reformulate(terms, response = "survival::Surv(time,outcome)")
  survival::coxph(form, data = q, ties = "efron", x = TRUE, y = TRUE, model = FALSE, na.action = stats::na.fail,
    singular.ok = TRUE)
}
standardize_mortality_ph_risk <- function(fit, q, knots, levels, times) {
  bn <- c("rcs1", "rcs2", "rcs3");
  b <- stats::coef(fit)[bn];
  if (any(!is.finite(b))) stop("Spline exposure coefficients are not estimable.")
  bh <- survival::basehaz(fit, centered = TRUE);
  if (!nrow(bh) || any(!is.finite(bh$hazard)) || any(diff(bh$hazard) < - 1e-10)) stop("Invalid cumulative baseline hazard.")
  h <- c(0, bh$hazard)[findInterval(times, bh$time) + 1L];
  lp0 <- as.numeric(fit$linear.predictors - mortality_spline_basis(q$exposure, knots) %*% b)
  if (any(!is.finite(lp0))) stop("Non-finite centered linear predictors.")
  out <- vapply(levels, function(z) {
    rr <- exp(lp0 + as.numeric(mortality_spline_basis(z, knots) %*% b));
    vapply(h, function(a) if (a == 0) 0 else mean( - expm1( - a * rr)), numeric(1))
  }, numeric(length(times)))
  if (any(!is.finite(out)) || any(out < - 1e-12 | out > 1 + 1e-12)) stop("Non-finite marginal standardized risks.")
  out
}
mortality_ph_once <- function(q, knots, levels, times) {
  fit <- fit_mortality_ph(q, knots);
  standardize_mortality_ph_risk(fit, q, knots, levels, times)
}
mortality_ph_bootstrap_worker <- function(i) {
  w <- get(".mortality_ph_work", envir = .GlobalEnv)
  assign(".Random.seed", w$streams[[i]], envir = .GlobalEnv)
  ans <- mortality_capture(function() {
    ii <- sample.int(nrow(w$q), nrow(w$q), replace = TRUE)
    mortality_ph_once(w$q[ii, , drop = FALSE], w$knots, w$levels, w$times)
  })
  ans$replicate <- as.integer(i)
  ans
}
run_mortality_ph_extensions <- function(base, death_censor, reference_result, cores = 30L, boot = 200L, seed =
  20260910L, horizon = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival before running mortality analysis.")
  for (nm in c("cores", "boot", "seed")) {
    v <- get(nm);
    if (length(v) != 1L || !is.finite(v) || v != as.integer(v)) stop(nm, " must be an integer.")
  }
  if (cores < 1L || cores > 30L || boot < 0L || seed < 1L || !is.finite(horizon) || length(horizon) != 1L ||
    horizon <= 0) stop("Use cores=1..30, boot>=0, positive seed and horizon.")
  cc <- tryCatch(as.integer(as.Date(death_censor)), error = function(e) NA_integer_);
  if (length(cc) != 1L || is.na(cc)) stop("Supply the exact death censor date used in the completed analysis.")
  need <- c("eid", "Z000", "blood_day", "death_day", mortality_covariates);
  if (length(m <- setdiff(need, names(base)))) stop("Missing baseline columns: ", paste(m, collapse = ", "))
  d <- as.data.frame(base)[as.data.frame(base)$blood_day < cc, , drop = FALSE];
  b <- d$blood_day;
  dd <- d$death_day;
  event <- !is.na(dd) & dd > b & dd <= cc;
  follow <- (ifelse(event, dd, cc) - b) / 365.25;
  k <- is.finite(follow) & follow > 0
  q <- data.frame(time = follow[k], outcome = as.integer(event[k]), exposure = d$Z000[k], d[k, mortality_covariates,
    drop = FALSE], check.names = FALSE)
  if (!nrow(q) || anyNA(q) || any(!is.finite(q$exposure)) || sum(q$outcome) < 20L || anyDuplicated(d$eid[k])) stop("Invalid mortality cohort: check complete cases, EIDs and death dates.")
  rr <- as.data.frame(reference_result);
  if ("exposure" %in% names(rr)) rr <- rr[rr$exposure == "Z000", , drop = FALSE];
  if ("landmark_years" %in% names(rr)) rr <- rr[rr$landmark_years == 0, , drop = FALSE]
  if (nrow(rr) != 1L || !all(c("N", "Events", "Beta") %in% names(rr))) stop("reference_result must contain one original integrated mortality result.")
  refnum <- as.numeric(unlist(rr[1, c("N", "Events", "Beta")], use.names = FALSE));
  if (any(!is.finite(refnum))) stop("The original integrated mortality result is not finite; review that model first.")
  if (nrow(q) != rr$N || sum(q$outcome) != rr$Events) stop(sprintf("Mortality cohort mismatch: reconstructed N=%d/deaths=%d; original N=%d/deaths=%d. Reuse the exact OOF data, covariates and censor date.",
    nrow(q), sum(q$outcome), rr$N, rr$Events))
  envnames <- c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OPENBLAS_DEFAULT_NUM_THREADS",
    "GOTO_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS");
  oldenv <- Sys.getenv(envnames, unset = NA_character_);
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(envnames)), envnames)))
  on.exit({
    isna <- is.na(oldenv);
    if (any(isna)) Sys.unsetenv(envnames[isna]);
    if (any(!isna)) do.call(Sys.setenv, as.list(oldenv[!isna]))
  }, add = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    oldthreads <- data.table::getDTthreads();
    data.table::setDTthreads(1L);
    on.exit(data.table::setDTthreads(oldthreads), add = TRUE)
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    oldblas <- RhpcBLASctl::blas_get_num_procs();
    RhpcBLASctl::blas_set_num_threads(1L);
    on.exit(RhpcBLASctl::blas_set_num_threads(oldblas), add = TRUE)
  }
  knots <- as.numeric(stats::quantile(q$exposure, c(.05, .35, .65, .95), type = 7));
  if (any(diff(knots) <= 0)) stop("Repeated exposure knots: the four-knot spline cannot be fitted.")
  linear <- mortality_capture(function() fit_mortality_ph(q, knots, "linear"));
  if (is.null(linear$value)) stop("Mortality validation model failed: ", linear$error)
  beta <- unname(stats::coef(linear$value)["exposure"]);
  if (!is.finite(beta) || abs(beta - rr$Beta) > 1e-6 * max(1, abs(rr$Beta))) stop(sprintf("Mortality coefficient mismatch: reconstructed Beta=%.9g; original Beta=%.9g. Do not use mismatched data or restandardize Z000.",
    beta, rr$Beta))
  support_n <- min(nrow(q), max(20L, ceiling(.05 * nrow(q))));
  supported <- sort(q$time, decreasing = TRUE)[support_n];
  used_horizon <- min(horizon, supported);
  if (used_horizon <= 0) stop("No supported follow-up horizon.")
  times <- seq(0, used_horizon, length.out = 81L);
  levels <- as.numeric(stats::quantile(q$exposure, c(.1, .5, .9), type = 7));
  names(levels) <- c("P10", "P50", "P90");
  zgrid <- sort(unique(c(seq(stats::quantile(q$exposure, .01), stats::quantile(q$exposure, .99), length.out =
    201L), 0)))
  main <- mortality_capture(function() {
    fit <- fit_mortality_ph(q, knots);
    b <- stats::coef(fit);
    nn <- c("rcs1", "rcs2", "rcs3");
    if (any(!is.finite(b[nn]))) stop("Spline coefficients are not estimable.")
    V <- stats::vcov(fit)[nn, nn, drop = FALSE];
    X <- sweep(mortality_spline_basis(zgrid, knots), 2, mortality_spline_basis(0, knots)[1, ], "-");
    eta <- as.numeric(X %*% b[nn]);
    se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
    if (any(!is.finite(c(eta, se)))) stop("Spline HR or uncertainty is non-finite.")
    nul <- fit_mortality_ph(q, knots, "null");
    LL <- c(nul$loglik[2], linear$value$loglik[2], fit$loglik[2]);
    DF <- c(sum(is.finite(stats::coef(nul))), sum(is.finite(stats::coef(linear$value))), sum(is.finite(b)))
    stat <- pmax(0, 2 * (LL[3] - LL[1:2]));
    df <- DF[3] - DF[1:2];
    if (any(df <= 0)) stop("Spline likelihood-ratio models are not properly nested.")
    tests <- data.frame(term = c("Overall association", "Nonlinearity"), statistic = stat, df = df, P = stats::pchisq(stat,
      df, lower.tail = FALSE));
    tests$P_Bonf <- stats::p.adjust(tests$P, "bonferroni", n = 2L);
    tests$family <- "New mortality spline: overall and nonlinearity";
    tests$m_tests <- 2L
    list(curve = data.frame(z = zgrid, HR = exp(eta), Lower95 = exp(eta - 1.96 * se), Upper95 = exp(eta +
      1.96 * se)), risk = standardize_mortality_ph_risk(fit, q, knots, levels, times), tests = tests)
  })
  if (is.null(main$value)) stop("Spline/marginal-risk computation failed: ", main$error, if (nzchar(main$warning)) paste0("; warnings: ",
    main$warning))
  logrow <- function(x, stage, replicate = NA_integer_) data.frame(stage = stage, replicate = replicate, status =
    x$status, warning = x$warning, error = x$error, stringsAsFactors = FALSE)
  logs <- rbind(logrow(linear, "original mortality model validation"), logrow(main, "spline and marginal risks"))
  if (used_horizon < horizon) logs <- rbind(logs, data.frame(stage = "follow-up support", replicate = NA_integer_,
    status = "warning", warning = sprintf("Requested %.2f years shortened to %.2f years to retain at least %d at risk.",
    horizon, used_horizon, support_n), error = ""))
  risk <- data.frame(time = rep(times, 3L), level = rep(names(levels), each = length(times)), z = rep(levels,
    each = length(times)), risk = as.numeric(main$value$risk), Lower95 = NA_real_, Upper95 = NA_real_)
  rngkind <- RNGkind();
  hadseed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE);
  if (hadseed) oldseed <- get(".Random.seed", envir = .GlobalEnv);
  on.exit({
    do.call(RNGkind, as.list(rngkind));
    if (hadseed) assign(".Random.seed", oldseed, envir = .GlobalEnv) else if (exists(".Random.seed", envir =
      .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG");
  set.seed(seed);
  streams <- vector("list", boot);
  s <- get(".Random.seed", envir = .GlobalEnv);
  if (boot) for (i in seq_len(boot)) {
    streams[[i]] <- s;
    s <- parallel::nextRNGStream(s)
  }
  work <- list(q = q, knots = knots, levels = levels, times = times, streams = streams);
  replicates <- vector("list", boot)
  if (boot) {
    ncore <- min(as.integer(cores), boot);
    message("Mortality cohort verified; ", boot, " participant bootstraps, ", ncore, " processes, one BLAS thread each.")
    cl <- parallel::makePSOCKcluster(ncore);
    on.exit(parallel::stopCluster(cl), add = TRUE)
    exports <- c("mortality_covariates", "mortality_capture", "mortality_spline_basis", "fit_mortality_ph",
      "standardize_mortality_ph_risk", "mortality_ph_once", "mortality_ph_bootstrap_worker");
    parallel::clusterExport(cl, exports, envir = environment(run_mortality_ph_extensions));
    parallel::clusterCall(cl, function(w, libs) {
      .libPaths(libs);
      assign(".mortality_ph_work", w, envir = .GlobalEnv);
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L);
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L);
      NULL
    }, work, .libPaths())
    batches <- split(seq_len(boot), ceiling(seq_len(boot) / max(1L, ncore * 2L)));
    for (ii in batches) {
      replicates[ii] <- parallel::parLapplyLB(cl, ii, mortality_ph_bootstrap_worker);
      message("Mortality bootstrap: ", max(ii), "/", boot, " completed.")
    }
    logs <- rbind(logs, do.call(rbind, lapply(replicates, function(x) logrow(x, "participant bootstrap", x$replicate))))
  }
  good <- if (boot) vapply(replicates, function(x) x$status %in% c("ok", "warning") && is.matrix(x$value) &&
    all(is.finite(x$value)), logical(1)) else logical();
  n_success <- sum(good);
  needed <- max(100L, ceiling(.8 * boot));
  ci_ok <- boot > 0L && n_success >= needed;
  bootvalues <- if (n_success) do.call(cbind, lapply(replicates[good], function(x) as.numeric(x$value))) else matrix(numeric(),
    nrow(risk), 0L)
  if (ci_ok) {
    ci <- t(apply(bootvalues, 1, stats::quantile, probs = c(.025, .975), names = FALSE, type = 7));
    risk$Lower95 <- ci[, 1];
    risk$Upper95 <- ci[, 2]
  } else logs <- rbind(logs, data.frame(stage = "bootstrap intervals", replicate = NA_integer_, status =
    "warning", warning = if (boot == 0L) "boot=0: point estimates only; marginal-risk confidence intervals were not computed." else sprintf("Only %d eligible bootstrap fits; %d required. All replicate results retained; risk confidence intervals omitted.",
    n_success, needed), error = ""))
  summary <- risk[risk$time == max(times), , drop = FALSE];
  summary$horizon <- used_horizon
  diffpoint <- main$value$risk[nrow(main$value$risk), 3] - main$value$risk[nrow(main$value$risk), 1];
  diffci <- c(NA_real_, NA_real_)
  if (ci_ok) diffci <- stats::quantile(bootvalues[nrow(risk), ] - bootvalues[length(times), ], c(.025, .975),
    names = FALSE, type = 7)
  contrast <- data.frame(comparison = "P90 minus P10", horizon = used_horizon, risk_difference = diffpoint,
    Lower95 = diffci[1], Upper95 = diffci[2])
  atrisk_times <- sort(unique(c(0, pretty(c(0, used_horizon), n = 5L), used_horizon)));
  atrisk_times <- atrisk_times[atrisk_times >= 0 & atrisk_times <= used_horizon];
  at_risk <- data.frame(time = atrisk_times, N = vapply(atrisk_times, function(t) sum(q$time >= t), integer(1)))
  manifest <- list(version = "mortality_ph_v1", N = nrow(q), Events = sum(q$outcome), death_censor = as.Date(cc,
    origin = "1970-01-01"), covariates = mortality_covariates, validation = list(original_beta = rr$Beta,
    reconstructed_beta = beta), knots = knots, knot_quantiles = c(.05, .35, .65, .95), exposure_levels =
    levels, HR_reference = 0, requested_horizon = horizon, horizon = used_horizon, minimum_at_risk = support_n,
    boot_requested = boot, boot_success = n_success, boot_required = needed, boot_CI_available = ci_ok, seed =
    seed, cores = cores, main_status = main$status, standardization = "Exact average over observed eligible cohort covariates; each bootstrap refits and averages over its resampled participants; spline knots and P10/P50/P90 fixed.",
    uncertainty = "Spline HR: pointwise coefficient Wald intervals. Absolute risks: participant bootstrap percentile intervals, conditional on the already trained OOF clocks and pooled acceleration transformation.",
    preparation = attr(base, "mortality_acceleration"), audit = attr(base, "mortality_audit"))
  list(curve = main$value$curve, risk = risk, risk_summary = summary, risk_contrast = contrast, at_risk =
    at_risk, tests = main$value$tests, logs = logs, manifest = manifest, exposure = data.frame(z = q$exposure),
    bootstrap_results = replicates)
}

# Disease-mortality similarity across the common eligible pathway profiles.
mortality_profile_eligibility <- function(d) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table.")
  d <- data.table::copy(data.table::as.data.table(d))
  needed <- c("exposure", "icd3", "Beta", "Effect", "Lower95", "Upper95", "P_Bonf", "status")
  missing <- setdiff(needed, names(d))
  if (length(missing)) stop("Model table missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(d[, c("icd3", "exposure"), with = FALSE])) stop("Duplicate disease/exposure models.")
  for (nm in c("Beta", "Effect", "Lower95", "Upper95", "P_Bonf")) if (!is.numeric(d[[nm]])) stop("Expected numeric model column: ",
    nm)
  if (any(is.finite(d$P_Bonf) & (d$P_Bonf < 0 | d$P_Bonf > 1))) stop("Invalid P_Bonf outside [0,1].")
  for (nm in c("warning", "error")) if (!nm %in% names(d)) data.table::set(d, j = nm, value = rep("", nrow(d)))
  status <- as.character(d$status)
  warning_text <- as.character(d$warning)
  error_text <- as.character(d$error)
  has_error <- !is.na(error_text) & nzchar(trimws(error_text)) & !trimws(error_text) %in% c("NA", "<NA>")
  accepted <- status %in% c("ok", "warning", "warning_review")
  eligible <- accepted & !has_error & is.finite(d$Beta)
  effect_ok <- eligible & is.finite(d$Effect) & d$Effect > 0
  ci_ok <- effect_ok & is.finite(d$Lower95) & is.finite(d$Upper95) & d$Lower95 > 0 & d$Upper95 >= d$Lower95 &
  d$Lower95 <= d$Effect & d$Upper95 >= d$Effect
  has_warning <- status %in% c("warning", "warning_review") | (!is.na(warning_text) & nzchar(trimws(warning_text)) &
    !trimws(warning_text) %in% c("NA", "<NA>"))
  reason <- ifelse(!accepted, paste0("model_status: ", status), ifelse(has_error, "recorded_error", ifelse(!is.finite(d$Beta),
    "nonfinite_Beta", "")))
  data.table::set(d, j = "trusted", value = status %in% c("ok", "warning") & is.finite(d$Beta))
  data.table::set(d, j = "review", value = status == "warning_review" & is.finite(d$Beta))
  data.table::set(d, j = "eligible", value = eligible)
  data.table::set(d, j = "effect_eligible", value = effect_ok)
  data.table::set(d, j = "ci_valid", value = ci_ok)
  data.table::set(d, j = "has_warning", value = has_warning)
  data.table::set(d, j = "retained_warning", value = has_warning & eligible)
  data.table::set(d, j = "exclusion_reason", value = reason)
  data.table::set(d, j = "sig", value = eligible & is.finite(d$P_Bonf) & d$P_Bonf < .05)
  d
}
mortality_disease_chapters <- c("Infection", "Neoplasms", "Blood/immune", "Endocrine/metabolic", "Mental/behavioural",
  "Nervous", "Eye", "Ear", "Circulatory", "Respiratory", "Digestive", "Skin", "Musculoskeletal", "Genitourinary",
  "Other")
mortality_disease_labels <- function(d, labels = NULL) {
  d <- data.table::copy(data.table::as.data.table(d))
  code <- toupper(as.character(d$icd3));
  letter <- substr(code, 1, 1)
  n <- suppressWarnings(as.integer(substr(code, 2, 3)))
  ch <- rep("Other", length(code))
  lut <- c(A = "Infection", B = "Infection", C = "Neoplasms", E = "Endocrine/metabolic", F = "Mental/behavioural",
    G = "Nervous", I = "Circulatory", J = "Respiratory", K = "Digestive", L = "Skin", M = "Musculoskeletal",
    N = "Genitourinary")
  j <- match(letter, names(lut));
  ch[!is.na(j)] <- unname(lut[j[!is.na(j)]])
  ch[which(letter == "D" & n <= 48)] <- "Neoplasms"
  ch[which(letter == "D" & n >= 50 & n <= 89)] <- "Blood/immune"
  ch[which(letter == "H" & n <= 59)] <- "Eye"
  ch[which(letter == "H" & n >= 60 & n <= 95)] <- "Ear"
  d[, `:=`(chapter = ch, chapter_short = ch, disease = code)]
  if (!is.null(labels)) {
    if (!all(c("icd3", "disease") %in% names(labels))) stop("Disease labels require icd3 and disease.")
    if (anyDuplicated(labels$icd3)) stop("Duplicate ICD3 in disease labels.")
    j <- match(code, toupper(as.character(labels$icd3)))
    good <- !is.na(j) & !is.na(labels$disease[j]) & nzchar(as.character(labels$disease[j]))
    data.table::set(d, i = which(good), j = "disease", value = as.character(labels$disease[j[good]]))
  }
  d
}
select_mortality_profile_diseases <- function(d, codes = NULL, n = 30L) {
  z <- mortality_disease_labels(mortality_profile_eligibility(d));
  z[, chapter_order := match(chapter_short, mortality_disease_chapters)];
  if (!is.null(codes)) {
    codes <- unique(toupper(trimws(as.character(codes))));
    if (!length(codes) || any(!codes %in% z$icd3)) stop('selected_diseases contains empty or unknown ICD3 codes');
    s <- z[match(codes, icd3)];
    if (any(!s$sig)) stop('Disease selection must contain eligible Bonferroni-significant diseases');
    s[, selection_reason := 'User-specified significant disease']
  } else {
    z <- z[sig == TRUE];
    if (!nrow(z)) {
      z[, `:=`(selection_reason = character(), selection_rank = integer())];
      return(z)
    };
    data.table::setorder(z, P_Bonf, icd3);
    first <- z[!duplicated(chapter_short)];
    first[, selection_reason := 'Chapter coverage: smallest original Bonferroni P'];
    if (nrow(first) > n) first <- first[seq_len(n)];
    extra <- head(z[!icd3 %in% first$icd3], max(0L, n - nrow(first)));
    extra[, selection_reason := 'Remaining places: smallest original Bonferroni P'];
    s <- data.table::rbindlist(list(first, extra), fill = TRUE)
  };
  data.table::setorder(s, chapter_order, P_Bonf, icd3);
  s[, selection_rank := seq_len(.N)];
  s
}
run_disease_mortality_profiles <- function(inc_path, mort_path, selected, disease_labels = NULL) {
  i <- mortality_profile_eligibility(inc_path);
  m <- mortality_profile_eligibility(mort_path);
  sel <- unique(toupper(trimws(as.character(selected$icd3))));
  i <- i[icd3 %in% sel]

  # Completely unavailable disease profiles are reported, without suppressing other diseases.
  available_sel <- sel[vapply(sel, function(cd) any(i$icd3 == cd), logical(1))]
  common_for <- function(flag) {
    if (!length(available_sel)) return(character())
    sets <- lapply(available_sel, function(cd) i$exposure[i$icd3 == cd & i[[flag]]])
    sort(Reduce(intersect, c(list(m$exposure[m[[flag]]]), sets)))
  }
  common <- common_for('eligible')
  profile <- function(cd, ids) {
    a <- i[icd3 == cd];
    x <- a$Beta[match(ids, a$exposure)];
    y <- m$Beta[match(ids, m$exposure)]
    nn <- sum(is.finite(x) & is.finite(y))
    st <- if (!nrow(a)) 'no_pathway_estimates_for_disease' else if (length(ids) < 10L) 'insufficient_common_pathways' else if (nn !=
      length(ids)) 'nonfinite_profile' else if (stats::sd(x) <= 0 || stats::sd(y) <= 0) 'zero_variance_profile' else 'Exploratory; no P value'
    list(n = if (nrow(a)) length(ids) else NA_integer_, n_matched = nn, rho = if (identical(st, 'Exploratory; no P value')) stats::cor(x,
      y, method = 'spearman') else NA_real_, status = st)
  }
  ans <- data.table::rbindlist(lapply(sel, function(cd) {
    a <- i[icd3 == cd];
    z <- profile(cd, common)
    data.table::data.table(icd3 = cd, n_pathways = z$n, n_matched_pathways = z$n_matched, n_incident_eligible =
      sum(a$eligible), n_mortality_eligible = sum(m$eligible), n_incident_retained_warning = sum(a$has_warning &
      a$eligible & a$exposure %in% common), n_mortality_retained_warning = sum(m$has_warning & m$eligible &
      m$exposure %in% common), rho = z$rho, status = z$status)
  }), fill = TRUE)
  if (!length(sel)) ans <- data.table::data.table(icd3 = character(), n_pathways = integer(), n_matched_pathways =
    integer(), n_incident_eligible = integer(), n_mortality_eligible = integer(), n_incident_retained_warning =
    integer(), n_mortality_retained_warning = integer(), rho = numeric(), status = character())
  ans <- mortality_disease_labels(ans, disease_labels);
  data.table::setorder(ans, rho, icd3)
  ans[, `:=`(comparison = "Incident disease vs all-cause mortality", profile = "Pathway log hazard ratio per SD acceleration")]
  source_rows <- data.table::rbindlist(list(data.table::copy(i)[, source := 'incident_pathway'], data.table::copy(m)[,
    source := 'mortality_pathway']), fill = TRUE)
  source_rows[, `:=`(used_in_common_profile = eligible & exposure %in% common, retained_warning = eligible &
    has_warning & exposure %in% common, omission_reason = ifelse(!is.finite(Beta), 'nonfinite_coefficient',
    ifelse(!eligible, exclusion_reason, ifelse(!exposure %in% common, 'not_shared_across_all_available_selected_diseases',
    ''))))]
  list(data = ans, common_pathways = common, source_models = source_rows)
}


# ==============================================================================
# 07 | All-protein comparator
# ==============================================================================

# All-protein comparator: manuscript Methods and supplementary Figures S3-S5.
# Protein values and median-imputed covariates are prepared upstream. This module
# retains the supplied model settings and performs no additional imputation.
all_protein_config <- function() {
  list(seed = 666L, folds = 5L, cores = 4L, learner_threads = 1L, association_cores = 8L, expected_proteins =
    2920L, expected_prev = 117L, expected_inc = 274L, min_cases = 200L, alpha = .05, disease_censor = "2022-05-31",
    death_censor = "2024-08-31", covariates = c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index",
    "smoking", "alcohol", "physical_activity", "bmi"), imputation = "none", preprocessing_provenance = "original supplied pre-imputed NPX matrix",
    covariate_missing = "complete_case", allow_failed_candidates = FALSE, fit_timeout = c(train = 3600, predict =
    Inf), mortality_knots_probs = c(.05, .35, .65, .95), mortality_curve_probs = c(.01, .99))
}

# Shared analytical utilities
sf_assert <- function(test, msg) {
  if (length(test) != 1L || is.na(test) || !test) stop(msg, call. = FALSE)
  invisible(TRUE)
}
sf_require <- function(pkgs) {
  bad <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(bad)) stop("Missing R packages: ", paste(bad, collapse = ", "), ". Install the required packages before running this analysis.",
    call. = FALSE)
}
sf_hash <- function(object) {
  sf_require("digest")
  digest::digest(object, algo = "xxhash64", serialize = TRUE)
}
sf_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w));
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  err <- if (inherits(value, "error")) conditionMessage(value) else ""
  list(value = if (nzchar(err)) NULL else value, warnings = unique(warnings), error = err, status = if (nzchar(err)) "error" else if (length(warnings)) "warning" else "ok")
}
`%||%` <- function(x, y) if (is.null(x)) y else x
sf_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- lapply(x, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[cols]
  })
  ans <- do.call(rbind, out);
  rownames(ans) <- NULL;
  ans
}
sf_metrics <- function(y, p) {
  sf_assert(length(y) == length(p) && length(y) > 1L && all(is.finite(y)) && all(is.finite(p)), "Metrics require aligned finite paired observations.")
  c(MAE = mean(abs(y - p)), RMSE = sqrt(mean((y - p) ^ 2)), PCC = if (stats::sd(y) > 0 && stats::sd(p) >
    0) stats::cor(y, p) else NA_real_)
}
sf_acceleration <- function(pred, age) {
  sf_assert(length(pred) == length(age) && all(is.finite(pred)) && all(is.finite(age)), "Age correction requires aligned finite values.")
  f <- stats::lm.fit(cbind(1, age), pred);
  r <- as.numeric(f$residuals);
  s <- stats::sd(r)
  sf_assert(is.finite(s) && s > 0 && all(is.finite(f$coefficients)), "Age acceleration is not estimable.")
  list(residual = r, z = (r - mean(r)) / s, intercept = unname(f$coefficients[1]), slope = unname(f$coefficients[2]),
    sd = s)
}
sf_cache <- function(key, signature, compute, cfg) {
  sf_assert(is.function(compute), "Compute argument must be a function.")
  compute()
}
sf_log <- function(cfg, context, capture) {
  if (!length(capture$warnings) && !nzchar(capture$error %||% "")) return(invisible(NULL))
  invisible(data.frame(context = as.character(context), status = capture$status, warning = paste(capture$warnings,
    collapse = " | "), error = capture$error %||% "", stringsAsFactors = FALSE))
}

# Prepared input objects and original learner specifications
sf_ids <- function(x, where) {
  z <- as.character(x)
  sf_assert(length(z) > 0L && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z), paste(where, "requires unique non-missing participant IDs."))
  z
}
sf_date_days <- function(x, where) {
  if (inherits(x, "Date")) return(as.integer(x))
  if (inherits(x, "POSIXt")) return(as.integer(as.Date(x)))
  if (is.numeric(x)) {
    sf_assert(all(is.na(x) | (x > - 40000 & x < 80000 & x == floor(x))), paste(where, "numeric dates must be days since1970-01-01."))
    return(as.integer(x))
  }
  q <- trimws(as.character(x));
  q[q %in% c("", "NA")] <- NA_character_
  out <- as.Date(q, format = "%Y-%m-%d")
  sf_assert(!any(!is.na(q) & is.na(out)), paste(where, "contains invalid ISO dates."))
  as.integer(out)
}
sf_configure_learner <- function(learner, cfg, saved_timeout = NULL) {
  sf_assert(inherits(learner, "Learner"), "Saved model specification did not produce an mlr3 Learner.")
  sf_assert(is.function(learner$encapsulate), "mlr3 >= 0.21.0 is required for the encapsulate() API.")

  # A real fallback object is required by mlr3 even when errors must propagate.
  # It is never accepted as an age model: sf_check_learner rejects all use.
  fb <- mlr3::lrn("regr.featureless")
  args <- list(method = "callr", fallback = fb)
  if ("when" %in% names(formals(learner$encapsulate))) args$when <- function(cond, stage) FALSE
  do.call(learner$encapsulate, args)
  timeout <- cfg$fit_timeout %||% saved_timeout %||% c(train = 600, predict = 600)
  if (is.numeric(timeout) && length(timeout) == 1L) timeout <- c(train = unname(timeout), predict = unname(timeout))
  sf_assert(is.numeric(timeout) && length(timeout) == 2L && setequal(names(timeout), c("train", "predict")) &&
    all(!is.na(timeout) & timeout > 0), "Learner timeout must be positive seconds, or named c(train=...,predict=...).")
  learner$timeout <- timeout[c("train", "predict")]
  if (identical(learner$id, "regr.lightgbm")) {
    sf_require("lightgbm")
    sf_assert(utils::packageVersion("lightgbm") >= "4.0.0", "LightGBM >= 4.0.0 is required for portable callr/RDS model serialization.")
    if ("serializable" %in% learner$param_set$ids()) learner$param_set$values$serializable <- TRUE
  }
  learner
}
sf_recover_learners <- function(specs, cfg) {
  sf_assert(is.list(specs) && length(specs) > 0L && !is.null(names(specs)) && !anyNA(names(specs)) && all(nzchar(names(specs))) &&
    !anyDuplicated(names(specs)), "Original saved learner specifications must be a uniquely named list.")
  lapply(specs, function(sp) {
    if (inherits(sp, "Learner")) {
      sf_require(sp$packages)
      return(sf_configure_learner(sp$clone(deep = TRUE), cfg, sp$timeout))
    }
    sf_assert(is.list(sp) && is.character(sp$id) && length(sp$id) == 1L, "Invalid saved learner descriptor: expected id, params and packages.")
    sf_require(sp$packages %||% character())
    learner <- do.call(mlr3::lrn, c(list(sp$id), sp$params %||% list()))
    sf_configure_learner(learner, cfg, sp$timeout)
  })
}
sf_region <- function(x) {
  x <- as.character(x)
  map <- list(Scotland = c("11004", "11005"), Wales = c("11003", "11022", "11023"), `NE/Yorkshire` = c("11009",
    "11017", "11027", "11010", "11014"), `North West` = c("11008", "11001", "11016", "10003", "11024", "11025"),
    Midlands = c("11021", "11013", "11006"), South = c("11011", "11028", "11002", "11007", "11026"), London =
    c("11012", "11020", "11018"))
  out <- x
  for (nm in names(map)) out[x %in% map[[nm]]] <- nm
  sf_assert(all(is.na(out) | out %in% names(map)), "Unknown assessment_centre code/region; edit exact mapping, not a guessed recode.")
  factor(out)
}
sf_prepare_outcomes <- function(ctx, oof, cfg, covariates, blood, outcomes, disease_catalog) {
  sf_require("data.table")
  oof <- as.data.frame(oof)
  sf_assert(all(c("eid", "age", "I", "G") %in% names(oof)), "Clinical input lacks paired full-sample OOF accelerations.")
  sf_ids(oof$eid, "Paired OOF")
  cov <- as.data.frame(covariates)
  ceid <- if ("eid" %in% names(cov)) as.character(cov$eid) else rownames(cov)
  sf_ids(ceid, "Covariates")
  jj <- match(as.character(oof$eid), ceid);
  sf_assert(!anyNA(jj), "Missing OOF IDs in covariates.")
  cov <- cov[jj, , drop = FALSE];
  cov$age <- oof$age
  sf_assert(all(cfg$covariates %in% names(cov)), paste("Missing covariates:", paste(setdiff(cfg$covariates,
    names(cov)), collapse = ", ")))
  cov$assessment_centre <- sf_region(cov$assessment_centre)
  for (nm in c("age", "sex", "town_index", "bmi")) {
    raw <- as.character(cov[[nm]]);
    z <- suppressWarnings(as.numeric(raw))
    sf_assert(!any(!is.na(raw) & nzchar(raw) & is.na(z)), paste("Nonnumeric", nm, "requires explicit recoding."))
    z[!is.finite(z)] <- NA_real_;
    cov[[nm]] <- z
  }
  sf_assert(all(is.na(cov$sex) | cov$sex %in% c(0, 1)), "UKB sex must be0female/1male.")
  for (nm in setdiff(cfg$covariates, c("age", "sex", "town_index", "bmi"))) cov[[nm]] <- factor(cov[[nm]])
  blood <- as.data.frame(blood);
  sf_ids(blood$eid, "Blood dates")
  sf_assert("blood_date" %in% names(blood), "Blood file requires blood_date.")
  ob <- outcomes
  sf_assert(all(c("hosp_events", "death_data") %in% names(ob)), "Outcome RDS needs hosp_events and death_data.")
  death <- as.data.frame(ob$death_data);
  sf_ids(death$eid, "Death dates")
  dcols <- intersect(c("40000-0.0", "40000-1.0"), names(death))
  ds <- lapply(dcols, function(nm) sf_date_days(death[[nm]], nm))
  if (length(ds)) {
    dd <- do.call(pmin, c(ds, na.rm = TRUE));
    dd[!is.finite(dd)] <- NA_real_
    if ("death_date" %in% names(death)) {
      alt <- sf_date_days(death$death_date, "death_date");
      dd[is.na(dd)] <- alt[is.na(dd)]
    }
  } else {
    sf_assert("death_date" %in% names(death), "No recognized death date column.")
    dd <- sf_date_days(death$death_date, "death_date")
  }
  ib <- match(as.character(oof$eid), as.character(blood$eid));
  id <- match(as.character(oof$eid), as.character(death$eid))
  sf_assert(!anyNA(ib) && !anyNA(id), "OOF participant absent from blood/death tables.")
  base <- data.frame(eid = as.character(oof$eid), cov[cfg$covariates], blood_day = sf_date_days(blood$blood_date,
    "blood_date")[ib], death_day = dd[id], I = oof$I, G = oof$G, check.names = FALSE)
  ok <- stats::complete.cases(base[c(cfg$covariates, "blood_day", "I", "G")]) & (is.na(base$death_day) |
    base$death_day > base$blood_day)
  sf_assert(identical(cfg$covariate_missing, "complete_case"), "Supported covariate policy is original complete-case code; supply already imputed covariates if used in final analyses.")
  base <- base[ok, , drop = FALSE];
  rownames(base) <- NULL
  ev <- data.table::as.data.table(ob$hosp_events)
  sf_assert(all(c("eid", "icd3", "event_date") %in% names(ev)), "Hospital events require eid,icd3,event_date.")
  ev <- ev[as.character(eid) %in% base$eid, c("eid", "icd3", "event_date"), with = FALSE]
  ev[, eid := as.character(eid)];
  ev[, icd3 := toupper(trimws(as.character(icd3)))];
  ev[, dx_day := sf_date_days(event_date, "event_date")]
  ev <- ev[!is.na(dx_day) & grepl("^[A-N][0-9]{2}$", icd3)]
  first <- as.data.frame(ev[, list(dx_day = min(dx_day)), by = c("eid", "icd3")])
  cat <- as.data.frame(disease_catalog)
  names(cat)[names(cat) == "official_disease"] <- "disease"
  names(cat)[names(cat) == "official_chapter_label"] <- "category"
  req <- c("icd3", "disease", "category", "scope", "eligible_prevalent", "eligible_incident")
  sf_assert(all(req %in% names(cat)), "Disease catalog lacks required canonical fields.")
  for (nm in c("eligible_prevalent", "eligible_incident")) {
    z <- tolower(as.character(cat[[nm]]));
    sf_assert(all(z %in% c("true", "false", "1", "0", "yes", "no")), paste("Invalid catalog flag", nm));
    cat[[nm]] <- z %in% c("true", "1", "yes")
  }
  sf_assert(!anyDuplicated(cat$icd3), "Disease catalog contains duplicated ICD3 codes.")
  sf_assert(sum(cat$eligible_prevalent) == cfg$expected_prev && sum(cat$eligible_incident) == cfg$expected_inc,
    "Frozen disease catalog does not match117prevalent/274incident. Verify the intended catalog.")
  list(base = base, first = first, catalog = cat, fingerprint = sf_hash(list(base = base, first = first, catalog =
    cat, censor = c(cfg$disease_censor, cfg$death_censor))))
}
sf_prepare_comparator_context <- function(original, protein_data, cfg) {
  sf_require(c("mlr3", "mlr3learners", "mlr3extralearners", "data.table", "digest", "callr"))
  sf_assert(is.list(original), "Supply the original fitted pathway-clock result.")
  O <- as.data.frame(original$full_cv$oof, check.names = FALSE)
  sf_assert(all(c("eid", "label", "immune_age", "fold") %in% names(O)), "Original result lacks full_cv$oof fields.")
  ids <- sf_ids(O$eid, "Original OOF");
  age <- as.numeric(O$label)
  sf_assert(all(is.finite(age)), "Original age is not finite.")
  raw <- protein_data
  sf_assert(is.data.frame(raw) || is.matrix(raw), "Protein input must be a data.frame/matrix with participant row names.")
  rid <- sf_ids(rownames(raw), "Protein input row names")
  if ("eid" %in% colnames(raw)) raw <- raw[, setdiff(colnames(raw), "eid"), drop = FALSE]
  if ("label" %in% colnames(raw)) raw <- raw[, setdiff(colnames(raw), "label"), drop = FALSE]
  sf_assert(!anyNA(colnames(raw)) && !anyDuplicated(colnames(raw)), "Protein feature names are missing/duplicated.")
  ii <- match(ids, rid);
  sf_assert(!anyNA(ii), "Original OOF participants are absent from protein input.")
  sf_assert(all(vapply(as.data.frame(raw), is.numeric, logical(1))), "Protein input contains nonnumeric columns. Supply protein-only NPX data.")
  X <- as.matrix(raw[ii, , drop = FALSE]);
  storage.mode(X) <- "double";
  rownames(X) <- ids
  rm(raw);
  invisible(gc())
  sf_assert(ncol(X) == cfg$expected_proteins, paste("Expected", cfg$expected_proteins, "proteins, found",
    ncol(X), "; verify input before editing expectation."))
  sf_assert(!any(is.infinite(X)), "Protein input contains Inf.")
  if (identical(cfg$imputation, "none")) sf_assert(all(is.finite(X)), "Original input should be pre-imputed; missing values found. Set training_median only for an explicitly documented new preprocessing analysis.")
  tr <- match(as.character(original$split$train_eid), ids)
  te <- match(as.character(original$split$test_eid), ids)
  sf_assert(!anyNA(c(tr, te)) && !anyDuplicated(c(tr, te)) && setequal(c(tr, te), seq_along(ids)), "Saved train/test IDs do not partition original OOF population.")
  dp <- as.data.frame(original$split$train_prediction)
  dev_fold <- dp$fold[match(ids[tr], as.character(dp$eid))]
  full_fold <- O$fold
  sf_assert(!anyNA(dev_fold) && !anyNA(full_fold) && length(unique(dev_fold)) >= 2L && length(unique(full_fold)) >=
    2L, "Saved CV folds are missing.")
  learners <- sf_recover_learners(original$learners$stage1, cfg)
  fp <- sf_hash(list(ids = ids, age = age, X = X, train = tr, test = te, dev_fold = dev_fold, full_fold =
    full_fold, learners = lapply(learners, function(l) list(id = l$id, params = l$param_set$values)), preprocess =
    cfg$preprocessing_provenance))
  list(X = X, age = age, ids = ids, original = original, train_idx = tr, test_idx = te, dev_fold = dev_fold,
    full_fold = full_fold, learners = learners, fingerprint = fp)
}

# General-proteome learner selection and fixed-fold prediction
sf_model_engine_version <- "2026-09-24.4"
sf_learner_spec <- function(learner) {
  sf_assert(inherits(learner, "Learner"), "Learner must be a recovered mlr3 Learner")
  list(id = learner$id, params = learner$param_set$values, timeout = learner$timeout, execution = "callr-strict-no-fallback",
    packages = setNames(vapply(learner$packages, function(p) as.character(utils::packageVersion(p)), character(1)),
    learner$packages))
}
sf_model_signature <- function(ctx, cfg, features, learner = NULL, extra = NULL) {
  list(engine = sf_model_engine_version, input = ctx$fingerprint, features = features, learner = if (!is.null(learner)) sf_learner_spec(learner),
    imputation = cfg$imputation, threads = cfg$learner_threads, extra = extra)
}
sf_check_learner <- function(learner, context = "model") {
  condition_text <- function(x) {
    if (inherits(x, "condition")) return(conditionMessage(x))
    if (is.list(x)) return(paste(vapply(x, condition_text, character(1)), collapse = " | "))
    paste(as.character(x), collapse = " | ")
  }
  errors <- tryCatch(learner$errors, error = function(e) NULL)
  bad <- if (is.numeric(errors)) any(errors > 0, na.rm = TRUE) else length(errors) > 0L
  details <- if (!is.numeric(errors) && length(errors)) condition_text(errors) else character()
  lg <- tryCatch(learner$log, error = function(e) NULL)
  if (!is.null(lg) && NROW(lg)) {
    if ("class" %in% names(lg)) {
      ix <- which(tolower(as.character(lg$class)) == "error")
      bad <- bad || length(ix) > 0L
      if (length(ix) && "condition" %in% names(lg)) details <- c(details, vapply(lg$condition[ix], condition_text,
        character(1)))
    }
    txt <- paste(vapply(lg, condition_text, character(1)), collapse = " ")

    # Pre-training a fallback is normal. Only actual fallback use is rejected.
    bad <- bad || grepl("using fallback learner", txt, ignore.case = TRUE)
  }
  bad <- bad || is.null(learner$model)
  details <- unique(details[!is.na(details) & nzchar(details)])
  sf_assert(!bad, paste0(context, " failed, timed out, or invoked a fallback", if (length(details)) paste0(": ",
    paste(details, collapse = " | ")) else "; see audit log"))
  invisible(learner)
}
sf_prepare_matrix <- function(X, policy = "none", imputer = NULL) {
  X <- as.matrix(X)
  sf_assert(is.numeric(X) && ncol(X) > 0L && !is.null(colnames(X)), "Protein inputs must be a named numeric matrix")
  sf_assert(!anyDuplicated(colnames(X)), "Duplicate protein columns")
  sf_assert(!any(is.infinite(X)), "Infinite protein measurements are not allowed")
  if (is.null(imputer)) {
    if (identical(policy, "training_median")) {
      imputer <- apply(X, 2L, stats::median, na.rm = TRUE)
      sf_assert(all(is.finite(imputer)), "A training protein is entirely missing")
    } else {
      sf_assert(identical(policy, "none"), "imputation must be 'none' or 'training_median'")
      imputer <- setNames(rep(NA_real_, ncol(X)), colnames(X))
    }
  }
  if (anyNA(X)) {
    sf_assert(identical(policy, "training_median"), "Missing protein values: provide original finite matrix or explicitly use training_median")
    for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- imputer[colnames(X)[j]]
  }
  sf_assert(all(is.finite(X)), "Invalid values after training-derived imputation")
  list(X = X, imputer = imputer)
}
sf_fit_single <- function(X, y, learner, seed, cfg) {
  sf_require(c("mlr3", "data.table"))
  sf_assert(nrow(X) == length(y) && all(is.finite(y)), "Misaligned or invalid ages")
  policy <- if (is.null(cfg$imputation)) "none" else cfg$imputation
  pre <- sf_prepare_matrix(X, policy)
  features <- colnames(pre$X)
  internal <- sprintf("F%05d", seq_along(features))

  # The reconstruction path can supply the original global feature-name map.
  if (!is.null(cfg$internal_feature_map)) {
    internal <- unname(cfg$internal_feature_map[features])
    sf_assert(!anyNA(internal) && !anyDuplicated(internal), "Incomplete internal feature mapping")
  }
  dat <- data.table::as.data.table(pre$X);
  data.table::setnames(dat, internal)
  dat[, label := as.numeric(y)]
  task <- mlr3::TaskRegr$new("supp_age", dat, target = "label")
  l <- learner$clone(deep = TRUE);
  l$reset()

  # Recovered encapsulation/timeouts are retained. Any fallback is detected above.
  if (!is.null(cfg$learner_threads)) {
    thread_param <- intersect(c("num_threads", "num.threads", "nthread", "n_jobs", "threads"), l$param_set$ids())
    for (nm in thread_param) l$param_set$values[[nm]] <- as.integer(cfg$learner_threads)
  }
  set.seed(seed)
  cap <- sf_capture({
    l$train(task);
    sf_check_learner(l, l$id);
    invisible(NULL)
  })
  sf_log(cfg, paste0("fit_single:", l$id, ":seed=", seed), cap)
  if (!is.null(cap$error) && nzchar(as.character(cap$error))) stop(cap$error, call. = FALSE)

  # LightGBM native pointers do not survive ordinary RDS reloads reliably.
  # Keep the documented portable model string and restore before prediction.
  # https://lightgbm.readthedocs.io/en/stable/R/reference/lgb.load.html
  booster_model_str <- NULL
  booster_best_iter <- NULL
  if (inherits(l$model, "lgb.Booster")) {
    sf_require("lightgbm")
    booster_model_str <- l$model$save_model_to_string(NULL)
    booster_best_iter <- l$model$best_iter
  }
  list(learner = l, learner_name = learner$id, feature_names = features, internal_features = internal, imputer =
    pre$imputer, imputation = policy, seed = seed, booster_model_str = booster_model_str, booster_best_iter =
    booster_best_iter, specification = sf_learner_spec(learner), engine = sf_model_engine_version)
}
sf_predict_single <- function(model, X, cfg) {
  sf_assert(all(model$feature_names %in% colnames(X)), "Prediction matrix lacks trained protein features")
  pre <- sf_prepare_matrix(X[, model$feature_names, drop = FALSE], model$imputation, model$imputer)
  dat <- data.table::as.data.table(pre$X);
  data.table::setnames(dat, model$internal_features)
  l <- model$learner
  if (!is.null(model$booster_model_str)) {

    # Populate portable raw bytes before callr serializes the restored booster.
    # A new lgb.load() object alone still contains a process-local native pointer.
    restored <- lightgbm::lgb.load(model_str = model$booster_model_str)
    if (!is.null(model$booster_best_iter)) restored$best_iter <- model$booster_best_iter
    lightgbm::lgb.make_serializable(restored)
    l$model <- restored
  }
  limit <- l$timeout[["predict"]]
  message("Predict ", l$id, ": ", nrow(dat), " samples x ", ncol(dat), " features; time limit = ", if (is.infinite(limit)) "unlimited" else paste(limit,
    "seconds"))
  started <- proc.time()[["elapsed"]]
  cap <- sf_capture({
    q <- l$predict_newdata(dat);
    sf_check_learner(l, l$id);
    q
  })
  sf_log(cfg, paste0("predict_single:", l$id), cap)
  if (!is.null(cap$error) && nzchar(as.character(cap$error))) stop(cap$error, call. = FALSE)
  q <- cap$value
  j <- match(seq_len(nrow(X)), q$row_ids)
  sf_assert(!anyNA(j) && !anyDuplicated(q$row_ids), "Prediction row IDs are misaligned")
  p <- as.numeric(q$response[j])
  sf_assert(length(p) == nrow(X) && all(is.finite(p)), "Invalid single-stage predictions")
  message("Predict ", l$id, " completed in ", format(round(proc.time()[["elapsed"]] - started, 1), trim =
    TRUE), " seconds.")
  p
}
sf_single_oof <- function(X, y, fold, learner, seed, cfg, identity = NULL) {
  sf_assert(length(fold) == length(y) && !anyNA(fold) && length(unique(fold)) >= 2L, "Invalid OOF fold vector")
  p <- rep(NA_real_, length(y))
  for (v in sort(unique(fold))) {
    tr <- which(fold != v);
    va <- which(fold == v)
    label <- paste0(learner$id, ", fold ", v, "/", length(unique(fold)), " (train=", length(tr), ", validation=",
      length(va), ", features=", ncol(X), ")")
    compute <- function() {
      message("OOF train: ", label)
      tryCatch({
        fit <- sf_fit_single(X[tr, , drop = FALSE], y[tr], learner, seed + v, cfg)
        sf_predict_single(fit, X[va, , drop = FALSE], cfg)
      }, error = function(e) stop(label, ": ", conditionMessage(e), call. = FALSE))
    }
    if (is.null(identity)) p[va] <- compute() else {
      sig <- list(identity = identity, fold = fold, heldout = v, seed = seed, learner = sf_learner_spec(learner),
        engine = sf_model_engine_version)
      values <- sf_cache(paste0("oof_", identity$key, "_", v), sig, compute, cfg)
      sf_assert(is.numeric(values) && length(values) == length(va) && all(is.finite(values)), paste("Invalid saved OOF predictions for",
        label))
      p[va] <- values
      message("OOF ready: ", label, ".")
    }
  }
  sf_assert(all(is.finite(p)), "OOF prediction vector is incomplete")
  p
}
sf_select_single <- function(X, y, fold, learners, seed, cfg, identity = NULL) {
  sf_assert(length(learners) && !is.null(names(learners)), "Named learner candidates required")

  # The caller may supply an already fingerprinted data/row/feature identity.
  # Bind the computation identity to the actual matrix and outcome.
  selection_identity <- list(input = if (is.null(identity)) sf_hash(list(X = X, y = y)) else identity, features =
    colnames(X), n = nrow(X), fold = fold, seed = seed, imputation = cfg$imputation, threads = cfg$learner_threads,
    internal_feature_map = cfg$internal_feature_map, engine = sf_model_engine_version)
  selection_key <- sf_hash(selection_identity)
  scores <- predictions <- list();
  failures <- character()

  # All candidates receive identical folds. No test values enter this function.
  for (j in seq_along(learners)) {
    nm <- names(learners)[j]
    candidate_identity <- list(selection = selection_identity, candidate = nm, learner = sf_learner_spec(learners[[j]]),
      seed = seed + 10000L * j)
    candidate_identity$key <- paste0("candidate_", sf_hash(candidate_identity))
    message("Candidate [", j, "/", length(learners), "]: ", nm)
    cap <- sf_capture(sf_single_oof(X, y, fold, learners[[j]], seed + 10000L * j, cfg, identity = candidate_identity))
    sf_log(cfg, paste0("selection:", nm), cap)
    if (!is.null(cap$error) && nzchar(as.character(cap$error))) {
      failures[nm] <- as.character(cap$error)
      scores[[nm]] <- data.frame(learner = nm, MAE = NA_real_, RMSE = NA_real_, PCC = NA_real_, status =
        "failed", error = as.character(cap$error))
    } else {
      predictions[[nm]] <- cap$value
      met <- sf_metrics(y, cap$value)
      scores[[nm]] <- data.frame(learner = nm, MAE = met[["MAE"]], RMSE = met[["RMSE"]], PCC = met[["PCC"]],
        status = "ok", error = "")
    }

    # Retain completed candidate scores in the returned leaderboard.
  }
  leaderboard <- sf_bind(scores)
  if (length(failures) && !isTRUE(cfg$allow_failed_candidates)) {
    stop("Candidate comparison incomplete (no silent learner substitution). ", "Resolve the failed candidate before rerunning. ",
      paste(names(failures), failures, sep = ": ", collapse = "; "), call. = FALSE)
  }
  good <- leaderboard[leaderboard$status == "ok", , drop = FALSE]
  sf_assert(nrow(good) > 0L, "No candidate completed development cross-validation")
  good <- good[order(good$MAE, good$RMSE, - good$PCC, good$learner), , drop = FALSE]
  best <- good$learner[1L]
  leaderboard$rank <- match(leaderboard$learner, good$learner)
  leaderboard$selected <- leaderboard$learner == best

  # A failed selected refit is an error, not a switch informed by test results.
  fit <- sf_fit_single(X, y, learners[[best]], seed + 900000L, cfg)
  fit$learner_name <- best
  list(learner_name = best, leaderboard = leaderboard, fit = fit, pred_oof = predictions[[best]], all_candidate_oof =
    predictions)
}
sf_align_original <- function(tab, ids, value) {
  sf_assert(is.data.frame(tab) && all(c("eid", value) %in% names(tab)), paste("Original result lacks eid/",
    value))
  sf_assert(!anyDuplicated(as.character(tab$eid)), "Duplicate IDs in original prediction table")
  j <- match(as.character(ids), as.character(tab$eid))
  sf_assert(!anyNA(j), "Original prediction table does not cover required participant IDs")
  p <- as.numeric(tab[[value]][j])
  sf_assert(all(is.finite(p)), "Original predictions are not all finite")
  p
}
sf_full_oof <- function(ctx, s1, cfg) {
  nm <- s1$general_learner_name
  sf_assert(nm %in% names(ctx$learners), "Selected all-protein learner is not available")
  learner <- ctx$learners[[nm]]
  sig <- sf_model_signature(ctx, cfg, colnames(ctx$X), learner, list(fold = ctx$full_fold, seed = cfg$seed +
    20000L, role = "post_selection_full_OOF"))
  general <- sf_cache("general_full_oof", sig, function() sf_single_oof(ctx$X, ctx$age, ctx$full_fold, learner,
    cfg$seed + 20000L, cfg, identity = list(key = "general_full", signature = sig)), cfg)
  immune <- sf_align_original(ctx$original$full_cv$oof, ctx$ids, "immune_age")
  ia <- sf_acceleration(immune, ctx$age);
  ga <- sf_acceleration(general, ctx$age)
  out <- data.frame(eid = ctx$ids, age = ctx$age, immune_age = immune, general_age = general, I = ia$z, G =
    ga$z, fold = ctx$full_fold)
  cal <- sf_bind(list(data.frame(clock = "I", intercept = ia$intercept, slope = ia$slope, sd = ia$sd), data.frame(clock =
    "G", intercept = ga$intercept, slope = ga$slope, sd = ga$sd)))
  out
}

# The all-protein model occupied position j=2 in the original comparison.
# Preserve that fixed seed and the original saved development folds.
sf_select_general <- function(ctx, cfg) {
  sf_assert(ncol(ctx$X) == cfg$expected_proteins, "The comparator requires all 2,920 retained proteins.")
  tr <- ctx$train_idx;
  te <- ctx$test_idx;
  feats <- colnames(ctx$X)
  sig <- list(engine = sf_model_engine_version, input = ctx$fingerprint, features = feats, train = ctx$ids[tr],
    test = ctx$ids[te], fold = ctx$dev_fold, candidates = lapply(ctx$learners, sf_learner_spec), seed = cfg$seed +
    10100L * 2L, imputation = cfg$imputation, threads = cfg$learner_threads)
  fit <- sf_select_single(ctx$X[tr, feats, drop = FALSE], ctx$age[tr], ctx$dev_fold, ctx$learners, cfg$seed +
    10100L * 2L, cfg, identity = sig[setdiff(names(sig), "candidates")])
  list(general_learner_name = fit$learner_name, selection = fit)
}

# Separate immune and all-protein disease associations
sf_dis_bool <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "yes", "1")
}
sf_dis_scope <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("all", "all participants", "both", "both sexes")) return("All")
  if (x %in% c("female", "female only", "women", "females")) return("Female")
  if (x %in% c("male", "male only", "men", "males")) return("Male")
  stop("Unknown disease sex scope: ", x, ". Preserve the original catalog scope.")
}
sf_dis_riskset <- function(base, first, code, scope, type, cfg) {
  disease_end <- as.integer(as.Date(cfg$disease_censor))
  death_end <- as.integer(as.Date(cfg$death_censor))
  sf_assert(all(is.finite(c(disease_end, death_end))), "Invalid censor dates.")
  scope <- sf_dis_scope(scope)
  d <- base[base$blood_day <= disease_end, , drop = FALSE]
  if (scope == "Female") d <- d[d$sex == 0, , drop = FALSE]
  if (scope == "Male") d <- d[d$sex == 1, , drop = FALSE]
  dx <- first$dx_day[match(d$eid, first$eid)]
  if (type == "prevalent") {
    d$.outcome <- as.integer(!is.na(dx) & dx <= d$blood_day)
    d$.time <- NA_real_
    d$.competing <- 0L
  } else {
    end <- min(disease_end, death_end)
    keep <- (is.na(dx) | dx > d$blood_day) & d$blood_day < end
    d <- d[keep, , drop = FALSE]
    dx <- dx[keep]
    event <- !is.na(dx) & dx > d$blood_day & dx <= end & (is.na(d$death_day) | dx <= d$death_day)
    competing <- !event & !is.na(d$death_day) & d$death_day > d$blood_day & d$death_day <= end
    exit <- ifelse(event, dx, ifelse(competing, d$death_day, end))
    d$.outcome <- as.integer(event)
    d$.time <- (exit - d$blood_day) / 365.25
    d$.competing <- as.integer(competing)
    d <- d[is.finite(d$.time) & d$.time > 0, , drop = FALSE]
  }
  rownames(d) <- NULL
  d
}
sf_dis_fit_matrix <- function(X, y, time, type, exposure, indices = NULL, detailed = TRUE) {
  warnings <- character()
  error <- ""
  n <- length(y)
  if (is.null(indices)) indices <- seq_len(n)
  empty <- data.frame(clock = exposure, beta = NA_real_, se = NA_real_, effect = NA_real_, lower = NA_real_,
    upper = NA_real_, p = NA_real_, logp = NA_real_, neglog10p = NA_real_)
  ans <- tryCatch(withCallingHandlers({
    if (!length(indices) || length(unique(y[indices])) != 2L) stop("Both outcome states are required for model estimation.")
    if (type == "prevalent") {

      # The point fit uses every row once, with unit frequency weights.
      w <- tabulate(indices, nbins = n)
      fit <- stats::glm.fit(x = X, y = y, weights = w, family = stats::binomial(), control = stats::glm.control(maxit =
        50L))
      if (!isTRUE(fit$converged)) warnings <- c(warnings, "glm did not converge")
      class(fit) <- c("glm", "lm")
      sm <- summary(fit)$coefficients
      beta <- fit$coefficients[exposure]
      se <- setNames(rep(NA_real_, length(exposure)), exposure)
      good <- exposure %in% rownames(sm)
      se[good] <- sm[exposure[good], "Std. Error"]
    } else {

      # Retain the original point-model matrix and Efron handling of ties.
      xx <- X[indices, , drop = FALSE]
      yy <- survival::Surv(time[indices], y[indices])
      fit <- survival::coxph.fit(x = xx, y = yy, strata = rep.int(0L, length(indices)), offset = rep(0, length(indices)),
        init = NULL, control = survival::coxph.control(iter.max = 50L), weights = rep(1, length(indices)),
        method = "efron", rownames = as.character(seq_along(indices)), resid = FALSE)
      beta <- fit$coefficients[exposure]
      se <- setNames(sqrt(pmax(diag(fit$var), 0)), colnames(X))[exposure]
    }
    ok <- is.finite(beta) & is.finite(se) & se > 0
    out <- empty
    out$beta[ok] <- beta[ok]
    out$se[ok] <- se[ok]
    out$effect[ok] <- exp(beta[ok])
    out$lower[ok] <- exp(beta[ok] - stats::qnorm(.975) * se[ok])
    out$upper[ok] <- exp(beta[ok] + stats::qnorm(.975) * se[ok])
    out$logp[ok] <- log(2) + stats::pnorm( - abs(beta[ok] / se[ok]), log.p = TRUE)
    out$p[ok] <- exp(out$logp[ok])
    out$neglog10p[ok] <- - out$logp[ok] / log(10)
    if (any(!ok)) warnings <- c(warnings, "At least one clock coefficient is not estimable")
    out
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    empty
  })
  list(table = ans, warnings = unique(warnings), error = error, status = if (nzchar(error)) "error" else if (length(warnings)) "warning" else "ok")
}
sf_dis_design <- function(d, covariates, type) {

  # Removing a constant covariate (e.g. sex in a female-only disease) does not
  # omit an estimable adjustment term. Both clocks always use this same design.
  if (!nrow(d)) {
    mm <- matrix(numeric(), 0L, 2L, dimnames = list(NULL, c("I", "G")))
    return(list(I = mm[, "I", drop = FALSE], G = mm[, "G", drop = FALSE], covariates_used = character(), covariates_constant =
      covariates))
  }
  used <- covariates[vapply(d[covariates], function(v) {
    length(unique(v[!is.na(v)])) > 1L
  }, logical(1))]
  dd <- droplevels(d[c("I", "G", used)])
  for (nm in names(dd)) if (is.character(dd[[nm]])) dd[[nm]] <- factor(dd[[nm]])
  mm <- stats::model.matrix(stats::reformulate(c("I", "G", used)), dd)
  if (type == "incident") mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  sf_assert(all(is.finite(mm)), "Disease model matrix contains nonfinite values.")
  sf_assert(all(c("I", "G") %in% colnames(mm)), "Missing clock design columns.")
  list(I = mm[, colnames(mm) != "G", drop = FALSE], G = mm[, colnames(mm) != "I", drop = FALSE], covariates_used =
    used, covariates_constant = setdiff(covariates, used))
}
sf_dis_endpoint <- function(d, meta, type, cfg, signature) {
  ident <- data.frame(code = meta$icd3, icd3 = meta$icd3, disease = meta$disease, category = meta$category,
    scope = sf_dis_scope(meta$scope), type = type)
  n <- nrow(d)
  cases <- sum(d$.outcome)
  py <- if (type == "incident") sum(d$.time) else NA_real_
  counts <- cbind(ident, data.frame(n = n, cases = cases, controls = n - cases, events = if (type == "incident") cases else NA_integer_,
    person_years = py, competing_deaths = if (type == "incident") sum(d$.competing) else NA_integer_, prevalence =
    if (type == "prevalent" && n > 0L) cases / n else NA_real_, incidence_rate_1000py = if (type == "incident" &&
    py > 0) cases / py * 1000 else NA_real_, crude_event_proportion = if (type == "incident" && n > 0L) cases /
    n else NA_real_, median_observed_followup_years = if (type == "incident" && n > 0L) stats::median(d$.time) else NA_real_))
  design <- sf_dis_design(d, cfg$covariates, type)
  fit_i <- sf_dis_fit_matrix(design$I, d$.outcome, d$.time, type, "I")
  fit_g <- sf_dis_fit_matrix(design$G, d$.outcome, d$.time, type, "G")
  assoc <- cbind(counts[rep(1L, 2L), , drop = FALSE], rbind(fit_i$table, fit_g$table))
  audit <- sf_bind(lapply(list(I = fit_i, G = fit_g), function(a) {
    data.frame(status = a$status, warnings = paste(a$warnings, collapse = " | "), error = a$error)
  }))
  audit$model <- c("I", "G")
  audit$code <- meta$icd3
  audit$type <- type
  audit$constant_covariates <- paste(design$covariates_constant, collapse = ";")
  audit$below_original_event_threshold <- cases < cfg$min_cases || (type == "prevalent" && n - cases < cfg$min_cases)
  list(associations = assoc, counts = counts, diagnostics = audit)
}
sf_dis_adjust <- function(z, family_n) {
  z$m_tests <- as.integer(family_n)
  z$p_bonf <- pmin(1, z$p * family_n)
  z$logp_bonf <- pmin(0, z$logp + log(family_n))
  z
}
sf_run_disease <- function(clinical, type = c("prevalent", "incident"), cfg) {
  type <- match.arg(type)
  sf_require("survival")
  need <- c("eid", "I", "G", "blood_day", "death_day", cfg$covariates)
  sf_assert(all(need %in% names(clinical$base)), "Clinical base lacks disease-analysis columns.")
  sf_assert(all(c("icd3", "disease", "category", "scope", paste0("eligible_", type)) %in% names(clinical$catalog)),
    "Frozen disease catalog is incomplete.")
  catalog <- clinical$catalog[sf_dis_bool(clinical$catalog[[paste0("eligible_", type)]]), , drop = FALSE]
  catalog <- catalog[order(catalog$icd3), , drop = FALSE]
  rownames(catalog) <- NULL
  sf_assert(!anyNA(catalog[c("icd3", "disease", "category", "scope")]) && !anyDuplicated(catalog$icd3), "Frozen disease catalog must have unique codes and complete annotations.")
  expected <- if (type == "prevalent") cfg$expected_prev else cfg$expected_inc
  sf_assert(nrow(catalog) == expected, paste0("Frozen ", type, " catalog has ", nrow(catalog), " endpoints; expected ",
    expected, "."))
  base <- as.data.frame(clinical$base)
  sf_assert(!anyDuplicated(base$eid) && !anyNA(base$eid), "Clinical base EIDs must be unique and nonmissing.")
  sf_assert(all(stats::na.omit(base$sex) %in% c(0, 1)), "Disease sex coding must be 0 female / 1 male.")
  required <- unique(c("I", "G", "blood_day", cfg$covariates))
  complete <- stats::complete.cases(base[required])
  for (nm in required[vapply(base[required], is.numeric, logical(1))]) complete <- complete & is.finite(base[[nm]])
  death_invalid <- !is.na(base$death_day) & base$death_day <= base$blood_day
  keep <- complete & !death_invalid
  base <- base[keep, , drop = FALSE]
  sf_assert(nrow(base) > 0L, "No participants remain for paired disease models.")
  first <- as.data.frame(clinical$first)
  sf_assert(all(c("eid", "icd3", "dx_day") %in% names(first)), "Missing earliest-diagnosis columns.")
  first <- first[first$icd3 %in% catalog$icd3 & first$eid %in% base$eid & !is.na(first$dx_day), , drop =
    FALSE]
  sf_assert(!anyDuplicated(paste(first$eid, first$icd3, sep = "|")), "clinical$first must contain earliest date per EID/code.")
  split_first <- split(first[c("eid", "dx_day")], first$icd3)
  empty_first <- first[FALSE, c("eid", "dx_day"), drop = FALSE]
  signature <- sf_hash(list(version = "disease-point-only-v2", clinical = clinical$fingerprint, base = sf_hash(base[c("eid",
    "I", "G", "blood_day", "death_day", cfg$covariates)]), first = sf_hash(first), catalog = catalog, type =
    type, settings = cfg[c("min_cases", "covariates", "disease_censor", "death_censor", "alpha")]))
  stage <- if (type == "prevalent") "Prevalent comparator" else "Incident comparator"
  message(stage, ": ", nrow(catalog), " frozen endpoints; separate immune and all-protein models with Wald inference.",
    " Association bootstrap is disabled.")
  output <- sf_assoc_parallel(seq_len(nrow(catalog)), function(k) {
    meta <- catalog[k, , drop = FALSE]
    key <- paste0("disease_", type, "_", meta$icd3)
    endpoint_sig <- sf_hash(list(signature, code = meta$icd3))
    value <- sf_cache(key, endpoint_sig, function() {
      dx <- split_first[[as.character(meta$icd3)]]
      if (is.null(dx)) dx <- empty_first
      d <- sf_dis_riskset(base, dx, meta$icd3, meta$scope, type, cfg)
      sf_dis_endpoint(d, meta, type, cfg, endpoint_sig)
    }, cfg)
    value
  }, cfg, labels = paste(type, catalog$icd3), stage = stage)

  # Worker failures are explicit errors; do not silently remove endpoints.
  sf_assert(all(vapply(output, function(z) is.list(z) && all(c("associations", "counts", "diagnostics") %in%
    names(z)), logical(1))), "At least one disease endpoint failed outside fitting; inspect the worker error and rerun.")
  associations <- sf_dis_adjust(sf_bind(lapply(output, `[[`, "associations")), nrow(catalog))
  diagnostics <- sf_bind(lapply(output, `[[`, "diagnostics"))
  analysis_audit <- data.frame(type = type, input_n = nrow(clinical$base), excluded_incomplete = sum(!complete),
    excluded_invalid_death = sum(complete & death_invalid), common_complete_n = nrow(base), frozen_endpoints =
    nrow(catalog))
  list(associations = associations, counts = sf_bind(lapply(output, `[[`, "counts")), diagnostics = diagnostics,
    analysis_audit = analysis_audit, metadata = list(type = type, family_per_clock = nrow(catalog), association_inference =
    "Separate immune and all-protein model coefficients, standard errors, Wald 95% CI and P values; no association bootstrap",
    significance = "Finite warning estimates retained; nonestimable coefficients remain NA", fingerprint =
    signature))
}

# Separate mortality spline models
sf_mort_or <- function(x, default) if (is.null(x)) default else x
sf_mort_basis <- function(z, knots) {
  as.matrix(splines::ns(z, knots = knots[2:3], Boundary.knots = knots[c(1, 4)], intercept = FALSE))
}
sf_mort_capture_row <- function(cap, context) {
  data.frame(context = context, status = if (isTRUE(nzchar(cap$error))) "error" else if (length(cap$warnings)) "warning" else "ok",
    detail = paste(c(as.character(cap$error), as.character(cap$warnings)), collapse = " | "), stringsAsFactors =
    FALSE)
}
sf_mort_curve <- function(fit, knots, grid, clock) {
  out <- data.frame(clock = clock, z = grid, hr = NA_real_, lower = NA_real_, upper = NA_real_)
  if (is.null(fit)) return(out)
  nm <- paste0(".s", seq_len(3L))
  b <- stats::coef(fit);
  V <- stats::vcov(fit)
  if (!all(nm %in% names(b)) || any(!is.finite(b[nm]))) return(out)
  delta <- sweep(sf_mort_basis(grid, knots), 2, as.numeric(sf_mort_basis(0, knots)), FUN = "-")
  vv <- V[nm, nm, drop = FALSE]
  if (any(!is.finite(vv))) return(out)
  eta <- as.vector(delta %*% b[nm])
  se <- sqrt(pmax(0, rowSums((delta %*% vv) * delta)))
  out$hr <- exp(eta)
  out$lower <- exp(eta - stats::qnorm(.975) * se)
  out$upper <- exp(eta + stats::qnorm(.975) * se)
  out
}
sf_mort_fit <- function(d, cov_cols, clock, knots, iter_max = 50L) {
  dd <- d
  ss <- sf_mort_basis(dd[[clock]], knots)
  colnames(ss) <- paste0(".s", seq_len(ncol(ss)))
  dd <- cbind(dd, as.data.frame(ss))
  exposure_cols <- colnames(ss)
  rhs <- c(exposure_cols, cov_cols)
  form <- stats::as.formula(paste("survival::Surv(.time, .event) ~", if (length(rhs)) paste(rhs, collapse =
    " + ") else "1"))
  survival::coxph(form, data = dd, ties = "efron", singular.ok = TRUE, na.action = stats::na.fail, x = TRUE,
    y = TRUE, model = FALSE, control = survival::coxph.control(iter.max = iter_max))
}
sf_run_mortality_comparison <- function(clinical, cfg) {
  sf_require(c("survival", "splines"))
  b <- clinical$base
  req <- unique(c("eid", "blood_day", "death_day", "I", "G", cfg$covariates))
  sf_assert(all(req %in% names(b)), paste("Mortality input is missing:", paste(setdiff(req, names(b)), collapse =
    ", ")))
  sf_assert(!anyDuplicated(b$eid), "Mortality base must have exactly one row per participant.")
  sf_assert(all(vapply(b[c("blood_day", "death_day")], function(x) is.numeric(x) || inherits(x, "Date"), logical(1))),
    "blood_day and death_day must be integer-day dates (not Excel serial strings).")
  censor <- as.numeric(as.Date(cfg$death_censor))
  sf_assert(is.finite(censor), "cfg$death_censor is not a valid YYYY-MM-DD date.")
  knot_probs <- sf_mort_or(cfg$mortality_knots_probs, c(.05, .35, .65, .95))
  sf_assert(length(knot_probs) == 4L && all(diff(knot_probs) > 0) && min(knot_probs) > 0 && max(knot_probs) <
    1, "Specify four increasing mortality knot probabilities strictly between 0 and 1.")
  curve_probs <- sf_mort_or(cfg$mortality_curve_probs, c(.01, .99))
  sf_assert(length(curve_probs) == 2L && all(diff(curve_probs) > 0) && all(curve_probs >= 0 & curve_probs <=
    1), "mortality_curve_probs requires two increasing probabilities.")
  iter_max <- as.integer(sf_mort_or(cfg$cox_iter_max, 50L))
  blood <- as.numeric(b$blood_day);
  death <- as.numeric(b$death_day)
  eligible <- is.finite(blood) & blood < censor & (!is.finite(death) | death > blood) & is.finite(b$I) &
  is.finite(b$G)
  sf_assert(any(eligible), "No eligible participants remain for mortality analysis.")
  diag_rows <- list(data.frame(context = "eligibility", status = "info", detail = sprintf("Input %d; eligible %d; excluded %d (dates or nonfinite OOF acceleration).",
    nrow(b), sum(eligible), sum(!eligible))))
  b <- b[eligible, , drop = FALSE]
  sf_assert(all(stats::complete.cases(b[, cfg$covariates, drop = FALSE])), "Mortality covariates contain missing values: apply the prespecified common preprocessing first.")

  # Build one shared numeric design matrix for both clocks and all models.
  for (v in cfg$covariates) if (is.character(b[[v]])) b[[v]] <- factor(b[[v]])
  varying <- vapply(b[, cfg$covariates, drop = FALSE], function(x) length(unique(x)) > 1L, logical(1))
  covars <- cfg$covariates[varying]
  if (any(!varying)) diag_rows[[length(diag_rows) + 1L]] <- data.frame(context = "covariate_design", status =
    "info", detail = paste("Constant covariates omitted:", paste(cfg$covariates[!varying], collapse = ", ")))
  mm <- stats::model.matrix(if (length(covars)) stats::reformulate(covars) else ~ 1, data = b)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  sf_assert(all(is.finite(mm)), "Mortality covariate design contains nonfinite values.")
  if (ncol(mm)) mm <- mm[, apply(mm, 2, function(x) length(unique(x)) > 1L), drop = FALSE]
  orig_cov_names <- colnames(mm)
  cov_cols <- sprintf(".c%03d", seq_len(ncol(mm)));
  colnames(mm) <- cov_cols
  death <- as.numeric(b$death_day);
  blood <- as.numeric(b$blood_day)
  end <- pmin(ifelse(is.finite(death), death, censor), censor)
  d <- data.frame(eid = as.character(b$eid), .time = (end - blood) / 365.25, .event = as.integer(is.finite(death) &
    death <= censor), I = as.numeric(b$I), G = as.numeric(b$G), stringsAsFactors = FALSE)
  if (ncol(mm)) d <- cbind(d, as.data.frame(mm))
  sf_assert(all(d$.time > 0), "Mortality follow-up must be positive.")
  sf_assert(sum(d$.event) > 0L, "No observed deaths in mortality analysis.")
  knots <- lapply(c("I", "G"), function(cl) as.numeric(stats::quantile(d[[cl]], knot_probs, names = FALSE)))
  names(knots) <- c("I", "G")
  sf_assert(all(vapply(knots, function(k) all(is.finite(k)) && all(diff(k) > 0), logical(1))), "Spline knots are not distinct: age acceleration requires a nondegenerate distribution.")
  clocks <- c("I", "G")
  fits <- sf_assoc_parallel(clocks, function(cl) sf_capture(sf_mort_fit(d, cov_cols, cl, knots[[cl]], iter_max =
    iter_max)), cfg, labels = clocks, stage = "Mortality spline comparison")
  names(fits) <- clocks
  spline_rows <- list()
  for (cl in clocks) {
    cap <- fits[[cl]]
    sf_log(cfg, paste("Mortality spline", cl), cap)
    diag_rows[[length(diag_rows) + 1L]] <- sf_mort_capture_row(cap, paste0(cl, "_spline"))
    grid_lim <- stats::quantile(d[[cl]], curve_probs, names = FALSE)
    grid <- sort(unique(c(seq(grid_lim[1], grid_lim[2], length.out = 151L), 0)))
    spline_rows[[cl]] <- sf_mort_curve(cap$value, knots[[cl]], grid, cl)
  }
  knot_table <- sf_bind(lapply(clocks, function(cl) data.frame(clock = cl, probability = knot_probs, z =
    knots[[cl]])))
  list(spline = sf_bind(spline_rows), knots = knot_table, models = lapply(fits, `[[`, "value"), diagnostics =
    sf_bind(diag_rows), n = nrow(d), deaths = sum(d$.event), covariate_design = data.frame(model_column =
    cov_cols, original_column = orig_cov_names))
}

# Parallel association execution
sf_assoc_parallel <- function(X, FUN, cfg, labels = as.character(X), stage = "Association") {
  if (!length(X)) return(list())
  available <- if (requireNamespace("parallelly", quietly = TRUE)) as.integer(parallelly::availableCores()) else parallel::detectCores(logical =
    FALSE)
  if (!is.finite(available) || available < 1L) available <- 1L
  requested <- cfg$association_cores %||% cfg$cores %||% 1L
  sf_assert(length(requested) == 1L && is.finite(requested) && requested >= 1, "association_cores must be a positive integer.")
  workers <- min(as.integer(requested), available, length(X))
  sf_assert(length(labels) == length(X), "Parallel progress labels must match jobs.")
  message(stage, ": ", length(X), " jobs, ", workers, " worker(s); no association resampling.")

  # Keep this environment small: do not serialize a PSOCK cluster/socket in FUN.
  worker_env <- new.env(parent = .GlobalEnv)
  worker_env$values <- X;
  worker_env$task <- FUN;
  worker_env$labels <- labels
  worker_env$prefix <- stage;
  worker_env$n <- length(X)
  worker <- function(k) {
    start <- proc.time()[[3L]]
    message(prefix, " [", k, "/", n, "] START ", labels[[k]])
    result <- tryCatch(task(values[[k]]), error = function(e) stop(prefix, " job ", labels[[k]], ": ", conditionMessage(e),
      call. = FALSE))
    message(prefix, " [", k, "/", n, "] DONE ", labels[[k]], " (", format(round(proc.time()[[3L]] - start,
      1), nsmall = 1), " s)")
    result
  }
  environment(worker) <- worker_env
  if (workers == 1L) return(lapply(seq_along(X), worker))
  use_psock <- .Platform$OS.type == "windows" || identical(cfg$association_backend, "psock")
  if (!use_psock) {
    ans <- parallel::mclapply(seq_along(X), worker, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed =
      FALSE)
  } else {
    cl <- parallel::makePSOCKcluster(workers, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export the already defined analytical functions to workers. The integrated
    # script has no external R source-file dependency.
    env <- environment(sf_assoc_parallel)
    defs <- mget(c(ls(env, pattern = "^sf_"), "%||%"), envir = env, inherits = TRUE)
    parallel::clusterCall(cl, function(defs, libs) {
      .libPaths(libs)
      list2env(defs, envir = .GlobalEnv)
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      TRUE
    }, defs, .libPaths())
    ans <- parallel::parLapplyLB(cl, seq_along(X), worker)
  }
  failed <- which(vapply(ans, inherits, logical(1), "try-error"))
  sf_assert(!length(failed), paste0(stage, " parallel job failed: ", paste(labels[failed], collapse = ", ")))
  ans
}

# Prespecified disease families
sf_frozen_catalog <- function() {
  chapter_labels <- c("Certain infectious and parasitic diseases", "Neoplasms", "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
    "Endocrine, nutritional and metabolic diseases", "Mental and behavioural disorders", "Diseases of the nervous system",
    "Diseases of the eye and adnexa", "Diseases of the ear and mastoid process", "Diseases of the circulatory system",
    "Diseases of the respiratory system", "Diseases of the digestive system", "Diseases of the skin and subcutaneous tissue",
    "Diseases of the musculoskeletal system and connective tissue", "Diseases of the genitourinary system")
  fixed <- utils::read.delim(text = "A04|1|A|1|1|Other bacterial intestinal infections
A08|1|A|0|1|Viral and other specified intestinal infections
A09|1|A|0|1|Diarrhoea and gastro-enteritis of presumed infectious origin
A41|1|A|0|1|Other septicaemia
A49|1|A|0|1|Bacterial infection of unspecified site
B34|1|A|0|1|Viral infection of unspecified site
B37|1|A|0|1|Candidiasis
B95|1|A|1|1|Streptococcus and staphylococcus as the cause of diseases classified to other chapters
B96|1|A|1|1|Other bacterial agents as the cause of diseases classified to other chapters
B97|1|A|0|1|Viral agents as the cause of diseases classified to other chapters
B98|1|A|0|1|Other specified infectious agents as the cause of diseases classified to other chapters
C18|2|A|0|1|Malignant neoplasm of colon
C20|2|A|0|1|Malignant neoplasm of rectum
C34|2|A|0|1|Malignant neoplasm of bronchus and lung
C43|2|A|0|1|Malignant melanoma of skin
C44|2|A|1|1|Other malignant neoplasms of skin
C50|2|A|1|1|Malignant neoplasm of breast
C61|2|M|1|1|Malignant neoplasm of prostate
C64|2|A|0|1|Malignant neoplasm of kidney, except renal pelvis
C67|2|A|0|1|Malignant neoplasm of bladder
C77|2|A|1|1|Secondary and unspecified malignant neoplasm of lymph nodes
C78|2|A|0|1|Secondary malignant neoplasm of respiratory and digestive organs
C79|2|A|0|1|Secondary malignant neoplasm of other sites
C80|2|A|0|1|Malignant neoplasm without specification of site
C97|2|A|0|1|Malignant neoplasms of independent (primary) multiple sites
D05|2|A|0|1|Carcinoma in situ of breast
D12|2|A|1|1|Benign neoplasm of colon, rectum, anus and anal canal
D17|2|A|1|1|Benign lipomatous neoplasm
D18|2|A|0|1|Haemangioma and lymphangioma, any site
D22|2|A|1|1|Melanocytic naevi
D23|2|A|1|1|Other benign neoplasms of skin
D25|2|F|1|1|Leiomyoma of uterus
D35|2|A|0|1|Benign neoplasm of other and unspecified endocrine glands
D37|2|A|0|1|Neoplasm of uncertain or unknown behaviour of oral cavity and digestive organs
D47|2|A|0|1|Other neoplasms of uncertain or unknown behaviour of lymphoid, haematopoietic and related tissue
D50|3|A|1|1|Iron deficiency anaemia
D51|3|A|0|1|Vitamin B12 deficiency anaemia
D63|3|A|0|1|Anaemia in chronic diseases classified elsewhere
D64|3|A|1|1|Other anaemias
D68|3|A|0|1|Other coagulation defects
D69|3|A|0|1|Purpura and other haemorrhagic conditions
D70|3|A|0|1|Agranulocytosis
E03|4|A|1|1|Other hypothyroidism
E04|4|A|0|1|Other non-toxic goitre
E05|4|A|0|1|Thyrotoxicosis [hyperthyroidism]
E10|4|A|1|1|Insulin-dependent diabetes mellitus
E11|4|A|1|1|Non-insulin-dependent diabetes mellitus
E14|4|A|0|1|Unspecified diabetes mellitus
E16|4|A|0|1|Other disorders of pancreatic internal secretion
E21|4|A|0|1|Hyperparathyroidism and other disorders of parathyroid gland
E53|4|A|0|1|Deficiency of other B group vitamins
E55|4|A|0|1|Vitamin D deficiency
E66|4|A|1|1|Obesity
E78|4|A|1|1|Disorders of lipoprotein metabolism and other lipidaemias
E83|4|A|0|1|Disorders of mineral metabolism
E86|4|A|0|1|Volume depletion
E87|4|A|0|1|Other disorders of fluid, electrolyte and acid-base balance
E89|4|A|0|1|Postprocedural endocrine and metabolic disorders, not elsewhere classified
F00|5|A|0|1|Dementia in Alzheimer's disease
F01|5|A|0|1|Vascular dementia
F02|5|A|0|1|Dementia in other diseases classified elsewhere
F03|5|A|0|1|Unspecified dementia
F05|5|A|0|1|Delirium, not induced by alcohol and other psychoactive substances
F10|5|A|1|1|Mental and behavioural disorders due to use of alcohol
F17|5|A|0|1|Mental and behavioural disorders due to use of tobacco
F32|5|A|1|1|Depressive episode
F40|5|A|0|1|Phobic anxiety disorders
F41|5|A|0|1|Other anxiety disorders
G12|6|A|0|1|Spinal muscular atrophy and related syndromes
G20|6|A|0|1|Parkinson's disease
G25|6|A|0|1|Other extrapyramidal and movement disorders
G30|6|A|0|1|Alzheimer's disease
G31|6|A|0|1|Other degenerative diseases of nervous system, not elsewhere classified
G35|6|A|1|1|Multiple sclerosis
G40|6|A|1|1|Epilepsy
G43|6|A|0|1|Migraine
G45|6|A|0|1|Transient cerebral ischaemic attacks and related syndromes
G47|6|A|1|1|Sleep disorders
G55|6|A|1|1|Nerve root and plexus compressions in diseases classified elsewhere
G56|6|A|1|1|Mononeuropathies of upper limb
G62|6|A|0|1|Other polyneuropathies
G63|6|A|0|1|Polyneuropathy in diseases classified elsewhere
G81|6|A|0|1|Hemiplegia
G93|6|A|0|1|Other disorders of brain
G99|6|A|0|1|Other disorders of nervous system in diseases classified elsewhere
H02|7|A|1|1|Other disorders of eyelid
H04|7|A|0|1|Disorders of lachrymal system
H18|7|A|0|1|Other disorders of cornea
H25|7|A|1|1|Senile cataract
H26|7|A|1|1|Other cataract
H33|7|A|1|1|Retinal detachments and breaks
H35|7|A|0|1|Other retinal disorders
H36|7|A|0|1|Retinal disorders in diseases classified elsewhere
H40|7|A|0|1|Glaucoma
H43|7|A|0|1|Disorders of vitreous body
H52|7|A|0|1|Disorders of refraction and accommodation
H53|7|A|1|1|Visual disturbances
H54|7|A|0|1|Blindness and low vision
H81|8|A|0|1|Disorders of vestibular function
H90|8|A|0|1|Conductive and sensorineural hearing loss
H91|8|A|0|1|Other hearing loss
I07|9|A|0|1|Rheumatic tricuspid valve diseases
I08|9|A|0|1|Multiple valve diseases
I10|9|A|1|1|Essential (primary) hypertension
I12|9|A|0|1|Hypertensive renal disease
I20|9|A|1|1|Angina pectoris
I21|9|A|1|1|Acute myocardial infarction
I24|9|A|0|1|Other acute ischaemic heart diseases
I25|9|A|1|1|Chronic ischaemic heart disease
I26|9|A|1|1|Pulmonary embolism
I27|9|A|0|1|Other pulmonary heart diseases
I31|9|A|0|1|Other diseases of pericardium
I34|9|A|0|1|Nonrheumatic mitral valve disorders
I35|9|A|0|1|Nonrheumatic aortic valve disorders
I42|9|A|0|1|Cardiomyopathy
I44|9|A|0|1|Atrioventricular and left bundle-branch block
I45|9|A|0|1|Other conduction disorders
I46|9|A|0|1|Cardiac arrest
I47|9|A|1|1|Paroxysmal tachycardia
I48|9|A|1|1|Atrial fibrillation and flutter
I49|9|A|0|1|Other cardiac arrhythmias
I50|9|A|1|1|Heart failure
I51|9|A|0|1|Complications and ill-defined descriptions of heart disease
I61|9|A|0|1|Intracerebral haemorrhage
I63|9|A|0|1|Cerebral infarction
I65|9|A|0|1|Occlusion and stenosis of precerebral arteries, not resulting in cerebral infarction
I67|9|A|0|1|Other cerebrovascular diseases
I69|9|A|0|1|Sequelae of cerebrovascular disease
I70|9|A|0|1|Atherosclerosis
I71|9|A|0|1|Aortic aneurysm and dissection
I73|9|A|1|1|Other peripheral vascular diseases
I74|9|A|0|1|Arterial embolism and thrombosis
I77|9|A|0|1|Other disorders of arteries and arterioles
I80|9|A|1|1|Phlebitis and thrombophlebitis
I83|9|A|1|1|Varicose veins of lower extremities
I84|9|A|1|1|Haemorrhoids
I89|9|A|0|1|Other non-infective disorders of lymphatic vessels and lymph nodes
I95|9|A|0|1|Hypotension
J10|10|A|0|1|Influenza due to identified influenza virus
J12|10|A|0|1|Viral pneumonia, not elsewhere classified
J18|10|A|1|1|Pneumonia, organism unspecified
J22|10|A|1|1|Unspecified acute lower respiratory infection
J30|10|A|0|1|Vasomotor and allergic rhinitis
J32|10|A|0|1|Chronic sinusitis
J33|10|A|1|1|Nasal polyp
J34|10|A|1|1|Other disorders of nose and nasal sinuses
J38|10|A|0|1|Diseases of vocal cords and larynx, not elsewhere classified
J43|10|A|0|1|Emphysema
J44|10|A|1|1|Other chronic obstructive pulmonary disease
J45|10|A|1|1|Asthma
J47|10|A|0|1|Bronchiectasis
J69|10|A|0|1|Pneumonitis due to solids and liquids
J84|10|A|0|1|Other interstitial pulmonary diseases
J90|10|A|1|1|Pleural effusion, not elsewhere classified
J92|10|A|0|1|Pleural plaque
J93|10|A|0|1|Pneumothorax
J96|10|A|0|1|Respiratory failure, not elsewhere classified
J98|10|A|0|1|Other respiratory disorders
K02|11|A|1|1|Dental caries
K04|11|A|1|1|Diseases of pulp and periapical tissues
K08|11|A|1|1|Other disorders of teeth and supporting structures
K12|11|A|0|1|Stomatitis and related lesions
K13|11|A|0|1|Other diseases of lip and oral mucosa
K14|11|A|0|1|Diseases of tongue
K20|11|A|1|1|Oesophagitis
K21|11|A|1|1|Gastro-oesophageal reflux disease
K22|11|A|1|1|Other diseases of oesophagus
K25|11|A|1|1|Gastric ulcer
K26|11|A|1|1|Duodenal ulcer
K29|11|A|1|1|Gastritis and duodenitis
K30|11|A|1|1|Dyspepsia
K31|11|A|1|1|Other diseases of stomach and duodenum
K35|11|A|1|1|Acute appendicitis
K40|11|A|1|1|Inguinal hernia
K42|11|A|1|1|Umbilical hernia
K43|11|A|1|1|Ventral hernia
K44|11|A|1|1|Diaphragmatic hernia
K51|11|A|1|1|Ulcerative colitis
K52|11|A|1|1|Other non-infective gastro-enteritis and colitis
K55|11|A|0|1|Vascular disorders of intestine
K56|11|A|1|1|Paralytic ileus and intestinal obstruction without hernia
K57|11|A|1|1|Diverticular disease of intestine
K58|11|A|1|1|Irritable bowel syndrome
K59|11|A|1|1|Other functional intestinal disorders
K60|11|A|1|1|Fissure and fistula of anal and rectal regions
K62|11|A|1|1|Other diseases of anus and rectum
K63|11|A|1|1|Other diseases of intestine
K64|11|A|0|1|Haemorrhoids and perianal venous thrombosis
K66|11|A|0|1|Other disorders of peritoneum
K74|11|A|0|1|Fibrosis and cirrhosis of liver
K76|11|A|0|1|Other diseases of liver
K80|11|A|1|1|Cholelithiasis
K81|11|A|1|1|Cholecystitis
K82|11|A|0|1|Other diseases of gallbladder
K83|11|A|0|1|Other diseases of biliary tract
K85|11|A|0|1|Acute pancreatitis
K86|11|A|0|1|Other diseases of pancreas
K90|11|A|0|1|Intestinal malabsorption
K91|11|A|0|1|Postprocedural disorders of digestive system, not elsewhere classified
K92|11|A|1|1|Other diseases of digestive system
L02|12|A|0|1|Cutaneous abscess, furuncle and carbuncle
L03|12|A|1|1|Cellulitis
L08|12|A|0|1|Other local infections of skin and subcutaneous tissue
L30|12|A|0|1|Other dermatitis
L40|12|A|0|1|Psoriasis
L57|12|A|0|1|Skin changes due to chronic exposure to nonionising radiation
L72|12|A|1|1|Follicular cysts of skin and subcutaneous tissue
L82|12|A|0|1|Seborrhoeic keratosis
L89|12|A|0|1|Decubitus ulcer
L90|12|A|1|1|Atrophic disorders of skin
L97|12|A|0|1|Ulcer of lower limb, not elsewhere classified
L98|12|A|1|1|Other disorders of skin and subcutaneous tissue, not elsewhere classified
M06|13|A|1|1|Other rheumatoid arthritis
M10|13|A|0|1|Gout
M13|13|A|1|1|Other arthritis
M15|13|A|0|1|Polyarthrosis
M16|13|A|1|1|Coxarthrosis [arthrosis of hip]
M17|13|A|1|1|Gonarthrosis [arthrosis of knee]
M18|13|A|0|1|Arthrosis of first carpometacarpal joint
M19|13|A|1|1|Other arthrosis
M20|13|A|1|1|Acquired deformities of fingers and toes
M21|13|A|0|1|Other acquired deformities of limbs
M23|13|A|1|1|Internal derangement of knee
M24|13|A|0|1|Other specific joint derangements
M25|13|A|1|1|Other joint disorders, not elsewhere classified
M32|13|A|1|1|Systemic lupus erythematosus
M35|13|A|0|1|Other systemic involvement of connective tissue
M41|13|A|0|1|Scoliosis
M43|13|A|0|1|Other deforming dorsopathies
M46|13|A|0|1|Other inflammatory spondylopathies
M47|13|A|1|1|Spondylosis
M48|13|A|0|1|Other spondylopathies
M50|13|A|0|1|Cervical disk disorders
M51|13|A|1|1|Other intervertebral disk disorders
M54|13|A|1|1|Dorsalgia
M62|13|A|0|1|Other disorders of muscle
M65|13|A|1|1|Synovitis and tenosynovitis
M67|13|A|1|1|Other disorders of synovium and tendon
M70|13|A|0|1|Soft tissue disorders related to use, overuse and pressure
M72|13|A|1|1|Fibroblastic disorders
M75|13|A|1|1|Shoulder lesions
M77|13|A|0|1|Other enthesopathies
M79|13|A|1|1|Other soft tissue disorders, not elsewhere classified
M80|13|A|0|1|Osteoporosis with pathological fracture
M81|13|A|1|1|Osteoporosis without pathological fracture
M84|13|A|0|1|Disorders of continuity of bone
M85|13|A|0|1|Other disorders of bone density and structure
M89|13|A|0|1|Other disorders of bone
N12|14|A|0|1|Tubulo-interstitial nephritis, not specified as acute or chronic
N13|14|A|0|1|Obstructive and reflux uropathy
N17|14|A|0|1|Acute renal failure
N18|14|A|0|1|Chronic renal failure
N19|14|A|0|1|Unspecified renal failure
N20|14|A|1|1|Calculus of kidney and ureter
N28|14|A|0|1|Other disorders of kidney and ureter, not elsewhere classified
N30|14|A|1|1|Cystitis
N32|14|A|1|1|Other disorders of bladder
N35|14|A|1|1|Urethral stricture
N39|14|A|1|1|Other disorders of urinary system
N40|14|M|1|1|Hyperplasia of prostate
N41|14|M|0|1|Inflammatory diseases of prostate
N42|14|M|0|1|Other disorders of prostate
N43|14|M|0|1|Hydrocele and spermatocele
N47|14|M|0|1|Redundant prepuce, phimosis and paraphimosis
N48|14|M|0|1|Other disorders of penis
N50|14|M|1|1|Other disorders of male genital organs
N73|14|F|1|1|Other female pelvic inflammatory diseases
N80|14|F|1|0|Endometriosis
N81|14|F|1|1|Female genital prolapse
N83|14|F|1|1|Noninflammatory disorders of ovary, Fallopian tube and broad ligament
N84|14|F|1|1|Polyp of female genital tract
N85|14|F|1|1|Other noninflammatory disorders of uterus, except cervix
N90|14|F|0|1|Other noninflammatory disorders of vulva and perineum
N92|14|F|1|1|Excessive, frequent and irregular menstruation
N93|14|F|1|0|Other abnormal uterine and vaginal bleeding
N94|14|F|1|0|Pain and other conditions associated with female genital organs and menstrual cycle
N95|14|F|1|1|Menopausal and other perimenopausal disorders",
    header = FALSE, sep = "|", quote = "", comment.char = "", stringsAsFactors = FALSE, col.names = c("icd3",
    "chapter_id", "scope_code", "prevalent", "incident", "disease"))
  data.frame(icd3 = fixed$icd3, disease = fixed$disease, category = chapter_labels[fixed$chapter_id], scope =
    unname(c(A = "All", M = "Male", F = "Female")[fixed$scope_code]), eligible_prevalent = as.logical(fixed$prevalent),
    eligible_incident = as.logical(fixed$incident), stringsAsFactors = FALSE)
}

# Public entry point. all_data contains the 2,920 prepared protein columns and
# the chronological-age label used by the primary model, with participant EIDs
# as row names. Supplied covariates must already have the manuscript's upstream
# median imputation; missing-value filtering below is retained from the source.
run_all_protein_comparator <- function(result, all_data, covariates, blood, outcome_dates, disease_catalog =
  sf_frozen_catalog(), config = all_protein_config()) {
  cfg <- config
  sf_require("data.table")
  Sys.setenv(OMP_NUM_THREADS = as.character(cfg$learner_threads), MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS =
    "1")
  data.table::setDTthreads(1L)
  set.seed(cfg$seed)
  ctx <- sf_prepare_comparator_context(result, all_data, cfg)
  general <- sf_select_general(ctx, cfg)
  oof <- sf_full_oof(ctx, general, cfg)
  clinical <- sf_prepare_outcomes(ctx, oof, cfg, covariates, blood, outcome_dates, disease_catalog)
  list(general_model = general, full_oof = oof, prevalent = sf_run_disease(clinical, "prevalent", cfg), incident =
    sf_run_disease(clinical, "incident", cfg), mortality = sf_run_mortality_comparison(clinical, cfg), config =
    cfg)
}


# ==============================================================================
# 08 | External validation
# ==============================================================================

# External validation in the Gong and Chia cohorts.
# Sources: the two-cohort full-training pipeline, all-age evaluation, and the
# later Chia-internal calibration and frozen-contribution analysis.
# All source helper names are local to this function to avoid collisions.
#
# external_cohorts contains named entries imm_of_aging and chia_als. Each
# entry has data (first column label; sample IDs in row names), aligned meta,
# and optional protein_map. The common panel is defined without external age.
# annotation: 58 pathway rows with reactome_id, pathway_name, class_id, class_name.
# ukb_contribution: original 58-pathway UKB contributions; the source checks
# its original held-out reference sample size of 10,602 and leaves it unchanged.
# Chia disease models use internal age correction and internal residual SD,
# superseding the earlier training-OOF-based Chia age-acceleration definition.
# This function returns analytical objects and does not create plots/workbooks.
run_external_validation <- function(all_data, immlist, L1, L2, selected_pipeline, annotation, ukb_contribution,
  external_cohorts, S = 666L, KIN = 5L, NC = 1L, TLIM = 600, FUS = c("EqualWeight", "MAEWeight"), POST_BOOT =
  2000L) {

  # Reuse the learner combination selected by the main UKB analysis. The shared
  # protein panel is refitted with that pair; it does not trigger model selection.
  if (!is.list(selected_pipeline) || !all(c("stage1", "stage2", "pipeline") %in% names(selected_pipeline)) ||
    any(!vapply(selected_pipeline[c("stage1", "stage2", "pipeline")], function(x) is.character(x) && length(x) ==
    1L && !is.na(x) && nzchar(x), logical(1)))) stop("selected_pipeline must contain the main analysis stage1, stage2 and pipeline names.")
  BEST1 <- selected_pipeline$stage1;
  BEST2 <- selected_pipeline$stage2;
  key <- selected_pipeline$pipeline
  if (!BEST1 %in% names(L1) || (!BEST2 %in% FUS && !BEST2 %in% names(L2))) stop("The main selected learner pair is absent from the supplied learner specifications.")
  L1 <- L1[BEST1]
  L2 <- if (BEST2 %in% FUS) list() else L2[BEST2]
  if (!is.list(external_cohorts) || length(external_cohorts) != 2L || anyDuplicated(names(external_cohorts)) ||
    !setequal(names(external_cohorts), c("imm_of_aging", "chia_als"))) stop("external_cohorts must contain imm_of_aging and chia_als.")

  # Keep the source cohort order because prediction seeds include cohort position.
  external_cohorts <- external_cohorts[c("imm_of_aging", "chia_als")]

  # Common-protein selection, full-training refit, and frozen predictions
  EXTERNAL_FEATURE_MODE <- "common_panel_two_cohorts"
  for (nm in c("all_data", "immlist", "L1", "L2")) if (!exists(nm, inherits = TRUE)) stop("Required training/configuration object is missing: ",
    nm)
  required_pkgs <- c("data.table", "mlr3", "callr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(mlr3))
  suppressPackageStartupMessages(library(parallel))
  if (requireNamespace("mlr3learners", quietly = TRUE)) suppressPackageStartupMessages(library(mlr3learners))
  if (requireNamespace("mlr3extralearners", quietly = TRUE)) suppressPackageStartupMessages(library(mlr3extralearners))
  if (!is.data.frame(all_data) || is.data.table(all_data)) stop("all_data must be a data.frame retaining participant IDs as row names.")
  if (!is.list(immlist) || !length(immlist) || is.null(names(immlist))) stop("immlist must be a named list of modules.")
  if (!is.list(L2)) stop("L2 must be a named learner list; use list() for fixed fusions only.")
  alias_table <- data.frame(old = c("FGFR1OP", "LRMP", "CD97", "DARS", "CASC4", "WISP2", "HARS", "FOPNL",
    "WISP1"), new = c("CEP43", "IRAG2", "ADGRE5", "DARS1", "GOLM2", "CCN5", "HARS1", "CEP20", "CCN4"))
  feature_key <- function(x) {
    z <- toupper(trimws(as.character(x)))
    j <- match(z, alias_table$old);
    hit <- which(!is.na(j))
    z[hit] <- alias_table$new[j[hit]]
    z
  }
  read_external <- function(cohort, dataset) {
    e <- cohort
    if (!is.list(e) || !all(c("data", "meta") %in% names(e))) stop("Each external cohort object must contain data and meta.")
    d <- as.data.frame(e$data, check.names = FALSE);
    m <- as.data.frame(e$meta, check.names = FALSE)
    if (!nrow(d) || !identical(names(d)[1L], "label") || anyDuplicated(names(d))) stop("data must be nonempty, have label as its first column, and have unique column names.")
    ids <- rownames(d)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("External data row names are invalid.")
    if (anyDuplicated(rownames(m)) || !setequal(ids, rownames(m))) stop("External data and meta sample IDs do not agree.")
    m <- m[match(ids, rownames(m)), , drop = FALSE]
    if ("sample_id" %in% names(m) && !identical(as.character(m$sample_id), ids)) stop("meta$sample_id does not agree with row names.")
    list(data = d, meta = m, protein_map = if ("protein_map" %in% names(e)) e$protein_map else NULL)
  }
  external_input <- setNames(lapply(names(external_cohorts), function(nm) tryCatch(read_external(external_cohorts[[nm]],
    nm), error = identity)), names(external_cohorts))
  input_status <- rbindlist(lapply(names(external_input), function(nm) {
    z <- external_input[[nm]]
    data.table(dataset = nm, read_ok = !inherits(z, "error"), message = if (inherits(z, "error")) conditionMessage(z) else "",
      input = nm)
  }))
  if (!all(input_status$read_ok)) stop("Both external cohort objects must be valid to define the three-way protein intersection; inspect input_status.")
  if (is.null(names(immlist)) || anyNA(names(immlist)) || any(!nzchar(trimws(names(immlist)))) || anyDuplicated(names(immlist))) stop("immlist module names must be nonempty and unique.")
  if (!nrow(all_data) || !identical(names(all_data)[1L], "label") || anyNA(names(all_data)) || anyDuplicated(names(all_data))) stop("all_data must contain samples, label in its first column, and valid unique column names.")
  original_immlist <- immlist
  training_columns <- names(all_data)[ - 1L]
  training_keys <- feature_key(training_columns)
  train_annotation_map <- rbindlist(lapply(names(immlist), function(m) {
    fs <- unique(as.character(immlist[[m]]));
    fs <- fs[!is.na(fs) & nzchar(trimws(fs))]
    rbindlist(lapply(fs, function(f) {
      j <- which(training_keys == feature_key(f))
      data.table(module = m, immlist_feature = f, symbol = feature_key(f), training_column = if (length(j) ==
        1L) training_columns[j] else NA_character_, status = if (!length(j)) "absent_from_training" else if (length(j) >
        1L) "ambiguous_training_symbol" else "matched")
    }))
  }), fill = TRUE)
  if (!nrow(train_annotation_map)) stop("immlist contains no valid protein annotations.")
  if (any(train_annotation_map$status == "ambiguous_training_symbol")) stop("A training symbol maps to multiple columns; inspect train_annotation_map and verify the assays.")
  original_modules <- setNames(lapply(names(immlist), function(m) unique(train_annotation_map[module == m &
    status == "matched", training_column])), names(immlist))
  original_features <- unique(unlist(original_modules, use.names = FALSE))
  if (!length(original_features)) stop("No proteins match between immlist and the training data.")
  if (anyDuplicated(feature_key(original_features))) stop("Training protein symbols are duplicated after mapping.")
  coverage_one <- function(z, nm, features) {
    if (inherits(z, "error")) return(data.table(dataset = nm, feature = features, source_column = NA_character_,
      status = "file_unavailable", n_observed = NA_integer_, n_samples = NA_integer_))
    cols <- setdiff(names(z$data), "label");
    keys <- feature_key(cols)
    rbindlist(lapply(features, function(f) {
      j <- which(keys == feature_key(f));
      status <- if (!length(j)) "missing_column" else if (length(j) > 1L) "ambiguous_column" else "matched"
      nobs <- NA_integer_;
      source <- if (length(j) == 1L) cols[j] else NA_character_
      if (status == "matched") {
        v <- z$data[[source]]
        if (!is.numeric(v)) status <- "nonnumeric_column" else {
          nobs <- sum(is.finite(v))
          if (any(!is.na(v) & !is.finite(v))) status <- "nonfinite_column"
          else if (nobs == 0L) status <- "all_missing_column"
        }
      }
      data.table(dataset = nm, feature = f, source_column = source, status = status, n_observed = nobs, n_samples =
        nrow(z$data))
    }))
  }
  feature_coverage <- rbindlist(lapply(names(external_input), function(nm) coverage_one(external_input[[nm]],
    nm, original_features)))
  common <- original_features[vapply(original_features, function(f) all(feature_coverage[feature == f, status] ==
    "matched"), logical(1))]
  model_immlist <- lapply(original_modules, intersect, common)
  common_proteins <- data.table(training_feature = common, symbol = feature_key(common))
  for (ds in names(external_input)) {
    tab <- feature_coverage[dataset == ds]
    common_proteins[[paste0(ds, "_column")]] <- tab$source_column[match(common, tab$feature)]
  }
  feature_coverage[, retained_for_model := feature %in% common]
  module_coverage <- rbindlist(lapply(names(original_modules), function(m) {
    data.table(module = m, original_features = length(original_modules[[m]]), retained_features = length(model_immlist[[m]]),
      excluded_for_external = paste(setdiff(original_modules[[m]], model_immlist[[m]]), collapse = ";"))
  }))
  if (any(!lengths(model_immlist))) stop("No usable common proteins remain in these modules: ", paste(names(model_immlist)[!lengths(model_immlist)],
    collapse = ", "), "; inspect module_coverage. Empty modules cannot be silently removed.")
  message("Usable proteins shared by training data and both external cohorts: ", length(common), " / ", length(original_features),
    "; excluded ", length(setdiff(original_features, common)), " proteins; all training participants retained: ",
    nrow(all_data), " participants.")
  print(module_coverage)

  ## ---------- audit, packages and mapped modules ----------
  normL <- function(x, s) {
    if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) || any(!nzchar(trimws(names(x)))) ||
      anyDuplicated(names(x))) stop(s, " must be a uniquely named non-empty list");
    z <- lapply(x, \(q) tryCatch({
      if (is.character(q)) {
        if (length(q) != 1L) stop("learner ID must have length 1");
        q <- lrn(if (startsWith(q, "regr.")) q else paste0("regr.", q))
      } else if (is.list(q) && !inherits(q, "Learner")) {
        id <- q$id;
        q$id <- NULL;
        if (!is.character(id) || length(id) != 1L) stop("learner specification needs one id");
        id <- if (startsWith(id, "regr.")) id else paste0("regr.", id);
        q <- do.call(lrn, c(list(id), q))
      };
      if (!inherits(q, "Learner") || q$task_type != "regr") stop("not a regression learner");
      q
    }, error = identity));
    ok <- !vapply(z, inherits, logical(1), "error");
    st <- data.table(name = names(x), registered = ok, registry_error = vapply(z, \(q) if (inherits(q, "error")) conditionMessage(q) else "",
      ""));
    if (any(!ok)) warning("Skipping unregistered ", s, ": ", paste(st$name[!ok], collapse = ", "));
    z <- z[ok];
    if (!length(z)) stop("No registered learner in ", s);
    list(learners = z, status = st)
  }
  cap <- function(z) {
    fb <- lrn("regr.featureless");
    if ("when" %in% names(formals(z$encapsulate))) z$encapsulate("callr", fallback = fb, when = \(cond, stage) FALSE) else z$encapsulate("callr",
      fallback = fb);
    z$timeout <- c(train = TLIM, predict = TLIM);
    z
  }
  a1 <- normL(L1, "L1");
  a2 <- if (length(L2)) normL(L2, "L2") else list(learners = list(), status = data.table(name = character(),
    registered = logical(), registry_error = character()));
  L1 <- lapply(a1$learners, cap);
  L2 <- lapply(a2$learners, cap);
  l1_registry_status <- a1$status;
  l2_registry_status <- a2$status;
  setnames(l1_registry_status, "name", "stage1");
  setnames(l2_registry_status, "name", "stage2");
  rm(a1, a2)
  if (any(grepl("__", c(names(L1), names(L2)), fixed = TRUE)) || any(FUS %in% names(L2))) stop("Learner names cannot contain '__' or duplicate mandatory fusion names")
  stopifnot(KIN >= 2L, NC >= 1L, TLIM > 0, names(all_data)[1L] == "label", is.numeric(all_data$label), !anyNA(names(all_data)),
    !anyDuplicated(names(all_data)), is.list(immlist), length(immlist) > 0L)
  if (is.null(names(immlist)) || anyNA(names(immlist)) || any(!nzchar(trimws(names(immlist)))) || anyDuplicated(names(immlist))) stop("Module names are missing, blank or duplicated")
  eid <- rownames(all_data);
  if (is.null(eid) || anyNA(eid) || any(!nzchar(eid)) || anyDuplicated(eid)) stop("EID is missing, blank or duplicated")
  ip <- rownames(installed.packages());
  if (!"callr" %in% ip) stop("Missing mandatory package: callr")
  l1miss <- setNames(lapply(L1, \(z) setdiff(z$packages, ip)), names(L1));
  l1_package_status <- data.table(stage1 = names(L1), available = !lengths(l1miss), missing = vapply(l1miss,
    paste, "", collapse = ", "));
  badL1 <- l1_package_status$stage1[!l1_package_status$available];
  if (length(badL1)) warning("Skipping unavailable L1: ", paste(badL1, collapse = ", "));
  L1 <- L1[l1_package_status$stage1[l1_package_status$available]];
  if (!length(L1)) stop("No available L1 learner")
  l2miss <- setNames(lapply(L2, \(z) setdiff(z$packages, ip)), names(L2));
  l2_package_status <- data.table(stage2 = names(L2), available = !lengths(l2miss), missing = vapply(l2miss,
    paste, "", collapse = ", "));
  badL2 <- l2_package_status$stage2[!l2_package_status$available];
  if (length(badL2)) warning("Skipping unavailable L2: ", paste(badL2, collapse = ", "));
  L2 <- L2[l2_package_status$stage2[l2_package_status$available]]
  ann <- lapply(model_immlist, \(z) unique(as.character(z)));
  mods0 <- lapply(ann, intersect, names(all_data)[ - 1L]);
  unmapped <- Map(setdiff, ann, mods0);
  if (any(!lengths(mods0))) stop("No mapped feature: ", paste(names(mods0)[!lengths(mods0)], collapse = ", "))
  used <- unique(unlist(mods0, use.names = FALSE));
  fi <- sprintf("F%05d", seq_along(used));
  mi <- sprintf("M%03d", seq_along(mods0));
  f2i <- setNames(fi, used);
  mods <- lapply(mods0, \(z) unname(f2i[z]));
  names(mods) <- mi
  dat <- all_data[, c("label", used), drop = FALSE];
  names(dat) <- c("label", fi);
  if (!all(vapply(dat, is.numeric, logical(1))) || !all(vapply(dat, \(z) all(is.finite(z)), logical(1)))) stop("label/features must be finite numeric values")
  constant <- used[vapply(dat[ - 1L], \(z) length(unique(z)) < 2L, logical(1))];
  if (length(constant)) warning("Constant features retained because feature selection is disabled: ", paste(constant,
    collapse = ", "))
  feature_name_map <- data.table(original = used, internal = fi);
  module_name_map <- data.table(original = names(mods0), internal = mi);
  module_feature_map <- rbindlist(Map(\(m, o, i) data.table(module = m, original = o, internal = i), names(mods0),
    mods0, mods));
  reserved <- c("eid", "label", "fold", "prediction", "immune_age", "age_gap_raw", "age_acceleration", "age_acceleration_z",
    "sample_id", "dataset", "eligible_age_validation");
  if (any(names(mods0) %in% reserved)) stop("Module names conflict with output columns: ", paste(intersect(names(mods0),
    reserved), collapse = ", "))
  mapping <- data.table(module = names(mods0), module_internal = mi, annotated = lengths(ann), mapped = lengths(mods0),
    excluded_unmapped = lengths(unmapped));
  y <- dat$label;
  n <- length(y)
  tasks <- Map(\(z, i) TaskRegr$new(paste0("module", i), as.data.table(dat[, c("label", z), drop = FALSE]),
    target = "label"), mods, seq_along(mods));
  names(tasks) <- names(mods);
  stopifnot(all(vapply(tasks, \(z) identical(z$row_ids, seq_len(n)), logical(1))))

  ## ---------- inner folds and fit helpers ----------
  mkfold <- function(a, k, seed, b = 10L) {
    if (k > length(a)) stop("folds > samples");
    set.seed(seed);
    q <- unique(quantile(a, seq(0, 1, length.out = b + 1L), type = 8));
    g <- if (length(q) > 2L) cut(a, q, include.lowest = TRUE) else factor(rep(1L, length(a)));
    f <- integer(length(a));
    for (ii in split(seq_along(a), g, drop = TRUE)) {
      v <- rep_len(seq_len(k), length(ii));
      f[ii] <- v[sample.int(length(ii))]
    };
    if (any(tabulate(f, k) == 0L)) {
      v <- rep_len(seq_len(k), length(a));
      f <- v[sample.int(length(a))]
    };
    f
  }
  restoreM <- function(M) {
    j <- match(colnames(M), module_name_map$internal);
    if (anyNA(j)) stop("Unknown internal module name");
    colnames(M) <- module_name_map$original[j];
    M
  }
  MAP <- function(x, f, ...) {
    if (NC == 1L || .Platform$OS.type == "windows") lapply(x, f, ...) else mclapply(x, f, ..., mc.cores =
      min(NC, length(x)), mc.set.seed = TRUE, mc.preschedule = FALSE)
  }
  badfit <- function(l) {
    e <- tryCatch(l$errors, error = \(e) NULL);
    if ((is.numeric(e) && any(e > 0)) || (!is.null(e) && !is.numeric(e) && length(e) > 0L)) return(TRUE);
    z <- tryCatch(l$log, error = \(e) NULL);
    if (is.null(z) || !NROW(z)) return(FALSE);
    if ("class" %in% names(z) && any(tolower(as.character(z$class)) == "error")) return(TRUE);
    x <- if ("msg" %in% names(z)) z$msg else if ("condition" %in% names(z)) vapply(z$condition, \(q) if (inherits(q,
      "condition")) conditionMessage(q) else paste(as.character(q), collapse = " "), "") else unlist(z, use.names =
      FALSE);
    any(grepl("using fallback learner", paste(x, collapse = " "), ignore.case = TRUE))
  }
  chk <- function(l) {
    if (is.null(l$model) || badfit(l)) stop("Learner failed or timed out: ", l$id, call. = FALSE);
    invisible(l)
  };
  trn <- function(l, t, id = NULL) {
    l$train(t, row_ids = id);
    chk(l)
  }
  pp <- function(l, t, id) {
    q <- l$predict(t, row_ids = id);
    chk(l);
    j <- match(id, q$row_ids);
    if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Invalid/misaligned stage-1 row_ids");
    p <- q$response[j];
    if (length(p) != length(id) || any(!is.finite(p))) stop("Invalid stage-1 predictions");
    p
  }

  ## ---------- stage 1: one shared learner is fitted separately to every module ----------
  oof1 <- function(task, tr, base, inner, seed) {
    p <- numeric(length(tr));
    for (v in seq_len(KIN)) {
      va <- which(inner == v);
      l <- base$clone(deep = TRUE);
      l$reset();
      set.seed(seed + v);
      trn(l, task, tr[inner != v]);
      p[va] <- pp(l, task, tr[va])
    };
    if (any(!is.finite(p))) stop("Invalid inner OOF predictions");
    p
  }

  ## ---------- selected stage-2 fusion helper ----------
  fuse2 <- function(z0, z1, a, method) {
    if (!identical(colnames(z0), colnames(z1)) || any(!is.finite(z0)) || any(!is.finite(z1))) stop("Invalid/misaligned module predictions");
    e <- colMeans(abs(sweep(z0, 1L, a, "-")));
    if (method == "EqualWeight") w <- rep(1 / ncol(z0), ncol(z0)) else {
      zero <- e <= sqrt(.Machine$double.eps);
      if (any(zero)) w <- zero / sum(zero) else {
        w <- 1 / (e + sqrt(.Machine$double.eps) * max(1, median(e)));
        w <- w / sum(w)
      }
    };
    names(w) <- colnames(z0);
    list(pred = drop(z1[, names(w), drop = FALSE] %*% w), weights = data.table(module_internal = names(w),
      module_MAE = e[names(w)], weight = as.numeric(w)))
  }

  ## ---------- final full-training fit; L2 uses inner-OOF L1 predictions ----------
  # This function never receives external data or external ages.
  refit_full <- function(stage1, stage2, seed) {
    inn <- mkfold(y, KIN, seed)
    z <- MAP(seq_along(tasks), function(m) tryCatch({
      p0 <- oof1(tasks[[m]], ids, L1[[stage1]], inn, seed + 1000L * m)
      l <- L1[[stage1]]$clone(deep = TRUE);
      l$reset();
      set.seed(seed + 90000L + m)
      trn(l, tasks[[m]], ids)
      list(oof = p0, model = l)
    }, error = identity))
    failed <- which(vapply(z, inherits, logical(1), "error"))
    if (length(failed)) stop("Full-training module ", names(tasks)[failed[1L]], ": ", conditionMessage(z[[failed[1L]]]))
    z0 <- do.call(cbind, lapply(z, `[[`, "oof"));
    colnames(z0) <- names(tasks)
    if (any(!is.finite(z0))) stop("Invalid full-training inner OOF predictions")
    if (stage2 %in% FUS) {
      q <- fuse2(z0, z0, y, stage2)
      model2 <- list(kind = "fusion", method = stage2, weights = q$weights)
    } else {
      t <- TaskRegr$new("stage2_full_training", as.data.table(cbind(label = y, z0)), target = "label")
      l <- L2[[stage2]]$clone(deep = TRUE);
      l$reset();
      set.seed(seed + 50000L);
      trn(l, t)
      model2 <- list(kind = "learner", learner = l)
    }
    list(stage1 = setNames(lapply(z, `[[`, "model"), names(tasks)), stage2 = model2, inner_oof = z0, inner_fold =
      inn, n_training = length(ids), stage1_name = stage1, stage2_name = stage2)
  }
  predict_new_checked <- function(l, nd) {
    q <- l$predict_newdata(as.data.table(nd));
    chk(l)
    j <- match(seq_len(nrow(nd)), q$row_ids)
    if (anyNA(j) || anyDuplicated(q$row_ids)) stop("External prediction row_ids are missing/duplicated")
    p <- q$response[j]
    if (length(p) != nrow(nd) || any(!is.finite(p))) stop("Invalid external predictions")
    p
  }
  predict_frozen <- function(model, x, seed) {
    M <- matrix(NA_real_, nrow(x), length(mods), dimnames = list(rownames(x), names(mods)))
    for (m in seq_along(mods)) {
      set.seed(seed + m)
      nd <- x[, tasks[[m]]$feature_names, drop = FALSE]
      M[, m] <- predict_new_checked(model$stage1[[m]], nd)
    }
    if (model$stage2$kind == "fusion") {
      w <- model$stage2$weights
      p <- drop(M[, w$module_internal, drop = FALSE] %*% w$weight)
    } else {
      set.seed(seed + 10000L)
      p <- predict_new_checked(model$stage2$learner, as.data.frame(M, check.names = FALSE))
    }
    list(pred = p, module = M)
  }
  age_external <- function(z) {
    raw <- trimws(as.character(z$data$label));
    exact <- !is.na(raw) & grepl("^[0-9]+(\\.[0-9]+)?$", raw)
    a <- rep(NA_real_, length(raw));
    a[exact] <- as.numeric(raw[exact])
    ok <- exact & is.finite(a) & a >= 0 & a <= 120
    if ("age_status" %in% names(z$meta)) ok <- ok & !is.na(z$meta$age_status) & z$meta$age_status == "exact"
    if ("eligible_age_validation" %in% names(z$meta)) ok <- ok & tolower(as.character(z$meta$eligible_age_validation)) %in%
    c("true", "1")
    ok[is.na(ok)] <- FALSE;
    a[!ok] <- NA_real_
    list(age = a, eligible = ok)
  }

  ## ---------- refit the main selected pair on all development participants ----------
  ids <- seq_len(n)
  if (n < KIN) stop("The training sample is too small for the specified inner folds.")
  if (!BEST1 %in% names(L1) || (!BEST2 %in% FUS && !BEST2 %in% names(L2))) stop("The main selected learner pair is unavailable; no alternative is substituted.")
  message("Refitting the main selected pipeline: ", BEST1, " + ", BEST2)
  final_fit <- tryCatch(refit_full(BEST1, BEST2, S + 9100L), error = identity)
  if (inherits(final_fit, "error")) stop("The selected pipeline failed during the full-training refit; no alternative was substituted: ",
    conditionMessage(final_fit))
  training_medians <- vapply(all_data[, used, drop = FALSE], median, numeric(1), na.rm = TRUE)
  if (any(!is.finite(training_medians))) stop("Invalid training medians")
  model_bundle <- list(schema_version = "full_training_two_external_common_v2", run_id = "external_validation",
    best = list(stage1 = BEST1, stage2 = BEST2, pipeline = key), models = final_fit, training_eid = eid, n_training =
    n, feature_name_map = feature_name_map, module_name_map = module_name_map, module_features = mods0, module_feature_map =
    module_feature_map, training_medians = training_medians, alias_table = alias_table, external_cohorts =
    names(external_cohorts), common_proteins = common_proteins, training_feature_mapping = train_annotation_map,
    feature_mode = EXTERNAL_FEATURE_MODE, imputation = "training data already finite; external sporadic NA use fixed training medians",
    input_training_already_imputed = "Upstream preprocessing was not provided and cannot be audited by this script",
    session = sessionInfo())

  ## ---------- external cohorts: frozen model, independent prediction ----------
  external_results <- setNames(vector("list", length(external_input)), names(external_input))
  external_status <- list()
  for (di in seq_along(external_input)) {
    nm <- names(external_input)[di];
    z <- external_input[[di]]
    message("External evaluation: ", nm)
    if (inherits(z, "error")) {
      external_status[[di]] <- data.table(dataset = nm, status = "input_error", n_samples = NA_integer_, n_evaluated =
        NA_integer_, reason = conditionMessage(z))
      next
    }
    cov <- feature_coverage[dataset == nm & feature %in% used]
    cov <- cov[match(used, feature)]
    bad <- cov[status != "matched"]
    if (nrow(bad)) {
      external_status[[di]] <- data.table(dataset = nm, status = "incomplete_features", n_samples = nrow(z$data),
        n_evaluated = 0L, reason = paste0(nrow(bad), " required features unavailable; inspect feature_coverage"))
      external_results[[nm]] <- list(status = "incomplete_features", feature_coverage = cov, meta = z$meta)
      message("  Full-model evaluation could not proceed: ", nrow(bad), " required features unavailable; inspect feature_coverage.")
      next
    }
    evaluated <- tryCatch({
      x <- z$data[, cov$source_column, drop = FALSE];
      names(x) <- used
      missing_counts <- colSums(is.na(x))
      for (j in seq_along(x)) {
        v <- as.numeric(x[[j]]);
        v[is.na(v)] <- training_medians[j];
        x[[j]] <- v
      }
      names(x) <- fi;
      rownames(x) <- rownames(z$data)
      if (any(!is.finite(as.matrix(x)))) stop("External input is not finite after fixed imputation")

      # No external label is passed to predict_frozen.
      q <- predict_frozen(final_fit, x, S + 1000000L + 10000L * di)
      ag <- age_external(z);
      a <- ag$age
      M <- restoreM(q$module)
      pred <- data.frame(sample_id = rownames(z$data), dataset = nm, label = a, eligible_age_validation =
        ag$eligible, prediction = q$pred, immune_age = q$pred, M, check.names = FALSE)
      rownames(pred) <- rownames(z$data)
      stopifnot(identical(rownames(pred), rownames(z$meta)))
      list(status = "ok", prediction = pred, meta = z$meta, feature_coverage = cov, imputation_counts = data.table(feature =
        used, n_imputed = as.integer(missing_counts)), protein_map = z$protein_map)
    }, error = identity)
    if (inherits(evaluated, "error")) {
      external_status[[di]] <- data.table(dataset = nm, status = "prediction_or_evaluation_error", n_samples =
        nrow(z$data), n_evaluated = 0L, reason = conditionMessage(evaluated))
      external_results[[nm]] <- list(status = "prediction_or_evaluation_error", message = conditionMessage(evaluated))
      message("  Evaluation failure: ", conditionMessage(evaluated));
      next
    }
    external_results[[nm]] <- evaluated
    external_status[[di]] <- data.table(dataset = nm, status = "ok", n_samples = nrow(z$data), n_evaluated =
      sum(evaluated$prediction$eligible_age_validation), reason = "")
  }
  external_status <- rbindlist(external_status, fill = TRUE)
  desc <- function(x) lapply(x, function(z) list(id = z$id, params = z$param_set$values, packages = z$packages,
    timeout = z$timeout, encapsulation = z$encapsulation))
  result <- list(run_id = "external_validation", best = list(stage1 = BEST1, stage2 = BEST2, pipeline = key),
    n_training = n, training_eid = eid, mapping = mapping, module_coverage = module_coverage, excluded_unmapped =
    unmapped, common_proteins = common_proteins, training_feature_mapping = train_annotation_map, feature_name_map =
    feature_name_map, module_name_map = module_name_map, module_feature_map = module_feature_map, module_features =
    mods0, registry_status = list(stage1 = l1_registry_status, stage2 = l2_registry_status), package_status =
    list(stage1 = l1_package_status, stage2 = l2_package_status), external_status = external_status, external =
    external_results, feature_coverage = feature_coverage, training_medians = training_medians, training_age_range =
    range(y), config = list(seed = S, training_fraction = 1, inner_folds = KIN, cores = NC, per_fit_timeout_seconds =
    TLIM, metrics = c("MAE", "RMSE", "PCC"), bootstrap = POST_BOOT, feature_mode = EXTERNAL_FEATURE_MODE,
    shared_stage1 = TRUE, final_stage2_input = "full-development inner-OOF stage1 predictions", final_models_retained =
    TRUE, stage1_selection = "locked to the main UKB selected learner pair", external_model_selection = "none",
    external_imputation = "fixed training medians for sporadic NA only; no missing-column imputation", external_scaling =
    "none added; source NPX unchanged, model-internal learned transformations remain frozen", calibration =
    "Chia internal pooled exact-age calibration follows frozen prediction; external accuracy uses raw prediction",
    source_training_imputation = "input required finite; upstream preprocessing not supplied", external_cohorts =
    names(external_cohorts), feature_definition = "immlist proteins present in training and usable in both external cohorts",
    selection_failure_policy = "stop if the main selected pair cannot refit; no switch to another algorithm",
    model_identity = "two-cohort common-protein model; retrained using all development participants", original_module_features =
    original_modules), learners = list(stage1 = desc(L1), stage2 = desc(L2)), session = sessionInfo())
  print(result$best);
  print(external_status)

  # All-age evaluation and long-form prediction records
  post_bind <- function(xs) {
    xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
    if (!length(xs)) return(data.frame())
    all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
    xs <- lapply(xs, function(x) {
      for (k in setdiff(all_names, names(x))) x[[k]] <- NA;
      x[, all_names, drop = FALSE]
    })
    out <- do.call(rbind, xs);
    rownames(out) <- NULL;
    out
  }
  post_bool <- function(x) !is.na(x) & tolower(trimws(as.character(x))) %in% c("true", "1")
  post_diag <- function(x) {
    s <- tolower(gsub("[^[:alnum:]]", "", trimws(as.character(x))))
    ans <- rep(NA_character_, length(s))
    ans[s %in% c("healthycontrol", "healthycontrols", "healthy")] <- "Healthy_control"
    ans[s %in% c("als", "amyotrophiclateralsclerosis")] <- "ALS"
    ans[s %in% c("neurologicalcontrol", "neurologicalcontrols")] <- "Neurological_control"
    ans
  }
  post_metrics <- function(a, P, unit, B, seed) {
    P <- as.matrix(P);
    storage.mode(P) <- "double";
    cn <- colnames(P);
    ks <- c("MAE", "RMSE", "PCC")
    one <- function(y, X) {
      e <- sweep(X, 1L, y, "-");
      yc <- y - mean(y);
      xc <- sweep(X, 2L, colMeans(X), "-")
      ss <- sum(yc ^ 2);
      den <- sqrt(ss * colSums(xc ^ 2))
      pcc <- rep(NA_real_, ncol(X));
      ok <- is.finite(den) & den > 0
      pcc[ok] <- drop(crossprod(yc, xc))[ok] / den[ok]
      cbind(MAE = colMeans(abs(e)), RMSE = sqrt(colMeans(e ^ 2)), PCC = pcc)
    }
    if (length(a)) stopifnot(nrow(P) == length(a), all(is.finite(a)), all(is.finite(P)), length(unit) ==
      length(a))
    point <- if (length(a)) one(a, P) else matrix(NA_real_, ncol(P), 3L, dimnames = list(NULL, ks))
    groups <- split(seq_along(a), unit);
    ng <- length(groups)
    boot <- array(NA_real_, dim = c(B, ncol(P), 3L))
    if (B > 0L && ng >= 2L) {
      set.seed(seed)
      for (b in seq_len(B)) {
        ix <- unlist(groups[sample.int(ng, ng, replace = TRUE)], use.names = FALSE)
        boot[b, , ] <- one(a[ix], P[ix, , drop = FALSE])
      }
    }
    post_bind(lapply(seq_len(ncol(P)), function(j) post_bind(lapply(seq_along(ks), function(k) {
      vals <- boot[, j, k];
      vals <- vals[is.finite(vals)]
      ci <- if (length(vals) > 1L) quantile(vals, c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
      data.frame(clock_id = cn[j], metric = ks[k], estimate = point[j, k], lower95 = ci[1L], upper95 = ci[2L],
        n_valid_bootstrap = length(vals), N = length(a), n_units = ng, stringsAsFactors = FALSE)
    }))))
  }
  postprocess_external <- function(r, B = 2000L, seed = 666L) {
    if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 0 || B != as.integer(B)) stop("POST_BOOT must be a nonnegative integer.")
    if (!is.list(r) || !setequal(r$config$external_cohorts, c("imm_of_aging", "chia_als"))) stop("The input must be the two-cohort common-protein model result.")
    expected <- c("imm_of_aging", "chia_als")
    if (!all(expected %in% names(r$external)) || any(!vapply(r$external[expected], function(z) identical(z$status,
      "ok"), logical(1)))) stop("Both external cohorts must have completed predictions.")
    pathways <- names(r$module_features)
    if (!length(pathways) || anyNA(pathways) || anyDuplicated(pathways)) stop("Invalid pathway identities.")
    train_range <- r$training_age_range
    clock_map <- data.frame(clock_id = c("overall", sprintf("P%03d", seq_along(pathways))), clock_name =
      c("Overall immune age", pathways), clock_type = c("overall", rep("pathway", length(pathways))), prediction_column =
      c("immune_age", pathways), stringsAsFactors = FALSE)
    perf <- audit <- external <- list()
    for (di in seq_along(expected)) {
      nm <- expected[di];
      z <- r$external[[nm]];
      p <- as.data.frame(z$prediction, check.names = FALSE);
      m <- as.data.frame(z$meta, check.names = FALSE)
      need <- c("sample_id", "label", clock_map$prediction_column)
      if (!all(need %in% names(p)) || anyDuplicated(names(p)) || anyNA(p$sample_id) || anyDuplicated(p$sample_id)) stop(nm,
        ": invalid prediction sample IDs or pathway columns.")
      sid <- as.character(p$sample_id)
      if (!setequal(sid, rownames(m)) || anyDuplicated(rownames(m))) stop(nm, ": meta and prediction samples do not agree.")
      m <- m[match(sid, rownames(m)), , drop = FALSE]
      if ("sample_id" %in% names(m) && !identical(as.character(m$sample_id), sid)) stop(nm, ": meta$sample_id does not agree with row names.")
      if (is.numeric(p$label)) {
        age <- as.numeric(p$label);
        exact <- is.finite(age)
      } else {
        age_text <- trimws(as.character(p$label));
        exact <- !is.na(age_text) & grepl("^[0-9]+(\\.[0-9]+)?$", age_text)
        age <- rep(NA_real_, nrow(p));
        age[exact] <- as.numeric(age_text[exact]);
        exact <- exact & is.finite(age)
      }
      if ("eligible_age_validation" %in% names(p)) exact <- exact & post_bool(p$eligible_age_validation)
      if ("age_status" %in% names(m)) exact <- exact & !is.na(m$age_status) & as.character(m$age_status) ==
      "exact"
      age[!exact] <- NA_real_;
      valid_age <- exact
      reason <- ifelse(valid_age, "included", "missing_or_nonexact_age")
      P <- as.matrix(p[, clock_map$prediction_column, drop = FALSE])
      if (!is.numeric(P) || any(!is.finite(P))) stop(nm, ": frozen predictions contain invalid values.")
      colnames(P) <- clock_map$clock_id
      has_subject <- rep(FALSE, nrow(p));
      subject <- rep(NA_character_, nrow(p))
      if ("subject_id" %in% names(m)) {
        subject <- trimws(as.character(m$subject_id));
        has_subject <- !is.na(subject) & nzchar(subject)
      }
      unit <- ifelse(has_subject, paste0("participant:", subject), paste0("sample:", sid))
      diag <- rep(NA_character_, nrow(p));
      dcols <- intersect(c("Diagnosis", "diagnosis"), names(m))
      if (length(dcols)) {
        ds <- lapply(dcols, function(k) post_diag(m[[k]]));
        diag <- ds[[1L]]
        if (length(ds) > 1L) for (d in ds[ - 1L]) {
          disagree <- !is.na(diag) & !is.na(d) & diag != d
          if (any(disagree)) stop(nm, ": Diagnosis and diagnosis columns disagree.")
          take <- is.na(diag);
          diag[take] <- d[take]
        }
      }
      p$label <- age;
      p$eligible_age_validation <- exact;
      p$age_range_filter_applied <- FALSE
      p$eligible_age_analysis <- valid_age;
      p$age_exclusion_reason <- reason
      p$training_age_min <- train_range[1L];
      p$training_age_max <- train_range[2L]
      long <- post_bind(lapply(seq_len(nrow(clock_map)), function(j) data.frame(sample_id = sid, dataset =
        nm, unit_id = unit, diagnosis = diag, age = age, clock_map[rep(j, length(sid)), c("clock_id", "clock_name",
        "clock_type"), drop = FALSE], predicted_age = P[, j], ageacc = NA_real_, ageacc_z = NA_real_, eligible_age_analysis =
        valid_age, age_exclusion_reason = reason, stringsAsFactors = FALSE, check.names = FALSE)))

      # Preserve the original group-specific bootstrap seeds (Gong index 1;
      # Chia healthy-control index 2), while computing only reported populations.
      groups <- if (nm == "imm_of_aging") list(all_samples = seq_len(nrow(p))) else list(Healthy_control =
        which(diag == "Healthy_control"))
      for (g in names(groups)) {
        gi <- if (nm == "imm_of_aging") 1L else 2L
        rows <- groups[[g]];
        keep <- rows[valid_age[rows]]
        q <- post_metrics(age[keep], P[keep, , drop = FALSE], unit[keep], B, seed + di * 10000L + gi)
        q$clock_name <- clock_map$clock_name[match(q$clock_id, clock_map$clock_id)]
        q$clock_type <- clock_map$clock_type[match(q$clock_id, clock_map$clock_id)]
        q$dataset <- nm;
        q$group <- g;
        q$n_predictions_group <- length(rows)
        q$n_missing_or_nonexact_age <- sum(!exact[rows])
        q$n_outside_training_range_included <- sum(exact[rows] & (age[rows] < train_range[1L] | age[rows] >
          train_range[2L]), na.rm = TRUE)
        q$age_range_filter_applied <- FALSE
        q$training_age_min <- train_range[1L];
        q$training_age_max <- train_range[2L]
        q$age_min_evaluated <- if (length(keep)) min(age[keep]) else NA_real_
        q$age_max_evaluated <- if (length(keep)) max(age[keep]) else NA_real_
        q$bootstrap_unit <- if (length(keep) && all(has_subject[keep])) "participant" else if (any(has_subject[keep])) "sample_or_mixed" else "sample"
        q$pipeline <- r$best$pipeline;
        q$analysis_set <- "all exact ages; no training age range restriction"
        perf[[paste(nm, g, sep = "::")]] <- q
      }
      audit[[nm]] <- data.frame(dataset = nm, n_prediction = nrow(p), n_exact_age = sum(exact), n_evaluated =
        sum(valid_age), n_missing_or_nonexact = sum(!exact), age_range_filter_applied = FALSE, n_below_training_min_included =
        sum(exact & age < train_range[1L], na.rm = TRUE), n_above_training_max_included = sum(exact & age >
        train_range[2L], na.rm = TRUE), n_with_subject_id = sum(has_subject), training_age_min = train_range[1L],
        training_age_max = train_range[2L], stringsAsFactors = FALSE)
      external[[nm]] <- list(prediction = p, meta = m, ageacc = long)
    }
    performance <- post_bind(perf);
    inclusion <- post_bind(audit)
    output <- list(schema_version = "external_all_ages_tube_v2", source_model_run = r$run_id, pipeline =
      r$best, training_age_range = train_range, clock_map = clock_map, external = external, performance =
      performance, inclusion = inclusion, config = list(bootstrap = B, seed = seed, covariates = c("age",
      "sex", "tube"), age_range = "unrestricted; all exact ages included", age_range_filter_applied = FALSE,
      ageacc = "Chia internal calibration is applied after collecting frozen predictions; no Gong age-acceleration analysis",
      unit = "sample unless a recorded subject_id enables participant clustering", downstream_model = "Chia diagnosis models are supplied by the subsequent internal-calibration analysis"),
      session = sessionInfo())
    invisible(output)
  }

  # External performance and cohort characteristics
  # External validation: reuse saved predictions and frozen age models.
  # This file defines functions only. It never sources or retrains the upstream clock.
  f7_bind <- function(x) {
    x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
    if (!length(x)) return(data.frame())
    nn <- unique(unlist(lapply(x, names), use.names = FALSE))
    x <- lapply(x, function(z) {
      for (k in setdiff(nn, names(z))) z[[k]] <- NA;
      z[, nn, drop = FALSE]
    })
    z <- do.call(rbind, x);
    rownames(z) <- NULL;
    z
  }
  f7_need <- function(x, columns, label) {
    miss <- setdiff(columns, names(x))
    if (length(miss)) stop(label, " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  f7_key <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
  f7_clean <- function(x) {
    z <- trimws(as.character(x))
    z[is.na(x) | f7_key(z) %in% c("", "na", "nan", "unknown", "missing", "notavailable", "notreported")] <-
    NA_character_
    z
  }
  f7_bool <- function(x) !is.na(x) & tolower(trimws(as.character(x))) %in% c("true", "1")
  f7_diag <- function(x) {
    z <- f7_key(f7_clean(x));
    ans <- rep(NA_character_, length(z))
    ans[z %in% c("healthy", "healthycontrol", "healthycontrols")] <- "Healthy_control"
    ans[z %in% c("als", "amyotrophiclateralsclerosis")] <- "ALS"
    ans[z %in% c("neurologicalcontrol", "neurologicalcontrols")] <- "Neurological_control"
    ans
  }

  # Match stable IDs or exact full names only; never use approximate text matching.
  f7_match_pathway <- function(x, annotation) {
    x <- as.character(x);
    idx <- match(x, annotation$pathway_name)
    direct <- match(x, annotation$reactome_id);
    idx[is.na(idx)] <- direct[is.na(idx)]
    has_id <- grepl("R-HSA-[0-9]+", x)
    id <- rep(NA_character_, length(x))
    id[has_id] <- sub(".*?(R-HSA-[0-9]+).*", "\\1", x[has_id], perl = TRUE)
    parsed <- match(id, annotation$reactome_id);
    idx[is.na(idx)] <- parsed[is.na(idx)]
    if (anyNA(idx)) stop("No exact Reactome name/ID match for: ", paste(x[is.na(idx)], collapse = "; "), ". Update the annotation mapping explicitly; labels are never shortened or guessed.",
      call. = FALSE)
    annotation$reactome_id[idx]
  }
  f7_attach_annotation <- function(x, annotation) {
    j <- match(x$reactome_id, annotation$reactome_id)
    for (nm in setdiff(names(annotation), c("reactome_id", "pathway_name"))) x[[nm]] <- annotation[[nm]][j]
    x$full_name <- annotation$pathway_name[j]
    overall <- is.na(x$reactome_id) | x$reactome_id == "overall"
    if (any(overall)) {
      x$full_name[overall] <- "Integrated immune age"
      x$class_id[overall] <- "integrated";
      x$class_name[overall] <- "Integrated clock"
    }
    x
  }

  # Expects ageacc and ageacc_z calibrated separately for each clock across all
  # eligible Chia exact-age samples, with all diagnosis groups pooled.
  f7_chia_associations <- function(long, meta) {
    f7_need(long, c("sample_id", "clock_id", "clock_name", "clock_type", "age", "ageacc", "ageacc_z", "eligible_age_analysis",
      "unit_id"), "Chia age acceleration")
    long <- as.data.frame(long, stringsAsFactors = FALSE);
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (anyDuplicated(long[c("sample_id", "clock_id")]) || anyDuplicated(rownames(meta)) || !all(long$sample_id %in%
      rownames(meta))) stop("Chia prediction/meta sample alignment failed.")
    audit <- list();
    warnings <- list()
    add <- function(type, contrast = NA_character_, clock_id = NA_character_, variable = NA_character_, detail =
      "") {
      audit[[length(audit) + 1L]] <<- data.frame(type = type, contrast = contrast, clock_id = clock_id, variable =
        variable, detail = detail, stringsAsFactors = FALSE)
    }
    resolve <- function(label, aliases, normalizer) {
      hits <- names(meta)[f7_key(names(meta)) %in% f7_key(aliases)]
      if (!length(hits)) {
        add("missing_covariate", variable = label, detail = "No matching metadata column");
        return(rep(NA_character_, nrow(meta)))
      }
      vals <- lapply(hits, function(nm) normalizer(meta[[nm]]))
      if (length(vals) > 1L && !all(vapply(vals[ - 1L], function(z) identical(z, vals[[1L]]), logical(1)))) stop("Conflicting metadata aliases for ",
        label, ": ", paste(hits, collapse = ", "))
      add("resolved_covariate", variable = label, detail = paste(hits, collapse = ";"));
      vals[[1L]]
    }
    sex_norm <- function(x) {
      z <- f7_key(f7_clean(x));
      z[z %in% c("m", "male")] <- "male";
      z[z %in% c("f", "female")] <- "female";
      z[is.na(f7_clean(x))] <- NA_character_;
      z
    }
    cat_norm <- function(x) {
      z <- f7_key(f7_clean(x));
      z[is.na(f7_clean(x))] <- NA_character_;
      z
    }
    diag <- resolve("diagnosis", c("Diagnosis", "diagnosis"), f7_diag)
    sex <- resolve("sex", c("Sex", "sex", "biological_sex", "subject.biologicalSex"), sex_norm)
    tube <- resolve("tube", c("Plasma_collection_tube_type", "plasma_collection_tube", "tube_type", "tube"),
      cat_norm)
    i <- match(long$sample_id, rownames(meta));
    d <- long
    d$diagnosis <- diag[i];
    d$sex <- sex[i];
    d$tube <- tube[i]
    d$unit_id <- as.character(d$unit_id)
    if (anyNA(d$unit_id) || any(!nzchar(trimws(d$unit_id)))) stop("Missing sample/participant unit IDs in Chia long predictions.")
    clocks <- unique(d[, c("clock_id", "clock_name", "clock_type")])
    if (anyDuplicated(clocks$clock_id)) stop("Clock identity conflict in Chia age acceleration.")
    if (sum(clocks$clock_type == "pathway") != 58L || sum(clocks$clock_type == "overall") != 1L) stop("Expected one integrated clock and 58 pathways in Chia predictions.")
    eligible <- f7_bool(d$eligible_age_analysis) & is.finite(d$age)
    dm <- d[eligible & !duplicated(d$sample_id), , drop = FALSE]
    contingency <- lapply(c("sex", "tube"), function(v) {
      tab <- as.data.frame(table(diagnosis = ifelse(is.na(dm$diagnosis), "Missing_or_unknown", dm$diagnosis),
        level = ifelse(is.na(dm[[v]]), "Missing", dm[[v]])), stringsAsFactors = FALSE)
      names(tab)[3L] <- "n";
      tab$covariate <- v;
      tab
    })
    contrasts <- c(ALS_vs_Healthy = "ALS", Neurological_control_vs_Healthy = "Neurological_control")
    rows <- list()
    for (contrast in names(contrasts)) for (j in seq_len(nrow(clocks))) {
      clock <- clocks$clock_id[j];
      cases <- unname(contrasts[contrast])
      z <- d[eligible & d$clock_id == clock & !is.na(d$diagnosis) & d$diagnosis %in% c("Healthy_control",
        cases), , drop = FALSE]
      z$y <- as.integer(z$diagnosis == cases)
      row <- data.frame(contrast = contrast, case_group = cases, reference_group = "Healthy_control", clock_id =
        clock, clock_name = clocks$clock_name[j], clock_type = clocks$clock_type[j], model_spec = "age_sex_tube",
        status = "not_run", n_eligible = nrow(z), n = 0L, n_case = 0L, n_control = 0L, n_units = 0L, n_case_units =
        0L, n_control_units = 0L, n_missing = nrow(z), beta_per_chia_sd = NA_real_, se = NA_real_, OR_per_chia_sd =
        NA_real_, CI95_low = NA_real_, CI95_high = NA_real_, p = NA_real_, p_bonf_global = NA_real_, se_method =
        NA_character_, covariates_used = NA_character_, formula = NA_character_, converged = NA, design_rank =
        NA_integer_, design_ncol = NA_integer_, fitted_probability_min = NA_real_, fitted_probability_max =
        NA_real_, warning = "", stringsAsFactors = FALSE)
      cc <- stats::complete.cases(z[, c("ageacc_z", "age", "sex", "tube"), drop = FALSE]) & is.finite(z$ageacc_z) &
      is.finite(z$age)
      zz <- z[cc, , drop = FALSE]
      row$n <- nrow(zz);
      row$n_case <- sum(zz$y == 1L);
      row$n_control <- sum(zz$y == 0L)
      row$n_missing <- nrow(z) - nrow(zz);
      row$n_units <- length(unique(zz$unit_id))
      row$n_case_units <- length(unique(zz$unit_id[zz$y == 1L]));
      row$n_control_units <- length(unique(zz$unit_id[zz$y == 0L]))
      used <- c("age", "sex", "tube")
      constants <- c("sex", "tube")[vapply(c("sex", "tube"), function(v) length(unique(zz[[v]])) <= 1L, logical(1))]
      for (v in constants) add("constant_covariate_omitted", contrast, clock, v, "No variation among complete-case rows")
      used <- setdiff(used, constants);
      row$covariates_used <- paste(used, collapse = ";")
      form <- stats::reformulate(c("ageacc_z", used), response = "y");
      row$formula <- paste(deparse(form), collapse = "")
      if (!row$n_case || !row$n_control) row$status <- "missing_case_or_control" else if (length(unique(zz$ageacc_z)) <
        2L) row$status <- "constant_ageacc" else {
        for (v in intersect(c("sex", "tube"), used)) zz[[v]] <- droplevels(factor(zz[[v]]))
        fw <- character()
        out <- tryCatch(withCallingHandlers({
          X <- stats::model.matrix(form, zz);
          row$design_rank <- qr(X)$rank;
          row$design_ncol <- ncol(X)
          if (nrow(X) <= ncol(X) || qr(X)$rank < ncol(X)) stop("rank_deficient_or_insufficient_rows")

          # There must be at least one tube level shared by cases and controls.
          # A perfect tube/diagnosis split cannot identify a tube-adjusted association.
          if (!length(intersect(as.character(zz$tube[zz$y == 1L]), as.character(zz$tube[zz$y == 0L])))) stop("no_shared_tube_level_between_cases_and_controls")
          fit <- stats::glm(form, data = zz, family = stats::binomial(), control = stats::glm.control(maxit =
            100L))
          row$converged <- isTRUE(fit$converged)
          if (!isTRUE(fit$converged)) stop("nonconverged_fit")
          if (any(!is.finite(stats::coef(fit)))) stop("nonfinite_coefficients")
          mu <- stats::fitted(fit);
          row$fitted_probability_min <- min(mu);
          row$fitted_probability_max <- max(mu)
          if (isTRUE(fit$boundary)) fw <- c(fw, "glm boundary flag; finite converged estimate retained")
          if (any(mu < 1e-8 | mu > 1 - 1e-8)) fw <- c(fw, "Extreme fitted probabilities; finite converged estimate retained")
          V <- stats::vcov(fit);
          se_method <- "model_based"
          if (anyDuplicated(zz$unit_id) > 0L) {
            G <- length(unique(zz$unit_id));
            N <- nrow(X);
            P <- ncol(X)
            if (G < 2L) stop("insufficient_independent_units")
            bread <- solve(crossprod(X, X * as.numeric(mu * (1 - mu))))
            scores <- rowsum(X * as.numeric(zz$y - mu), zz$unit_id, reorder = FALSE)
            V <- G / (G - 1) * (N - 1) / (N - P) * bread %*% crossprod(scores) %*% bread
            se_method <- "participant_cluster_sandwich_HC1"
          }
          b <- unname(stats::coef(fit)["ageacc_z"]);
          se <- sqrt(V["ageacc_z", "ageacc_z"])
          val <- c(b, se, exp(b), exp(b - stats::qnorm(.975) * se), exp(b + stats::qnorm(.975) * se))
          if (any(!is.finite(val)) || se <= 0 || any(val[3:5] <= 0)) stop("nonfinite_or_unestimable_inference")
          list(beta = b, se = se, OR = val[3L], lower = val[4L], upper = val[5L], p = 2 * stats::pnorm( -
            abs(b / se)), se_method = se_method)
        }, warning = function(w) {
          fw <<- c(fw, conditionMessage(w));
          invokeRestart("muffleWarning")
        }), error = function(e) e)
        if (inherits(out, "error")) {
          row$status <- "non_estimable";
          fw <- c(conditionMessage(out), fw)
        } else {
          row$status <- if (length(fw)) "ok_with_warning" else "ok"
          row$beta_per_chia_sd <- out$beta;
          row$se <- out$se;
          row$OR_per_chia_sd <- out$OR
          row$CI95_low <- out$lower;
          row$CI95_high <- out$upper;
          row$p <- out$p;
          row$se_method <- out$se_method
        }
        row$warning <- paste(unique(fw), collapse = " | ")
      }
      if (nzchar(row$warning)) warnings[[length(warnings) + 1L]] <- data.frame(contrast = contrast, clock_id =
        clock, message = row$warning, stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- row
    }
    associations <- f7_bind(rows)
    pathway <- associations$clock_type == "pathway";
    overall <- associations$clock_type == "overall"
    associations$p_bonf_global[pathway] <- pmin(1, associations$p[pathway] * 116L)
    associations$p_bonf_global[overall] <- pmin(1, associations$p[overall] * 2L)
    associations$bonferroni_family_size <- ifelse(pathway, 116L, 2L)
    associations$significant <- is.finite(associations$p_bonf_global) & associations$p_bonf_global < .05
    associations$estimable <- is.finite(associations$OR_per_chia_sd) & is.finite(associations$CI95_low) &
    is.finite(associations$CI95_high)
    list(associations = associations, distribution = d, audit = f7_bind(audit), warnings = f7_bind(warnings),
      contingency = f7_bind(contingency))
  }

  # Accuracy estimates are read from the completed external evaluation.
  # Point estimates are recalculated only to verify the saved records, never to
  # select/retrain a clock or replace the saved bootstrap confidence intervals.
  f7_prepare <- function(post, raw, annotation) {
    chia_group <- "Healthy_control"
    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
    f7_need(annotation, c("reactome_id", "pathway_name", "class_id", "class_name"), "Reactome annotation")
    if (nrow(annotation) != 58L || anyDuplicated(annotation$reactome_id) || anyDuplicated(annotation$pathway_name)) stop("Annotation must identify exactly 58 unique pathways.")
    run <- c(as.character(raw$run_id), as.character(post$source_model_run))
    pipeline <- c(as.character(raw$best$pipeline), as.character(post$pipeline$pipeline))
    if (length(run) != 2L || anyNA(run) || length(unique(run)) != 1L) stop("Raw and postprocess results are from different or unidentified runs.")
    if (length(pipeline) != 2L || anyNA(pipeline) || length(unique(pipeline)) != 1L) stop("Saved pipeline identity mismatch.")
    if (is.null(post$clock_map) || is.null(post$performance)) stop("Supply the completed external postprocessing result.")
    clock_map <- as.data.frame(post$clock_map, stringsAsFactors = FALSE)
    f7_need(clock_map, c("clock_id", "clock_name", "clock_type", "prediction_column"), "Postprocess clock map")
    if (nrow(clock_map) != 59L || sum(clock_map$clock_type == "pathway") != 58L || sum(clock_map$clock_type ==
      "overall") != 1L || anyDuplicated(clock_map$clock_id)) stop("Expected one integrated clock and 58 pathways.")
    clock_map$reactome_id <- "overall"
    ip <- clock_map$clock_type == "pathway"
    clock_map$reactome_id[ip] <- f7_match_pathway(clock_map$clock_name[ip], annotation)
    if (anyDuplicated(clock_map$reactome_id) || !setequal(clock_map$reactome_id, c("overall", annotation$reactome_id))) stop("Clock-to-Reactome mapping is incomplete or duplicated.")
    scopes <- c(imm_of_aging = "all_samples", chia_als = chia_group)
    cohorts <- c(imm_of_aging = "Gong", chia_als = "Chia")
    perf <- pred <- verification <- summaries <- list()
    raw_perf <- as.data.frame(post$performance, stringsAsFactors = FALSE)
    f7_need(raw_perf, c("dataset", "group", "clock_id", "metric", "estimate", "lower95", "upper95", "N"),
      "Saved performance")
    same_num <- function(a, b, tol = 1e-8) {
      if (length(a) != length(b)) return(FALSE)
      if (!identical(is.na(a), is.na(b))) return(FALSE)
      ii <- !is.na(a);
      if (!all(is.finite(a[ii]) == is.finite(b[ii]))) return(FALSE)
      if (!all(a[ii][!is.finite(a[ii])] == b[ii][!is.finite(b[ii])])) return(FALSE)
      kk <- ii & is.finite(a) & is.finite(b)
      all(abs(a[kk] - b[kk]) <= tol)
    }
    meta_values <- function(m, aliases) {
      h <- names(m)[f7_key(names(m)) %in% f7_key(aliases)]
      if (!length(h)) return(rep(NA_character_, nrow(m)))
      val <- lapply(h, function(nm) f7_clean(m[[nm]]))
      if (length(h) > 1L && !all(vapply(val[ - 1L], function(x) identical(f7_key(x), f7_key(val[[1L]])), logical(1)))) stop("Conflicting metadata aliases: ",
        paste(h, collapse = "; "))
      val[[1L]]
    }
    for (dataset in names(scopes)) {
      z <- post$external[[dataset]];
      rz <- raw$external[[dataset]]
      if (is.null(z$ageacc) || is.null(z$prediction) || is.null(rz$prediction) || is.null(z$meta)) stop("Missing saved predictions/metadata for ",
        dataset)
      p <- as.data.frame(z$prediction, check.names = FALSE)
      rp <- as.data.frame(rz$prediction, check.names = FALSE)
      f7_need(p, c("sample_id", "label", "eligible_age_analysis", clock_map$prediction_column), paste(dataset,
        "predictions"))
      f7_need(rp, c("sample_id", clock_map$prediction_column), paste(dataset, "raw predictions"))
      if (anyDuplicated(p$sample_id) || anyDuplicated(rp$sample_id) || !setequal(p$sample_id, rp$sample_id)) stop(dataset,
        " raw/post sample identity mismatch.")
      rp <- rp[match(p$sample_id, rp$sample_id), , drop = FALSE]
      for (nm in clock_map$prediction_column) if (!same_num(p[[nm]], rp[[nm]], 1e-12)) stop(dataset, " raw/post prediction mismatch for ",
        nm)
      long <- as.data.frame(z$ageacc, stringsAsFactors = FALSE)
      f7_need(long, c("sample_id", "clock_id", "clock_type", "age", "predicted_age", "eligible_age_analysis",
        "unit_id", "diagnosis"), paste(dataset, "age acceleration"))
      if (anyDuplicated(long[c("sample_id", "clock_id")]) || !setequal(unique(long$clock_id), clock_map$clock_id) ||
        nrow(long) != nrow(p) * 59L) stop(dataset, " incomplete or duplicated 59-clock prediction records.")
      if (anyNA(long$unit_id) || any(!nzchar(trimws(as.character(long$unit_id))))) stop(dataset, " missing participant/sample unit identifiers.")
      overall <- long[long$clock_type == "overall", , drop = FALSE]
      if (anyDuplicated(overall$sample_id) || !setequal(overall$sample_id, p$sample_id)) stop(dataset, " integrated sample alignment failed.")
      overall <- overall[match(p$sample_id, overall$sample_id), , drop = FALSE]
      group_keep <- rep(TRUE, nrow(p))
      if (scopes[[dataset]] != "all_samples") group_keep <- !is.na(overall$diagnosis) & overall$diagnosis ==
      scopes[[dataset]]
      keep <- f7_bool(p$eligible_age_analysis) & is.finite(p$label) & group_keep
      q <- raw_perf[raw_perf$dataset == dataset & raw_perf$group == scopes[[dataset]] & raw_perf$metric %in%
        c("MAE", "PCC", "RMSE"), , drop = FALSE]
      if (nrow(q) != 177L || anyDuplicated(q[c("clock_id", "metric")]) || !setequal(q$clock_id, clock_map$clock_id)) stop(dataset,
        " expected saved MAE/PCC/RMSE and bootstrap CIs for all59 clocks.")
      if (!"n_units" %in% names(q)) q$n_units <- NA_integer_
      for (j in seq_len(nrow(clock_map))) {
        lj <- long[long$clock_id == clock_map$clock_id[j], , drop = FALSE]
        if (!setequal(lj$sample_id, p$sample_id)) stop(dataset, " incomplete sample set for clock ", clock_map$clock_id[j])
        lj <- lj[match(p$sample_id, lj$sample_id), , drop = FALSE]
        pr <- p[[clock_map$prediction_column[j]]]
        if (!same_num(lj$predicted_age, pr, 1e-12) || !same_num(lj$age, p$label, 1e-12) || !identical(f7_bool(lj$eligible_age_analysis),
          f7_bool(p$eligible_age_analysis))) stop(dataset, " long/wide prediction or age mismatch for ", clock_map$clock_id[j])
        if (!identical(as.character(lj$unit_id), as.character(overall$unit_id))) stop(dataset, " unit IDs vary across clocks.")
        use <- keep & is.finite(pr)
        y <- p$label[use];
        est <- pr[use]
        if (length(y) < 2L) stop(dataset, " fewer than two finite observations for clock ", clock_map$clock_id[j])
        calculated <- c(MAE = mean(abs(est - y)), RMSE = sqrt(mean((est - y) ^ 2)), PCC = if (stats::sd(y) >
          0 && stats::sd(est) > 0) stats::cor(y, est) else NA_real_)
        ii <- which(q$clock_id == clock_map$clock_id[j])
        expected <- calculated[as.character(q$metric[ii])]
        if (!same_num(q$estimate[ii], as.numeric(expected), 1e-8) || any(q$N[ii] != length(y))) stop(dataset,
          " saved performance does not match selected predictions for ", clock_map$clock_id[j])
        nu <- length(unique(lj$unit_id[use]))
        if (any(!is.na(q$n_units[ii]) & q$n_units[ii] != nu)) stop(dataset, " saved bootstrap unit count disagrees with predictions.")
        q$n_units[ii] <- nu
        verification[[length(verification) + 1L]] <- data.frame(dataset = dataset, clock_id = clock_map$clock_id[j],
          N = length(y), n_units = nu, metric = q$metric[ii], absolute_difference = abs(q$estimate[ii] -
          as.numeric(expected)), stringsAsFactors = FALSE)
      }
      q$cohort <- unname(cohorts[dataset])
      q$cohort_label <- if (dataset == "imm_of_aging") "Gong cohort" else "Chia: healthy controls"
      q$reactome_id <- clock_map$reactome_id[match(q$clock_id, clock_map$clock_id)]
      q$clock_type <- clock_map$clock_type[match(q$clock_id, clock_map$clock_id)]
      perf[[dataset]] <- q
      v <- overall[keep & is.finite(overall$predicted_age), , drop = FALSE]
      v$cohort <- unname(cohorts[dataset]);
      v$cohort_label <- q$cohort_label[1L];
      v$dataset <- dataset;
      v$group <- unname(scopes[dataset])
      pred[[dataset]] <- v

      # Cohort characteristics describe the exact-age population, with available
      # sample counts recorded separately; they are not model complete-case counts.
      m <- as.data.frame(z$meta, stringsAsFactors = FALSE)
      if (anyDuplicated(rownames(m)) || !all(p$sample_id %in% rownames(m))) stop(dataset, " metadata/sample mismatch.")
      m <- m[match(p$sample_id, rownames(m)), , drop = FALSE]
      sex <- f7_key(meta_values(m, c("Sex", "sex", "biological_sex", "subject.biologicalSex")))
      tube <- meta_values(m, c("Plasma_collection_tube_type", "plasma_collection_tube", "tube_type", "tube"))
      groups <- if (dataset == "chia_als") c("Healthy_control", "ALS", "Neurological_control") else "all_samples"
      for (g in groups) {
        ix <- if (g == "all_samples") rep(TRUE, nrow(p)) else !is.na(overall$diagnosis) & overall$diagnosis ==
        g
        ex <- ix & f7_bool(p$eligible_age_analysis) & is.finite(p$label)
        ag <- p$label[ex]
        tt <- sort(table(ifelse(is.na(tube[ex]), "Missing", tube[ex])), decreasing = TRUE)
        sex_known <- sex[ex] %in% c("m", "male", "f", "female")
        summaries[[length(summaries) + 1L]] <- data.frame(cohort = unname(cohorts[dataset]), group = g, n_samples_available =
          sum(ix), n_units_available = length(unique(overall$unit_id[ix])), n_exact_age = sum(ex), n_exact_age_units =
          length(unique(overall$unit_id[ex])), n_in_prediction_analysis = sum(ex & group_keep & is.finite(overall$predicted_age)),
          age_mean = if (length(ag)) mean(ag) else NA_real_, age_sd = if (length(ag) > 1L) stats::sd(ag) else NA_real_,
          age_min = if (length(ag)) min(ag) else NA_real_, age_max = if (length(ag)) max(ag) else NA_real_,
          female_n = sum(sex[ex] %in% c("f", "female")), male_n = sum(sex[ex] %in% c("m", "male")), sex_other_or_missing_n =
          sum(!sex_known), collection_tube_counts = paste(paste(names(tt), as.integer(tt), sep = ": "), collapse =
          "; "), stringsAsFactors = FALSE)
      }
    }
    performance <- f7_attach_annotation(f7_bind(perf), annotation)

    # Wide accuracy table: one record per cohort and clock, all CI values unchanged.
    wide <- list()
    for (cohort in c("Gong", "Chia")) {
      u <- performance[performance$cohort == cohort & performance$metric == "PCC", , drop = FALSE]
      u <- u[match(clock_map$clock_id, u$clock_id), , drop = FALSE]
      out <- u[, c("cohort", "cohort_label", "group", "clock_id", "reactome_id", "full_name", "clock_type",
        "class_id", "class_name", "N", "n_units"), drop = FALSE]
      for (metric in c("MAE", "RMSE", "PCC")) {
        m <- performance[performance$cohort == cohort & performance$metric == metric, , drop = FALSE]
        m <- m[match(out$clock_id, m$clock_id), , drop = FALSE]
        out[[metric]] <- m$estimate;
        out[[paste0(metric, "_lower95")]] <- m$lower95;
        out[[paste0(metric, "_upper95")]] <- m$upper95
      }
      wide[[cohort]] <- out
    }
    accuracy <- f7_bind(wide)
    cohort_summary <- f7_bind(summaries)
    audit <- list(prediction_verification = f7_bind(verification))
    settings <- list(run_id = raw$run_id, pipeline = raw$best$pipeline, chia_prediction_group = chia_group,
      accuracy = "Saved prediction point estimates and 95% percentile bootstrap intervals; MAE and RMSE in years; PCC is Pearson correlation with chronological age.",
      disease_contrasts = c("ALS vs healthy controls", "Other neurological disease controls vs healthy controls"),
      age_acceleration = "Each clock uses OLS predicted_age ~ chronological age across all eligible Chia exact-age samples, pooling healthy controls, ALS and other diagnoses. Age acceleration is the residual in years; residuals are centered and divided by their Chia-cohort SD. Calibration does not include diagnosis or fit groups separately. Adjusted OR per 1 Chia-cohort SD.",
      correction = "Integrated: Bonferroni 2. Pathways: Bonferroni 116 across the two planned contrasts.",
      counts = "Model counts use complete cases for age acceleration, age, sex and collection tube. Repeated samples use participant-cluster sandwich HC1 standard errors when participant IDs are available.",
      warning_policy = "Finite, converged, estimable warning fits retained. Failed/non-estimable models remain as NA rows; no pathway is silently removed.")
    clean_accuracy <- accuracy[, c("cohort", "group", "reactome_id", "full_name", "clock_type", "class_name",
      "N", "n_units", "MAE", "MAE_lower95", "MAE_upper95", "RMSE", "RMSE_lower95", "RMSE_upper95", "PCC",
      "PCC_lower95", "PCC_upper95"), drop = FALSE]
    tables <- list(Cohort_characteristics = cohort_summary, Prediction_accuracy = clean_accuracy)
    ans <- list(performance = performance, prediction = f7_bind(pred), cohort_summary = cohort_summary, clock_map =
      f7_attach_annotation(clock_map, annotation), annotation = annotation, accuracy = accuracy, tables =
      tables, audit = audit, settings = settings)
    ans
  }

  # Frozen affine pathway contributions and UKB comparison
  f7_affine_coefficients <- function(predict_stage, Xg, allx, tolerance = 1e-6) {
    mu <- colMeans(Xg)

    # Evaluate a basis through the frozen prediction API, preserving all saved
    # preprocessing. This extracts raw-pathway-age coefficients without fitting.
    m <- ncol(Xg);
    step <- pmax(apply(Xg, 2L, stats::sd), 1)
    centre <- matrix(mu, nrow = 1L, dimnames = list(NULL, colnames(Xg)))
    plus <- minus <- matrix(rep(mu, each = m), nrow = m, dimnames = list(NULL, colnames(Xg)))
    diag(plus) <- diag(plus) + step;
    diag(minus) <- diag(minus) - step
    q <- predict_stage(rbind(centre, plus, minus))
    beta <- (q[1L + seq_len(m)] - q[1L + m + seq_len(m)]) / (2 * step)
    intercept <- q[1L] - sum(beta * mu)
    centred <- sweep(allx, 2L, mu, "-")

    # Numerical affine checks at observed and mixed inputs; no age labels enter.
    probes <- rbind(centre, plus, minus, sweep(centred[seq_len(min(30L, nrow(centred))), , drop = FALSE] *
      .37, 2L, mu, "+"), sweep(centred[seq_len(min(30L, nrow(centred))), , drop = FALSE] * - 1.2, 2L, mu,
      "+"))
    probe_error <- max(abs(predict_stage(probes) - (intercept + drop(probes %*% beta))))
    if (!is.finite(probe_error) || probe_error > tolerance) stop("The stage-2 model fails the affine identity check; approximate contributions were not generated.")
    list(beta = beta, intercept = intercept, maximum_probe_error_years = probe_error)
  }
  f7_add_contribution <- function(prepared, post, raw, bundle, ukb, tolerance = 1e-6) {
    if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) stop("Invalid attribution verification tolerance.")
    if (!identical(bundle$schema_version, "full_training_two_external_common_v2") || !identical(bundle$run_id,
      raw$run_id) || !identical(bundle$best, raw$best) || !identical(post$source_model_run, raw$run_id) ||
      !identical(post$pipeline, raw$best) || !identical(prepared$settings$run_id, raw$run_id)) stop("Prepared data, frozen model and upstream results must identify the same run.")
    ann <- prepared$annotation
    mp <- as.data.frame(bundle$module_name_map, stringsAsFactors = FALSE)
    f7_need(mp, c("original", "internal"), "Frozen pathway mapping")
    if (nrow(mp) != 58L || anyNA(mp) || anyDuplicated(mp$original) || anyDuplicated(mp$internal)) stop("Frozen model must map 58 unique pathways.")
    mp$reactome_id <- f7_match_pathway(mp$original, ann)
    if (anyDuplicated(mp$reactome_id) || !setequal(mp$reactome_id, ann$reactome_id)) stop("Frozen model pathway IDs disagree with annotation.")
    f7_need(ukb, c("reactome_id", "n", "rank", "coefficient", "contribution_years", "contribution_pct"), "Original UKB contribution reference")
    if (nrow(ukb) != 58L || anyDuplicated(ukb$reactome_id) || !setequal(ukb$reactome_id, mp$reactome_id)) stop("UKB contribution reference must contain all 58 pathways exactly once.")
    for (v in c("n", "rank", "coefficient", "contribution_years", "contribution_pct")) if (!is.numeric(ukb[[v]]) ||
      any(!is.finite(ukb[[v]]))) stop("Invalid numeric UKB column: ", v)
    if (any(ukb$n != 10602L) || any(ukb$contribution_years < 0) || any(ukb$contribution_pct < 0) || sum(ukb$contribution_years) <=
      0 || abs(sum(ukb$contribution_pct) - 100) > 1e-6 || max(abs(ukb$contribution_pct - 100 * ukb$contribution_years /
      sum(ukb$contribution_years))) > 1e-6 || !setequal(ukb$rank, seq_len(58L)) || anyDuplicated(ukb$rank) ||
      any(diff(ukb$contribution_years[order(ukb$rank)]) > 1e-10)) stop("UKB reference values/ranks fail consistency checks.")
    mat <- pred <- list()
    for (ds in c("imm_of_aging", "chia_als")) {
      p <- as.data.frame(post$external[[ds]]$prediction, check.names = FALSE)
      rp <- as.data.frame(raw$external[[ds]]$prediction, check.names = FALSE)
      f7_need(p, c("sample_id", "immune_age", "label", "eligible_age_analysis", mp$original), paste(ds, "post predictions"))
      f7_need(rp, c("sample_id", "immune_age", mp$original), paste(ds, "raw predictions"))
      if (anyNA(p$sample_id) || anyDuplicated(p$sample_id) || anyDuplicated(rp$sample_id) || !setequal(p$sample_id,
        rp$sample_id)) stop("Raw/post participant IDs disagree: ", ds)
      rp <- rp[match(p$sample_id, rp$sample_id), , drop = FALSE]
      for (nm in c(mp$original, "immune_age")) {
        if (!is.numeric(p[[nm]]) || any(!is.finite(p[[nm]])) || !is.numeric(rp[[nm]]) || any(!is.finite(rp[[nm]])) ||
          max(abs(p[[nm]] - rp[[nm]])) > tolerance) stop("Saved raw/post predictions disagree: ", ds, " / ",
          nm)
      }
      X <- as.matrix(p[, mp$original, drop = FALSE]);
      colnames(X) <- mp$internal
      mat[[ds]] <- X;
      pred[[ds]] <- p
    }
    gong <- pred$imm_of_aging
    keep <- f7_bool(gong$eligible_age_analysis) & is.finite(gong$label)
    if (sum(keep) < 3L) stop("Insufficient exact-age Gong participants.")
    gp <- prepared$prediction[prepared$prediction$cohort %in% c("Gong", "Gong cohort"), , drop = FALSE]
    if (anyDuplicated(gp$sample_id) || !setequal(gp$sample_id, gong$sample_id[keep])) stop("Contribution participants must match the Gong accuracy population.")
    ix <- match(gong$sample_id[keep], gp$sample_id)
    if (max(abs(gp$predicted_age[ix] - gong$immune_age[keep])) > tolerance || max(abs(gp$age[ix] - gong$label[keep])) >
      tolerance) stop("Accuracy and contribution prediction/age pairs differ.")
    Xg <- mat$imm_of_aging[keep, , drop = FALSE]
    mu <- colMeans(Xg)
    stage <- bundle$models$stage2
    log_rows <- list()
    predict_stage <- if (identical(stage$kind, "fusion")) {
      w <- as.data.frame(stage$weights)
      f7_need(w, c("module_internal", "weight"), "Frozen fusion weights")
      if (nrow(w) != 58L || anyDuplicated(w$module_internal) || !setequal(w$module_internal, mp$internal) ||
        any(!is.finite(w$weight))) stop("Frozen fusion weights do not cover all pathways.")
      w <- w$weight[match(mp$internal, w$module_internal)]
      function(x) drop(x[, mp$internal, drop = FALSE] %*% w)
    } else if (identical(stage$kind, "learner")) {
      l <- stage$learner
      pkgs <- unique(c("data.table", "mlr3", as.character(l$packages)))
      absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(absent)) stop("Missing dependencies required by the frozen model: ", paste(absent, collapse =
        ", "), call. = FALSE)
      function(x) {
        if (is.null(l$model)) stop("Frozen learner has no fitted model.")
        q <- withCallingHandlers(l$predict_newdata(data.table::as.data.table(x)), warning = function(w) {
          log_rows[[length(log_rows) + 1L]] <<- data.frame(message = conditionMessage(w))
        })
        errors <- tryCatch(l$errors, error = function(e) NULL)
        if ((is.numeric(errors) && any(errors > 0)) || (!is.null(errors) && !is.numeric(errors) && length(errors) >
          0L)) stop("Frozen learner returned an error.")
        lg <- tryCatch(l$log, error = function(e) NULL)
        if (!is.null(lg) && NROW(lg)) {
          if ("class" %in% names(lg) && any(tolower(as.character(lg$class)) == "error")) stop("Frozen learner reported an error.")
          msg <- if ("msg" %in% names(lg)) lg$msg else if ("condition" %in% names(lg)) vapply(lg$condition,
            function(z) if (inherits(z, "condition")) conditionMessage(z) else paste(as.character(z), collapse =
            " "), character(1)) else unlist(lg, use.names = FALSE)
          if (any(grepl("using fallback learner", paste(msg, collapse = " "), ignore.case = TRUE))) stop("Frozen learner used a fallback prediction.")
        }
        j <- match(seq_len(nrow(x)), q$row_ids)
        if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Frozen prediction row IDs are missing or duplicated.")
        z <- q$response[j]
        if (length(z) != nrow(x) || any(!is.finite(z))) stop("Nonfinite frozen predictions.")
        as.numeric(z)
      }
    } else stop("Unsupported frozen stage-2 model kind.")
    af <- f7_affine_coefficients(predict_stage, Xg, do.call(rbind, mat), tolerance)
    beta <- af$beta;
    intercept <- af$intercept;
    probe_error <- af$maximum_probe_error_years
    check_rows <- lapply(names(mat), function(ds) {
      calc <- intercept + drop(mat[[ds]] %*% beta)
      frozen <- predict_stage(mat[[ds]])
      e <- max(abs(calc - pred[[ds]]$immune_age));
      ef <- max(abs(frozen - pred[[ds]]$immune_age))
      if (!is.finite(e) || !is.finite(ef) || e > tolerance || ef > tolerance) stop("The frozen model cannot reproduce the saved integrated predictions: ",
        ds, ". Verify the frozen model bundle.")
      data.frame(dataset = ds, n = nrow(mat[[ds]]), maximum_affine_error_years = e, maximum_frozen_prediction_error_years =
        ef, maximum_probe_error_years = probe_error)
    })
    cohorts <- c(imm_of_aging = "Gong", chia_als = "Chia")
    scopes <- c(imm_of_aging = "all_samples", chia_als = prepared$settings$chia_prediction_group)
    if (length(scopes) != 2L || !identical(scopes[["chia_als"]], "Healthy_control")) stop("Missing/invalid Chia contribution population.")
    all_contributions <- list();
    sum_errors <- numeric();
    overlaps <- integer()
    for (ds in names(cohorts)) {
      pp <- pred[[ds]]
      use <- f7_bool(pp$eligible_age_analysis) & is.finite(pp$label)
      if (ds == "chia_als") {
        long <- as.data.frame(post$external[[ds]]$ageacc, stringsAsFactors = FALSE)
        f7_need(long, c("sample_id", "clock_type", "diagnosis"), "Chia diagnosis for contribution population")
        overall <- long[long$clock_type == "overall", , drop = FALSE]
        if (anyDuplicated(overall$sample_id) || !setequal(overall$sample_id, pp$sample_id)) stop("Chia diagnosis/sample IDs do not match predictions.")
        diagnosis <- overall$diagnosis[match(pp$sample_id, overall$sample_id)]
        use <- use & !is.na(diagnosis) & diagnosis == scopes[[ds]]
      }
      plotted <- prepared$prediction[prepared$prediction$cohort == cohorts[[ds]], , drop = FALSE]
      if (sum(use) < 3L || anyDuplicated(plotted$sample_id) || !setequal(plotted$sample_id, pp$sample_id[use])) stop("Contribution population differs from the accuracy population: ",
        ds)
      ii <- match(pp$sample_id[use], plotted$sample_id)
      if (max(abs(plotted$predicted_age[ii] - pp$immune_age[use])) > tolerance || max(abs(plotted$age[ii] -
        pp$label[use])) > tolerance) stop("Accuracy prediction pairs disagree: ", ds)
      X <- mat[[ds]][use, , drop = FALSE];
      mu_cohort <- colMeans(X)
      phi <- sweep(sweep(X, 2L, mu_cohort, "-"), 2L, beta, "*")
      sum_error <- max(abs(rowSums(phi) - (pp$immune_age[use] - mean(pp$immune_age[use]))))
      if (sum_error > 2 * tolerance) stop("Participant-level centered additive identity failed: ", ds)
      imp <- colMeans(abs(phi))
      if (any(!is.finite(imp)) || sum(imp) <= 0) stop("Invalid or all-zero contributions: ", ds)
      d <- data.frame(cohort = cohorts[[ds]], group = scopes[[ds]], reactome_id = mp$reactome_id, N = nrow(X),
        coefficient = beta, background_age = mu_cohort, contribution_years = imp, contribution_pct = 100 *
        imp / sum(imp), stringsAsFactors = FALSE)
      d <- d[order( - d$contribution_years, d$reactome_id), , drop = FALSE]
      d$contribution_rank <- seq_len(nrow(d));
      rownames(d) <- NULL
      j <- match(d$reactome_id, ukb$reactome_id)
      d$UKB_N <- ukb$n[j];
      d$UKB_coefficient <- ukb$coefficient[j]
      d$UKB_contribution_years <- ukb$contribution_years[j];
      d$UKB_contribution_pct <- ukb$contribution_pct[j]
      d$UKB_contribution_rank <- ukb$rank[j]
      all_contributions[[ds]] <- d;
      sum_errors[ds] <- sum_error
      overlaps[ds] <- sum(d$UKB_contribution_rank[seq_len(10L)] <= 10L)
    }
    d <- do.call(rbind, all_contributions);
    rownames(d) <- NULL
    prepared$contributions <- d
    clean <- f7_attach_annotation(d, ann)
    prepared$tables$Pathway_contribution <- clean[, c("cohort", "group", "reactome_id", "full_name", "class_name",
      "N", "coefficient", "background_age", "contribution_years", "contribution_pct", "contribution_rank",
      "UKB_N", "UKB_coefficient", "UKB_contribution_years", "UKB_contribution_pct", "UKB_contribution_rank")]
    prepared$audit$contribution_verification <- do.call(rbind, check_rows)
    prepared$audit$contribution_verification$maximum_centered_sum_error_years <- unname(sum_errors[prepared$audit$contribution_verification$dataset])
    prepared$audit$contribution_warnings <- if (length(log_rows)) unique(do.call(rbind, log_rows)) else data.frame(message =
      character())
    prepared$settings$contribution_version <- "two_cohorts_centered_additive_v1"
    prepared$settings$contribution <- "Gong and Chia participants included in age-prediction accuracy; frozen common-protein integrator; mean(abs(beta_j * (pathway_age_ij - cohort_mean_j))). Each cohort is centered and normalized separately across all 58 pathways. UKB values/ranks reused unchanged from original held-out S04."
    prepared$settings$ranking <- "Rank all 58 pathways by decreasing within-cohort mean absolute centered additive contribution; exact ties use Reactome stable ID."
    prepared$settings$contribution_model_comparison <- "Both external cohorts use the completed common-protein external model; UKB uses the original main model. Comparison reflects both model and cohort differences."
    prepared$settings$contribution_intercept <- intercept
    prepared$settings$contribution_centered_sum_error_years <- sum_errors
    prepared$settings$contribution_top10_overlap <- overlaps
    prepared
  }

  # Chia-internal age correction and standardization
  # Chia-internal age-bias correction of predictions from the frozen UKB clocks.
  # This postprocessing function never fits or alters an age-prediction clock.
  # Every eligible exact-age Chia sample enters its clock's calibration, regardless
  # of diagnosis or availability of disease-model covariates. The same coefficients
  # and residual scale are then used for every disease group.
  # Run before choosing either a cached or freshly prepared External validation result.
  # The upstream wide predictions are independently checked against the frozen
  # raw predictions by f7_add_contribution() on both execution paths.
  f7_validate_chia_frozen_records <- function(long, wide, clock_map) {
    need <- function(x, columns, label) {
      if (!is.data.frame(x) || anyDuplicated(names(x)) || !all(columns %in% names(x))) stop("Chia frozen-record validation: invalid or incomplete ",
        label, ".", call. = FALSE)
    }
    need(clock_map, c("clock_id", "clock_name", "clock_type", "prediction_column"), "clock map")
    need(long, c("sample_id", "clock_id", "clock_name", "clock_type", "predicted_age", "age", "eligible_age_analysis"),
      "long records")
    need(wide, c("sample_id", "label", "eligible_age_analysis", as.character(clock_map$prediction_column)),
      "wide predictions")
    valid_ids <- function(x) !anyNA(x) && all(nzchar(trimws(as.character(x))))
    if (nrow(clock_map) != 59L || !valid_ids(clock_map$clock_id) || !valid_ids(clock_map$clock_name) || !valid_ids(clock_map$clock_type) ||
      !valid_ids(clock_map$prediction_column) || anyDuplicated(clock_map$clock_id) || anyDuplicated(clock_map$prediction_column) ||
      sum(clock_map$clock_type == "overall") != 1L || sum(clock_map$clock_type == "pathway") != 58L) stop("Chia frozen-record validation requires 59 uniquely mapped clocks.",
      call. = FALSE)
    if (!nrow(wide) || !valid_ids(wide$sample_id) || anyDuplicated(wide$sample_id) || !valid_ids(long$sample_id) ||
      !valid_ids(long$clock_id) || anyDuplicated(long[, c("sample_id", "clock_id"), drop = FALSE]) || nrow(long) !=
      59L * nrow(wide) || !setequal(long$clock_id, clock_map$clock_id)) stop("Chia frozen-record validation: incomplete or duplicated sample/clock grid.",
      call. = FALSE)
    same_num <- function(a, b) {
      if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b) || !identical(is.na(a), is.na(b))) return(FALSE)
      observed <- !is.na(a)
      if (!identical(is.finite(a[observed]), is.finite(b[observed]))) return(FALSE)
      infinite <- observed & !is.finite(a)
      finite <- observed & is.finite(a)
      all(a[infinite] == b[infinite]) && all(abs(a[finite] - b[finite]) <= 1e-12)
    }
    as_eligible <- function(x) {
      value <- tolower(trimws(as.character(x)))
      if (any(!is.na(value) & !value %in% c("true", "false", "1", "0"))) stop("Chia frozen-record validation: invalid eligibility flag.",
        call. = FALSE)
      !is.na(value) & value %in% c("true", "1")
    }
    wide_eligible <- as_eligible(wide$eligible_age_analysis)
    for (j in seq_len(nrow(clock_map))) {
      z <- long[as.character(long$clock_id) == as.character(clock_map$clock_id[j]), , drop = FALSE]
      if (!setequal(z$sample_id, wide$sample_id)) stop("Chia frozen-record validation: incomplete sample set for ",
        clock_map$clock_id[j], ".", call. = FALSE)
      z <- z[match(wide$sample_id, z$sample_id), , drop = FALSE]
      if (anyNA(z$clock_name) || anyNA(z$clock_type) || !all(as.character(z$clock_name) == as.character(clock_map$clock_name[j])) ||
        !all(as.character(z$clock_type) == as.character(clock_map$clock_type[j])) || !same_num(z$predicted_age,
        wide[[as.character(clock_map$prediction_column[j])]]) || !same_num(z$age, wide$label) || !identical(as_eligible(z$eligible_age_analysis),
        wide_eligible)) stop("Chia long/wide frozen prediction, age, eligibility or clock identity mismatch for ",
        clock_map$clock_id[j], ".", call. = FALSE)
    }
    invisible(TRUE)
  }
  f7_chia_internal_calibration <- function(long) {
    required <- c("sample_id", "clock_id", "clock_name", "clock_type", "age", "predicted_age", "eligible_age_analysis",
      "unit_id", "ageacc", "ageacc_z")
    if (!is.data.frame(long)) stop("Chia calibration requires a data frame.", call. = FALSE)
    missing <- setdiff(required, names(long))
    if (length(missing)) stop("Chia calibration is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    if (anyDuplicated(names(long))) stop("Chia calibration has duplicated column names.", call. = FALSE)
    if (!is.numeric(long$age) || !is.numeric(long$predicted_age)) stop("Chia calibration requires numeric age and predicted_age columns.",
      call. = FALSE)
    if (!nrow(long)) stop("Chia calibration requires all 59 clocks and their samples.", call. = FALSE)

    # Keep identifiers and eligibility as local working vectors. Input columns,
    # classes, attributes and row order are preserved in the returned long object.
    sample_id <- as.character(long$sample_id)
    clock_id <- as.character(long$clock_id)
    unit_id <- as.character(long$unit_id)
    invalid_id <- function(x) anyNA(x) || any(!nzchar(trimws(x)))
    if (invalid_id(sample_id) || invalid_id(clock_id) || invalid_id(unit_id)) stop("Chia calibration has missing or blank sample, clock or unit IDs.",
      call. = FALSE)
    if (anyDuplicated(long[, c("sample_id", "clock_id"), drop = FALSE])) stop("Chia calibration has duplicate sample/clock records.",
      call. = FALSE)
    clocks <- unique(long[, c("clock_id", "clock_name", "clock_type"), drop = FALSE])
    if (invalid_id(as.character(clocks$clock_name)) || invalid_id(as.character(clocks$clock_type)) || anyDuplicated(clocks$clock_id)) stop("Chia calibration has conflicting or missing clock labels/types.",
      call. = FALSE)
    if (nrow(clocks) != 59L || sum(clocks$clock_type == "overall") != 1L || sum(clocks$clock_type == "pathway") !=
      58L) stop("Chia calibration requires one integrated clock and 58 pathway clocks.", call. = FALSE)
    rownames(clocks) <- NULL
    eligibility_text <- tolower(trimws(as.character(long$eligible_age_analysis)))
    if (any(!is.na(eligibility_text) & !eligibility_text %in% c("true", "false", "1", "0"))) stop("Chia calibration has an invalid eligible_age_analysis value.",
      call. = FALSE)
    eligible <- !is.na(eligibility_text) & eligibility_text %in% c("true", "1")
    sample_set <- unique(sample_id)
    reference <- match(sample_set, sample_id)
    clock_rows <- lapply(as.character(clocks$clock_id), function(id) which(clock_id == id))
    for (j in seq_along(clock_rows)) {
      ii <- clock_rows[[j]]
      if (length(ii) != length(sample_set) || !setequal(sample_id[ii], sample_set)) stop("Chia calibration has an incomplete sample set for ",
        clocks$clock_id[j], ".", call. = FALSE)
      aligned <- ii[match(sample_set, sample_id[ii])]
      if (!identical(unname(long$age[aligned]), unname(long$age[reference])) || !identical(unname(eligible[aligned]),
        unname(eligible[reference])) || !identical(unname(unit_id[aligned]), unname(unit_id[reference]))) stop("Chia calibration has inconsistent age, eligibility or unit IDs across clocks.",
        call. = FALSE)
    }

    # Deep-copy also protects callers using data-frame subclasses with references.
    result_long <- unserialize(serialize(long, NULL))
    result_long$ageacc <- rep(NA_real_, nrow(long))
    result_long$ageacc_z <- rep(NA_real_, nrow(long))
    sd_tolerance <- sqrt(.Machine$double.eps)
    moment_tolerance <- 1e-8
    calibration_rows <- vector("list", 59L)
    audit_rows <- vector("list", 59L)
    for (j in seq_len(nrow(clocks))) {
      ii <- clock_rows[[j]]
      keep <- ii[eligible[ii] & is.finite(long$age[ii]) & is.finite(long$predicted_age[ii])]
      n <- length(keep)
      calibration <- data.frame(clocks[j, , drop = FALSE], n_calibration = n, n_calibration_units = length(unique(unit_id[keep])),
        intercept = NA_real_, slope = NA_real_, residual_mean_years = NA_real_, residual_sd_years = NA_real_,
        stringsAsFactors = FALSE, check.names = FALSE)
      audit <- data.frame(clocks[j, , drop = FALSE], status = "not_run", n_samples = length(ii), n_eligible_age =
        sum(eligible[ii]), n_calibration = n, n_calibration_units = length(unique(unit_id[keep])), n_ineligible_age =
        sum(!eligible[ii]), n_eligible_nonfinite_age = sum(eligible[ii] & !is.finite(long$age[ii])), n_eligible_nonfinite_prediction =
        sum(eligible[ii] & is.finite(long$age[ii]) & !is.finite(long$predicted_age[ii])), n_ageacc = 0L, n_ageacc_z =
        0L, fit_rank = NA_integer_, ageacc_available = FALSE, z_available = FALSE, residual_sd_near_zero =
        NA, mean_ageacc_years = NA_real_, sd_ageacc_years = NA_real_, mean_ageacc_z = NA_real_, sd_ageacc_z =
        NA_real_, corr_ageacc_age = NA_real_, corr_ageacc_z_age = NA_real_, z_mean_near_zero = NA, z_sd_near_one =
        NA, age_correlation_near_zero = NA, warning = "", error = "", stringsAsFactors = FALSE, check.names =
        FALSE)
      local_warnings <- character()
      if (n < 2L) {
        audit$status <- "insufficient_calibration_pairs"
        audit$error <- "At least two finite eligible age/prediction pairs are needed."
      } else {
        fit_result <- tryCatch(withCallingHandlers({
          fit_data <- data.frame(predicted_age = long$predicted_age[keep], age = long$age[keep])
          fit <- stats::lm(predicted_age ~ age, data = fit_data, na.action = stats::na.fail)
          audit$fit_rank <- fit$rank
          if (fit$rank != 2L) stop("Age-bias regression is rank deficient.")
          coefficients <- unname(stats::coef(fit))
          if (length(coefficients) != 2L || any(!is.finite(coefficients))) stop("Age-bias regression coefficients are nonfinite.")
          aa <- fit_data$predicted_age - (coefficients[1L] + coefficients[2L] * fit_data$age)
          if (any(!is.finite(aa))) stop("Age-bias regression produced nonfinite residuals.")
          residual_mean <- mean(aa)
          residual_sd <- stats::sd(aa)
          if (!is.finite(residual_mean) || !is.finite(residual_sd)) stop("Age-bias residual moments are nonfinite.")
          list(coefficients = coefficients, ageacc = aa, residual_mean = residual_mean, residual_sd = residual_sd)
        }, warning = function(w) {
          local_warnings <<- c(local_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
        if (inherits(fit_result, "error")) {
          audit$status <- "calibration_fit_failed"
          audit$error <- conditionMessage(fit_result)
        } else {
          calibration$intercept <- fit_result$coefficients[1L]
          calibration$slope <- fit_result$coefficients[2L]
          calibration$residual_mean_years <- fit_result$residual_mean
          calibration$residual_sd_years <- fit_result$residual_sd
          result_long$ageacc[keep] <- fit_result$ageacc
          audit$n_ageacc <- n
          audit$ageacc_available <- TRUE
          audit$mean_ageacc_years <- fit_result$residual_mean
          audit$sd_ageacc_years <- fit_result$residual_sd
          audit$residual_sd_near_zero <- fit_result$residual_sd <= sd_tolerance
          if (audit$residual_sd_near_zero) {

            # A valid linear fit can have zero residual variance. Its AA values
            # remain available, but division by that scale is not meaningful.
            audit$status <- "residual_sd_near_zero"
            local_warnings <- c(local_warnings, "Residual SD is at or below sqrt(machine epsilon) years; ageacc_z is unavailable.")
          } else {
            az <- (fit_result$ageacc - fit_result$residual_mean) / fit_result$residual_sd
            if (any(!is.finite(az))) {
              audit$status <- "standardization_failed"
              audit$error <- "Chia-internal standardization produced nonfinite values."
            } else {
              result_long$ageacc_z[keep] <- az
              audit$n_ageacc_z <- n
              audit$z_available <- TRUE
              audit$mean_ageacc_z <- mean(az)
              audit$sd_ageacc_z <- stats::sd(az)
              audit$corr_ageacc_age <- stats::cor(fit_result$ageacc, long$age[keep])
              audit$corr_ageacc_z_age <- stats::cor(az, long$age[keep])
              audit$z_mean_near_zero <- abs(audit$mean_ageacc_z) <= moment_tolerance
              audit$z_sd_near_one <- abs(audit$sd_ageacc_z - 1) <= moment_tolerance
              audit$age_correlation_near_zero <- abs(audit$corr_ageacc_age) <= moment_tolerance
              if (!isTRUE(audit$z_mean_near_zero) || !isTRUE(audit$z_sd_near_one) || !isTRUE(audit$age_correlation_near_zero)) local_warnings <-
              c(local_warnings, "Calibration moment audit exceeds the numerical tolerance.")
              audit$status <- if (length(local_warnings)) "ok_with_warning" else "ok"
            }
          }
        }
      }
      audit$warning <- paste(unique(local_warnings), collapse = " | ")
      calibration_rows[[j]] <- calibration
      audit_rows[[j]] <- audit
    }
    calibration <- do.call(rbind, calibration_rows)
    audit <- do.call(rbind, audit_rows)
    rownames(calibration) <- rownames(audit) <- NULL
    preserved <- setdiff(names(long), c("ageacc", "ageacc_z"))
    if (!identical(result_long[, preserved, drop = FALSE], long[, preserved, drop = FALSE])) stop("Chia calibration changed protected input columns.",
      call. = FALSE)
    settings <- list(cohort = "Chia", calibration_population = "All Chia eligible exact-age samples with finite age and frozen predicted age, jointly across diagnoses",
      calibration_timing = "Before disease-group restriction or disease-model complete-case selection", age_clock =
      "Frozen UKB age-clock predictions unchanged; no age-clock retraining", calibration_formula = "predicted_age ~ age",
      ageacc_formula = "predicted_age - (intercept + slope * age)", standardization_formula = "(ageacc - mean(ageacc in the calibration population)) / sd(ageacc in the calibration population)",
      residual_sd_definition = "Sample standard deviation with denominator n_calibration - 1", ageacc_reference =
      "Chia internal all exact-age participants", standardization_reference = "Chia internal all exact-age participants",
      group_specific_calibration = FALSE, diagnosis_in_calibration = FALSE, ukb_calibration_used = FALSE,
      weighting = "Unweighted OLS over eligible sample rows; unique participant/unit count reported separately",
      n_clocks = 59L, n_samples_per_clock = length(sample_set), residual_sd_tolerance_years = sd_tolerance,
      audit_moment_tolerance = moment_tolerance, noneligible_or_nonfinite_pairs = "Retained with NA ageacc and ageacc_z")
    list(long = result_long, calibration = calibration, audit = audit, settings = settings)
  }

  # Synchronize Chia analytical results
  # Replace all Chia age-acceleration consumers together after internal calibration.
  f7_update_chia_results <- function(prepared, internal, meta) {
    if (!is.list(internal) || !all(c("long", "calibration", "audit", "settings") %in% names(internal))) stop("Chia internal-calibration output is incomplete.")
    ann <- prepared$annotation;
    map <- prepared$clock_map
    f7_need(map, c("clock_id", "reactome_id"), "Prepared clock identity map")
    if (nrow(map) != 59L || anyDuplicated(map$clock_id) || anyDuplicated(map$reactome_id)) stop("Expected 59 unique clock IDs before updating Chia results.")
    disease <- f7_chia_associations(internal$long, meta)
    for (nm in c("associations", "distribution")) {
      j <- match(disease[[nm]]$clock_id, map$clock_id)
      if (anyNA(j)) stop("Unknown Chia clock ID in ", nm)
      disease[[nm]]$reactome_id <- map$reactome_id[j]
      disease[[nm]] <- f7_attach_annotation(disease[[nm]], ann)
    }
    if (nrow(disease$associations) != 118L) stop("Chia update must retain both contrasts for all 59 clocks.")
    counts <- lapply(unique(disease$associations$contrast), function(contrast) {
      z <- disease$associations[disease$associations$contrast == contrast, , drop = FALSE]
      row <- data.frame(contrast = contrast, case_group = z$case_group[1L], reference_group = z$reference_group[1L],
        n_clocks = nrow(z), stringsAsFactors = FALSE)
      for (nm in c("n", "n_case", "n_control", "n_units", "n_case_units", "n_control_units")) {
        row[[paste0(nm, "_min")]] <- min(z[[nm]])
        row[[paste0(nm, "_max")]] <- max(z[[nm]])
      }
      row
    })
    prepared$associations <- disease$associations
    prepared$contrast_counts <- f7_bind(counts)
    prepared$chia_internal_ageacc <- disease$distribution
    prepared$tables$Chia_disease_associations <- disease$associations[, c("contrast", "reactome_id", "full_name",
      "clock_type", "class_name", "n", "n_case", "n_control", "n_units", "n_case_units", "n_control_units",
      "OR_per_chia_sd", "CI95_low", "CI95_high", "p", "p_bonf_global"), drop = FALSE]
    calibration <- as.data.frame(internal$calibration, stringsAsFactors = FALSE)
    j <- match(calibration$clock_id, map$clock_id)
    if (nrow(calibration) != 59L || anyNA(j) || anyDuplicated(j)) stop("Calibration must identify all 59 clocks exactly once.")
    calibration$reactome_id <- map$reactome_id[j]
    calibration <- f7_attach_annotation(calibration, ann);
    calibration$cohort <- "Chia"
    prepared$tables$Chia_ageacc_calibration <- calibration[, c("cohort", "reactome_id", "full_name", "clock_type",
      "class_name", "n_calibration", "n_calibration_units", "intercept", "slope", "residual_mean_years", "residual_sd_years"),
      drop = FALSE]
    prepared$audit$chia_internal_calibration <- internal$audit
    prepared$audit$disease_covariates <- disease$audit
    prepared$audit$disease_warnings <- disease$warnings
    prepared$audit$diagnosis_covariate_counts <- disease$contingency
    prepared$audit$disease_status <- disease$associations[, c("contrast", "clock_id", "reactome_id", "status",
      "n_eligible", "n", "n_case", "n_control", "n_units", "n_case_units", "n_control_units", "n_missing",
      "converged", "se_method", "covariates_used", "formula", "design_rank", "design_ncol", "fitted_probability_min",
      "fitted_probability_max", "warning"), drop = FALSE]

    # Keep the original ages/predictions, while updating any attached AA.
    ch <- which(prepared$prediction$cohort == "Chia")
    overall <- internal$long[internal$long$clock_type == "overall", , drop = FALSE]
    j <- match(prepared$prediction$sample_id[ch], overall$sample_id)
    if (anyNA(j)) stop("Selected Chia samples are absent from internal calibration results.")
    for (nm in intersect(c("ageacc", "ageacc_z"), names(prepared$prediction))) prepared$prediction[[nm]][ch] <-
    overall[[nm]][j]
    for (nm in names(internal$settings)) prepared$settings[[nm]] <- internal$settings[[nm]]
    prepared$settings$age_acceleration_reference_version <- "chia_internal_all_exact_v1"
    prepared$settings$age_acceleration <- paste("Chia: for each of 59 clocks, fit predicted age ~ chronological age using all eligible exact-age Chia samples together;",
      "AA is the regression residual; Z = (AA - mean_Chia(AA)) / SD_Chia(AA).", "The same Chia parameters apply to healthy controls, ALS and other neurological disease samples.",
      "UKB-trained age predictions are unchanged. Gong age acceleration is outside the Chia disease analysis.")
    prepared$settings$disease_effect_scale <- "Adjusted odds ratio per one Chia-cohort residual SD"
    prepared$settings$chia_calibration_population <- "All eligible exact-age Chia samples with finite clock predictions, before disease-model covariate complete-case filtering; not stratified by diagnosis."
    prepared
  }

  # All-age performance uses the unchanged frozen age predictions. Chia's
  # disease-analysis age acceleration is then recalibrated within Chia.
  post <- postprocess_external(result, B = POST_BOOT, seed = S)
  f7_validate_chia_frozen_records(post$external$chia_als$ageacc, post$external$chia_als$prediction, post$clock_map)
  internal <- f7_chia_internal_calibration(post$external$chia_als$ageacc)
  post_for_analysis <- post
  post_for_analysis$external$chia_als$ageacc <- internal$long
  wide <- as.data.frame(post$external$chia_als$prediction, check.names = FALSE)
  overall <- internal$long[internal$long$clock_type == "overall", , drop = FALSE]
  i <- match(wide$sample_id, overall$sample_id)
  if (anyNA(i)) stop("Chia wide/long sample IDs disagree.")
  wide$age_acceleration <- overall$ageacc[i]
  wide$age_acceleration_z <- overall$ageacc_z[i]
  post_for_analysis$external$chia_als$prediction <- wide
  post_for_analysis$chia_internal_calibration <- internal$calibration
  post_for_analysis$config$ageacc <- "Chia: pooled eligible-age internal calibration and internal residual SD."
  prepared <- f7_prepare(post_for_analysis, result, annotation)
  prepared$annotation <- annotation
  prepared <- f7_add_contribution(prepared, post, result, model_bundle, ukb_contribution)
  prepared <- f7_update_chia_results(prepared, internal, post$external$chia_als$meta)

  # The original source reports contribution ranks and top-10 overlap, without
  # fitting a separate rank-correlation test. Contributions use predicted pathway
  # ages directly, not the internally recalibrated Chia age accelerations.
  list(raw = result, model_bundle = model_bundle, post = post_for_analysis, analysis = prepared, chia_calibration =
    internal)
}
