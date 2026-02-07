#!/usr/bin/env Rscript
# Step11_GLM_Ranch_AllModels.R
# Ranch-level exploratory models:
#   1) Merge community and environment
#   2) Alpha diversity (richness, Shannon, evenness)
#   3) NMDS + envfit
#   4) Simple GLMs (alpha and NMDS1 vs env predictors)
#   5) ManyGLM (optional, if mvabund is installed)
#   6) Mantel test between Bray-Curtis and env distance
#
# Inputs:
#   --comm   : taxa_matrix_by_ranch_SUM_TRANS.csv (rows = ranches, cols = taxa, plus sample_id)
#   --env    : ranch_env_HIRES.csv (one row per ranch, from Step11b_Aggregate_HIRES_to_Ranch)
#   --outdir : output directory
#
# Notes:
#   - This is exploratory. N = 4 ranches, so all inferential results must be interpreted cautiously.
#   - Script will correct the El_Mszquite sample_id typo if present.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a) > 0) a else b
cli_i  <- function(...) cat("[INFO]", sprintf(...), "\n")
cli_ok <- function(...) cat("[OK]",   sprintf(...), "\n")
cli_w  <- function(...) cat("[WARN]", sprintf(...), "\n")

# ---------------------------
# CLI options
# ---------------------------
opt_list <- list(
  make_option(c("--comm"),   type = "character",
              help = "Ranch community matrix CSV (rows = ranches, cols = taxa, plus sample_id)."),
  make_option(c("--env"),    type = "character",
              help = "Ranch environment CSV (ranch_env_HIRES.csv)."),
  make_option(c("--outdir"), type = "character",
              help = "Output directory."),
  make_option(c("--id-col"), type = "character", default = "sample_id",
              help = "Sample ID column name (default = sample_id)."),
  make_option(c("--env-cols"), type = "character",
              default = "ElevDEM,slope,roughness,topo_res,EVI_amp,EVI_area,EVI_greenup,TmeanK,TmaxK",
              help = "Comma separated list of environmental variables to use.")
)

opt <- parse_args(OptionParser(option_list = opt_list))

ID_COL   <- opt$`id-col`
ENV_COLS <- strsplit(opt$`env-cols`, ",", fixed = TRUE)[[1]] |> trimws()

# If user does not pass paths, you can set defaults here (edit if needed)
COMM_CSV <- opt$comm %||% "taxa_matrix_by_ranch_SUM_TRANS.csv"
ENV_CSV  <- opt$env  %||% "Step11a_HIRES_RANCH/ranch_env_HIRES.csv"
OUTDIR   <- opt$outdir %||% "Step11_AllModels_RANCH"

if (!file.exists(COMM_CSV)) stop("Community file not found: ", COMM_CSV)
if (!file.exists(ENV_CSV))  stop("Env file not found: ", ENV_CSV)
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

cli_i("comm   : %s", COMM_CSV)
cli_i("env    : %s", ENV_CSV)
cli_i("outdir : %s", OUTDIR)
cli_i("ID_COL : %s", ID_COL)
cli_i("ENV_COLS: %s", paste(ENV_COLS, collapse = ", "))

# ---------------------------
# Load data
# ---------------------------
comm <- read_csv(COMM_CSV, show_col_types = FALSE)
env  <- read_csv(ENV_CSV,  show_col_types = FALSE)

if (!ID_COL %in% names(comm)) stop("Community file must have ID column: ", ID_COL)
if (!ID_COL %in% names(env))  stop("Env file must have ID column: ", ID_COL)

# Clean IDs and fix known typo
comm[[ID_COL]] <- gsub("\\s+", "_", comm[[ID_COL]])
env[[ID_COL]]  <- gsub("\\s+", "_", env[[ID_COL]])

if ("El_Mszquite" %in% comm[[ID_COL]]) {
  cli_w("Fixing sample_id typo: El_Mszquite -> El_Mesquite in community table")
  comm[[ID_COL]][comm[[ID_COL]] == "El_Mszquite"] <- "El_Mesquite"
}

# Align to shared ranch IDs
ids_comm <- unique(comm[[ID_COL]])
ids_env  <- unique(env[[ID_COL]])
ids_keep <- intersect(ids_comm, ids_env)

cli_i("Ranches in community: %s", paste(ids_comm, collapse = ", "))
cli_i("Ranches in env      : %s", paste(ids_env, collapse = ", "))
cli_i("Overlap ranch IDs   : %s", paste(ids_keep, collapse = ", "))

if (length(ids_keep) < 4) {
  cli_w("Fewer than 4 overlapping ranches. Continuing, but results will be very limited.")
}

comm <- comm[comm[[ID_COL]] %in% ids_keep, ]
env  <- env[ env[[ID_COL]] %in% ids_keep, ]

# ---------------------------
# Community matrix and alpha diversity
# ---------------------------
taxa_cols <- setdiff(names(comm), ID_COL)
if (length(taxa_cols) < 2) stop("Need at least 2 taxa columns in community file.")

comm_taxa <- comm[, taxa_cols]
comm_taxa[] <- lapply(comm_taxa, function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
})

mat <- as.matrix(comm_taxa)
rownames(mat) <- comm[[ID_COL]]

# Richness, Shannon, evenness
S        <- rowSums(mat > 0)
Shannon  <- diversity(mat, index = "shannon")
Evenness <- ifelse(S > 0, Shannon / log(S), NA_real_)

alpha_df <- tibble(
  !!ID_COL := comm[[ID_COL]],
  Richness = S,
  Shannon  = Shannon,
  Evenness = Evenness
)

alpha_out <- file.path(OUTDIR, "Ranch_alpha_diversity.csv")
write_csv(alpha_df, alpha_out)
cli_ok("Wrote alpha diversity table: %s", alpha_out)

# ---------------------------
# NMDS and envfit
# ---------------------------
cli_i("Running NMDS (k = 2) on Bray Curtis at ranch level...")
set.seed(42)

# Run NMDS on Bray–Curtis
nmds <- metaMDS(
  mat,
  distance = "bray",
  k        = 2,
  trymax   = 50,
  autotransform = FALSE,
  trace    = FALSE
)

# Extract *site* scores only (one row per ranch)
scr <- scores(nmds, display = "sites")

# Turn into data.frame, but DO NOT assume axis names like "MDS1"/"MDS2"
nmds_scores <- as.data.frame(scr)

# Capture whatever the axis names actually are (e.g., "NMDS1", "NMDS2" or "MDS1", "MDS2")
axis_names <- colnames(nmds_scores)

if (length(axis_names) < 2) {
  stop("NMDS returned fewer than 2 axes – cannot build NMDS1/NMDS2 table.")
}

# Add the ranch ID column from rownames
nmds_scores[[ID_COL]] <- rownames(nmds_scores)

# Build a clean NMDS scores tibble:
#   sample_id, NMDS1, NMDS2
# We explicitly pull the first two axes whatever they are called.
nmds_df <- tibble(
  !!ID_COL := nmds_scores[[ID_COL]],
  NMDS1    = nmds_scores[[axis_names[1]]],
  NMDS2    = nmds_scores[[axis_names[2]]]
)

# Save NMDS scores
nmds_out <- file.path(OUTDIR, "Ranch_NMDS_scores.csv")
write_csv(nmds_df, nmds_out)
cli_ok("Wrote NMDS scores: %s", nmds_out)

# Prepare env for envfit: subset to ENV_COLS
env_sub <- env %>%
  select(all_of(c(ID_COL, ENV_COLS)))

# Order env rows to match NMDS rownames
env_sub <- env_sub[match(rownames(mat), env_sub[[ID_COL]]), ]
env_mat <- as.data.frame(env_sub[, ENV_COLS, drop = FALSE])

fit <- envfit(nmds, env_mat, permutations = 999)

# Save envfit summary
ef_out <- file.path(OUTDIR, "Envfit_summary.txt")
sink(ef_out)
cat("Envfit results (ranch level):\n\n")
print(fit)
sink()
cli_ok("Wrote envfit summary: %s", ef_out)

# NMDS plot with envfit vectors
pdf(file.path(OUTDIR, "Fig_Ranch_NMDS_envfit.pdf"), width = 6.5, height = 5.5)
plot(nmds, type = "n")
points(nmds, display = "sites", pch = 19, col = "black")
text(nmds, display = "sites", labels = rownames(mat), pos = 3)
plot(fit, p.max = 0.1, col = "red")  # show only env vars with p <= 0.1
title("Ranch NMDS with envfit vectors (exploratory)")
dev.off()
cli_ok("Wrote NMDS envfit figure.")

# ---------------------------
# Merge env + alpha + NMDS into one table
# ---------------------------
merged <- env_sub %>%
  left_join(alpha_df, by = ID_COL) %>%
  left_join(nmds_df, by = ID_COL)

merged_out <- file.path(OUTDIR, "Ranch_Env_Alpha_NMDS.csv")
write_csv(merged, merged_out)
cli_ok("Wrote merged ranch table: %s", merged_out)

# ---------------------------
# Simple GLMs (very small N, fully exploratory)
# ---------------------------
glm_results <- list()

fit_if_possible <- function(formula_str, data, name) {
  form <- as.formula(formula_str)
  vars <- all.vars(form)
  if (!all(vars %in% names(data))) {
    cli_w("Skipping model '%s' - missing vars: %s",
          name, paste(setdiff(vars, names(data)), collapse = ", "))
    return(NULL)
  }
  dsub <- data[, vars, drop = FALSE]
  dsub <- dsub[complete.cases(dsub), ]
  if (nrow(dsub) < 4) {
    cli_w("Skipping model '%s' - fewer than 4 complete cases.", name)
    return(NULL)
  }
  m <- lm(form, data = dsub)
  m
}

# Choose simple one-predictor models to avoid overfitting
glm_results[["Shannon_vs_Elev"]]   <- fit_if_possible("Shannon ~ scale(ElevDEM)",      merged, "Shannon_vs_Elev")
glm_results[["Shannon_vs_EVIarea"]] <- fit_if_possible("Shannon ~ scale(EVI_area)",    merged, "Shannon_vs_EVIarea")
glm_results[["Richness_vs_Elev"]]  <- fit_if_possible("Richness ~ scale(ElevDEM)",     merged, "Richness_vs_Elev")
glm_results[["NMDS1_vs_Elev"]]     <- fit_if_possible("NMDS1 ~ scale(ElevDEM)",        merged, "NMDS1_vs_Elev")
glm_results[["NMDS1_vs_EVIarea"]]  <- fit_if_possible("NMDS1 ~ scale(EVI_area)",       merged, "NMDS1_vs_EVIarea")

summary_rows <- list()

for (nm in names(glm_results)) {
  m <- glm_results[[nm]]
  if (is.null(m)) next

  sum_file <- file.path(OUTDIR, paste0("GLM_", nm, "_summary.txt"))
  sink(sum_file)
  cat("Model:", nm, "(ranch level, exploratory)\n\n")
  print(summary(m))
  sink()

  s <- summary(m)
  # simple approximate R2
  r2 <- 1 - sum(s$residuals^2) / sum((s$model[[1]] - mean(s$model[[1]]))^2)

  coefs <- as.data.frame(coef(summary(m)))
  coefs$term <- rownames(coefs)

  for (i in seq_len(nrow(coefs))) {
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      model    = nm,
      term     = coefs$term[i],
      estimate = coefs$Estimate[i],
      se       = coefs$`Std. Error`[i],
      t_value  = coefs$`t value`[i],
      p_value  = coefs$`Pr(>|t|)`[i],
      n        = length(s$residuals),
      r2       = r2
    )
  }
}

if (length(summary_rows) > 0) {
  glm_summary <- bind_rows(summary_rows)
  glm_out <- file.path(OUTDIR, "GLM_summary_table.csv")
  write_csv(glm_summary, glm_out)
  cli_ok("Wrote GLM summary table: %s", glm_out)
} else {
  cli_w("No GLMs were successfully fitted.")
}

# ---------------------------
# ManyGLM (optional, mvabund)
# ---------------------------
if (requireNamespace("mvabund", quietly = TRUE)) {
  cli_i("mvabund found - running ManyGLM (exploratory, ElevDEM as predictor)...")
  library(mvabund)

  elev_vec <- env_sub$ElevDEM[match(rownames(mat), env_sub[[ID_COL]])]
  if (any(is.na(elev_vec))) {
    cli_w("Missing ElevDEM for some ranches, skipping ManyGLM.")
  } else {
    y_mv <- mvabund(mat)
    dat_mv <- data.frame(ElevDEM = scale(elev_vec))

    mg <- manyglm(y_mv ~ ElevDEM, data = dat_mv, family = "poisson")
    mg_sum <- summary(mg, test = "wald")

    sink(file.path(OUTDIR, "ManyGLM_ElevDEM_summary.txt"))
    cat("ManyGLM ~ ElevDEM (ranch level, exploratory)\n\n")
    print(mg_sum)
    sink()
    cli_ok("Wrote ManyGLM summary for ElevDEM.")
  }
} else {
  cli_w("Package mvabund not installed - skipping ManyGLM step.")
}

# ---------------------------
# Mantel test: Bray Curtis vs env distance
# ---------------------------
cli_i("Running Mantel test (Bray Curtis vs env distance on scaled predictors)...")

# Bray Curtis dissimilarity on community
bc <- vegdist(mat, method = "bray")

# Euclidean distance on scaled env predictors
env_scaled <- merged %>%
  select(all_of(ENV_COLS)) %>%
  scale()

env_dist <- dist(env_scaled)

mant <- mantel(bc, env_dist, permutations = 999)

mant_out <- file.path(OUTDIR, "Mantel_summary.txt")
sink(mant_out)
cat("Mantel test: Bray Curtis vs env distance (ranch level, exploratory)\n\n")
print(mant)
sink()
cli_ok("Wrote Mantel summary: %s", mant_out)

cli_ok("Ranch-level all models complete. Outputs in: %s", OUTDIR)