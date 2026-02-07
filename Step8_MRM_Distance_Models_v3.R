#!/usr/bin/env Rscript
# Step8_MRM_Distance_Models_v3.R
# Multiple Regression on distance matrices (MRM) with robust checks + QC plots

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ecodist)
  library(ggplot2)
})

opt_list <- list(
  make_option("--dist_dir",  type="character", help="Folder with Step7 distance CSVs"),
  make_option("--outdir",    type="character", default="Step8_MRM", help="Output folder"),
  make_option("--perms",     type="integer",   default=9999, help="MRM permutations"),
  make_option("--standardize", action="store_true", default=FALSE,
              help="Z-score predictors (vectorized distances) before MRM"),
  make_option("--method",    type="character", default="linear",
              help="MRM method: 'linear' (default) or 'logistic'")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cli <- function(t, ...) cat(sprintf("[%s] %s\n", t, sprintf(...)))

read_mat <- function(path) {
  df <- suppressWarnings(read.csv(path, check.names=FALSE, row.names=1))
  as.matrix(df)
}
upper_vec <- function(M) M[upper.tri(M, diag = FALSE)]
safe_sd <- function(x) { s <- sd(x); if (is.na(s) || s == 0) 1 else s }

if (is.null(opt$dist_dir) || !dir.exists(opt$dist_dir)) {
  stop("dist_dir not found: ", opt$dist_dir)
}
expected <- c(
  "D_bray.csv",
  "D_geo_km.csv",
  "D_env_scaled.csv",
  "D_cost_res_slope_global.csv",
  "D_cost_res_rough_global.csv",
  "D_cost_combined_global.csv"
)
paths <- file.path(opt$dist_dir, expected)
names(paths) <- sub("\\.csv$", "", expected)
missing <- names(paths)[!file.exists(paths)]
if (length(missing)) stop("Missing distance matrices in dist_dir: ", paste(missing, collapse=", "))
mats <- lapply(paths, read_mat)

rn0 <- rownames(mats[[1]]); cn0 <- colnames(mats[[1]])
for (nm in names(mats)) {
  if (!identical(rownames(mats[[nm]]), rn0) || !identical(colnames(mats[[nm]]), cn0)) {
    stop("Rownames/colnames mismatch in matrix: ", nm)
  }
}
sites <- rn0
cli("OK", "Sites (k=%d): %s", length(sites), paste(sites, collapse=", "))

D_comm <- mats[["D_bray"]]
X <- list()
for (nm in names(mats)) if (nm != "D_bray") X[[nm]] <- upper_vec(mats[[nm]])
y <- upper_vec(D_comm)
if (opt$standardize) X <- lapply(X, function(v) (v - mean(v)) / safe_sd(v))

if (!(opt$method %in% c("linear","logistic"))) stop("Use --method linear or logistic.")

run_mrm1 <- function(yM, xM, nperm = 9999, method = "linear") {
  fit <- ecodist::MRM(as.dist(yM) ~ as.dist(xM),
                      nperm  = nperm,
                      method = method)
  cf <- as.data.frame(fit$coef)
  cf$term <- rownames(fit$coef)
  names(cf) <- tolower(gsub("\\W+", "", names(cf)))
  coef_candidates <- c("coef", "estimate", "b", "beta")
  p_candidates    <- c("pval", "p", "pvalue", "prgtf")
  coef_col <- coef_candidates[coef_candidates %in% names(cf)][1]
  p_col    <- p_candidates[p_candidates %in% names(cf)][1]
  if (is.na(coef_col) || is.na(p_col)) {
    num_cols <- names(cf)[vapply(cf, is.numeric, logical(1))]
    if (length(num_cols) >= 2) {
      coef_col <- num_cols[1]
      p_col    <- num_cols[2]
    } else {
      stop("Could not identify coefficient and p-value columns in MRM output.")
    }
  }
  cf <- cf %>% dplyr::filter(term != "(Intercept)")
  if (nrow(cf) != 1L) {
    warning("Single-predictor MRM returned an unexpected number of terms; using first non-intercept row.")
    cf <- cf[1, , drop = FALSE]
  }
  cf[[coef_col]] <- suppressWarnings(as.numeric(cf[[coef_col]]))
  cf[[p_col]]    <- suppressWarnings(as.numeric(cf[[p_col]]))
  tibble::tibble(
    term = cf$term,
    coef = cf[[coef_col]],
    R2   = fit$r.squared,
    p    = cf[[p_col]]
  )
}

single_results <- list()
k <- length(sites)
for (nm in names(X)) {
  xM <- matrix(0, nrow = k, ncol = k,
               dimnames = list(sites, sites))
  xM[upper.tri(xM)] <- X[[nm]]
  xM <- xM + t(xM)
  res <- run_mrm1(D_comm, xM, nperm = opt$perms, method = opt$method)
  res$predictor <- nm
  single_results[[nm]] <- res
}
single_df <- dplyr::bind_rows(single_results) %>%
  dplyr::select(predictor, term, coef, R2, p)

Xd <- as.data.frame(X, check.names = FALSE)

run_mrm_multi <- function(yM, Xdf, formula_str, nperm, method="linear") {
  terms_rhs <- all.vars(delete.response(terms(as.formula(formula_str))))
  Xsub <- Xdf[, terms_rhs, drop=FALSE]
  fit <- ecodist::MRM(as.dist(yM) ~ ., data = Xsub, nperm = nperm, method = method)
  cf <- as.data.frame(fit$coef)
  cf$term <- rownames(fit$coef)
  names(cf) <- tolower(gsub("\\W+", "", names(cf)))
  coef_col <- c("coef","estimate","b","beta")
  p_col    <- c("pval","p","pvalue","prgtf")
  coef_col <- coef_col[coef_col %in% names(cf)][1]
  p_col    <- p_col[p_col %in% names(cf)][1]
  if (is.na(coef_col) || is.na(p_col)) {
    nums <- names(cf)[vapply(cf, is.numeric, logical(1))]
    if (length(nums) >= 2) { coef_col <- nums[1]; p_col <- nums[2] }
  }
  cf <- cf %>% dplyr::filter(term != "(Intercept)")
  if (nrow(cf) == 0) {
    terms_txt <- "(no terms)"
  } else {
    cf[[coef_col]] <- suppressWarnings(as.numeric(cf[[coef_col]]))
    cf[[p_col]]    <- suppressWarnings(as.numeric(cf[[p_col]]))
    terms_txt <- paste(sprintf("%s (coef=%.3g, p=%.3g)", cf$term, cf[[coef_col]], cf[[p_col]]),
                       collapse = " | ")
  }
  tibble::tibble(model = formula_str, R2 = fit$r.squared, terms = terms_txt)
}

multi_models <- list(
  "as.dist(D_comm) ~ D_geo_km",
  "as.dist(D_comm) ~ D_env_scaled",
  "as.dist(D_comm) ~ D_cost_res_slope_global",
  "as.dist(D_comm) ~ D_cost_res_rough_global",
  "as.dist(D_comm) ~ D_cost_combined_global",
  "as.dist(D_comm) ~ D_geo_km + D_env_scaled",
  "as.dist(D_comm) ~ D_geo_km + D_cost_combined_global",
  "as.dist(D_comm) ~ D_geo_km + D_env_scaled + D_cost_combined_global"
)
multi_df <- dplyr::bind_rows(lapply(multi_models, function(fml) {
  run_mrm_multi(D_comm, Xd, fml, nperm=opt$perms, method=opt$method)
}))

write.csv(single_df, file.path(opt$outdir, "MRM_single_predictor.csv"), row.names = FALSE)
write.csv(multi_df,  file.path(opt$outdir, "MRM_multi_predictor.csv"),  row.names = FALSE)

qc_dir <- file.path(opt$outdir, "QC"); dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
yvec <- upper_vec(D_comm)
for (nm in names(X)) {
  df <- data.frame(D_comm = yvec, predictor = X[[nm]])
  p <- ggplot(df, aes(predictor, D_comm)) +
    geom_point(size=3, alpha=0.8) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.8, color="black") +
    theme_minimal(base_size = 14) +
    labs(title = paste("D_bray vs", nm), x = nm, y = "D_bray")
  ggsave(file.path(qc_dir, paste0("scatter_Dbray_vs_", nm, ".png")), p,
         width=6.5, height=5.5, dpi=220)
}

sink(file.path(opt$outdir, "SUMMARY.txt"))
cat("Step 8 — MRM results\n")
cat(sprintf("- Sites (k=%d): %s\n\n", length(sites), paste(sites, collapse=", ")))
cat("Single-predictor MRM (D_bray ~ X)\n"); print(single_df)
cat("\nMulti-predictor MRM\n"); print(multi_df)
sink()

cli("OK", "MRM complete. Results in %s", opt$outdir)