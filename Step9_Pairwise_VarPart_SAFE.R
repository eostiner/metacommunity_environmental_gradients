#!/usr/bin/env Rscript
# Step9_Pairwise_VarPart_SAFE.R
# Minimal Step 9 for n=4 sites with robust R handling:
# - Single-block adj-R² (E, G, optional R)
# - All 2-way varparts (E–G, E–R, G–R) with Cailliez correction
# - Explicit file paths for E/G/R
# - Optional NA filling in R: keep NA, fill with high cost, or constant

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr)
  library(vegan); library(ape)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
cli <- function(level, ...) cat(sprintf("[%s] %s\n", level, sprintf(...)))

opt <- parse_args(OptionParser(option_list=list(
  make_option("--comm_file", type="character", help="CSV square dissimilarity (e.g., D_bray.csv)"),
  make_option("--dist_dir",  type="character", help="Directory (used only for fallback to D_cost_combined_global.csv)"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--method",    type="character", default="bray"),
  make_option("--env_file",  type="character", help="FULL PATH to environmental distance CSV"),
  make_option("--geo_file",  type="character", help="FULL PATH to geographic distance CSV"),
  make_option("--res_file",  type="character", default=NULL, help="(optional) FULL PATH to resistance CSV"),
  make_option("--env_threshold", type="double", default=0.60, help="PCoA cumulative variance threshold per block"),
  make_option("--kmax",          type="integer", default=1, help="Max PCoA axes per block (use 1 for n=4)"),
  make_option("--r_fill", type="character", default="none", help="How to handle NA in R: none | high | const"),
  make_option("--r_fill_mult", type="double", default=1000, help="Multiplier for 'high' mode"),
  make_option("--r_fill_const", type="double", default=1e6, help="Constant for 'const' mode")
)))
stopifnot(!is.null(opt$comm_file), !is.null(opt$outdir), !is.null(opt$env_file), !is.null(opt$geo_file))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
read_square <- function(path){
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  nm <- names(df)
  if (ncol(df) > 1 && (ncol(df)-1) == nrow(df)) {
    rn <- df[[1]]; mat <- as.matrix(df[,-1,drop=FALSE])
    rownames(mat) <- as.character(rn); colnames(mat) <- nm[-1]; return(mat)
  }
  if (nrow(df) == ncol(df)) { mat <- as.matrix(df); rownames(mat) <- colnames(mat); return(mat) }
  stop(sprintf("Not a square distance matrix: %s", path))
}
as_dist <- function(M){ M <- as.matrix(M); diag(M) <- 0; stats::as.dist(M) }
sanity_dist <- function(M, label){
  if (any(is.na(M))) cli("WARN","%s contains NA", label)
  if (!isTRUE(all.equal(M, t(M), tolerance=1e-10)))
    stop(sprintf("%s is not symmetric", label))
  diag(M) <- 0
  M
}
same_ids <- function(A, B, tagA, tagB){
  if (!setequal(rownames(A), rownames(B)))
    stop(sprintf("ID mismatch between %s and %s.\n%s IDs: %s\n%s IDs: %s",
                 tagA, tagB, tagA, paste(rownames(A), collapse=", "),
                 tagB, paste(rownames(B), collapse=", ")))
}
pcoa_block <- function(D, tag, thr=0.60, kmax=1){
  D <- as.matrix(D); diag(D) <- 0
  off <- D[upper.tri(D)]
  if (!any(is.finite(off)) || length(unique(off[is.finite(off)])) < 2)
    stop(sprintf("Degenerate %s (no variation)", tag))
  pc  <- ape::pcoa(as_dist(D), correction = "cailliez")
  eig <- pc$values$Relative_eig; eig[is.na(eig)] <- 0
  k   <- min(kmax, max(1, (which(cumsum(eig) >= thr)[1] %||% length(eig))))
  Ax  <- as.data.frame(pc$vectors[, seq_len(k), drop=FALSE])
  colnames(Ax) <- sprintf("%s_PC%02d", tag, seq_len(k))
  Ax
}
single_adjR2 <- function(D, X){
  if (is.null(X) || ncol(X) == 0) return(NA_real_)
  m <- vegan::capscale(as_dist(D) ~ ., data = X, add = TRUE)
  vegan::RsquareAdj(m)$adj.r.squared
}
pair_varpart <- function(D, X, Y, xname, yname){
  vp <- vegan::varpart(as_dist(D), X, Y)
  data.frame(
    pair        = paste0(xname,"-",yname),
    unique_X    = vp$part$indfract$Adj.R.squared[1],
    unique_Y    = vp$part$indfract$Adj.R.squared[2],
    shared      = vp$part$indfract$Adj.R.squared[3],
    total_adjR2 = vp$part$fract$Adj.R.squared,
    X = xname, Y = yname,
    stringsAsFactors = FALSE
  )
}
first_or_na <- function(x){ if (length(x)==0) NA_real_ else x[1] }

fill_R_na <- function(R, mode = "none", mult = 1000, const = 1e6){
  mode <- tolower(mode)
  n_na <- sum(is.na(R))
  if (n_na == 0 || mode == "none") {
    if (n_na > 0) cli("INFO","R has %d NA; mode=none (keeping NA).", n_na)
    return(list(R=R, n_filled=0, fill_value=NA_real_))
  }
  fin <- R[is.finite(R)]
  if (!length(fin)) {
    cli("WARN","R has no finite values; cannot fill meaningfully. Leaving as NA.")
    return(list(R=R, n_filled=0, fill_value=NA_real_))
  }
  if (mode == "high") {
    fv <- max(fin, na.rm=TRUE) * mult
  } else if (mode == "const") {
    fv <- const
  } else stop("Unknown --r_fill mode. Use: none | high | const")
  R2 <- R; R2[is.na(R2)] <- fv
  cli("INFO","Filled %d NA in R with %s value = %.4g", n_na, mode, fv)
  list(R=R2, n_filled=n_na, fill_value=fv)
}

# ---------- load inputs ----------
cli("INFO","Reading response: %s", opt$comm_file)
D_comm <- read_square(opt$comm_file)
D_comm <- sanity_dist(D_comm, "D_comm")

if (!file.exists(opt$env_file)) stop(sprintf("env_file not found: %s", opt$env_file))
if (!file.exists(opt$geo_file)) stop(sprintf("geo_file not found: %s", opt$geo_file))

E <- read_square(opt$env_file)
G <- read_square(opt$geo_file)
same_ids(D_comm, E, "D_comm", "E"); same_ids(D_comm, G, "D_comm", "G")
E <- sanity_dist(E[rownames(D_comm), rownames(D_comm)], "E")
G <- sanity_dist(G[rownames(D_comm), rownames(D_comm)], "G")

R <- NULL
if (!is.null(opt$res_file) && nzchar(opt$res_file)) {
  if (!file.exists(opt$res_file)) {
    cli("WARN","res_file not found: %s (continuing without R)", opt$res_file)
  } else {
    R <- read_square(opt$res_file)
  }
} else if (!is.null(opt$dist_dir)) {
  cand <- file.path(opt$dist_dir, "D_cost_combined_global.csv")
  if (file.exists(cand)) R <- read_square(cand)
}

R_fill_info <- list(n_filled=0, fill_value=NA_real_, mode=opt$r_fill)
if (!is.null(R)) {
  same_ids(D_comm, R, "D_comm", "R")
  R <- sanity_dist(R[rownames(D_comm), rownames(D_comm)], "R")
  if (any(is.na(R))) {
    Rfi <- fill_R_na(R, mode=opt$r_fill, mult=opt$r_fill_mult, const=opt$r_fill_const)
    R <- Rfi$R; R_fill_info <- c(R_fill_info, Rfi)
  }
}

cli("INFO","Sites: %d; rule-of-thumb max usable predictors per model: %d",
    nrow(D_comm), max(0, nrow(D_comm)-1))

# ---------- PCoA blocks ----------
cli("INFO","Building PCoA blocks (thr=%.2f, kmax=%d)", opt$env_threshold, opt$kmax)
E_pc <- pcoa_block(E, "E", thr=opt$env_threshold, kmax=opt$kmax)
G_pc <- pcoa_block(G, "G", thr=opt$env_threshold, kmax=opt$kmax)
R_pc <- NULL
if (!is.null(R)) {
  R_pc <- try(pcoa_block(R, "R", thr=opt$env_threshold, kmax=opt$kmax), silent=TRUE)
  if (inherits(R_pc, "try-error")) { cli("WARN","Skipping R (degenerate/no variation)"); R_pc <- NULL }
}

# ---------- single-block ----------
single_tbl <- tibble::tibble(
  model  = c("baseline_null","E_only","G_only","R_only"),
  adj_R2 = c(0,
             single_adjR2(D_comm, E_pc),
             single_adjR2(D_comm, G_pc),
             if (!is.null(R_pc)) single_adjR2(D_comm, R_pc) else NA_real_)
)
single_tbl_out <- dplyr::mutate(single_tbl, adj_R2_round = round(adj_R2, 4))
readr::write_csv(single_tbl_out, file.path(opt$outdir, sprintf("single_block_adjR2_%s.csv", opt$method)))

# ---------- pairwise varparts ----------
res_list <- list()
res_list[["E-G"]] <- pair_varpart(D_comm, E_pc, G_pc, "E","G")
if (!is.null(R_pc)) {
  res_list[["E-R"]] <- pair_varpart(D_comm, E_pc, R_pc, "E","R")
  res_list[["G-R"]] <- pair_varpart(D_comm, G_pc, R_pc, "G","R")
}
resdf <- dplyr::bind_rows(res_list)
resdf_out <- dplyr::mutate(resdf, dplyr::across(where(is.numeric), ~round(.x, 4)))
readr::write_csv(resdf_out, file.path(opt$outdir, sprintf("pairwise_varpart_%s.csv", opt$method)))

# ---------- compact consensus ----------
UE <- mean(c(resdf$unique_X[resdf$X=="E"], resdf$unique_Y[resdf$Y=="E"]), na.rm=TRUE)
UG <- mean(c(resdf$unique_X[resdf$X=="G"], resdf$unique_Y[resdf$Y=="G"]), na.rm=TRUE)
UR <- mean(c(resdf$unique_X[resdf$X=="R"], resdf$unique_Y[resdf$Y=="R"]), na.rm=TRUE)

SEG <- first_or_na(resdf$shared[resdf$pair=="E-G"])
SER <- first_or_na(resdf$shared[resdf$pair=="E-R"])
SGR <- first_or_na(resdf$shared[resdf$pair=="G-R"])

summary_tbl <- tibble::tibble(
  U_E=UE, U_R=UR, U_G=UG,
  S_EG=SEG, S_ER=SER, S_GR=SGR,
  total_EG = first_or_na(resdf$total_adjR2[resdf$pair=="E-G"]),
  total_ER = first_or_na(resdf$total_adjR2[resdf$pair=="E-R"]),
  total_GR = first_or_na(resdf$total_adjR2[resdf$pair=="G-R"])
)
summary_tbl_out <- dplyr::mutate(summary_tbl, dplyr::across(where(is.numeric), ~round(.x, 4)))
readr::write_csv(summary_tbl_out, file.path(opt$outdir, sprintf("pairwise_consensus_%s.csv", opt$method)))

# ---------- run metadata ----------
writeLines(c(
  sprintf("method,%s", opt$method),
  sprintf("comm_file,%s", opt$comm_file),
  sprintf("env_file,%s", opt$env_file),
  sprintf("geo_file,%s", opt$geo_file),
  sprintf("res_file,%s", opt$res_file %||% ""),
  sprintf("env_threshold,%.2f", opt$env_threshold),
  sprintf("kmax,%d", opt$kmax),
  sprintf("r_fill,%s", opt$r_fill),
  sprintf("r_fill_mult,%.4g", opt$r_fill_mult),
  sprintf("r_fill_const,%.4g", opt$r_fill_const)
), con=file.path(opt$outdir, "RUN_METADATA.txt"))

cli("OK","SAFE Step 9 complete. Outputs in: %s", opt$outdir)