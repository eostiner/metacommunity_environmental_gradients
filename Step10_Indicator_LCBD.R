#!/usr/bin/env Rscript
# Step10_Indicator_LCBD.R (cleaned + robust)
# Indicator taxa (IndVal) and site contributions (LCBD/SCBD)
# Works with:
#   (A) taxa-as-rows:  first column = Taxon, remaining = samples (your hellinger_matrix.csv)
#   (B) samples-as-rows: row = sample, first/explicit column = sample ID

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(vegan)
  library(indicspecies)  # multipatt
  library(adespatial)    # beta.div for LCBD/SCBD
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

cli <- function(level, ...) cat(sprintf("[%s] %s\n", level, sprintf(...)))

# ---------- CLI ----------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--comm_file",  type="character",
              help="Community matrix CSV. Either taxa-as-rows (Taxon + samples) or samples-as-rows."),
  make_option("--meta_file",  type="character",
              help="Metadata CSV (one row per sample)."),
  make_option("--group_col",  type="character",
              help="Metadata column defining groups (e.g., ranch)."),
  make_option("--id_col",     type="character", default="sample_id",
              help="Metadata column with sample IDs [default: %default]."),
  make_option("--transform",  type="character", default="hellinger",
              help="One of: hellinger | pa | raw (applied to community before LCBD/SCBD)."),
  make_option("--indval_mode", type="character", default="pa",
              help="IndVal mode: 'pa' (presence/absence) or 'r.g' (abundance)."),
  make_option("--top_k",      type="integer", default=25,
              help="Top-K indicator taxa to write separately."),
  make_option("--outdir",     type="character",
              help="Output directory.")
)))

if (is.null(opt$comm_file) || is.null(opt$meta_file) ||
    is.null(opt$group_col) || is.null(opt$outdir)) {
  stop("Required: --comm_file, --meta_file, --group_col, --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
norm_ids <- function(x){
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

as_pa <- function(M){
  M[] <- ifelse(M > 0, 1, 0)
  M
}

# ---------- Load community (auto-detect orientation) ----------
cli("INFO", "Reading community matrix: %s", opt$comm_file)
comm_raw <- readr::read_csv(opt$comm_file, show_col_types = FALSE)
comm_raw <- as.data.frame(comm_raw, check.names = FALSE)

cols_lower <- tolower(names(comm_raw))

if ("taxon" %in% cols_lower) {
  # ---- Case A: taxa-as-rows (your hellinger_matrix.csv format) ----
  tax_col_name <- names(comm_raw)[which(cols_lower == "taxon")[1]]
  cli("INFO", "Detected taxa-as-rows format using '%s' as Taxon column.", tax_col_name)

  taxa <- as.character(comm_raw[[tax_col_name]])
  mat  <- comm_raw[, setdiff(names(comm_raw), tax_col_name), drop = FALSE]

  # numeric matrix: taxa x samples
  mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
  mat[is.na(mat)] <- 0

  # sample IDs from column names → normalized
  sample_ids <- norm_ids(colnames(mat))

  # transpose → samples x taxa
  mat_t <- t(as.matrix(mat))
  colnames(mat_t) <- taxa
  rownames(mat_t) <- sample_ids

  comm <- as.data.frame(mat_t, check.names = FALSE)

} else {
  # ---- Case B: samples-as-rows ----
  cli("INFO", "No 'Taxon' column found; assuming samples-as-rows community matrix.")

  id_col_in_comm <- NULL

  if (opt$id_col %in% names(comm_raw)) {
    id_col_in_comm <- opt$id_col
  } else if ("sample_id" %in% names(comm_raw)) {
    id_col_in_comm <- "sample_id"
  } else {
    # fall back to first column if it looks non-numeric
    if (!is.numeric(comm_raw[[1]])) {
      id_col_in_comm <- names(comm_raw)[1]
    }
  }

  if (is.null(id_col_in_comm)) {
    stop("Could not identify a sample ID column in the community CSV. ",
         "Please include a 'sample_id' column or specify --id_col that exists there.")
  }

  cli("INFO", "Detected samples-as-rows format with sample IDs in '%s'.", id_col_in_comm)

  # normalize ID column name to 'sample_id'
  if (id_col_in_comm != "sample_id") {
    names(comm_raw)[names(comm_raw) == id_col_in_comm] <- "sample_id"
  }

  rn <- norm_ids(as.character(comm_raw[["sample_id"]]))
  comm <- comm_raw[, setdiff(names(comm_raw), "sample_id"), drop = FALSE]
  comm[] <- lapply(comm, function(x) suppressWarnings(as.numeric(x)))
  comm[is.na(comm)] <- 0
  rownames(comm) <- rn
}

cli("INFO", "Community matrix: %d samples × %d taxa", nrow(comm), ncol(comm))

# ---------- Load metadata ----------
cli("INFO", "Reading metadata: %s", opt$meta_file)
meta <- readr::read_csv(opt$meta_file, show_col_types = FALSE)

if (!(opt$id_col %in% names(meta))) {
  stop(sprintf("id_col '%s' not found in metadata.", opt$id_col))
}
if (!(opt$group_col %in% names(meta))) {
  stop(sprintf("group_col '%s' not found in metadata.", opt$group_col))
}

meta[[opt$id_col]]    <- norm_ids(as.character(meta[[opt$id_col]]))
meta[[opt$group_col]] <- as.factor(meta[[opt$group_col]])

# ---------- Align sample IDs ----------
ov <- intersect(rownames(comm), meta[[opt$id_col]])

cli("INFO", "Samples in community: %d; in metadata: %d; overlap: %d",
    nrow(comm), nrow(meta), length(ov))

if (length(ov) < 4) {
  cat("\nExamples not in metadata (first 20):\n")
  print(head(setdiff(rownames(comm), meta[[opt$id_col]]), 20))
  cat("\nExamples not in community (first 20):\n")
  print(head(setdiff(meta[[opt$id_col]], rownames(comm)), 20))
  stop("Fewer than 4 overlapping samples between community and metadata.")
}

# reorder / subset
comm <- comm[ov, , drop = FALSE]
meta <- meta[match(ov, meta[[opt$id_col]]), , drop = FALSE]
stopifnot(identical(rownames(comm), meta[[opt$id_col]]))

# ---------- Transform for analysis (for LCBD/SCBD) ----------
if (tolower(opt$transform) == "hellinger") {
  cli("INFO", "Applying Hellinger transform.")
  comm_tr <- vegan::decostand(comm, method = "hellinger")
} else if (tolower(opt$transform) == "pa") {
  cli("INFO", "Using presence/absence for LCBD/SCBD.")
  comm_tr <- as_pa(comm)
} else if (tolower(opt$transform) == "raw") {
  cli("INFO", "Using raw (untransformed) abundances for LCBD/SCBD.")
  comm_tr <- comm
} else {
  stop("transform must be one of: hellinger | pa | raw")
}

# Matrix for IndVal (PA by default)
X_ind <- if (tolower(opt$indval_mode) == "pa") {
  cli("INFO", "Using presence/absence for IndVal.")
  as_pa(comm)
} else {
  cli("INFO", "Using abundance-based IndVal (r.g).")
  comm
}

group <- droplevels(meta[[opt$group_col]])
if (nlevels(group) < 2) {
  stop("group_col must have at least 2 groups.")
}

# ---------- Indicator taxa (IndVal) ----------
cli("INFO", "Running IndVal (mode=%s) for %d taxa across %d groups",
    opt$indval_mode, ncol(X_ind), nlevels(group))

set.seed(42)
ind_fun <- if (tolower(opt$indval_mode) == "pa") "IndVal.g" else "r.g"

ind <- indicspecies::multipatt(X_ind,
                               cluster = group,
                               func    = ind_fun,
                               control = how(nperm = 999))

sign_tab <- as.data.frame(ind$sign)
sign_tab <- tibble::rownames_to_column(sign_tab, var = "taxon")

# Try to pull a p-value column if present
guess_p <- sign_tab[["p.value"]] %||%
           sign_tab[["p.value."]] %||%
           sign_tab[["p_value"]]

res_out <- sign_tab %>%
  mutate(p_value = guess_p %||% NA_real_) %>%
  arrange(p_value)

readr::write_csv(res_out,
                 file.path(opt$outdir, "indicators_full.csv"))
readr::write_csv(head(res_out, opt$top_k),
                 file.path(opt$outdir, "indicators_topK.csv"))

# ---------- LCBD / SCBD ----------
cli("INFO", "Computing LCBD/SCBD on '%s' matrix", tolower(opt$transform))

bd <- adespatial::beta.div(comm_tr, method = "hellinger")
LCBD <- tibble(sample_id = rownames(comm_tr),
               LCBD      = bd$LCBD)
SCBD <- tibble(taxon = colnames(comm_tr),
               SCBD  = bd$SCBD)

readr::write_csv(LCBD, file.path(opt$outdir, "LCBD.csv"))
readr::write_csv(SCBD, file.path(opt$outdir, "SCBD.csv"))

# ---------- Run metadata ----------
writeLines(c(
  sprintf("comm_file,%s",  opt$comm_file),
  sprintf("meta_file,%s",  opt$meta_file),
  sprintf("id_col,%s",     opt$id_col),
  sprintf("group_col,%s",  opt$group_col),
  sprintf("transform,%s",  opt$transform),
  sprintf("indval_mode,%s",opt$indval_mode),
  sprintf("top_k,%d",      opt$top_k)
), con = file.path(opt$outdir, "RUN_METADATA.txt"))

cli("OK", "Step 10 complete. Outputs in: %s", opt$outdir)