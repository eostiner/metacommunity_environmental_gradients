#!/usr/bin/env Rscript
# Step5_Preflight_FixIDs.R
# Purpose: Make presence/hellinger/metadata line up on SampleID before Step5 NMDS.
# What it does:
#   • Normalizes ID tokens (space→_, lower vs upper, hyphens/periods)
#   • Detects if matrices have samples as COLUMNS and transposes to rows = samples
#   • Tries light heuristics (strip common prefixes like "Trap_", remove parentheses)
#   • If still no overlap, proposes & applies a suggested ID map via base::adist
#   • Writes fixed files + a report of overlaps
#
# Inputs:
#   --presence, --hellinger, --metadata, --outdir
#
# Outputs (in --outdir):
#   presence_fixed.csv, hellinger_fixed.csv, metadata_fixed.csv
#   suggested_id_map.csv  (only if mapping was used)
#   preflight_report.txt  (counts & the first few mismatches)

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(tidyr); library(stringr)
})

# ---------- CLI ----------
opt <- list(
  make_option("--presence",  type="character"),
  make_option("--hellinger", type="character"),
  make_option("--metadata",  type="character"),
  make_option("--outdir",    type="character")
)
args <- parse_args(OptionParser(option_list = opt))
if (any(sapply(args[c("presence","hellinger","metadata","outdir")], is.null))) {
  stop("Required: --presence --hellinger --metadata --outdir")
}
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helpers ----------
read_auto <- function(p){
  df <- tryCatch(read.csv(p, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) == 1) df <- read.csv(p, sep = "\t", check.names = FALSE)
  df
}

norm_token <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[.\\-]+", "_", x)
  x <- gsub("\\(+[^)]*\\)", "", x)       # remove parenthetical suffixes
  x <- gsub("__+", "_", x)
  x <- gsub("^Trap[_-]?", "", x, ignore.case = TRUE)
  x <- gsub("^Sample[_-]?", "", x, ignore.case = TRUE)
  x <- gsub("^ID[_-]?", "", x, ignore.case = TRUE)
  x
}

norm_ids_col <- function(df){
  df <- as.data.frame(df)
  candidates <- c("SampleID","sample_id","Sample","sample","SampleName","sample_name","Sample_ID","sampleID","ID","id")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) {
    names(df)[names(df) == hit[1]] <- "SampleID"
  } else if (!is.null(rownames(df)) && !all(rownames(df) %in% as.character(seq_len(nrow(df))))) {
    df <- cbind(SampleID = rownames(df), df); rownames(df) <- NULL
  } else {
    names(df)[1] <- "SampleID"
  }
  df$SampleID <- norm_token(df$SampleID)
  df
}

# If samples are columns, transpose to rows=Samples
reorient_if_needed <- function(mat_df, known_ids){
  df <- as.data.frame(mat_df)
  # Try: are many column names “like” known ids after normalization?
  cn <- norm_token(colnames(df))
  overlap_cols <- sum(cn %in% known_ids)
  if (overlap_cols >= max(2, floor(nrow(df)*0.1))) {
    # samples in columns → transpose
    # use first column as a taxa label if it isn’t a sample
    # keep numeric only beyond the first id-like column set
    rownames(df) <- NULL
    colnames(df) <- cn
    # Find sample columns (in known_ids)
    sample_cols <- which(colnames(df) %in% known_ids)
    X <- as.matrix(df[, sample_cols, drop = FALSE])
    X <- t(X)
    out <- as.data.frame(X, check.names = FALSE)
    out$SampleID <- rownames(out); rownames(out) <- NULL
    out <- out[, c("SampleID", setdiff(names(out), "SampleID"))]
  } else {
    # assume samples already in rows
    out <- norm_ids_col(df)
  }
  # numeric coercion
  taxa_cols <- setdiff(names(out), "SampleID")
  for (c in taxa_cols) out[[c]] <- suppressWarnings(as.numeric(out[[c]]))
  out[is.na(out)] <- 0
  out
}

best_map_from_adist <- function(from, to){
  M <- adist(from, to, partial = TRUE, ignore.case = TRUE)
  rownames(M) <- from; colnames(M) <- to
  picks <- apply(M, 1, function(v) to[which.min(v)])
  distv <- apply(M, 1, min)
  data.frame(from = names(picks), to = unname(picks), edit_dist = as.integer(distv), row.names = NULL, stringsAsFactors = FALSE)
}

# ---------- Load ----------
PA_raw <- read_auto(args$presence)
HE_raw <- read_auto(args$hellinger)
ME_raw <- read_auto(args$metadata)

# normalize metadata, get known IDs
meta <- norm_ids_col(ME_raw)
if (!("ranch" %in% names(meta))) {
  cand <- names(meta)[grepl("ranch|site|location", names(meta), ignore.case = TRUE)][1]
  if (is.na(cand)) stop("Metadata needs a ranch/site/location column.")
  names(meta)[names(meta) == cand] <- "ranch"
}
meta$ranch <- norm_token(meta$ranch)
known_ids <- unique(meta$SampleID)

# Reorient matrices against known_ids
PA <- reorient_if_needed(PA_raw, known_ids)
HE <- reorient_if_needed(HE_raw, known_ids)

# First overlap check
sp <- intersect(PA$SampleID, known_ids)
sh <- intersect(HE$SampleID, known_ids)

# If either is zero, try a column-name normalization pass
if (length(sp) == 0 || length(sh) == 0) {
  # Re-normalize SampleID tokens in PA/HE a second time (safety)
  PA$SampleID <- norm_token(PA$SampleID)
  HE$SampleID <- norm_token(HE$SampleID)
  sp <- intersect(PA$SampleID, known_ids)
  sh <- intersect(HE$SampleID, known_ids)
}

used_map <- FALSE
map_df  <- NULL

# If still zero overlap, propose & apply a mapping via edit distance
if (length(sp) == 0 || length(sh) == 0) {
  mapP <- best_map_from_adist(PA$SampleID, known_ids)
  mapH <- best_map_from_adist(HE$SampleID, known_ids)
  # Accept only close matches (edit distance <= 3) to avoid bad mappings
  mapP_ok <- mapP %>% filter(edit_dist <= 3)
  mapH_ok <- mapH %>% filter(edit_dist <= 3)

  if (nrow(mapP_ok) > 0) {
    PA$SampleID <- ifelse(PA$SampleID %in% mapP_ok$from,
                          mapP_ok$to[match(PA$SampleID, mapP_ok$from)],
                          PA$SampleID)
  }
  if (nrow(mapH_ok) > 0) {
    HE$SampleID <- ifelse(HE$SampleID %in% mapH_ok$from,
                          mapH_ok$to[match(HE$SampleID, mapH_ok$from)],
                          HE$SampleID)
  }

  sp <- intersect(PA$SampleID, known_ids)
  sh <- intersect(HE$SampleID, known_ids)

  if (length(sp) > 0 || length(sh) > 0) {
    used_map <- TRUE
    map_df <- full_join(mapP_ok %>% mutate(source="presence"),
                        mapH_ok %>% mutate(source="hellinger"),
                        by = c("from","to","edit_dist"),
                        multiple = "all")
  }
}

# Final shared set (need intersection across ALL THREE)
shared <- Reduce(intersect, list(PA$SampleID, HE$SampleID, known_ids))

# Report + write
rep_txt <- file.path(args$outdir, "preflight_report.txt")
con <- file(rep_txt, "wt")
writeLines(sprintf("[Preflight] presence n=%d, hellinger n=%d, metadata n=%d", nrow(PA), nrow(HE), nrow(meta)), con)
writeLines(sprintf("[Preflight] overlap with metadata: presence=%d, hellinger=%d, shared(all three)=%d",
                   length(intersect(PA$SampleID, known_ids)),
                   length(intersect(HE$SampleID, known_ids)),
                   length(shared)), con)
if (length(shared) == 0) {
  # print first few non-overlaps to help eyeball
  missP <- setdiff(PA$SampleID, known_ids)[1:min(10, length(setdiff(PA$SampleID, known_ids)))]
  missH <- setdiff(HE$SampleID, known_ids)[1:min(10, length(setdiff(HE$SampleID, known_ids)))]
  writeLines("\n[Preflight] Examples of non-matching SampleID (presence):", con)
  writeLines(paste0("  - ", missP), con)
  writeLines("\n[Preflight] Examples of non-matching SampleID (hellinger):", con)
  writeLines(paste0("  - ", missH), con)
}
close(con)

if (used_map && !is.null(map_df) && nrow(map_df) > 0) {
  readr::write_csv(map_df, file.path(args$outdir, "suggested_id_map.csv"))
}

# Keep only shared samples, ordered the same way in all three
if (length(shared) == 0) stop("Still no overlapping SampleID after normalization/mapping.")

PA_out <- PA %>% filter(SampleID %in% shared) %>% arrange(SampleID)
HE_out <- HE %>% filter(SampleID %in% shared) %>% arrange(SampleID)
ME_out <- meta %>% filter(SampleID %in% shared) %>% arrange(SampleID)

readr::write_csv(PA_out, file.path(args$outdir, "presence_fixed.csv"))
readr::write_csv(HE_out, file.path(args$outdir, "hellinger_fixed.csv"))
readr::write_csv(ME_out, file.path(args$outdir, "metadata_fixed.csv"))

cat("[OK] Preflight fix complete. Wrote aligned files + report to:", args$outdir, "\n")