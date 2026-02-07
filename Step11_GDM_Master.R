#!/usr/bin/env Rscript
# Step11_GDM_Master.R  — dynamic CLI version
# GDM with permutation importance + hardened bootstrap CIs.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(gdm)
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a) > 0) a else b
cli_ok <- function(...) cat("[OK]",   sprintf(...), "\n")
cli_i  <- function(...) cat("[INFO]", sprintf(...), "\n")
cli_w  <- function(...) cat("[WARN]", sprintf(...), "\n")
cli_e  <- function(...) cat("[ERR]",  sprintf(...), "\n")

# =========================
# CLI (dynamic paths and settings)
# =========================
opt_list <- list(
  make_option(c("--comm"),   type = "character", help = "Community CSV (wide matrix: samples x taxa, or taxa x samples)."),
  make_option(c("--meta"),   type = "character", help = "Metadata CSV (ID + lon/lat + env cols)."),
  make_option(c("--outdir"), type = "character", help = "Output directory."),
  make_option(c("--id-col"),  type = "character", default = "sample_id"),
  make_option(c("--lon-col"), type = "character", default = "lon"),
  make_option(c("--lat-col"), type = "character", default = "lat"),
  make_option(c("--env-cols"), type = "character",
              default = "ElevDEM,slope,roughness,topo_res,EVI_amp,EVI_area,EVI_greenup,TmeanK,TmaxK"),
  make_option(c("--boot"),    type = "integer",  default = 500),
  make_option(c("--splines"), type = "integer",  default = 3),
  make_option(c("--nperm"),   type = "integer",  default = 250),
  make_option(c("--seed"),    type = "integer",  default = 42)
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Nogales defaults if flags omitted
COMM_CSV <- opt$comm %||% "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/USE_Analysis_Pseudo_4_ranches/Step11_Sensitivity_FINAL/MASTER_transposed_for_within.csv"
META_CSV <- opt$meta %||% "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/Step11a_HIRES/samples_env_HIRES.csv"
OUTDIR   <- opt$outdir %||% "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/USE_Analysis_Pseudo_4_ranches/Step11_Sensitivity_FINAL"

ID_COL   <- opt$`id-col`
LON_COL  <- opt$`lon-col`
LAT_COL  <- opt$`lat-col`
ENV_COLS <- strsplit(opt$`env-cols`, ",", fixed = TRUE)[[1]] |> trimws()
BOOT     <- opt$boot
SPLINES  <- opt$splines
N_PERM   <- opt$nperm
SEED     <- opt$seed

set.seed(SEED)
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

cli_i("Inputs:")
cli_i("  comm   : %s", COMM_CSV)
cli_i("  meta   : %s", META_CSV)
cli_i("  outdir : %s", OUTDIR)
cli_i("Cols : id=%s, lon=%s, lat=%s", ID_COL, LON_COL, LAT_COL)
cli_i("ENV  : %s", paste(ENV_COLS, collapse = ", "))
cli_i("Params: boot=%d, splines=%d, nPerm=%d, seed=%d", BOOT, SPLINES, N_PERM, SEED)

# =========================
# HELPERS
# =========================
to_numeric <- function(x) suppressWarnings(as.numeric(x))

# -------------------------
# NEW: helper to auto orient the community matrix
# -------------------------
auto_orient_comm <- function(comm_raw, id_col) {
  # If ID column already present, assume rows = samples and return as is
  if (id_col %in% names(comm_raw)) {
    cli_i("Community table already has ID column '%s'; assuming rows = samples.", id_col)
    return(comm_raw)
  }

  # Try a couple of common alternatives for the ID column name
  alt_ids <- c("SampleID", "sample", "site_id", "SiteID")
  alt_hit <- alt_ids[alt_ids %in% names(comm_raw)]
  if (length(alt_hit) > 0) {
    old <- alt_hit[1]
    cli_w("ID column '%s' not found. Using '%s' as ID and renaming to '%s'.",
          id_col, old, id_col)
    comm2 <- comm_raw
    names(comm2)[names(comm2) == old] <- id_col
    return(comm2)
  }

  # If we reach here, there is no obvious ID column.
  # We now assume the format is: taxa in rows, samples in columns, first column = taxon name.
  cli_w("ID column '%s' not found. Trying taxa x samples format (first column = taxon, remaining columns = samples).", id_col)

  if (ncol(comm_raw) < 2) {
    stop("Cannot auto orient community matrix: fewer than 2 columns.")
  }

  taxon_col <- names(comm_raw)[1]
  cli_i("Assuming first column '%s' is taxon name; all other columns are samples.", taxon_col)

  # Extract taxa and sample matrix
  taxa <- as.character(comm_raw[[taxon_col]])
  mat  <- as.data.frame(comm_raw[, -1, drop = FALSE], stringsAsFactors = FALSE)

  # Make sure columns represent samples
  sample_ids <- names(mat)
  cli_i("Detected %d samples from column names.", length(sample_ids))

  # Coerce to numeric
  mat_num <- as.matrix(apply(mat, 2, to_numeric))
  rownames(mat_num) <- taxa

  # Transpose so rows = samples, columns = taxa
  mat_t <- t(mat_num)
  mat_t <- as.data.frame(mat_t, stringsAsFactors = FALSE)
  mat_t <- tibble::rownames_to_column(mat_t, var = id_col)

  cli_ok("Reoriented community matrix: now %d samples x %d taxa.", nrow(mat_t), ncol(mat_t) - 1)
  mat_t
}

# =========================
# LOAD DATA
# =========================
cli_i("Reading community table: %s", COMM_CSV)
comm_raw <- suppressMessages(read_csv(COMM_CSV, show_col_types = FALSE))

# Use new helper to orient the community table to rows = samples, columns = taxa
comm <- auto_orient_comm(comm_raw, ID_COL)

if (!ID_COL %in% names(comm)) {
  stop(sprintf("After orientation, ID column '%s' is still missing from community table.", ID_COL))
}

cli_i("Reading metadata table: %s", META_CSV)
meta <- suppressMessages(read_csv(META_CSV, show_col_types = FALSE))

needed_meta <- c(ID_COL, LON_COL, LAT_COL)
if (!all(needed_meta %in% names(meta))) {
  stop(sprintf("Metadata must contain: %s", paste(needed_meta, collapse = ", ")))
}

missing_env <- setdiff(ENV_COLS, names(meta))
if (length(missing_env) > 0) {
  cli_w("Missing env columns in metadata: %s", paste(missing_env, collapse = ", "))
  ENV_COLS <- setdiff(ENV_COLS, missing_env)
}
if (length(ENV_COLS) == 0) stop("No valid ENV_COLS found in metadata.")

# =========================
# SANITY and MERGE
# =========================
# Clean spaces from IDs so joins behave
comm[[ID_COL]] <- gsub("\\s+", "_", comm[[ID_COL]])
meta[[ID_COL]] <- gsub("\\s+", "_", meta[[ID_COL]])

n_comm <- length(unique(comm[[ID_COL]]))
n_meta <- length(unique(meta[[ID_COL]]))
n_olap <- sum(unique(comm[[ID_COL]]) %in% unique(meta[[ID_COL]]))
cli_i("Samples: comm=%d, meta=%d, overlap=%d", n_comm, n_meta, n_olap)
if (n_olap < 4) stop("Insufficient sample ID overlap (< 4).")

# Identify taxa columns
taxa_cols <- setdiff(names(comm), c(ID_COL, "ranch", "trap"))
comm[taxa_cols] <- lapply(comm[taxa_cols], function(x) to_numeric(replace(x, is.na(x), 0)))

# Drop all zero taxa
nonzero_taxa <- vapply(comm[taxa_cols], function(x) sum(x, na.rm = TRUE) > 0, logical(1))
dropped <- setdiff(taxa_cols, taxa_cols[nonzero_taxa])
if (length(dropped)) {
  cli_w("Dropping %d all-zero taxa (e.g., %s)", length(dropped), paste(head(dropped, 5), collapse = ", "))
}
taxa_cols <- taxa_cols[nonzero_taxa]

# Merge with metadata
DF <- comm %>%
  dplyr::select(all_of(c(ID_COL, taxa_cols))) %>%
  dplyr::left_join(meta %>% dplyr::select(all_of(c(ID_COL, LON_COL, LAT_COL, ENV_COLS))), by = ID_COL) %>%
  dplyr::filter(!is.na(.data[[LON_COL]]), !is.na(.data[[LAT_COL]]))

if (nrow(DF) < 4) stop("Need >= 4 samples with coordinates after merge.")
for (v in ENV_COLS) DF[[v]] <- to_numeric(DF[[v]])

# =========================
# GDM INPUT (bioFormat = 1)
# =========================
sp_mat_core <- DF %>%
  dplyr::select(all_of(taxa_cols)) %>%
  as.data.frame(stringsAsFactors = FALSE)

sp_mat <- data.frame(
  site = DF[[ID_COL]],
  sp_mat_core,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

site_df <- DF %>%
  dplyr::select(all_of(c(ID_COL, LON_COL, LAT_COL, ENV_COLS))) %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(site_df) <- c("site", "lon", "lat", ENV_COLS)
site_df$site <- as.character(site_df$site)
site_df$lon  <- suppressWarnings(as.numeric(site_df$lon))
site_df$lat  <- suppressWarnings(as.numeric(site_df$lat))
for (v in ENV_COLS) site_df[[v]] <- suppressWarnings(as.numeric(site_df[[v]]))

sp_table <- gdm::formatsitepair(
  bioData    = as.data.frame(sp_mat, stringsAsFactors = FALSE),
  bioFormat  = 1,
  siteColumn = "site",
  XColumn    = "lon",
  YColumn    = "lat",
  predData   = as.data.frame(site_df, stringsAsFactors = FALSE)
)

# =========================
# FIT GDM and VARIMP
# =========================
n_env  <- length(ENV_COLS)
n_pred <- n_env + 1
splines_vec <- rep(SPLINES, n_pred)

cli_i("Fitting GDM (geo=TRUE, splines=%s) with %d predictors …",
      paste(splines_vec, collapse = ","), n_pred)

gdm_mod <- gdm(sp_table, geo = TRUE, splines = splines_vec)

saveRDS(gdm_mod, file.path(OUTDIR, "GDM_model.rds"))
write_csv(tibble(deviance_explained = gdm_mod$explained),
          file.path(OUTDIR, "GDM_deviance_explained.csv"))
cli_ok("Deviance explained: %.3f", gdm_mod$explained)

# --- Permutation importance ---
cli_i("Permutation importance (nPerm=%d) …", N_PERM)
varImp <- gdm.varImp(sp_table, geo = TRUE, nPerm = N_PERM, parallel = FALSE)
readr::write_csv(tibble::as_tibble(varImp$summary),
                 file.path(OUTDIR, "GDM_varImp_permutation.csv"))

# =========================
# BOOTSTRAP CIs
# =========================
extract_ispline_heights <- function(mod, envs) {
  spl <- isplineExtract(mod)
  out <- sapply(envs, function(v) {
    if (v %in% colnames(spl$y)) {
      max(spl$y[, v], na.rm = TRUE)
    } else {
      NA_real_
    }
  })
  as.numeric(out)
}

MIN_UNIQUE_SITES <- 4
MIN_TAXA         <- 2

cli_i("Bootstrapping (B=%d) by resampling sites with replacement …", BOOT)
site_ids <- unique(DF[[ID_COL]])
boot_heights <- matrix(NA_real_, nrow = BOOT, ncol = length(ENV_COLS))
colnames(boot_heights) <- ENV_COLS
boot_dev <- rep(NA_real_, BOOT)

for (b in seq_len(BOOT)) {
  samp <- sample(site_ids, length(site_ids), replace = TRUE)
  if (length(unique(samp)) < MIN_UNIQUE_SITES) next
  B <- DF[DF[[ID_COL]] %in% samp, , drop = FALSE]

  sp_b <- B %>%
    dplyr::select(all_of(taxa_cols)) %>%
    as.data.frame()
  sp_b[] <- lapply(sp_b, function(x) {
    x[is.na(x)] <- 0
    x
  })
  nz_b <- vapply(sp_b, function(col) sum(as.numeric(col), na.rm = TRUE) > 0, logical(1))
  if (sum(nz_b) < MIN_TAXA) next
  sp_b <- sp_b[, nz_b, drop = FALSE]

  site_b <- B %>%
    dplyr::select(all_of(c(ID_COL, LON_COL, LAT_COL, ENV_COLS))) %>%
    as.data.frame(stringsAsFactors = FALSE)
  colnames(site_b) <- c("site", "lon", "lat", ENV_COLS)
  site_b$site <- as.character(site_b$site)
  site_b$lon  <- suppressWarnings(as.numeric(site_b$lon))
  site_b$lat  <- suppressWarnings(as.numeric(site_b$lat))
  for (v in ENV_COLS) site_b[[v]] <- suppressWarnings(as.numeric(site_b[[v]]))

  bio_b <- data.frame(site = site_b$site, sp_b, check.names = FALSE)
  sp_tab_b <- try(
    gdm::formatsitepair(
      bioData    = bio_b,
      bioFormat  = 1,
      siteColumn = "site",
      XColumn    = "lon",
      YColumn    = "lat",
      predData   = site_b
    ),
    silent = TRUE
  )
  if (inherits(sp_tab_b, "try-error")) next
  if (nrow(sp_tab_b) < 1) next

  mod_b <- try(
    gdm(sp_tab_b, geo = TRUE, splines = rep(SPLINES, length(ENV_COLS) + 1)),
    silent = TRUE
  )
  if (inherits(mod_b, "try-error")) next
  if (is.null(mod_b) || is.null(mod_b$explained)) next

  boot_dev[b] <- mod_b$explained
  ih <- try(extract_ispline_heights(mod_b, ENV_COLS), silent = TRUE)
  if (!inherits(ih, "try-error")) boot_heights[b, ] <- ih

  if (b %% max(1, floor(BOOT / 10)) == 0) {
    cli_i("  … %d/%d", b, BOOT)
  }
}

bh <- as.data.frame(boot_heights)
bh$iter <- seq_len(BOOT)
readr::write_csv(bh, file.path(OUTDIR, "GDM_bootstrap_ispline_heights.csv"))
readr::write_csv(tibble(iter = seq_len(BOOT), deviance_explained = boot_dev),
                 file.path(OUTDIR, "GDM_bootstrap_deviance.csv"))

valid_cols <- vapply(as.data.frame(boot_heights), function(col) any(!is.na(col)), logical(1))
if (!any(valid_cols)) {
  cli_w("All predictors lack valid bootstrap values; skipping CI summary and importance plot.")
} else {
  env_for_ci <- ENV_COLS[valid_cols]
  summ <- lapply(env_for_ci, function(v) {
    x <- boot_heights[, v]
    x <- x[!is.na(x)]
    tibble(
      predictor = v,
      mean      = mean(x),
      sd        = sd(x),
      q025      = quantile(x, 0.025),
      q500      = quantile(x, 0.5),
      q975      = quantile(x, 0.975)
    )
  }) %>%
    bind_rows()
  readr::write_csv(summ, file.path(OUTDIR, "GDM_ispline_height_CI.csv"))
}

# =========================
# PLOTS
# =========================
spl <- isplineExtract(gdm_mod)
pdf(file.path(OUTDIR, "Fig_GDM_Isplines.pdf"), width = 8.0, height = 6.5)
par(mfrow = c(ceiling(length(ENV_COLS) / 2), 2), mar = c(4, 4, 2, 1))
for (v in ENV_COLS) {
  if (!v %in% colnames(spl$y)) next
  plot(
    spl$x[, v], spl$y[, v],
    type = "l",
    xlab = v,
    ylab = "I-spline partial ecological distance",
    main = v
  )
  grid()
}
dev.off()

ci_path <- file.path(OUTDIR, "GDM_ispline_height_CI.csv")
if (file.exists(ci_path)) {
  ci <- readr::read_csv(ci_path, show_col_types = FALSE)
  ci$predictor <- factor(ci$predictor, levels = ci$predictor[order(ci$mean, decreasing = TRUE)])
  p <- ggplot(ci, aes(x = predictor, y = mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2) +
    labs(
      x = NULL,
      y = "I-spline height (importance)",
      title = "GDM predictor importance (bootstrap)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUTDIR, "Fig_GDM_Importance_Bootstrap.pdf"), p, width = 7.5, height = 5.0)
} else {
  cli_w("CI file not found at %s; skipping importance plot.", ci_path)
}

pred_ecodist <- tryCatch(gdm.predict(gdm_mod, sp_table), error = function(e) NULL)
if (!is.null(pred_ecodist)) {
  dd <- tibble(geoDist = sp_table$dist, ecoDist = pred_ecodist)
  readr::write_csv(dd, file.path(OUTDIR, "GDM_distance_decay.csv"))
  pdf(file.path(OUTDIR, "Fig_GDM_distance_decay.pdf"), width = 6.5, height = 5.0)
  plot(
    dd$geoDist, dd$ecoDist,
    pch = 19, cex = 0.6,
    xlab = "Geographic distance (km)",
    ylab = "Predicted ecological distance",
    main = "GDM distance–decay"
  )
  grid()
  dev.off()
}

# =========================
# SUMMARY and SESSION INFO
# =========================
summary_lines <- c(
  sprintf("Deviance explained: %.3f", gdm_mod$explained),
  sprintf("Predictors used: %s", paste(ENV_COLS, collapse = ", ")),
  sprintf("Bootstrap iterations: %d", BOOT),
  sprintf("Permutation nPerm: %d", N_PERM)
)
writeLines(summary_lines, con = file.path(OUTDIR, "GDM_summary.txt"))
writeLines(c(capture.output(sessionInfo())), con = file.path(OUTDIR, "sessionInfo.txt"))

cli_ok("All outputs written to: %s", OUTDIR)