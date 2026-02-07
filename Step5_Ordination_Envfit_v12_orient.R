#!/usr/bin/env Rscript

# Step5_Ordination_Envfit_v12_orient.R
# NMDS ordination + envfit + PERMANOVA / PERMDISP for ranch communities.
# This version is simplified to assume:
#   - Hellinger matrix has taxa as rows and samples as columns
#     (first column = taxon name, remaining columns = SampleIDs)
#   - metadata_aligned.csv contains at least: SampleID, ranch, elev_m
# It automatically re-orients the Hellinger matrix so rows = samples
# and columns = taxa, then aligns on SampleID before running analyses.

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
})

ensure_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---------- ranch name normalization ----------

ranch_norm <- function(x, alias = list()) {
  if (length(alias)) {
    for (nm in names(alias)) {
      x[tolower(x) %in% tolower(alias[[nm]])] <- nm
    }
  }
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x <- ifelse(tolower(x) %in% c("el mezquital","mezquital","el mesquite","mesquite"),
              "El Mesquite", x)
  x <- ifelse(tolower(x) %in% c("casa viejas","casavieja","casa vieja"),
              "Casa Vieja", x)
  x <- ifelse(tolower(x) %in% c("la venada","lavenada"),
              "La Venada", x)
  x <- ifelse(tolower(x) %in% c("el tule","eltule"),
              "El Tule", x)
  x
}

# ---------- I/O helpers ----------

read_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv","txt","tsv")) {
    first <- readr::read_lines(path, n_max = 1)
    sep <- if (stringr::str_count(first, "\t") > stringr::str_count(first, ",")) "\t" else ","
    readr::read_delim(path, delim = sep, col_types = cols())
  } else if (ext %in% c("xls","xlsx")) {
    suppressPackageStartupMessages(library(readxl))
    readxl::read_excel(path)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

# Re-orient any taxa x samples table to samples x taxa with SampleID rownames.
# Assumes first column is a taxon / feature ID, remaining columns are samples.
reorient_to_samples_as_rows <- function(df) {
  df <- as.data.frame(df)

  # If already in sample-by-taxon form with a SampleID column, just use it.
  if ("SampleID" %in% names(df)) {
    rn <- df$SampleID
    taxa_cols <- setdiff(names(df), "SampleID")
    comm <- as.matrix(df[, taxa_cols, drop = FALSE])
    rownames(comm) <- rn
    return(comm)
  }

  if (ncol(df) < 2) {
    stop("Hellinger matrix must have at least 2 columns (taxon + samples).")
  }

  taxa_col  <- names(df)[1]
  taxa_ids  <- df[[taxa_col]]
  mat       <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- taxa_ids

  # transpose so rows = samples, cols = taxa
  comm <- t(mat)
  # sample IDs come from original column names (excluding taxon col)
  rownames(comm) <- colnames(df)[-1]
  comm
}

# Align community matrix and metadata on overlapping SampleID
align_two <- function(comm, meta) {
  if (!("SampleID" %in% names(meta))) {
    sid <- names(meta)[grepl("sampleid|sample_id|sample", names(meta), ignore.case = TRUE)]
    if (!length(sid)) stop("Metadata has no SampleID column.")
    meta <- meta |> dplyr::rename(SampleID = dplyr::all_of(sid[1]))
  }

  meta$SampleID   <- as.character(meta$SampleID)
  rownames(comm)  <- as.character(rownames(comm))

  common <- intersect(rownames(comm), meta$SampleID)
  if (!length(common)) {
    stop("No overlapping SampleID across hellinger and metadata.")
  }

  comm2 <- comm[common, , drop = FALSE]
  meta2 <- meta[match(common, meta$SampleID), , drop = FALSE]
  list(comm = comm2, meta = meta2)
}

# ---------- CLI ----------

option_list <- list(
  make_option(c("--hellinger"), type = "character",
              help = "Hellinger matrix (taxa x samples; first col = taxon, others = SampleIDs)."),
  make_option(c("--metadata"),  type = "character",
              help = "Aligned metadata (metadata_aligned.csv)."),
  make_option(c("--outdir"),    type = "character",
              help = "Output directory."),
  make_option(c("--ranch_order"), type = "character",
              default = "La Venada,El Mesquite,Casa Vieja,El Tule",
              help = "Comma-separated ranch order for plotting [default: %default]."),
  make_option(c("--alias_map"), type = "character", default = NULL,
              help = "Optional JSON mapping of ranch synonyms.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$hellinger) || is.null(opt$metadata) || is.null(opt$outdir)) {
  stop("You must provide --hellinger, --metadata, and --outdir.")
}

ensure_dir(opt$outdir)

# ---------- read inputs ----------

HE_raw <- read_auto(opt$hellinger)
META   <- read_auto(opt$metadata)

# clean metadata columns
if (!("SampleID" %in% names(META))) {
  sid <- names(META)[grepl("sampleid|sample_id|sample", names(META), ignore.case = TRUE)]
  if (!length(sid)) stop("Metadata must have a SampleID column or a recognizable equivalent.")
  META <- META |> dplyr::rename(SampleID = dplyr::all_of(sid[1]))
}
if (!("ranch" %in% names(META))) {
  rc <- names(META)[grepl("ranch", names(META), ignore.case = TRUE)]
  if (!length(rc)) stop("Metadata must contain a ranch column.")
  META <- META |> dplyr::rename(ranch = dplyr::all_of(rc[1]))
}
if (!("elev_m" %in% names(META))) {
  ec <- names(META)[grepl("elev|elevation", names(META), ignore.case = TRUE)]
  if (!length(ec)) stop("Metadata must contain an elevation column (elev_m/elevation).")
  META <- META |> dplyr::rename(elev_m = dplyr::all_of(ec[1]))
}

# normalize ranches and order factor
alias <- if (!is.null(opt$alias_map)) {
  jsonlite::fromJSON(opt$alias_map)
} else {
  list()
}
META$ranch <- ranch_norm(META$ranch, alias)

order_vec <- strsplit(opt$ranch_order, ",")[[1]] |> trimws()
META$ranch <- factor(META$ranch, levels = order_vec, ordered = TRUE)

# ---------- reorient hellinger and align ----------

HE_comm <- reorient_to_samples_as_rows(HE_raw)
pair    <- align_two(HE_comm, META)
HE      <- as.matrix(pair$comm)
META    <- pair$meta

# drop samples / taxa with all zero
zero_cols <- which(colSums(HE, na.rm = TRUE) == 0)
zero_rows <- which(rowSums(HE, na.rm = TRUE) == 0)
if (length(zero_cols)) {
  HE <- HE[, -zero_cols, drop = FALSE]
}
if (length(zero_rows)) {
  keep_ids <- rownames(HE)
  META <- META[META$SampleID %in% keep_ids, , drop = FALSE]
  HE   <- HE[keep_ids, , drop = FALSE]
}

if (nrow(HE) < 2) {
  stop("Not enough samples after filtering to run NMDS.")
}

set.seed(123)
nmds <- metaMDS(HE, distance = "bray", k = 2, trymax = 200,
                autotransform = FALSE, trace = FALSE)
S    <- scores(nmds, display = "sites")

# ---------- envfit for elevation ----------

env  <- envfit(nmds, META$elev_m, permutations = 999)
vec  <- scores(env, display = "vectors")

if (!is.null(vec) && nrow(vec) > 0) {
  v <- as.numeric(vec[1, 1:2])
  r2    <- env$r[1]^2
  pval  <- env$pvals[1]
} else {
  v    <- c(NA_real_, NA_real_)
  r2   <- NA_real_
  pval <- NA_real_
}

# ---------- orient so elevation arrow is sensible ----------

S_disp <- S
v_disp <- v

if (all(is.finite(v_disp))) {
  ang <- atan2(v_disp[2], v_disp[1])
  rot <- matrix(c(cos(-ang), -sin(-ang),
                  sin(-ang),  cos(-ang)),
                nrow = 2, byrow = TRUE)
  S_rot <- as.matrix(S) %*% t(rot)
  v_rot <- as.numeric(rot %*% matrix(v_disp, ncol = 1))

  elev_cor <- suppressWarnings(cor(S_rot[,1], META$elev_m, use = "complete.obs"))
  if (!is.na(elev_cor) && elev_cor > 0) {
    S_rot[,1] <- -S_rot[,1]
    v_rot[1]  <- -v_rot[1]
  }

  S_disp <- S_rot
  v_disp <- v_rot
}

# ---------- centroids ----------

df_scores <- data.frame(
  SampleID = rownames(S_disp),
  NMDS1    = S_disp[,1],
  NMDS2    = S_disp[,2]
) |>
  dplyr::left_join(META[, c("SampleID","ranch","elev_m")], by = "SampleID")

centroids <- df_scores |>
  dplyr::group_by(ranch) |>
  dplyr::summarise(
    x      = mean(NMDS1),
    y      = mean(NMDS2),
    elev_m = mean(elev_m, na.rm = TRUE),
    n      = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(ranch)

# ---------- PERMANOVA & PERMDISP ----------

dist_mat <- vegdist(HE, method = "bray")
adon     <- adonis2(dist_mat ~ ranch, data = META, permutations = 999)
bet      <- betadisper(dist_mat, META$ranch)
bet_perm <- permutest(bet, permutations = 999)

# ---------- save tables ----------

write_csv(df_scores, file.path(opt$outdir, "NMDS_scores_samples.csv"))
write_csv(centroids, file.path(opt$outdir, "NMDS_centroids_by_ranch.csv"))

# PERMANOVA
adon_tab <- data.frame(
  term = rownames(adon),
  F    = adon$F,
  R2   = adon$R2,
  p    = adon$`Pr(>F)`
)
write_csv(adon_tab, file.path(opt$outdir, "PERMANOVA_ranch.csv"))

# PERMDISP
bet_tab <- data.frame(
  ranch      = bet$group,          # use factor directly
  dispersion = bet$distances       # distance to group centroid
)

bet_p <- data.frame(
  F = bet_perm$tab[1, "F"],
  p = bet_perm$tab[1, "Pr(>F)"]
)

write_csv(bet_tab, file.path(opt$outdir, "PERMDISP_group_dispersion.csv"))
write_csv(bet_p,   file.path(opt$outdir, "PERMDISP_overall.csv"))

# ---------------------------------------------------
# envfit table (fully robust to missing values)
# ---------------------------------------------------
env_path <- file.path(opt$outdir, "envfit_elevation.csv")

valid_vector <- isTRUE(!is.null(v_disp)) &&
                isTRUE(length(v_disp) >= 2) &&
                isTRUE(is.finite(r2)) &&
                isTRUE(is.finite(pval))

if (valid_vector) {
  # Normal case
  env_tab <- data.frame(
    variable = "elev_m",
    r2       = r2,
    p        = pval,
    NMDS1    = v_disp[1],
    NMDS2    = v_disp[2]
  )
} else {
  # Fallback safe output
  env_tab <- data.frame(
    variable = "elev_m",
    r2       = NA_real_,
    p        = NA_real_,
    NMDS1    = NA_real_,
    NMDS2    = NA_real_
  )
  message("[WARN] envfit elevation vector missing; writing NA fallback row.")
}

write_csv(env_tab, env_path)

# ---------- plotting ----------

pal <- c("La Venada" = "#1f77b4",
         "El Mesquite" = "#2ca02c",
         "Casa Vieja"  = "#9467bd",
         "El Tule"     = "#d62728")

p <- ggplot(df_scores, aes(x = NMDS1, y = NMDS2, color = ranch)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_segment(
    data = df_scores,
    aes(
      x    = centroids$x[match(ranch, centroids$ranch)],
      y    = centroids$y[match(ranch, centroids$ranch)],
      xend = NMDS1,
      yend = NMDS2
    ),
    color       = "grey85",
    size        = 0.4,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = centroids,
    aes(x = x, y = y, fill = ranch),
    color       = "black",
    size        = 4,
    shape       = 21,
    stroke      = 0.7,
    inherit.aes = FALSE
  ) +
ggrepel::geom_text_repel(
  data = centroids,
  aes(x = x, y = y, label = ranch),
  size = 4.5,
  fontface = "bold",
  segment.color = "grey40",
  max.overlaps = Inf,
  min.segment.length = 0,
  fill = NA
 ) +
  scale_color_manual(values = pal, drop = FALSE) +
  scale_fill_manual(values = pal, drop = FALSE) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.title   = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = sprintf("NMDS (Bray–Curtis) | stress = %.3f", nmds$stress),
    color = "Ranch",
    fill  = "Ranch"
  )

if (all(is.finite(v_disp))) {
  arrow_scale <- 0.8 * max(diff(range(df_scores$NMDS1)),
                            diff(range(df_scores$NMDS2)))
  p <- p +
    geom_segment(
      aes(
        x    = 0,
        y    = 0,
        xend = v_disp[1] * arrow_scale,
        yend = v_disp[2] * arrow_scale
      ),
      arrow       = arrow(length = unit(0.25, "cm")),
      color       = "black",
      linewidth   = 0.7,
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x      = v_disp[1] * arrow_scale * 1.05,
      y      = v_disp[2] * arrow_scale * 1.05,
      label  = sprintf("elev_m (envfit)\nR² = %.2f, p = %.3f", r2, pval),
      hjust  = 0,
      vjust  = 0,
      size   = 3
    )
}

ggsave(
  file.path(opt$outdir, "NMDS_ordination_envfit_THESIS.png"),
  plot   = p,
  width  = 7.2,
  height = 5.2,
  dpi    = 400
)

ggsave(
  file.path(opt$outdir, "NMDS_ordination_envfit_THESIS.pdf"),
  plot   = p,
  width  = 7.2,
  height = 5.2,
  dpi    = 400
)

message("[OK] Wrote NMDS ordination, envfit, PERMANOVA, and PERMDISP outputs to: ", opt$outdir)