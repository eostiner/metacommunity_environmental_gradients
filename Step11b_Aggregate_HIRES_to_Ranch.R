#!/usr/bin/env Rscript
# Step11b_Aggregate_HIRES_to_Ranch.R
#
# WHAT THIS DOES
#   Takes the hi-res environmental table at the trap level
#   (samples_env_HIRES.csv from Step11a) and collapses it to
#   one row per ranch, using unweighted means.
#
#   This produces ranch_env_HIRES.csv, which is the metadata
#   used by Step11_GDM_Master.R alongside the ranch-level
#   community matrix (taxa_matrix_by_ranch_*_TRANS.csv).
#
# INPUTS (hard-coded paths for Nogales METACOMMUNITY):
#   - samples_env_HIRES.csv  (trap-level)
#   - samples.csv            (to get ranch IDs for each trap)
#
# OUTPUT:
#   - ranch_env_HIRES.csv    (one row per ranch)
#
# RATIONALE:
#   - Traps are equal-effort spatial replicates within each ranch.
#     For GDM at ranch scale, we want the typical environmental
#     context of each ranch, so we use unweighted means across
#     traps per ranch for ElevDEM, slope, roughness, topo_res,
#     EVI metrics, and temperature.
#
#   - The resulting file has:
#       sample_id, lon, lat,
#       ElevDEM, slope, roughness, topo_res,
#       EVI_amp, EVI_area, EVI_greenup,
#       TmeanK, TmaxK
#     with sample_id matching the ranch IDs used in the
#     ranch-level community matrix.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -------------------------------------------------------------------
# 1) Define paths (Nogales-specific, but can be parameterized later)
# -------------------------------------------------------------------
hires_csv <- "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/Step11a_HIRES_RANCH/samples_env_HIRES.csv"

samples_csv <- "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/samples.csv"

out_csv <- "/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY/Step11a_HIRES_RANCH/ranch_env_HIRES.csv"

message("[INFO] HIRES env (trap level): ", hires_csv)
message("[INFO] samples with ranch IDs: ", samples_csv)
message("[INFO] Output (ranch level): ", out_csv)

if (!file.exists(hires_csv)) {
  stop("HIRES CSV not found at: ", hires_csv)
}
if (!file.exists(samples_csv)) {
  stop("samples.csv not found at: ", samples_csv)
}

# -------------------------------------------------------------------
# 2) Read trap-level HIRES env and original samples table
# -------------------------------------------------------------------
hires <- read_csv(hires_csv, show_col_types = FALSE)
samples <- read_csv(samples_csv, show_col_types = FALSE)

if (!("sample_id" %in% names(hires))) {
  stop("samples_env_HIRES.csv must contain 'sample_id' column.")
}
if (!all(c("sample_id", "ranch") %in% names(samples))) {
  stop("samples.csv must contain 'sample_id' and 'ranch' columns.")
}

# -------------------------------------------------------------------
# 3) Attach ranch IDs to each trap row
# -------------------------------------------------------------------
hires2 <- hires %>%
  left_join(samples %>% select(sample_id, ranch), by = "sample_id")

if (any(is.na(hires2$ranch))) {
  warning("Some HIRES rows have no matching ranch in samples.csv.")
}

# -------------------------------------------------------------------
# 4) Collapse to ranch: unweighted mean across traps per ranch
# -------------------------------------------------------------------
ranch_env <- hires2 %>%
  group_by(ranch) %>%
  summarise(
    lon         = mean(lon,        na.rm = TRUE),
    lat         = mean(lat,        na.rm = TRUE),
    ElevDEM     = mean(ElevDEM,    na.rm = TRUE),
    slope       = mean(slope,      na.rm = TRUE),
    roughness   = mean(roughness,  na.rm = TRUE),
    topo_res    = mean(topo_res,   na.rm = TRUE),
    EVI_amp     = mean(EVI_amp,    na.rm = TRUE),
    EVI_area    = mean(EVI_area,   na.rm = TRUE),
    EVI_greenup = mean(EVI_greenup,na.rm = TRUE),
    TmeanK      = mean(TmeanK,     na.rm = TRUE),
    TmaxK       = mean(TmaxK,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # GDM expects sample_id column to match the community matrix
  rename(sample_id = ranch)

# -------------------------------------------------------------------
# 5) Write output and print a small preview
# -------------------------------------------------------------------
write_csv(ranch_env, out_csv)
cat("[OK] Wrote ranch-level HIRES env to:\n   ", out_csv, "\n", sep = "")
print(ranch_env)