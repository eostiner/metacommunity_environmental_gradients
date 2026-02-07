#!/usr/bin/env Rscript
# Step11a_Prepare_HiRes_Rasters_and_Metadata.R  (DYNAMIC / CLI VERSION)
# ------------------------------------------------------------------------------
# WHAT THIS DOES
#   1) Reads high-resolution rasters (DEM required; EVI optional; temps coarse OK).
#   2) Derives slope, roughness, topo_res from DEM at target resolution.
#   3) Aligns everything to a common projected CRS + grid resolution.
#   4) Extracts values at your samples and writes a new samples_env_HIRES.csv.
#   5) Saves per-layer GeoTIFFs + one stacked raster.
#
# HOW TO RUN (examples; flags explained below)
#   Rscript Step11a_Prepare_HiRes_Rasters_and_Metadata.R \
#     --samples "/path/to/samples_env.csv" \
#     --outdir  "/path/to/HIRES_Prepped" \
#     --dem     "/path/to/DEM_30m.tif" \
#     --tmean   "/path/to/TmeanK.tif" \
#     --tmax    "/path/to/TmaxK.tif" \
#     --evi-amp "/path/to/EVI_amp_10m.tif" \
#     --evi-area "/path/to/EVI_area_10m.tif" \
#     --evi-greenup "/path/to/EVI_greenup_10m.tif" \
#     --crs "EPSG:32612" \
#     --res 30 \
#     --roughwin 5 \
#     --id-col sample_id --lon-col lon --lat-col lat
#
# MINIMUM RUN (if you don’t have EVI hires yet; falls back to temps + DEM-derived preds)
#   Rscript Step11a_Prepare_HiRes_Rasters_and_Metadata.R \
#     --samples "/path/to/samples_env.csv" \
#     --outdir  "/path/to/HIRES_Prepped" \
#     --dem     "/path/to/DEM_30m.tif" \
#     --tmean   "/path/to/TmeanK.tif" \
#     --tmax    "/path/to/TmaxK.tif"
#
# NOTES
#   • If --crs is omitted, we AUTO-DETECT a UTM zone from your sample lon/lat.
#   • If --res is omitted, we use the native DEM resolution (rounded).
#   • Outputs are compatible with Step11_GDM_Master.R. Point META_CSV to the new CSV.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse) # <-- CLI parser
  library(terra)    # raster ops
  library(sf)       # vector ops
  library(dplyr)
  library(readr)
})

# ---------------------------
# CLI: define flags / options
# ---------------------------
opt_list <- list(
  make_option(c("--samples"), type = "character", help = "Path to samples_env CSV (must include ID/LON/LAT columns)."),
  make_option(c("--outdir"),  type = "character", help = "Output directory for HIRES results."),
  make_option(c("--dem"),     type = "character", help = "High-resolution DEM GeoTIFF (required)."),

  make_option(c("--tmean"),   type = "character", help = "Mean temperature raster (coarse ok).", default = NA_character_),
  make_option(c("--tmax"),    type = "character", help = "Max temperature raster (coarse ok).",  default = NA_character_),

  make_option(c("--evi-amp"),     type = "character", help = "EVI amplitude raster (optional).",  default = NA_character_),
  make_option(c("--evi-area"),    type = "character", help = "EVI area raster (optional).",       default = NA_character_),
  make_option(c("--evi-greenup"), type = "character", help = "EVI greenup raster (optional).",    default = NA_character_),

  make_option(c("--crs"),     type = "character", help = "Target projected CRS (e.g., 'EPSG:32612'). Auto-detect if missing.", default = NA_character_),
  make_option(c("--res"),     type = "double",    help = "Target grid resolution in meters. Default = DEM native.", default = NA_real_),
  make_option(c("--roughwin"),type = "integer",   help = "Roughness window size (pixels). Default=5.", default = 5L),

  make_option(c("--id-col"),  type = "character", help = "ID column name in samples CSV.", default = "sample_id"),
  make_option(c("--lon-col"), type = "character", help = "Longitude column name.",        default = "lon"),
  make_option(c("--lat-col"), type = "character", help = "Latitude column name.",         default = "lat"),

  make_option(c("--env-cols"), type = "character",
              help = "Comma-sep env names for GDM order. Default='ElevDEM,slope,roughness,topo_res,EVI_amp,EVI_area,EVI_greenup,TmeanK,TmaxK'",
              default = "ElevDEM,slope,roughness,topo_res,EVI_amp,EVI_area,EVI_greenup,TmeanK,TmaxK")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# ---------------------------
# Input validation & printing
# ---------------------------
stopifnot(!is.null(opt$samples), !is.null(opt$outdir), !is.null(opt$dem))
if (!file.exists(opt$samples)) stop("--samples not found: ", opt$samples)
if (!file.exists(opt$dem))     stop("--dem not found: ", opt$dem)
if (!is.na(opt$tmean) && !file.exists(opt$tmean)) stop("--tmean not found: ", opt$tmean)
if (!is.na(opt$tmax)  && !file.exists(opt$tmax))  stop("--tmax not found: ",  opt$tmax)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

ENV_COLS <- strsplit(opt$`env-cols`, ",", fixed = TRUE)[[1]] |> trimws()

message("[INFO] inputs:")
message("  samples   : ", opt$samples)
message("  outdir    : ", opt$outdir)
message("  dem       : ", opt$dem)
message("  tmean     : ", ifelse(is.na(opt$tmean), "(none)", opt$tmean))
message("  tmax      : ", ifelse(is.na(opt$tmax),  "(none)", opt$tmax))
message("  evi_amp   : ", ifelse(is.na(opt$`evi-amp`), "(none)", opt$`evi-amp`))
message("  evi_area  : ", ifelse(is.na(opt$`evi-area`), "(none)", opt$`evi-area`))
message("  evi_green : ", ifelse(is.na(opt$`evi-greenup`), "(none)", opt$`evi-greenup`))
message("  CRS       : ", ifelse(is.na(opt$crs), "(auto-detect UTM)", opt$crs))
message("  res (m)   : ", ifelse(is.na(opt$res), "(DEM native)", opt$res))
message("  roughwin  : ", opt$roughwin)
message("  env_cols  : ", paste(ENV_COLS, collapse = ", "))

ID_COL  <- opt$`id-col`; LON_COL <- opt$`lon-col`; LAT_COL <- opt$`lat-col`

# ---------------------------
# Helpers
# ---------------------------
project_to <- function(r, crs) {
  # project only if needed
  if (is.na(terra::crs(r)) || terra::crs(r) != crs) terra::project(r, crs) else r
}

auto_utm <- function(lon, lat) {
  # crude UTM zone chooser; works fine for small AOIs
  zone <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  hemi <- if (mean(lat, na.rm = TRUE) >= 0) "326" else "327"  # 326xx = north, 327xx = south
  sprintf("EPSG:%s%02d", hemi, zone)
}

# ---------------------------
# Read samples & decide CRS/res
# ---------------------------
samples <- readr::read_csv(opt$samples, show_col_types = FALSE)
needed <- c(ID_COL, LON_COL, LAT_COL)
if (!all(needed %in% names(samples))) stop("samples CSV missing columns: ", paste(needed, collapse = ", "))

# Set CRS
target_crs <- if (is.na(opt$crs)) {
  utm <- auto_utm(samples[[LON_COL]], samples[[LAT_COL]])
  message("[INFO] auto-detected UTM CRS: ", utm)
  utm
} else opt$crs

# ---------------------------
# Read DEM and determine target resolution  ✅ (terra-safe)
# ---------------------------
dem0 <- rast(opt$dem)
dem0 <- project_to(dem0, target_crs)

# Native DEM resolution (in projected meters); pick target_res accordingly
dem_res_native <- mean(res(dem0))
target_res <- if (is.na(opt$res) || opt$res <= 0) round(dem_res_native) else opt$res
message("[INFO] target grid resolution (m): ", target_res)

# ---------------------------
# Resample DEM to target grid  ✅ (use terra::ext(), not raster::extent())
# ---------------------------
# Build a template grid using the DEM’s extent, desired resolution, and target CRS
tmpl <- rast(ext = ext(dem0), res = target_res, crs = target_crs)

# Align the DEM to the template (resample with bilinear)
dem  <- resample(dem0, tmpl, method = "bilinear")
names(dem) <- "ElevDEM"

# Samples -> sf in WGS84 -> target CRS
pts_sf <- st_as_sf(samples, coords = c(LON_COL, LAT_COL), crs = 4326) |> st_transform(target_crs)

# ---------------------------
# DEM-derived predictors
# ---------------------------
ElevDEM <- dem
slope <- terrain(ElevDEM, v = "slope", unit = "degrees"); names(slope) <- "slope"
roughness <- terrain(ElevDEM, v = "roughness");           names(roughness) <- "roughness"

# local relief as topo_res
w <- matrix(1, nrow = opt$roughwin, ncol = opt$roughwin)
local_max <- focal(ElevDEM, w, fun = max, na.policy = "omit")
local_min <- focal(ElevDEM, w, fun = min, na.policy = "omit")
topo_res  <- local_max - local_min; names(topo_res) <- "topo_res"

# ---------------------------
# Optional EVI rasters
# ---------------------------
# get_and_align(): project to target_crs and resample to DEM template (ElevDEM).
get_and_align <- function(path, template, crs) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    message("[WARN] Raster not found (skipping): ", path)
    return(NULL)
  }
  rast(path) |> project_to(crs) |> resample(template, method = "bilinear")
}

# Try to read each EVI layer; if any is missing, it stays NULL (safe).
EVI_amp     <- get_and_align(opt$`evi-amp`,     ElevDEM, target_crs)
if (!is.null(EVI_amp))     names(EVI_amp)     <- "EVI_amp"

EVI_area    <- get_and_align(opt$`evi-area`,    ElevDEM, target_crs)
if (!is.null(EVI_area))    names(EVI_area)    <- "EVI_area"

EVI_greenup <- get_and_align(opt$`evi-greenup`, ElevDEM, target_crs)
if (!is.null(EVI_greenup)) names(EVI_greenup) <- "EVI_greenup"

# ---------------------------
# Temperatures (optional but recommended)
# ---------------------------
TmeanK <- if (!is.na(opt$tmean)) { get_and_align(opt$tmean, ElevDEM, target_crs) } else NULL
TmaxK  <- if (!is.na(opt$tmax))  { get_and_align(opt$tmax,  ElevDEM, target_crs) } else NULL
if (!is.null(TmeanK)) names(TmeanK) <- "TmeanK"
if (!is.null(TmaxK))  names(TmaxK)  <- "TmaxK"

# ---------------------------
# STACK & SAVE RASTERS
# ---------------------------
layers <- list(
  ElevDEM,
  slope,
  roughness,
  topo_res,
  EVI_amp,
  EVI_area,
  EVI_greenup,
  TmeanK,
  TmaxK
)

# Keep only non-NULL layers
avail <- Filter(Negate(is.null), layers)
if (length(avail) == 0) {
  stop("No raster layers available to stack (at least DEM-derived layers should exist).")
}

r_stack <- rast(avail)

# Reorder to the requested ENV_COLS if present
names_req <- ENV_COLS
keep <- intersect(names_req, names(r_stack))  # <- define `keep` here
if (length(keep) == 0) {
  # fall back to whatever we actually have
  keep <- names(r_stack)
  message("[WARN] None of ENV_COLS matched the stack; keeping all available layers: ",
          paste(keep, collapse = ", "))
}
r_stack <- r_stack[[keep]]

# Write each layer + one stack
for (nm in names(r_stack)) {
  outp <- file.path(opt$outdir, sprintf("%s_%sm.tif", nm, target_res))
  writeRaster(r_stack[[nm]], outp, overwrite = TRUE)
}
writeRaster(r_stack, file.path(opt$outdir, sprintf("ENV_HIRES_STACK_%sm.tif", target_res)),
            overwrite = TRUE)
            
            # ---------------------------
# EXTRACT TO POINTS -> NEW samples_env_HIRES.csv
# ---------------------------
# safety: make sure `keep` exists even if code was reorganized
if (!exists("keep")) keep <- names(r_stack)

# 1) extract raster values at sample locations (in target CRS)
# terra::extract returns a data.frame with an ID column first
vals <- terra::extract(r_stack, vect(pts_sf))
if (is.null(dim(vals))) {
  # single-layer extract can return a vector; coerce to data.frame
  vals <- data.frame(vals)
} else {
  vals <- as.data.frame(vals[, -1, drop = FALSE])  # drop the first ID column
}
colnames(vals) <- keep

# 2) put lon/lat back in WGS84 so your CSV matches other files
pts_ll <- sf::st_transform(pts_sf, 4326)
coords_ll <- sf::st_coordinates(pts_ll)

# 3) assemble output table: ID + lon/lat + extracted predictors (only those we kept)
out_df <- samples %>%
  dplyr::select(all_of(ID_COL)) %>%
  dplyr::mutate(
    !!LON_COL := coords_ll[, 1],
    !!LAT_COL := coords_ll[, 2]
  ) %>%
  dplyr::bind_cols(vals) %>%
  dplyr::select(all_of(c(ID_COL, LON_COL, LAT_COL, keep)))

# 4) write it
out_csv <- file.path(opt$outdir, "samples_env_HIRES.csv")
readr::write_csv(out_df, out_csv)
message("[OK] Wrote new high-res metadata CSV: ", out_csv)

# quick preview
print(utils::head(out_df, 5))

# ---------------------------
# README & done
# ---------------------------
readme <- paste0(
  "High-resolution env prep complete.\n",
  "CRS: ", target_crs, "\n",
  "Resolution (m): ", target_res, "\n",
  "Variables written: ", paste(keep, collapse = ", "), "\n",
  "Metadata CSV: ", out_csv, "\n"
)
writeLines(readme, con = file.path(opt$outdir, "HIRES_README.txt"))
message("[OK] Done. Layers: ", paste(keep, collapse = ", "))