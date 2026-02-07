#!/usr/bin/env Rscript
# Usage:
# Rscript Step11_extract_env_from_geotiffs.R samples.csv out.csv buffer_meters raster1.tif[:alias] [raster2.tif[:alias]] ...
suppressPackageStartupMessages({
  if (!requireNamespace("terra", quietly=TRUE)) {
    stop("Please install the 'terra' package: install.packages('terra')", call.=FALSE)
  }
  library(terra)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage:\n  Rscript Step11_extract_env_from_geotiffs.R <samples.csv> <out.csv> <buffer_meters> <r1.tif[:name]> [r2.tif[:name]] ...\n", call.=FALSE)
}
in_csv  <- args[1]
out_csv <- args[2]
buf_m   <- as.numeric(args[3])
rast_specs <- args[-(1:3)]

pts <- read.csv(in_csv, check.names = FALSE)
stopifnot(all(c("sample_id","lat","lon") %in% names(pts)))

v <- vect(pts, geom=c("lon","lat"), crs="EPSG:4326")  # carries all columns (incl. sample_id)

extract_one <- function(spec) {
  if (grepl(":", spec) && !grepl("^https?://", spec)) {
    file <- sub(":.*$", "", spec)
    alias <- sub("^.*?:", "", spec)
  } else {
    file <- spec; alias <- NULL
  }
  r <- rast(file)
  v2 <- project(v, crs(r))
  if (!is.na(buf_m) && buf_m > 0) {
    if (is.lonlat(r)) {
      deg <- buf_m / 111320
      vv <- buffer(v2, width = deg)
    } else {
      vv <- buffer(v2, width = buf_m)
    }
    e <- terra::extract(r, vv, fun=mean, na.rm=TRUE)
  } else {
    e <- terra::extract(r, v2)
  }
  out <- as_tibble(e[,-1, drop=FALSE])
  if (!is.null(alias)) {
    names(out) <- paste0(alias, if (ncol(out)>1) paste0("_b", seq_len(ncol(out))) else "")
  } else {
    names(out) <- make.names(names(out), unique=TRUE)
  }
  out
}

vals_list <- lapply(rast_specs, extract_one)
vals <- bind_cols(vals_list)
res <- bind_cols(pts, vals)
write.csv(res, out_csv, row.names = FALSE)
cat("[OK] wrote", out_csv, "\n")
