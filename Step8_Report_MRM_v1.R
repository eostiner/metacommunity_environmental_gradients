#!/usr/bin/env Rscript
# Step8_Report_MRM_v1.R — Simple bars; two R² values outside at bar end
# Panel A: per-predictor bars (colored). Label outside: "total R²=.. | individual R²=.."
# Panel B: per-model bars (neutral color). Label outside: "total R²=.. | best individual R²=.."
# Notes:
#   individual R² (per predictor) = single-predictor fit from MRM_single_predictor.csv
#   total R² (per bar)            = model R² from the corresponding CSV
#   best individual R² (Panel B)  = max of the included predictors’ individual R²

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(stringr)
  library(tidyr); library(ggplot2); library(cowplot)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--mrm_dir", type="character", help="Step8_MRM directory with MRM_*.csv"),
  make_option("--outdir",  type="character", default="Step8_MRM/REPORT"),
  make_option("--best_predictor", type="character", default="D_cost_combined_global"),
  make_option("--title", type="character", default="Distance–Dissimilarity Models (MRM)"),
  make_option("--dpi", type="integer", default=300)
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] mrm_dir=%s\n", opt$mrm_dir))

# ---------- Files ----------
single_fp <- file.path(opt$mrm_dir, "MRM_single_predictor.csv")
multi_fp  <- file.path(opt$mrm_dir, "MRM_multi_predictor.csv")
if (!file.exists(single_fp)) stop("Missing MRM_single_predictor.csv at: ", single_fp)
if (!file.exists(multi_fp))  stop("Missing MRM_multi_predictor.csv at: ", multi_fp)

single_raw <- suppressMessages(readr::read_csv(single_fp, show_col_types = FALSE))
multi_raw  <- suppressMessages(readr::read_csv(multi_fp,  show_col_types = FALSE))

# ---------- Normalize ----------
norm_single <- function(df){
  if (!"model" %in% names(df)) if ("predictor" %in% names(df)) df <- dplyr::rename(df, model = predictor)
  if (!"model" %in% names(df)) stop("[single] Need 'model' or 'predictor'. Have: ", paste(names(df), collapse=", "))
  if (!"R2" %in% names(df))    stop("[single] Missing 'R2'.")
  if (!"terms" %in% names(df)) df$terms <- ""
  df %>% select(model, R2, terms) %>% mutate(R2 = suppressWarnings(as.numeric(R2)))
}
norm_multi <- function(df){
  need <- c("model","R2","terms")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("[multi] Missing: ", paste(miss, collapse=", "))
  df %>% select(model, R2, terms) %>% mutate(R2 = suppressWarnings(as.numeric(R2)))
}
single <- norm_single(single_raw)
multi  <- norm_multi(multi_raw)

# Drop non-finite R2
drop_bad <- function(df, name){
  bad <- df %>% filter(!is.finite(R2) | is.na(R2))
  if (nrow(bad)) readr::write_csv(bad, file.path(opt$outdir, paste0("DROPPED_", name, "_missing_R2.csv")))
  df %>% filter(is.finite(R2))
}
single <- drop_bad(single, "single")
multi  <- drop_bad(multi,  "multi")

# ---------- Pretty names & parsing ----------
clean_name <- function(x){
  x %>%
    str_replace("^D_geo_km$", "Geographic distance") %>%
    str_replace("^D_env_scaled$", "Environmental distance") %>%
    str_replace("^D_cost_combined_global$", "Combined resistance") %>%
    str_replace("^D_cost_res_slope_global$", "Slope resistance") %>%
    str_replace("^D_cost_res_rough_global$", "Roughness resistance")
}
pretty_multi <- function(model){
  s <- model %>% str_replace_all("as\\.dist\\(|\\)", "") %>% str_replace("^D_comm\\s*~\\s*", "")
  parts <- str_split(s, "\\+")[[1]] %>% str_trim() %>% clean_name()
  paste(parts, collapse = " + ")
}
get_parts <- function(model_clean_label){
  str_split(model_clean_label, "\\s*\\+\\s*")[[1]] %>% str_trim()
}

# ---------- Panel A stats (individual & total per predictor) ----------
# If a predictor appears twice in single CSV, treat:
#   total R²      = max(R2)
#   individual R² = min(R2)
single_stats <- single %>%
  transmute(pred = clean_name(model), R2 = pmin(pmax(R2, 0), 1)) %>%
  group_by(pred) %>%
  summarise(total_R2 = max(R2, na.rm = TRUE),
            individual_R2 = min(R2, na.rm = TRUE),
            .groups = "drop")

# ---------- Panel B data (total per model; best individual among included) ----------
ind_map <- setNames(single_stats$individual_R2, single_stats$pred)

multi_plot <- multi %>%
  transmute(model_clean = vapply(model, pretty_multi, FUN.VALUE = character(1)),
            R2 = pmin(pmax(R2, 0), 1)) %>%
  group_by(model_clean) %>%
  summarise(total_R2 = max(R2, na.rm=TRUE), .groups="drop") %>%
  mutate(best_individual_R2 = vapply(model_clean, function(lbl){
           parts <- get_parts(lbl)
           vals <- as.numeric(ind_map[parts])
           vals <- vals[is.finite(vals)]
           ifelse(length(vals)==0, NA_real_, max(vals, na.rm=TRUE))
         }, FUN.VALUE = numeric(1))) %>%
  arrange(desc(total_R2))

# ---------- Colors ----------
concepts <- c("Geographic distance","Environmental distance","Combined resistance","Slope resistance","Roughness resistance")
if (requireNamespace("RColorBrewer", quietly = TRUE)) {
  base_pal <- RColorBrewer::brewer.pal(5, "Set2")
} else {
  base_pal <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")
}
cols <- setNames(base_pal, concepts)

# ---------- Theme & layout helpers ----------
fmt <- function(x) sprintf("%.2f", x)

base_theme <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        plot.margin = margin(6, 46, 6, 10))  # big right margin for outside labels

# Add right headroom for text and disable clipping
x_upper <- max(1, max(single_stats$total_R2, na.rm=TRUE), max(multi_plot$total_R2, na.rm=TRUE)) + 0.26
x_upper <- min(x_upper, 1.50)

# ---------- Panel A (simple colored bars; two values outside) ----------
p_single <- ggplot(single_stats, aes(x = reorder(pred, total_R2), y = total_R2, fill = pred)) +
  geom_col(width = 0.70, color = "white") +
  geom_text(aes(y = total_R2 + 0.035,
                label = paste0("total R²=", fmt(total_R2), "  |  individual R²=", fmt(individual_R2))),
            hjust = 0, size = 3.2) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = cols, breaks = concepts, drop = FALSE) +
  scale_y_continuous(limits = c(0, x_upper), expand = expansion(mult = c(0, 0.02))) +
  labs(title = paste0(opt$title, ": Single predictor"),
       x = NULL, y = "R²") +
  base_theme

# ---------- Panel B (simple neutral bars; two values outside) ----------
p_multi <- ggplot(multi_plot, aes(x = reorder(model_clean, total_R2), y = total_R2)) +
  geom_col(width = 0.70, fill = "grey70", color = "white") +
  geom_text(aes(y = total_R2 + 0.035,
                label = paste0("total R²=", fmt(total_R2),
                               "  |  best individual R²=", fmt(best_individual_R2))),
            hjust = 0, size = 3.0) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, x_upper), expand = expansion(mult = c(0, 0.02))) +
  labs(title = paste0(opt$title, ": Multi predictor"),
       x = NULL, y = "R²") +
  base_theme

# ---------- Layout & save ----------
h_single <- max(3.8, 0.60 * nrow(single_stats) + 2.0)
h_multi  <- max(4.8, 0.60 * nrow(multi_plot)  + 2.8)
total_h  <- h_single + h_multi + 1.3
page_w   <- 12.8  # wider page to ensure labels never run off

fig <- cowplot::plot_grid(p_single, p_multi, ncol = 1, labels = c("A","B"),
                          rel_heights = c(h_single, h_multi))

ggsave(file.path(opt$outdir, "MRM_model_R2_comparison.png"), fig,
       width = page_w, height = total_h, dpi = opt$dpi)
ggsave(file.path(opt$outdir, "MRM_model_R2_comparison.pdf"),  fig,
       width = page_w, height = total_h)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(opt$outdir, "MRM_model_R2_comparison.svg"), fig,
         width = page_w, height = total_h, device = svglite::svglite)
}

# ---------- Export tables ----------
write_csv(single_stats, file.path(opt$outdir, "MRM_single_predictor_individual_and_total_R2.csv"))
write_csv(multi_plot,  file.path(opt$outdir, "MRM_multi_predictor_total_and_best_individual_R2.csv"))

cat(sprintf("[OK] Wrote report to: %s\n", opt$outdir))