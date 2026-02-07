#!/usr/bin/env Rscript
# ============================================================
# Step10_SCBD_Figures.R  (clean classes + full & "meaningful" figs)
# ------------------------------------------------------------
# Input : SCBD CSV with columns: class, order, family, genus, species, scbd
# Output: Fig1/2 (PNG+SVG) + SCBD_family_summary.csv + SCBD_genus_summary.csv
#         + cleaned SCBD CSV with harmonized 'class'
#         + optional "meaningful" versions of both figures
# Usage :
#   Rscript Step10_SCBD_Figures.R <SCBD_in.csv> <outdir>
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel", repos="https://cloud.r-project.org")
  library(ggrepel)
  if (!requireNamespace("magrittr", quietly=TRUE)) install.packages("magrittr", repos="https://cloud.r-project.org")
  library(magrittr)   # for %>%
})

# -------------------- knobs for "meaningful" filter --------------------
cum_share_threshold <- 0.80   # keep families until this cumulative share of total SCBD
min_abs_scbd        <- 0.001  # and also keep any family >= this absolute SCBD
min_n_genera        <- 2      # require ≥ this many genera per family (set to 1 to disable)

# -------------------- helpers --------------------
msg <- function(...) cat(sprintf(...), "\n")

# within-facet reordering for x labels
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  scale_x_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

# minimal cheat sheet where we *know* the class regardless of file hiccups
family_to_class_cheats <- c(
  # Arachnids
  "Philodromidae"="Arachnida","Salticidae"="Arachnida","Gnaphosidae"="Arachnida",
  "Thomisidae"="Arachnida","Lycosidae"="Arachnida","Scorpionidae"="Arachnida",
  "Buthidae"="Arachnida","Ixodidae"="Arachnida",
  # Insects (common families in your set)
  "Chironomidae"="Insecta","Dolichopodidae"="Insecta","Muscidae"="Insecta",
  "Sarcophagidae"="Insecta","Calliphoridae"="Insecta","Syrphidae"="Insecta",
  "Miridae"="Insecta","Cicadellidae"="Insecta","Aleyrodidae"="Insecta",
  "Formicidae"="Insecta","Pteromalidae"="Insecta","Encyrtidae"="Insecta",
  "Bruchinae"="Insecta","Chrysomelidae"="Insecta","Carabidae"="Insecta",
  "Tortricidae"="Insecta","Gracillariidae"="Insecta","Zygaenidae"="Insecta",
  # Collembola
  "Entomobryidae"="Collembola","Hypogastruridae"="Collembola"
)

harmonize_class <- function(df) {
  df <- df %>%
    mutate(across(c(class, order, family, genus, species), ~trimws(as.character(.)))) %>%
    mutate(
      class  = ifelse(is.na(class)  | !nzchar(class),  NA, class),
      order  = ifelse(is.na(order)  | !nzchar(order),  NA, order),
      family = ifelse(is.na(family) | !nzchar(family), NA, family),
      genus  = ifelse(is.na(genus)  | !nzchar(genus),  NA, genus),
      species= ifelse(is.na(species)| !nzchar(species), "sp.", species)
    )

  # 1) cheat-sheet override
  cheat_df <- tibble(family = names(family_to_class_cheats),
                     class_fix = unname(family_to_class_cheats))
  df <- df %>%
    left_join(cheat_df, by = "family") %>%
    mutate(class = ifelse(!is.na(class_fix), class_fix, class)) %>%
    select(-class_fix)

  # 2) majority vote class per family (resolve conflicts)
  fam_vote <- df %>%
    filter(!is.na(family)) %>%
    count(family, class, name = "n") %>%
    group_by(family) %>%
    slice_max(n, with_ties = TRUE) %>%
    ungroup() %>%
    group_by(family) %>%
    mutate(tie = n > 0 & n == max(n) & dplyr::n() > 1) %>%
    ungroup()

  tied <- fam_vote %>% filter(tie) %>% distinct(family)
  if (nrow(tied) > 0) {
    msg("[WARN] Families with tied majority class (left unchanged): %s",
        paste(tied$family, collapse = ", "))
  }

  fam_majority <- fam_vote %>%
    filter(!tie) %>%
    group_by(family) %>% slice_max(n, with_ties = FALSE) %>% ungroup() %>%
    select(family, class_majority = class)

  df %>%
    left_join(fam_majority, by = "family") %>%
    mutate(class = ifelse(is.na(class) & !is.na(class_majority), class_majority, class)) %>%
    select(-class_majority) %>%
    mutate(
      class  = ifelse(is.na(class),  "Unknown_class",  class),
      order  = ifelse(is.na(order),  "Unknown_order",  order),
      family = ifelse(is.na(family), "Unknown_family", family),
      genus  = ifelse(is.na(genus),  "Unknown_genus",  genus)
    )
}

# -------------------- args --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript Step10_SCBD_Figures.R <SCBD_in.csv> <outdir>", call. = FALSE)

infile <- args[1]; outdir <- args[2]
msg("[INFO] Input file: %s", infile)
msg("[INFO] Output directory: %s", outdir)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------- read --------------------
df <- read_csv(infile, show_col_types = FALSE)
names(df) <- tolower(names(df))
need <- c("class","order","family","genus","species","scbd")
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop(sprintf("Missing required columns: %s", paste(miss, collapse=", ")), call.=FALSE)

df <- df %>% mutate(scbd = suppressWarnings(as.numeric(scbd))) %>% filter(!is.na(scbd))

# -------------------- fix class inconsistencies --------------------
df <- harmonize_class(df)

# write cleaned CSV
clean_csv <- file.path(outdir, "SCBD_HIGHER_TAXA_classfixed.csv")
write_csv(df, clean_csv)
msg("[INFO] Wrote cleaned CSV: %s", clean_csv)

# -------------------- summaries --------------------
fam_tot <- df %>%
  group_by(class, order, family) %>%
  summarise(total_scbd = sum(scbd, na.rm = TRUE),
            n_genera   = n_distinct(genus),
            .groups = "drop") %>%
  arrange(desc(total_scbd))

fam_order <- fam_tot %>% arrange(desc(total_scbd)) %>% pull(family) %>% unique()

gen_tot <- df %>%
  group_by(class, order, family, genus) %>%
  summarise(genus_scbd = sum(scbd, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(genus_scbd))

write_csv(fam_tot, file.path(outdir, "SCBD_family_summary.csv"))
write_csv(gen_tot, file.path(outdir, "SCBD_genus_summary.csv"))

# -------------------- Figure 1: All families --------------------
K_label <- 25
fam_for_plot <- fam_tot %>% mutate(family = factor(family, levels = fam_order))

p_fam <- ggplot(fam_for_plot, aes(x = total_scbd, y = family, fill = order)) +
  geom_col(width = 0.75) +
  geom_text(
    data = fam_for_plot %>% slice_max(total_scbd, n = K_label),
    aes(label = sprintf("%.3g", total_scbd)),
    hjust = -0.1, size = 2.8, color = "black"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Total SCBD by Family", x = "Total SCBD Contribution", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "Fig1_SCBD_by_family.png"), p_fam, width = 8, height = 14, dpi = 300)
ggsave(file.path(outdir, "Fig1_SCBD_by_family.svg"), p_fam, width = 8, height = 14)

# -------------------- Figure 2: Top Genera within Top Families (ALL) --------------------
set.seed(42)
J_fam <- 15
G_gen  <- 8

top_fams_all <- fam_tot %>% slice_max(total_scbd, n = J_fam) %>% pull(family)

gen_top_all <- gen_tot %>%
  filter(family %in% top_fams_all) %>%
  group_by(family) %>% slice_max(genus_scbd, n = G_gen, with_ties = FALSE) %>% ungroup() %>%
  distinct(family, genus, .keep_all = TRUE) %>%
  mutate(
    family = factor(family, levels = top_fams_all),
    genus  = reorder_within(genus, genus_scbd, family),
    lab    = gsub("___.*$", "", genus)
  )

p_gen_all <- ggplot(gen_top_all, aes(x = genus, y = genus_scbd, color = order)) +
  geom_segment(aes(x = genus, xend = genus, y = 0, yend = genus_scbd),
               linewidth = 0.5, alpha = 0.55) +
  geom_point(size = 2.1) +
  ggrepel::geom_text_repel(
    aes(label = lab),
    size = 2.4,
    direction = "y",
    box.padding = 0.3,
    point.padding = 0.15,
    force = 2.0,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.alpha = 0.7,
    nudge_y = 0.0002,
    hjust = 1.02,
    nudge_x = -0.25
  ) +
  facet_wrap(~ family, scales = "free_y", ncol = 3) +
  scale_x_reordered() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.35))) +
  labs(title = "Top Genera Contributions within Top Families",
       x = NULL, y = "Genus-level SCBD") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(outdir, "Fig2_TopGenera_within_TopFamilies.png"), p_gen_all, width = 12, height = 10, dpi = 300)
ggsave(file.path(outdir, "Fig2_TopGenera_within_TopFamilies.svg"), p_gen_all, width = 12, height = 10)

# ---------# -------------------- "Meaningful" subset (optional, on by default) --------------------
fam_ranked <- fam_tot %>%
  arrange(desc(total_scbd)) %>%
  mutate(cum_scbd  = cumsum(total_scbd),
         cum_share = cum_scbd / sum(total_scbd, na.rm = TRUE))

meaningful_fams <- fam_ranked %>%
  filter(cum_share <= cum_share_threshold | total_scbd >= min_abs_scbd) %>%
  filter(n_genera >= min_n_genera) %>%
  arrange(desc(total_scbd))

if (nrow(meaningful_fams) == 0) {
  meaningful_fams <- fam_ranked %>% slice_head(n = 20)
}

# --- DEDUP for plotting: keep one row per family (the highest-SCBD instance) ---
meaningful_fams_plot <- meaningful_fams %>%
  arrange(desc(total_scbd)) %>%
  distinct(family, .keep_all = TRUE)

# Figure 1 (meaningful only)
meaningful_levels <- rev(meaningful_fams_plot$family)
p_fam_meaningful <- ggplot(
  meaningful_fams_plot %>% mutate(family = factor(family, levels = meaningful_levels)),
  aes(x = total_scbd, y = family, fill = order)
) +
  geom_col(width = 0.7, color = "black", alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", total_scbd)),
            hjust = -0.18, size = 3.0, color = "black", fontface = "italic") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15)),
    limits = c(0, max(meaningful_fams_plot$total_scbd) * 1.15)
  ) +
  labs(title = "Families with Meaningful SCBD Contribution",
       subtitle = paste0("Rules: ≤", round(cum_share_threshold*100),
                         "% cumulative OR ≥", min_abs_scbd,
                         " SCBD; ≥", min_n_genera, " genera"),
       x = "Total SCBD Contribution", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10, color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(outdir, "Fig1_SCBD_by_family_meaningful.png"), p_fam_meaningful, width = 8, height = 9, dpi = 400)
ggsave(file.path(outdir, "Fig1_SCBD_by_family_meaningful.svg"), p_fam_meaningful, width = 8, height = 9)

# Figure 2 (meaningful only; same labeling)
top_fams_mean <- as.character(meaningful_fams_plot$family)
gen_tot_mean  <- gen_tot %>% filter(family %in% top_fams_mean)

gen_top_mean <- gen_tot_mean %>%
  group_by(family) %>% slice_max(genus_scbd, n = G_gen, with_ties = FALSE) %>% ungroup() %>%
  distinct(family, genus, .keep_all = TRUE) %>%
  mutate(
    family = factor(family, levels = top_fams_mean),
    genus  = reorder_within(genus, genus_scbd, family),
    lab    = gsub("___.*$", "", genus)
  )

p_gen_mean <- ggplot(gen_top_mean, aes(x = genus, y = genus_scbd, color = order)) +
  geom_segment(aes(x = genus, xend = genus, y = 0, yend = genus_scbd),
               linewidth = 0.5, alpha = 0.55) +
  geom_point(size = 2.1) +
  ggrepel::geom_text_repel(
    aes(label = lab),
    size = 2.4,
    direction = "y",
    box.padding = 0.3,
    point.padding = 0.15,
    force = 2.0,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.alpha = 0.7,
    nudge_y = 0.0002,
    hjust = 1.02,
    nudge_x = -0.25
  ) +
  facet_wrap(~ family, scales = "free_y", ncol = 3) +
  scale_x_reordered() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.35))) +
  labs(title = "Top Genera within Families with Meaningful SCBD",
       x = NULL, y = "Genus-level SCBD") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(outdir, "Fig2_TopGenera_within_MeaningfulFamilies.png"), p_gen_mean, width = 12, height = 10, dpi = 300)
ggsave(file.path(outdir, "Fig2_TopGenera_within_MeaningfulFamilies.svg"), p_gen_mean, width = 12, height = 10)

# -------------------- done --------------------
msg("[OK] Wrote:")
msg(" - %s", file.path(outdir, "Fig1_SCBD_by_family.png"))
msg(" - %s", file.path(outdir, "Fig1_SCBD_by_family.svg"))
msg(" - %s", file.path(outdir, "Fig2_TopGenera_within_TopFamilies.png"))
msg(" - %s", file.path(outdir, "Fig2_TopGenera_within_TopFamilies.svg"))
msg(" - %s", file.path(outdir, "Fig1_SCBD_by_family_meaningful.png"))
msg(" - %s", file.path(outdir, "Fig1_SCBD_by_family_meaningful.svg"))
msg(" - %s", file.path(outdir, "Fig2_TopGenera_within_MeaningfulFamilies.png"))
msg(" - %s", file.path(outdir, "Fig2_TopGenera_within_MeaningfulFamilies.svg"))
msg(" - %s", file.path(outdir, "SCBD_family_summary.csv"))
msg(" - %s", file.path(outdir, "SCBD_genus_summary.csv"))
msg(" - %s", clean_csv)