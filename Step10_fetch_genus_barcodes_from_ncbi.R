#!/usr/bin/env Rscript
# Fetch short barcode sequences from NCBI for a list of genera.
# Uses rentrez only (CRAN). Outputs per-genus FASTA and a combined FASTA + mapping CSV.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(rentrez)
})

opt_list <- list(
  make_option("--csv",       type="character", help="Input CSV with a genus column."),
  make_option("--genus_col", type="character", default="genus", help="Column name containing genus names."),
  make_option("--marker",    type="character", default="COI",
              help='Target marker. One of: COI, 16S, 18S, 28S.'),
  make_option("--minlen",    type="integer", default=150, help="Minimum sequence length."),
  make_option("--maxlen",    type="integer", default=900, help="Maximum sequence length."),
  make_option("--retmax",    type="integer", default=100, help="Max records to fetch per genus."),
  make_option("--email",     type="character", default=Sys.getenv("NCBI_EMAIL"), help="NCBI email (or set NCBI_EMAIL)."),
  make_option("--api_key",   type="character", default=Sys.getenv("NCBI_API_KEY"), help="NCBI API key (or set NCBI_API_KEY)."),
  make_option("--outdir",    type="character", help="Output directory.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$csv), !is.null(opt$outdir))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- config NCBI ----
if (nzchar(opt$email)) {
  # Try both old and new rentrez styles
  if ("rentrez_env" %in% ls("package:rentrez")) {
    # Newer rentrez versions store email in global options
    options(rentrez.email = opt$email)
  } else {
    # Older versions had a helper
    if (exists("set_entrez_email", where = asNamespace("rentrez"), mode = "function"))
      get("set_entrez_email", asNamespace("rentrez"))(opt$email)
  }
  message("[INFO] Using email: ", opt$email)
}

if (nzchar(opt$api_key)) {
  if ("rentrez_env" %in% ls("package:rentrez")) {
    options(rentrez.key = opt$api_key)
  } else {
    if (exists("set_entrez_key", where = asNamespace("rentrez"), mode = "function"))
      get("set_entrez_key", asNamespace("rentrez"))(opt$api_key)
  }
  message("[INFO] Using NCBI API key.")
}

# ---- marker query terms ----
marker_terms <- list(
  COI = '(cox1[Gene] OR coxI[Gene] OR COI[All Fields] OR "cytochrome oxidase subunit I"[Title])',
  `16S` = '(16S[All Fields] OR rrnS[Gene] OR "small subunit ribosomal RNA"[Title])',
  `18S` = '(18S[All Fields] OR SSU[All Fields] OR "18S ribosomal RNA"[Title])',
  `28S` = '(28S[All Fields] OR LSU[All Fields] OR "28S ribosomal RNA"[Title])'
)
if (!(opt$marker %in% names(marker_terms))) {
  stop("Unknown --marker. Use one of: COI, 16S, 18S, 28S")
}
marker_q <- marker_terms[[opt$marker]]

# ---- read CSV safely ----
dat <- suppressMessages(read_csv(opt$csv, show_col_types = FALSE, guess_max = 1e6))
cn <- names(dat)
# try to normalize probable genus column if needed
if (!(opt$genus_col %in% cn)) {
  alt <- cn[grepl("genus", tolower(cn))]
  if (length(alt) == 1) {
    message("[WARN] --genus_col not found. Using detected column: ", alt)
    opt$genus_col <- alt
  } else {
    stop("genus_col not found in the CSV. Columns are: ", paste(cn, collapse=", "))
  }
}

gens <- dat[[opt$genus_col]] %>%
  as.character() %>%
  trimws() %>%
  .[nchar(.) > 0] %>%
  unique()

# Keep only the 1st word (genus token), drop anything clearly not a genus
gens <- sub("^([A-Za-z]+).*", "\\1", gens)
gens <- gens[grepl("^[A-Za-z]+$", gens)]
gens <- unique(gens)

if (length(gens) == 0) stop("No parsable genus names found.")

# ---- helpers ----
fasta_wrap <- function(seq, width=70) {
  paste(strwrap(seq, width=width), collapse="\n")
}

fetch_one_genus <- function(genus) {
  org_q <- paste0(genus, "[Organism]")
  len_q <- sprintf("(%d:%d[SLEN])", opt$minlen, opt$maxlen)
  mt    <- marker_q

  # prefer mitochondria for COI
  mito_q <- if (opt$marker == "COI") '(mitochondrion[Filter] OR mitochondria[All Fields])' else ""
  q <- paste(org_q, mt, len_q, mito_q)
  q <- gsub("  +", " ", q)

  # search
  srch <- tryCatch(
    entrez_search(db="nuccore", term=q, retmax=opt$retmax, use_history=TRUE),
    error=function(e) NULL
  )
  if (is.null(srch) || srch$count == 0) return(list(records=NULL, fasta=NULL))

  # fetch fasta
  fasta_txt <- tryCatch(
    entrez_fetch(db="nuccore", web_history=srch$web_history, rettype="fasta", retmax=opt$retmax),
    error=function(e) NULL
  )
  if (is.null(fasta_txt) || !nzchar(fasta_txt)) return(list(records=NULL, fasta=NULL))

  # parse headers minimally
  lines <- strsplit(fasta_txt, "\n")[[1]]
  hdr_i <- grep("^>", lines)
  if (length(hdr_i) == 0) return(list(records=NULL, fasta=NULL))

  recs <- vector("list", length(hdr_i))
  out_fa <- character(length(hdr_i))
  for (k in seq_along(hdr_i)) {
    i <- hdr_i[k]
    j <- if (k < length(hdr_i)) hdr_i[k+1] - 1 else length(lines)
    hdr <- sub("^>", "", lines[i])
    seq <- paste(lines[(i+1):j], collapse="")
    acc <- sub(" .*", "", hdr)
    len <- nchar(gsub("\\s+", "", seq))
    recs[[k]] <- data.frame(
      genus   = genus,
      accession = acc,
      length  = len,
      header  = hdr,
      stringsAsFactors = FALSE
    )
    out_fa[k] <- paste0(">", acc, " | genus=", genus, " | len=", len, "\n", fasta_wrap(seq))
  }
  list(records = bind_rows(recs), fasta = paste(out_fa, collapse="\n"))
}

# ---- main loop ----
all_records <- list()
all_fasta   <- list()
missed      <- character(0)

message("[INFO] Querying NCBI for ", length(gens), " genera, marker=", opt$marker,
        ", len=[", opt$minlen, ",", opt$maxlen, "], retmax=", opt$retmax)

for (g in gens) {
  cat("[INFO] ", g, " … ")
  res <- fetch_one_genus(g)
  if (is.null(res$fasta)) {
    cat("0 hits\n")
    missed <- c(missed, g)
  } else {
    cat("ok\n")
    all_records[[g]] <- res$records
    all_fasta[[g]]   <- res$fasta
    # per-genus FASTA
    writeLines(res$fasta, file.path(opt$outdir, paste0("barcode_", opt$marker, "_", g, ".fasta")))
  }
}

if (length(all_records)) {
  REC <- bind_rows(all_records)
  write_csv(REC, file.path(opt$outdir, paste0("barcode_", opt$marker, "_mapping.csv")))
  writeLines(paste(unlist(all_fasta), collapse="\n"),
             file.path(opt$outdir, paste0("barcode_", opt$marker, "_ALL.fasta")))

  # pick a representative per genus: closest to 650bp for COI, else median length
  target <- if (opt$marker == "COI") 650L else round((opt$minlen + opt$maxlen)/2)
  reps <- REC %>%
    group_by(genus) %>%
    slice_min(abs(length - target), with_ties = FALSE) %>%
    ungroup()
  # subset FASTA for representatives
  fasta_all <- paste(unlist(all_fasta), collapse="\n")
  # crude filter by accession lines
  keep <- reps$accession
  blocks <- unlist(strsplit(fasta_all, "(?=>)", perl=TRUE))
  blocks <- blocks[nzchar(blocks)]
  keep_blocks <- blocks[vapply(keep, function(acc) any(grepl(paste0("^>", acc, "\\b"), blocks)), logical(1L))]
  writeLines(paste(keep_blocks, collapse=""), file.path(opt$outdir, paste0("barcode_", opt$marker, "_REPRESENTATIVES.fasta")))
} else {
  warning("No sequences retrieved for any genus.")
}

if (length(missed)) {
  writeLines(missed, file.path(opt$outdir, "no_hits_genera.txt"))
  message("[WARN] No hits for genera: ", paste(missed, collapse=", "))
}

message("[OK] Done. Outputs in: ", opt$outdir)