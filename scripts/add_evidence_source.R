#!/usr/bin/env Rscript
## add_evidence_source.R -- register a new evidence source in the AD loci table.
##
## Converts a ColocBoost-style coloc summary .bed (one row per variant, with
## semicolon-delimited event lists) into per-context toploci files, and registers
## them in metadata_analysis.csv / contexts_metadata.csv.
##
##   Rscript add_evidence_source.R --name transmap \
##       --bed trans_xQTL_only_colocalization_summary_table.bed
##
## Then rerun build_AD_locus_table.R. The build auto-detects any Method in
## metadata_analysis.csv that points at per-context toploci files and has no
## dedicated loader, so registering the rows is the only step needed. Safe to
## re-run: conversion is skipped if outputs exist, and metadata rows for the
## same --name are replaced, not duplicated.
suppressPackageStartupMessages({library(data.table)})

ROOT <- Sys.getenv("AD_LOCI_ROOT", unset = getwd())
.f <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
## This script lives in scripts/; the repository root is its parent directory.
REPO <- if (length(.f)) dirname(normalizePath(dirname(.f[1]))) else getwd()
STG <- Sys.getenv("AD_LOCI_CONFIG", unset = file.path(REPO, "config"))
a <- commandArgs(TRUE); getarg <- function(k,d=NA){i<-match(k,a); if(is.na(i)) d else a[i+1]}
NAME <- getarg("--name"); BED <- getarg("--bed"); DIR <- getarg("--dir", paste0("trans_toploci_",NAME))
STAMP <- format(Sys.Date(), "%Y%m%d")
if (is.na(NAME) || is.na(BED)) stop("need --name and --bed")
setwd(ROOT)

## ---- the converter (was toploci_explode.awk) -----------------------------
## Source columns: 1 #chr 2 start 3 end 4 a1 5 a2 6 variant_ID 7 region_ID
##                 8 event_ID(list) 9 cos_ID 10 vcp 11 cos_npc 14 zscore(list)
## Emitted toploci columns: #chr start end a1 a2 variant_ID region_ID event_ID
##                 cs_coverage_0.95(<-cosN) purity(<-cos_npc) PIP(<-vcp) z(<-zscore[i])
AWK_PROG <- c(
 'BEGIN{ FS=OFS="\\t" }',
 'NR==1 { next }',
 '{',
 '  cosn = 0',
 '  if (match($9, /cos[0-9]+/)) cosn = substr($9, RSTART+3, RLENGTH-3)+0',
 '  if (cosn == 0) next',
 '  n  = split($8,  ev, /; */)',
 '  nz = split($14, zs, /; */)',
 '  for (i=1; i<=n; i++) {',
 '    e = ev[i]; if (e == "") continue',
 '    ctx = e; sub(/_(cis|trans)_.*$/, "", ctx)',
 '    if (ctx == "" || ctx == e) continue',
 '    gsub(/[^A-Za-z0-9_.-]/, "", ctx)',
 '    z = (i <= nz) ? zs[i] : ""',
 '    f = D ctx "." NAME ".exported.toploci.bed"',
 '    if (!(ctx in seen)) {',
 '      print "#chr","start","end","a1","a2","variant_ID","region_ID","event_ID","cs_coverage_0.95","purity","PIP","z" > f',
 '      seen[ctx] = 1 }',
 '    print $1,$2,$3,$4,$5,$6,$7,e,cosn,$11,$10,z >> f',
 '  }',
 '}')

if (!dir.exists(DIR) || !length(list.files(DIR, pattern="\\.bed\\.gz$"))) {
  dir.create(DIR, showWarnings=FALSE)
  prog <- tempfile(fileext=".awk"); writeLines(AWK_PROG, prog)
  cat("[1] converting", BED, "->", DIR, "(several minutes)\n")
  st <- system2("awk", c("-v", paste0("D=",DIR,"/"), "-v", paste0("NAME=",NAME), "-f", prog, BED))
  if (st != 0) stop("awk failed with status ", st)
  for (f in list.files(DIR, pattern="\\.bed$", full.names=TRUE)) system2("gzip", c("-f", f))
} else cat("[1] per-context files already present in", DIR, "- skipping conversion\n")

files <- list.files(DIR, pattern="\\.bed\\.gz$")
stopifnot(length(files) > 0)
ctxs <- sub("\\..*$", "", files)
cat("    contexts:", length(ctxs), "->", paste(ctxs, collapse=", "), "\n")

## ---- reuse a context label the pipeline knows, else mint one -------------
ctxmeta <- fread(file.path(STG,"contexts_metadata.csv"), colClasses="character")
known <- vapply(ctxs, function(c){h <- grep(paste0("^",c,"_"), ctxmeta$context, value=TRUE)
                                 if (length(h)) h[1] else paste0(c,"_",NAME)}, "")
cat("[2] reusing", sum(known %in% ctxmeta$context), "existing labels;",
    sum(!known %in% ctxmeta$context), "new\n")

## ---- register ------------------------------------------------------------
mfile <- file.path(STG,"metadata_analysis.csv")
b <- paste0(mfile,".bak_pre_",NAME,"_",STAMP); if (!file.exists(b)) file.copy(mfile, b)
mtd <- fread(mfile, colClasses="character"); n0 <- nrow(mtd); mtd <- mtd[Method != NAME]
rows <- data.table(context=unname(known), `Data Type`="", Cohort="", Modality="", Method=NAME,
                   Path=file.path(DIR, files), context_broad="", summary_file="",
                   variant_level_method="TRUE", summary_file_ad="")
mtd2 <- rbind(mtd, rows[, names(mtd), with=FALSE]); fwrite(mtd2, mfile)
cat("[3] metadata_analysis.csv:", n0, "->", nrow(mtd2), "rows (", nrow(rows), NAME, ")\n")

need <- setdiff(unname(known), ctxmeta$context)
if (length(need)) {
  cf <- file.path(STG,"contexts_metadata.csv")
  cb <- paste0(cf,".bak_pre_",NAME,"_",STAMP); if (!file.exists(cb)) file.copy(cf, cb)
  add <- data.table(context_short=need, context_hi=need, context=need, context_trans=need,
                    context_coloc="", context_qr="", context_broad="", context_snsQTL="",
                    fill_color="#e6e6e6", font_color="black")
  fwrite(rbind(ctxmeta, add[, names(ctxmeta), with=FALSE]), cf)
  cat("[4] contexts_metadata.csv: added", length(need), "\n")
} else cat("[4] contexts_metadata.csv: no new contexts needed\n")
cat("\nDONE. Registered '",NAME,"' in metadata_analysis.csv; rerun build_AD_locus_table.R.\n", sep="")
