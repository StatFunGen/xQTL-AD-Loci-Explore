#!/usr/bin/env Rscript
## Check a release directory against the expectations of the 188-loci build.
## usage: Rscript validate_outputs.R <release_dir>
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript validate_outputs.R <release_dir>", call. = FALSE)
release <- args[1]
if (!dir.exists(release)) stop("no such release directory: ", release, call. = FALSE)

fail <- character()
note <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fail <<- c(fail, msg)
}

cat("release:", release, "\n\n")

## ---- expected outputs ------------------------------------------------------
expected <- c("res_AD_variants_xQTL.csv.gz",
              "res_all_single_gwas_finemapping_cs50orgreater.csv.gz",
              "res_coloc_AD_xQTL_unified.csv.gz",
              "res_AD_XWAS.csv.gz",
              "res_AD_cTWAS_pip075.csv.gz")
for (f in expected) note(file.exists(file.path(release, f)), paste("present:", f))

xl <- list.files(release, pattern = "^unified_AD_loci_xQTL_summary.*\\.xlsx$")
note(length(xl) > 0, "present: unified_AD_loci_xQTL_summary*.xlsx")

## ---- content checks --------------------------------------------------------
vf <- file.path(release, "res_AD_variants_xQTL.csv.gz")
if (file.exists(vf)) {
  v <- fread(vf, select = c("locus_index", "gene_name", "top_confidence"))
  n_loci <- uniqueN(v[!is.na(locus_index)]$locus_index)
  note(n_loci == 188, sprintf("locus count is 188 (found %d)", n_loci))

  t <- v[!is.na(top_confidence) & top_confidence != "" &
         !is.na(gene_name) & gene_name != ""]
  g <- t[, .(tier = min(top_confidence)), by = .(gene = gene_name)]
  note(nrow(g) == 508, sprintf("tiered genes is 508 (found %d)", nrow(g)))

  tb <- table(g$tier)
  cat("\n  tier distribution:", paste(sprintf("%s=%d", names(tb), tb), collapse = " "), "\n")
  note("T6" %in% names(tb) && tb[["T6"]] > 0,
       "T6 is populated (genes with gene-level evidence only)")
  cat("\n")
}

if (length(fail)) {
  cat(sprintf("FAILED %d check(s)\n", length(fail)))
  quit(status = 1)
}
cat("all checks passed\n")
