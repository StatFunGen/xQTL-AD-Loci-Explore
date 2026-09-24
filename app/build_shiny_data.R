#!/usr/bin/env Rscript
# build_shiny_data.R - refresh the AD Loci Explorer table from a pipeline release.
#
#   Rscript build_shiny_data.R <release_dir> [out.csv]
#
# ---------------------------------------------------------------------------
# WHY THIS IS A JOIN AND NOT A REBUILD
#
# The explorer's table mixes three kinds of column, with three different origins:
#
#   (A) LOCUS / VARIANT   - ADlocus, variant_ID, rsID, chr, pos, cV2F, p-values,
#                           GWAS sources. Produced by the integration pipeline and
#                           refreshed from the release. SAFE TO REGENERATE.
#
#   (B) TIER              - top_confidence (T1-T5). NOT produced by the pipeline.
#                           Tier assignment is a downstream gene-prioritization step
#                           the locus-level build script that
#                           no release carries. Joined from a preserved reference.
#
#   (C) GENE-LEVEL        - trans_*, ct_*, context, n_contexts, TWAS/MR/cTWAS flags.
#                           Also from the downstream script. Carried forward by gene.
#
# Regenerating this table from the release alone would silently blank every (B) and
# (C) column and drop any gene whose only support is gene-level rather than localized
# (EPHX2 is the worked example). So (A) refreshes, (B) and (C) are joined, and every
# row records where its evidence came from.
#
# A FULL refresh of (C) requires re-running the downstream prioritization script on
# the new release. This script does not attempt that; it reports the gap instead.
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(openxlsx)})

args    <- commandArgs(trailingOnly = TRUE)
release <- if (length(args) >= 1) args[1] else stop("usage: build_shiny_data.R <release_dir> [out.csv]")
outfile <- if (length(args) >= 2) args[2] else "data_refreshed.csv"

here      <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
prev_data <- file.path(here, "data.csv")

stopifnot(dir.exists(release), file.exists(prev_data))

# ---- (A) locus / variant evidence, from the release ------------------------
xl <- list.files(release, pattern = "^unified_AD_loci_xQTL_summary.*\\.xlsx$", full.names = TRUE)
if (!length(xl)) stop("no unified_AD_loci_xQTL_summary*.xlsx in ", release)
message("reading ", basename(xl[1]))
# The sheet carries a merged title banner above the header, and which row that
# lands on has varied between builds. Probe for the row that actually holds the
# column names rather than assuming row 1.
read_with_header <- function(path, probe = "Variant.ID", max_skip = 5) {
  for (sr in seq_len(max_skip)) {
    d <- tryCatch(read.xlsx(path, sheet = 1, startRow = sr), error = function(e) NULL)
    if (!is.null(d) && probe %in% names(d)) {
      message("  header found on row ", sr)
      return(as.data.table(d))
    }
  }
  stop("could not locate a header row containing '", probe, "' in ", basename(path))
}
new <- read_with_header(xl[1])

# release column -> explorer column
# The xlsx is a PRESENTATION export: a title banner, a row of merged group
# headers, then display labels on row 3 - not machine column names. Map those
# labels to the names the explorer expects.
map <- c(
  "ADlocus"                             = "ADlocus",
  "Variant.ID"                          = "variant_ID",
  "Rsid"                                = "rsid",
  "Chr"                                 = "chr",
  "Pos"                                 = "pos",
  "Effect.allele"                       = "effect_allele",
  "cV2F.score"                          = "cv2f_score",
  "maximum.inclusion.score"             = "max_inclusion",
  "maximum.inclusion.score.method"      = "max_inclusion_method",
  "cV2F.rank"                           = "cv2f_rank",
  "Maximum.zscore"                      = "max_zscore",
  "Min.p-value"                         = "min_pval",
  "Significance"                        = "significance",
  "log10(pval)"                         = "log10pval",
  "GWAS.associated.to.this.variant"     = "gwas_assoc",
  "xQTL.target.gene"                    = "gene"
)
missing <- setdiff(names(map), names(new))
if (length(missing)) stop("release xlsx is missing expected columns: ", paste(missing, collapse = ", "))

A <- new[, names(map), with = FALSE]
setnames(A, names(map), unname(map))
message(sprintf("  (A) refreshed: %d rows, %d loci", nrow(A), uniqueN(A$ADlocus)))

# ---- (B) tier assignment ----------------------------------------------------
# Tiers come from the release itself. gene_prio_utils.R assigns top_confidence
# (T1-T6) per row of res_AD_variants_xQTL.csv.gz during the locus-level build,
# and that table also carries gene_name, so no external gene reference or
# previously published tier file is needed. A gene is reported at its strongest
# tier.
tier_file <- file.path(release, "res_AD_variants_xQTL.csv.gz")
if (!file.exists(tier_file))
  stop("release is missing res_AD_variants_xQTL.csv.gz: ", tier_file)
tt <- fread(tier_file, select = c("gene_name", "top_confidence"))
tt <- tt[!is.na(top_confidence) & top_confidence != "" &
         !is.na(gene_name) & gene_name != ""]
tiers <- tt[, .(top_confidence = min(top_confidence)), by = .(gene = gene_name)]
tier_src <- basename(tier_file); tier_key <- "gene"
message(sprintf("  (B) tiers from %s, keyed on %s: %d genes (%s)",
                tier_src, tier_key, nrow(tiers),
                paste(sprintf("%s=%d", names(table(tiers$top_confidence)),
                              table(tiers$top_confidence)), collapse = " ")))


# ---- (C) gene-level columns, carried forward -------------------------------
prev <- fread(prev_data)
carry <- grep("^(trans_|ct_)|^(gene_id|context|n_contexts|ordered_contexts|dist_tss|dist_tes|max_twas_z|max_twas_ctx|twas_sig|mr_sig|ctwas_sig|has_trans|xqtl_max_inclusion|variant_rank)$",
              names(prev), value = TRUE)
C <- unique(prev[, c("variant_ID", carry), with = FALSE], by = "variant_ID")
message(sprintf("  (C) carried forward: %d columns for %d variants", length(carry), nrow(C)))

# ---- assemble ---------------------------------------------------------------
out <- merge(A, C, by = "variant_ID", all.x = TRUE)
out <- merge(out, tiers, by = tier_key, all.x = TRUE)

out[, evidence_locus := "release"]
out[, evidence_gene  := fifelse(is.na(context) & is.na(has_trans), "missing", "202605")]
## --- why there is NO T6 backfill here ---------------------------------
## An earlier version of this script promoted genes to T6 when twas_sig /
## mr_sig / ctwas_sig were TRUE but no tier had been assigned. That was
## wrong. Those flags are carried forward from the PREVIOUS data.csv and
## tagged evidence_gene = "202605"; they are not produced by the release
## being built. Checked against out_20260917_fix: of the 36 genes that
## backfill promoted, 0 appear in this release's own
## res_AD_XWAS_MR_filtered_TWAS_sig_overlapADloci.csv.gz or
## res_AD_cTWAS_pip075_overlapADloci.csv.gz, while all 36 are present in
## res_adxub. So the tier chain did see them and correctly found no
## gene-level evidence to award T6 on.
## Genes with this release's TWAS/MR/cTWAS evidence (251 of them) are all
## already tiered, so the pipeline is internally consistent and the
## explorer now matches the workbook. Do not reinstate this without first
## confirming the evidence comes from the release being built.

## --- drop trans-CCA ---------------------------------------------------
## The CCA trans-program columns are carried through from the release
## directory but are not shown anywhere in the explorer; drop them so the
## delivered table matches what the app actually uses.
.cca <- intersect(c("trans_cca_n","trans_cca_programs","trans_cca_contexts"), names(out))
if (length(.cca)) {
  out[, (.cca) := NULL]
  message("[cca] dropped columns: ", paste(.cca, collapse=", "))
}
out[, evidence_tier  := fifelse(is.na(top_confidence), "untiered", tier_src)]


# ---- report -----------------------------------------------------------------
prev_genes <- unique(na.omit(prev$gene)); new_genes <- unique(na.omit(out$gene))
message("\n--- refresh report ---")
message(sprintf("rows            : %d  (was %d)", nrow(out), nrow(prev)))
message(sprintf("loci            : %d  (was %d)", uniqueN(out$ADlocus), uniqueN(prev$ADlocus)))
message(sprintf("genes           : %d  (was %d)", length(new_genes), length(prev_genes)))
message(sprintf("genes lost      : %d  %s", length(setdiff(prev_genes, new_genes)),
                paste(utils::head(setdiff(prev_genes, new_genes), 8), collapse = ", ")))
message(sprintf("genes new       : %d", length(setdiff(new_genes, prev_genes))))
message(sprintf("tiered rows     : %d / %d", sum(out$evidence_tier != "untiered"), nrow(out)))
message(sprintf("gene-level gaps : %d rows have no 202605 gene-level evidence", sum(out$evidence_gene == "missing")))
message("\nNOTE: rows flagged evidence_gene='missing' belong to loci or variants that did")
message("not exist in the 202605 build. Their trans/cell-type columns stay empty until the")
message("downstream prioritization script is re-run on this release.")

## Tier labels: the carried-forward gene-level columns still spell the
## confidence level "CL#"; the pipeline now emits "T#". Normalise so the
## Tier column and the contexts string agree.
for (.c in intersect(c("ordered_contexts","xQTL_effects"), names(out)))
  out[, (.c) := gsub("\\(CL([0-9])", "(T\\1", get(.c))]

fwrite(out, outfile)

# ---- build provenance ---------------------------------------------------
# Written beside the data so the explorer can state on screen exactly which
# release and tier source it is serving. Plain key,value CSV, no extra deps.
## NOTE: data.table(key=) is a reserved arg (sets the table key), not a column
## named "key" - build as a data.frame first, then convert.
prov <- as.data.table(data.frame(stringsAsFactors = FALSE,
  key = c("built_at","built_by","release","release_path","source_xlsx",
          "tier_source","n_rows","n_loci","n_genes","builder"),
  value = as.character(c(
    format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"),
    Sys.info()[["user"]],
    basename(release),
    normalizePath(release, mustWork = FALSE),
    basename(xl[1]),
    tier_src,  ## the tier file actually used for tiering, not the published reference
    nrow(out), uniqueN(out$ADlocus), length(unique(na.omit(out$gene))),
    "build_shiny_data.R"))))
provfile <- file.path(dirname(outfile), "build_provenance.csv")
fwrite(prov, provfile)
message("wrote ", provfile)

message("\nwrote ", outfile)
