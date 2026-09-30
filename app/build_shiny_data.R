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

## The release is required. The previous data.csv is not: it only supplies the
## columns that still have no release source. Without it the build still runs,
## and reports which columns come out empty instead of refusing to start.
stopifnot(dir.exists(release))
.have_prev <- file.exists(prev_data)
if (!.have_prev)
  warning("[carry] no previous data.csv at ", prev_data,
          "; building from the release alone. Columns with no release source ",
          "will be empty.", call. = FALSE, immediate. = TRUE)

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

## Context and the ordered-context list also come from the workbook: they sit at
## the same one-row-per-variant-gene granularity as the rest of section (A).
## Carrying them forward from the previous data.csv left them without any context
## added in this release, so newly admitted genes showed no evidence at all.
## Matched by pattern so that a reworded workbook header warns instead of
## silently falling back to the stale carried value.
## gene_id joins them: the app holds 189 blanks and the workbook supplies an ID
## for them, with no row where the two hold different IDs, so this only fills gaps.
.ctx_src <- c(context          = "^Context$",
              gene_id          = "^gene\\.ID$",
              ordered_contexts = "^Ordered\\.contexts")
for (.nm in names(.ctx_src)) {
  .col <- grep(.ctx_src[[.nm]], names(new), value = TRUE)
  if (length(.col) == 1) A[, (.nm) := new[[.col]]]
  else warning("[context] no single workbook column matches ", .ctx_src[[.nm]],
               "; ", .nm, " stays carried", call. = FALSE, immediate. = TRUE)
}
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
prev <- if (.have_prev) fread(prev_data) else data.table(variant_ID = character(0), gene = character(0))
## ---- trans block: take it from the release, not from the previous data.csv ----
## The release derives the whole trans block in one pass (gene lists, context
## lists, their counts, and credible-set coverage). Carrying these forward pinned
## them to an older build: the lists disagreed with the release on about two
## thirds of variant-gene pairs, and the counts disagreed with their own lists.
## Release column names differ from the app's, so map them explicitly; anything
## the release does not provide stays carried as before.
.trans_map <- c(
  trans_genes            = "trans_genes",
  trans_contexts         = "trans_contexts",
  trans_coverage         = "trans_coverage",
  trans_n_genes          = "n_trans_genes",
  trans_n_contexts       = "n_trans_contexts",
  trans_gpQTL_genes      = "trans_genes_gpQTL",
  trans_gpQTL_contexts   = "trans_contexts_gpQTL",
  trans_gpQTL_n_genes    = "n_trans_genes_gpQTL",
  trans_pQTL_genes       = "trans_genes_pQTL",
  trans_pQTL_contexts    = "trans_contexts_pQTL",
  trans_pQTL_n_genes     = "n_trans_genes_pQTL",
  trans_snRNA_genes      = "trans_genes_snRNA",
  trans_snRNA_contexts   = "trans_contexts_snRNA",
  trans_snRNA_n_genes    = "n_trans_genes_snRNA",
  trans_hotspot_programs = "trans_hotspot_programs",
  trans_hotspot_contexts = "trans_contexts_hotspot_programs",
  trans_hotspot_n_genes  = "n_trans_hotspot_programs")
.rel_names  <- names(fread(tier_file, nrows = 1))
.trans_have <- .trans_map[.trans_map %in% .rel_names]
if (length(.trans_have) < length(.trans_map))
  warning("[trans] release is missing: ",
          paste(setdiff(.trans_map, .rel_names), collapse = ", "),
          " -- those stay carried", call. = FALSE, immediate. = TRUE)

carry <- grep("^(trans_|ct_)|^(gene_id|context|n_contexts|ordered_contexts|dist_tss|dist_tes|max_twas_z|max_twas_ctx|twas_sig|mr_sig|ctwas_sig|has_trans|xqtl_max_inclusion|variant_rank)$",
              names(prev), value = TRUE)
carry <- setdiff(carry, names(.trans_have))

## ---- gene-level columns that the release recomputes every run ----------
## Same failure mode as the trans block: these were carried forward from the
## previous data.csv, so they stayed pinned to whatever build first produced
## them. Checked against the release on 4,195 app rows, the carried copies
## disagreed on roughly two thirds of variant-gene pairs, by up to 6.5 Mb for
## the TSS/TES distances. The release recomputes all of them, so it is their
## source. ordered_contexts is deliberately NOT here: the app holds a rendered
## display string (e.g. "Inh eQTL- (T5,n=2)") that no single release column
## supplies, so it stays carried until its derivation is reproduced here.
.rel_map <- c(dist_tss           = "distance_from_tss",
              dist_tes           = "distance_from_tes",
              twas_sig           = "TWAS_signif",
              mr_sig             = "MR_signif",
              ctwas_sig          = "cTWAS_signif",
              max_twas_z         = "twas_z_gene_max",
              max_twas_ctx       = "twas_z_gene_max_context",
              variant_rank       = "variant_rank_xqtl",
              xqtl_max_inclusion = "max_variant_inclusion_probability",
              n_contexts         = "n_contexts")
.rel_names <- names(fread(tier_file, nrows = 0))
.rel_have  <- .rel_map[.rel_map %in% .rel_names]
.rel_miss  <- .rel_map[!.rel_map %in% .rel_names]
if (length(.rel_miss))
  warning("[release] not found in ", basename(tier_file),
          ", these stay carried forward and may be stale: ",
          paste(names(.rel_miss), collapse = ", "),
          call. = FALSE, immediate. = TRUE)
carry <- setdiff(carry, names(.rel_have))
carry <- setdiff(carry, "has_trans")
carry <- setdiff(carry, grep("^ct_", carry, value = TRUE))
carry <- setdiff(carry, names(A))   # section (A) wins over the carried copy
## Keyed on the variant AND the gene: these are gene-level columns, and 721
## variants carry more than one gene, so keying on the variant alone gave
## every gene at a variant the first gene's values.
C <- unique(prev[, c("variant_ID", "gene", carry), with = FALSE],
            by = c("variant_ID", "gene"))
message(sprintf("  (C) carried forward: %d columns for %d variant-gene pairs%s",
                length(carry), nrow(C),
                if (length(carry)) paste0(" -- ", paste(sort(carry), collapse = ", ")) else ""))
## Everything named above still comes from the previous data.csv and can go
## stale the way the trans block did. As of this patch that is ordered_contexts
## (a rendered display string with no single release column) and the ct_*_xQTL
## cell-type flags (the release'"'"'s celltypes column is empty, and no other
## column was confirmed to mean the same thing). Both need a decision about
## what they should be derived from before they can move to the release.

# ---- assemble ---------------------------------------------------------------
out <- merge(A, C, by = c("variant_ID", "gene"), all.x = TRUE)
out <- merge(out, tiers, by = tier_key, all.x = TRUE)

## (C2) refresh the trans block from the release, keyed on the variant AND the gene
if (length(.trans_have)) {
  .tr <- fread(tier_file, select = c("variant_ID", "gene_name", unname(.trans_have)))
  setnames(.tr, c("gene_name", unname(.trans_have)), c("gene", names(.trans_have)))
  .tr <- unique(.tr, by = c("variant_ID", "gene"))
  out <- merge(out, .tr, by = c("variant_ID", "gene"), all.x = TRUE)
  .k <- names(.trans_have)[1]
  message(sprintf("  (C2) trans block from release: %d columns, %d of %d pairs matched",
                  length(.trans_have), sum(!is.na(out[[.k]]) & out[[.k]] != ""), nrow(out)))
}

## (C3) refresh the gene-level block from the release, keyed on variant AND gene
if (length(.rel_have)) {
  .rv <- fread(tier_file, select = c("variant_ID", "gene_name", unname(.rel_have)))
  setnames(.rv, c("gene_name", unname(.rel_have)), c("gene", names(.rel_have)))
  .nonblank <- function(x) !(is.na(x) | trimws(as.character(x)) == "")
  .rv[, .nb := Reduce(`+`, lapply(.SD, function(x) as.integer(.nonblank(x)))),
      .SDcols = names(.rel_have)]
  setorderv(.rv, c("variant_ID", "gene", ".nb"), c(1L, 1L, -1L))
  .rv <- unique(.rv, by = c("variant_ID", "gene"))
  .rv[, .nb := NULL]
  out <- merge(out, .rv, by = c("variant_ID", "gene"), all.x = TRUE)
  message(sprintf("  (C3) gene-level from release: %d columns; %s",
                  length(.rel_have),
                  paste(sprintf("%s=%d", names(.rel_have),
                                vapply(names(.rel_have),
                                       function(k) sum(.nonblank(out[[k]])), integer(1))),
                        collapse = " ")))
}

out[, evidence_locus := "release"]
## (C5) cell-type flags, derived from the release rather than carried forward.
##
## Follows how the upstream code assigns cell types: it does not keep a map in
## the script, it joins contexts_metadata.csv. complete_ADlocus_level_summary.R
## merges sn-sQTL cell types against that table's context_snsQTL column, and
## build_AD_locus_table.R groups bulk monocyte, macrophage and microglia into
## one immune group. Both rules are reproduced here from the same table, so a
## context added to the config is picked up without editing this file.
##
## The context list per variant-gene pair is xQTL_contexts, the release's copy
## of gene_prio_utils.R's XQTL_contexts. The app's own `context` column holds
## at most one context per row and is not a substitute.
##
## This changes what the flags say. The carried values cannot be reproduced
## from any context list in the release: of 1,898 rows carrying ct_Exc_xQTL,
## 594 have no contexts at all and 549 more have contexts naming no excitatory
## type. They are residue from an older build, so roughly a third to a half of
## the current ticks disappear here. That is the correction, not a regression.
.cmeta <- file.path(release, "_used_contexts_metadata.csv")
if (!file.exists(.cmeta) && nzchar(Sys.getenv("AD_LOCI_CONFIG")))
  .cmeta <- file.path(Sys.getenv("AD_LOCI_CONFIG"), "contexts_metadata.csv")
if (file.exists(.cmeta)) {
  .cm <- fread(.cmeta)
  .flag_of <- function(b) fcase(
    b %in% c("bulk_macrophage_eQTL", "bulk_microglia_eQTL", "bulk_monocyte_eQTL"),
                            "ct_Bulk_Immune_xQTL",
    grepl("^Ast_", b),      "ct_Ast_xQTL",
    grepl("^Exc_", b),      "ct_Exc_xQTL",
    grepl("^Inh_", b),      "ct_Inh_xQTL",
    grepl("^Mic_", b),      "ct_Microglia_xQTL",
    grepl("^Oli_", b),      "ct_Oli_xQTL",
    grepl("^OPC_", b),      "ct_OPC_xQTL",
    grepl("^bulk_brain_", b), "ct_Brain_xQTL",
    default = NA_character_)
  .cm[, .flag := .flag_of(context_broad)]
  .cmap <- setNames(.cm$.flag, .cm$context)
  .cmap <- .cmap[!is.na(.cmap) & !is.na(names(.cmap)) & names(.cmap) != ""]
  .cx <- fread(tier_file, select = c("variant_ID", "gene_name", "xQTL_contexts"))
  setnames(.cx, "gene_name", "gene")
  .cx <- unique(.cx[!is.na(xQTL_contexts) & xQTL_contexts != ""],
                by = c("variant_ID", "gene"))
  .cstr <- .cx[out[, .(variant_ID, gene)], on = .(variant_ID, gene), xQTL_contexts]
  .ctxs <- strsplit(ifelse(is.na(.cstr), "", as.character(.cstr)), "[|;]")
  .none <- is.na(.cstr)
  .seen <- setdiff(unique(trimws(unlist(.ctxs))), "")
  .unknown <- setdiff(.seen, .cm$context)
  if (length(.unknown))
    warning("[ct] contexts absent from contexts_metadata, contributing to no flag: ",
            paste(.unknown, collapse = ", "), call. = FALSE, immediate. = TRUE)
  .nobucket <- setdiff(intersect(.seen, .cm$context), names(.cmap))
  if (length(.nobucket))
    message(sprintf("  (C5) %d contexts have no cell-type flag by design (%s)",
                    length(.nobucket), paste(utils::head(.nobucket, 6), collapse = ", ")))
  for (.b in sort(unique(unname(.cmap)))) {
    .keys <- names(.cmap)[.cmap == .b]
    .v <- vapply(.ctxs, function(v) any(trimws(v) %in% .keys), logical(1))
    .v[.none] <- NA
    set(out, j = .b, value = .v)
  }
  message(sprintf("  (C5) cell-type flags from contexts_metadata: %s",
                  paste(sprintf("%s=%d", sort(unique(unname(.cmap))),
                                vapply(sort(unique(unname(.cmap))),
                                       function(k) sum(out[[k]] %in% TRUE), integer(1))),
                        collapse = " ")))
} else {
  warning("[ct] no contexts_metadata found; cell-type flags left as carried",
          call. = FALSE, immediate. = TRUE)
}

## (C4) has_trans is the union of the trans count columns, per Q3 of the
## flagship section-5 notebook:
##   trans_cols <- c("# trans genes", "# pQTL trans genes", "# gpQTL trans genes",
##                   "# snRNA trans genes", "# trans cca programs",
##                   "# trans hotspots programs genes")
##   has_trans := rowSums(...) > 0
## Five of those six are app columns, refreshed from the release in (C2). The
## sixth, the CCA program count, has no app column of its own and is read here
## rather than added, so the drift check does not see an unreferenced column.
##
## Not derived from Method == 'trans_finemapping': that is how gene_prio_utils.R
## builds it upstream, but the exported table drops the trans rows and keeps
## only their summaries, so applying it here yields FALSE everywhere.
##
## Against the carried column this agrees on 3,789 of 3,821 rows. All 32
## disagreements run the same way -- carried TRUE where every count is zero --
## which is the carry-forward staleness, not a second definition. Rows that
## were previously blank resolve to FALSE, which is what the definition gives.
.tcc <- intersect(c("trans_n_genes", "trans_pQTL_n_genes", "trans_gpQTL_n_genes",
                    "trans_snRNA_n_genes", "trans_hotspot_n_genes"), names(out))
.tot <- rowSums(sapply(.tcc, function(k) suppressWarnings(as.numeric(out[[k]]))),
                na.rm = TRUE)
if ("n_trans_cca_programs" %in% .rel_names) {
  .cca <- fread(tier_file, select = c("variant_ID", "gene_name", "n_trans_cca_programs"))
  setnames(.cca, "gene_name", "gene")
  .cca <- .cca[, .(cca_n = suppressWarnings(max(as.numeric(n_trans_cca_programs),
                                                na.rm = TRUE))),
               by = .(variant_ID, gene)]
  .cca[!is.finite(cca_n), cca_n := 0]
  .j <- .cca[out[, .(variant_ID, gene)], on = .(variant_ID, gene), cca_n]
  .tot <- .tot + ifelse(is.na(.j), 0, .j)
} else {
  warning("[release] no n_trans_cca_programs; has_trans built from ",
          length(.tcc), " of 6 columns", call. = FALSE, immediate. = TRUE)
}
out[, has_trans := .tot > 0]
message(sprintf("  (C4) has_trans from %d trans count columns: %d TRUE, %d FALSE",
                length(.tcc) + as.integer("n_trans_cca_programs" %in% .rel_names),
                sum(out$has_trans), sum(!out$has_trans)))

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

## ---- delivery drift check -------------------------------------------------
## Every column written here should be one the app actually reads. A column that
## is delivered but never referenced is invisible to users while still appearing
## in the released tables - that is how whole loci went missing from the explorer
## while remaining in the workbook. Warn rather than fail, so a deliberate
## addition can still ship, but it has to be noticed.
## The app sources sit beside this script, not beside the output file. Deriving
## the directory from outfile meant that writing the table anywhere else left
## .src empty and skipped the check without saying so.
.app_dir <- here
.src <- c(list.files(file.path(.app_dir, "modules"), pattern = "[.]R$", full.names = TRUE),
          list.files(.app_dir, pattern = "^app[.]R$", full.names = TRUE))
if (!length(.src))
  warning("[drift] no app sources found under ", .app_dir,
          " -- delivery drift was NOT checked", call. = FALSE, immediate. = TRUE)
if (length(.src)) {
  .code  <- paste(unlist(lapply(.src, readLines, warn = FALSE)), collapse = "\n")
  .never <- names(out)[!vapply(names(out), function(cc) grepl(cc, .code, fixed = TRUE), logical(1))]
  ## provenance columns are intentionally not surfaced in the UI
  .never <- setdiff(.never, c("evidence_locus", "evidence_gene", "evidence_tier"))
  if (length(.never)) {
    warning("[drift] delivered but never read by the app: ", paste(.never, collapse = ", "),
            " -- surface them in the app or stop delivering them.", call. = FALSE, immediate. = TRUE)
  } else {
    message("[drift] every delivered column is referenced by the app")
  }
}

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
