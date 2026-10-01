## ---- (T) trans gene rows ------------------------------------------------
## The gene column is filled from cis results only, so a locus whose only
## evidence is trans shows no gene at all. This reads the paired trans
## fine-mapping table from the release and adds one row per locus-gene,
## tagged cis_trans = "trans" so nothing cis-derived is attached to it.
## Sourced from build_shiny_data.R; expects `out`, `release` and `here`.

out[, cis_trans := "cis"]

.tf <- file.path(release,
  "res_all_transgene_single_context_finemapping_cs50orgreater_overlapADloci.csv.gz")

.gref <- Sys.getenv("AD_LOCI_GENE_REF", unset = "")
if (!nzchar(.gref)) {
  .rn   <- "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list"
  .cand <- c(file.path(dirname(here), "resource", "references", .rn),
             file.path(dirname(dirname(here)), "AD_loci_xQTL", "resource", "references", .rn))
  .gref <- c(.cand[file.exists(.cand)], "")[1]
}

if (!file.exists(.tf)) {
  message("  (T) trans fine-mapping table not in the release; no trans rows added")
} else {
  .td <- fread(.tf, select = c("variant_ID","gene_ID","context","ADlocus","ADlocusID","PIP"))
  ## some trans effects are not assigned to an AD locus; they have no home
  ## in a locus-keyed table, so drop them rather than invent a 189th locus.
  .td <- .td[!is.na(ADlocus) & trimws(ADlocus) != ""]
  .td[, gid := tstrsplit(gene_ID, ".", fixed = TRUE)[[1]]]
  if (nzchar(.gref) && file.exists(.gref)) {
    .gr  <- fread(.gref, header = FALSE)
    .gr[, gid := tstrsplit(V4, ".", fixed = TRUE)[[1]]]
    .map <- unique(.gr[, .(gid, sym = V5)], by = "gid")
    .td  <- merge(.td, .map, by = "gid", all.x = TRUE)
  } else {
    message("  (T) no gene reference found; trans genes stay as Ensembl IDs")
    .td[, sym := NA_character_]
  }
  .td[, gene := fifelse(is.na(sym) | sym == "", gene_ID, sym)]

  ## one row per locus-gene; the representative variant is the highest-PIP one
  .tp <- .td[order(-PIP)][, .(variant_ID    = variant_ID[1],
                              context       = context[1],
                              max_inclusion = PIP[1],
                              trans_n_contexts_gene = uniqueN(context),
                              trans_source_contexts = paste(unique(context), collapse = "|")),
                          by = .(ADlocus, ADlocusID, gene, gene_id = gene_ID)]

  ## Variant-level fields belong to the variant, so inherit them. Gene-level
  ## ones (tier, TWAS/MR/cTWAS, distances, context flags) are all computed on
  ## cis pairs and must stay empty on a trans row.
  .vcols <- intersect(c("rsid","chr","pos","effect_allele","cv2f_score","cv2f_rank",
                        "max_zscore","min_pval","significance","log10pval","gwas_assoc"),
                      names(out))
  if (length(.vcols)) {
    .vv <- unique(out[, c("variant_ID", .vcols), with = FALSE], by = "variant_ID")
    .tp <- merge(.tp, .vv, by = "variant_ID", all.x = TRUE)
  }

  .tp[, `:=`(cis_trans            = "trans",
             max_inclusion_method = "trans_finemapping",
             evidence_tier        = "untiered",
             evidence_locus       = "release",
             evidence_gene        = "release",
             has_trans            = TRUE)]
  .tp <- .tp[!is.na(gene) & gene != ""]

  message(sprintf("  (T) trans rows added: %d across %d loci, %d genes",
                  nrow(.tp), uniqueN(.tp$ADlocus), uniqueN(.tp$gene)))
  out <- rbindlist(list(out, .tp), use.names = TRUE, fill = TRUE)
  rm(.td, .tp)
}

## ---- (T2) transmap pairs ------------------------------------------------
## Transmap is trans evidence keyed on events rather than genes. Where the
## event names an Ensembl gene the target can be placed on the genome and
## drawn on the Trans tab; glycan, metabolite and hotspot-module targets
## have no position and are left out, the same rule hotspot already follows.
.tmf <- file.path(release,
  "res_all_transmap_single_context_finemapping_cs95_overlapADloci.csv.gz")
.gr2 <- if (nzchar(.gref) && file.exists(.gref)) fread(.gref, header = FALSE) else NULL
if (!is.null(.gr2)) .gr2[, gid := tstrsplit(V4, ".", fixed = TRUE)[[1]]]

if (!file.exists(.tmf)) {
  message("  (T2) transmap table not in the release; no transmap pairs written")
} else {
  .tm <- fread(.tmf, select = c("variant_ID", "context", "event_ID", "ADlocus", "PIP"))
  .tm <- .tm[!is.na(ADlocus) & nzchar(trimws(ADlocus))]
  .tm[, gid := regmatches(event_ID, regexpr("ENSG[0-9]+", event_ID))[1], by = event_ID]
  .n_all <- nrow(.tm)
  .tm <- .tm[!is.na(gid) & grepl("^ENSG", gid)]
  if (!is.null(.gr2)) {
    .tm <- merge(.tm, unique(.gr2[, .(gid, sym = V5)], by = "gid"), by = "gid", all.x = TRUE)
  } else .tm[, sym := NA_character_]
  .tm[, gene := fifelse(is.na(sym) | sym == "", gid, sym)]
  ## one row per locus-gene; the highest-PIP variant stands in for the pair
  .tp2 <- .tm[order(-PIP)][, .(variant_ID = variant_ID[1],
                               context    = context[1],
                               max_pip    = PIP[1],
                               n_ctx      = uniqueN(context),
                               contexts   = paste(unique(context), collapse = "|")),
                           by = .(ADlocus, gene, gene_id = gid)]
  fwrite(.tp2, file.path(here, "transmap_pairs.csv"))
  message(sprintf("  (T2) transmap pairs: %d locus-gene pairs across %d loci (%d of %d rows had a gene target)",
                  nrow(.tp2), uniqueN(.tp2$ADlocus), nrow(.tm), .n_all))

  ## gene_positions.csv must cover every trans target or its arc silently
  ## vanishes. Top it up here from the reference so it stops being a hand
  ## step. Chromosomes are stored bare (21, not chr21) to match the file.
  .gpf <- file.path(here, "gene_positions.csv")
  if (file.exists(.gpf) && !is.null(.gr2)) {
    .gp   <- fread(.gpf)
    .need <- setdiff(unique(c(out[cis_trans == "trans", gene], .tp2$gene)), .gp$gene)
    if (length(.need)) {
      .add <- unique(.gr2[V5 %in% .need, .(gene = V5, chr = sub("^chr", "", V1), start = V2, end = V3)], by = "gene")
      fwrite(rbind(.gp, .add[, names(.gp), with = FALSE]), .gpf)
      message(sprintf("  (T2) gene_positions.csv topped up with %d genes (%d still unplaceable)",
                      nrow(.add), length(.need) - nrow(.add)))
    }
  }
  rm(.tm, .tp2)
}
