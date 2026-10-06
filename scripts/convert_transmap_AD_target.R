#!/usr/bin/env Rscript
## convert_transmap_AD_target.R -- turn the transmap AD-target colocalization table
## (one row per colocalization set, variants and outcomes packed in columns) into the
## ColocBoost export layout the build reads (one row per variant).
##
## Rules:
## - Every set keeps its AD (focal) signal, so it can define or support an AD locus.
## - cis outcomes become gene evidence, renamed to the registry context of their tissue.
## - trans outcomes are not gene evidence here; they already reach the table as transmap
##   annotation. A trans-only set is written with the AD outcome alone.
## - cis gene x context pairs already present in an existing ColocBoost export are dropped,
##   so no colocalization is counted twice.
## - sets repeated across overlapping regions are written once.
##
## Rscript convert_transmap_AD_target.R <in.tsv> <out.bed> <existing_export.bed> [...]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) >= 2)
x <- fread(a[1]); out <- a[2]; prev <- a[-(1:2)]
ctx_of <- c(AC='AC_DeJager_eQTL', DLPFC='DLPFC_DeJager_eQTL', PCC='PCC_DeJager_eQTL',
            Exc='Exc_DeJager_eQTL', Inh='Inh_DeJager_eQTL', Oli='Oli_DeJager_eQTL',
            OPC='OPC_DeJager_eQTL', Ast='Ast_DeJager_eQTL', Mic='Mic_DeJager_eQTL')
seen <- character(0)
for (f in prev) { p <- fread(f); for (i in seq_len(nrow(p))) { e <- trimws(strsplit(p$event_ID[i], ';')[[1]]); c <- sub('_ENSG.*$', '', e); seen <- c(seen, paste(p$region_ID[i], c[!grepl('^AD_', c)])) } }
seen <- unique(seen)
rows <- list(); n_drop <- 0L; keys <- character(0)
for (r in seq_len(nrow(x))) {
  v <- trimws(strsplit(x$colocalized_variables[r], ';')[[1]]); w <- trimws(strsplit(x$colocalized_variables_vcp[r], ';')[[1]])
  stopifnot(length(v) == length(w))
  ev <- trimws(strsplit(x$colocalized_outcomes[r], ';')[[1]]); ev <- ev[ev != '']
  ad <- ev[grepl('^AD_', ev)]; if (!length(ad)) ad <- x$focal_outcome[r]
  cis <- ev[grepl('_cis_', ev)]
  key <- paste(paste(sort(v), collapse=','), paste(sort(ev), collapse=','))
  if (key %in% keys) next; keys <- c(keys, key)
  tis <- sub('_cis_.*$', '', cis); g <- regmatches(cis, regexpr('ENSG[0-9]+', cis)); cx <- unname(ctx_of[tis])
  stopifnot(!anyNA(cx))
  ok <- !(paste(g, cx) %in% seen); n_drop <- n_drop + sum(!ok); g <- g[ok]; cx <- cx[ok]
  base <- data.table(variant_ID=v, vcp=as.numeric(w), cos_npc=x$cos_npc[r], min_npc_outcome=x$min_npc_outcome[r])
  sid <- gsub('[^A-Za-z0-9]', '_', x$cos_id[r])
  if (length(g)) {
    for (gg in unique(g)) { e <- c(paste0(unique(cx[g == gg]), '_', gg), ad)
      rows[[length(rows)+1]] <- base[, `:=`(region_ID=gg, event_ID=paste(e, collapse='; '), cos_ID=paste0(gg, ':tmAD', r, '_', sid), coef=paste(rep('NA', length(e)), collapse=';'))][] ; base <- copy(base[, 1:4]) }
  } else {
    rows[[length(rows)+1]] <- base[, `:=`(region_ID=paste0('tmAD_', x$region[r]), event_ID=paste(ad, collapse='; '), cos_ID=paste0('tmAD', r, '_', sid), coef=paste(rep('NA', length(ad)), collapse=';'))][]
  }
}
y <- rbindlist(rows)
sp <- tstrsplit(y$variant_ID, ':', fixed=TRUE)
y[, `:=`(`#chr`=sp[[1]], start=as.integer(sp[[2]]), end=as.integer(sp[[2]]), a1=sp[[4]], a2=sp[[3]])]
setcolorder(y, c('#chr','start','end','a1','a2','variant_ID','region_ID','event_ID','cos_ID','vcp','cos_npc','min_npc_outcome','coef'))
fwrite(y, out, sep='\t', quote=FALSE)
message(sprintf('[tmAD] %d input sets, %d written once, %d rows, %d genes, %d duplicated cis pairs dropped, %d sets written as AD-only', nrow(x), length(keys), nrow(y), uniqueN(y[!grepl('^tmAD_', region_ID)]$region_ID), n_drop, uniqueN(y[grepl('^tmAD_', region_ID)]$cos_ID)))
