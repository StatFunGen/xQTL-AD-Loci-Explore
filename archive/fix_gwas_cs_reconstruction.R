## Rebuild res_all_single_gwas_finemapping_..._converged.csv.gz from the converged
## extension checkpoint, WITHOUT discarding original credible sets.
##
## Bug: resfp was built only from gwas_new (the extension output) and then merged to
## res_gw with all.x=TRUE. Any genuine original credible set that never entered
## gwas_new was silently dropped (2355 of 3355 CS variants; 368 regions -> 88),
## including APOE. Fix: re-add those variants carrying their ORIGINAL cs assignment.
##
## Usage: Rscript fix_gwas_cs_reconstruction.R <out_dir>
suppressMessages({library(data.table); library(stringr)})
args <- commandArgs(trailingOnly=TRUE)
OUT  <- if (length(args)) args[1] else stop("need out dir")
PR <- Sys.getenv("AD_LOCI_ROOT", unset = getwd())
STG  <- file.path(PR,"repro/main_text/5_AD_xQTL_genes_cis_trans/staging/gene_priorization_table")
fp   <- function(d,f) file.path(d,f)
ckpt <- fp(OUT,'all_adlocis_extended_any0.8ANDmin0.5union_withfsusie.csv.gz')
dest <- fp(OUT,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz')
stopifnot(file.exists(ckpt))

gwas_new  <- fread(ckpt)
gwas_newf <- gwas_new[source=='AD_GWAS_finemapping']
gwas_newf[, study := strsplit(locuscontext_id,'_chr')[[1]][1], by='locuscontext_id']
gwas_newf[, coverage := ifelse(str_detect(locuscontext_id,'cs95'),'cs_coverage_0.95',
                        ifelse(str_detect(locuscontext_id,'cs70'),'cs_coverage_0.7','cs_coverage_0.5'))]
gwas_newf[, cs_num := as.numeric(str_extract(locuscontext_id,'[0-9]+$'))]
resfp <- dcast(gwas_newf, study+region+variant_ID~coverage, value.var='cs_num',
               fun.aggregate=function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm=TRUE))
for (m in c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5'))
  if (!m %in% names(resfp)) resfp[, (m) := NA_real_]

mtd <- fread(fp(STG,"metadata_analysis.csv"), header=TRUE, select=1:6)
res_gw <- rbindlist(lapply(file.path(PR, mtd[Method=='AD_GWAS_finemapping']$Path), fread), fill=TRUE)
for (cc in c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5'))
  if (is.character(res_gw[[cc]])) res_gw[, (cc) := as.integer(sub('.*_','',get(cc)))]
res_gw[, study := event_ID]

resfp <- merge(resfp, res_gw[, .(variant_ID,PIP,z,study,
                 cs_coverage_0.95_original=cs_coverage_0.95,
                 cs_coverage_0.7_original =cs_coverage_0.7,
                 cs_coverage_0.5_original =cs_coverage_0.5)],
               by=c('study','variant_ID'), all.x=TRUE)

## ---- FIX: recover original credible sets dropped by the extension step ----
gw_cs <- res_gw[cs_coverage_0.95!=0 | cs_coverage_0.7!=0 | cs_coverage_0.5!=0]
miss  <- gw_cs[!resfp, on=c('study','variant_ID')]
message('[fix] CS variants in res_gw: ', nrow(gw_cs),
        ' | dropped by extension, re-added: ', nrow(miss))
if (nrow(miss)) {
  add <- miss[, .(study, region=gene_ID, variant_ID,
                  cs_coverage_0.95, cs_coverage_0.7, cs_coverage_0.5, PIP, z,
                  cs_coverage_0.95_original=cs_coverage_0.95,
                  cs_coverage_0.7_original =cs_coverage_0.7,
                  cs_coverage_0.5_original =cs_coverage_0.5)]
  resfp <- rbind(resfp, add, fill=TRUE)
}
## --------------------------------------------------------------------------

resfp[, event_ID := study]
for (k in c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5',
            'cs_coverage_0.95_original','cs_coverage_0.7_original','cs_coverage_0.5_original'))
  resfp[is.na(get(k)) | is.infinite(get(k)), (k) := 0]
resfp[, n.extension := sum(`cs_coverage_0.95`!=`cs_coverage_0.95_original`), by=.(study,region,`cs_coverage_0.95`)]
resfp[, n.extension := sum(`cs_coverage_0.7` !=`cs_coverage_0.7_original`),  by=.(study,region,`cs_coverage_0.7`)]
resfp <- resfp[!is.na(PIP)]

message('[fix] rows: ', nrow(resfp), '  regions: ', uniqueN(resfp$region),
        '  studies: ', uniqueN(resfp$study))
apoe <- merge(resfp, unique(res_gw[,.(variant_ID,chr,pos)]), by='variant_ID', all.x=TRUE)
message('[fix] APOE rows (chr19:44905791-44909967): ',
        apoe[chr=='chr19' & pos>=44905791 & pos<=44909967, .N])
fwrite(resfp, dest)
message('[fix] wrote ', dest)
