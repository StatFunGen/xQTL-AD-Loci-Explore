
#Gene locus summary
#Input needed:
#- metadata_analysis.csv containing all exported tables paths
#- contexts_metadata.csv containing contexts (datasets) label mapping between the different analysis/exported tables
#- columns_metadata.tsv containing all columns to keep for the excel sheets, with associated metadata (column width, coloring..)
#- excel_metadata.tsv containing some meta information for the excel sheet construction
#- pattern_coloring.tsv containing all pattern to color in the excel sheet
#- gp_coordinates_clean.txt containing the glyco protein identifier to gene link
#- cv2f_score_dir: path to directory containing the .cv2f file for each chromosme
#- gene_names: path to table with gene_id gene_name link
#- ld_meta_file: path to ADSP_R4_EUR_LD ld_meta_file.tsv
#- gene_gtf: path to genes gtf
#1 metadata is optional
#- long_table_columns_selection.csv to generate a long table with selected columns from the table xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz generated in step III 
#(in this table each row is a variant-ADlocus-Method-context-gene_name information, and so facilitate querying informations ) 

#install.packages(c('openxlsx'))
## ---- project root -------------------------------------------------------
## Override with AD_LOCI_ROOT; otherwise derived from this script's location.
.args <- commandArgs(FALSE)
.file <- sub('^--file=', '', .args[grep('^--file=', .args)])
## This script lives in scripts/; the repository root is its parent directory.
SCRIPTS <- if (length(.file)) normalizePath(dirname(.file[1])) else file.path(getwd(), 'scripts')
REPO <- if (length(.file)) dirname(SCRIPTS) else getwd()
PROJECT_ROOT <- Sys.getenv('AD_LOCI_ROOT', unset = REPO)
ROOTP   <- function(...) file.path(PROJECT_ROOT, ...)
## Repository layout: config/ holds the metadata tables, gene_prio_utils.R sits
## beside this script. STAGING holds three precomputed tables that are too large
## to distribute with the code; set AD_LOCI_STAGING to wherever they are.
CONFIG <- file.path(REPO, 'config')
STAGING <- Sys.getenv('AD_LOCI_STAGING',
  unset = ROOTP('repro/main_text/5_AD_xQTL_genes_cis_trans/staging/gene_priorization_table'))
cat('PROJECT_ROOT:', PROJECT_ROOT, '
')
source(file.path(SCRIPTS, 'gene_prio_utils.R'))

## ---- preflight -----------------------------------------------------------
## Fail here with a readable list rather than part-way through a long run.
.need <- c(
  file.path(CONFIG, c('metadata_analysis.csv','contexts_metadata.csv',
                      'columns_metadata.tsv','excel_metadata.tsv','pattern_coloring.tsv')),
  file.path(STAGING, c('gwas_variants_cor0.5.csv.gz',
                       'res_APOE_interaction_summ.csv.gz',
                       'res_msex_interaction_summ.csv.gz')))
.gone <- .need[!file.exists(.need)]
if (length(.gone)) {
  stop('missing required inputs:\n  ', paste(.gone, collapse='\n  '),
       '\n\nConfig files ship with this repository. The three staging tables are\n',
       'available separately: the LD table from Synapse syn75082260, the two interaction\n       summaries on request. Point AD_LOCI_STAGING at the directory holding them.',
       call. = FALSE)
}
if (dir.exists(STAGING)) setwd(STAGING)
#INPUT FILES
metadata_analysis <- file.path(CONFIG, 'metadata_analysis.csv')
contexts_metadata <- file.path(CONFIG, 'contexts_metadata.csv')
columns_metadata <- file.path(CONFIG, 'columns_metadata.tsv')
excel_metadata <- file.path(CONFIG, 'excel_metadata.tsv')
pattern_coloring <- file.path(CONFIG, 'pattern_coloring.tsv')
gp_coordinates_clean=ROOTP('gpQTL/gp_coordinates_clean.txt')
cv2f_score_dir=ROOTP('analysis_result/cv2f/score/xqtl_feature_max_allcv2f/')
gene_names=ROOTP('resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list')
ld_meta_file=ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv')
gene_gtf=ROOTP('resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf')



out <- Sys.getenv('AD_LOCI_OUT', unset = file.path(PROJECT_ROOT, paste0('out_', format(Sys.Date(), '%Y%m%d'))))
STAMP <- sub('^out_', '', basename(out))   # version stamp used in output filenames
dir.create(out)


mtd<-fread(metadata_analysis,header = T)
contexts<-fread(contexts_metadata)

#0)prep integration####
#harmonization of contexts names
# contexts<-fread(contexts_metadata)
#merge with coloc context
# coloccont<-fread('colocname_to_contexts.tsv')
# coloccont[context=='DLPFC_Klein_gpQTL',context:='DLPFC_Klein_gpQTL_unadjusted']
# fwrite(coloccont,'colocname_to_contexts.tsv')
# contexts<-merge(contexts,coloccont,all.x = T,by='context')
# setdiff(coloccont$context,contexts$context)
# fwrite(contexts,contexts_metadata)
# #merge with Qr context name
# #on excel

#I) Prep/harmionized the exported tables for each methods/resource####
#for each, long the table at variant context, gene level and create an unique locus-event identifier


# setnames(contexts,'context_hi','context_broad')
# mtd<-merge(mtd,unique(contexts),all.x = T,by='context')
# fwrite(mtd[,.(`Data Type`,Cohort,Modality,context,context_broad,Method,Path)],'all_analysis_summary_tables_metadata.csv')
table(mtd$Method)
# AD_GWAS_finemapping           APOE interaction                      Coloc                 ColocBoost                      ctwas 
# 8                         23                          2                          3                          1 
# fSuSiE_finemapping                         LR           msex interaction  multi_context_finemapping     multi_gene_finemapping 
# 10                         35                         35                          2                         21 
# QR single_context_finemapping          trans_finemapping                       twas 
# 22                         48                         36                          1 

# 1) the finemapping tables####
#  Finemapping: singlecontext, multicontexts, fsusie, multi gene, transqtl
# - singlecontexts####
res_sc<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='single_context_finemapping']$Path),function(f)fread(f)),fill = T)
res_sc[!is.na(cat)]
res_sc[,context:=str_extract(event_ID[1], '^.+?(?:_chr|_ENSG|_gp_)')%>%
         gsub('(_chr|_ENSG|_gp_)', '', .),by='event_ID']
table(res_sc$context)
#create locus_context_id for cs95>cs70>cs50
res_sc[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_scf<-res_sc[(hit)]
res_scf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_scf[,locuscontext_id:=paste(event_ID,credibleset,sep='_')]
res_scf[,region:=gene_ID]
fwrite(res_scf,fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz'))
res_scf<-fread(fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz'))
# unique(res_scf,by='variant_ID')[alt(variant_ID)==a2]
table(res_scf$context)

mtd[Method=='single_context_finemapping',summary_file:=fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz')]

# - multicontexts####
res_mc<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='multi_context_finemapping']$Path),function(f)fread(f)),fill = T)
res_mc[lfsr=='']
#first rm thos without lfsr
res_mcf<-res_mc[!lfsr=='']

#one row per context
multi<-res_mcf[,.(context=strsplit(event_ID,';')[[1]],
                  conditional_effect=strsplit(conditional_effect,';')[[1]]|>as.numeric(),
                  lfsr=strsplit(lfsr,';')[[1]]|>as.numeric()),
               by=c('event_ID','gene_ID')]|>unique()
res_mcf<-merge(res_mcf[,-c('conditional_effect','lfsr')],multi,by=c('event_ID','gene_ID'),allow.cartesian = T)

#create locus_context_id for cs95>cs70>cs50
res_mcf[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_mcf<-res_mcf[(hit)]
res_mcf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_mcf[,locuscontext_id:=paste(event_ID,gene_ID,credibleset,sep='_')]
res_mcf[,region:=gene_ID]
max(res_mcf$lfsr)
fwrite(res_mcf,fp(out,'res_all_multi_context_finemapping_cs50orgreater_lfsr0.01.csv.gz'))

mtd[Method=='multi_context_finemapping',summary_file:=fp(out,'res_all_multi_context_finemapping_cs50orgreater_lfsr0.01.csv.gz')]

# - fsusie##### 
file.exists(file.path(PROJECT_ROOT,mtd[Method=='fSuSiE_finemapping']$Path))
res_fs<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='fSuSiE_finemapping']$Path),function(f)fread(f)),fill = T)

#create locus_context_id for cs95>cs70>cs50
res_fs[,hit:=cs_coverage_0.95!=0]
res_fs[(hit)]#everything is cs95
res_fsf<-res_fs
# res_fsf<-res_fs[(hit)]
# res_fsf[,
#         credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
#                             ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
#                                    paste0('cs50_',cs_coverage_0.5)))]
res_fsf[,locuscontext_id:=cs_id]
res_fsf[,region:=region_ID]
res_fsf[,context:=event_ID]
res_fsf[,chr:=seqid(variant_ID)]
res_fsf[,pos:=pos(variant_ID)]

setdiff(res_fsf$context,mtd$context)
res_fsf<-res_fsf[,-'grid_position']
colnames(res_fsf)
res_fsf[,top_effect_coord:=strsplit(epi_mark_effects[1],split = ';')[[1]]|>as.numeric()|>abs()|>which.max(),by='locuscontext_id']
res_fsf[,top_effect:=strsplit(epi_mark_effects[1],split = ';')[[1]][top_effect_coord]|>as.numeric(),by='locuscontext_id']
res_fsf[,n_epi_marks:=strsplit(epi_mark_effects[1],split = ';')[[1]]|>as.numeric()|>length(),by='locuscontext_id']
res_fsf[,-'grid_positions']
summary(res_fsf$n_epi_marks)
summary(res_fsf[str_detect(context,'ATAC')]$n_epi_marks)
table(unique(res_fsf,by='locuscontext_id')$context)
res_fsf[str_detect(context,'ATAC')][,.(epi_mark_names)]
fwrite(res_fsf,fp(out,'res_all_fsusie_finemapping_cs95.csv.gz'))
res_fsf<-fread(fp(out,'res_all_fsusie_finemapping_cs95.csv.gz'))
mtd[Method=='fSuSiE_finemapping',summary_file:=fp(out,'res_all_fsusie_finemapping_cs95.csv.gz')]

# res_fsf<-UnifyLoci(res_fsf,variant_col = 'variant_ID',locus_col = 'locuscontext_id',group.by = 'region_ID')

# - multi gene####
res_mg<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='multi_gene_finemapping']$Path),function(f)fread(f)[,Path:=f]),fill = T)
res_mg[lfsr=='']
#first rm thos without lfsr
res_mgf<-res_mg[!lfsr=='']
res_mgf[,region:=gene_ID]
#one row per gene
multi<-res_mgf[,.(gene_ID=strsplit(event_ID,';')[[1]]|>str_extract('ENSG[0-9]+'),
                  conditional_effect=strsplit(conditional_effect,';')[[1]]|>as.numeric(),
                  lfsr=strsplit(lfsr,';')[[1]]|>as.numeric()),
               by=c('event_ID','region')]|>unique()
res_mgf<-merge(res_mgf[,-c('conditional_effect','lfsr','gene_ID')],multi,by=c('event_ID','region'),allow.cartesian = T)

#create locus_context_id for cs95>cs70>cs50
res_mgf[,hit:=cs_coverage_0.95!=0]

res_mgf<-res_mgf[(hit)]
res_mgf[,
        credibleset:=paste0('cs95_',cs_coverage_0.95)]
#need add context
res_mgf[,context:=strsplit(basename(Path[1]),'\\.')[[1]][1],by='Path']
setdiff(res_mgf$context,contexts$context) 

res_mgf[str_detect(context,'Oli|Exc|Inh|Ast|Mic|OPC|monocyte|DLPFC|AC|PCC'),context:=str_remove(context,'^ROSMAP_')]
res_mgf[str_detect(context,'MSBB'),context:=paste0('BM_',str_extract(context,'[0-9]+'),'_MSBB_eQTL')]
setdiff(res_mgf$context,contexts$context) 

res_mgf[,locuscontext_id:=paste(context,'multigene',event_ID,credibleset,sep='_')]

fwrite(res_mgf,fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz'))
res_mgf<-fread(fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz'))

mtd[Method=='multi_gene_finemapping',summary_file:=fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz')]


# - the transQTL####

contexts<-fread(contexts_metadata)

# Finemap
res_ts<-fread(file.path(PROJECT_ROOT,mtd[Method=='trans_finemapping']$Path))

table(res_ts$resource)
cat(sort(unique(res_ts$resource)),sep = '\n')
res_ts<-merge(res_ts,contexts[,.(resource=context_trans,context)],all.x = TRUE)
res_ts[is.na(context)]#ok

table(res_ts$context)
#create locus_context_id for cs95>cs70>cs50
res_ts[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_tsf<-res_ts[(hit)]
res_tsf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_tsf[,region:=gene_ID]

res_tsf[,locuscontext_id:=paste(event_ID,region,'trans',credibleset,sep='_')]
res_tsf[,gene_ID:=str_extract(event_ID,'ENSG[0-9]+')]
res_tsf[is.na(gene_ID),.(resource,variant_ID,event_ID,region)]#for gpqtl need mapping
gpmap<-fread(gp_coordinates_clean)
gpmap[str_detect(ensembl_id,';'),ensembl_id:=str_extract(ensembl_id,'ENSG[0-9]+')]
res_tsf[is.na(gene_ID),gp_ID:=str_extract(event_ID,'gp_[0-9]+')]

res_tsf[is.na(gene_ID),gene_ID:=gpmap[gp_ID,on='ID']$ensembl_id]
res_tsf[!is.na(gp_ID)]
res_tsf<-res_tsf[,-c('gp_ID','resource')]
fwrite(res_tsf,fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz'))
res_tsf<-fread(fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz'))

mtd[Method=='trans_finemapping',summary_file:=fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz')]


#snuc####
if (nrow(mtd[Method=='trans_single_gene_snRNA'])) {
res_ts2<-mtd[Method=='trans_single_gene_snRNA',fread(file.path(PROJECT_ROOT,Path)),by='context']
res_ts2[is.na(context )]
setdiff(res_ts2$context,contexts$context)#OK

res_ts2[,gene_ID:=str_extract(event_ID,'ENSG[0-9]+')]
#create locus_context_id 
res_ts2[,hit:=cs_coverage_0.95!=0,by=.(event_ID)]

res_ts2f<-res_ts2[(hit)]
res_ts2f[,credibleset:=paste0('cs95_',cs_coverage_0.95)]

res_ts2f[,locuscontext_id:=paste(region_ID,event_ID,'transNoHF',credibleset,sep='_')]
fwrite(res_ts2f,fp(out,'res_snuc_transgene_NoHF_single_context_finemapping_cs95.csv.gz'))
res_ts2f<-fread(fp(out,'res_snuc_transgene_NoHF_single_context_finemapping_cs95.csv.gz'))

mtd[Method=='trans_single_gene_snRNA',summary_file:=fp(out,'res_snuc_transgene_NoHF_single_context_finemapping_cs95.csv.gz')]
} else cat('[skip] trans_single_gene_snRNA - no rows registered in metadata_analysis.csv\n')

#trans pQTL
if (nrow(mtd[Method=='trans_pQTL'])) {
res_tspqtl<-mtd[Method=='trans_pQTL',fread(file.path(PROJECT_ROOT,Path)),by='context']
setdiff(res_tspqtl$context,contexts$context)#OK

res_tspqtl[,gene_ID:=str_extract(event_ID,'ENSG[0-9]+')]
#create locus_context_id 
res_tspqtl[,hit:=cs_coverage_0.95!=0,by=.(event_ID)]

res_tspqtlf<-res_tspqtl[(hit)]
res_tspqtlf[,
         credibleset:=paste0('cs95_',cs_coverage_0.95)]
res_tspqtlf[,locuscontext_id:=paste(region_ID,event_ID,credibleset,sep='_')]
fwrite(res_tspqtlf,fp(out,'res_transpQTL_NoHF_single_context_finemapping_cs95.csv.gz'))
res_tspqtlf<-fread(fp(out,'res_transpQTL_NoHF_single_context_finemapping_cs95.csv.gz'))

mtd[Method=='trans_pQTL',summary_file:=fp(out,'res_transpQTL_NoHF_single_context_finemapping_cs95.csv.gz')]
} else cat('[skip] trans_pQTL - no rows registered in metadata_analysis.csv\n')

#trans gpQTL
if (nrow(mtd[Method=='trans_gpQTL'])) {
res_tsgpqtl<-mtd[Method=='trans_gpQTL',fread(file.path(PROJECT_ROOT,Path)),by='context']
setdiff(res_tsgpqtl$context,contexts$context)#OK
table(res_tsgpqtl$context)
#add gene
gpmap<-fread(gp_coordinates_clean)
gpmap[str_detect(ensembl_id,';'),ensembl_id:=str_extract(ensembl_id,'ENSG[0-9]+')]
res_tsgpqtl[,gp_ID:=str_extract(event_ID,'gp_[0-9]+')]

res_tsgpqtl[,gene_ID:=gpmap[gp_ID,on='ID']$ensembl_id]
res_tsgpqtl[!is.na(gp_ID)]
res_tsgpqtl<-res_tsgpqtl[,-c('gp_ID')]

#create locus_context_id 
res_tsgpqtl[,hit:=cs_coverage_0.95!=0,by=.(event_ID)]

res_tsgpqtlf<-res_tsgpqtl[(hit)]
res_tsgpqtlf[,
            credibleset:=paste0('cs95_',cs_coverage_0.95)]
res_tsgpqtlf[,locuscontext_id:=paste(region_ID,event_ID,credibleset,sep='_')]
fwrite(res_tsgpqtlf,fp(out,'res_transgpQTL_NoHF_single_context_finemapping_cs95.csv.gz'))
res_tspqtlf<-fread(fp(out,'res_transgpQTL_NoHF_single_context_finemapping_cs95.csv.gz'))

mtd[Method=='trans_gpQTL',summary_file:=fp(out,'res_transgpQTL_NoHF_single_context_finemapping_cs95.csv.gz')]
} else cat('[skip] trans_gpQTL - no rows registered in metadata_analysis.csv\n')

#trans cca
if (nrow(mtd[Method=='trans_cca'])) {
res_tscca<-mtd[Method=='trans_cca',fread(file.path(PROJECT_ROOT,Path)),by='context']

#create locus_context_id 
res_tscca[,hit:=cs_coverage_0.95!=0,by=.(event_ID)]

res_tsccaf<-res_tscca[(hit)]
res_tsccaf[,
           credibleset:=paste0('cs95_',cs_coverage_0.95)]
res_tsccaf[,locuscontext_id:=paste(event_ID,'transCCA',credibleset,sep='_')]
fwrite(res_tsccaf,fp(out,'res_all_transCCA_single_context_finemapping_cs95.csv.gz'))
res_tsccaf<-fread(fp(out,'res_all_transCCA_single_context_finemapping_cs95.csv.gz'))

mtd[Method=='trans_cca',summary_file:=fp(out,'res_all_transCCA_single_context_finemapping_cs95.csv.gz')]
} else cat('[skip] trans_cca - no rows registered in metadata_analysis.csv\n')

#### single-context toploci evidence sources ################################
## To add a new source: run add_evidence_source.R (registers it in
## metadata_analysis.csv), then add one entry below. No other edit is needed.
## Any Method in metadata_analysis.csv that points at per-context toploci files
## and has no dedicated loader above is picked up automatically -- registering
## rows in metadata_analysis.csv is the ONLY step needed to add a new source.
TOPLOCI_HANDLED <- c('single_context_finemapping','multi_context_finemapping',
  'multi_gene_finemapping','fSuSiE_finemapping','AD_GWAS_finemapping','trans_finemapping',
  'trans_pQTL','trans_cca','trans_single_gene_snRNA','trans_gpQTL','ColocBoost','Coloc','twas','ctwas',
  'cTWAS','LR','msex interaction','APOE interaction','QR','AD_xQTL_colocalization','KNIGHT')
TOPLOCI_SOURCES <- setdiff(
  unique(mtd[grepl('top_?loci\\.bed(\\.gz)?$', Path) & Path != '']$Method),
  TOPLOCI_HANDLED)
cat('[toploci] auto-detected sources:', paste(TOPLOCI_SOURCES, collapse=', '), '\n')

for (.m in TOPLOCI_SOURCES) {
  if (!nrow(mtd[Method == .m])) { cat('[toploci]', .m, '- no rows registered, skipped\n'); next }
  .r <- mtd[Method == .m, fread(file.path(PROJECT_ROOT, Path)), by = 'context']
  .unknown <- setdiff(.r$context, contexts$context)
  if (length(.unknown)) warning('unknown contexts for ', .m, ': ', paste(.unknown, collapse=', '))
  .r[, hit := cs_coverage_0.95 != 0, by = .(event_ID)]
  .r <- .r[(hit)]
  .r[, credibleset := paste0('cs95_', cs_coverage_0.95)]
  .r[, locuscontext_id := paste(region_ID, event_ID, credibleset, sep = '_')]
  .f <- fp(out, sprintf('res_all_%s_single_context_finemapping_cs95.csv.gz', gsub('_', '', .m)))
  fwrite(.r, .f)
  mtd[Method == .m, summary_file := .f]
  cat('[toploci]', .m, '->', basename(.f), nrow(.r), 'rows,', uniqueN(.r$context), 'contexts\n')
}
fwrite(mtd, metadata_analysis)
############################################################################



# the gwas results: consolidate gwas locus with min corr 0.5 ####
res_gw<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='AD_GWAS_finemapping']$Path),function(f)fread(f)),fill = T)
## FIX (2026-09-17): the 8-study GWAS fine-mapping files encode the credible-set index as
## text ('susie_rss_0' = no CS, 'susie_rss_1', ...). The script expects integers (0 = no CS);
## with text, `x != 0` and `x > 0` are TRUE for every row (string comparison), so every
## variant became a 'hit'/'cs95' and loci were inflated. Extract the trailing integer.
for (.cc in intersect(c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5'), names(res_gw))) {
  if (is.character(res_gw[[.cc]])) res_gw[, (.cc) := as.integer(sub('.*_', '', get(.cc)))]
}
res_gw[!is.na(cat)]
res_gw[,context:=event_ID]
res_gw[,region:=gene_ID]
res_gw[,study:=event_ID]
#
#consolidate gwas locus with min corr 0.5 res_gwf<-res_gw[(hit)]
res_gw[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]
res_gwf<-res_gw[(hit)]
res_gwf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_gwf[,region:=gene_ID]

res_gwf[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]

fwrite(res_gwf,fp(out,'res_all_single_gwas_finemapping_cs50orgreater.csv.gz'))

# 2) The non linear tables####
# the QR table####
#two pvalue: bonferonni wich is adjusted for each genes, while qvalue is adjusted by number of variants tested within a gene
#for 1
#res<-fread(ROOTP('analysis_result/marginal_significant_qtl/cis_association/KNIGHT/eQTL/Brain/QR/Knight_quantile_eQTL.cis_pairs.significant_qtl.filtered_bonferroni_BH_adjusted.tsv.gz'))  ## disabled: exploratory read, unused downstream
#res[,.(molecular_trait_object_id,variant_id,pvalue,qvalue,p_bonferroni_adj,coef_heter)][coef_heter>-0.6931472]$molecular_trait_object_id|>unique()|>length()
#res[,.(molecular_trait_object_id,variant_id,pvalue,qvalue,p_bonferroni_adj,coef_heter)]$molecular_trait_object_id|>unique()|>length()


#need add LR/QR correlation test to separate real non linear effect from linear
#Pending Anjing update
#update in analysis_result/quantile_qtl/export/pure_qr.qtl.bonf_fdr_add_xi.classification.tsv.gz
#edit mtd
mtd[Method=='QR',Path:='analysis_result/quantile_qtl/export/pure_qr.qtl.bonf_fdr_add_xi.classification.tsv.gz']
fwrite(mtd,metadata_analysis)

res_qr<-rbindlist(lapply(file.path(PROJECT_ROOT,unique(mtd[Method=='QR']$Path)),function(f){
  
  fread(f)
  
  
}))
res_qr<-res_qr[p_bonferroni_adj<0.05][,.(chrom,pos,ref,alt,molecular_trait_object_id,variant_id,maf,cis_start,cis_end,feature_tss,feature_tes,pvalue ,qvalue,coef_heter,context,n_variants,p_bonferroni_adj,xi,xi_pval,classification)]

res_qr[,gene_ID:=str_extract(molecular_trait_object_id,'ENSG[0-9]+')]
res_qr[,variant_ID:=str_replace_all(variant_id,'_',':')]
setnames(res_qr,'context','context_qr')
res_qr<-merge(res_qr,contexts[,.(context,context_qr,context_hi)],by = 'context_qr',all.x = T)
res_qr[is.na(context)]
#n genes
unique(res_qr$gene_ID)

fwrite(res_qr,fp(out,'res_qr_summ.csv.gz'))
res_qr<-fread(fp(out,'res_qr_summ.csv.gz'))

mtd[Method=='QR',summary_file:=fp(out,'res_qr_summ.csv.gz')]
mtd<-unique(mtd,by=c('Path'))
mtd[Method=='QR',Cohort:='ROSMAP, MSBB, KNIGHT']
mtd[Method=='QR',`Data Type`:='eQTL & pQTL']
mtd[Method=='QR',`Modality`:='-']

fwrite(mtd,metadata_analysis)

# the interactions tables####
#have been corrected 
# #for 1
# interactions<-fread(ROOTP('analysis_result/marginal_significant_qtl/cis_association/KNIGHT/eQTL/Brain/interaction/msex/xqtl_protocol_data.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.bed.bed.cis_pairs.significant_qtl.q_bonferroni_min_adjusted_events_qvalue.tsv.gz'))
# max(interactions$qvalue_interaction)
#for all
#msex
#msex_summ<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='msex interaction']$Path),function(f){
#  file=list.files(f,pattern='q_bonferroni_min_adjusted_events_qvalue.tsv.gz',full.names = T)
#  if(length(file)==1){
#    return(fread(file)[,Path:=str_remove(f,ROOTP(''))])
#    
#  }else{
#    return(data.table(Path=str_remove(f,ROOTP(''))))
#  }
#  
#  
#}),fill = T)
#msex_summ[,variant_id:=str_replace_all(variant_id,'_',':')]
#msex_summ[,gene_ID:=str_extract(molecular_trait_id,'ENSG[0-9]+')]
#msex_summ[,variant_ID:=str_replace_all(variant_id,'_',':')]
#msex_summ<-merge(msex_summ,mtd[,.(Path,context,context_broad)],by = 'Path')
#msex_summ<-msex_summ[!is.na(variant_id)]
#fwrite(msex_summ,fp(out,'res_msex_interaction_summ.csv.gz'))
msex_summ<-fread(file.path(STAGING,'res_msex_interaction_summ.csv.gz'))  ## inputs absent: precomputed summary shipped in staging
fwrite(msex_summ,fp(out,'res_msex_interaction_summ.csv.gz'))  ## copy into out/ for provenance

mtd[Method=='msex interaction',summary_file:=fp(out,'res_msex_interaction_summ.csv.gz')]


#define locus: lead SNP and variants signif in its block with r>0.1
#NOT RUN for now, if needed later
#for one block
# block<-'chr17_42087601_45383525'
# ldblock<-load_LD_matrix(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'),
#                         region = data.frame(chr=seqid(block),
#                                             start=start(block),
#                                             end=end(block)))$combined_LD_matrix
# msex_summf<-msex_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
# 
# table(msex_summf$context)
# loci<-getLocusID(msex_summf[context=='DLPFC_Bennett_pQTL'],ldblock,order.by = 'pvalue_msex_interaction',group.by = 'gene_ID')
# table(loci$locus_id)
# 
# 
# #for all
# library(parallel)
# msex_summ<-msex_summ[!is.na(variant_id)]
# ldmeta<-fread(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'))
# blocks<-bed_inter(unique(msex_summ[,.(seqid(variant_id),pos(variant_id)-1,pos(variant_id),variant_id)]),ldmeta,select = 8)[[1]]|>unique()
# ldmetaf<-ldmeta[path%in%blocks]
# ldmetaf[,block:=paste(`#chrom`,start,end,sep='_')]
# msex_summ<-rbindlist(mclapply(unique(ldmetaf$block),function(block){
#   message(block)
#   ldblock<-load_LD_matrix(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'),
#                           region = data.frame(chr=seqid(block),
#                                               start=start(block),
#                                               end=end(block)))$combined_LD_matrix
#   summf<-msex_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
#   
#   if(nrow(summf)>0){
#     summf<-getLocusID(summf,ldblock,order.by = 'pvalue_msex_interaction',group.by = c('gene_ID','context'))
#     
#   }else{
#     summf<-data.table()
#   }
# 
#   return(summf)
# },mc.cores = 1),fill = T)

#	APOE interaction
#apoe_summ<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='APOE interaction']$Path),function(f){
#  file=list.files(f,pattern='q_bonferroni_min_adjusted_events_qvalue.tsv.gz',full.names = T)
#  if(length(file)==1){
#    return(fread(file)[,Path:=str_remove(f,ROOTP(''))])
#    
#  }else{
#    return(data.table(Path=str_remove(f,ROOTP(''))))
#  }
#  
#  
#}),fill = T)
#apoe_summ
#apoe_summ<-merge(apoe_summ,mtd[,.(Path,context,context_broad)],by = 'Path')
#apoe_summ[,gene_ID:=str_extract(molecular_trait_id,'ENSG[0-9]+')]
#apoe_summ[,variant_ID:=str_replace_all(variant_id,'_',':')]
#
#apoe_summ<-apoe_summ[!is.na(variant_id)]
#fwrite(apoe_summ,fp(out,'res_APOE_interaction_summ.csv.gz'))
apoe_summ<-fread(file.path(STAGING,'res_APOE_interaction_summ.csv.gz'))  ## inputs absent: precomputed summary shipped in staging
fwrite(apoe_summ,fp(out,'res_APOE_interaction_summ.csv.gz'))  ## copy into out/ for provenance
## (removed duplicate re-read of res_APOE_interaction_summ.csv.gz)
mtd[Method=='APOE interaction',summary_file:=fp(out,'res_APOE_interaction_summ.csv.gz')]

#define locus: lead SNP and variants signif in its block with r>0.1
#NOT RUN for now, if needed later
# ldmeta<-fread(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'))
# blocks<-bed_inter(unique(apoe_summ[,.(seqid(variant_id),pos(variant_id)-1,pos(variant_id),variant_id)]),ldmeta,select = 8)[[1]]|>unique()
# ldmetaf<-ldmeta[path%in%blocks]
# ldmetaf[,block:=paste(`#chrom`,start,end,sep='_')]
# apoe_summ<-rbindlist(mclapply(unique(ldmetaf$block),function(block){
#   message(block)
#   ldblock<-load_LD_matrix(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'),
#                           region = data.frame(chr=seqid(block),
#                                               start=start(block),
#                                               end=end(block)))$combined_LD_matrix
#   summf<-apoe_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
#   if(nrow(summf)>0){
#     summf<-getLocusID(summf,ldblock,order.by = 'pvalue_msex_interaction',group.by = c('gene_ID','context'))
#     
#   }else{
#     summf<-data.table()
#   }
#   
#   return(summf)
# },mc.cores = 2))



#the sn_sQTL ####
ressqtl<-mtd[Method=='LR'&`Data Type`=='sQTL',fread(file.path(PROJECT_ROOT,Path),col.names = c('splice_site', 'variant_ID', 'pval', 'beta', 'se', 'FDR', 'if_significant', 'spliced_genes', 'gene_type', 'splice_type', 'cell_type'))]

ressqtl
table(ressqtl$splice_type)
table(ressqtl$cell_type)

ressqtl[,gene_name:=strsplit(spliced_genes,',')[[1]][1],by=.(spliced_genes)]
unique(ressqtl$gene_name)
length(unique(ressqtl[is.na(gene_name)]$splice_site))
length(unique(ressqtl$splice_site)) #6895/33512 without gene identified for the splice site
#rm them
ressqtlf<-ressqtl[!is.na(gene_name)]

#get gene id
trans<-fread(gene_names)
ressqtlf<-merge(ressqtlf,unique(trans[,.(gene_ID=gene_id,gene_name)]),all.x = TRUE)
ressqtlf[is.na(gene_ID)]$gene_name|>unique() #121/ 
ressqtlf$gene_name|>unique()|>length() #121/ 7238 genes without matching gene name on our reference

ressqtlf<-ressqtlf[!is.na(gene_ID)]

#add harmonized context name
contexts<-fread(contexts_metadata)

ressqtlf<-merge(ressqtlf,unique(contexts[,.(cell_type=context_snsQTL ,context,context_hi)]),by='cell_type')

ressqtlf<-ressqtlf[,-'cell_type']

fwrite(ressqtlf,fp(out,'res_sn_sqtl_summ.csv.gz'))
ressqtlf<-fread(fp(out,'res_sn_sqtl_summ.csv.gz'))

mtd[Method=='LR'&`Data Type`=='sQTL',summary_file:=fp(out,'res_sn_sqtl_summ.csv.gz')]

fwrite(mtd,metadata_analysis)



# 3)the colocs##### 
table(mtd$Method)
# - xQTLxADs####
res_c<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='ColocBoost'][Modality=='AD_xQTL_colocalization']$Path),
                        function(f)fread(f)),fill=TRUE)
res_c[,gene_ID:=region_ID]
res_c[,gwas_source:=str_extract(event_ID,'AD_[A-Za-z0-9_]+$')]
res_c$gwas_source|>table()
res_c[gwas_source=='AD_Bellenguez',gwas_source:='AD_Bellenguez_2022']

#transform in long format with column gene_id, study, variant_id, vcp,
res_c2<-merge(res_c[,.(
  variant_ID,
  vcp,
  cos_npc,
  min_npc_outcome),
  by=.(gene_ID,cos_ID,gwas_source)],
  res_c[,.(
    event_ID=strsplit(event_ID,'\\; ')[[1]],
    coef=as.numeric(strsplit(coef,'\\;')[[1]])),
    by=.(gene_ID,cos_ID,gwas_source)],
  allow.cartesian=TRUE)
res_c2[event_ID=='AD_Bellenguez',event_ID:='AD_Bellenguez_2022']

res_c2[,top_variant:=vcp==max(vcp),by=.(gene_ID,cos_ID,gwas_source)]
unique(res_c2[(top_variant)],by=c('cos_ID','gene_ID'))
#res_c<-merge(res_c,res_cRu[,.(cos_id,gene_id=gene,top_variable,purity)],by=c('gene_id','cos_id'))
res_c2[,locuscontext_id:=paste(cos_ID,gwas_source,gene_ID,sep='_')]

res_c2[,n.variant:=length(unique(variant_ID)),by=.(locuscontext_id)]
#add the matching context name with fp
res_c2[,context_coloc:=str_remove(event_ID,gene_ID[1])|>str_remove('(_|:)$')|>str_remove('adjusted_gp_[0-9]+|P[0-9]+')|>str_remove('_\\|[A-Z0-9]+$')|>str_remove('chr[0-9]+__[A-Z0-9]+')|>str_remove('_chr[0-9]+:[0-9]+:[0-9]+:clu_[0-9]+_[+-]:[A-Z]+')|>str_remove('chr[0-9]+__')|>str_remove('_\\|')|>str_remove('_$'),by='gene_ID']
unique(res_c2$context_coloc)|>sort()|>cat(sep = '\n')
setdiff(res_c2$context_coloc,contexts$context_coloc)#OK

res_c<-merge(res_c2[,-c('context')],unique(contexts[context_coloc!=''&!str_detect(context,'_u_|_p_')][,.(context_coloc,context)]),by='context_coloc',all.x = T)

res_c[,chr:=seqid(variant_ID)]
res_c[,pos:=pos(variant_ID)]
res_c[is.na(context)]#ok
table(res_c[is.na(context)]$context_coloc)


fwrite(res_c,fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))


# - ADs #####
res_cad<-fread(file.path(PROJECT_ROOT,mtd[Method=='ColocBoost'][Modality=='AD_meta_colocalization']$Path))
#transform in long format with column gene_id, study, variant_id, vcp,
res_cad2<-merge(res_cad[,.(
  variant_ID,
  vcp,
  cos_npc,
  min_npc_outcome
),
by=.(cos_ID)],
res_cad[,.(
  event_ID=strsplit(event_ID,'\\; ')[[1]],
  coef=as.numeric(strsplit(coef,'\\;')[[1]])),
  
  by=.(cos_ID)],
allow.cartesian=TRUE)
table(res_cad2$event_ID)
setdiff(unique(res_cad2$event_ID),contexts$context)
res_cad2[event_ID=='AD_Bellenguez',event_ID:='AD_Bellenguez_2022']

res_cad2[event_ID=='AD_Bellenguez_EADB',event_ID:='AD_Bellenguez_EADB_2022']
res_cad2[,context:=event_ID]
res_cad2[,gwas_source:=event_ID]

res_cad2[,top_variant:=vcp==max(vcp),by=.(cos_ID,gwas_source)]
res_cad2[,locuscontext_id:=cos_ID]


fwrite(res_cad2,fp(out,'res_coloc_meta_AD.csv.gz'))

# - fsusie Coloc####
res_cfs<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='Coloc']$Path),function(f)fread(f)),fill=TRUE,use.names = T)

res_cfs[event_ID=='ROSMAP_mQTL',event_ID:='ROSMAP_DLPFC_mQTL']
res_cfs[event_ID=='ROSMAP_haQTL',event_ID:='ROSMAP_DLPFC_haQTL']
res_cfs[,context:=event_ID]
unique(res_cfs$context)
unique(res_cfs$ge)

unique(res_cfs$CoS)|>length()
res_cfs[is.na(CoS)]
setdiff(res_cfs$context,contexts$context)
res_cfs[is.na(context)]
res_cfs[,gwas_source:=AD]
res_cfs$gwas_source|>table()
#keep only those in CoS
res_cfsf<-res_cfs[!is.na(CoS)]
setnames(res_cfsf,'CoS','cos_ID')
#setnames(res_cfsf,'variant_id','variant_ID')

res_cfsf[,top_variant:=SNP_PPH4==max(SNP_PPH4),by=.(cos_ID,gwas_source)]
unique(res_cfsf[(top_variant)],by=c('cos_ID','gwas_source'))
#res_c<-merge(res_c,res_cRu[,.(cos_id,gene_id=gene,top_variable,purity)],by=c('gene_id','cos_id'))
res_cfsf[,locuscontext_id:=paste(TADB_region,gwas_source,context,coloc_index,sep='_')]

res_cfsf[,n.variant:=length(unique(variant_ID)),by=.(locuscontext_id)]
res_cfsf[!str_detect(variant_ID,'^chr'),variant_ID:=paste0('chr',variant_ID)]

res_cfsf[,chr:=seqid(variant_ID)]

fwrite(res_cfsf,fp(out,'res_coloc_AD_epiQTL.csv.gz'))


#  1-3) AD locis Extension/harmonization with FP and ADmeta coloc ####

update_gwas_cs_harmo=FALSE
if(update_gwas_cs_harmo|!file.exists(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))){
  #export table: RUN export_mr.R
  library(pecotmr)
  
  bed_inter<- function(a, b, opt1="-wa",
                       opt2="-wb",out_dir=".",
                       select=NULL, col.names=NULL){
    require(data.table)
    l<-list(a,b)
    #order chr and pos
    l<-lapply(l, function(x){
      chr_cols<-colnames(x)[1:2]
      setorderv(x,chr_cols)
    })
    
    files_to_rm<-c(FALSE,FALSE)
    file_paths<-sapply(1:2, function(i){
      x<-l[[i]]
      if(is.data.frame(x)){
        files_to_rm[i]<-TRUE
        file_path<-file.path(out_dir,paste0("temp",i,".bed"))
        fwrite(x,file_path,sep="\t",col.names = FALSE,scipen = 999)
        
      }else{
        file_path<-x
      }
      return(file_path)
    })
    
    out_file<-file.path(out_dir,"temp_inter.bed")
    cmd<-paste("bedtools intersect -a",file_paths[1],"-b",file_paths[2], opt1, opt2,">",out_file)
    message("run in shell : ",cmd)
    system(cmd)
    message("done.")
    
    if(!is.null(col.names)){
      dt<-fread(out_file,select = select,col.names = col.names)
      file.remove(out_file)
      file.remove(file_paths[files_to_rm])
      return(dt)
    }else{
      dt<-fread(out_file,select = select)
      file.remove(out_file)
      file.remove(file_paths[files_to_rm])
      return(dt)
    }
  }
  
  
  #if r > 0.8, AND the union of the two CS have min abs(r) > 0.5
  #first get long table variant1 variant2 r for abs(r)>0.5
  
  res_gwf <- fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater.csv.gz'))
  res_gwf[,source:='AD_GWAS_finemapping']
  #bind with ADxQTL, ADAD and fsusie colocs
  res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))
  res_c[,chr:=seqid(variant_ID)]
  res_c[,pos:=pos(variant_ID)]
  res_cad<-fread(fp(out,'res_coloc_meta_AD.csv.gz'))
  res_cad[,chr:=seqid(variant_ID)]
  res_cad[,pos:=pos(variant_ID)]
  
  res_cfsf<-fread(fp(out,'res_coloc_AD_epiQTL.csv.gz'))
  res_cfsf[,.(variant_ID,cos_ID,context)]
  res_cfsf[,chr:=seqid(variant_ID)]
  res_cfsf[,pos:=pos(variant_ID)]
  
  res_cosad<-rbindlist(list(res_c[str_detect(context,'^AD')][,source:='AD_xQTL_colocalization'],
                            res_cad[,source:='AD_meta_colocalization'],
                            res_cfsf[,source:='AD_epiQTL_colocalization']),fill=TRUE)
  all_variants<-Reduce(union,list(res_gwf$variant_ID,res_c$variant_ID,res_cad$variant_ID,res_cfsf$variant_ID))
  sum(duplicated(seqid(all_variants))&duplicated(pos(all_variants)))#ok
  
  #for 1 block
  # b1<-res_gwf$region[1]
  # ldblock<-load_LD_matrix(ROOTP('resource/ADSP_R4_EUR_LD/ld_meta_file.tsv'),
  #                         region = data.frame(chr=seqid(b1),
  #                                             start=start(b1),
  #                                             end=end(b1)))$combined_LD_matrix
  # variants<-intersect(rownames(ldblock),res_gwf$variant_id)
  # variants_dt<-data.table(variant_id=variants)
  # variants_dt
  # variants_cors<-variants_dt[,{
  #   vars_cor<-rownames(ldblock)[abs(ldblock[,variant_id])>0.5]
  #   
  #   list(variant2=vars_cor,
  #        r=ldblock[vars_cor,variant_id])
  # },by=.(variant_id)]
  # variants_cors
  
  #for all block
  ldmeta<-fread(ld_meta_file)
  
  regions=bed_inter(data.table(chr=seqid(all_variants),
                               start=pos(all_variants)-1,
                               end=pos(all_variants),
                               variant_ID=all_variants)[!is.na(start)],
                    ldmeta[,.(`#chrom`,start,end,block=paste(`#chrom`,start,end,sep='_'))])[[8]]|>unique()
  
  
## LD-based variant correlations ------------------------------------------
## load_LD_matrix() is absent from the installed pecotmr, so the compute block
## (load_LD_matrix over each region -> pairwise |r| > 0.5) is disabled and the
## precomputed correlation table shipped in staging is used instead.
variants_cors<-fread(file.path(STAGING,'gwas_variants_cor0.5.csv.gz'))
fwrite(variants_cors,fp(out,'gwas_variants_cor0.5.csv.gz'))  ## copy into out/ for provenance
cat('[LD] precomputed variant correlations:', nrow(variants_cors), 'pairs,',
    uniqueN(c(variants_cors$variant_ID, variants_cors$variant_ID2)), 'variants\n')
  
  #if r > 0.8, AND the union of the two CS have min abs(r) > 0.5
  #per region, iterate in all CS of a GWAS finding if correlation with other GWAS CS meet this condition, if it meet,
  #extend the CS 
  variants_corsf<-variants_cors[variant_ID%in%all_variants&variant_ID2%in%all_variants]
  #iteration 1
  i<-1
  gwas_new<-rbindlist(lapply(regions,function(r){
    message(r)
    res_gwff<-res_gwf[region==r]
    #get long table of the CS
    res_gwff_cs<-melt(res_gwff,measure.vars=c('cs_coverage_0.95',
                                              'cs_coverage_0.7',
                                              'cs_coverage_0.5'),
                      variable.name='coverage',value.name = 'cs_num')
    res_gwff_cs<-res_gwff_cs[cs_num!=0]
    
    #add the coloc sets
    res_cf<-res_cosad[chr==seqid(r)&pos>=start(r)&pos<=end(r)]
    
    res_gwff_cs<-rbind(res_gwff_cs,res_cf,fill=TRUE)
    
    #iterate per CS
    
    res_gwff_cs_ext<-rbindlist(lapply(unique(res_gwff_cs$locuscontext_id),function(l){
      message(l)
      res_gwff_cs_study<-res_gwff_cs[locuscontext_id==l]
      res_gwff_cs_others<-res_gwff_cs[locuscontext_id!=l]
      res_gwff_cs_study_ext<-res_gwff_cs_study[,{
        vars_cs<-variant_ID #the CS
        #check if this CS have any cor >0.8 with CS from others GWAS and mincorr 0.8
        #if yes, extract these variants
        to_add<-res_gwff_cs_others[,to_merge:={
          other_cs<-variant_ID
          mincorr50_vars<-variants_corsf[variant_ID%in%vars_cs&variant_ID2%in%other_cs]
          if(nrow(mincorr50_vars)>0){
            #any 0.8 between those CS
            have_highLD<-mincorr50_vars[,any(abs(r)>0.8)]
            #all variants pairs have r>0.5
            have_all50<-nrow(mincorr50_vars)==length(other_cs)*length(vars_cs)
            rep(have_highLD&have_all50,.N) #rep is to return TRUE or FALSE for each Variants of thise CS
          }else{
            rep(FALSE,.N)
          }
          
        },
        by=.(locuscontext_id)][(to_merge)]$variant_ID
        
        new<-setdiff(to_add,vars_cs)
        if(length(new)>0){
          message(length(new),' variants added')
          
        }
        
        merged_cs<-union(vars_cs,to_add)
        
        data.table(variant_ID=merged_cs)[,original_cs:=variant_ID%in%vars_cs]
        
      },by=.(locuscontext_id,source)]
      
      return(res_gwff_cs_study_ext)
    }),fill = TRUE)
    return(res_gwff_cs_ext[,region:=r])
  }),fill = TRUE)
  
  #iteration 2 and +
  change<-TRUE
  while(change){
    i<-i+1
    message('iteration ', i)
    change=FALSE
    
    gwas_toploci_old<-copy(gwas_new)
    gwas_new<-rbindlist(lapply(regions,function(r){
      message(r)
      res_gwff_cs<-gwas_new[region==r]
      
      #iterate by study
      res_gwff_cs_ext<-rbindlist(lapply(unique(res_gwff_cs$locuscontext_id),function(l){
        message(l)
        res_gwff_cs_study<-res_gwff_cs[locuscontext_id==l]
        res_gwff_cs_others<-res_gwff_cs[locuscontext_id!=l]
        #iterate per CS
        res_gwff_cs_study_ext<-res_gwff_cs_study[,{
          vars_cs<-variant_ID #the CS
          vars_original<-variant_ID[(original_cs)]
          #check if this CS have any cor >0.8 with CS from others GWAS and mincorr 0.8
          #if yes, extract these variants
          to_add<-res_gwff_cs_others[,to_merge:={
            other_cs<-variant_ID
            mincorr50_vars<-variants_corsf[variant_ID%in%vars_cs&variant_ID2%in%other_cs]
            if(nrow(mincorr50_vars)>0){
              #any 0.8 between those CS
              have_highLD<-mincorr50_vars[,any(abs(r)>0.8)]
              #all variants pairs have r>0.5
              have_all50<-nrow(mincorr50_vars)==length(other_cs)*length(vars_cs)
              rep(have_highLD&have_all50,.N) #rep is to return TRUE or FALSE for each Variants of thise CS
            }else{
              rep(FALSE,.N)
            }
            
          },
          by=.(locuscontext_id)][(to_merge)]$variant_ID
          new<-setdiff(to_add,vars_cs)
          if(length(new)>0){
            message(length(new),' variants added')
            change<<-TRUE
            
          }
          merged_cs<-union(vars_cs,to_add)
          
          data.table(variant_ID=merged_cs)[,original_cs:=variant_ID%in%vars_original]
          
        },by=.(locuscontext_id,source)]
        return(res_gwff_cs_study_ext)
      }),fill = TRUE)
      return(res_gwff_cs_ext[,region:=r])
    }),fill = TRUE)
    
    #gwas_new[,cs_num_ext:=ifelse(gwas_toploci_old,cs_num,-cs_num)]
    ntotadded=nrow(gwas_new[!(original_cs)])
    message(ntotadded,' total variants added')
    
  }
  
  
  message('converged after ',i,' iteration')
  gwas_new
  #recreate top_loci table
  fwrite(gwas_new,fp(out,'all_adlocis_extended_any0.8ANDmin0.5union_withfsusie.csv.gz'))
  gwas_new<-fread(fp(out,'all_adlocis_extended_any0.8ANDmin0.5union_withfsusie.csv.gz'))
  
  #saved back the different methods
  #FINEMAP
  table(gwas_new$source)
  gwas_newf<-gwas_new[source=='AD_GWAS_finemapping']
  gwas_newf[,study:=strsplit(locuscontext_id,'_chr')[[1]][1],by='locuscontext_id']
  gwas_newf[,coverage:=ifelse(str_detect(locuscontext_id,'cs95'),'cs_coverage_0.95',
                              ifelse(str_detect(locuscontext_id,'cs70'),'cs_coverage_0.7',
                                     'cs_coverage_0.5'))]
  gwas_newf[,cs_num:=as.numeric(str_extract(locuscontext_id,'[0-9]+$'))]
  
  resfp<-dcast(gwas_newf,
               study+region+variant_ID~coverage,
               value.var = 'cs_num',
               fun.aggregate = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  for (.miss_col in c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5')) {
    if (!.miss_col %in% names(resfp)) resfp[, (.miss_col) := NA_real_]
  }
  #add pip, z
  mtd<-fread(metadata_analysis,header = T,select = 1:6)
  
  res_gw<-rbindlist(lapply(file.path(PROJECT_ROOT,mtd[Method=='AD_GWAS_finemapping']$Path),function(f)fread(f)),fill = T)
## FIX (2026-09-17): the 8-study GWAS fine-mapping files encode the credible-set index as
## text ('susie_rss_0' = no CS, 'susie_rss_1', ...). The script expects integers (0 = no CS);
## with text, `x != 0` and `x > 0` are TRUE for every row (string comparison), so every
## variant became a 'hit'/'cs95' and loci were inflated. Extract the trailing integer.
for (.cc in intersect(c('cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5'), names(res_gw))) {
  if (is.character(res_gw[[.cc]])) res_gw[, (.cc) := as.integer(sub('.*_', '', get(.cc)))]
}
  res_gw[,study:=event_ID]
  resfp<-merge(resfp,res_gw[,.(variant_ID,PIP,z,study,cs_coverage_0.95_original=cs_coverage_0.95,
                               cs_coverage_0.7_original=cs_coverage_0.7, cs_coverage_0.5_original=cs_coverage_0.5)],
               by=c('study','variant_ID'),all.x = T)

    ## ---- FIX (2026-09-17): do not discard original credible sets ----------
    ## resfp is built only from gwas_new (the extension output); the merge above is
    ## all.x=TRUE, so any genuine ORIGINAL credible set that never entered gwas_new
    ## was silently dropped (2355 of 3355 CS variants; 368 regions -> 88, incl. APOE).
    ## Re-add those variants carrying their original cs assignment (no extension).
    .gw_cs <- res_gw[`cs_coverage_0.95`!=0 | `cs_coverage_0.7`!=0 | `cs_coverage_0.5`!=0]
    .miss  <- .gw_cs[!resfp, on=c('study','variant_ID')]
    message('[fix] CS variants: ', nrow(.gw_cs), ' | dropped by extension, re-added: ', nrow(.miss))
    if (nrow(.miss)) {
      resfp <- rbind(resfp, .miss[, .(study, region=gene_ID, variant_ID,
                       `cs_coverage_0.95`, `cs_coverage_0.7`, `cs_coverage_0.5`, PIP, z,
                       `cs_coverage_0.95_original`=`cs_coverage_0.95`,
                       `cs_coverage_0.7_original` =`cs_coverage_0.7`,
                       `cs_coverage_0.5_original` =`cs_coverage_0.5`)], fill=TRUE)
    }
    ## ----------------------------------------------------------------------
  resfp[,event_ID:=study]
  resfp[is.na(cs_coverage_0.95),cs_coverage_0.95:=0]
  resfp[is.na(cs_coverage_0.7),cs_coverage_0.7:=0]
  resfp[is.na(cs_coverage_0.5),cs_coverage_0.5:=0]
  resfp[is.na(cs_coverage_0.95_original),cs_coverage_0.95_original:=0]
  resfp[is.na(cs_coverage_0.7_original),cs_coverage_0.7_original:=0]
  resfp[is.na(cs_coverage_0.5_original),cs_coverage_0.5_original:=0]
  
  #distrib n.extension 
  resfp[,n.extension:=sum(cs_coverage_0.95!=cs_coverage_0.95_original),by=.(study,region,cs_coverage_0.95)]
  summary(unique(resfp[cs_coverage_0.95>0],by=c('study','region','cs_coverage_0.95'))$n.extension)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.000   0.000   2.000   7.559   6.000 144.000 
  
  resfp[,n.extension:=sum(cs_coverage_0.7!=cs_coverage_0.7_original),by=.(study,region,cs_coverage_0.7)]
  summary(unique(resfp[cs_coverage_0.7>0],by=c('study','region','cs_coverage_0.95'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.00    0.00    5.00   18.99   23.00  144.00 
  
  resfp[!is.na(PIP)]
  
  #export
  #fwrite(resfp,ROOTP('interactive_analysis/<user>/export/AD_GWAS_finemapping_109_blocks_top_loci_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
  
  fwrite(resfp,fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
  
  #COLOCS
  table(gwas_new$source)
  
  #ADmeta coloc
  res_cad<-fread(fp(out,'res_coloc_meta_AD.csv.gz'))
  res_cad[,chr:=seqid(variant_ID)]
  res_cad[,pos:=pos(variant_ID)]
  
  
  res_cadnew<-gwas_new[source=='AD_meta_colocalization']
  
  res_cadnew<-merge(res_cadnew,res_cad,all.x = T,
                    by=c('locuscontext_id','variant_ID'))
  res_cadnew[,chr:=seqid(variant_ID)]
  res_cadnew[,pos:=pos(variant_ID)]
  
  res_cadnew<-rbindlist(lapply(unique(res_cad$event_ID),function(g){
    loci<-unique(res_cosad[event_ID==g]$locuscontext_id)
    
    rbind(res_cadnew[locuscontext_id%in%loci&!is.na(event_ID)],res_cadnew[locuscontext_id%in%loci&is.na(event_ID),event_ID:=g][])
    
  }))
  res_cadnew[,context:=event_ID]
  res_cadnew[,gwas_source:=event_ID]
  res_cadnew<-unique(res_cadnew)
  res_cadnew[is.na(locuscontext_id)]
  table(res_cadnew$locuscontext_id)
  unique(res_cadnew$locuscontext_id)
  
  #extension
  res_cadnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
  summary(unique(res_cadnew,by=c('locuscontext_id'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   0.0     0.0     2.0    21.4    11.0  1036.0 
  
  
  fwrite(res_cadnew,fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
  
  #ADxQTL coloc
  res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))
  
  table(gwas_new$source)
  table(gwas_new$source)
  
  res_cnew<-gwas_new[source=='AD_xQTL_colocalization']
  unique(res_cnew$locuscontext_id)
  res_cnew[,context:=paste0('AD_',strsplit(locuscontext_id,'AD_|_ENSG')[[1]][2]),by='locuscontext_id']
  unique(res_cnew$context)
  
  res_cnew<-merge(res_cnew,res_c,all = T,
                  by=c('locuscontext_id','context','variant_ID'))
  
  # res_cnew[,event_ID:=event_ID[!is.na(event_ID)][1],by=.(locuscontext_id,context)]
  res_cnew[,context_coloc:=context_coloc[!is.na(context_coloc)][1],by=.(locuscontext_id,context)]
  #res_cnew[,context:=context[!is.na(context)][1],by=.(locuscontext_id,context)]
  
  res_cnew[,gwas_source:=gwas_source[!is.na(gwas_source)][1],by=.(locuscontext_id)]
  
  
  #extension
  res_cnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
  summary(unique(res_cnew,by=c('locuscontext_id'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.0     3.0    10.0   135.2    49.0  2639.0 
  
  res_cnew[is.na(locuscontext_id)]
  fwrite(res_cnew,fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
  
  #fsusie Coloc
  res_cfsf<-fread(fp(out,'res_coloc_AD_epiQTL.csv.gz'))
  table(res_cfsf$context)
  
  table(gwas_new$source)
  
  res_cfnew<-gwas_new[source=='AD_epiQTL_colocalization']
  res_cfnew<-merge(res_cfnew,unique(res_cfsf[,.(locuscontext_id,gwas_source)]),by='locuscontext_id')
  
  res_cfnew<-merge(res_cfnew,res_cfsf,all = T,
                   by=c('locuscontext_id','gwas_source','variant_ID'))
  
  # res_cfnew[is.na(context)]
  # res_cfnew[,context:=context[!is.na(context)][1],by=.(locuscontext_id)]
  
  #extension
  res_cfnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
  summary(unique(res_cfnew,by=c('locuscontext_id'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.00    0.00    3.00   28.31   36.00  284.00 
  
  res_cfnew[is.na(locuscontext_id)]
  fwrite(res_cfnew,fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
}

#Finemapping####
res_gwf<-fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#create locus_context_id for cs95>cs70>cs50

res_gwf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_gwf[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]

#unified locus
res_gwf[,chr:=seqid(variant_ID)]
res_gwf<-UnifyLoci(res_gwf,variant_col = 'variant_ID',
                   locus_col = 'locuscontext_id',group.by = 'chr')
unique(res_gwf$uni.locus_id)#186
unique(res_gwf[str_detect(credibleset,'cs95')]$uni.locus_id)|>length()#92

res_gwf[,context:=event_ID]
res_gwf[,gwas_source:=event_ID]

fwrite(res_gwf,fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz'))
res_gwf<-fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz'))

mtd[Method=='AD_GWAS_finemapping',summary_file:=fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz')]

# #QC: how many loci without extension?
# res_gwfo<-res_gwf[!is.na(cs_coverage_0.95_original)][cs_coverage_0.5_original!=0|cs_coverage_0.7_original!=0|cs_coverage_0.95_original!=0][,-c('uni.locus_id')]
# res_gwfo[,
#         credibleset:=ifelse(cs_coverage_0.95_original>0,paste0('cs95_',cs_coverage_0.95_original),
#                             ifelse(cs_coverage_0.7_original>0,paste0('cs70_',cs_coverage_0.7_original),
#                                    paste0('cs50_',cs_coverage_0.5_original)))]
# res_gwfo[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]
# 
# res_gwfo<-UnifyLoci(res_gwfo,variant_col = 'variant_ID',
#                    locus_col = 'locuscontext_id',group.by = 'chr')
# unique(res_gwfo$uni.locus_id)|>length()#174/195
# unique(res_gwfo[str_detect(credibleset,'cs95')]$uni.locus_id)|>length()#92/91
# 


#ADxQTL coloc####
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#add the coloc locus overlap between studies
res_c[,chr:=seqid(variant_ID)]
res_c[,pos:=pos(variant_ID)]
res_c<-UnifyLoci(res_c[,-'uni.locus_id'],variant_col = 'variant_ID',
                 locus_col = 'locuscontext_id',group.by = c('chr'),rm_overlap = T)
unique(res_c$uni.locus_id)#146
table(res_c$context)
setdiff(res_c$context,contexts$context)
fwrite(res_c,fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

mtd[Method=='ColocBoost'&Modality=='AD_xQTL_colocalization',summary_file:=fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]
mtd[Method=='ColocBoost'&Modality=='AD_xQTL_colocalization']


#ADAD colocs####
res_cad2<-fread(fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

mtd[Method=='ColocBoost'&Modality=='AD_meta_colocalization',summary_file:=fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]


#fsusie Coloc####
res_cf<-fread(fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#add the coloc locus overlap between studies
res_cf[,chr:=seqid(variant_ID)]
res_cf[,pos:=pos(variant_ID)]
res_cf<-UnifyLoci(res_cf[,-'uni.locus_id'],variant_col = 'variant_ID',
                  locus_col = 'locuscontext_id',group.by = c('chr'),rm_overlap = T)
unique(res_cf$uni.locus_id)#24
table(res_cf$context)#
res_cf[context=='']#the extension
fwrite(res_cf,fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

mtd[Method=='Coloc',summary_file:=fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]
fwrite(mtd,metadata_analysis)

#  4) the TWAS/MR####
table(mtd$Method)

res_twas<-fread(file.path(PROJECT_ROOT,mtd[Method=='twas']$Path))
setdiff(res_twas$context,contexts$context)
res_twas[,gene_ID:=str_extract(molecular_id,'ENSG[0-9]+')]
res_twas[,event_ID:=context]
res_twas[str_detect(context,'_chr|_un'),context:=strsplit(context[1],'_chr|_gp_')[[1]][1],by='context']
setdiff(res_twas$context,contexts$context)#ok

res_twas[,gwas_source:=paste0('AD_',gwas_study)]
res_twas<-res_twas[,-'gwas_study']

#TWAS signif
res_twas<-unique(res_twas[order(gene_ID,context,gwas_source,twas_pval)],by=c('gene_ID','context','gwas_source','method'))

#TWAS sig: (pvalue<2.5e-6 for >50% methods OR (pvalue<2.5e-6 AND the best method))
nmethods<-unique(res_twas$method)|>length()
res_twas[!is.na(twas_pval),TWAS_signif:=(twas_pval<2.5e-6&(is_selected_method))|(sum(twas_pval<2.5e-6)>=nmethods/2),
         by=c('gene_ID','block','context','gwas_source')]
res_twas[(TWAS_signif)]$gene_ID|>unique()|>length()#598 genes

fwrite(res_twas,fp(out,'res_AD_XWAS.csv.gz'))
res_twas<-fread(fp(out,'res_AD_XWAS.csv.gz'))


#    - MR####
# if(!any(mtd$Method=='MR')){
#   if(mtd[Method=='MR',!file.exists(file.path(PROJECT_ROOT,Path))]){
#     #export table: RUN export_mr.R
#     #for all gene_id
#     #get file paths
#     #format: /data/analysis_result/twas/[context].[block].mr_results.tsv.gz
#     system('ls /data/analysis_result/twas/AD >files_twas_folder.txt')
#     files<-fread('files_twas_folder.txt',col.names = 'file',header = F)[str_detect(file,'mr_result.tsv.gz$')]
#     files[,context:=strsplit(file,'\\.')[[1]][1],by='file']
#     files[,block:=strsplit(file,'\\.')[[1]][2],by='file']
#     files[,result_type:=strsplit(file,'\\.')[[1]][3],by='file']
#     
#     
#     table(files$context)
#     filesfa<-merge(files,unique(res_twas[,.(block,gene_ID)]),by = 'block',allow.cartesian = T)
#     table(filesfa$context)
#     
#     mr_res<-rbindlist(lapply(unique(filesfa$file),function(f){
#       i<-which(unique(filesfa$file)==f)
#       if(i%%100==0){
#         message(i,'/',length(unique(filesfa$file)))
#       }
#       genes=filesfa[file==f]$gene_ID
#       p<-fp(ROOTP('analysis_result/twas/'),f)
#       mr<-fread(p)
#       if(nrow(mr)>0){
#         return(mr[!is.na(meta_pval)][gene_name%in%genes])
#         
#       }else{
#         return(mr)
#         
#       }
#       
#       
#     }),fill=T)
#     
#     setnames(mr_res,'gene_id','gene_ID')
#     
#     mr_resf<-mr_res[!is.na(gene_ID)]
#     
#     #add chr pos of twas
#     mr_resf<-merge(mr_resf,res_twas[,.(`#chr`,start,end,gene_ID,TADB_start,TADB_end,context)],
#                    by=c('context','gene_ID'),all.x = T)
#     
#     #export
#     fwrite(mr_resf[order(`#chr`,start)],ROOTP('interactive_analysis$AD_LOCI_ROOT/interactive_analysis/<user>/export/FunGen_mr.exported.bed.gz'))
#     mr_resf<-fread(ROOTP('interactive_analysis$AD_LOCI_ROOT/interactive_analysis/<user>/export/FunGen_mr.exported.bed.gz'))
#     
#     #add to mtd
#     mtd<-fread(metadata_analysis)
#     mtd<-rbind(mtd[Method!='MR'],data.table(
#       'Data Type'='Gene & GWAS',
#       Cohort='ROSMAP & MSBB & AD',
#       Method='MR',
#       Path='interactive_analysis$AD_LOCI_ROOT/interactive_analysis/<user>/export/FunGen_mr.exported.bed.gz'),
#       fill=TRUE)
#     mtd[Method=='MR']
#     fwrite(mtd,metadata_analysis)
#   }
#  
# }

res_mr<-mtd[Method=='MR',fread(file.path(PROJECT_ROOT,Path))]
res_mr[str_detect(context,'DLPFC_Bennett_pQTL'),context:='DLPFC_Bennett_pQTL']
res_mr[!is.na(Q_pval)][duplicated(paste(gene_ID,context,gwas_study))]

#harmonize gwas source /study
res_mr[,gwas_source:=paste0('AD_',gwas_study)]
table(res_mr$gwas_source)

# add to TWAS and call MR signif 
res_twas<-fread(fp(out,'res_AD_XWAS.csv.gz'))

res_mrtwas<-merge(unique(res_mr[order(gene_ID,context,gwas_source,-num_CS,I2)],by=c('gene_ID','context','gwas_source')),
                  unique(res_twas[!is.na(twas_pval)][order(gene_ID,context,gwas_source,twas_pval)],by=c('gene_ID','context','gwas_source')),all=TRUE)
res_mrtwas[,MR_signif:=TWAS_signif&cpip>0.5&num_CS>=2&I2<0.5]

res_mrtwas[(TWAS_signif)]
res_mrtwas[(MR_signif)]
fwrite(res_mrtwas,fp(out,'res_AD_XWAS_MR.csv.gz'))
fwrite(res_mrtwas[(TWAS_signif)],fp(out,'res_AD_XWAS_MR_filtered_TWAS_sig.csv.gz'))

mtd[Method%in%c('twas','MR'),summary_file:=fp(out,'res_AD_XWAS_MR_filtered_TWAS_sig.csv.gz')]
fwrite(mtd,metadata_analysis)

#cTWAS####
#$cs column not left blank (having values such as L1 L2 L3... 
#although the $pip can be small) and also check $pip column values being greater than 0.75
contexts<-fread(contexts_metadata)  ## was fp(out,..): contexts_metadata.csv lives in staging (cwd), never written to out/
ROOTP('analysis_result/ctwas/export/')
res_ctwas<-fread(file.path(PROJECT_ROOT,mtd[Method=='cTWAS']$Path))

res_ctwas[cs!=''&susie_pip>0.75]$gene_name|>unique() #68
res_ctwas[cs!=''&susie_pip>0.75]$context|>unique() #11

setdiff(res_ctwas$context,contexts$context)#OK
setnames(res_ctwas,'gene_id','gene_ID')

res_ctwasf<-res_ctwas[cs!=''&susie_pip>0.75]
res_ctwasf<-res_ctwasf[,-c('group','type','gene_name')]
res_ctwasf[,cTWAS_signif:=TRUE]
res_ctwasf[(cTWAS_signif)]

fwrite(res_ctwasf,fp(out,'res_AD_cTWAS_pip075.csv.gz'))
res_ctwasf<-fread(fp(out,'res_AD_cTWAS_pip075.csv.gz'))

mtd[Method%in%c('cTWAS'),summary_file:=fp(out,'res_AD_cTWAS_pip075.csv.gz')]
fwrite(mtd,metadata_analysis)


#Metadata of harmonized analysis check and saving####
unique(mtd[!is.na(summary_file)][,.(Method,summary_file)])


#confirm all have variant_ID, gene_ID, context, locuscontext_id
cols_to_check<-c('variant_ID', 'gene_ID', 'context','locuscontext_id')

unique(mtd[!is.na(summary_file)&summary_file!=''][,.(Method,summary_file)])[,{
  message(summary_file[1])
  missing<-setdiff(cols_to_check,colnames(fread(summary_file[1],nrows = 1)))
  if(length(missing)>0){
    message('missing ',paste(missing,collapse = ' '))
    
  }else{
    message('OK')
  }
},by='summary_file']

#ok except TWAS/MR/cTWAS for variant_ID and locuscontext_id, interaction/QR lacking locuscontext_id, and epiQTL lacking gene_id but normal
mtd[,variant_level_method:=!Method%in%c('MR','twas','cTWAS')]
fwrite(mtd,metadata_analysis)


#II) get unified 'AD associated' loci ####
#outputs:ADlocus variant ADlocus_event ADmethod
#here we want to integrate GWAS locus found in i) single gwas finemapping, ii) ADxQTL coloc, iii) ADxAD colocs
mtd<-fread(metadata_analysis)
mtd[summary_file=='']

res_gwf<-fread(mtd[Method=='AD_GWAS_finemapping']$summary_file[1])

res_cad<-fread(mtd[Modality=='AD_xQTL_colocalization']$summary_file[1])[str_detect(context,'^AD')]

res_cad2<-fread(mtd[Modality=='AD_meta_colocalization']$summary_file[1])
res_cad3<-fread(mtd[Method=='Coloc']$summary_file[1])


res_ad<-rbindlist(list(res_gwf[,.(variant_ID,locuscontext_id,AD_method='single_gwas_finemapping',gwas_source=event_ID,credibleset,PIP)],
                       res_cad[,.(variant_ID,locuscontext_id,AD_method='ADxQTL_coloc',gwas_source=gwas_source,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad2[,.(variant_ID,locuscontext_id,AD_method='ADxAD_coloc',gwas_source=event_ID,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad3[,.(variant_ID,locuscontext_id,AD_method='ADFsusie_coloc',gwas_source=gwas_source,SNP_PPH4,L_PP.H4.abf,xQTL_L,AD_L)]),fill = T)
res_ad[,chr:=seqid(variant_ID)]


res_ad<-UnifyLoci(res_ad,variant_col = 'variant_ID',locus_col = 'locuscontext_id',
                  group.by = 'chr')
setnames(res_ad,'uni.locus_id','ADlocus')
message('UNIFYLOCI_NLOCI: ',uniqueN(res_ad$ADlocus))
fwrite(res_ad,fp(out,'checkpoint_res_ad_unified.csv.gz'))


#add gwas zscore
res_ad<-rbindlist(lapply(unique(res_ad$gwas_source),function(g){
  message(g)
  file=ps(ROOTP('analysis_result/summary_stats_qced/AD_GWAS_RSS_QC_RAISS_imputed_concatenate_sumstat/'),g,'_RSS_QC_RAISS_imputed.tsv.gz')
  if(file.exists(file)){
    gwas<-fread(file)
    merge(res_ad[gwas_source==g],gwas[,.(variant_ID=paste0('chr',variant_id),
                                         gwas_zscore=z)],
          all.x = T,by='variant_ID')
  }else{
    message(g,' zscores not found')
    res_ad[gwas_source==g]
  }
  
}),fill = T)

#ncs95 or ADxAD coloc
res_ad[str_detect(credibleset,'cs95')]$ADlocus|>unique()|>length()#88
res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc']$ADlocus|>unique()|>length()#116
res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc'|AD_method=='ADxAD_coloc'|AD_method=='ADxQTL_coloc']$ADlocus|>unique()|>length()#191
#res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc'|AD_method=='ADxAD_coloc'|AD_method=='ADxQTL_coloc'|AD_method=='TWAS']$ADlocus|>unique()|>length()#191
res_ad$ADlocus|>unique()|>length()#249

#add more meaningful ADLocus name: chr-[start]-[end]
res_ad[,is.cs95:=any(str_detect(credibleset,'cs95'),na.rm = T),by='ADlocus']
res_ad[,ADlocusID:=paste(seqid(ADlocus),min(pos(variant_ID)),max(pos(variant_ID)),sep = '_'),by='ADlocus']

res_ad[,ADlocus.n.variant.union:=length(unique(variant_ID)),by='ADlocus']
setcolorder(res_ad,c('ADlocusID','ADlocus.n.variant.union','is.cs95'))

table(unique(res_ad,by=c('ADlocus','variant_ID'))$ADlocusID)

#save
unique(res_ad[(is.cs95)]$ADlocusID)
fwrite(res_ad,fp(out,'AD_loci_unified_cs95orColocs.csv.gz'))
res_ad<-fread(fp(out,'AD_loci_unified_cs95orColocs.csv.gz'))


#one variant per line
#with  chr, pos, effect_allele, GWAS methods,
res_ad[,variant_inclusion_probability:=ifelse(is.na(vcp),PIP,vcp)]
res_ad[is.na(variant_inclusion_probability),variant_inclusion_probability:=SNP_PPH4]

res_ad[,max_variant_inclusion_probability:=max(variant_inclusion_probability,na.rm = T),by='variant_ID']

res_ad[,max_variant_inclusion_probability_method:=AD_method[which.max(variant_inclusion_probability)],by='variant_ID']
res_ad[,n.variants:=length(unique(variant_ID)),by='ADlocus']

res_ad[,chr:=seqid(variant_ID,only_num =T )]
res_ad[,pos:=pos(variant_ID)]
res_ad[,max_zscore:=gwas_zscore[which.max(abs(gwas_zscore))],by='variant_ID']
res_ad[,min_pval:=getPval(max_zscore)]
res_ad[,susie_coverage:=sort(unique(str_extract(credibleset,'cs[0-9]+')),decreasing = T)[1],by='ADlocusID']

res_ad[order(variant_ID,-variant_inclusion_probability),GWAS_methods:=paste(unique(AD_method),collapse =  '|'),by='variant_ID']

res_ad[order(variant_ID,-variant_inclusion_probability),gwas_sources:=paste(unique(gwas_source),collapse =  '|'),by='variant_ID']
#res_ad[,gwas_source_effect:=NULL]
res_ad[!is.na(gwas_zscore),gwas_source_effect:=paste0(str_extract(gwas_source,'[BWJK]'),ifelse(gwas_zscore>0,'+','-'))]
res_ad[is.na(gwas_zscore),gwas_source_effect:=str_extract(gwas_source,'[BWJK]')]
res_ad[order(variant_ID,-variant_inclusion_probability),gwas_sources_effects:=paste(gwas_source_effect[!is.na(gwas_source_effect)&!duplicated(str_extract(gwas_source_effect,'[BWJK]'))],collapse =  '|'),by='variant_ID']


# remove variants with min_pval>1e-5 or in cs70 or cs50 only
res_adf<-res_ad[ADlocus%in%ADlocus[min_pval<1e-5]]
unique(res_adf$ADlocusID) #183
#only cs95 or coloc
res_adf<-res_adf[ADlocusID%in%ADlocusID[susie_coverage=='cs95'|str_detect(GWAS_methods,'coloc')]]
length(unique(res_adf$ADlocusID)) #183
unique(res_adf[min_pval<5e-8]$ADlocusID)|>length() #90
unique(res_adf[susie_coverage=='cs95']$ADlocusID)|>length()  #88
res_adf[is.na(is.cs95),is.cs95:=FALSE]


res_adfv<-unique(res_adf[order(variant_ID,-abs(gwas_zscore),-variant_inclusion_probability)][,.(ADlocusID,
                                                                                                n.variants,
                                                                                                chr,
                                                                                                pos,
                                                                                                variant_ID,
                                                                                                GWAS_methods,
                                                                                                max_variant_inclusion_probability,
                                                                                                max_variant_inclusion_probability_method,
                                                                                                gwas_sources,
                                                                                                gwas_sources_effects,
                                                                                                is.cs95,
                                                                                                max_zscore,
                                                                                                min_pval,
                                                                                                susie_coverage,
                                                                                                ADlocus)],
                 by='variant_ID')

res_adfv[,max_variant_inclusion_probability_rank:=rank(-max_variant_inclusion_probability),by='ADlocus']
res_adfv[,n.variants_0.1:=sum(max_variant_inclusion_probability>0.1),by='ADlocus']

res_adfv[order(chr,ADlocus,-max_variant_inclusion_probability)]
res_adfv[is.na(susie_coverage)]
table(res_adfv$ADlocusID)
#unique(res_adfv[,.(ADlocus,GWAS_methods)][order(-str_length(GWAS_methods))],by='ADlocus')[101:194]
res_adfv[,median_pos:=median(pos),by='ADlocus']
#add locus index
adloci<-unique(res_adfv[order(chr,median_pos,-max_variant_inclusion_probability)],by='ADlocus')
adloci[,locus_index:=1:.N]
res_adfv<-merge(res_adfv[,-'locus_index'],adloci[,.(ADlocus,locus_index)])
setcolorder(res_adfv,'locus_index')
unique(res_adfv[str_detect(GWAS_methods,'Fsusie')],by='ADlocusID')|>nrow()#18
unique(res_adfv[GWAS_methods=='ADFsusie_coloc'],by='ADlocusID')|>nrow()#2

#add v2f score
files<-list.files(cv2f_score_dir,pattern = '.cv2f$',full.names = TRUE)
v2f<-rbindlist(lapply(files,function(f)fread(f)))
v2f[!str_detect(SNP,'rs')]
res_adfv<-merge(res_adfv,v2f[,.(chr=CHR,pos=BP,cV2F,SNP)],
                by=c('chr','pos'),all.x=TRUE)


res_adfv[is.na(cV2F)]#867/7600
res_adfv[!is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]
res_adfv[is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]
res_adfv[is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]$max_variant_inclusion_probability|>summary()


res_adfv[!is.na(cV2F),cV2F_rank:=rank(-cV2F),by='ADlocus']

fwrite(res_adfv[order(locus_index)],fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))
res_adfv<-fread(fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))

#export
fwrite(res_adfv[order(locus_index)],'../../../../../xqtl-resources/data/genes/AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz')




#III) Get the merged long variant-level table  AD xQTL overlap: each row a variant-ADlocus-Method-context-gene####
#add those ADlocus annot to all variant level summ table
mtd<-fread(metadata_analysis)  ## was fp(out,..): config file lives in staging (cwd), never written to out/
res_adfv<-fread(fp(out,'/AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))
update_summary_ad=FALSE

mtd[file.exists(summary_file),summary_file_ad:=paste0(tools::file_path_sans_ext(summary_file[1],compression = T),'_overlapADloci.csv.gz'),by='summary_file']

mtd[!file.exists(summary_file)&summary_file!='']

#first check if need flipping to harmonize variant_ID
mismatches<-unique(mtd,by='summary_file')[file.exists(summary_file)&(variant_level_method),{
  message(Method[1])
  
  res<-fread(summary_file[1],tmpdir = ROOTP('repro/tmp/'))
  res[,chr:=seqid(variant_ID,only_num = T)]
  res[,pos:=pos(variant_ID)]
  res[,variant_id:=variant_ID]
  
  res<-merge(res[,-c('ADlocus','locus_index','ADlocusID','variant_ID')],
             res_adfv[,.(variant_ID,chr,pos,ADlocus,locus_index,ADlocusID)],all.x=T,by=c('chr','pos'))
  print(table(res[!is.na(variant_ID)][variant_id!=variant_ID]$context))
  res[!is.na(variant_ID)][variant_id!=variant_ID][,.(variant_id,variant_ID,context)]
},by=c('Method','Modality')]
mismatches<-unique(mismatches)
unique(mismatches$variant_ID)#213/5318 variants mismatch
unique(mismatches[,.(variant_id,variant_ID)])

effect_cols<-c('conditional_effect','coef','effect_direction','top_effect','beta','z')

mismatches[,to_flip:=ref(variant_id)==alt(variant_ID)&alt(variant_id)==ref(variant_ID)]
mismatches[(to_flip)]
mismatches[ref(variant_id)==alt(variant_ID)&alt(variant_id)!=ref(variant_ID)]
table(mismatches[(to_flip)]$context)
table(mismatches[(to_flip)]$Method)
table(mismatches[(to_flip)][,paste(Method,context)])
fwrite(mismatches,'variants_to_flip.csv.gz')
mismatches<-fread('variants_to_flip.csv.gz')
toflip<-mismatches[(to_flip)]$variant_id

#Annotate each variant level Methods table  with the AD Loci

mtd[file.exists(summary_file)&(variant_level_method),{
  message(Method[1])
  
  if(!file.exists(summary_file_ad[1])|update_summary_ad){
    res<-fread(summary_file[1],tmpdir = ROOTP('repro/tmp/'))
    res[,chr:=seqid(variant_ID,only_num = T)]
    res[,pos:=pos(variant_ID)]
    res[,variant_id:=variant_ID]
    
    if(any(unique(res$variant_id)%in%toflip)){
      #flip variant_ID
      res<-rbind(res[!variant_id%in%toflip],
                 merge(res[variant_id%in%toflip][,-'variant_ID'],unique(mismatches[,.(variant_id,variant_ID)]),by='variant_id'))
      #flip effect
      effect_colsf<-intersect(effect_cols,colnames(res))
      if(length(effect_colsf)>0){

          res<-rbind(res[!variant_id%in%toflip],
                     res[variant_id%in%toflip,(effect_colsf):=lapply(.SD,function(x)-x),.SDcols = effect_colsf]
                     )
        
        
      }
    }
    
    if('locuscontext_id'%in%colnames(res)){
      #n.variants
      res[,n.variant.locus:=.N,by='locuscontext_id']
    }
    res<-merge(res[,-c('ADlocus','locus_index','ADlocusID','chr','pos')],
               res_adfv[,.(variant_ID,chr,pos,ADlocus,locus_index,ADlocusID)],
               all.x=T,by=c('variant_ID'))
    
    res[,is.in.ad.locus:=!is.na(ADlocus)]
    
    res[is.na(variant_ID),variant_ID:=variant_id]
    
    if('locuscontext_id'%in%colnames(res)){
      
      res[,n.variant.overlap.ad:=sum(is.in.ad.locus),by=.(locuscontext_id)]
      res[,overlap.ad:=any(is.in.ad.locus),by=.(locuscontext_id)]
      resf<-res[(overlap.ad)]
      
    }else{
      resf<-res[(is.in.ad.locus)]
      
    }
    
    fwrite(res,summary_file[1])
    
    
    fwrite(resf,summary_file_ad[1])
  }
  
  summary_file_ad[1]
},by=c('summary_file')]


fwrite(mtd,metadata_analysis)
#get the long AD overlap table####

mtd<-fread(metadata_analysis)
res_adfv<-fread(fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))


res_adx<-rbindlist(lapply(unique(mtd[file.exists(summary_file_ad)&(variant_level_method)]$summary_file_ad),function(f){
  message(f)
  m=mtd[summary_file_ad==f]$Method[1]
  message(m)
  if(m=='ColocBoost'){
    m=mtd[summary_file_ad==f]$Modality[1]
    summary_file_ad=mtd[Modality==m]$summary_file_ad[1]
    
  }else if(m=='LR'){
    m='sn_sQTL'
    summary_file_ad=mtd[Method=='LR'&`Data Type`=='sQTL']$summary_file_ad[1]
    
  }else {
    summary_file_ad=mtd[Method==m]$summary_file_ad[1]
    
  }
  res<-fread(summary_file_ad,tmpdir = ROOTP('repro/tmp/'))[,Method:=m]
  message('found associated with ',length(unique(res$ADlocus)),' AD locus')
  return(res)
}),fill = T)
res_adx
res_adx[!str_detect(gene_ID,'^ENSG')]$context|>table() #OK

res_adx[context=='']$Method|>table() #Coloc and Meta but normal
res_adx[context=='']$context|>table() #Coloc and Meta but normal





#populate all locuscontext_id lvl with the AD locus overlapping it
res_adx[is.na(locus_index)]$Method|>table()

res_adx[order(-n.variant.overlap.ad),locus_index:=ifelse(is.na(locus_index),locus_index[!is.na(locus_index)][1],locus_index),
        by='locuscontext_id']
res_adx[,locus_index:=ifelse(is.na(locus_index),locus_index[!is.na(locus_index)][1],locus_index),
        by='variant_ID']

res_adx<-merge(res_adx[,-c('ADlocus','ADlocusID')],unique(res_adfv[,.(locus_index,ADlocus,ADlocusID)]))


#rm variants from extension in the AD methods
res_adx<-res_adx[!(is.na(vcp)&str_detect(Method,'coloc'))]
res_adx<-res_adx[!(is.na(PIP)&str_detect(Method,'finemapping'))]

#add GWAS infos : 
dupcols<-setdiff(intersect(colnames(res_adfv),colnames(res_adx)),c('variant_ID','ADlocus','ADlocusID','locus_index'))
res_adx<-merge(res_adx[,.SD,.SDcols = !dupcols],res_adfv,
               by=c('variant_ID','ADlocus','ADlocusID','locus_index'),all = T)
res_adx[locus_index==1]

res_adx[,only_by_proxi:=!any(context%in%c('AD_Bellenguez_EADB_2022','AD_Bellenguez_EADI_2022',
                                          'AD_Wightman_ExcludingUKBand23andME_2021',
                                          'AD_Kunkle_Stage1_2019')),by='locus_index']
res_adx[,APOE_region:=chr==19&pos>43905790 &pos<45905791]
unique(res_adx[(APOE_region)]$locus_index)


res_adx[is.na(locus_index)]$Method|>table()


res_adx[variant_ID=='chr16:19721950:T:TAATC']

#trans[gene_name=="APOE"][,.(start-1e6,end+1e6)]

# #add gwaszscore for gwas method
# res_adx<-merge(res_adx,unique(res_ad[,.(variant_ID,context=gwas_source,gwas_zscore)]),all.x = T,by=c('variant_ID','context'))
res_adx[(is.in.ad.locus)]|>nrow()

#bind to TWAS/MR
#for twas, get ADlocus-gene connection to annotate 
summary_file=mtd[Method=='twas']$summary_file[1]
summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file,compression = T),'_overlapADloci.csv.gz')
update_twasad<-TRUE
if(!file.exists(summary_file_ad)|update_twasad){
  res_adgenes<-unique(res_adx[,.(ADlocus,gene_ID,ADlocusID,locus_index)])
  
  
  res_tw<-fread(summary_file,tmpdir = ROOTP('repro/tmp/'))
  res_twad<-merge(res_tw[,-c('ADlocus','ADlocusID','locus_index')],res_adgenes,by='gene_ID',all.x = T,allow.cartesian = T)
  res_twadf<-res_twad[!is.na(ADlocus)][(TWAS_signif)]
  unique(res_twadf$gene_ID)
  fwrite(res_twadf,summary_file_ad)
  mtd[Method%in%c('twas','MR'),summary_file_ad:=..summary_file_ad]
  fwrite(mtd,metadata_analysis)
}


res_twadf<-fread(summary_file_ad,tmpdir = ROOTP('repro/tmp/'))
res_twadf[(MR_signif)]
res_twadf[is.na(MR_signif),MR_signif:=FALSE]

#bind it
res_adx<-rbind(res_adx,res_twadf[,Method:='TWAS/MR'],fill=T)

#add cTWAS
summary_file=mtd[Method=='cTWAS']$summary_file[1]
summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file,compression = T),'_overlapADloci.csv.gz')
update_ctwasad<-TRUE
if(!file.exists(summary_file_ad)|update_ctwasad){
  res_adgenes<-unique(res_adx[,.(ADlocus,gene_ID,ADlocusID,locus_index)])
  
  
  res_ctw<-fread(summary_file,tmpdir = ROOTP('repro/tmp/'))
  res_ctwad<-merge(res_ctw[,-c('ADlocus','ADlocusID','locus_index')],res_adgenes,by='gene_ID',all.x = T,allow.cartesian = T)
  res_ctwadf<-res_ctwad[!is.na(ADlocus)]
  unique(res_ctwadf$gene_ID)
  fwrite(res_ctwadf,summary_file_ad)
  mtd[Method%in%c('cTWAS'),summary_file_ad:=..summary_file_ad]
  fwrite(mtd,metadata_analysis)
}

# #lacking AD genes are of interest ?
# trans<-fread(ROOTP('resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list'))
# trans[setdiff(res_ctwad$gene_ID,res_ctwadf$gene_ID),on='gene_id']$gene_name
# # [1] "EPDR1"      "MAP3K1"     "GTPBP1"     "RBX1"       "ZFYVE21"    "RIPK2"      "PLA2G12A"   "HS3ST3B1"   "UBA2"      
# # [10] "COL5A1"     "RPS15A"     "PLXNC1"     "ABCA8"      "ZCCHC24"    "SMG8"       "TRANK1"     "PIK3CD"     "SLC16A11"  
# # [19] "PACS2"      "GPR141"     "XPNPEP3"    "IGHG2"      "AP000295.1" "IFNAR2"     "EEF1G"   
# msig<-fread(ROOTP('interactive_analysis/<user>/github.com/AlexandrePelletier/xqtl-paper/resources/all_CPandGOs_gene_and_genesets.csv.gz'))
# resor<-OR3(trans[unique(res_ctwad$gene_ID),on='gene_id']$gene_name,split(msig$gene,msig$pathway),background =unique(msig$gene) )
# resor #yes participate to e.g. GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY (PIK3CD,MAP3K1,IGHG2); GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS (RBX1,RIPK2)
# 
# #does the CS of these genes overlap with our AD loci?
# res_ctwadnf<-res_ctwad[!gene_ID%in%res_ctwadf$gene_ID]
# res_ctwasal<-fread(ROOTP('analysis_result/ctwas/archive/export/summary/ctwas_single_fmprs_merged.tsv.gz'))
# res_ctwasalf<-merge(res_ctwasal,unique(res_ctwadnf[,.(region_id,cs)]))
# res_ctwasalf[type=='SNP',variant_ID:=paste0('chr',id)]
# unique(res_ctwasalf,by=c('region_id','cs'))
# res_ctwasalfadloci<-merge(res_ctwasalf,res_adfv,by='variant_ID')
# res_ctwadnf[cs=='L1'&region_id=='17_57489969_60570445']
# trans['ENSG00000167447',on='gene_id']#SMG8
# res_adx[locus_index==144&cV2F_rank==1][,.(Method,gene_ID,context)]
# 


#bind it
res_ctwadf<-fread(summary_file_ad,tmpdir = ROOTP('repro/tmp/'))

res_adx<-rbind(res_adx,res_ctwadf[,Method:='cTWAS'],fill=T)

#populate TWAS info at variant level
res_adx[,TWAS_signif:=any(TWAS_signif),by=.(ADlocus,gene_ID,context,gwas_source)]
res_adx[,MR_signif:=any(MR_signif),by=.(ADlocus,gene_ID,context,gwas_source)]

res_adx[,cTWAS_signif:=any(Method=='cTWAS'),by=.(ADlocus,gene_ID,context)]


table(unique(res_adx,by=c('locuscontext_id','ADlocus','gene_ID','context'))$Method)
# AD_GWAS_finemapping     AD_meta_colocalization     AD_xQTL_colocalization           APOE interaction 
# 341                        249                       2458                         39 
# Coloc         fSuSiE_finemapping           msex interaction  multi_context_finemapping 
# 93                        137                          2                        747 
# multi_gene_finemapping                         QR single_context_finemapping          trans_finemapping 
# 527                        848                       2333                          2 
# TWAS/MR 
# 661 
#add lacking contexts
res_adx[is.na(context)]$locuscontext_id


#split contexts FOR sQTL: if 'UP' put as u-sQTL (Unproductive sQTL), if 'PR' as p-sQTL (productive)
res_adx[str_detect(event_ID,':UP:')&!str_detect(context,'u_sQTL'),context:=str_replace(context,'sQTL','u_sQTL')]
res_adx[str_detect(event_ID,':PR:')&!str_detect(context,'p_sQTL'),context:=str_replace(context,'sQTL','p_sQTL')]
res_adx[str_detect(event_ID,':UP:')]$context


#add the context broad and short
contexts<-fread(contexts_metadata)
setdiff(res_adx$context,contexts$context) 
setdiff(contexts$context,res_adx$context)

res_adx<-merge(res_adx[,-c('context_broad','context_short')],
               unique(contexts[,.(context,context_broad,context_short)]),all.x = T,by='context')
table(res_adx[,.(Method,context)])

table(res_adx[,.(Method,context_short)])


#add GENE info tss /tes
if(!file.exists(file.path(out,'genes_infos.csv.gz'))){
  trans<-fread(gene_names)
  
  gtf<-fread(gene_gtf,
             select = c(1,3,4,5,7,9),col.names = c('chr','biotype','start','end','sens','info'))
  gtf[,gene_id:=str_extract(info,'ENSG[0-9]+')]
  trans<-merge(trans,gtf[biotype=='gene'][,.(gene_id,tss=ifelse(sens=='+',start,end),tes=ifelse(sens=='+',end,start))])
  setnames(trans,'gene_id','gene_ID')
  trans<-unique(trans,by='gene_ID')
  fwrite(trans,file.path(out,'genes_infos.csv.gz'))
}else{
  trans<-fread(file.path(out,'genes_infos.csv.gz'),tmpdir = ROOTP('repro/tmp/'))
}

res_adx<-merge(res_adx[,-c('tss','tes','gene_name')],trans[,-c('start','end','#chr')],all.x = T,by='gene_ID')

res_adx[,chr:=seqid(variant_ID[1],only_num = T),by='variant_ID']
res_adx[,pos:=pos(variant_ID[1]),by='variant_ID']

res_adx[,effect_allele:=alt(variant_ID[1]),by='variant_ID']
res_adx[,distance_from_tss:=pos-tss]
res_adx[,distance_from_tes:=pos-tes]
res_adx[str_detect(gene_ID,';')]$gene_ID|>unique()
res_adx[str_detect(gene_ID,';')]$Method|>table()
res_adx[str_detect(gene_ID,'ENSG')&is.na(gene_name)]

#associate to a gene each epiQTL: populate gene name for epiQTL if geneQTL overlap the locus
epi_methods<-c('fSuSiE_finemapping','Coloc')
geneqtl_methods<-c('single_context_finemapping',
                   'multi_context_finemapping',
                   'AD_xQTL_colocalization')

res_adx<-rbind(res_adx[!Method%in%epi_methods|is.na(Method)],
               res_adx[Method%in%epi_methods]|>split(by=c('locus_index','context'))|>lapply(function(epi){
                 epivars<-epi$variant_ID
                 genes_overlapping_epi<-res_adx[Method%in%geneqtl_methods][variant_ID%in%epivars]$gene_name|>unique()
                 
                 epi_genes<-rbindlist(lapply(genes_overlapping_epi,function(g)epi[,gene_name:=g]|>unique()))
                 return(epi_genes)
               })|>rbindlist())
res_adx[,gene_ID:=gene_ID[!is.na(gene_ID)][1],by='gene_name']

res_adx[locus_index==1]
table(res_adx$context_broad)
table(res_adx$context_short)


#add broader context and qtl type
res_adx[,context_broad2:=ifelse(context_broad%in%c('bulk_monocyte_eQTL','bulk_macrophage_eQTL','bulk_microglia_eQTL'),'Immune (bulk)',
                                ifelse(str_detect(context_broad,'bulk_brain'),'Brain',
                                       str_remove(context_broad,'_eQTL|_snATAC|_brain|_sQTL')))]
table(res_adx$context_broad2)

res_adx[,qtl_type:=str_extract(context_broad,'(u_|p_)?[a-z]+QTL')]

res_adx[str_detect(context_broad,'snATAC'),qtl_type:='caQTL']
table(res_adx$qtl_type)

#Fix some numerical value
res_adx[twas_z==Inf,twas_z:=max(abs(res_adx[!is.infinite(twas_z)&!is.na(twas_z)][['twas_z']]),na.rm = TRUE)]
res_adx[twas_z==-Inf,twas_z:=-max(abs(res_adx[!is.infinite(twas_z)&!is.na(twas_z)][['twas_z']]),na.rm = TRUE)]


#Summarize table creating column at variant gene -main context level
#here we create the confidence score per variant-gene-context short
unique(res_adx[,.(context,context_short)])

## FIX: SummarizeTable() (gene_prio_utils.R) references many columns by bare name; an
## evidence type with zero rows this run leaves its columns absent -> "Object 'X' not found".
## Pre-create only the columns SummarizeTable READS. Never pre-create ones it assigns to.
.SummarizeTable_orig <- SummarizeTable
## Columns SummarizeTable ASSIGNS to must NOT be pre-created: a bare NA column is
## logical, and data.table coerces later assignments INTO that type ('CL6' -> NA,
## 3 -> TRUE) rather than upgrading it. Only pre-create columns it merely READS.
.ST_assigned <- local({
  txt <- readLines('gene_prio_utils.R', warn = FALSE)
  st  <- grep('^\\s*SummarizeTable\\s*<-\\s*function', txt)
  blk <- txt[st[1]:length(txt)]
  nxt <- grep('^[A-Za-z.][A-Za-z0-9._]*\\s*<-\\s*function', blk); nxt <- nxt[nxt > 1]
  if (length(nxt)) blk <- blk[1:(nxt[1] - 1)]
  unique(sub('[[:space:]]*:=$', '',
    unlist(regmatches(blk, gregexpr('[A-Za-z._][A-Za-z0-9._]*[[:space:]]*:=', blk)))))
})
message('guard: ', length(.ST_assigned), ' cols are assigned by SummarizeTable; leaving those to it')
SummarizeTable <- function(res_adx, group.by, ...) {
  syms <- unique(all.names(body(.SummarizeTable_orig)))
  syms <- syms[grepl('^[A-Za-z_][A-Za-z0-9._]*$', syms)]
  skip <- c(names(formals(.SummarizeTable_orig)), 'res_adxv', 'T', 'F', names(res_adx), .ST_assigned)
  syms <- setdiff(syms, skip)
  syms <- syms[!vapply(syms, function(s) exists(s, envir = globalenv()), logical(1))]
  if (length(syms)) {
    message('guard: read-only cols absent this run, added as NA: ', paste(syms, collapse = ', '))
    for (s in syms) res_adx[, (s) := NA]
  }
  for (.att in 1:25) {
    .r <- tryCatch(.SummarizeTable_orig(res_adx, group.by = group.by, ...), error = function(e) e)
    if (!inherits(.r, 'error')) return(.r)
    .msg <- conditionMessage(.r)
    .mm <- regmatches(.msg, regexpr("bject '[^']+' not found", .msg))
    if (!length(.mm)) stop(.r)
    .miss <- sub("' not found$", '', sub("^bject '", '', .mm))
    message('guard retry: creating missing read-only column ', .miss)
    res_adx[, (.miss) := NA]
  }
  stop('guard: too many missing-column retries')
}

res_adx<-SummarizeTable(res_adx,group.by = 'context_short')
colnames(res_adx)
res_adx[variant_ID=='chr7:143414109:CCA:C'&gene_name=='ZYX'][,.(TWAS_signif_gene,xQTL_effects)]

unique(res_adx[,.(context,context_short)])

fwrite(res_adx,fp(out,'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz'))
table(res_adx$Method)
#no Trans
res_adx<-fread(fp(out,'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz'))

fwrite(res_adx[!str_detect(Method,'trans')],fp(out,'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5_noTrans.csv.gz'))


# #OPTIONAL
# #keep only columns of interest
# # #filter cols: keep only AD locus, gene, method, context, variant_id and gwas infos,  and key method specific infos
# colnames(res_adx)|>cat(sep='\n')
# 
# if(file.exists('long_table_columns_selection.csv')){
#   
#   cols<-fread('long_table_columns_selection.csv')$column_name
#   setdiff(cols,colnames(res_adx))
#   setdiff(colnames(res_adx),cols)
#   
#   colsf<-intersect(cols,colnames(res_adx))
#   res_adxf<-res_adx[,..colsf]
#   res_adxf<-res_adxf[order(locus_index)]
#   res_adxf
#   
#   fwrite(res_adxf,fp(out,'res_allanalysis_ADloci_overlap_selectedcols.csv.gz'))
#   
# }





#IV) WIDE TABLE CREATION  ####
res_adx<-fread(fp(out,'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz'),tmpdir = ROOTP('repro/tmp/'))

res_adxub<-WideTable(res_adx,split.by=c('context_broad2','qtl_type'))
res_adxub[, xQTL_effects := str_replace_all(xQTL_effects, '\\(CL([0-9])', '(T\\1')]
res_adxub[locus_index==152][gene_name=='ARL17A']$xQTL_effects|>unique()
res_adxub[gene_name=='KNOP1'][is.na(locus_index)]|>nrow()
res_adxub[gene_name=='KNOP1'][,.(variant_ID,locus_index,xQTL_effects,gene_name,Method,coverage_xqtl)]
res_adxub[variant_ID=='chr17:45667624:C:T'&gene_name=='ARL17A'][,.(TWAS_signif_gene,xQTL_effects)]

fwrite(res_adxub,fp(out,'res_AD_variants_xQTL.csv.gz'))
res_adxub<-fread(fp(out,'res_AD_variants_xQTL.csv.gz'))

## --- T6 backfill (pipeline-native) ------------------------------------
## T6 means "TWAS/MR evidence only, with no localized AD-xQTL support".
## The tier chain in gene_prio_utils.R is an if/else over rows of the xQTL
## overlap table with T6 as its final else, so a gene with no localized
## xQTL support never reaches that branch and comes back untiered.
## The gene-level XWAS/MR and cTWAS evidence that DEFINES T6 is produced by
## this same release, so read it here. (An earlier attempt keyed off
## TWAS_signif_gene / MR_signif_gene in res_adxub: those columns exist by
## name but are empty at this stage, because the gene-level flags are not
## joined until build_shiny_data.R. That run matched 0 genes.)
.blankt <- function(v) is.na(v) | v == ''
.t6src <- c(fp(out, 'res_AD_XWAS_MR_filtered_TWAS_sig_overlapADloci.csv.gz'),
            fp(out, 'res_AD_cTWAS_pip075_overlapADloci.csv.gz'))
.t6ids <- unique(unlist(lapply(.t6src, function(f) {
  if (!file.exists(f)) { cat('[T6 backfill] MISSING:', f, '\n'); return(character(0)) }
  x  <- fread(f)
  gc <- intersect(c('gene_ID','gene_id'), names(x))
  if (!length(gc)) { cat('[T6 backfill] no gene column in', f, '\n'); return(character(0)) }
  sub('\\..*$', '', as.character(x[[gc[1]]]))
})))
cat('[T6 backfill] TWAS/MR/cTWAS genes in this release:', length(.t6ids), '\n')
.key    <- sub('\\..*$', '', as.character(res_adxub$gene_ID))
.tiered <- unique(.key[!.blankt(res_adxub$top_confidence)])
.t6hit  <- setdiff(intersect(.t6ids, unique(.key)), .tiered)
cat('[T6 backfill] genes given T6:', length(.t6hit), '\n')
if (length(.t6hit)) {
  res_adxub[.key %in% .t6hit & .blankt(top_confidence), top_confidence := 'T6']
  cat('[T6 backfill] rows updated:',
      sum(.key %in% .t6hit & res_adxub$top_confidence == 'T6'), '\n')
}


#FILTER: keep only variants with GWAS PIP/VCP > 0.1, for maximum of 5
# if non of the variant has GWAS PIP/VCP > 0.1 we just show top one based on GWAS PIP/VCP and on xQTLPIP/VCP
res_adxub[,top_variants:=((max_variant_inclusion_probability>=0.1)|cV2F_rank<=5)|((max_variant_inclusion_probability_rank==1|variant_rank_xqtl==1|(rank(pval)==1&!is.na(pval)))),by='locus_index']
res_adxubf<-res_adxub[(top_variants)][!is.na(locus_index)][variant_ID!='']

## ---- FIX 2026-09-17 keep every gene that has a tier ----------------------
## top_variants is ranked per locus_index, so a gene that doesn't own its
## locus's top variant is dropped from the workbook entirely (e.g. EPHX2 at
## locus 66, outranked by a neighbouring gene). Add back one row -- the gene's
## own best variant -- for every gene that carries a confidence level but has
## no surviving row, so the sheet covers the same genes as the tier table.
.blank2 <- function(v) is.na(v) | v==''
.cand <- res_adxub[!.blank2(top_confidence) &
                   !is.na(locus_index) & !(gene_ID %in% unique(res_adxubf$gene_ID))]
if (nrow(.cand)) {
  .cand[, .pick := -fifelse(is.na(max_variant_inclusion_probability), -Inf,
                            as.numeric(max_variant_inclusion_probability))]
  .cand[, .hasvar := as.integer(.blank2(variant_ID))]
  setorder(.cand, .hasvar, top_confidence, .pick, na.last = TRUE)
  .add <- .cand[, .SD[1], by = .(locus_index, gene_ID)]
  .add[, c('.pick','.hasvar') := NULL]
  message('[sheet fix] genes re-added to workbook table: ', uniqueN(.add$gene_ID),
          ' (rows: ', nrow(.add), ')')
  res_adxubf <- rbind(res_adxubf, .add, fill = TRUE)
}
## -------------------------------------------------------------------------

#### gene-level confidence tiers (folded in from repro_tier.R) ##############
## res_adxubf already carries top_confidence, verified byte-identical to
## str_replace(str_extract(xQTL_effects,'(?:CL|T)[0-9](?=,n=)'),'CL','T') on the 20260911 build (0 disagreements
## over 89,087 non-NA rows).  A gene's tier is the best (lowest) CL it reaches.
gene_tiers <- res_adxub[!is.na(top_confidence) & top_confidence != '',
                         .(tier = min(as.integer(str_extract(top_confidence, '[0-9]')))),
                         by = .(gene_ID, gene_name)]
gene_tiers[, tier := paste0('T', tier)]
setorder(gene_tiers, tier, gene_name)
fwrite(gene_tiers, fp(out, 'gene_tier_assignment.csv'))
cat('[tiers] genes with a tier:', nrow(gene_tiers), '  ->', fp(out, 'gene_tier_assignment.csv'), '\n')
print(gene_tiers[, .N, by = tier][order(tier)])
############################################################################

nrow(res_adxubf)#4161
unique(res_adxubf$locus_index)
res_adxubf[variant_ID=='chr16:19711044:T:C'][gene_name=='KNOP1']$locus_index|>table()
res_adxubf[gene_name=='KNOP1'][,.(variant_ID,locus_index,xQTL_effects,gene_name,Method,variant_rank_xqtl)]

#some stats
unique(res_adxubf$gene_name)|>length()#511
unique(res_adxubf[min_pval<5e-8]$gene_name)#259
res_adxubf$top_confidence|>table()
# CL1  CL2  CL3  CL4  CL5 
# 198   72  233  716 1266 


unique(res_adxubf[order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()

unique(res_adxubf[gwas_significance=='ns'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()


unique(res_adxubf[gwas_significance!='ns'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()
unique(res_adxubf[gwas_significance=='genome wide'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()

res_loc<-unique(res_adxubf[!is.na(locus_index)][order(locus_index,top_confidence)],
                by='locus_index')
res_loc[,gwas_significance:=ifelse(gwas_significance=='ns','p<1e-5',ifelse(gwas_significance=='suggestive','p<1e-6','p<5e-8'))]
ggplot(res_loc)+geom_bar(aes(x=top_confidence,fill=gwas_significance))+theme_bw()

res_loc[gwas_significance=='p<5e-8']|>nrow()

#Main Excel Sheet creation ####
#get the columns metadata ready

cols<-fread(columns_metadata)[(keep==1)]  ## was fp(out,..): config file lives in staging (cwd), never written to out/
colsmtd<-fread(excel_metadata)  ## was fp(out,..): config file lives in staging (cwd), never written to out/
colorsmtd<-fread(pattern_coloring)  ## was fp(out,..): config file lives in staging (cwd), never written to out/

cols<-PrepColsMtd(cols,colsmtd,res_adxubf)
unique(cols[,.(parent_column,grandparent_column)])|>tail(100)
unique(cols[,.(parent_column,grandparent_column)])|>tail(100)

setdiff(colnames(res_adxubf),cols$r_name)
#get the main sheet
wb<-CreateExcelFormat(res_adxubf,columns_mtd =cols,
                      colors = colorsmtd)

saveWorkbook(wb, fp(out, paste0('unified_AD_loci_xQTL_summary_', STAMP, '.xlsx')), overwrite = TRUE)

#One sheet per broad context Creation #####
#split per context keeping central information 

res_adxubf_list<-split(res_adxubf,by = 'context_broad2')

for(cont in colsmtd[wildcard=='context_broad2']$r_name){
  message(cont)
  res_adxc<-res_adx[context_broad2==cont]
  
  
  res_adxc<-SummarizeTable(res_adxc,group.by = 'qtl_type')
  
  res_adxcub<-WideTable(res_adxc)
  
  
  #FILTER: keep only variants with GWAS PIP/VCP > 0.1, for maximum of 5
  # if non of the variant has GWAS PIP/VCP > 0.1 we just show top one based on GWAS PIP/VCP and on xQTLPIP/VCP
  
  res_adxcub[,top_variants:=((max_variant_inclusion_probability>=0.1)|cV2F_rank<=5)|((max_variant_inclusion_probability_rank==1|variant_rank_xqtl==1|(rank(pval)==1&!is.na(pval)))),by='ADlocus']
  res_adxcubf<-res_adxcub[(top_variants)]
  
  #add supplemental cols
  res_adxcubf[,gwas_significance:=ifelse(min_pval<5e-8,'genome wide',
                                         ifelse(min_pval<1e-6,
                                                'suggestive',
                                                'ns'))]
  res_adxcubf[,mlog10pval:=-log10(min_pval)]
  
  #variant inclusion top confidence level
  res_adxcubf[,top_confidence:=str_replace(str_extract(xQTL_effects,'(?:CL|T)[0-9](?=,n=)'),'CL','T')]
  res_adxcubf$top_confidence|>unique()
  
  wb<-CreateExcelFormat(res_adxcubf,columns_mtd =cols,colors = colorsmtd,
                        wb = wb,sheet_name = cont)
  
  
}


saveWorkbook(wb, fp(out, paste0('unified_AD_loci_xQTL_summary_', STAMP, '.xlsx')), overwrite = TRUE)


