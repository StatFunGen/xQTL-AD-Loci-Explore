
fp<-function(...)file.path(...)
ps<-function(...)paste0(...)
library(data.table)
setDTthreads(threads = 0)
message(getDTthreads(),' threads available for data.table')

library(stringr)
library(ggplot2)
library(tidyverse)

get_gene_id <- function(gene_name){
  gene_ref <- read.table('/data/resource/references//Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list', header = F)
  gene_name = gene_name
  gene_id = gene_ref %>% filter(V5 == gene_name) %>% pull(V4)
  return(gene_id)
}


get_gene_info <- function(gene_id = NULL, gene_name = NULL){
  if(is.null(gene_id)) gene_id <- get_gene_id(gene_name)
  FunGen_gene <- fread("/data/interactive_analysis/rf2872/resource/Gene_Block.tsv")
  gene_info <- FunGen_gene %>% filter(region_id == gene_id) %>% mutate(gene_name = gene_name)
  target_LD_ids <- gene_info%>% pull(LD_matrix_id) %>% str_split(., ',', simplify = T) %>% unlist
  target_TAD_ids <- gene_info%>% pull(TADB_id) %>% str_split(., ',', simplify = T) %>% unlist
  target_sums_ids <- gene_info %>% pull(LD_sumstats_id) %>% str_split(., ',', simplify = T) %>% unlist
  gene_region <- paste0(gene_info[['#chr']], ":", gene_info[['start']], "-", gene_info[['end']])
  return(list(gene_info = gene_info, target_LD_ids = target_LD_ids, target_sums_ids = target_sums_ids, gene_region = gene_region, target_TAD_ids = target_TAD_ids))
}

getPval<-function(zscore){
  log_p_one_tail <- pnorm(abs(zscore), lower.tail = FALSE, log.p = TRUE)
  log_p_two_tail <- log(2) + log_p_one_tail  # log(2 * p)
  p_value_two_tail <- exp(log_p_two_tail)
  return(p_value_two_tail)
}

#genomics coordinates manipulation ####
start<-function(x,start_pos=2)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][start_pos]))
end<-function(x,end_pos=3)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][end_pos]))

seqid<-function(x,only_num=FALSE){
  if(only_num){
    str_extract(x,'[0-9]+')|>as.numeric()
  }else{
    sapply(x,function(x)strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][1])
    
  }
}

pos<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][2]))


ref<-function(x)sapply(x,function(x){
  vec=strsplit(x,"-|:|_|,|\\[|\\]")[[1]]
  return(vec[length(vec)-1])
})

alt<-function(x)sapply(x,function(x){
  vec=strsplit(x,"-|:|_|,|\\[|\\]")[[1]]
  return(vec[length(vec)])
})

#get Unified locis
#loci: data.table with 'variant_id' and $locus_col
#locus_col: column in data.table containing the locus id to overlap.
#overlap_col: if already computed, column in data.table containing the locus_id
#group.by: by which column(s) to separate unification of loci
UnifyLoci<-function(loci,variant_col='variant_id',
                    locus_col='locuscontext_id',group.by=NULL,
                    order.by='pos',order.by.ascending=TRUE,
                    overlap_col=NULL,rm_overlap=TRUE){
  if(variant_col!='variant_id'){
    loci$variant_id<-loci[[variant_col]]
    
  }
  
  if(locus_col!='locus_id'){
    loci$locus_id<-loci[[locus_col]]
    
  }
  
  
  
  if(is.null(overlap_col)){
    loci[,overlapping_locus.var:=paste(unique(locus_id),collapse ='|'),by='variant_id']
    
    loci[,overlapping_locus:=paste(setdiff(unlist(strsplit(overlapping_locus.var,'\\|')),locus_id),collapse = '|'),by=locus_col]
    
    
  }else{
    loci$overlapping_locus<-loci[[overlap_col]]
    
  }
  loci_u<-unique(loci,by=locus_col)
  if(order.by=='pos'){
    loci_u<-loci_u[order(seqid(variant_id),pos(variant_id))]
  }else{
    setorderv(loci_u,cols = order.by,order = as.numeric(order.by.ascending))
    
  }
  
  
  if(!is.null(group.by)){
    loci_u[,group:=apply(.SD,1,function(x)paste(x,collapse = '.')),.SDcols=group.by]
    
  }
  
  loci_u[,uni.locus_id:={
    if(!is.null(group.by)){
      pref.loc=unique(group)
    }else{
      pref.loc='loc'
    }
    
    #get the different sets list
    #set number back to the locus
    #for each gwasoverlap, if intersect, i create the union 
    overlaps_list<-strsplit(overlapping_locus,'\\|')
    overlaps_list<-lapply(1:.N, function(i)c(overlaps_list[[i]],locus_id[i]))
    overlap_set_num<-1:.N
    already_set<-c()
    for(i in 1:.N){
      if(!i%in%already_set){
        overlap_set_num[i]<-i
        already_set<-c(already_set,i)
      }
      if(.N>1){
        for(j in 2:.N){
          
          if(length(intersect(overlaps_list[[overlap_set_num[i]]],
                              overlaps_list[[overlap_set_num[j]]]))>0){
            overlaps_list[[overlap_set_num[i]]]<-union(overlaps_list[[overlap_set_num[i]]],overlaps_list[[overlap_set_num[j]]])
            
            overlap_set_num[j]<-overlap_set_num[i]
            already_set<-c(already_set,j)
          }
          
        }
      }
      
    }
    paste(pref.loc,overlap_set_num,sep='_')
    
  },by=group.by]
  
  loci<-merge(loci,loci_u[,.SD,.SDcols=c(locus_col,group.by,'uni.locus_id')],by=c(locus_col,group.by))
  if(rm_overlap){
    loci<-loci[,-c('overlapping_locus.var','overlapping_locus')]
    if(!is.null(overlap_col)){
      loci[[overlap_col]]<-NULL
    }
  }
  return(loci)
}



#VariantsLevelSummary####
#annotation AD locus overlapped table with confidence score and summarize 
SummarizeTable<-function(res_adx,
                         group.by='context_short',
                         xqtl_methods=c('fSuSiE_finemapping','single_context_finemapping',
                                        'multi_context_finemapping','AD_xQTL_colocalization',
                                        'TWAS/MR','cTWAS','Coloc')){
  #res_adx<-res_adx[!str_detect(context,'^AD')]
  
  
  #COntexts summary
  ##for each locus, add confidence level of each broad context
  #C1: cTWAS/MR, cs95 single context fine-mapping overlap
  #C2: cTWAS/MR, colocalization
  #C3: TWAS and (cs95 fine-mapping overlap OR colocalization)
  #C4: cs95 single context fine-mapping overlap
  #C5: colocalization OR any fine-mapping overlap (multicontext, cs50, cs70..)
  #C6: TWAS only
  message('summarizing per ',group.by,' the xQTL evidence for each variant')
  res_adx[!str_detect(context,'^AD'),confidence_lvl:={
    if((any(MR_signif,na.rm = T)|any(cTWAS_signif,na.rm = T))&(any(Method%in%c('single_context_finemapping','fSuSiE_finemapping')&str_detect(credibleset,'cs95')&susie_coverage%in%c('cs95','cs70'),na.rm = T))){
      'T1'
    }else if((any(MR_signif,na.rm = T)|any(cTWAS_signif,na.rm = T))&any(Method=='AD_xQTL_colocalization',na.rm = T)){
      'T2'
    }else if(any(TWAS_signif,na.rm = T)&(any(Method%in%c('single_context_finemapping','fSuSiE_finemapping')&str_detect(credibleset,'cs95'),na.rm = T)|any(Method=='AD_xQTL_colocalization',na.rm = T))){
      'T3'
    }else if(any(Method%in%c('single_context_finemapping','fSuSiE_finemapping')&str_detect(credibleset,'cs95')&susie_coverage%in%c('cs95','cs70'),na.rm = T)){
      'T4'
    }else if(any(Method%in%c('AD_xQTL_colocalization','multi_context_finemapping',
                             'single_context_finemapping','fSuSiE_finemapping','sn_sQTL','Coloc'),na.rm = T)){
      'T5'
    }else{
      'T6'
    }
  },
  by=c('variant_ID','gene_ID',group.by)]
  
  res_adx[,confidence_lvl:=factor(confidence_lvl,levels = paste0('T',1:6))]
  
  if('xQTL_effect'%in%colnames(res_adx)){
    res_adx[,xQTL_effect:=NULL]
    
  }
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_effect:=paste0(get(group.by),'.')]
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_effect:=ifelse(is.na(conditional_effect),xQTL_effect,ifelse(conditional_effect>0,paste0(get(group.by),'+'),paste0(get(group.by),'-')))[order(is.na(conditional_effect))][1],by=c('variant_ID','gene_ID',group.by)]
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_effect:=ifelse(is.na(coef),xQTL_effect,paste0(xQTL_effect,ifelse(coef>0,'+','-')))[order(is.na(coef))][1],by=c('variant_ID','gene_ID',group.by)]
  
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),n_study:=length(unique(context)),by=c('ADlocus','gene_ID',group.by)]
  
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_effect:=paste0(xQTL_effect,' (',confidence_lvl,',n=',n_study,')'),
          by=c('variant_ID','gene_ID',group.by)]
  res_adx[order(variant_ID,confidence_lvl,-abs(twas_z),-cos_npc),xQTL_effects:=paste(unique(xQTL_effect[!is.na(xQTL_effect)]),collapse =  '|'),by=c('variant_ID','gene_ID')]
  
  #add ncontexts and methods by variant
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_methods:=paste(unique(Method[order(-ifelse(is.na(vcp),PIP,vcp))]),collapse =  '|'),by=c('variant_ID','gene_ID')]
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),xQTL_contexts:=paste(unique(context[order(-ifelse(is.na(vcp),PIP,vcp))]),collapse =  '|'),by=c('variant_ID','gene_ID')]
  
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),n_contexts:=length(unique(context)),by=c('variant_ID','gene_ID')]
  
  res_adx[!str_detect(context,'^AD')&Method%in%c(xqtl_methods,'sn_sQTL'),n_methods:=length(unique(Method)),by=c('variant_ID','gene_ID')]
  
  
  #Add top twas zscore per gene
  res_adx[!str_detect(context,'^AD'),TWAS_signif_gene:=any(TWAS_signif,na.rm = T),by=c('gene_ID')]
  res_adx[!str_detect(context,'^AD'),MR_signif_gene:=any(MR_signif,na.rm = T),by=c('gene_ID')]
  res_adx[!str_detect(context,'^AD'),cTWAS_signif_gene:=any(cTWAS_signif,na.rm = T),by=c('gene_ID')]
  
  
  res_adx[!str_detect(context,'^AD'),twas_pval_gene_min:=twas_pval[which.min(twas_pval)][1],by=.(gene_ID)]
  res_adx[!str_detect(context,'^AD'),twas_z_gene_max:=twas_z[which.min(twas_pval)][1],by=.(gene_ID)]
  res_adx[!str_detect(context,'^AD'),twas_z_gene_max_context:=context[which.min(twas_pval)][1],by=.(gene_ID)]
  
  #xQTL level variant prioritization
  #MAX_VIP for xQTL like for GWAS 
  #i.e. max_vip, vip_rank, top_method, top_context pER GENE
  
  res_adx[,variant_inclusion_probability_xqtl:=max(ifelse(is.na(vcp),PIP,vcp)[Method%in%xqtl_methods&!str_detect(context,'^AD')],na.rm = T),by=c('variant_ID','gene_ID',group.by)]
  res_adx[is.infinite(variant_inclusion_probability_xqtl),variant_inclusion_probability_xqtl:=NA]
  #add variant rank
  res_adxv<-unique(res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD')][order(confidence_lvl,-variant_inclusion_probability_xqtl)],by=c('variant_ID','gene_ID',group.by))
  res_adxv[,variant_rank_xqtl:=rank(-variant_inclusion_probability_xqtl),
           by=c('ADlocus','gene_ID',group.by)]
  
  
  res_adx<-merge(res_adx[,-c('variant_rank_xqtl')],res_adxv[,.SD,.SDcols=c('variant_ID','ADlocus','gene_ID',group.by,'variant_rank_xqtl')],
                 by=c('ADlocus','gene_ID',group.by,'variant_ID'),all.x=T)
  #add top context associated
  res_adx[,method_xqtl:=Method[which.max(variant_inclusion_probability_xqtl)],by=c('variant_ID','gene_ID',group.by)]
  res_adx[,context_xqtl:=context[which.max(variant_inclusion_probability_xqtl)],by=c('variant_ID','gene_ID',group.by)]
  res_adx[,n.variant_xqtl:=n.variant.locus[which.max(variant_inclusion_probability_xqtl)],by=c('variant_ID','gene_ID',group.by)]
  res_adx[,effect_xqtl:=ifelse(str_detect(Method,'fSuSiE'),top_effect,ifelse(str_detect(Method,'finemapping'),conditional_effect,coef))[which.max(variant_inclusion_probability_xqtl)],by=c('variant_ID','gene_ID',group.by)]
  
  res_adx[,coverage_xqtl:=str_extract(credibleset,'cs[0-9]+')[which.max(PIP)],by=c('variant_ID','gene_ID',group.by)]
  
  
  #add transQTL infos
  #standard
  res_adx[,have_trans_effect:=any(str_detect(Method,'trans')),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),n_trans_genes:=length(unique(gene_name[Method=='trans_finemapping'&gene_name!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_genes:=paste(unique(gene_name[Method=='trans_finemapping'&gene_name!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),n_trans_contexts:=length(unique(context[Method=='trans_finemapping'&context!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts:=paste(unique(context[Method=='trans_finemapping'&context!='']),collapse = '|'),by=.(variant_ID)]
  
  
   #others trans
  #snuc
  res_adx[order(-have_trans_effect,-PIP),n_trans_genes_snRNA:=length(unique(gene_name[Method=='trans_single_gene_snRNA'&gene_name!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_genes_snRNA:=paste(unique(gene_name[Method=='trans_single_gene_snRNA'&gene_name!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts_snRNA:=paste(unique(context[Method=='trans_single_gene_snRNA'&context!='']),collapse = '|'),by=.(variant_ID)]
  
  #pqtl
  res_adx[order(-have_trans_effect,-PIP),n_trans_genes_pQTL:=length(unique(gene_name[Method=='trans_pQTL'&gene_name!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_genes_pQTL:=paste(unique(gene_name[Method=='trans_pQTL'&gene_name!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts_pQTL:=paste(unique(context[Method=='trans_pQTL'&context!='']),collapse = '|'),by=.(variant_ID)]
  
  #gpqtl
  res_adx[order(-have_trans_effect,-PIP),n_trans_genes_gpQTL:=length(unique(gene_name[Method=='trans_gpQTL'&gene_name!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_genes_gpQTL:=paste(unique(gene_name[Method=='trans_gpQTL'&gene_name!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts_gpQTL:=paste(unique(context[Method=='trans_gpQTL'&context!='']),collapse = '|'),by=.(variant_ID)]
  
  #programs based
  res_adx[Method=='trans_cca']$event_ID
  res_adx[Method=='trans_cca',program:=str_extract(event_ID,'Program_[0-9]+')]
  res_adx[order(-have_trans_effect,-PIP),n_trans_cca_programs:=length(unique(program[Method=='trans_cca'&program!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_cca_programs:=paste(unique(program[Method=='trans_cca'&program!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts_cca_programs:=paste(unique(context[Method=='trans_cca'&context!='']),collapse = '|'),by=.(variant_ID)]
  
  res_adx[Method=='trans_hotspot']$event_ID
  res_adx[Method=='trans_hotspot',program:=paste(str_extract(event_ID,'minNumGene[0-9]+'),str_extract(event_ID,'module_[0-9]+_WPC[0-9]+'),sep='_')]
  res_adx[order(-have_trans_effect,-PIP),n_trans_hotspot_programs:=length(unique(program[Method=='trans_hotspot'&program!=''])),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_hotspot_programs:=paste(unique(program[Method=='trans_hotspot'&program!='']),collapse = '|'),by=.(variant_ID)]
  res_adx[order(-have_trans_effect,-PIP),trans_contexts_hotspot_programs:=paste(unique(context[Method=='trans_hotspot'&context!='']),collapse = '|'),by=.(variant_ID)]

  
  #mv genes
  res_adx[!str_detect(context,'^AD'),mv_genes:=paste(unique(gene_name[Method=='multi_gene_finemapping'&!gene_name=='']),collapse = '|'),by=.(variant_ID,ADlocus)]
  res_adx[!str_detect(context,'^AD'),n.mv_genes:=length(unique(gene_name[Method=='multi_gene_finemapping'&!gene_name==''])),by=.(variant_ID,ADlocus)]
  
  #QR 
  res_adx[order(p_bonferroni_adj),qr_contexts:=paste(unique(context[Method=='QR'&context!='']),collapse = '|'),by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_qr_context:=context[which.min(p_bonferroni_adj)][1],by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_qr_classification:=classification[which.min(p_bonferroni_adj)][1],by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_qr_p_bonferroni_adj:=min(p_bonferroni_adj,na.rm = T),by=.(variant_ID,gene_ID)]
  
  
  #interactions
  res_adx[order(pvalue_APOE_interaction),apoe4_interactions:=paste(unique(context[Method=='APOE interaction'&context!='']),collapse = '|'),by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_apoe4_context:=context[which.min(pvalue_APOE_interaction)][1],by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_apoe4_p_value:=min(pvalue_APOE_interaction,na.rm = T),by=.(variant_ID,gene_ID)]
  
  res_adx[order(pvalue_msex_interaction),msex_interactions:=paste(unique(context[Method=='msex interaction'&context!='']),collapse = '|'),by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_msex_context:=context[which.min(pvalue_msex_interaction)][1],by=.(variant_ID,gene_ID)]
  res_adx[!str_detect(context,'^AD'),top_msex_p_value:=min(pvalue_msex_interaction,na.rm = T),by=.(variant_ID,gene_ID)]
  
  return(res_adx)
  
  
}

#WideTable####
#From AD locus overlapped long table of  variant_ID level Method/context/gene_ID summ stats, create a wide table at variant-gene level adding columns information per 'split.by' level 

WideTable<-function(res_adx,split.by=c('context_broad2','qtl_type')){
  
  #rm MEthpd= colocnooost and context = GWAS 
  res_adx<-res_adx[!(Method=='AD_xQTL_colocalization'&str_detect(context,'AD'))]
  #   Add per broad context infos
  # separate by main context (cell type or tissue)
  message('generating columns per ',paste(split.by,collapse = ' and '))
  
  res_adxl<-split(res_adx[!context_broad2%in%c('','GWAS')&!is.na(context_broad2)],by = split.by)
  
  #for each add inclusion prob, Fp, coloc, multicontext,datasets, caQTL, QR, Interaction
  
  res_adxlu<-lapply(res_adxl, function(res){
    message(unique(res[,paste(get(split.by[1]),get(split.by[2]))]))
    #add for top inclusion score, context, fp pip, fp pip effect,  vcp, coloc effect, multicontext pip, TWAS Z, MR T/F, QR Pval, APOE4 Pval, Sex Pval, 
    res[,incl_prob:=ifelse(is.na(vcp),PIP,vcp)]
    res[order(-incl_prob),max_dataset:=context[1],by=c('variant_ID','gene_ID')]
    
    res[,max_inclusion_probability:=max(incl_prob,na.rm = T),by=c('variant_ID','gene_ID')]
    res[is.infinite(max_inclusion_probability),max_inclusion_probability:=NA]
    
    
    
    res[,max_finemap_pip:=PIP[Method%in%c('single_context_finemapping')][context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    res[,max_finemap_effect:=conditional_effect[Method%in%c('single_context_finemapping')][context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    
    
    res[,max_coloc_vcp:=vcp[Method=='AD_xQTL_colocalization'][context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    res[,max_coloc_effect:=coef[Method=='AD_xQTL_colocalization'][context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    
    
    res[,max_multicontext_pip:=PIP[Method=='multi_context_finemapping'][context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    res[,max_multicontext_effect:=conditional_effect[Method=='multi_context_finemapping'&context==max_dataset[1]][1],by=c('variant_ID','gene_ID')]
    
    
    res[,max_TWAS_Z:=twas_z[which.max(abs(twas_z[context==max_dataset[1]]))],by=c('gene_ID')]
    
    res[,max_MR:=any(MR_signif[context==max_dataset[1]]),by=c('gene_ID')]
    
    res[,max_QR_padj:=p_bonferroni_adj[Method=='QR'&context==max_dataset[1]]|>min(na.rm = T),by=c('variant_ID','gene_ID')]
    res[is.infinite(max_QR_padj),max_QR_padj:=NA]
    
    #res[,QR_type:=classification[!is.na(classification)][1],by=c('variant_ID','gene_ID')]
    
    res[,max_APOE4_interaction_pval:=pvalue_APOE_interaction[Method=='APOE interaction'&context==max_dataset[1]]|>min(na.rm = T),by=c('variant_ID','gene_ID')]
    res[,max_msex_interaction_pval:=pvalue_msex_interaction[Method=='msex interaction'&context==max_dataset[1]]|>min(na.rm = T),by=c('variant_ID','gene_ID')]
    res[is.infinite(max_APOE4_interaction_pval),max_APOE4_interaction_pval:=NA]
    res[is.infinite(max_msex_interaction_pval),max_msex_interaction_pval:=NA]
    
    
    #add QTL type specific info
    if('sn_sQTL'%in%unique(res$Method)){
      #sn sQTL spe info beta, pval, splice site
      #FIXME IF more than one context to aggregate
      res[,nsplice_sites:=length(unique(splice_site)),by=c('variant_ID','gene_ID')]
      res[order(pval,-abs(beta)),
          max_splice_sites:=paste(unique(splice_site),collapse = '|'),by=c('variant_ID','gene_ID')]
      
      res[order(pval,-abs(beta)),
          max_splice_types:=paste(unique(splice_type[!duplicated(splice_site)]),collapse = '|'),by=c('variant_ID','gene_ID')]
      
      res[order(pval,-abs(beta)),max_sn_splice_pval:=pval[!duplicated(splice_site)][1],
          by=c('variant_ID','gene_ID')]
      res[order(pval,-abs(beta)),max_sn_splice_beta:=beta[!duplicated(splice_site)][1],
          by=c('variant_ID','gene_ID')]
      
      res<-res[,-c('max_finemap_pip','max_inclusion_probability','max_finemap_effect',
                   'max_coloc_vcp','max_coloc_effect',
                   'max_multicontext_pip','max_multicontext_effect',
                   'max_TWAS_Z','max_MR','max_QR_pval',
                   'max_APOE4_interaction_pval','max_msex_interaction_pval')]
      
      
    }
    
    
    if(any(str_detect(unique(res$context_short),'bulk (u-|s-)?sQTL'))){
      res[,max_njunctions:=length(unique(event_ID[context==max_dataset[1]][1])),by=c('variant_ID','gene_ID')]
      res[order(-ifelse(is.na(vcp),PIP,vcp),-ifelse(is.na(coef),conditional_effect,coef)),
          max_junctions:=paste(unique(str_extract(event_ID,'chr[0-9]+:[0-9]+:[0-9]+')[context==max_dataset[1]][1]),collapse = '|'),by=c('variant_ID','gene_ID')]
      res[order(-ifelse(is.na(vcp),PIP,vcp),-ifelse(is.na(coef),abs(conditional_effect),abs(coef))),
          max_junctions_type:=paste(str_extract(event_ID,'(?<=:)NE|UP|PR(?=:)')[context==max_dataset[1]][1][!duplicated(str_extract(event_ID,'chr[0-9]+:[0-9]+:[0-9]+'))],collapse = '|'),by=c('variant_ID','gene_ID')]
      
    }
    
    
    if(unique(res$qtl_type)=='gpQTL'){
      res[order(-ifelse(is.na(vcp),PIP,vcp)),glyco_changes:=paste(unique(str_extract(event_ID[context==max_dataset[1]][1],'(un)adjusted_gp_[0-9]+')),collapse = '|'),by=c('variant_ID','gene_ID')]
      
    }
    
    if('fSuSiE_finemapping'%in%unique(res$Method)){
      
      #peak position
      res[,top_epi_mark:=strsplit(epi_mark_names[context==max_dataset[1]][1],';')[[1]][top_effect_coord[context==max_dataset[1]][1][1]],
          by=c('variant_ID','gene_ID')]
      res[,top_epi_mark_position:=strsplit(epi_mark_positions[context==max_dataset[1]][1],';')[[1]][top_effect_coord[context==max_dataset[1]][1][1]],
          by=c('variant_ID','gene_ID')]
      
      
      res<-res[,-c('max_coloc_vcp','max_coloc_effect',
                   'max_multicontext_pip','max_multicontext_effect',
                   'max_TWAS_Z','max_MR','max_QR_pval',
                   'max_APOE4_interaction_pval','max_msex_interaction_pval')]
    }else{
      res<-res[,-c('effect_peak_start','effect_peak_end','effect_peak_index','top_effect')]
    }
    
    
    #then n tot datasets, datasets.
    res[,n_datasets:=length(unique(context)),by=c('variant_ID','gene_ID')]
    
    res[order(-ifelse(is.na(vcp),PIP,vcp),-ifelse(is.na(coef),conditional_effect,coef)),datasets:=paste(unique(context),collapse = '|'),by=c('variant_ID','gene_ID')]
    
    
    
    return(unique(res,by=c('variant_ID','gene_ID')))
    
  })
  #add context suffixes
  contexts_cols=c('variant_ID','gene_ID',setdiff(Reduce(union,lapply(res_adxlu,colnames)),colnames(res_adx)))
  
  res_adx[variant_ID=='chr7:143414019:C:G'][str_detect(Method,'fSuSiE')]
  
  res_adxlu<-lapply(names(res_adxlu), function(c){
    cols_to_keep<-intersect(contexts_cols,colnames(res_adxlu[[c]]))
    to_rename<-setdiff(cols_to_keep,c('variant_ID','gene_ID'))
    newnames<-paste0(to_rename,'.',c)
    res<-res_adxlu[[c]][,..cols_to_keep]
    setnames(res,to_rename,newnames)
    return(res)
  })
  
  
  #Get the variant level summary table####
  #ORDER: Order locus per gene and variant importance, keeping the top context for each variant
  
  res_adxu<-unique(res_adx[!str_detect(Method,'trans')][order(locus_index,
                                                              confidence_lvl,n_study,
                                                              -abs(twas_z_gene_max),
                                                              cV2F_rank,variant_rank_xqtl)],
                   by=c('locus_index','gene_ID','variant_ID'))
  #remove row if variant without gene while variant  present with a gene
  res_adxu<-res_adxu[order(locus_index,-gene_ID)][!(duplicated(variant_ID)&(gene_ID==''|is.na(gene_ID)))]
  #merge with the context specific columns
  
  res_adxub<-Reduce(function(x,y)merge(x,y,by=c('variant_ID','gene_ID'),all=T),c(list(res_adxu),
                                                                                 res_adxlu))[order(locus_index,
                                                                                                   confidence_lvl,n_study,
                                                                                                   -abs(twas_z_gene_max),
                                                                                                   cV2F_rank,variant_rank_xqtl)]
  
  #add supplemental cols
  res_adxub[,gwas_significance:=ifelse(min_pval<5e-8,'genome wide',
                                       ifelse(min_pval<1e-6,
                                              'suggestive',
                                              'ns'))]
  res_adxub[,mlog10pval:=-log10(min_pval)]
  
  #variant inclusion top confidence level
  res_adxub[,top_confidence:=str_replace(str_extract(xQTL_effects,'(?:CL|T)[0-9](?=,n=)'),'CL','T')]

  ## ---- FIX 2026-09-17 gene-level confidence -------------------------------
  ## Gene-level evidence (TWAS/MR/cTWAS, trans) carries no variant_ID, so the
  ## merge above leaves it on its own row (variant_ID==''), while the gene's real
  ## variant rows have xQTL_effects/top_confidence NA. Downstream filters drop the
  ## blank-variant row, silently losing the gene's confidence (e.g. APOE CL6).
  ## Carry the gene's best xQTL_effects/top_confidence onto its variant rows in
  ## the same locus when those are blank. Rows that already have a CL keep it.
  .blank <- function(v) is.na(v) | v==''
  res_adxub[, `:=`(
      .g_eff = { k <- which(!.blank(top_confidence))
                 if (length(k)) { j <- k[which.min(as.integer(str_extract(top_confidence[k],'[0-9]')))]
                                  xQTL_effects[j] } else NA_character_ },
      .g_tc  = { k <- which(!.blank(top_confidence))
                 if (length(k)) top_confidence[k][which.min(as.integer(str_extract(top_confidence[k],'[0-9]')))]
                 else NA_character_ }),
    by = .(locus_index, gene_ID)]
  .n_fill <- res_adxub[.blank(top_confidence) & !is.na(.g_tc) & !.blank(variant_ID), .N]
  res_adxub[.blank(top_confidence) & !is.na(.g_tc) & !.blank(variant_ID),
            `:=`(xQTL_effects = .g_eff, top_confidence = .g_tc)]
  res_adxub[, c('.g_eff','.g_tc') := NULL]
  message('[WideTable fix] variant rows given gene-level confidence: ', .n_fill)
  ## -------------------------------------------------------------------------
  
  
  return(res_adxub)
}

#PrepColsMtd: use colsmtd to duplicate the columns metadata of cols according to the wildcard in parent and grandparentcolumn
PrepColsMtd<-function(cols,colsmtd,res){
  #replicate columns for each $variable adding the varnames
  variable_gp<-str_extract(cols$grandparent_column,'\\$[A-Za-z0-9_.]+')|>na.omit()|>str_remove('\\$')|>as.vector()
  
  variable_p<-str_extract(cols$parent_column,'\\$[A-Za-z0-9_.]+')|>na.omit()|>str_remove('\\$')|>as.vector()
  
  variable_gp_vec<-colsmtd[wildcard==variable_gp]$r_name
  
  variable_p_vec<-colsmtd[wildcard==variable_p]$r_name
  
  cols_extended<-rbindlist(lapply(variable_gp_vec,function(v1){
    color=colsmtd[wildcard==variable_gp][r_name==v1]$color
    grandparentname<-colsmtd[wildcard==variable_gp][r_name==v1]$excel_name
    
    cols_ext<-rbindlist(lapply(variable_p_vec,function(v2){
      name=paste(v1,v2,sep = '.')
      
      return(cols[expanded==1][,r_name:=paste(r_name,name,sep = '.')][,parent_column:=v2])
      
    }))
    cols_ext[,parent_fill:=color]
    cols_ext[,grandparent_column:=grandparentname]
    
    return(cols_ext)
  }))
  cols=rbind(cols[expanded==0],cols_extended)
  
  
  
  return(cols)
}



#create formatted Excel Sheet####
#res_summ: the summary table to transform as formatted excel sheet
#columns_mtd: metadata for the columns to keep on this table, need at least 'r_name', 'excel_name'	'parent_column'	'grandparent_column'	'parent_fill'	'width'	'to_divide'	'wrap_text'
#colors: patterns to colors cells or text. must have 'pattern','fill_color' and 'font_color'.
CreateExcelFormat<-function(res_summ,columns_mtd,colors,wb=NULL,
                            file=NULL,
                            sheet_name="Gene Locus table",
                            start_rheader=3){
  require(openxlsx)
  # #rm NA for epi QTL
  # cols_with_nachar<-colnames(res_summ)[sapply(res_summ, function(x)any(is.character(x)&str_detect(x,'NA')))]
  # if(length(cols_with_nachar)>0){
  #   message('NA detected to remove in ',paste(cols_with_nachar,collapse = ', '))
  #   res_summ[,(cols_with_nachar):=lapply(.SD,function(x)str_remove_all(x,'NA')),.SDcols=cols_with_nachar]
  #   
  # }
  # 
  #remove columns all empty
  emptycols<-colnames(res_summ)[sapply(res_summ, function(x)all(is.na(x)|x==''))]
  if(length(emptycols)>0){
    message(paste(emptycols,sep = ', '),' is empty, removing')
    res_summ<-res_summ[,.SD,.SDcols=!emptycols]
  }
  
  #keep only columns in metadata and in table
  
  cols_to_keep<-intersect(columns_mtd$r_name,colnames(res_summ))
  cols_absent<-setdiff(columns_mtd$r_name,colnames(res_summ))
  if(length(cols_absent)>0){
    columns_mtd<-columns_mtd[r_name%in%cols_to_keep]
    
    if(length(cols_absent)>6){
      message(paste(head(cols_absent),collapse = ', '),' and ',length(cols_absent)-6 ,' others are not present in the summary table')
      
    }else{
      message(paste(cols_absent,collapse = ', '),' are not present in the summary table')
      
    }
    
  }
  
  res_summ<-copy(res_summ)[,..cols_to_keep][]|>unique()
  
  
  #split columns with '|' separator according to metadata
  cols_to_divide<-columns_mtd[(to_divide==1)]$r_name
  message('dividing columns: ',paste(cols_to_divide,collapse = ', '))
  
  cols_to_create_list<-lapply(cols_to_divide, function(col){
    n_col_to_add<-res_summ[[col]]|>str_count('\\|')|>max(na.rm = T)
    message('dividing ',col,' into ',n_col_to_add+1,' columns')
    if(n_col_to_add>0){
      c(col,paste0(col,'_',1:n_col_to_add))
      
    }else{
      col
    }
  })
  names(cols_to_create_list)<-cols_to_divide
  # create the splitted cols in the metadata and the table
  for(col in names(cols_to_create_list)[sapply(cols_to_create_list,function(x)length(x)>1)]){
    
    tocreate=cols_to_create_list[[col]]
    res_summ[,(tocreate):=str_split(.SD[[1]],'\\|',simplify = T)|>as.data.table(),.SDcols=col]
    
    columns_mtd<-rbindlist(list(columns_mtd[1:which(columns_mtd$r_name==col) ],
                                rbindlist(lapply(cols_to_create_list[[col]],function(newcol)columns_mtd[col,on='r_name' ][,r_name:=newcol])),
                                columns_mtd[(which(columns_mtd$r_name==col)+1):nrow(columns_mtd) ]))
    
  }
  columns_mtd<-unique(columns_mtd)[!is.na(r_name)&r_name!=''&!duplicated(r_name)]
  
  #columns ordering accroding to mtd
  setcolorder(res_summ,neworder =columns_mtd$r_name)
  
  
  
  
  #transform '|' as ', ' separator for columns to wrap text
  for(col in columns_mtd[wrap_text==1]$r_name){
    res_summ[,(col):=str_replace_all(.SD[[1]],'\\|',', '),.SDcols=col]
  }
  
  
  
  setdiff(columns_mtd$r_name,colnames(res_summ))
  if(is.null(wb)){
    wb <- createWorkbook()
  }
  
  #add Sheet
  addWorksheet(wb,sheet_name)
  
  writeData(wb, sheet_name, res_summ, startCol = 1, startRow = start_rheader)
  
  #merge cells column name if from the columns to divide
  for(col in names(cols_to_create_list)){
    cols_pos<-which(colnames(res_summ)%in%cols_to_create_list[[col]])
    mergeCells(wb,sheet_name, cols = cols_pos, rows = start_rheader) 
  }
  
  #add conditional formatting
  for(context in unique(colors$pattern)){
    for(col in names(cols_to_create_list)){
      cols_pos<-which(colnames(res_summ)%in%cols_to_create_list[[col]])
      conditionalFormatting(wb, sheet_name,
                            cols = cols_pos,
                            rows = (start_rheader+1):(nrow(res_summ)+(start_rheader+1)), 
                            rule = context, 
                            type=ifelse('type'%in%colnames(colors),colors[pattern==context]$type[1],'beginsWith'),
                            style = createStyle(bgFill = colors[pattern==context]$fill_color[1],
                                                fontColour = colors[pattern==context]$font_color[1],
                                                fontSize = 10,
                                                border='TopBottomLeftRight'))  
      
    }
  }
  
  # Define style for odd numbers 
  oddStyle <- createStyle(fontColour = "#FFFFFF", 
                          bgFill = "#7b848a",
                          fontSize = 10,
                          border='TopBottomLeftRight')  
  
  
  #odd number
  conditionalFormatting(
    wb, sheet_name,
    cols = which(colnames(res_summ)=='locus_index'),
    rows = (start_rheader+1):(nrow(res_summ)+(start_rheader+1)),
    rule = "MOD(A4,2)=1",
    style = oddStyle,
    type = "expression"
  )
  
  #COLORING of content
  #color in green TRUE for colones with logical
  greenStyle <- createStyle(fontColour = "#0d944a", bgFill = "#cef5df")  # Green
  
  logicals<-colnames(res_summ)[sapply(res_summ,function(x)is.logical(x)|any(str_detect(x,'TRUE'),na.rm = T))]
  
  for(col in logicals){
    colnum<-which(colnames(res_summ)==col)
    col_letter<-int2col(colnum)
    conditionalFormatting(
      wb, sheet_name,
      cols =colnum,
      rows = (start_rheader+1):(nrow(res_summ)+(start_rheader+1)),
      rule = paste0(col_letter,"4=TRUE"),
      style = greenStyle,
      type = "expression"
    )
  }
  
  
  #add color scale
  cols_to_colorscales<-which(colnames(res_summ)%in%  columns_mtd[to_colorscale==1]$r_name)
  
  # Apply color gradient on Score (column E)
  for(col in cols_to_colorscales){
    if(all(c(1,-1)%in%unique(sign(res_summ[[col]])))){
      
      conditionalFormatting(
        wb, sheet_name,
        cols = col, 
        rows = 2:(nrow(res_summ) + 1),
        type = "colourScale",
        style = c(columns_mtd[col,]$colorlow,columns_mtd[col,]$colormid,columns_mtd[col,]$colorhigh),
        rule = c(min(res_summ[[col]],na.rm = TRUE), 0, max(res_summ[[col]],na.rm = TRUE)) #
      )
    }else{
      conditionalFormatting(
        wb, sheet_name,
        cols = col, 
        rows = 2:(nrow(res_summ) + 1),
        type = "colourScale",
        style = c(columns_mtd[col,]$colorlow,columns_mtd[col,]$colormid,columns_mtd[col,]$colorhigh)
      )
    }
    
  }
  
  #add parent columns
  writeData(wb, sheet_name, str_remove(t(columns_mtd)[c(4,3,2),],'\\.')|>matrix(nrow = 3),
            startCol = 1, startRow = 1,
            colNames = F,rowNames = F)
  
  #merge cells
  sepas<-c(which(columns_mtd$grandparent_column!=''&!duplicated(columns_mtd$grandparent_column)),nrow(columns_mtd)+1)
  for(i in 1:(length(sepas)-1) ){
    
    mergeCells(wb,sheet_name, cols = c(sepas[i],sepas[i+1]-1),
               rows = 1) 
  }
  
  sepas<-c(which(columns_mtd$parent_column!=''&(!duplicated(columns_mtd[,paste(grandparent_column,parent_column)])|columns_mtd$parent_column=='.')),nrow(columns_mtd)+1)
  for(i in 1:(length(sepas)-1) ){
    
    mergeCells(wb,sheet_name, cols = c(sepas[i],sepas[i+1]-1),
               rows = 2) 
  }
  
  
  #Adjust Column Width and Row Height
  setRowHeights(wb, sheet_name, rows = 1:3,
                heights = c(40, 35,60))
  setRowHeights(wb, sheet_name, rows = (start_rheader+1):(nrow(res_summ)+(start_rheader+1)),
                heights = 37)
  
  setColWidths(wb,sheet_name,
               cols=which(columns_mtd$width!='auto'),
               widths =as.numeric(columns_mtd[width!='auto']$width),
               hidden = as.numeric(columns_mtd[width!='auto']$width)==0)
  
  setColWidths(wb,sheet_name,
               cols=which(columns_mtd$width=='auto'),widths ='auto')
  
  #format the table content
  #for numeric h center, for character left
  to_center<-which(sapply(res_summ,function(x)is.numeric(x)|all(str_length(x)<=5)))
  addStyle(wb, sheet_name, style = createStyle(wrapText = F,
                                               valign = 'center',
                                               halign = 'left',
                                               border ='TopBottomLeftRight' ),
           rows = 4:(nrow(res_summ)+4),
           cols=setdiff(1:ncol(res_summ),to_center),
           gridExpand = TRUE)
  
  addStyle(wb, sheet_name, style = createStyle(wrapText = F,
                                               valign = 'center',
                                               halign = 'center',
                                               border ='TopBottomLeftRight' ),
           rows = 4:(nrow(res_summ)+4),
           cols=to_center, gridExpand = TRUE)
  
  
  #format headers
  addStyle(wb, sheet_name, style = createStyle(wrapText = TRUE,
                                               valign = 'center',
                                               halign = 'center',
                                               textDecoration = 'bold',
                                               border ='TopBottomLeftRight'),
           rows = 3,cols=1:ncol(res_summ), 
           gridExpand = TRUE)
  
  sepas<-c(which(columns_mtd$parent_fill!=''),nrow(columns_mtd)+1)
  for(i in 1:(length(sepas)-1) ){
    
    addStyle(wb, sheet_name, style = createStyle(wrapText = TRUE,
                                                 fontSize = 14,
                                                 valign = 'center',
                                                 halign = 'center',
                                                 textDecoration = 'bold',
                                                 border ='TopBottomLeftRight',
                                                 fgFill =columns_mtd[sepas[i],]$parent_fill ),
             rows = 1,cols=sepas[i]:(sepas[i+1]-1))
    addStyle(wb, sheet_name, style = createStyle(wrapText = TRUE,
                                                 fontSize = 12,
                                                 valign = 'center',
                                                 halign = 'center',
                                                 textDecoration = 'bold',
                                                 border ='TopBottomLeftRight',
                                                 fgFill =columns_mtd[sepas[i],]$parent_fill ),
             rows = 2,
             cols=sepas[i]:(sepas[i+1]-1))
  }
  
  
  #reformat specific columns
  addStyle(wb, sheet_name, style = createStyle(wrapText = T,
                                               valign = 'center',
                                               halign = 'center',
                                               border ='TopBottomLeftRight' ),
           rows = 4:(nrow(res_summ)+4),
           cols=which(columns_mtd$wrap_text==1),
           gridExpand = TRUE)
  
  if(is.null(file)){
    return(wb)
    
  }else{
    #save
    saveWorkbook(wb, file, overwrite = TRUE)
    
  }
  
}
