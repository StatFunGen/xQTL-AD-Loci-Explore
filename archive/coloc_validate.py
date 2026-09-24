import gzip,csv,collections
R221='releases/loci221_20260911T191256Z'; R182='releases/loci182_20260911T163540Z'
VAR='AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'
COL='res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged_overlapADloci.csv.gz'
MET='res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged_overlapADloci.csv.gz'
def read_loci(p):
    iv={}
    with gzip.open(p,'rt') as f:
        rd=csv.DictReader(f); cs=rd.fieldnames
        cc=[c for c in cs if c.lower() in ('#chr','chr','chrom')][0]
        pc=[c for c in cs if c.lower() in ('pos','start','position')][0]
        for r in rd:
            k=r['ADlocusID']
            try: q=int(float(r[pc]))
            except: continue
            ch=str(r[cc]).replace('chr','')
            if k in iv:
                a=iv[k]; a[1]=min(a[1],q); a[2]=max(a[2],q)
            else: iv[k]=[ch,q,q]
    return iv
old=read_loci(R182+'/'+VAR); new=read_loci(R221+'/'+VAR)
print('old loci',len(old),' new loci',len(new))
byc=collections.defaultdict(list)
for k,(c,lo,hi) in old.items(): byc[c].append((k,lo,hi))
n2o={k:[ok for ok,olo,ohi in byc[c] if olo<=hi and lo<=ohi] for k,(c,lo,hi) in new.items()}
o2n=collections.defaultdict(list)
for nk,oks in n2o.items():
    for ok in oks: o2n[ok].append(nk)
cls={}
for nk,oks in n2o.items():
    if len(oks)==0: cls[nk]='novel'
    elif len(oks)>1: cls[nk]='refined'
    else: cls[nk]='refined' if len(o2n[oks[0]])>1 else 'unchanged'
print(collections.Counter(cls.values()))
def load_coloc(p):
    mv=collections.defaultdict(float); gn=collections.defaultdict(set); cs=collections.defaultdict(set); nrow=0
    with gzip.open(p,'rt') as f:
        for r in csv.DictReader(f):
            k=r.get('ADlocusID','')
            if not k or k=='NA': continue
            nrow+=1
            try: v=float(r.get('vcp') or 0)
            except: v=0.0
            if v>mv[k]: mv[k]=v
            if r.get('gene_ID'): gn[k].add(r['gene_ID'])
            if r.get('cos_ID'): cs[k].add(r['cos_ID'])
    print('  rows with ADlocusID:',nrow,' loci touched:',len(mv))
    return mv,gn,cs
print('\n== xQTL coloc =='); mvx,gnx,csx=load_coloc(R221+'/'+COL)
print('== meta coloc ==');  mvm,gnm,csm=load_coloc(R221+'/'+MET)
grp={'INCREASED (refined+novel)':[k for k,c in cls.items() if c!='unchanged'],'  - refined':[k for k,c in cls.items() if c=='refined'],'  - novel':[k for k,c in cls.items() if c=='novel'],'UNCHANGED':[k for k,c in cls.items() if c=='unchanged']}
print('\n%-28s %5s %9s %9s %9s %9s %9s'%('group','n','anyCoS','vcp>=.5','vcp>=.8','vcp>=.95','medGenes'))
for nm,ks in grp.items():
    n=len(ks) or 1
    a=sum(1 for k in ks if k in mvx)
    f5=sum(1 for k in ks if mvx.get(k,0)>=0.5)
    f8=sum(1 for k in ks if mvx.get(k,0)>=0.8)
    f95=sum(1 for k in ks if mvx.get(k,0)>=0.95)
    g=sorted(len(gnx.get(k,())) for k in ks); med=g[len(g)//2] if g else 0
    print('%-28s %5d %8.1f%% %8.1f%% %8.1f%% %8.1f%% %9d'%(nm,len(ks),100*a/n,100*f5/n,100*f8/n,100*f95/n,med))
print('\n-- same, using GWAS-meta coloc file --')
for nm,ks in grp.items():
    n=len(ks) or 1
    a=sum(1 for k in ks if k in mvm)
    f8=sum(1 for k in ks if mvm.get(k,0)>=0.8)
    u=sum(1 for k in ks if (k in mvx or k in mvm))
    print('%-28s n=%-4d anyCoS %5.1f%%  vcp>=.8 %5.1f%%  union(xQTL|meta) %5.1f%%'%(nm,len(ks),100*a/n,100*f8/n,100*u/n))
import statistics
print('\n=== locus WIDTH by class (bp) ===')
w={'unchanged':[],'refined':[],'novel':[]}
for k,c in cls.items():
    ch,lo,hi=new[k]; w[c].append(hi-lo)
for c in ('unchanged','refined','novel'):
    v=sorted(w[c]); n=len(v)
    print('%-10s n=%-4d median %9d   mean %9d   IQR %9d - %9d'%(c,n,statistics.median(v),statistics.mean(v),v[n//4],v[3*n//4]))
inc=w['refined']+w['novel']
print('%-10s n=%-4d median %9d   mean %9d'%('INCREASED',len(inc),statistics.median(inc),statistics.mean(inc)))
