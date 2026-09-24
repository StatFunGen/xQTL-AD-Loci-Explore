#!/usr/bin/env python3
"""Assign AD locus-to-gene evidence tiers T1-T5 from a pipeline release.

Tier definitions follow the FunGen-xQTL manuscript Methods:
  T1  95% single-context CS overlap + (MR or causal-TWAS)
  T2  95% AD-xQTL CoS               + (MR or causal-TWAS)
  T3  (95% CS overlap or CoS)       + TWAS
  T4  95% CS overlap alone
  T5  CoS alone, or CS overlap at 50%/70% coverage
A gene is reported at its highest supported tier.
"""
import gzip, csv, collections, sys

RELEASE = sys.argv[1] if len(sys.argv) > 1 else 'releases/loci221_20260911T191256Z'

VARIANTS = 'res_AD_variants_xQTL.csv.gz'
# T5 counts localized evidence from ANY AD-xQTL colocalization source, so all three
# ColocBoost outputs are read. Using only the xQTL file silently drops genes whose
# sole support is epigenomic or meta-GWAS colocalization.
COLOC_FILES = [
    'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged_overlapADloci.csv.gz',
    'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged_overlapADloci.csv.gz',
    'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged_overlapADloci.csv.gz',
]
CTWAS    = 'res_AD_cTWAS_pip075_overlapADloci.csv.gz'
MR       = 'res_AD_XWAS_MR_filtered_TWAS_sig_overlapADloci.csv.gz'

def rows(name):
    with gzip.open(f'{RELEASE}/{name}', 'rt') as f:
        yield from csv.DictReader(f)

def in_credible_set(value):
    """cs_coverage_* holds a credible-set INDEX. '' means none; '0' is a
    placeholder, not membership. Only a non-zero integer is real CS membership."""
    try:
        return int(str(value).strip()) != 0
    except (TypeError, ValueError):
        return False

# --- localized evidence -------------------------------------------------
# Methods: localized evidence must resolve BOTH the posterior variant set and the
# molecular context. Single-context fine-mapping supplies context directly; multi-context
# components reach the tiers only via CoS. So CS-overlap evidence is restricted here.
SINGLE_CONTEXT = 'single_context_finemapping'

cs95, cs_low = set(), set()          # (locus, gene) with 95% / 50-70% CS overlap
for r in rows(VARIANTS):
    key = (r.get('ADlocusID'), r.get('gene_ID'))
    if not all(key):
        continue
    if r.get('Method') != SINGLE_CONTEXT:
        continue
    if in_credible_set(r.get('cs_coverage_0.95')):
        cs95.add(key)
    elif in_credible_set(r.get('cs_coverage_0.7')) or in_credible_set(r.get('cs_coverage_0.5')):
        cs_low.add(key)

cos = set()
for _f in COLOC_FILES:
    try:
        n0 = len(cos)
        for r in rows(_f):
            if r.get('ADlocusID') and r.get('gene_ID'):
                cos.add((r['ADlocusID'], r['gene_ID']))
        print(f'  coloc {_f.split("_unified")[0][10:]:<12} +{len(cos)-n0} new pairs')
    except FileNotFoundError:
        print(f'  coloc {_f[:40]:<40} MISSING')

# --- gene-level support -------------------------------------------------
ctwas_genes = {r['gene_ID'] for r in rows(CTWAS) if r.get('cTWAS_signif') == 'TRUE'}
mr_genes, twas_genes = set(), set()
for r in rows(MR):
    if r.get('MR_signif') == 'TRUE':
        mr_genes.add(r['gene_ID'])
    if r.get('TWAS_signif') == 'TRUE':
        twas_genes.add(r['gene_ID'])

# --- tier assignment ----------------------------------------------------
def tier(locus, gene):
    strong = gene in mr_genes or gene in ctwas_genes
    has_cs, has_cos = (locus, gene) in cs95, (locus, gene) in cos
    if has_cs and strong:                      return 'T1'
    if has_cos and strong:                     return 'T2'
    if (has_cs or has_cos) and gene in twas_genes: return 'T3'
    if has_cs:                                 return 'T4'
    if has_cos or (locus, gene) in cs_low:     return 'T5'
    return None

# Published tiering restricts genes to those with a gene_name on the GRCh38.103
# reference (complete_ADlocus_level_summary_202605.R merges this map, then drops
# unmatched rows). Without it, unmapped Ensembl IDs inflate the stand-alone-CS tier.
GENE_REF = 'alexandre/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list'
import os, re
ref_genes = set()
if os.path.exists(GENE_REF):
    with open(GENE_REF) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 5 and parts[3].startswith('ENSG') and parts[4].strip():
                ref_genes.add(parts[3].split('.')[0])
    print(f'  gene reference : {len(ref_genes)} genes with a gene_name')

def base(g):
    return str(g).split('.')[0]

ORDER = ['T1', 'T2', 'T3', 'T4', 'T5']
best = {}
dropped = set()
for locus, gene in cs95 | cos | cs_low:
    if ref_genes and base(gene) not in ref_genes:
        dropped.add(gene)
        continue
    t = tier(locus, gene)
    if t and (gene not in best or ORDER.index(t) < ORDER.index(best[gene])):
        best[gene] = t

counts = collections.Counter(best.values())
print(f'  dropped (no gene_name on reference): {len(dropped)} genes')
print(f'release: {RELEASE}\n')
print(f'  evidence pairs : 95%CS {len(cs95)}   CoS {len(cos)}   low-cov CS {len(cs_low)}')
print(f'  gene support   : MR {len(mr_genes)}   cTWAS {len(ctwas_genes)}   TWAS {len(twas_genes)}\n')
for t in ORDER:
    print(f'  {t}  {counts.get(t, 0):5d} genes')
t14 = sum(counts.get(t, 0) for t in ORDER[:4])
print(f'\n  T1-T2 {sum(counts.get(t,0) for t in ORDER[:2]):5d}')
print(f'  T3-T4 {sum(counts.get(t,0) for t in ORDER[2:4]):5d}')
print(f'  T1-T4 {t14:5d}   <- stringent set')
print(f'  T1-T5 {t14 + counts.get("T5", 0):5d}   <- all tiers')
