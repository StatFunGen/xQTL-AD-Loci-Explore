#!/usr/bin/env python3
"""Assign AD locus-to-gene confidence levels, reproducing Alexandre's classifier.

Source of truth: alexandre/gene_prio_utils.R, SummarizeTable(), lines 170-196.
His labels are CL1-CL6; the manuscript renames CL1-CL5 to T1-T5 and keeps CL6
(TWAS evidence only, no localization) as a separate untiered reference set.

  CL1/T1  (MR or cTWAS) AND cs95 single-context/fSuSiE fine-mapping overlap
  CL2/T2  (MR or cTWAS) AND AD-xQTL colocalization
  CL3/T3  TWAS AND (cs95 fine-mapping overlap OR colocalization)
  CL4/T4  cs95 single-context/fSuSiE fine-mapping overlap alone
  CL5/T5  colocalization OR any fine-mapping overlap (multi-context, cs50, cs70)
  CL6/T6  TWAS only

Classification is per (variant_ID, gene_ID, context_short), matching his
by= clause; a gene's reported level is the strongest it reaches anywhere.
Rows whose context starts with 'AD' are excluded, as in his filter.
"""
import gzip, csv, collections, sys

RELEASE = sys.argv[1] if len(sys.argv) > 1 else 'releases/loci182_20260911T163540Z'
SRC     = 'res_AD_variants_xQTL.csv.gz'
FINEMAP = {'single_context_finemapping', 'fSuSiE_finemapping'}
CL5_MTD = {'AD_xQTL_colocalization', 'multi_context_finemapping',
           'single_context_finemapping', 'fSuSiE_finemapping', 'sn_sQTL', 'Coloc'}

def truthy(v):
    return str(v).strip().upper() == 'TRUE'

groups = collections.defaultdict(lambda: {
    'mr': False, 'ctwas': False, 'twas': False,
    'cs95_fm': False, 'cs95_fm_cov': False, 'cs95_fm_strict': False,
    'coloc': False, 'cl5': False})
genes_seen, genes_ad_only = set(), set()

with gzip.open(f'{RELEASE}/{SRC}', 'rt') as f:
    for r in csv.DictReader(f):
        gene = r.get('gene_ID')
        if not gene:
            continue
        genes_seen.add(gene)
        if str(r.get('context', '')).startswith('AD'):
            genes_ad_only.add(gene)
            continue                      # his !str_detect(context,'^AD') filter
        g = groups[(r.get('variant_ID'), gene, r.get('context_short'))]
        mtd, cset = r.get('Method'), str(r.get('credibleset') or '')
        cov = str(r.get('susie_coverage') or '')
        if truthy(r.get('MR_signif')):    g['mr'] = True
        if truthy(r.get('cTWAS_signif')): g['ctwas'] = True
        if truthy(r.get('TWAS_signif')):  g['twas'] = True
        if mtd in FINEMAP and 'cs95' in cset:
            g['cs95_fm'] = True
            if cov in ('cs95', 'cs70'):
                g['cs95_fm_cov'] = True          # as written by Alexandre
            if cov == 'cs95':
                g['cs95_fm_strict'] = True       # cs95 only, per the Methods text
        if mtd == 'AD_xQTL_colocalization': g['coloc'] = True
        if mtd in CL5_MTD:                  g['cl5'] = True

def level(g, strict=False):
    key = 'cs95_fm_strict' if strict else 'cs95_fm_cov'
    strong = g['mr'] or g['ctwas']
    if strong and g[key]:                              return 1
    if strong and g['coloc']:                          return 2
    if g['twas'] and (g['cs95_fm'] or g['coloc']):     return 3
    if g[key]:                                         return 4
    if g['cl5']:                                       return 5
    return 6

best, best_strict = {}, {}
for (_v, gene, _c), g in groups.items():
    for tgt, st in ((best, False), (best_strict, True)):
        lv = level(g, st)
        if gene not in tgt or lv < tgt[gene]:
            tgt[gene] = lv

unclassified = genes_seen - set(best)          # genes with only AD-context rows
counts = collections.Counter(best.values())

print(f'release: {RELEASE}')
print(f'  variant x gene x context groups : {len(groups)}')
print(f'  genes in source                 : {len(genes_seen)}')
print(f'  genes classified                : {len(best)}')
print(f'  genes UNCLASSIFIED (AD-context only): {len(unclassified)}')
print()
for i in range(1, 7):
    print(f'  T{i} (CL{i})  {counts.get(i, 0):5d} genes')
t14 = sum(counts.get(i, 0) for i in (1, 2, 3, 4))
print()
print(f'  T1-T2 {counts.get(1,0)+counts.get(2,0):5d}')
print(f'  T3-T4 {counts.get(3,0)+counts.get(4,0):5d}')
print(f'  T1-T4 {t14:5d}   <- stringent set (published: 145)')
print(f'  T1-T5 {t14+counts.get(5,0):5d}   <- tiered (published: 410)')
print(f'  T6    {counts.get(6,0):5d}   <- TWAS only, untiered reference set')
print(f'  TOTAL {len(best)+len(unclassified):5d}   <- every gene accounted for')

# ---- sensitivity: does accepting cs70 into CL1/CL4 change the stringent set? ----
cs = collections.Counter(best_strict.values())
t14s = sum(cs.get(i, 0) for i in (1, 2, 3, 4))
moved = {g: (best[g], best_strict[g]) for g in best if best[g] != best_strict[g]}
print()
print('=== cs95-only variant (CL1/CL4 require susie_coverage == cs95) ===')
for i in range(1, 7):
    print(f'  T{i}  {cs.get(i,0):5d}   (as-written: {counts.get(i,0)})')
print(f'\n  T1-T4 {t14s:5d}   (as-written: {t14})')
print(f'  T1-T5 {t14s+cs.get(5,0):5d}   (as-written: {t14+counts.get(5,0)})')
print(f'\n  genes that change level: {len(moved)}')
dem = collections.Counter((a, b) for a, b in moved.values())
for (a, b), n in sorted(dem.items()):
    print(f'    T{a} -> T{b} : {n}')
out = [g for g, (a, b) in moved.items() if a <= 4 < b]
print(f'\n  genes dropping OUT of the stringent T1-T4 set: {len(out)}')
print('   ', ', '.join(sorted(out)[:15]))

# ---- export per-gene assignments and the tier summary ----------------------
import os
tag = os.path.basename(RELEASE.rstrip('/'))
RULE = {1: 'CL1  (MR or cTWAS) + cs95 fine-mapping overlap',
        2: 'CL2  (MR or cTWAS) + colocalization',
        3: 'CL3  TWAS + (cs95 fine-mapping overlap or colocalization)',
        4: 'CL4  cs95 fine-mapping overlap alone',
        5: 'CL5  colocalization or any fine-mapping overlap (multi-context, cs50/cs70)',
        6: 'CL6  TWAS only, no localized evidence'}
with open(f'gene_tier_{tag}.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['gene_ID', 'tier', 'tier_cs95_only', 'rule'])
    for g in sorted(best, key=lambda x: (best[x], x)):
        w.writerow([g, f'T{best[g]}', f'T{best_strict[g]}', RULE[best[g]]])
with open(f'tier_summary_{tag}.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['tier', 'genes', 'genes_cs95_only', 'definition'])
    for i in range(1, 7):
        w.writerow([f'T{i}', counts.get(i, 0), cs.get(i, 0), RULE[i]])
    w.writerow(['T1-T4', t14, t14s, 'stringent set'])
    w.writerow(['T1-T5', t14 + counts.get(5, 0), t14s + cs.get(5, 0), 'all tiers'])
    w.writerow(['TOTAL', len(best), len(best_strict), 'every gene classified'])
print(f'\\nwrote gene_tier_{tag}.csv and tier_summary_{tag}.csv')
