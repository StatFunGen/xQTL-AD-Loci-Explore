#!/usr/bin/env python3
"""Pre-deploy release check for the AD Loci Explorer.

Recomputes every headline number from the data files and compares it with
the recorded manifest. Run with --init to record the current values as the
new expected truth; run with no arguments to check. Exit code 1 on any FAIL.
"""
import csv, json, os, re, sys, gzip, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ROOT, 'scripts', 'release_manifest.json')
EXPECTED_MODALITIES = ['caQTL','eQTL','gpQTL','haQTL','mQTL','pQTL','sQTL']
LEGACY_NAMES = ['Kellis','DeJager']

def rows(path):
    with open(os.path.join(ROOT, path), newline='') as fh:
        return list(csv.DictReader(fh))

def measure():
    d = rows('data.csv')
    m = {}
    m['data_rows'] = len(d)
    m['loci'] = len({r['ADlocus'] for r in d if r['ADlocus']})
    m['genes'] = len({r['gene'] for r in d if r['gene']})
    tiers = collections.Counter(r.get('top_confidence','') for r in d)
    m['tier_counts'] = {k: v for k, v in sorted(tiers.items()) if k}
    m['genes_t1_t5'] = len({r['gene'] for r in d
                            if r.get('top_confidence') in ('T1','T2','T3','T4','T5') and r['gene']})
    m['genes_t1_t6'] = len({r['gene'] for r in d
                            if re.fullmatch(r'T[1-6]', r.get('top_confidence') or '') and r['gene']})
    sig = collections.Counter(r.get('significance','') for r in d)
    m['significance'] = {k: v for k, v in sorted(sig.items())}
    blob = open(os.path.join(ROOT, 'data.csv')).read()
    mods = collections.Counter(re.findall(r'[A-Za-z]+QTL', blob))
    for junk in ('ADxQTL', 'xQTL'):
        mods.pop(junk, None)
    m['modalities'] = sorted(mods)
    m['legacy_name_hits'] = sum(blob.count(n) for n in LEGACY_NAMES)
    m['chroms'] = sorted({r['chr'] for r in d if r['chr']},
                         key=lambda c: (len(c), c))
    rx = re.compile(r'^(\S+)\s+([A-Za-z-]+QTL)([+.-]*)\s*\(T([0-9]),n=([0-9]+)\)$')
    tot = ok = na = 0
    for r in d:
        for tok in (r.get('ordered_contexts') or '').split(';'):
            t = tok.strip()
            if not t: continue
            tot += 1
            if rx.match(t): ok += 1
            elif t.startswith('NA.'): na += 1
    m['ordered_tokens'] = tot
    m['ordered_tokens_parsed'] = ok
    m['ordered_tokens_na'] = na
    ctxs = sorted({(r.get('context') or '').strip() for r in d if (r.get('context') or '').strip()})
    m['contexts'] = len(ctxs)
    src = open(os.path.join(ROOT, 'modules', '45_views.R')).read()
    keys = set(re.findall(r'^  "([^"]+)" = c\(', src, re.M))
    m['ctx_map_entries'] = len(keys)
    m['contexts_not_in_map'] = sorted(set(ctxs) - keys)
    m['locus_summary_rows'] = len(rows('downloads/AD_locus_summary_release.csv'))
    m['gene_tier_rows'] = len(rows('downloads/gene_tier_assignment_release.csv'))
    with gzip.open(os.path.join(ROOT,'downloads/AD_locus_variants_release.csv.gz'),'rt',newline='') as fh:
        m['variant_rows'] = sum(1 for _ in csv.DictReader(fh))
    return m, d

def ui_literals():
    src = open(os.path.join(ROOT, 'modules', '80_ui.R')).read()
    out = {}
    hit = re.findall(r'([0-9]{1,3}(?:,[0-9]{3})+) cis xQTL credible sets', src)
    out['credible_sets_prose'] = hit[0] if hit else None
    hit = re.findall(r'dc-count-n dc-mono", "([0-9,]+)"', src)
    out['credible_sets_kpi'] = hit[0] if hit else None
    hit = re.findall(r'span (\w+) molecular modalities', src)
    out['modality_word'] = hit[0] if hit else None
    return out

def main():
    m, d = measure()
    lit = ui_literals()
    init = '--init' in sys.argv
    if init:
        json.dump({'measured': m, 'ui': lit}, open(MANIFEST,'w'), indent=2, sort_keys=True)
        print('manifest written to', MANIFEST)
        return 0
    if not os.path.exists(MANIFEST):
        print('FAIL  no manifest; run with --init first')
        return 1
    exp = json.load(open(MANIFEST))
    fails = []
    def check(name, got, want):
        ok = got == want
        print(('PASS  ' if ok else 'FAIL  ') + name + ': ' + repr(got) +
              ('' if ok else '  expected ' + repr(want)))
        if not ok: fails.append(name)
    for k in sorted(exp['measured']):
        check(k, m.get(k), exp['measured'][k])
    # invariants that must hold regardless of the manifest
    tot = sum(m['significance'].values())
    check('significance classes sum to data_rows', tot, m['data_rows'])
    check('every ordered_contexts token is parsed or a null marker',
          m['ordered_tokens_parsed'] + m['ordered_tokens_na'], m['ordered_tokens'])
    check('every context value is in CTX_MAP', m['contexts_not_in_map'], [])
    check('modality set matches release', m['modalities'], EXPECTED_MODALITIES)
    check('no legacy dataset names', m['legacy_name_hits'], 0)
    check('gene tier file reconciles with T1-T6 gene set',
          m['gene_tier_rows'], m['genes_t1_t6'])
    bad = [c for c in m['chroms'] if c not in [str(i) for i in range(1,23)]]
    check('all chromosomes in 1-22', bad, [])
    # UI literals that are hardcoded rather than derived
    for k in sorted(exp['ui']):
        check('ui literal ' + k, lit.get(k), exp['ui'][k])
    if lit['credible_sets_kpi'] != lit['credible_sets_prose']:
        print('FAIL  credible-set figure disagrees between KPI and prose')
        fails.append('credible_sets_consistency')
    print('')
    print('UNVERIFIED  cis-xQTL credible sets (' + str(lit['credible_sets_kpi']) +
          ') is a hardcoded literal with no backing file in this app.')
    print('')
    print(('FAILED: ' + str(len(fails))) if fails else 'All checks passed.')
    return 1 if fails else 0

if __name__ == '__main__':
    sys.exit(main())
