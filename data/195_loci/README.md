# Released tables - 188-loci build

Locus-level results for the 195 candidate Alzheimer's disease loci in this
release, produced by build_AD_locus_table.R. See the repository README for how
to regenerate them.

Every table here comes from the same build. Verified locus and row counts:

| file | rows | loci | contents |
|---|---|---|---|
| unified_AD_loci_xQTL_summary_20261002.xlsx | - | 195 | unified summary workbook: locus, gene and variant evidence across fine-mapping, colocalization, TWAS/MR and cTWAS, with confidence tiers T1-T6 |
| AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz | 7,618 | 195 | one row per variant-locus: variants entering a locus through a 95% credible set or a colocalization and reaching p < 1e-5, with per-variant GWAS summary columns |

| ADlocus_name_map_188_to_195.csv | 195 | - | maps each locus to its name in the 188-locus build: same, boundary widened, or new |

## Not included here

xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz is the
long-form table with one row per variant-locus-method-context-gene. At 133 MB it exceeds
what a git repository can hold, as does the noTrans variant at 111 MB.
Both are available on request: open an issue on this repository or email
jaempawi@bu.edu.

The unfiltered AD_loci_unified_cs95orColocs.csv.gz spans 450 candidate regions
before the p-value filter, so it does not describe this release and is not
published here.

## Provenance

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.

GWAS credible sets in this build were fine-mapped with SuSiE-RSS EB-mix. Colocalization evidence includes the tpQTL and CUIMC1 eQTL ColocBoost export, trans-xQTL-only colocalization and transmap AD-target colocalization. Trans fine-mapping is limited to the selected single-context set (FunGen_xQTL.trans.exported.toploci); genome-wide trans-eGene, trans-hotspot, trans pQTL and trans gpQTL fine-mapping are not included. See `../README.md` for what differs between the builds.