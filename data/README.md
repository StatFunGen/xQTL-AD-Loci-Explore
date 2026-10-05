# Released tables

Three builds are published here. They are kept in separate directories because the
same workbook filename has been used across builds, and the `202605` / `202609`
suffix is the build version of the summary script, not the number of loci.

| directory | build | workbook | rows | loci |
|---|---|---|---|---|
| 183_loci/ | earlier release | unified_AD_loci_xQTL_summary_202605.xlsx | 4,114 | 183 |
| 188_loci/ | earlier release | unified_AD_loci_xQTL_summary_202609.xlsx | 4,195 | 188 |
| 195_loci/ | current release, shown in the Explorer | unified_AD_loci_xQTL_summary_20261002.xlsx | 5,016 | 195 |

Counts are measured from the `Gene Locus table` sheet of each workbook and agree
on both `locus index` and `ADlocus`.

## What differs between 183_loci and 188_loci

Both builds draw on the same eight AD GWAS studies. They differ in how those
studies were fine-mapped, and in one added evidence source.

| | 183_loci (202605) | 188_loci (202609) |
|---|---|---|
| GWAS studies | the same eight | the same eight, in the current format |
| GWAS fine-mapping | SuSiE-RSS version 1 | SuSiE-RSS EB-mix |
| colocalization | as in that release | adds transmap colocalization |

The credible sets were recomputed with a different method rather than simply
extended, so the two locus sets are not nested and neither contains the other.
Comparing the two workbooks on variant IDs: 143 loci carry the same identifier,
8 loci in the 188-loci build share no variant with any locus in the 183-loci
build, and 8 loci in the 183-loci build share no variant with any locus in the
188-loci build; the remainder shift their boundaries. Of the 1,182 variants
carrying a maximum inclusion score in both, 869 are identical, 120 are higher in
the 188-loci build and 193 are lower.

That some scores are *lower* is expected here and is worth understanding. Within
a single build the inclusion score is a maximum across methods and sources, so
adding a source can only raise it. It can fall only when the underlying
probabilities are themselves recomputed, which is exactly what changing the
fine-mapping method does. A drop between these two builds is therefore a
property of the method change, not evidence that something was removed.

Locus numbers are ordinals within one build, so the same number does not denote
the same locus across the two. Always compare on variant IDs.

## What differs between 188_loci and 195_loci

Both builds use the same eight AD GWAS studies, fine-mapped the same way. The 195-loci build changes the xQTL evidence:

| | 188_loci (202609) | 195_loci (20261002) |
|---|---|---|
| colocalization | ColocBoost exports and transmap colocalization | adds tpQTL and CUIMC1 eQTL colocalization with AD, and transmap AD-target colocalization |
| trans evidence | all trans fine-mapping results | the selected single-context trans fine-mapping set plus trans-xQTL-only and AD-target colocalization; genome-wide trans-eGene, trans-hotspot, trans pQTL and trans gpQTL fine-mapping are not included |

Loci are defined by a 95% credible set or an AD colocalization, so the added AD colocalizations define 7 new loci, and 4 existing loci gain a wider boundary. All 188 earlier loci are kept. Locus names are ordinals and are renumbered; `195_loci/ADlocus_name_map_188_to_195.csv` maps each locus to its earlier name. Gene tiers use cis evidence only.

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
