# xQTL AD Loci Explorer

Code and released tables behind the **AD Loci Explorer**, a browsable view of the
ADSP FunGen-xQTL release across 188 candidate Alzheimer's disease loci.

Live app: https://jenny-empawi.shinyapps.io/xQTL-AD-loci-Explore/

## Layout

```
scripts/
  build_AD_locus_table.R   locus-level build: evidence integration and tier assignment
  gene_prio_utils.R        helper functions, including the T1-T6 tier rules
  add_evidence_source.R    register a new evidence source (see 'Adding data')
  validate_outputs.R       check a release against the expectations of this build
  run_pipeline.qsub        end-to-end job script for the SCC scheduler
config/                    metadata tables the build reads
data/                      released tables, one directory per build
app/                       Shiny application, runnable from a clone
archive/                   superseded and ad-hoc scripts; not part of the pipeline
DEPENDENCIES.md            R version and package versions
```

## Requirements

R 4.5.3 with the packages listed in DEPENDENCIES.md:

```r
install.packages(c("data.table","stringr","openxlsx","tidyverse","ggplot2",
                   "scales","circlize","shiny","bslib","DT","plotly",
                   "shinycssloaders"))
remotes::install_github("StatFunGen/pecotmr")   # not on CRAN
```

## Running it

There are two ways in, depending on whether you want to reproduce the published
release or add your own data to it.

### Option 1: reproduce the release from scratch

Rebuilds the locus table from the registered inputs. This needs access to the
tables listed in `config/metadata_analysis.csv`, which are controlled-access; see
Inputs for how to obtain them.

```bash
export AD_LOCI_ROOT=/path/to/AD_loci_xQTL      # where the input trees live
export AD_LOCI_STAGING=/path/to/staging        # three precomputed tables, see Inputs
export AD_LOCI_OUT=$AD_LOCI_ROOT/out_$(date +%Y%m%d)    # optional
export AD_LOCI_CONFIG=/path/to/private/config          # optional, see Inputs

Rscript scripts/build_AD_locus_table.R                 # -> $AD_LOCI_OUT
Rscript scripts/validate_outputs.R  $AD_LOCI_OUT       # expected tables, 188 loci, 508 genes
Rscript app/build_shiny_data.R $AD_LOCI_OUT app/data.csv
```

On the SCC the same three steps run as one job. From the repository root:

```bash
qsub scripts/run_pipeline.qsub
```

The job records the commit, host, time and the config it used in
`_provenance.csv` beside the outputs, then validates the release it produced.

The build checks every required input before starting and stops with a list of
whatever is missing, rather than failing part-way through a long run.

### Option 2: add your own evidence to the published build

Register new results as an evidence source and rebuild with them included.
`scripts/add_evidence_source.R` converts a ColocBoost-style colocalization
summary and adds its rows to the registry; see Adding data for the input format,
and for what a source in any other format needs.

```bash
Rscript scripts/add_evidence_source.R --name mysource --bed coloc_summary.bed
Rscript scripts/build_AD_locus_table.R
```

Adding a source is a configuration change, not a code change. Be aware that the
rebuild still reads the whole registry, so this route needs the same input
access as option 1; it saves you from modifying the pipeline, not from having
the data.

If all you want is the Explorer, neither option is needed: `app/data.csv` ships
with the repository. See Running the app.

## Inputs

`config/metadata_analysis.csv` is the registry of every input the build reads.
One row per exported analysis table:

| column | meaning |
|---|---|
| context | dataset label, joined to config/contexts_metadata.csv |
| Data Type | eQTL, sQTL, pQTL, ... |
| Cohort | ROSMAP, MSBB, KNIGHT, ... |
| Modality | assay or analysis modality |
| Method | how the table was produced; selects the loader in the build |
| Path | location of the table, relative to `$AD_LOCI_ROOT` |
| summary_file, summary_file_ad | written by the build at run time; leave blank |
| variant_level_method | whether the method contributes variant-level evidence |

Paths under `$AD_LOCI_ROOT` refer to the FunGen-xQTL release, Synapse
`syn68872650`. The layout on Synapse does not necessarily mirror these relative
paths; email jaempawi@bu.edu if you need help locating a particular table. Some
rows carry `<user>` and `<collaborator>` placeholders where
tables were exported from per-user analysis directories; replace these with the
directory names in your own copy.

Three precomputed tables are too large to distribute with the code and are read
from `$AD_LOCI_STAGING`:

| file | size | source |
|---|---|---|
| gwas_variants_cor0.5.csv.gz | 295 MB | Synapse `syn75082260` |
| res_APOE_interaction_summ.csv.gz | 45 MB | on request |
| res_msex_interaction_summ.csv.gz | 1.6 MB | on request |

The two interaction summaries are not redistributed publicly. To request them,
open an issue on this repository or email jaempawi@bu.edu.

The correlation table drives the LD-based credible-set extension, so the
released variant sets cannot be reproduced without it.

## Adding data

To add a new evidence source you register it in `config/metadata_analysis.csv`.
The build auto-detects any Method there that points at per-context toploci files
and has no dedicated loader, so registration is the only step.

For a ColocBoost-style colocalization summary, `scripts/add_evidence_source.R` does the
conversion and registration for you:

```bash
Rscript scripts/add_evidence_source.R --name transmap \
        --bed trans_xQTL_only_colocalization_summary_table.bed
Rscript scripts/build_AD_locus_table.R
```

It is safe to re-run: conversion is skipped when outputs already exist, and rows
for the same `--name` are replaced rather than duplicated.

### Input format

The `--bed` file is tab-delimited with a header, one row per variant, and event
lists delimited by `;`:

| # | column | notes |
|---|---|---|
| 1 | #chr | |
| 2 | start | |
| 3 | end | |
| 4 | a1 | |
| 5 | a2 | |
| 6 | variant_ID | chr:pos:ref:alt |
| 7 | region_ID | |
| 8 | event_ID | semicolon-delimited list; context is the part before `_cis_` or `_trans_` |
| 9 | cos_ID | must contain `cosN`; rows with `cos0` are skipped |
| 10 | vcp | becomes PIP in the emitted file |
| 11 | cos_npc | becomes purity |
| 14 | zscore | semicolon-delimited, positionally matched to event_ID |

It emits one toploci file per context with the columns the build expects:
`#chr start end a1 a2 variant_ID region_ID event_ID cs_coverage_0.95 purity PIP z`.

A source in any other format needs its own loader in `scripts/build_AD_locus_table.R`;
follow one of the existing `Method ==` blocks.

## Methodology

This section describes what the build does, step by step, and the thresholds it
applies. Two terms recur. A *credible set* is the small group of variants that
fine-mapping says most likely contains the causal one; 95% coverage (cs95) means
the set is built to have a 95% chance of containing it, and cs70 and cs50 are
weaker versions of the same idea. *Colocalization* asks whether an AD signal and
a molecular signal in the same place are driven by the same variant rather than
by two neighbouring ones.

### Step 1: the AD locus set

The locus set is imported, not derived here: GWAS fine-mapping is run upstream
and is not repeated by this pipeline. The credible sets behind the current
build come from eight AD GWAS studies fine-mapped with SuSiE-RSS EB-mix; the
earlier 183-loci release used SuSiE-RSS version 1 on the same studies, which is
why the two are not nested (see `data/README.md`).
At startup the build scans the registry
for the analysis methods that point at per-context top-loci files and reports
what it found, so the source list is discovered from the configuration rather
than hardcoded. The job log records it as `[toploci] auto-detected sources`.

### Step 2: harmonization

Nothing is merged until the identifiers agree.

*Contexts.* Dataset and assay labels are mapped to one harmonized context name
through `config/contexts_metadata.csv`. Single-nucleus tables that carry a cell
type instead of a context are mapped by cell type.

*Variants.* Identifiers arrive in several spellings. Underscores are converted
to colons and a `chr` prefix is added where it is missing, but the merge itself
is done on chromosome and position rather than on the identifier string, and the
AD locus table's spelling is then adopted for the merged row. This is what
absorbs allele-order differences: where two sources name the same position with
the alleles written the other way round, the merge still succeeds, and the build
counts the disagreeing spellings instead of dropping those variants (213 of
5,318 in this build).

*Events and genes.* Event identifiers have trailing `_chr`, `_ENSG` and `_gp`
fragments stripped. Genes are matched on the Ensembl ID with the version suffix
removed, so `ENSG00000141510.16` and `ENSG00000141510` are treated as the same
gene.

*GWAS studies.* Source and study names are harmonized so that one GWAS is not
counted twice under two labels.

### Step 3: evidence families and their significance rules

| evidence | what makes it count |
|---|---|
| fine-mapping | membership in a credible set; the coverage level is recorded on each row as cs95, cs70 or cs50, and the strongest level available is kept per locus |
| colocalization | a variant-level colocalization probability (VCP), or the SNP-level PPH4 where that is what the method reports |
| TWAS | gene-level p < 2.5e-6, either in the selected best method or in more than half of the methods run, evaluated per gene, block, context and GWAS source |
| MR | TWAS-significant, and cPIP > 0.5, and at least two credible sets, and I2 < 0.5 |
| cTWAS | a non-empty credible set and SuSiE PIP > 0.75 |

A row that names a method but lacks that method's statistic is dropped, so a
colocalization row with no VCP and a fine-mapping row with no PIP never enter
the table.

### Step 4: one number per variant

Each surviving row is reduced to a single `variant_inclusion_probability` by
precedence rather than by averaging: the colocalization VCP when present,
otherwise the fine-mapping PIP, otherwise the SNP-level PPH4. The maximum across
every method and source is then taken per variant, and the method that produced
that maximum is recorded beside it. Because it is a maximum, adding a source can
only raise a variant's score and never lower it; a score that falls between two
builds therefore means the inputs or the filtering changed, not merely that
something was added. The GWAS methods, sources and effect directions behind each
variant are collapsed into parallel `|`-separated fields on the same row.

### Step 5: assembling and filtering loci

A locus identifier is built from its chromosome and the first and last positions
of its member variants, and the loci are then numbered in order. That number is
an ordinal within a single build, so the same number does not denote the same
locus in another build and cross-build comparisons have to be made on variant
IDs.

Three filters decide what survives. A locus is kept only if it carries a GWAS
credible set at 95% coverage or GWAS colocalization, so a locus supported only
by a cs70 or cs50 set is dropped. A locus is kept only if its best GWAS p-value
is below 1e-5. Each surviving locus is then represented by up to five variants
whose GWAS PIP or VCP exceeds 0.1; where no variant clears 0.1, the single best
one is shown, ranked on the GWAS probability first and the xQTL probability
second. GWAS significance is reported in bands: p < 1e-5, p < 1e-6, p < 5e-8.

Two flags are carried per locus: whether it is supported only by proxy-based
GWAS, and whether it falls in the APOE region (chr19:43,905,790-45,905,791),
which is reported separately because the extreme linkage disequilibrium there
makes fine-mapping and colocalization hard to read.

### Step 6: tiers, outputs and provenance

Tiers are assigned per variant, gene and context group by the rules in the next
section, and each gene is reported at the strongest tier it reaches anywhere.
The build writes the locus-level summary, the variant-level membership table and
the unified workbook to `$AD_LOCI_OUT`, together with a `_provenance.csv`
recording the commit, host, time and the exact configuration the run used, and a
copy of the two config tables as they stood. `scripts/validate_outputs.R` then
checks that release: the expected tables are present, 188 loci, 508 genes, and
the published tier distribution.

The output directory is deliberately created empty. The GWAS credible-set
harmonization recomputes from its inputs rather than reusing a previous
release's artifacts, so a rerun cannot silently inherit stale intermediate
files.

## Confidence tiers

| tier | evidence |
|---|---|
| T1 | cTWAS or MR, plus a 95% credible set overlap from single-context or fSuSiE fine-mapping |
| T2 | cTWAS or MR, plus colocalization |
| T3 | TWAS, plus either a 95% credible set overlap or colocalization |
| T4 | a 95% credible set overlap from single-context or fSuSiE fine-mapping |
| T5 | colocalization, or a fine-mapping overlap that is not a single-context 95% set (multi-context, cs50, cs70) |
| T6 | gene-level TWAS, MR, or cTWAS only, with no localized xQTL signal |

`scripts/gene_prio_utils.R` assigns `top_confidence` per row during the build. T1-T5
describe genes with localized AD-xQTL support, ordered by strength of evidence.
**T6** covers genes supported only by gene-level TWAS/MR or cTWAS evidence, with
no localized xQTL signal; because the tier chain is evaluated over rows of the
xQTL overlap table, such a gene never reaches the final branch, so the build
assigns T6 directly from the gene-level tables produced earlier in the same run.

A gene is reported at its strongest tier. The app takes its tiers from the
release, so the table and the Explorer cannot diverge.

## Running the app

```r
shiny::runApp("app")
```

`app/data.csv` ships with the repository, so the Explorer runs without the
pipeline. Note that `app/build_shiny_data.R` carries a number of gene-level
columns forward from the existing `app/data.csv`, so that file is an input as
well as an output.

## Data

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
