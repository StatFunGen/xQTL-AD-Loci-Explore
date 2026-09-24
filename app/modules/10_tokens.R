# ---- palette ---------------------------------------------------------------
TWAS_ONLY <- "T6"                # a real tier: TWAS/MR evidence only
NO_TIER   <- "No tier assigned"  # gene is in the locus, pipeline gave it no tier
NO_GENE <- "No gene mapped"
NO_P    <- "no p-value"   # variant rows the GWAS table carries no p-value for

NAVY <- "#ffffff"   # navbar/page chrome
TEAL <- "#2a78d6"   # the single accent
INK  <- "#0b0b0b"   # primary ink
MUTE <- "#52514e"   # secondary ink
LINE <- "#e1e0d9"   # hairline separator
PAPER<- "#fcfcfb"   # page / chart surface
na_fill <- function(x, fill) { x[is.na(x)] <- fill; x }   # base replacement for tidyr::replace_na
SIG_COL <- c("genome wide" = "#0ca30c", "suggestive" = "#fab219", "ns" = "#898781")
SIG_COL[NO_P] <- "#2a78d6"

# ---- cell-type config ------------------------------------------------------
ct_cols   <- c("ct_Brain_xQTL","ct_Exc_xQTL","ct_Inh_xQTL","ct_Oli_xQTL",
               "ct_OPC_xQTL","ct_Ast_xQTL","ct_Microglia_xQTL","ct_Bulk_Immune_xQTL")
ct_labels <- c("Brain","Excitatory","Inhibitory","Oligodendrocyte",
               "OPC","Astrocyte","Microglia","Bulk Immune")
ct_colors <- c("#2a78d6","#eb6834","#1baf7a","#eda100",
               "#e87ba4","#008300","#4a3aa7","#e34948")
names(ct_colors) <- ct_labels

trans_types <- c(
  "cis-eGene"   = "trans_n_genes",      # (kept generic label; counts below)
  "snRNA"       = "trans_snRNA_n_genes",
  "pQTL"        = "trans_pQTL_n_genes",
  "gpQTL"       = "trans_gpQTL_n_genes",
  "Hotspot"     = "trans_hotspot_n_genes")

# ---- citation identifiers for this release ---------------------------------
# Fill synapse_id or zenodo_doi once the release is deposited. The citation
# panel, the BibTeX export and the RIS export read these directly, so nothing
# else needs editing. Set publisher to the repository that issues the record.
REL_CITE <- list(
  publisher  = "Synapse, AD Knowledge Portal",
  synapse_id = "syn68872650",
  zenodo_doi = "",
  grant      = "U01AG072572"
)
