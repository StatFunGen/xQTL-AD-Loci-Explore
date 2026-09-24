# ---- load & coerce ---------------------------------------------------------
dat <- read_csv("data.csv", show_col_types = FALSE)
dat <- dat %>% mutate(
  chr            = suppressWarnings(as.integer(chr)),
  # 214 rows carry no GWAS p-value. Without an explicit level they become NA
  # and the significance filter drops them even with every box ticked - the
  # same trap the tier column had. "no p-value" is not a claim of ns.
  significance   = factor(
    ifelse(is.na(significance) | !nzchar(trimws(as.character(significance))),
           NO_P, as.character(significance)),
    levels = c("genome wide","suggestive","ns", NO_P)),
  # Every gene the upstream classifier sees lands in T1-T6. Only variant rows
  # with no xQTL target gene at all lack a tier; give them an explicit level
  # so factor() does not turn them into NA and the tier filter keeps them.
  top_confidence = factor(
    # T6 is a real tier (CL6 in the upstream classifier: TWAS-only support).
    # Only rows with no xQTL target gene at all are left without a tier;
    # those get NO_GENE, which is not a tier label.
    ifelse(is.na(top_confidence) | !nzchar(trimws(as.character(top_confidence))),
           ifelse(is.na(gene) | !nzchar(trimws(as.character(gene))), NO_GENE, NO_TIER),
           as.character(top_confidence)),
    levels = c("T1","T2","T3","T4","T5","T6", NO_TIER, NO_GENE)),
  across(c(twas_sig, mr_sig, ctwas_sig, has_trans, all_of(ct_cols)), as.logical)
)

# precompute the cell-type dot column once (static per row)
# top_confidence is defined as the gene's MAXIMUM confidence level, but the source
# table populates it only on the row where the tier was determined, leaving the
# gene's other rows blank. Without propagating it, a T1 gene such as CR1 or BIN1
# displays as untiered on a subset of its variants. Carry each gene's best tier
# across all of its rows; genes with no tier anywhere keep the NO_TIER label.
.lv <- levels(dat$top_confidence)
.best <- dat %>%
  filter(!is.na(gene), !top_confidence %in% c(NO_TIER, NO_GENE)) %>%
  group_by(gene) %>%
  summarise(.best_tier = .lv[min(as.integer(top_confidence))], .groups = "drop")
dat <- dat %>%
  left_join(.best, by = "gene") %>%
  mutate(top_confidence = factor(
           ifelse(is.na(.best_tier), as.character(top_confidence), .best_tier),
           levels = .lv)) %>%
  select(-.best_tier)

chr <- intToUtf8
q <- function() intToUtf8(34)
ct_short <- abbreviate(ct_labels, minlength = 3L, named = FALSE)
dot_html_row <- function(i) {
  on <- which(vapply(ct_cols, function(cc) isTRUE(dat[[cc]][i]), logical(1)))
  if (!length(on)) return(paste0(chr(60), "span aria-label=", q(), "No cell type recorded", q(),
    " style=", q(), "color:#cbd5e1", q(), chr(62), "&mdash;", chr(60), "/span", chr(62)))
  paste(sprintf(
    paste0(chr(60), "span class=", q(), "ctpill", q(), " title=", q(), "%s", q(),
           " aria-label=", q(), "%s", q(), " style=", q(),
           "background:%s22;border:1px solid %s;color:#233947", q(), chr(62), "%s",
           chr(60), "/span", chr(62)),
    ct_labels[on], ct_labels[on], ct_colors[on], ct_colors[on], ct_short[on]), collapse = "")
}
dat$ct_dots <- vapply(seq_len(nrow(dat)), dot_html_row, character(1))

# choices
loci  <- sort(unique(dat$ADlocus))
genes <- sort(unique(na.omit(dat$gene)))

# ---- headline counts, derived from the data --------------------------------
# These used to be hardcoded May-2026 figures (183 loci / 410 genes / 145 T1-T4
# and a fixed tier table) and silently went stale when the data was rebuilt for
# the 221-locus release. Derive them from `dat` for the same reason the build
# stamp is derived: the narrative must not be able to contradict the table.
.gt      <- unique(dat[!is.na(dat$gene) & dat$gene != "", c("gene","top_confidence")])
tier_n   <- table(factor(as.character(.gt$top_confidence),
                         levels = c("T1","T2","T3","T4","T5","T6")))
TIERN    <- function(k) format(as.integer(tier_n[[k]]), big.mark = ",")
n_t1_t5  <- sum(as.integer(tier_n[c("T1","T2","T3","T4","T5")]))
n_t1_t6  <- sum(as.integer(tier_n[c("T1","T2","T3","T4","T5","T6")]))
n_t6    <- sum(as.integer(tier_n[c("T6")]))
n_t1_t4  <- sum(as.integer(tier_n[c("T1","T2","T3","T4")]))
n_loci_total <- length(loci)
n_loci_qtl   <- length(unique(dat$ADlocus[!is.na(dat$gene) & dat$gene != ""]))
