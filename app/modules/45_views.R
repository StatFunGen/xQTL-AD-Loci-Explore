# ---- shared view helpers: cell types, trans pairs, batch lookup, quality ----

CT_COLS <- c(ct_Exc_xQTL = "Excitatory neuron",
             ct_Inh_xQTL = "Inhibitory neuron",
             ct_Ast_xQTL = "Astrocyte",
             ct_Microglia_xQTL = "Microglia",
             ct_Oli_xQTL = "Oligodendrocyte",
             ct_OPC_xQTL = "OPC",
             ct_Bulk_Immune_xQTL = "Bulk immune",
             ct_Brain_xQTL = "Bulk brain")
CT_COLS <- CT_COLS[names(CT_COLS) %in% names(dat)]
CT_CHOICES <- stats::setNames(names(CT_COLS), unname(CT_COLS))

# ---- explicit context lookup -----------------------------------------------
# The raw `context` strings use several inconsistent token orders (cell type
# first, dataset first, modality in the middle, or no modality at all), so
# cell type, dataset and modality cannot be read off token position. This
# table is exhaustive over the distinct context values in the release.
# Each entry is c(label, dataset, modality, cell-type bucket).
CTX_MAP <- list(
  "AC_CUIMC1_eQTL" = c("Anterior cingulate", "CUIMC1", "eQTL", "ct_Brain_xQTL"),
  "Ast_CUIMC1_eQTL" = c("Astrocyte", "CUIMC1", "eQTL", "ct_Ast_xQTL"),
  "Ast_MIT_eQTL" = c("Astrocyte", "MIT", "eQTL", "ct_Ast_xQTL"),
  "Ast_mega_eQTL" = c("Astrocyte", "mega", "eQTL", "ct_Ast_xQTL"),
  "BM_10_MSBB_eQTL" = c("BM 10", "MSBB", "eQTL", "ct_Brain_xQTL"),
  "BM_22_MSBB_eQTL" = c("BM 22", "MSBB", "eQTL", "ct_Brain_xQTL"),
  "BM_36_MSBB_eQTL" = c("BM 36", "MSBB", "eQTL", "ct_Brain_xQTL"),
  "BM_44_MSBB_eQTL" = c("BM 44", "MSBB", "eQTL", "ct_Brain_xQTL"),
  "DLPFC_Bennett_pQTL" = c("DLPFC", "Bennett", "pQTL", "ct_Brain_xQTL"),
  "DLPFC_CUIMC1_eQTL" = c("DLPFC", "CUIMC1", "eQTL", "ct_Brain_xQTL"),
  "DLPFC_Klein_gpQTL_adjusted" = c("DLPFC", "Klein", "gpQTL", "ct_Brain_xQTL"),
  "DLPFC_Klein_gpQTL_unadjusted" = c("DLPFC", "Klein", "gpQTL", "ct_Brain_xQTL"),
  "Exc_CUIMC1_eQTL" = c("Excitatory", "CUIMC1", "eQTL", "ct_Exc_xQTL"),
  "Exc_MIT_eQTL" = c("Excitatory", "MIT", "eQTL", "ct_Exc_xQTL"),
  "Exc_mega_eQTL" = c("Excitatory", "mega", "eQTL", "ct_Exc_xQTL"),
  "Inh_CUIMC1_eQTL" = c("Inhibitory", "CUIMC1", "eQTL", "ct_Inh_xQTL"),
  "Inh_MIT_eQTL" = c("Inhibitory", "MIT", "eQTL", "ct_Inh_xQTL"),
  "Inh_mega_eQTL" = c("Inhibitory", "mega", "eQTL", "ct_Inh_xQTL"),
  "KnightADRC_mQTL" = c("Brain", "KnightADRC", "mQTL", "ct_Brain_xQTL"),
  "Knight_eQTL_brain" = c("Brain", "Knight", "eQTL", "ct_Brain_xQTL"),
  "Knight_pQTL_brain" = c("Brain", "Knight", "pQTL", "ct_Brain_xQTL"),
  "MSBB_BM36_pQTL" = c("BM 36", "MSBB", "pQTL", "ct_Brain_xQTL"),
  "MSBB_mQTL" = c("Brain", "MSBB", "mQTL", "ct_Brain_xQTL"),
  "Metabrain_Cerebellum" = c("Cerebellum", "Metabrain", "eQTL", "ct_Brain_xQTL"),
  "Metabrain_Cortex" = c("Cortex", "Metabrain", "eQTL", "ct_Brain_xQTL"),
  "Metabrain_Hippocampus" = c("Hippocampus", "Metabrain", "eQTL", "ct_Brain_xQTL"),
  "MiGA_GFM_eQTL" = c("GFM", "MiGA", "eQTL", "ct_Brain_xQTL"),
  "MiGA_GTS_eQTL" = c("GTS", "MiGA", "eQTL", "ct_Brain_xQTL"),
  "MiGA_SVZ_eQTL" = c("SVZ", "MiGA", "eQTL", "ct_Brain_xQTL"),
  "MiGA_THA_eQTL" = c("THA", "MiGA", "eQTL", "ct_Brain_xQTL"),
  "Mic_CUIMC1_eQTL" = c("Microglia", "CUIMC1", "eQTL", "ct_Microglia_xQTL"),
  "Mic_MIT_eQTL" = c("Microglia", "MIT", "eQTL", "ct_Microglia_xQTL"),
  "Mic_mega_eQTL" = c("Microglia", "mega", "eQTL", "ct_Microglia_xQTL"),
  "OPC_CUIMC1_eQTL" = c("OPC", "CUIMC1", "eQTL", "ct_OPC_xQTL"),
  "OPC_MIT_eQTL" = c("OPC", "MIT", "eQTL", "ct_OPC_xQTL"),
  "OPC_mega_eQTL" = c("OPC", "mega", "eQTL", "ct_OPC_xQTL"),
  "Oli_CUIMC1_eQTL" = c("Oligodendrocyte", "CUIMC1", "eQTL", "ct_Oli_xQTL"),
  "Oli_MIT_eQTL" = c("Oligodendrocyte", "MIT", "eQTL", "ct_Oli_xQTL"),
  "Oli_mega_eQTL" = c("Oligodendrocyte", "mega", "eQTL", "ct_Oli_xQTL"),
  "PCC_CUIMC1_eQTL" = c("PCC", "CUIMC1", "eQTL", "ct_Brain_xQTL"),
  "ROSMAP_AC_p_sQTL" = c("Anterior cingulate", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_AC_sQTL" = c("Anterior cingulate", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_AC_u_sQTL" = c("Anterior cingulate", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_DLPFC_haQTL" = c("DLPFC", "ROSMAP", "haQTL", "ct_Brain_xQTL"),
  "ROSMAP_DLPFC_mQTL" = c("DLPFC", "ROSMAP", "mQTL", "ct_Brain_xQTL"),
  "ROSMAP_DLPFC_p_sQTL" = c("DLPFC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_DLPFC_sQTL" = c("DLPFC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_DLPFC_u_sQTL" = c("DLPFC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_PCC_p_sQTL" = c("PCC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_PCC_sQTL" = c("PCC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "ROSMAP_PCC_u_sQTL" = c("PCC", "ROSMAP", "sQTL", "ct_Brain_xQTL"),
  "STARNET_eQTL_Mac" = c("Macrophage", "STARNET", "eQTL", "ct_Bulk_Immune_xQTL"),
  "monocyte_ROSMAP_eQTL" = c("Monocyte", "ROSMAP", "eQTL", "ct_Bulk_Immune_xQTL")
)

.ctx_bucket <- vapply(as.character(dat$context), function(v) {
  if (is.na(v) || !nzchar(v)) return(NA_character_)
  m <- CTX_MAP[[v]]
  if (is.null(m)) NA_character_ else m[[4]]
}, character(1), USE.NAMES = FALSE)

# rows whose own assay context belongs to one cell type or region group
ct_rows <- function(col) {
  if (is.null(col) || !nzchar(col)) return(dat[0, ])
  dat[!is.na(.ctx_bucket) & .ctx_bucket == col, , drop = FALSE]
}

# ---- trans pairs, one row per source-target ----
.mk_trans_tbl <- function() {
  mods <- c(snRNA = "trans_snRNA", pQTL = "trans_pQTL",
            gpQTL = "trans_gpQTL", Hotspot = "trans_hotspot")
  res <- list()
  for (m in names(mods)) {
    p  <- unname(mods[[m]])
    tc <- if (m == "Hotspot") paste0(p, "_programs") else paste0(p, "_genes")
    cc <- paste0(p, "_contexts")
    nc <- paste0(p, "_n_genes")
    if (!tc %in% names(dat)) next
    keep <- !is.na(dat[[tc]]) & nzchar(dat[[tc]])
    if (!any(keep)) next
    d   <- dat[keep, , drop = FALSE]
    lst <- strsplit(d[[tc]], "\\s*,\\s*")
    n   <- lengths(lst)
    res[[m]] <- data.frame(
      locus    = rep(d$ADlocus, n),
      source   = rep(d$gene, n),
      rsid     = rep(d$rsid, n),
      target   = trimws(unlist(lst)),
      modality = m,
      context  = if (cc %in% names(d)) rep(d[[cc]], n) else NA_character_,
      n_genes  = if (nc %in% names(d)) rep(d[[nc]], n) else NA_integer_,
      tier     = rep(as.character(d$top_confidence), n),
      stringsAsFactors = FALSE)
  }
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res)
  out <- out[!is.na(out$target) & nzchar(out$target), , drop = FALSE]
  rownames(out) <- NULL
  unique(out)
}
trans_all <- .mk_trans_tbl()

# ---- batch lookup: paste genes, rsIDs or regions ----
batch_tokens <- function(txt, limit = 300) {
  if (is.null(txt)) return(character(0))
  x <- unlist(strsplit(txt, "[[:space:],;]+"))
  x <- trimws(x)
  x <- x[nzchar(x)]
  utils::head(unique(x), limit)
}

# data.csv gene_id does not line up with the gene column, so Ensembl IDs are
# resolved through gid_map (the release symbol-to-Ensembl map) instead.
ens_to_symbol <- function(tk) {
  base <- sub("\\..*$", "", toupper(trimws(tk)))
  if (!exists("gid_map") || !length(gid_map)) return(character(0))
  hit <- names(gid_map)[!is.na(gid_map) & toupper(sub("\\..*$", "", gid_map)) == base]
  unique(hit)
}

batch_kind <- function(tok) {
  ifelse(grepl("^ENSG[0-9]+", tok, ignore.case = TRUE), "Ensembl ID",
  ifelse(grepl("^rs[0-9]+$", tok, ignore.case = TRUE), "rsID",
  ifelse(grepl("^chr[0-9XYM]+[:_-][0-9]+([-_][0-9]+)?$", tok, ignore.case = TRUE), "region",
         "gene")))
}

.region_parts <- function(tok) {
  p <- regmatches(tok, regexec("^chr([0-9XYM]+)[:_-]([0-9]+)(?:[-_]([0-9]+))?$", tok, ignore.case = TRUE))[[1]]
  if (!length(p)) return(NULL)
  ch <- suppressWarnings(as.integer(p[2]))
  a  <- suppressWarnings(as.numeric(p[3]))
  b  <- if (nzchar(p[4])) suppressWarnings(as.numeric(p[4])) else a
  if (is.na(ch) || is.na(a)) return(NULL)
  list(chr = ch, start = min(a, b), end = max(a, b))
}

# one row per query token, with the matching data rows attached
batch_match <- function(tokens) {
  if (!length(tokens)) return(NULL)
  out <- list()
  for (tk in tokens) {
    k <- batch_kind(tk)
    hit <- dat[0, ]
    if (k == "gene") {
      hit <- dat[!is.na(dat$gene) & toupper(dat$gene) == toupper(tk), ]
    } else if (k == "Ensembl ID") {
      sym <- ens_to_symbol(tk)
      hit <- if (length(sym)) dat[!is.na(dat$gene) & toupper(dat$gene) %in% toupper(sym), ] else dat[0, ]
    } else if (k == "rsID") {
      hit <- dat[!is.na(dat$rsid) & tolower(dat$rsid) == tolower(tk), ]
    } else if (k == "region") {
      rp <- .region_parts(tk)
      if (!is.null(rp)) hit <- dat[!is.na(dat$chr) & dat$chr == rp$chr &
                                   !is.na(dat$pos) & dat$pos >= rp$start & dat$pos <= rp$end, ]
    }
    out[[length(out) + 1]] <- data.frame(
      query = tk, kind = k,
      matched = nrow(hit) > 0,
      loci  = if (nrow(hit)) paste(sort(unique(hit$ADlocus)), collapse = ", ") else NA_character_,
      n_loci = length(unique(hit$ADlocus)),
      genes = if (nrow(hit)) paste(sort(unique(hit$gene)), collapse = ", ") else NA_character_,
      n_genes = length(unique(hit$gene[nzchar(hit$gene)])),
      best_tier = if (nrow(hit) && any(!is.na(hit$top_confidence)))
                    as.character(sort(hit$top_confidence)[1]) else NA_character_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# the underlying data rows for a set of tokens, for the detail table and download
batch_rows <- function(tokens) {
  if (!length(tokens)) return(dat[0, ])
  idx <- logical(nrow(dat))
  for (tk in tokens) {
    k <- batch_kind(tk)
    if (k == "gene") {
      idx <- idx | (!is.na(dat$gene) & toupper(dat$gene) == toupper(tk))
    } else if (k == "Ensembl ID") {
      sym <- ens_to_symbol(tk)
      if (length(sym)) idx <- idx | (!is.na(dat$gene) & toupper(dat$gene) %in% toupper(sym))
    } else if (k == "rsID") {
      idx <- idx | (!is.na(dat$rsid) & tolower(dat$rsid) == tolower(tk))
    } else if (k == "region") {
      rp <- .region_parts(tk)
      if (!is.null(rp)) idx <- idx | (!is.na(dat$chr) & dat$chr == rp$chr &
                                      !is.na(dat$pos) & dat$pos >= rp$start & dat$pos <= rp$end)
    }
  }
  dat[idx, , drop = FALSE]
}

# ---- Cell types panel ----
panel_cell_types <- function() {
  nav_panel(
    "Cell types", icon = icon("vials"),
    div(class = "tc-top",
        h1(class = "doc-h", style = "font-size:22px;margin:0", "Evidence by cell type"),
        div(class = "dc-sub", style = "max-width:820px",
            "Everything the release records for one brain or immune cell type: which AD loci ",
            "speak to it, which genes are named in it, and how strong that evidence is. ",
            "The same data as the Loci tab, read down the cell-type axis instead.")),
    div(class = "tc-bar",
        div(class = "tc-f",
            selectInput("cty_ct", info("Cell type",
              "Cell types come from the cell-type resolved xQTL columns. Bulk brain and bulk immune are tissue-level, not single cell."),
              choices = CT_CHOICES, width = "100%")),
        div(class = "tc-f tc-f2",
            selectizeInput("cty_tier", info("Evidence tier",
              "Optional. Leave empty for every tier. Works on its own or together with the search box."),
              choices = .lv, multiple = TRUE, width = "100%",
              options = list(placeholder = "Any tier"))),
        div(class = "tc-f tc-f3",
            textInput("cty_q", info("Find",
              "Matches gene symbol, AD locus or rsID. Works on its own or on top of the tier filter."),
              value = "", width = "100%", placeholder = "Gene, locus or rsID"))),
    uiOutput("cty_active"),
    uiOutput("cty_stats"),
    div(class = "panel", style = "margin-top:18px",
        div(class = "dc-sechead",
            div(class = "dc-figt", "Evidence tiers in this cell type"),
            div(class = "dc-figc", "How many gene records reach each tier here.")),
        withSpinner(plotOutput("cty_tier_plot", height = "210px"),
                    type = 8, color = "#2a78d6", size = 0.6)),
    div(class = "panel", style = "margin-top:20px",
        div(class = "dc-sechead",
            div(class = "dc-figt", "Gene records"),
            div(class = "dc-figc", "One row per gene and AD locus with evidence in this cell type.")),
        withSpinner(DTOutput("cty_tbl"), type = 8, color = "#2a78d6", size = 0.6),
        div(class = "dc-exp", style = "margin-top:10px",
            downloadLink("dl_cty", "CSV", class = "dc-btn"))))
}

# ---- Batch panel: paste a list of genes, rsIDs or regions ----
panel_batch <- function() {
  nav_panel(
    "Batch", icon = icon("clipboard-list"),
    div(class = "tc-top",
        h1(class = "doc-h", style = "font-size:22px;margin:0", "Check a list against this release"),
        div(class = "dc-sub", style = "max-width:820px",
            "Paste a candidate list from your own screen and see which entries this release ",
            "supports. Gene symbols, Ensembl IDs, rsIDs and genomic regions can be mixed in ",
            "one list, separated by spaces, commas or new lines.")),
    div(class = "bt-wrap",
        div(class = "bt-in",
            textAreaInput("bt_in", info("Your list",
              "Up to 300 entries. A region is written chr2:127088929-127118906. Entries that match nothing are reported back rather than dropped."),
              value = "", rows = 7, width = "100%",
              placeholder = "BIN1\nAPOE\nrs429358\nchr2:127088929-127118906\nENSG00000130203"),
            div(class = "bt-btns",
                actionButton("bt_go", "Match list", class = "dc-btn dc-btn-go"),
                actionLink("bt_demo", "Load an example", class = "gp-clear"),
                actionLink("bt_clear", "Clear", class = "gp-clear"))),
        div(class = "bt-side", uiOutput("bt_summary"))),
    uiOutput("bt_body"))
}

# body is rendered server-side so the page stays empty until a list is matched
batch_body_ui <- function() {
  tagList(
    div(class = "panel", style = "margin-top:20px",
        div(class = "dc-sechead",
            div(class = "dc-figt", "Your list, entry by entry"),
            div(class = "dc-figc", "What each entry resolved to, including the ones that matched nothing.")),
        DTOutput("bt_query_tbl")),
    div(class = "panel", style = "margin-top:20px",
        div(class = "dc-sechead",
            div(class = "dc-figt", "Evidence tiers across the matched entries"),
            div(class = "dc-figc", "Each matched gene record counted once, at the tier it reached.")),
        plotOutput("bt_tier_plot", height = "210px")),
    div(class = "panel", style = "margin-top:20px",
        div(class = "dc-sechead",
            div(class = "dc-figt", "Matched gene records"),
            div(class = "dc-figc", "Every release row behind the summary above.")),
        DTOutput("bt_rows_tbl"),
        div(class = "dc-exp", style = "margin-top:10px",
            downloadLink("dl_batch", "CSV", class = "dc-btn"))))
}

# ---- context strings split into their three parts ----------------------
# Contexts are written like Inh_DeJager_eQTL or BM_10_MSBB_eQTL: the assay is
# the token ending in QTL, the dataset is the token before it, and whatever
# comes first is the cell type or brain region.
# the release writes these datasets by first author; the site shows the site name
DSET_FULL <- c(Kellis = "MIT", DeJager = "CUIMC1")

.ctx_one <- function(v) {
  if (is.na(v) || !nzchar(v)) return(c(NA_character_, NA_character_, NA_character_))
  m <- CTX_MAP[[v]]
  if (!is.null(m)) return(m[1:3])
  # unknown context: show it whole rather than splitting it on a guess
  c(v, NA_character_, NA_character_)
}

ctx_cols <- function(v) {
  v <- as.character(v)
  pick <- function(one, slot) {
    if (is.na(one) || !nzchar(one)) return(NA_character_)
    parts <- trimws(unlist(strsplit(one, ",")))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NA_character_)
    u <- vapply(parts, function(p) .ctx_one(p)[slot], character(1), USE.NAMES = FALSE)
    u <- unique(u[!is.na(u) & nzchar(u)])
    if (!length(u)) NA_character_ else paste(u, collapse = ", ")
  }
  list(ctx  = vapply(v, pick, character(1), 1L, USE.NAMES = FALSE),
       dset = vapply(v, pick, character(1), 2L, USE.NAMES = FALSE),
       mod  = vapply(v, pick, character(1), 3L, USE.NAMES = FALSE))
}

# rewrite the dataset token inside a raw context string, for the places that
# still show them whole
ds_relabel <- function(x) {
  x <- as.character(x)
  for (k in names(DSET_FULL))
    x <- gsub(paste0("_", k, "_"), paste0("_", unname(DSET_FULL[[k]]), "_"), x, fixed = TRUE)
  x
}
