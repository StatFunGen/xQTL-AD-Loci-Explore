# ---- badge helpers (HTML) --------------------------------------------------
sig_badge <- function(x) {
  lbl <- c("genome wide"="GW","suggestive"="Sug","ns"="NS")
  lbl[NO_P] <- "no p"   # otherwise this level renders as a bare NA badge
  ifelse(is.na(x), "—", sprintf(
    '<span style="background:%s22;color:%s;border:1px solid %s55;padding:1px 7px;border-radius:10px;font-size:11px;font-weight:600">%s</span>',
    SIG_COL[as.character(x)], SIG_COL[as.character(x)], SIG_COL[as.character(x)], lbl[as.character(x)]))
}

# Pipeline strings arrive machine-formatted (pipe separated, underscore joined).
# Render them as prose for display; the stored value is untouched.
pretty_list <- function(x, sep = ", ") {
  x <- as.character(x)
  if (length(x) != 1 || is.na(x) || !nzchar(x)) return(x)
  v <- trimws(unlist(strsplit(x, "[|;]")))
  v <- gsub("_", " ", v)
  v <- gsub("\\bgwas\\b", "GWAS", v, ignore.case = TRUE)
  v <- gsub("\\bqtl\\b", "QTL", v, ignore.case = TRUE)
  v <- gsub("finemapping", "fine-mapping", v, ignore.case = TRUE)
  v <- v[nzchar(v)]
  if (!length(v)) return("")
  paste(unique(v), collapse = sep)
}

# Explicit evidence states. A blank is not a negative result, so the three
# cases are rendered distinctly and each carries its own tooltip.
ev_mark <- function(x) {
  x <- as.logical(x)
  vapply(seq_along(x), function(i) {
    v <- x[i]
    if (is.na(v))
      return('<span class="ev ev-na" role="img" aria-label="Not available or not tested" title="Not available: this was not assessed for this entry, which is not the same as tested and found negative">&mdash;</span>')
    if (isTRUE(v))
      return('<span class="ev ev-yes" role="img" aria-label="Detected" title="Detected: evidence reported for this gene">&#10003;</span>')
    '<span class="ev ev-ns" role="img" aria-label="Tested but not statistically significant" title="Tested, not significant">NS</span>'
  }, character(1), USE.NAMES = FALSE)
}

CTX_COL <- c(Ast = "#2a78d6", Exc = "#1f7a5c", Inh = "#7b5ea7", Mic = "#da532c",
             OPC = "#a8791b", Oli = "#0f8f9e", bMono = "#b0446f", bMac = "#5e7326",
             bMic = "#8a5fa8", bulk = "#6b7b85")
CTX_FULL <- c(Ast = "Astrocyte", Exc = "Excitatory", Inh = "Inhibitory",
              Mic = "Microglia", OPC = "OPC", Oli = "Oligodendrocyte",
              bMono = "Monocyte", bMac = "Macrophage", bMic = "Microglia (blood)",
              bulk = "Bulk")
CTX_ORD <- c("Ast", "Exc", "Inh", "Mic", "OPC", "Oli", "bMono", "bMac", "bMic", "bulk")
MOD_ORD <- c("eQTL", "sQTL", "pQTL", "mQTL", "caQTL", "haQTL", "gpQTL")
MOD_SHP <- c(eQTL = 21, sQTL = 22, pQTL = 23, mQTL = 24, caQTL = 3, haQTL = 124, gpQTL = 4)
GWS <- -log10(5e-8)

parse_ctx_tokens <- function(d) {
  rx <- "^(\\S+)[[:space:]]+([A-Za-z-]+QTL)([+.-]*)[[:space:]]*\\(T([0-9]),n=([0-9]+)\\)$"
  d <- d[!is.na(d$gene) & d$gene != "" & !is.na(d$ordered_contexts) &
         nzchar(trimws(as.character(d$ordered_contexts))), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    tk <- trimws(unlist(strsplit(as.character(d$ordered_contexts[i]), ";")))
    tk <- tk[nzchar(tk)]
    m  <- regmatches(tk, regexec(rx, tk))
    ok <- vapply(m, length, integer(1)) == 6
    if (!any(ok)) return(NULL)
    m <- m[ok]
    data.frame(gene = as.character(d$gene[i]),
               ctx  = vapply(m, function(z) z[2], character(1)),
               mod  = sub("^[pu]-", "", vapply(m, function(z) z[3], character(1))),
               tier = as.integer(vapply(m, function(z) z[5], character(1))),
               n    = as.integer(vapply(m, function(z) z[6], character(1))),
               stringsAsFactors = FALSE)
  }))
  if (is.null(out) || !nrow(out)) return(NULL)
  out %>% dplyr::group_by(gene, ctx, mod) %>%
    dplyr::summarise(tier = min(tier), n = max(n), .groups = "drop") %>%
    as.data.frame()
}

dl_row <- function(id, file, desc, rows, fmt, size, cols) {
  tags$tr(
    tags$td(class = "f",
      div(class = "fn", file),
      div(class = "de", desc),
      div(class = "cl dc-mono", paste0("Columns: ", cols))),
    tags$td(class = "n dc-mono", rows),
    tags$td(class = "fm", fmt),
    tags$td(class = "n dc-mono", size),
    tags$td(class = "a", downloadLink(id, "Download", class = "dc-btn")))
}

dl_ext <- function(file, desc, fmt, size, href) {
  tags$tr(class = "ext",
    tags$td(class = "f",
      div(class = "fn", file),
      div(class = "de", desc)),
    tags$td(class = "n dc-mono", "—"),
    tags$td(class = "fm", fmt),
    tags$td(class = "n dc-mono", size),
    tags$td(class = "a",
      tags$a(href = href, target = "_blank", rel = "noopener",
             class = "dc-btn", "Catalogue →")))
}

dl_group <- function(title, note, rows) {
  div(class = "dl-grp",
    div(class = "dl-grp-h", span(class = "t", title), span(class = "n", note)),
    tags$table(class = "dl-tbl", tags$tbody(rows)))
}
