# ---- locus registry & search index ---------------------------------------
# Built once at startup. Every path that turns a biological entity into a
# locus goes through these, so the UI no longer depends on the ordinal
# pipeline IDs (chr17_57) that are regenerated with each build.
# Genome build is GRCh38 (stated in the manuscript and in the GRCh38.103
# gene region list); the data files do not carry a build column.
GENOME_BUILD <- "GRCh38"
# Stated, not derived: no data file carries a release or build column.
DATA_RELEASE <- "2026-09"

.chr_norm <- function(v) paste0("chr", sub("^chr", "", as.character(v)))

locus_registry <- dat %>%
  dplyr::filter(!is.na(pos)) %>%
  dplyr::group_by(ADlocus) %>%
  dplyr::summarise(
    chr       = .chr_norm(chr[1]),
    start     = min(pos, na.rm = TRUE),
    end       = max(pos, na.rm = TRUE),
    n_genes   = dplyr::n_distinct(gene[!is.na(gene) & gene != ""]),
    lead_rsid = {
      i <- which.max(ifelse(is.na(log10pval), -Inf, log10pval))
      if (length(i)) as.character(rsid[i]) else NA_character_
    },
    .groups = "drop") %>%
  dplyr::mutate(region = sprintf("%s:%.0f-%.0f", chr, start, end))

gname_map <- local({
  f <- "data/gene_ensembl_map.csv"
  if (!file.exists(f)) return(character(0))
  m <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (is.null(m$GENENAME)) return(character(0))
  m <- m[!is.na(m$GENENAME) & nzchar(m$GENENAME) & !duplicated(m$SYMBOL), ]
  out <- as.character(m$GENENAME); names(out) <- as.character(m$SYMBOL); out
})

gid_map <- local({
  f <- "data/gene_ensembl_map.csv"
  if (!file.exists(f)) return(character(0))
  m <- utils::read.csv(f, stringsAsFactors = FALSE)
  m <- m[!is.na(m$ENSEMBL) & nzchar(m$ENSEMBL) & !duplicated(m$SYMBOL), ]
  out <- as.character(m$ENSEMBL); names(out) <- as.character(m$SYMBOL); out
})

search_index <- local({
  bl <- function(v) is.na(v) | v == ""
  g <- dat %>% dplyr::filter(!bl(gene)) %>%
    dplyr::mutate(
      .tr = suppressWarnings(as.integer(sub("^T", "", as.character(top_confidence)))),
      .tr = ifelse(is.na(.tr), 99L, .tr),
      .lp = ifelse(is.na(log10pval), -Inf, log10pval)) %>%
    dplyr::group_by(term = gene) %>%
    dplyr::arrange(.tr, dplyr::desc(.lp), .by_group = TRUE) %>%
    dplyr::slice(1) %>% dplyr::ungroup() %>%
    dplyr::transmute(term, ADlocus, tier = as.character(top_confidence),
                     kind = "gene")
  v <- dat %>% dplyr::filter(!bl(rsid)) %>%
    dplyr::distinct(term = rsid, ADlocus) %>%
    dplyr::mutate(kind = "variant", tier = NA_character_)
  dplyr::bind_rows(g, v) %>%
    dplyr::left_join(dplyr::select(locus_registry, ADlocus, region), by = "ADlocus") %>%
    dplyr::group_by(term, kind) %>% dplyr::mutate(n_hits = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      value = paste(kind, term, ADlocus, sep = "|"),
      label = term,
      gid = ifelse(kind == "gene" & term %in% names(gid_map),
                   unname(gid_map[term]), "")) %>%
    dplyr::arrange(kind, term) %>%
    as.data.frame()
})

parse_region_str <- function(x) {
  x <- gsub("[ ,]", "", as.character(x))
  m <- regmatches(x, regexec("^(?:chr)?([0-9XYMxym]+):([0-9]+)-([0-9]+)$", x))[[1]]
  if (length(m) != 4) return(NULL)
  list(chr = paste0("chr", m[2]), start = as.numeric(m[3]), end = as.numeric(m[4]))
}

# Returns list(locus, via, n) or NULL. `n` is how many loci the entity matched
# in THIS build, so the caller can tell the user when a saved link is no
# longer unambiguous after a boundary redraw.
resolve_locus <- function(q) {
  hit <- function(loc, via, n) list(locus = loc, via = via, n = as.integer(n))
  if (!is.null(q$gene) && nzchar(q$gene)) {
    h <- search_index[search_index$kind == "gene" &
                        toupper(search_index$term) == toupper(q$gene), ]
    if (nrow(h)) return(hit(h$ADlocus[1], "gene", nrow(h)))
  }
  if (!is.null(q$variant) && nzchar(q$variant)) {
    h <- search_index[search_index$kind == "variant" &
                        search_index$term == q$variant, ]
    if (nrow(h)) return(hit(h$ADlocus[1], "variant", nrow(h)))
  }
  if (!is.null(q$region) && nzchar(q$region)) {
    r <- parse_region_str(q$region)
    if (!is.null(r)) {
      ov <- locus_registry$chr == r$chr &
            locus_registry$start <= r$end & locus_registry$end >= r$start
      if (any(ov)) {
        exact <- ov & locus_registry$start == r$start & locus_registry$end == r$end
        if (any(exact))
          return(hit(locus_registry$ADlocus[which(exact)[1]], "region", sum(ov)))
        w <- pmin(locus_registry$end[ov], r$end) -
             pmax(locus_registry$start[ov], r$start)
        return(hit(locus_registry$ADlocus[which(ov)[which.max(w)]],
                   "region", sum(ov)))
      }
    }
  }
  lp <- if (!is.null(q$locus_pick)) q$locus_pick else q$locus_id
  if (!is.null(lp) && lp %in% locus_registry$ADlocus) return(hit(lp, "legacy", 1))
  NULL
}

locus_to_query <- function(locus, gene = NULL) {
  i <- match(locus, locus_registry$ADlocus)
  if (is.na(i)) return("?")
  p <- c(region = locus_registry$region[i])
  if (!is.na(locus_registry$lead_rsid[i]))
    p["variant"] <- locus_registry$lead_rsid[i]
  if (!is.null(gene) && nzchar(gene)) p["gene"] <- gene
  paste0("?", paste(names(p),
                    vapply(p, utils::URLencode, "", reserved = TRUE),
                    sep = "=", collapse = "&"))
}

# Ordinal blue ramp, light surface. Steps 250-650 of the sequential hue:
# the light end clears 2:1 against the page, and every chip carries its
# label, which is the secondary encoding the adjacent steps require.
conf_pal <- c(T1 = "#78170f", T2 = "#c9452a", T3 = "#f2a184", T4 = "#93c3ef", T5 = "#2a72c4")
conf_txt <- c(T1 = "#ffffff", T2 = "#ffffff", T3 = "#233947", T4 = "#233947", T5 = "#ffffff")
conf_pal["T6"] <- "#10375e"   # dark blue end of the T1-T6 ramp
conf_txt["T6"] <- "#ffffff"
conf_pal[NO_TIER] <- "#f0efec"
conf_txt[NO_TIER] <- "#52514e"
conf_pal[NO_GENE] <- "#f7f7f5"
conf_txt[NO_GENE] <- "#898781"

# Rows with no xQTL gene get a quiet em dash, not a full-width tier chip -
# the category name lives in the filter, where it needs to be readable.
conf_badge <- function(x) ifelse(is.na(x) | as.character(x) == NO_GENE,
  "<span style=\"color:#9aa4b2\" title=\"no xQTL target gene\">\u2014</span>", sprintf(
  '<span style="background:%s;color:%s;padding:1px 7px;border-radius:4px;font-size:11px;font-weight:700;font-family:monospace">%s</span>',
  conf_pal[as.character(x)], conf_txt[as.character(x)], as.character(x)))
bool_badge <- function(x) ifelse(isTRUE(x),
  '<span style="color:#0b6e4f;font-weight:700">●</span>',
  '<span style="color:#e2e8f0">○</span>')

info <- function(label, tip)  # ScotPHO-style inline definition marker
  span(label, tags$button(type = "button", class = "ibub",
    `aria-label` = paste("Explain", if (is.character(label)) label else "this item"),
    onclick = paste0("event.preventDefault();event.stopPropagation();",
                      "this.classList.toggle(", intToUtf8(39), "ibub-open", intToUtf8(39), ");"),
    span(`aria-hidden` = "true", "\u24d8"),
    span(class = "ibub-t", role = "tooltip", tip)))

# Stable URL serialisation for the multi-locus filter. Ordinal locus IDs are
# build-specific, so the address bar carries coordinates instead; ordinals are
# still accepted on the way in so older links keep working.
loci_to_regions <- function(x) {
  i <- match(as.character(x), locus_registry$ADlocus)
  ifelse(is.na(i), as.character(x), locus_registry$region[i])
}

regions_to_loci <- function(x) {
  x <- as.character(x); out <- x
  for (j in seq_along(x)) {
    r <- parse_region_str(x[j])
    if (is.null(r)) next
    k <- which(locus_registry$chr == r$chr &
               locus_registry$start == r$start &
               locus_registry$end == r$end)
    if (!length(k))
      k <- which(locus_registry$chr == r$chr &
                 locus_registry$start <= r$end & locus_registry$end >= r$start)
    if (length(k)) out[j] <- locus_registry$ADlocus[k[1]]
  }
  out
}

# ---- Trans map -------------------------------------------------------------
# Source AD locus and distal target gene, both in genome coordinates, so the
# question "which chromosome acts on which chromosome" can be read directly.
CHR_ORD <- as.character(1:22)
CHR_LEN <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
             159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
             114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
             58617616, 64444167, 46709983, 50818468, 156040895)
CHR_LEN <- CHR_LEN[seq_along(CHR_ORD)]  # autosomes only; no chrX in this release
chr_offsets <- data.frame(
  chr = CHR_ORD, len = CHR_LEN,
  off = c(0, cumsum(as.numeric(CHR_LEN))[-length(CHR_LEN)]),
  stringsAsFactors = FALSE)

gene_pos <- local({
  f <- "gene_positions.csv"
  if (!file.exists(f)) return(NULL)
  g <- readr::read_csv(f, show_col_types = FALSE)
  g$chr <- as.character(g$chr)
  g[!duplicated(g$gene), ]
})

trans_pairs <- local({
  if (is.null(gene_pos)) return(NULL)
  d <- dat[!is.na(dat$trans_genes) & nzchar(trimws(as.character(dat$trans_genes))) &
           !is.na(dat$gene) & dat$gene != "" & !is.na(dat$pos), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  mods <- c(snRNA = "trans_snRNA_n_genes", pQTL = "trans_pQTL_n_genes",
            gpQTL = "trans_gpQTL_n_genes", Hotspot = "trans_hotspot_n_genes")
  out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    tg <- unique(trimws(unlist(strsplit(as.character(d$trans_genes[i]), "[;,|]"))))
    tg <- tg[nzchar(tg)]
    if (!length(tg)) return(NULL)
    gsets <- list(snRNA = "trans_snRNA_genes", pQTL = "trans_pQTL_genes",
                  gpQTL = "trans_gpQTL_genes")
    mm <- character(0); gg <- character(0)
    for (mname in names(gsets)) {
      cc <- gsets[[mname]]
      if (!cc %in% names(d)) next
      vv <- as.character(d[[cc]][i])
      if (is.na(vv) || !nzchar(vv)) next
      g2 <- unique(trimws(unlist(strsplit(vv, "[;,|]"))))
      g2 <- g2[nzchar(g2)]
      if (!length(g2)) next
      gg <- c(gg, g2); mm <- c(mm, rep(mname, length(g2)))
    }
    if (!length(gg)) return(NULL)
    tg <- gg
    data.frame(ADlocus = as.character(d$ADlocus[i]),
               src = as.character(d$gene[i]),
               src_chr = sub("^chr", "", as.character(d$chr[i])),
               src_pos = as.numeric(d$pos[i]),
               tgt = tg,
               modality = mm,
               stringsAsFactors = FALSE)
  }))
  if (is.null(out) || !nrow(out)) return(NULL)
  out <- unique(out)
  i <- match(out$tgt, gene_pos$gene)
  out$tgt_chr <- gene_pos$chr[i]
  out$tgt_pos <- (as.numeric(gene_pos$start[i]) + as.numeric(gene_pos$end[i])) / 2
  out <- out[!is.na(out$tgt_chr) & out$src_chr %in% CHR_ORD & out$tgt_chr %in% CHR_ORD, ]
  if (!nrow(out)) return(NULL)
  out$x <- out$src_pos + chr_offsets$off[match(out$src_chr, chr_offsets$chr)]
  out$y <- out$tgt_pos + chr_offsets$off[match(out$tgt_chr, chr_offsets$chr)]
  out$modality <- factor(out$modality,
                         levels = intersect(c("snRNA", "pQTL", "gpQTL", "Hotspot"),
                                            unique(out$modality)))
  out
})

MOD_COL <- c(snRNA = "#2a78d6", pQTL = "#1baf7a", gpQTL = "#eda100",
             Hotspot = "#da532c", unspecified = "#898781")

trans_map_plot <- function(d, note_empty = "No trans pairs with mapped coordinates.") {
  if (is.null(d) || !nrow(d)) return(NULL)
  o <- chr_offsets
  ggplot(d, aes(x = x, y = y)) +
    geom_vline(xintercept = o$off, colour = LINE, linewidth = 0.25) +
    geom_hline(yintercept = o$off, colour = LINE, linewidth = 0.25) +
    geom_point(aes(colour = modality), size = 2.1, alpha = 0.85) +
    scale_colour_manual(values = MOD_COL, name = NULL, drop = TRUE) +
    scale_x_continuous(breaks = o$off + o$len / 2, labels = o$chr,
                       limits = c(0, max(o$off + o$len)), expand = expansion(mult = 0.005)) +
    scale_y_continuous(breaks = o$off + o$len / 2, labels = o$chr,
                       limits = c(0, max(o$off + o$len)), expand = expansion(mult = 0.005)) +
    labs(x = "AD locus position", y = "Trans target gene position") + g_theme +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          axis.text = element_text(size = 8, colour = MUTE))
}

trans_chr_matrix_plot <- function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  m <- d %>% dplyr::count(src_chr, tgt_chr, name = "n")
  m$src_chr <- factor(m$src_chr, levels = rev(CHR_ORD))
  m$tgt_chr <- factor(m$tgt_chr, levels = CHR_ORD)
  ggplot(m, aes(x = tgt_chr, y = src_chr, fill = n)) +
    geom_tile(colour = PAPER, linewidth = 0.8) +
    scale_fill_gradient(low = "#eaf1fb", high = "#1c5cab", name = "Trans pairs") +
    labs(x = "Target chromosome", y = "AD locus chromosome") + g_theme +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          axis.text = element_text(size = 9, colour = INK))
}

mbp <- function(p) sprintf("%.1f Mb", as.numeric(p) / 1e6)

trans_illu_one <- function(src_chr, src_pos, src_lab, tgt_chr, tgt_pos, tgt_gene,
                           modality, src_gene) {
  W <- 860; PADL <- 54; PADR <- 26
  ls <- as.numeric(CHR_LEN[match(src_chr, CHR_ORD)])
  lt <- as.numeric(CHR_LEN[match(tgt_chr, CHR_ORD)])
  if (is.na(ls) || is.na(lt)) return(NULL)
  span <- W - PADL - PADR
  xs <- PADL + span * min(1, max(0, src_pos / ls))
  xt <- PADL + span * min(1, max(0, tgt_pos / lt))
  col <- unname(MOD_COL[as.character(modality)])
  if (is.na(col) || !nzchar(col)) col <- unname("#95a3ad")
  mid <- (xs + xt) / 2
  lab_anchor <- if (xs > W - 220) "end" else if (xs < 220) "start" else "middle"
  lab_x <- min(W - PADR, max(PADL, xs))
  paste0(
    '<svg class="ti" viewBox="0 0 860 148" preserveAspectRatio="xMidYMid meet" role="img">',
  '<text class="ti-lab" x="', sprintf('%.0f', lab_x), '" y="20" text-anchor="', lab_anchor, '">',
      htmltools::htmlEscape(paste0(src_gene, " \u00b7 ", src_lab)), '</text>',
    '<text class="ti-chr" x="8" y="45">chr', src_chr, '</text>',
    '<rect class="ti-bar" x="', PADL, '" y="34" width="', span, '" height="13" rx="6.5"/>',
    '<rect class="ti-hit ti-src" x="', sprintf('%.1f', xs - 3.5), '" y="31" width="7" height="19" rx="3"/>',
    '<circle class="ti-dot ti-src" cx="', sprintf('%.1f', xs), '" cy="27" r="4.6"/>',
    '<path class="ti-link" style="stroke:', col, '" d="M', sprintf('%.1f', xs), ',53 C',
      sprintf('%.1f', xs), ',70 ', sprintf('%.1f', xt), ',80 ', sprintf('%.1f', xt), ',95"/>',
    '<text class="ti-tag" style="fill:', col, '" x="', sprintf('%.0f', PADL),
      '" y="78" text-anchor="start">trans association',
      if (!is.na(modality) && modality != 'assay not recorded')
        paste0(' &#183; ', htmltools::htmlEscape(as.character(modality))) else '', '</text>',
    '<text class="ti-chr" x="8" y="115">chr', tgt_chr, '</text>',
    '<rect class="ti-bar" x="', PADL, '" y="104" width="', span, '" height="13" rx="6.5"/>',
    '<rect class="ti-hit" style="fill:', col, '" x="', sprintf('%.1f', xt - 3.5),
      '" y="101" width="7" height="19" rx="3"/>',
    '<text class="ti-gene" style="fill:', col, '" x="', sprintf('%.0f', xt),
      '" y="138" text-anchor="middle">', htmltools::htmlEscape(tgt_gene), '</text>',
    '<text class="ti-pos" x="', sprintf('%.0f', xt), '" y="138" text-anchor="middle" dy="0">',
      '</text>',
    '</svg>')
}

CEN_FRAC <- c(0.503,0.375,0.455,0.264,0.268,0.354,0.377,0.313,0.318,0.297,0.395,
              0.268,0.155,0.166,0.186,0.408,0.269,0.221,0.437,0.436,0.269,0.283,0.30)
names(CEN_FRAC) <- CHR_ORD

trans_card_svg <- function(src_chr, src_pos, src_lab, tg) {
  W <- 1000; X0 <- 156; XW <- 820
  L1 <- as.numeric(CHR_LEN[match("1", CHR_ORD)])
  chrs <- unique(tg$tgt_chr)
  chrs <- chrs[order(match(chrs, CHR_ORD))]
  H <- 96 + length(chrs) * 84
  bw <- function(ch) XW * as.numeric(CHR_LEN[match(ch, CHR_ORD)]) / L1
  bx <- function(ch, p) X0 + XW * (as.numeric(p) / L1)
  bar <- function(ch, y, lab, sub) {
    w <- bw(ch); cf <- CEN_FRAC[[ch]]
    paste0(
      '<text class="tc-chr" x="', X0 - 14, '" y="', y + 13, '" text-anchor="end">chr', ch, '</text>',
      '<text class="tc-role" x="', X0 - 14, '" y="', y + 26, '" text-anchor="end">', sub, '</text>',
      '<rect class="tc-bar" x="', X0, '" y="', y, '" width="', round(w, 1),
      '" height="20" rx="10"/>',
      '<rect class="tc-cen" x="', round(X0 + w * cf - 2.5, 1), '" y="', y,
      '" width="5" height="20"/>')
  }
  ys <- 46
  out <- bar(src_chr, ys, src_lab, "AD locus")
  xs <- bx(src_chr, src_pos)
  out <- paste0(out,
    '<rect class="tc-mark tc-src" x="', round(xs - 1.5, 1), '" y="', ys - 5,
    '" width="3" height="30"/>',
    '<text class="tc-lab tc-srclab" x="', round(xs, 1), '" y="', ys - 12,
    '" text-anchor="middle">', htmltools::htmlEscape(src_lab), '</text>')
  for (i in seq_along(chrs)) {
    ch <- chrs[i]; yy <- 96 + (i - 1) * 84 + 10
    out <- paste0(out, bar(ch, yy, ch, "distal target"))
    sub <- tg[tg$tgt_chr == ch, , drop = FALSE]
    sub <- sub[order(sub$tgt_pos), , drop = FALSE]
    for (j in seq_len(nrow(sub))) {
      col <- unname(MOD_COL[as.character(sub$modality[j])])
      if (is.na(col) || !nzchar(col)) col <- "#6b7b85"
      xt <- bx(ch, sub$tgt_pos[j])
      out <- paste0(out,
        '<path class="tc-link" style="stroke:', col, '" d="M', round(xs, 1), ',', ys + 27,
        ' C', round(xs, 1), ',', yy - 26, ' ', round(xt, 1), ',', yy - 34, ' ',
        round(xt, 1), ',', yy - 6, '"/>',
        '<rect class="tc-mark" style="fill:', col, '" x="', round(xt - 1.5, 1),
        '" y="', yy - 5, '" width="3" height="30"/>',
        '<text class="tc-gene" style="fill:', col, '" x="', round(xt, 1), '" y="', yy + 38,
        '" text-anchor="middle">', htmltools::htmlEscape(sub$tgt[j]), '</text>',
        '<text class="tc-pos" x="', round(xt, 1), '" y="', yy + 50,
        '" text-anchor="middle">', mbp(sub$tgt_pos[j]), '</text>')
    }
  }
  paste0('<svg class="tc-svg" viewBox="0 0 ', W, ' ', H,
         '" preserveAspectRatio="xMidYMid meet" role="img">', out, '</svg>')
}

TIER_SEQ <- c("T1", "T2", "T3", "T4", "T5", "T6")

locus_table <- local({
  rows <- lapply(seq_len(nrow(locus_registry)), function(i) {
    lo <- as.character(locus_registry$ADlocus[i])
    d  <- dat[dat$ADlocus == lo, , drop = FALSE]
    g  <- d[!is.na(d$gene) & d$gene != "", , drop = FALSE]
    e  <- parse_ctx_tokens(d)
    gu <- if (nrow(g)) g[!duplicated(g$gene), , drop = FALSE] else g
    tv <- as.character(gu$top_confidence)
    k  <- match(tv, TIER_SEQ)
    best <- if (any(!is.na(k))) TIER_SEQ[min(k, na.rm = TRUE)] else NA_character_
    top <- NA_character_
    if (nrow(gu)) {
      ord <- order(ifelse(is.na(k), 99L, k),
                   -ifelse(is.na(gu$log10pval), -Inf, gu$log10pval))
      top <- as.character(gu$gene[ord[1]])
    }
    cts <- if (is.null(e)) character(0) else intersect(CTX_ORD, unique(e$ctx))
    mds <- if (is.null(e)) character(0) else intersect(MOD_ORD, unique(e$mod))
    tc  <- vapply(TIER_SEQ, function(tt) sum(tv == tt, na.rm = TRUE), integer(1))
    data.frame(
      ADlocus = lo,
      region  = as.character(locus_registry$region[i]),
      chr     = as.character(locus_registry$chr[i]),
      lead    = as.character(locus_registry$lead_rsid[i]),
      best    = best,
      top_gene = top,
      cells   = paste(cts, collapse = ","),
      mods    = paste(mds, collapse = ", "),
      n_genes = nrow(gu),
      n_rows  = nrow(d),
      t1 = tc[[1]], t2 = tc[[2]], t3 = tc[[3]],
      t4 = tc[[4]], t5 = tc[[5]], t6 = tc[[6]],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$best_k <- ifelse(is.na(match(out$best, TIER_SEQ)), 99L, match(out$best, TIER_SEQ))
  out[order(out$best_k, -out$n_genes), , drop = FALSE]
})

lt_cells <- function(s) {
  if (is.na(s) || !nzchar(s)) return("<span class=\"lt-na\">&mdash;</span>")
  v <- strsplit(s, ",")[[1]]
  paste0('<span class="lt-dots">',
    paste0(vapply(v, function(k) {
      col <- unname(CTX_COL[k]); if (is.na(col)) col <- "#6b7b85"
      nm  <- if (k %in% names(CTX_FULL)) CTX_FULL[[k]] else k
      sprintf('<span class="lt-dot" style="background:%s" title="%s"></span>', col, nm)
    }, character(1)), collapse = ""), '</span>')
}

lt_bar <- function(r) {
  v <- as.numeric(r[paste0("t", 1:6)])
  tot <- sum(v)
  if (!is.finite(tot) || tot <= 0) return("<span class=\"lt-na\">&mdash;</span>")
  seg <- vapply(seq_len(6), function(i) {
    if (v[i] <= 0) return("")
    sprintf('<span style="width:%.2f%%;background:%s" title="%s: %d genes"></span>',
            100 * v[i] / tot, unname(conf_pal[[TIER_SEQ[i]]]), TIER_SEQ[i], as.integer(v[i]))
  }, character(1))
  paste0('<span class="lt-ebar">', paste0(seg, collapse = ""), '</span>')
}

ctx_tier_counts <- local({
  e <- parse_ctx_tokens(dat)
  if (is.null(e) || !nrow(e)) return(NULL)
  g <- e %>% dplyr::group_by(ctx, gene) %>%
    dplyr::summarise(tier = min(tier), .groups = "drop") %>% as.data.frame()
  m <- matrix(0L, nrow = length(CTX_ORD), ncol = 6,
              dimnames = list(CTX_ORD, TIER_SEQ))
  for (i in seq_len(nrow(g))) {
    r <- match(g$ctx[i], CTX_ORD); k <- as.integer(g$tier[i])
    if (!is.na(r) && !is.na(k) && k >= 1 && k <= 6) m[r, k] <- m[r, k] + 1L
  }
  m[order(-rowSums(m)), , drop = FALSE]
})

locus_tier_counts <- local({
  b <- locus_table$best
  v <- vapply(TIER_SEQ, function(k) sum(!is.na(b) & b == k), integer(1))
  c(v, "No tier" = sum(is.na(b)))
})

mod_counts <- local({
  e <- parse_ctx_tokens(dat)
  if (is.null(e) || !nrow(e)) return(NULL)
  v <- vapply(MOD_ORD, function(m) sum(e$mod == m), integer(1))
  v <- v[v > 0]
  v[order(-v)]
})

trans_chr_counts <- local({
  if (is.null(trans_pairs) || !nrow(trans_pairs)) return(NULL)
  t <- unique(trans_pairs[, c("ADlocus", "tgt", "tgt_chr")])
  v <- table(factor(t$tgt_chr, levels = CHR_ORD))
  v <- as.integer(v); names(v) <- CHR_ORD
  v[v > 0]
})

trans_hub_counts <- local({
  if (is.null(trans_pairs) || !nrow(trans_pairs)) return(NULL)
  t <- unique(trans_pairs[, c("ADlocus", "tgt")])
  v <- table(t$ADlocus)
  v <- sort(v, decreasing = TRUE)
  v <- v[seq_len(min(12, length(v)))]
  out <- as.integer(v); names(out) <- names(v)
  out
})


gene_ctx_breadth <- local({
  e <- parse_ctx_tokens(dat)
  if (is.null(e) || !nrow(e)) return(NULL)
  g <- unique(e[, c("gene", "ctx")])
  v <- table(table(g$gene))
  out <- as.integer(v)
  names(out) <- ifelse(names(v) == "1", "1 cell type", paste(names(v), "cell types"))
  out
})

ctx_unique_counts <- local({
  e <- parse_ctx_tokens(dat)
  if (is.null(e) || !nrow(e)) return(NULL)
  g <- unique(e[, c("gene", "ctx")])
  n <- table(g$gene)
  u <- g[g$gene %in% names(n)[n == 1], , drop = FALSE]
  v <- table(factor(u$ctx, levels = CTX_ORD))
  out <- as.integer(v); names(out) <- names(v)
  out[out > 0]
})
