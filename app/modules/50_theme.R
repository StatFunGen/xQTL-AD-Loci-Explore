# ---- ggplot theme ----------------------------------------------------------
g_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major  = element_line(color = "#232a33", linewidth = 0.3),
        plot.background   = element_rect(fill = "#151b23", color = NA),
        panel.background  = element_rect(fill = "#151b23", color = NA),
        legend.background = element_rect(fill = "#151b23", color = NA),
        legend.key        = element_rect(fill = "#151b23", color = NA),
        legend.text       = element_text(color = "#c3c2b7"),
        legend.title      = element_text(color = MUTE),
        plot.title = element_text(color = INK, face = "bold", size = 13),
        axis.text  = element_text(color = MUTE),
        axis.title = element_text(color = MUTE))

# Compact row panel: no axes, category label in ink, value at the bar end.
# Used by the four supporting dashboard panels so they read as a readout
# rather than four competing charts.
rowbar <- function(d, fillcol = "#3987e5", pal = NULL) {
  d <- as.data.frame(d)
  if (!nrow(d)) return(NULL)
  d$lab <- factor(as.character(d$lab), levels = rev(as.character(d$lab)))
  p <- ggplot(d, aes(n, lab))
  p <- if (is.null(pal)) p + geom_col(fill = fillcol, width = 0.55)
       else p + geom_col(aes(fill = lab), width = 0.55) +
              scale_fill_manual(values = pal, guide = "none")
  p + geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.25,
                size = 3.1, colour = MUTE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(x = NULL, y = NULL) + g_theme +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(colour = INK, size = 9.5),
          panel.grid   = element_blank(),
          plot.margin  = margin(4, 10, 2, 2))
}


# ---- bslib theme -----------------------------------------------------------

dk_theme <- function(on) {
  if (!isTRUE(on)) return(NULL)
  ggplot2::theme(
    plot.background   = ggplot2::element_rect(fill = "#15202a", colour = NA),
    panel.background  = ggplot2::element_rect(fill = "#15202a", colour = NA),
    legend.background = ggplot2::element_rect(fill = "#15202a", colour = NA),
    legend.key        = ggplot2::element_rect(fill = "#15202a", colour = NA),
    legend.box.background = ggplot2::element_rect(fill = "#15202a", colour = NA),
    text        = ggplot2::element_text(colour = "#e8eef2"),
    axis.text   = ggplot2::element_text(colour = "#95a5af"),
    axis.title  = ggplot2::element_text(colour = "#95a5af"),
    legend.text = ggplot2::element_text(colour = "#95a5af"),
    axis.line.x = ggplot2::element_line(colour = "#2b3a45"))
}

# ---- Home river (cell type -> assay -> evidence tier) ----------------------
.sk_data <- local({
  e <- parse_ctx_tokens(dat)
  if (is.null(e) || !nrow(e)) return(NULL)
  e <- e[!is.na(e$ctx) & !is.na(e$mod) & !is.na(e$tier), , drop = FALSE]
  e$ctxf <- unname(CTX_FULL[e$ctx])
  keep <- c("eQTL", "sQTL", "pQTL", "mQTL")
  e$modg <- ifelse(e$mod %in% keep, e$mod, "other assays")
  e$tierl <- paste0("T", e$tier)
  e
})

.g2l <- local({
  if (is.null(dat$gene) || is.null(dat$ADlocus)) return(NULL)
  d <- unique(data.frame(gene = as.character(dat$gene),
                         loc  = as.character(dat$ADlocus), stringsAsFactors = FALSE))
  d[!is.na(d$gene) & !is.na(d$loc), , drop = FALSE]
})

.sk_links <- function(e, scol, tcol) {
  key <- paste(e[[scol]], e[[tcol]], e$tierl, sep = "\r")
  sp <- split(seq_len(nrow(e)), key)
  do.call(rbind, lapply(names(sp), function(kk) {
    ix <- sp[[kk]]; pr <- strsplit(kk, "\r", fixed = TRUE)[[1]]
    gs <- unique(e$gene[ix])
    nl <- if (is.null(.g2l)) NA_integer_ else length(unique(.g2l$loc[.g2l$gene %in% gs]))
    data.frame(s = pr[1], t = pr[2], b = pr[3], Freq = length(ix),
               ng = length(gs), nl = nl, stringsAsFactors = FALSE)
  }))
}

home_river <- function(dark = FALSE) {
  e <- .sk_data
  if (is.null(e) || !nrow(e)) return(NULL)
  bg  <- "#ffffff"; ink <- "#233947"; sub <- "#5c6b75"
  ctxs <- unname(CTX_FULL[CTX_ORD]); ctxs <- ctxs[ctxs %in% e$ctxf]
  mods <- c("eQTL", "sQTL", "pQTL", "mQTL", "other assays")
  mods <- mods[mods %in% e$modg]
  tiers <- TIER_SEQ[TIER_SEQ %in% e$tierl]
  nodes <- c(ctxs, mods, tiers)
  idx <- stats::setNames(seq_along(nodes) - 1L, nodes)
  tot <- c(vapply(ctxs,  function(x) sum(e$ctxf == x),  1L),
           vapply(mods,  function(x) sum(e$modg == x),  1L),
           vapply(tiers, function(x) sum(e$tierl == x), 1L))
  lab <- nodes
  ncol <- c(unname(CTX_COL[CTX_ORD])[unname(CTX_FULL[CTX_ORD]) %in% ctxs],
            rep("#8fa0ac", length(mods)), unname(conf_pal[tiers]))
  L <- rbind(.sk_links(e, "ctxf", "modg"), .sk_links(e, "modg", "tierl"))
  L <- L[L$Freq > 0, , drop = FALSE]
  hex <- unname(conf_pal[L$b])
  rgba <- vapply(hex, function(h) {
    v <- grDevices::col2rgb(h)
    sprintf("rgba(%d,%d,%d,0.46)", v[1], v[2], v[3])
  }, character(1))
  llab <- sprintf("%s to %s · %s gene%s across %s %s in %s",
                  L$s, L$t,
                  format(L$ng, big.mark = ",", trim = TRUE),
                  ifelse(L$ng == 1, "", "s"),
                  format(L$nl, big.mark = ",", trim = TRUE),
                  ifelse(L$nl == 1, "locus", "loci"), L$b)
  ann <- function(x, txt, anc) list(x = x, y = -0.06, xref = "paper", yref = "paper",
    text = txt, showarrow = FALSE, font = list(size = 12, color = sub), xanchor = anc)
  plotly::plot_ly(type = "sankey", orientation = "h", arrangement = "snap",
    node = list(label = lab, pad = 16, thickness = 15, color = ncol,
                line = list(color = bg, width = 0.5),
                customdata = unname(format(tot, big.mark = ",", trim = TRUE)),
                hovertemplate = "%{label}<br>%{customdata} gene records<extra></extra>"),
    link = list(source = unname(idx[L$s]), target = unname(idx[L$t]),
                value = L$Freq, color = unname(rgba), label = llab,
                hovertemplate = "%{label}<extra></extra>")) %>%
    plotly::layout(font = list(family = "Inter, sans-serif", size = 12.5, color = ink),
      margin = list(l = 8, r = 8, t = 8, b = 30),
      paper_bgcolor = bg, plot_bgcolor = bg,
      annotations = list(ann(0, "Cell type", "left"), ann(0.5, "Assay", "center"),
                         ann(1, "Evidence tier", "right"))) %>%
    plotly::config(displayModeBar = FALSE)
}

# ---- Trans circos ----------------------------------------------------------
TMOD_COL <- c(snRNA = "#2a78d6", pQTL = "#da532c", gpQTL = "#17868f",
              Hotspot = "#a8791b", `assay not recorded` = "#95a3ad")

trans_loci_choices <- local({
  if (is.null(trans_pairs) || !nrow(trans_pairs)) return(character(0))
  t <- unique(trans_pairs[, c("ADlocus", "tgt")])
  v <- sort(table(t$ADlocus), decreasing = TRUE)
  stats::setNames(names(v), sprintf("%s  (%d distal gene%s)", names(v), as.integer(v), ifelse(as.integer(v) == 1, "", "s")))
})

trans_chord <- function(loc, top = 40) {
  lp <- trans_pairs[trans_pairs$ADlocus == loc, , drop = FALSE]
  if (!is.null(top) && nrow(lp)) {
    tn <- sort(table(as.character(lp$tgt)), decreasing = TRUE)
    if (length(tn) > top) {
      keep <- names(tn)[seq_len(top)]
      lp <- lp[match(as.character(lp$tgt), keep, nomatch = 0L) > 0L, , drop = FALSE]
    }
  }
  circos.clear()
  if (!nrow(lp)) { plot.new(); text(.5, .5, "No distal targets for this locus"); return(invisible()) }
  ch <- chr_offsets[match(CHR_ORD, as.character(chr_offsets$chr)), , drop = FALSE]
  src_ch <- as.character(lp$src_chr[1]); src_pos <- lp$src_pos[1]
  tg <- unique(lp[, c("tgt", "tgt_chr", "tgt_pos")])
  tg <- tg[!is.na(tg$tgt) & nzchar(as.character(tg$tgt)), , drop = FALSE]
  par(mar = c(0.4, 0.4, 1.4, 0.4))
  circos.par(gap.degree = 2.2, start.degree = 90, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.004, 0.004), points.overflow.warning = FALSE)
  circos.initialize(factors = factor(CHR_ORD, levels = CHR_ORD),
                    xlim = cbind(rep(0, nrow(ch)), ch$len))
  labs <- rbind(
    data.frame(s = as.character(tg$tgt_chr), x = tg$tgt_pos,
               l = as.character(tg$tgt), stringsAsFactors = FALSE),
    data.frame(s = src_ch, x = src_pos, l = as.character(lp$src[1]), stringsAsFactors = FALSE))
  labs <- labs[labs$s %in% CHR_ORD, , drop = FALSE]
  if (exists("circos.labels") && nrow(labs)) {
    circos.labels(sectors = labs$s, x = labs$x, labels = labs$l, cex = 0.52,
                  side = "outside", niceFacing = TRUE, col = "#3c4d59", line_col = "#c9d1d7")
  }
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.055,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = if (s == src_ch) "#f7e3da" else "#eceff2",
                  border = if (s == src_ch) "#da532c" else "#d3d9de", lwd = 0.7)
      circos.text(CELL_META$xcenter, 1.9, s, cex = 0.56, col = "#5c6b75",
                  facing = "inside", niceFacing = TRUE)
      if (s == src_ch) {
        hw <- 6e6
        circos.rect(max(0, src_pos - hw), 0, min(CELL_META$xlim[2], src_pos + hw), 1,
                    col = "#da532c", border = NA)
      }
    })
  for (i in seq_len(nrow(lp))) {
    if (!(as.character(lp$src_chr[i]) %in% CHR_ORD)) next
    if (!(as.character(lp$tgt_chr[i]) %in% CHR_ORD)) next
    md <- as.character(lp$modality[i])
    cl <- if (md %in% names(TMOD_COL)) TMOD_COL[[md]] else "#95a3ad"
    circos.link(as.character(lp$src_chr[i]), lp$src_pos[i],
                as.character(lp$tgt_chr[i]), lp$tgt_pos[i],
                col = grDevices::adjustcolor(cl, alpha.f = 0.5), lwd = 1.3)
  }
  title(main = sprintf("%s reaches %d gene%s on other chromosomes", loc, nrow(tg),
                       ifelse(nrow(tg) == 1, "", "s")),
        cex.main = 1.0, col.main = "#233947", font.main = 1, line = 0.1)
  present <- intersect(names(TMOD_COL), unique(as.character(lp$modality)))
  if (length(present))
    legend("bottomleft", legend = present, fill = unname(TMOD_COL[present]),
           border = NA, bty = "n", cex = 0.78, text.col = "#5c6b75")
  circos.clear()
}
