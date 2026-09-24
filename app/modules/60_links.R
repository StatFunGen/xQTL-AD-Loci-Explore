# ---- external variant resources ---------------------------------------
# variant_ID arrives as "chr1:17117119:G:A"; gnomAD and Open Targets want
# positional forms, so parse once and build each link from the parts.
variant_parts <- function(vid) {
  p <- strsplit(gsub("^chr", "", as.character(vid)), ":", fixed = TRUE)[[1]]
  if (length(p) < 4) return(NULL)
  list(chr = p[1], pos = p[2], ref = p[3], alt = p[4])
}

ext_links <- function(vid, rsid = NA) {
  v <- variant_parts(vid)
  out <- list()
  if (!is.null(v)) {
    out <- c(out, list(
      tags$a(href = sprintf("https://gnomad.broadinstitute.org/variant/%s-%s-%s-%s?dataset=gnomad_r4",
                            v$chr, v$pos, v$ref, v$alt),
             target = "_blank", rel = "noopener", class = "ext-link", "gnomAD"),
      tags$a(href = sprintf("https://genetics.opentargets.org/variant/%s_%s_%s_%s",
                            v$chr, v$pos, v$ref, v$alt),
             target = "_blank", rel = "noopener", class = "ext-link", "Open Targets"),
      tags$a(href = sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr%s:%s-%s",
                            v$chr, max(as.numeric(v$pos) - 500, 1), as.numeric(v$pos) + 500),
             target = "_blank", rel = "noopener", class = "ext-link", "UCSC")))
  }
  if (!is.na(rsid) && nzchar(rsid) && grepl("^rs", rsid)) {
    out <- c(out, list(
      tags$a(href = paste0("https://www.ncbi.nlm.nih.gov/snp/", rsid),
             target = "_blank", rel = "noopener", class = "ext-link", "dbSNP")))
  }
  if (!length(out)) return(NULL)
  div(class = "ext-links", tags$span(class = "ext-label", "Variant resources:"), out)
}

gene_links <- function(sym, ensg = NA) {
  out <- list()
  if (!is.na(sym) && nzchar(sym))
    out <- c(out, list(tags$a(href = paste0("https://gnomad.broadinstitute.org/gene/", sym, "?dataset=gnomad_r4"),
                              target = "_blank", rel = "noopener", class = "ext-link", "gnomAD gene")))
  if (!is.na(ensg) && nzchar(ensg))
    out <- c(out, list(tags$a(href = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", ensg),
                              target = "_blank", rel = "noopener", class = "ext-link", "Ensembl")))
  if (!length(out)) return(NULL)
  div(class = "ext-links", tags$span(class = "ext-label", "Gene resources:"), out)
}

# ---- build provenance ---------------------------------------------------
# Read at start-up from build_provenance.csv, which build_shiny_data.R writes
# beside the data. Nothing here is hardcoded, so the on-screen build stamp
# cannot silently go stale when the data is refreshed.
prov <- local({
  f <- "build_provenance.csv"
  if (!file.exists(f)) return(setNames(character(0), character(0)))
  x <- utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
  setNames(x$value, x$key)
})
pv <- function(k, default = "unknown") if (k %in% names(prov)) prov[[k]] else default
prov_chip <- sprintf("%s \u00b7 %s loci \u00b7 %s genes \u00b7 built %s",
                     pv("release"), pv("n_loci"), pv("n_genes"), pv("built_at"))
prov_tip  <- sprintf("tiers: %s | source table: %s", pv("tier_source"), pv("source_xlsx"))

# Light chart surfaces. Additive so the base g_theme definition is untouched.
g_theme <- g_theme + theme(
  plot.background  = element_rect(fill = PAPER, colour = NA),
  panel.background = element_rect(fill = PAPER, colour = NA),
  panel.grid.major = element_line(colour = LINE, linewidth = 0.3),
  panel.grid.minor = element_blank(),
  axis.text    = element_text(colour = MUTE),
  axis.title   = element_text(colour = MUTE),
  plot.title   = element_text(colour = INK),
  legend.text  = element_text(colour = MUTE),
  legend.title = element_text(colour = INK),
  legend.key            = element_blank(),
  legend.background     = element_rect(fill = PAPER, colour = NA),
  legend.box.background = element_rect(fill = PAPER, colour = NA),
  strip.background      = element_rect(fill = PAPER, colour = NA),
  strip.text            = element_text(colour = MUTE))

# chr2:127805598-128187321 -> chr2:127,805,598-128,187,321
pretty_region <- function(r) {
  m <- regmatches(r, regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", r))[[1]]
  if (length(m) != 4) return(r)
  a <- as.numeric(m[3]); b <- as.numeric(m[4])
  if (isTRUE(a == b)) return(sprintf("%s:%s", m[2], format(a, big.mark = ",")))
  sprintf("%s:%s\u2013%s", m[2], format(a, big.mark = ","), format(b, big.mark = ","))
}

# Compact keys shown next to the tables, so the tier scale and the evidence
# states can be read without opening the methods page.
TIER_DEFS <- c(
  T1 = "Localized xQTL support plus orthogonal gene-level evidence",
  T2 = "Multi-trait colocalization plus orthogonal evidence",
  T3 = "Localized support plus TWAS",
  T4 = "Localized support on its own",
  T5 = "Additional localized evidence at relaxed coverage",
  T6 = "TWAS or MR evidence only, with no localized AD-xQTL support")

tier_key <- function() {
  div(class = "tierkey",
    span(class = "keylab", "Evidence tier"),
    lapply(names(TIER_DEFS), function(k)
      span(class = "tk", title = TIER_DEFS[[k]],
           span(class = "tchip",
                style = sprintf("background:%s;color:%s", conf_pal[[k]], conf_txt[[k]]),
                k))))
}

ev_key <- function() {
  div(class = "evkey",
    span(class = "keylab", "Evidence status"),
    span(HTML('<span class="ev ev-yes">&#10003;</span>'), " detected"),
    span(HTML('<span class="ev ev-ns">NS</span>'), " tested, not significant"),
    span(HTML('<span class="ev ev-na">&mdash;</span>'), " not available"),
    span(class = "keynote",
         "A dash means the evidence is unknown for that entry, not that it was tested and found negative."))
}
hero_banner <- local({
  set.seed(7)
  circ <- function(cx, cy, R, id, sw1, sw2, nd) {
    ang <- seq(0, 2 * pi, length.out = 200)
    strand <- function(ph) {
      x <- seq(cx - R, cx + R, length.out = 90)
      y <- cy + R * 0.42 * sin((x - cx) / (R * 0.30) + ph)
      paste0("M", paste(sprintf("%.1f,%.1f", x, y), collapse = " L"))
    }
    xr <- seq(cx - R * 0.94, cx + R * 0.94, length.out = 18)
    y1 <- cy + R * 0.42 * sin((xr - cx) / (R * 0.30))
    y2 <- cy + R * 0.42 * sin((xr - cx) / (R * 0.30) + pi)
    rungs <- paste0(sprintf(
      '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#7ba3cc" stroke-width="1.1" opacity="0.7"/>',
      xr, y1, xr, y2), collapse = "")
    th <- runif(nd, 0, 2 * pi); rr <- R * runif(nd, 0.30, 1.02)
    nx <- cx + rr * cos(th); ny <- cy + rr * sin(th)
    den <- paste0(sprintf(
      '<path d="M%.1f,%.1f Q%.1f,%.1f %.1f,%.1f" fill="none" stroke="#5b87b5" stroke-width="%.1f"/>',
      cx, cy, cx + rr * 0.55 * cos(th + 0.35), cy + rr * 0.55 * sin(th + 0.35), nx, ny, sw2),
      collapse = "")
    nodes <- paste0(sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#5b87b5" opacity="0.75"/>',
                            nx, ny, runif(nd, 2.2, 4.6)), collapse = "")
    paste0(
      '<clipPath id="', id, '"><circle cx="', cx, '" cy="', cy, '" r="', R, '"/></clipPath>',
      '<circle cx="', cx, '" cy="', cy, '" r="', R,
      '" fill="none" stroke="#5b87b5" stroke-width="', sw1, '"/>',
      '<g clip-path="url(#', id, ')" opacity="0.7">',
      '<path d="', strand(0), '" fill="none" stroke="#7ba3cc" stroke-width="2"/>',
      '<path d="', strand(pi), '" fill="none" stroke="#7ba3cc" stroke-width="2"/>',
      rungs, '</g>',
      '<g clip-path="url(#', id, ')">', den, nodes,
      '<circle cx="', cx, '" cy="', cy, '" r="', round(R / 19), '" fill="#5b87b5"/></g>')
  }
  HTML(paste0(
    '<svg viewBox="0 0 1440 520" preserveAspectRatio="xMidYMid slice" aria-hidden="true" focusable="false">',
    '<g opacity="0.42">', circ(300, 250, 190, "em-a", 2.4, 1.5, 27), '</g>',
    '<g opacity="0.24">', circ(1185, 150, 108, "em-b", 2.0, 1.2, 18), '</g>',
    '</svg>'))
})
