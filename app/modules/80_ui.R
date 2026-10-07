# ---- shared filter sidebar (Explore) ---------------------------------------
explore_sidebar <- sidebar(
  width = 300, title = "Filters",
  textInput("q", info("Filter rows by gene or rsID", "Gene symbol, rsID, position (chr8:27369273), variant (chr8:27369273:A:C) or range (chr8:27300000-27400000)."), ""),
  selectizeInput("locus", "AD locus", choices = c("All" = "", loci), multiple = TRUE,
                 options = list(placeholder = "All loci")),
  checkboxGroupInput("sig", info("Significance", "Minimum GWAS p-value. Genome wide below 5e-8, suggestive below 1e-6, otherwise NS. Rows without a GWAS p-value form their own group."),
                     choices = c("genome wide","suggestive","ns", NO_P),
                     selected = c("genome wide","suggestive","ns", NO_P)),
  checkboxGroupInput("tier", info("Evidence tier", "Best evidence tier of the target gene, from T1 (strongest) to T6. T6 means TWAS or MR evidence without a fine-mapped AD-xQTL variant, and it does not exclude the gene. No tier assigned marks a gene in the locus without a tier; No gene mapped marks a variant without a target gene."),
                     choices = c(paste0("T", 1:6), NO_TIER, NO_GENE), selected = c(paste0("T", 1:6), NO_TIER, NO_GENE), inline = TRUE),
  checkboxGroupInput("ct", info("Context", "Target gene has xQTL evidence in the selected contexts."),
                     choices = setNames(ct_cols, ct_labels)),
  sliderInput("minlp", info("Min \u2212log10(p)", "Minimum GWAS significance of the variant."),
              min = 0, max = ceiling(max(dat$log10pval, na.rm = TRUE)), value = 0, step = 0.5),
  checkboxInput("transonly", "Trans evidence only", FALSE),
  hr(),
  downloadButton("dl", "Download filtered data", class = "btn-sm btn-outline-primary"),
  tags$button(type = "button", class = "btn btn-sm btn-outline-secondary link-btn",
    onclick = paste0("navigator.clipboard.writeText(window.location.href);",
                     "var t=this.innerHTML;this.innerHTML='Link copied';",
                     "var b=this;setTimeout(function(){b.innerHTML=t},1500);"),
    "Copy link to this view"),
  div(class = "note",
      "The address bar records these filters, so a bookmarked or shared link ",
      "reopens the same view.")
)

# ---- UI --------------------------------------------------------------------
COLTIPS <- c(
  Locus = "AD locus of this row. Labels are numbered within this release and can change between releases.",
  Gene = "Target gene of the variant at this locus. Opens Open Targets in a new tab.",
  rsID = "dbSNP identifier of the variant. Opens dbSNP in a new tab.",
  Pos = "Chromosome and position of the variant, GRCh38.",
  Sig = "GWAS significance of the variant. Genome wide below 5e-8, suggestive below 1e-6.",
  Tier = "Best evidence tier of the target gene across the whole release, so it can differ from the evidence in this row.",
  cV2F = "Consensus variant-to-function score. Higher values mean the variant is more likely to affect function.",
  TWAS = "Transcriptome-wide association result for the gene, gene-level and not tied to a variant. A filled mark means tested and significant, an open mark tested and not significant. An empty cell means no result is recorded in this release, which is different from a zero.",
  MR = "Mendelian randomization result for the gene. A filled mark means tested and significant, an open mark tested and not significant. An empty cell means no result is recorded in this release, which is different from a zero.",
  cTWAS = "Causal TWAS result for the gene. A filled mark means tested and significant, an open mark tested and not significant. An empty cell means no result is recorded in this release, which is different from a zero.",
  Trans = "Whether the variant has trans-xQTL evidence for genes outside this locus.",
  `Cell types` = "Contexts (single-nucleus cell types or bulk tissue) with xQTL evidence for this variant and gene. An empty cell means no result is recorded for that context, which is different from a tested null.",
  `Max xQTL PIP` = "Largest xQTL fine-mapping PIP across the variants recorded for this gene at this locus.",
  `GWAS PIP` = "Largest posterior inclusion probability (PIP) of the variant in an AD GWAS credible set. It is on a different scale from the xQTL PIP.",
  `xQTL PIP` = "Largest xQTL fine-mapping PIP of the variant across methods and contexts.",
  `TWAS z` = "Largest absolute TWAS z-score for the gene.",
  `#Ctx` = "Number of contexts with xQTL evidence for this variant and gene.",
  `-log10p` = "GWAS minimum -log10 p-value for the variant.",
  `Ordered contexts` = "Contexts ranked by strength of evidence, strongest first."
)

js_gene_render <- DT::JS(paste0(
  "function(data, type, row){ if (type !== 'display' || !data) return data; ",
  "return '<a class=ext-link target=_blank rel=noopener href=https://platform.opentargets.org/search?q=' ",
  "+ encodeURIComponent(data) + '>' + data + '</a>'; }"))

js_rsid_render <- DT::JS(paste0(
  "function(data, type, row){ if (type !== 'display' || !data) return data; ",
  "if (!/^rs[0-9]+$/.test(data)) return data; ",
  "return '<a class=ext-link target=_blank rel=noopener href=https://www.ncbi.nlm.nih.gov/snp/' ",
  "+ data + '>' + data + '</a>'; }"))

hdr_js <- DT::JS(paste0(
  "function(thead){ var t = ", jsonlite::toJSON(as.list(COLTIPS), auto_unbox = TRUE), "; ",
  "$(thead).find('th').each(function(){ var k = $(this).text().trim(); ",
  "if (t[k]) { $(this).attr('title', t[k]); } }); }"))

link_defs <- function(gene_idx = NULL, rsid_idx = NULL) {
  out <- list()
  if (!is.null(gene_idx) && length(gene_idx))
    out <- c(out, list(list(targets = gene_idx, render = js_gene_render)))
  if (!is.null(rsid_idx) && length(rsid_idx))
    out <- c(out, list(list(targets = rsid_idx, render = js_rsid_render)))
  out
}

igv_head <- tagList(
  tags$script(src = "https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js", defer = NA),
  ## Tell the wanggroup.org page the first view has drawn, so it can drop its loading screen.
  tags$script(HTML("$(document).one('shiny:idle', function(){ try { if (window.parent !== window) window.parent.postMessage('adx-ready', 'https://wanggroup.org'); } catch (e) {} }); (function(){ if (window.parent === window) return; function send(){ try { window.parent.postMessage({ adxq: location.search }, 'https://wanggroup.org'); } catch (e) {} } ['pushState','replaceState'].forEach(function(k){ var o = history[k]; history[k] = function(){ var r = o.apply(this, arguments); send(); return r; }; }); })();")),
  tags$script(HTML(paste0(
    "var adIgvBrowser = null; var adIgvPending = null; var adIgvTries = 0;",
    "function adIgvGo(loc){",
    "  if (!loc) return; adIgvPending = loc;",
    "  var el = document.getElementById('igv_div'); if (!el) return;",
    "  if (typeof igv === 'undefined' || !el.offsetWidth) {",
    "    if (adIgvTries++ > 200) return;",
    "    setTimeout(function(){ adIgvGo(adIgvPending); }, 300); return;",
    "  }",
    "  adIgvTries = 0;",
    "  if (adIgvBrowser) {",
    "    try { adIgvBrowser.search(loc); adIgvBrowser.resize(); } catch(e) {}",
    "    return;",
    "  }",
    "  el.innerHTML = '';",
    "  igv.createBrowser(el, { genome: 'hg38', locus: loc }).then(function(b){",
    "    adIgvBrowser = b;",
    "    setTimeout(function(){ try { b.resize(); } catch(e) {} }, 250);",
    "  }).catch(function(){",
    "    el.innerHTML = '<div class=\\\"note\\\">The genome browser could not be loaded. The rest of the locus evidence on this page is unaffected.</div>';",
    "  });",
    "}",
    "window.addEventListener('resize', function(){",
    "  if (adIgvBrowser) { try { adIgvBrowser.resize(); } catch(e) {} }",
    "});",
    "$(document).on('shiny:connected', function(){",
    "  Shiny.addCustomMessageHandler('igv_locus', function(msg){ adIgvGo(msg.locus); });",
    "});")))
)

extra_css <- tags$style(HTML('
.statline{display:flex;flex-wrap:wrap;gap:26px;align-items:baseline;margin:2px 0 22px 0;color:#0b0b0b;font-size:15px}
.statline .sk{color:#52514e;font-size:11px;letter-spacing:.08em;text-transform:uppercase;margin-right:8px}
.statline .shiny-text-output{display:inline;font-variant-numeric:tabular-nums}
.hero{max-width:740px;margin:60px auto 0 auto;padding:0 16px;text-align:center}
.hero-title{font-size:34px;font-weight:600;letter-spacing:-.01em;color:#0b0b0b;margin-bottom:10px}
.hero-sub{font-size:16px;color:#52514e;margin-bottom:30px;line-height:1.55}
.hero-search{margin:0 auto 22px auto;max-width:560px;text-align:left}
.hero-search .selectize-input{padding:11px 14px;font-size:16px}
.hero-stats{display:flex;justify-content:center;gap:30px;flex-wrap:wrap;color:#52514e;font-size:14px;margin-bottom:26px}
.hero-stats strong{color:#0b0b0b;font-weight:600;font-variant-numeric:tabular-nums}
.hero-links{display:flex;justify-content:center;gap:28px;font-size:14px}
.hero-links a{color:#2a78d6;text-decoration:none}
.hero-links a:hover{text-decoration:underline}
.sheet{max-width:1080px;margin:54px auto 0 auto;padding:0 16px 40px 16px}
.sheet-head{display:flex;align-items:baseline;justify-content:space-between;gap:24px;flex-wrap:wrap;border-top:1px solid #e1e0d9;padding-top:22px;margin-bottom:28px}
.sheet-title{font-size:13px;letter-spacing:.08em;text-transform:uppercase;color:#52514e}
.sheet .statline{margin:0}
.panel{margin-bottom:34px}
.panel-title{font-size:13px;color:#52514e;margin-bottom:12px}
.locus-bar{display:flex;gap:44px;align-items:flex-start;flex-wrap:wrap;border-top:1px solid #e1e0d9;padding-top:22px;margin-bottom:32px}
.locus-pick{min-width:240px}
.locus-bar .locus-record{flex:1 1 380px}
.tierkey,.evkey{display:flex;gap:14px;flex-wrap:wrap;align-items:center;font-size:12px;color:#52514e;margin:0 0 12px 0}
.keylab{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:#898781;margin-right:4px}
.tierkey .tk{display:inline-flex;align-items:center}
.tchip{display:inline-block;min-width:24px;text-align:center;border-radius:3px;padding:2px 6px;font-size:11px;font-weight:600}
.ev{font-size:12px;font-weight:600}
.ev-yes{color:#0ca30c}
.ev-ns{color:#898781;font-weight:500}
.ev-na{color:#b9b7b0}
.keynote{color:#898781}
.viewbar{margin-bottom:12px}
.viewbar .radio-inline{font-size:13px;margin-right:18px}
.doc{max-width:1160px;margin:34px auto 64px auto;padding:0 16px;display:flex;gap:56px;align-items:flex-start}
.doc-toc{position:sticky;top:70px;flex:0 0 190px;font-size:13px;line-height:1.95}
.doc-toc a{color:#52514e;text-decoration:none;display:block}
.doc-toc a:hover{color:#2a78d6;text-decoration:underline}
.doc-toc .keylab{margin-bottom:8px}
.doc-body{flex:1 1 auto;min-width:0;max-width:760px}
.doc-sec{border-top:1px solid #e1e0d9;padding:30px 0 6px 0;background:transparent;border-radius:0;box-shadow:none}
.doc-sec:first-child{border-top:none;padding-top:0}
.doc-h{font-size:20px;font-weight:600;color:#0b0b0b;margin:0 0 14px 0}
.doc-body p{line-height:1.68;margin:0 0 15px 0}
.doc-body .note{background:transparent;border-left:2px solid #e1e0d9;padding:2px 0 2px 14px;color:#52514e}
.doc-body table{width:100%;font-size:13px;margin:6px 0 18px 0;border-collapse:collapse}
.doc-body th{text-align:left;font-weight:600;color:#52514e;font-size:11px;letter-spacing:.06em;text-transform:uppercase;border-bottom:1px solid #e1e0d9;padding:6px 12px 6px 0}
.doc-body td{border-bottom:1px solid #f0efec;padding:7px 12px 7px 0;vertical-align:top}
.doc-body .accordion{border:none}
@media (max-width: 900px){.doc{flex-direction:column;gap:20px}.doc-toc{position:static;flex:none}}
.start{max-width:900px;margin:34px auto 60px auto;padding:26px 16px 0 16px;border-top:1px solid #e1e0d9}
.start-grid{display:flex;gap:34px;flex-wrap:wrap;margin:14px 0 20px 0}
.start-step{flex:1 1 230px;display:flex;gap:12px;align-items:flex-start}
.start-n{flex:0 0 22px;height:22px;border-radius:50%;background:#f0efec;color:#52514e;font-size:12px;font-weight:600;display:flex;align-items:center;justify-content:center;margin-top:2px}
.start-step b{font-size:14px;color:#0b0b0b}
.start-step p{font-size:13px;line-height:1.6;color:#52514e;margin:4px 0 0 0}
.start-ex{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin-bottom:18px}
.start-ex a{font-size:13px;color:#2a78d6;text-decoration:none;border:1px solid #e1e0d9;border-radius:3px;padding:3px 10px}
.start-ex a:hover{border-color:#2a78d6}
.start-note{font-size:12px;color:#898781;line-height:1.6;max-width:70ch}
.tn-row{border-top:1px solid #f0efec;padding:12px 0}
.tn-row:first-child{border-top:none}
.tn-src{font-size:15px;font-weight:600;color:#0b0b0b;margin-bottom:5px}
.tn-count{font-weight:400;font-size:12px;color:#898781;margin-left:10px}
.tn-mods{margin-bottom:7px}
.tn-mod{display:inline-block;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:#52514e;background:#f0efec;border-radius:3px;padding:2px 7px;margin-right:6px}
.tn-tg{display:flex;flex-wrap:wrap;gap:6px}
.tn-chip{font-size:12px;color:#2a78d6;border:1px solid #e1e0d9;border-radius:3px;padding:2px 7px}
.brand-logo{height:34px;width:auto;margin-right:12px;vertical-align:middle}
.orgs{max-width:900px;margin:0 auto 64px auto;padding:24px 16px 0 16px;border-top:1px solid #e1e0d9;display:flex;gap:44px;align-items:center;flex-wrap:wrap}
.orgs img{height:34px;width:auto}
.orgs .keylab{flex:1 1 100%;margin-bottom:4px}
.panel-note{font-size:12.5px;line-height:1.6;color:#898781;max-width:70ch;margin:-4px 0 12px 0}
.brand-b,.hero-title,.locus-title,.doc-h{color:#0d1b2a}
.navbar{border-bottom:2px solid #0d1b2a}
.site-footer{margin-top:56px;padding:34px 16px 44px 16px;background:#fbfbf9}
.foot-in{max-width:1080px;margin:0 auto;display:flex;gap:48px;flex-wrap:wrap}
.foot-col{flex:1 1 210px;min-width:190px}
.foot-h{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:#898781;margin-bottom:9px}
.foot-t{font-size:12.5px;line-height:1.65;color:#52514e;margin-bottom:3px}
.foot-col a{display:block;font-size:12.5px;color:#2a78d6;text-decoration:none;margin-bottom:3px}
.foot-col a:hover{text-decoration:underline}
.foot-logos{max-width:1080px;margin:28px auto 0 auto;display:flex;gap:34px;align-items:center;flex-wrap:wrap}
.foot-logos img{height:30px;width:auto}
.hero-ex{display:flex;gap:10px;align-items:center;flex-wrap:wrap;justify-content:center;margin:-6px 0 24px 0;font-size:13px}
.hero-ex-lab{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:#898781;margin-right:2px}
.hero-ex a{color:#2a78d6;text-decoration:none;border-bottom:1px solid #cfe0f5}
.hero-ex a:hover{border-bottom-color:#2a78d6}
.hero-ex-k{color:#898781;font-size:11.5px;margin-right:12px}
.tl-head{display:flex;gap:18px;align-items:baseline;flex-wrap:wrap;margin-bottom:12px}
.tl-h-src{font-size:15px;font-weight:600;color:#233947}
.tl-h-n{font-size:12px;color:#898781}
.tl-list{border-top:1px solid #f0efec}
.tl-row{display:flex;gap:16px;align-items:baseline;flex-wrap:wrap;border-bottom:1px solid #f0efec;padding:8px 0}
.tl-coord{flex:0 0 150px;font-size:13px;color:#233947;font-variant-numeric:tabular-nums}
.tl-gene{flex:0 0 130px;font-size:13px;font-weight:600;color:#2a78d6}
.tl-mod{flex:0 0 110px;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:#52514e}
.tl-src{font-size:12px;color:#898781}
.hero{position:relative}
.hero > *{position:relative;z-index:1}
.hero-banner{position:absolute;left:0;right:0;top:-18px;height:170px;overflow:hidden;pointer-events:none;z-index:0}
.hero-banner svg{width:100%;height:170px;display:block}
.bn-helix{fill:none;stroke:#233947;stroke-width:1.1;opacity:.12}
.bn-curve{fill:none;stroke:#2a78d6;stroke-width:1.1;opacity:.14;stroke-dasharray:4 5}
.bn-dots circle{fill:#2a78d6;opacity:.15}
.hero-wrap{position:relative;overflow:hidden;padding-bottom:28px;background:linear-gradient(180deg,#eef4fb 0%,#f7fafd 52%,#ffffff 100%)}
.hero-wrap > .hero,.hero-wrap > .start{position:relative;z-index:1}
.hero-bg{position:absolute;left:0;right:0;bottom:0;height:124px;z-index:0;pointer-events:none}
.hero-bg svg{width:100%;height:124px;display:block}
.kary rect{fill:#233947;opacity:.055}
.kary rect.kb-hi{fill:#2a78d6;opacity:.13}
.kary rect.kb-cen{fill:#ffffff;opacity:.6}
.hero{margin-top:0;padding-top:64px}
.start{border-top:1px solid rgba(35,57,71,.10)}
.ti-head{font-size:13px;color:#4a5b66;margin:0 0 10px 2px}
.ti-wrap{display:flex;flex-direction:column;gap:12px}
.ti-card{border:1px solid #e6e4e0;border-radius:10px;background:#fff;padding:8px 12px 4px}
.ti{width:100%;height:auto;display:block;max-height:180px}
.ti-bar{fill:#eceef0;stroke:#dcdfe2;stroke-width:.8}
.ti-hit{fill:#233947}
.ti-src{fill:#da532c}
.ti-dot{stroke:#fff;stroke-width:1.4}
.ti-link{fill:none;stroke-width:2;stroke-linecap:round;opacity:.8}
.ti-lab{font-size:12.5px;font-weight:600;fill:#233947}
.ti-chr{font-size:11.5px;fill:#6b7a85}
.ti-tag{font-size:11px;font-weight:500;letter-spacing:.2px;paint-order:stroke;stroke:#ffffff;stroke-width:3.5px;stroke-linejoin:round}
.ti-gene{font-size:12.5px;font-weight:600}
.ti-sub{font-size:11.5px;color:#7b8790;padding:0 2px 6px;text-align:center}
.rv-bar{display:flex;gap:16px;align-items:flex-end;flex-wrap:wrap;margin:10px 0 2px;padding:12px 14px;background:#f7f8fa;border:1px solid #e9e7e3;border-radius:10px}
.rv-f{flex:1 1 190px;min-width:170px}
.rv-f .form-group{margin-bottom:0}
.rv-f label{font-size:11.5px;letter-spacing:.3px;text-transform:uppercase;color:#6b7a85;font-weight:600;margin-bottom:4px}
.rv-b{flex:0 0 auto;min-width:0}
.rv-b .btn{margin-bottom:1px}
.rv-stat{font-size:12.5px;color:#4a5b66;margin:10px 2px 2px}
.rv-dl{display:flex;gap:18px;flex-wrap:wrap;font-size:12.5px;margin:6px 2px 0}
.rv-dl a{color:#2a78d6;text-decoration:none}
.rv-dl a:hover{text-decoration:underline}
/* ---- design canvas language ---- */
body,.navbar,.nav-link,h1,h2,h3{font-family:Inter,system-ui,-apple-system,sans-serif}
.dc-mono,.dc-chip,.dc-count-n,.dc-num,code,kbd{font-family:"IBM Plex Mono",ui-monospace,monospace}
.card,.btn,.form-control,.form-select,.selectize-input,.selectize-dropdown,.panel,.sheet,.badge,.nav-link,.dataTables_wrapper .dataTable,.shiny-input-container input{border-radius:0 !important}
a{color:var(--text);text-decoration:none}
a:hover{color:var(--cta)}
.navbar{border-bottom:1px solid var(--line);min-height:62px}
.navbar .nav-link{font-size:13.5px;font-weight:500;letter-spacing:.01em}
.navbar .nav-link.active{box-shadow:inset 0 -2px 0 var(--text);color:var(--text)}
.dc-eyebrow{font-size:11.5px;letter-spacing:.1em;text-transform:uppercase;color:var(--muted)}
.dc-kicker{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted)}
.dc-hero{position:relative;background-color:#fff;background-image:url(hero-bg.jpg);background-repeat:no-repeat;background-position:right center;background-size:cover;overflow:hidden;padding:74px 40px 60px;border-bottom:1px solid var(--line)}@media (prefers-contrast:more),(forced-colors:active){.dc-hero{background-image:none}}
.dc-hero-bg,.dc-hero-wash{position:absolute;left:0;right:0;top:0;bottom:0}
.dc-hero-bg svg{width:100%;height:100%;display:block}
.dc-hero-wash{background:radial-gradient(ellipse 58% 70% at 44% 48%,rgba(255,255,255,.90) 0%,rgba(255,255,255,.58) 50%,rgba(255,255,255,.06) 100%)}
.dc-hero-in{position:relative;max-width:940px;margin:0 auto;text-align:center}
.dc-h1{margin:14px 0 0;font-size:44px;line-height:1.12;font-weight:600;letter-spacing:-.025em;color:var(--text)}
.dc-lede{margin:16px auto 0;max-width:680px;font-size:15px;line-height:1.6;color:var(--muted)}
.dc-searchwrap{display:flex;margin:34px auto 0;max-width:760px;border:1px solid var(--text);background:#fff;align-items:stretch}
.dc-searchwrap .form-group,.dc-searchwrap .shiny-input-container{margin:0;flex:1 1 auto;min-width:0}
.dc-searchwrap .selectize-input{border:0 !important;box-shadow:none !important;min-height:56px;padding:16px 18px;font-size:16px;background:#fff}
.dc-searchwrap .selectize-input input{font-size:16px}
.dc-searchbtn{display:flex;align-items:center;padding:0 30px;background:var(--text);color:#fff;font-size:15px;font-weight:500;white-space:nowrap}
.dc-try{display:flex;justify-content:center;align-items:center;gap:10px;margin-top:16px;font-size:13px;color:var(--muted);flex-wrap:wrap}
.dc-chip{border:1px solid var(--line);padding:5px 11px;background:#fff;font-size:12.5px;color:var(--text)}
.dc-chip:hover{border-color:var(--text)}
.dc-counts{display:flex;flex-wrap:nowrap;align-items:flex-start;justify-content:center;gap:clamp(16px,2.6vw,42px);margin-top:44px;padding-top:6px}
.dc-counts>div{flex:0 1 auto;min-width:0;text-align:center}
.dc-count-n{font-size:clamp(21px,2.3vw,34px);font-weight:600;letter-spacing:-.02em;line-height:1;white-space:nowrap}
.dc-count-l{font-size:12.5px;color:var(--muted);margin-top:4px}
.dc-band{background:#f4f5f6;border-bottom:1px solid var(--line);padding:36px 40px 40px}
.dc-sec{background:#fff;border-bottom:1px solid var(--line);padding:32px 40px 34px}
.dc-sec:last-child{border-bottom:0}
.dc-sechead{display:flex;align-items:flex-end;justify-content:space-between;border-bottom:1px solid var(--line);padding-bottom:12px;gap:24px;flex-wrap:wrap}
.dc-h2{margin:0;font-size:18px;font-weight:600;letter-spacing:-.015em}
.dc-sub{font-size:13px;color:var(--muted);margin-top:5px;line-height:1.55}
.dc-grid3{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:1px;background:var(--line);border:1px solid var(--line);margin-top:20px}
.dc-cell{background:#fff;padding:24px 24px 22px;display:flex;flex-direction:column;gap:9px;min-height:132px}
.dc-cell-t{font-size:16px;font-weight:600;letter-spacing:-.01em}
.dc-cell-b{font-size:13px;line-height:1.6;color:var(--muted)}
.dc-cell-c{margin-top:auto;font-size:13px;color:var(--cta)}
.dc-btn{font-size:12.5px;border:1px solid var(--line);padding:6px 12px;background:#fff;color:var(--text);display:inline-block}
.dc-btn:hover{border-color:var(--text);color:var(--text)}
.dc-exp{display:flex;align-items:center;gap:8px}
.dc-figt{font-size:15px;font-weight:600}
.dc-figc{font-size:12.5px;color:var(--muted);margin-top:4px;line-height:1.55}
.dc-note{font-size:12.5px;color:var(--muted);line-height:1.6}
.dc-rule{border-top:1px solid var(--rule)}
.dc-tierdefs{display:grid;grid-template-columns:repeat(3,1fr);gap:10px 32px;margin-top:16px;padding-top:14px;border-top:1px solid var(--rule);font-size:12.5px;line-height:1.55;color:var(--muted)}
.dc-tierdefs b{color:var(--text);font-weight:600}
.dc-legend{display:flex;flex-direction:column;gap:7px}
.dc-legend span.k{display:inline-flex;align-items:center;gap:9px;font-size:12.5px}
.dc-dot{width:10px;height:10px;display:inline-block;flex:0 0 auto}
.dc-howto{margin-top:16px;padding-top:14px;border-top:1px solid var(--rule)}
.dc-howto-g{display:grid;grid-template-columns:auto 1fr;gap:6px 10px;font-size:12.5px;line-height:1.5;color:var(--muted);margin-top:9px;max-width:760px}
.dc-howto-g .n{color:var(--text);font-family:"IBM Plex Mono",ui-monospace,monospace}
.navbar .nav-link.active{box-shadow:none;border-bottom:2px solid var(--cta);font-weight:600}
.rv-bar{background:#fff;border:0;border-top:1px solid var(--rule);border-bottom:1px solid var(--rule);padding:14px 0 12px;margin:0 0 4px;border-radius:0}
.rv-stat{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12.5px}
.doc-wrap{display:flex;align-items:stretch;gap:0;background:#fff;border:1px solid var(--line)}
.doc-side{flex:0 0 260px;border-right:1px solid var(--line);padding:26px 22px 34px;background:#fff}
.doc-nav{display:flex;flex-direction:column;gap:9px;margin-top:12px;font-size:12.5px}
.doc-nav a{color:var(--muted)}
.doc-nav a.on{color:var(--text);font-weight:600}
.doc-nav a:hover{color:var(--cta)}
.doc-side-note{margin-top:22px;padding-top:16px;border-top:1px solid var(--rule);font-size:12px;line-height:1.6;color:var(--muted)}
.doc-main{flex:1 1 auto;min-width:0;padding:30px 44px 46px}
.doc-h1{margin:10px 0 0;font-size:22px;font-weight:600;letter-spacing:-.015em;scroll-margin-top:80px}
.doc-h2{margin:34px 0 0;padding-top:22px;border-top:1px solid var(--line);font-size:18px;font-weight:600;letter-spacing:-.015em;scroll-margin-top:80px}
.doc-p{margin:12px 0 0;font-size:14px;line-height:1.7;max-width:760px;color:var(--text)}
.doc-tiers{display:flex;flex-direction:column;gap:1px;background:var(--line);border:1px solid var(--line);margin-top:22px;max-width:900px}
.doc-tier{background:#fff;padding:14px 18px;display:flex;gap:14px;align-items:flex-start}
.doc-tier .sw{width:13px;height:13px;flex:0 0 auto;margin-top:2px;display:inline-block}
.doc-tier .tk{flex:0 0 34px;font-weight:600;font-size:13.5px}
.doc-tier .df{flex:1 1 auto;font-size:13px;line-height:1.6}
.doc-steps{display:grid;grid-template-columns:auto 1fr;gap:10px 12px;font-size:13.5px;line-height:1.6;max-width:720px;margin-top:16px}
.doc-steps .n{color:var(--cta);font-weight:600}
.doc-callout{margin-top:16px;border:1px solid var(--cta);padding:18px 20px;display:flex;gap:16px;align-items:flex-start;max-width:860px}
.doc-callout .bar{width:10px;height:10px;background:var(--cta);flex:0 0 auto;margin-top:6px;display:inline-block}
.doc-callout-t{font-size:14px;font-weight:600}
.doc-tbl{width:100%;border-collapse:collapse;font-size:13px;margin-top:16px;max-width:900px}
.doc-tbl th{text-align:center;padding:10px 8px;font-weight:500;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);border-bottom:1px solid var(--line)}
.doc-tbl th:first-child{text-align:left;padding-left:0}
.doc-tbl td{padding:7px 8px;border-bottom:1px solid var(--rule);text-align:center;font-family:"IBM Plex Mono",ui-monospace,monospace}
.doc-tbl td:first-child{text-align:left;padding-left:0;font-family:Inter,system-ui,sans-serif}
.doc-tbl td.z{color:var(--faint)}
.dc-sechead{flex-wrap:nowrap}
.dc-sechead > *{min-width:0}
.dc-exp{flex:0 0 auto;white-space:nowrap}
@media (max-width:900px){.dc-sechead{flex-wrap:wrap}}
.doc{display:flex;align-items:flex-start;gap:0;background:#fff;border:1px solid var(--line)}
.doc-toc{flex:0 0 268px;border-right:1px solid var(--line);padding:28px 22px 34px;position:sticky;top:70px;display:flex;flex-direction:column;gap:9px;font-size:12.5px}
.doc-toc a{color:var(--muted);line-height:1.45}
.doc-toc a:hover{color:var(--cta)}
.doc-toc > b,.doc-toc .keylab{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);font-weight:600;margin-bottom:4px}
.doc-body{flex:1 1 auto;min-width:0;padding:30px 44px 46px;background:#fff}
.doc-sec{margin-top:34px;padding-top:22px;border-top:1px solid var(--line);scroll-margin-top:80px}
.doc-sec:first-child{margin-top:0;padding-top:0;border-top:0}
.doc-h{margin:0 0 10px;font-size:18px;font-weight:600;letter-spacing:-.015em;color:var(--text)}
.doc-sec:first-child .doc-h{font-size:22px}
.doc-body p{margin:12px 0 0;font-size:14px;line-height:1.7;max-width:780px;color:var(--text)}
.doc-body p.note{font-size:12.5px;color:var(--muted);line-height:1.6}
.doc-body a{color:var(--text);border-bottom:1px solid var(--line)}
.doc-body a:hover{color:var(--cta)}
.dl-top{background:#fff;border:1px solid var(--line);border-bottom:0;padding:26px 40px 24px}
.dl-top .dc-sub{margin-top:7px}
.dl-wrap{display:flex;gap:24px;align-items:flex-start;margin-top:24px}
.dl-main{flex:1 1 auto;min-width:0;background:#fff;border:1px solid var(--line);padding:22px 26px 26px}
.dl-side{flex:0 0 330px;display:flex;flex-direction:column;gap:20px}
.dl-card{background:#fff;border:1px solid var(--line);padding:20px 22px 22px}
.dl-card-t{font-size:14px;font-weight:600;margin-bottom:10px}
.dl-card-b{font-size:12.5px;line-height:1.6;color:var(--muted)}
.dl-code{margin:12px 0 0;background:#f4f5f6;border:1px solid var(--rule);padding:11px 12px;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:11.5px;line-height:1.6;color:var(--text);white-space:pre-wrap;overflow-x:auto;border-radius:0}
.dl-tiers{display:grid;grid-template-columns:auto 1fr;gap:8px 12px;font-size:12px;line-height:1.5;color:var(--muted)}
.dl-tiers .k{font-weight:600}
.dl-cite{margin-top:12px;font-size:12px;line-height:1.6;border-left:2px solid var(--cta);padding-left:11px;color:var(--text)}
.dl-grp{margin-top:22px}
.dl-grp-h{display:flex;align-items:baseline;gap:12px;border-bottom:1px solid var(--rule);padding-bottom:8px}
.dl-grp-h .t{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--text);font-weight:600}
.dl-grp-h .n{font-size:12px;color:var(--muted)}
.dl-tbl{width:100%;border-collapse:collapse;font-size:13px}
.dl-tbl td{padding:12px 8px;vertical-align:top;border-bottom:1px solid var(--rule)}
.dl-tbl td.f{padding-left:0}
.dl-tbl .fn{font-weight:500}
.dl-tbl .de{color:var(--muted);font-size:12px;margin-top:3px;line-height:1.5}
.dl-tbl .cl{color:var(--faint);font-size:11.5px;margin-top:5px;line-height:1.5}
.dl-tbl td.n{text-align:right;color:var(--muted);width:112px;white-space:nowrap}
.dl-tbl td.fm{color:var(--muted);width:90px}
.dl-tbl td.a{text-align:right;width:160px;padding-right:0;white-space:nowrap}
.dl-tbl tr.ext .fn{color:var(--muted)}
.dl-foot{margin-top:20px;padding-top:14px;border-top:1px solid var(--line);display:flex;align-items:center;gap:16px;font-size:12.5px;color:var(--muted);flex-wrap:wrap}
.dl-foot a{border-bottom:1px solid var(--line)}
.dl-foot a:nth-of-type(1){margin-left:auto}
@media (max-width:1100px){.dl-wrap{flex-wrap:wrap}.dl-side{flex:1 1 100%}}
.tc-top{background:#fff;border:1px solid var(--line);border-bottom:0;padding:26px 40px 20px}
.tc-top .dc-sub{margin-top:7px}
.tc-bar{display:flex;gap:26px;align-items:flex-end;flex-wrap:wrap;margin-top:18px;padding-top:16px;border-top:1px solid var(--rule)}
.tc-f{flex:0 0 260px}
.tc-f2{flex:1 1 320px}
.tc-f3{flex:0 0 auto;margin-left:auto}
.tc-f .form-group{margin-bottom:0}
.tc-f label,.tc-f .control-label{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);font-weight:600;margin-bottom:6px}
.tc-f .checkbox-inline{font-size:12.5px;margin-right:14px}
.tc-wrap{background:#f4f5f6;border:1px solid var(--line);border-top:0;padding:24px 40px 32px}
.tc-card{background:#fff;border:1px solid var(--line);padding:22px 26px 20px}
.tc-head{display:flex;align-items:flex-start;justify-content:space-between;border-bottom:1px solid var(--rule);padding-bottom:12px;margin-bottom:8px}
.tc-title{font-size:16px;font-weight:600;letter-spacing:-.01em;font-family:"IBM Plex Mono",ui-monospace,monospace}
.tc-sub{font-size:12.5px;color:var(--muted);margin-top:4px}
.tc-svg{width:100%;height:auto;display:block;overflow:visible}
.tc-bar-r{fill:#eceef0}
.tc-bar rect{}
.tc-svg .tc-bar{fill:#e6eaed;stroke:#d3d9de;stroke-width:.8}
.tc-svg .tc-cen{fill:#ffffff;opacity:.95}
.tc-svg .tc-mark{fill:#da532c}
.tc-svg .tc-src{fill:#da532c}
.tc-svg .tc-link{fill:none;stroke-width:1.3;stroke-dasharray:5 4;opacity:.8}
.tc-svg .tc-chr{font-size:12.5px;font-weight:500;fill:#233947;font-family:"IBM Plex Mono",ui-monospace,monospace}
.tc-svg .tc-role{font-size:11px;fill:#5c6b75}
.tc-svg .tc-lab{font-size:12px;font-family:"IBM Plex Mono",ui-monospace,monospace}
.tc-svg .tc-srclab{fill:#da532c}
.tc-svg .tc-gene{font-size:12.5px;font-weight:600}
.tc-svg .tc-pos{font-size:11px;fill:#5c6b75;font-family:"IBM Plex Mono",ui-monospace,monospace}
.tc-legend{display:flex;gap:22px;flex-wrap:wrap;margin-top:18px;padding-top:14px;border-top:1px solid var(--rule);font-size:12px;color:var(--muted);align-items:center}
.tc-legend .k{display:inline-flex;align-items:center;gap:7px}
.tc-legend .sw{width:4px;height:14px;display:inline-block}
.tc-legend .dash{width:22px;height:0;border-top:1.3px dashed #da532c;display:inline-block}
.lt-top{background:#fff;border:1px solid var(--line);border-bottom:0;padding:26px 40px 20px}
.lt-top .dc-sub{margin-top:7px}
.lt-bar{display:flex;gap:20px;align-items:flex-end;flex-wrap:wrap;margin-top:18px;padding-top:16px;border-top:1px solid var(--rule)}
.lt-f{flex:0 0 210px}
.lt-f .form-group{margin-bottom:0}
.lt-f label,.lt-f .control-label{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);font-weight:600;margin-bottom:6px}
.lt-stat{flex:1 1 auto;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12.5px;color:var(--muted);padding-bottom:8px}
.lt-e{flex:0 0 auto;margin-left:auto}
.lt-wrap{background:#f4f5f6;border:1px solid var(--line);border-top:0;padding:24px 40px 32px}
.lt-k{font-weight:600;font-size:13.5px}
.lt-co{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:11.5px;color:var(--muted);margin-top:2px}
.lt-dots{display:inline-flex;gap:4px;align-items:center}
.lt-dot{width:9px;height:9px;border-radius:50%;display:inline-block}
.lt-na{color:var(--faint)}
.lt-bar-w{display:flex;width:150px}
.lt-ebar{display:flex;width:150px;height:8px;background:#f0f2f4}
.lt-ebar span{display:block;height:8px}
.lt-wrap table.dataTable td{vertical-align:middle}
.lt-wrap table.dataTable tbody tr{cursor:pointer}
.navbar .navbar-nav{flex-wrap:nowrap}
.navbar .nav-link{font-size:13px;padding:.35rem .5rem}
.navbar .navbar-brand{margin-right:10px}
.brand-logo{height:30px}
.brand-b{font-size:15px}
@media (min-width:1200px){.lt-bar{flex-wrap:nowrap}.lt-stat{flex:0 1 auto;white-space:nowrap}}
.lt-stat{flex:1 1 140px}

.cb-wrap{margin-top:16px;overflow-x:auto}
.cb-grid{display:grid;grid-template-columns:168px minmax(160px,1fr) repeat(6,46px) 64px;gap:0 10px;align-items:center;min-width:700px}
.cb-h{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);padding-bottom:8px}
.cb-n{text-align:right;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px}
.cb-v{color:var(--muted);padding:6px 0}
.cb-tot{color:var(--text);font-weight:600;font-size:13px}
.cb-lab{font-size:13px;display:flex;align-items:center;gap:9px;padding:6px 8px 6px 0}
.cb-dot{width:9px;height:9px;border-radius:50%;display:inline-block;flex:0 0 auto}
.cb-bar{display:flex;height:20px;border-left:1px solid #c2cace}
.cb-bar span{display:block;height:20px}
.cb-side{flex:0 0 280px;border-left:1px solid var(--rule);padding-left:22px}
.cb-legend{display:grid;grid-template-columns:auto auto 1fr;gap:9px 10px;margin-top:12px;align-items:flex-start}
.cb-legend .sw{width:11px;height:11px;display:inline-block;margin-top:3px}
.cb-legend .tk{font-weight:600;font-size:12.5px}
.cb-legend .df{font-size:12.5px;line-height:1.4;color:var(--muted)}
.gp-top{background:#fff;border:1px solid var(--line);border-bottom:0;padding:22px 40px 20px;display:flex;align-items:flex-start;gap:28px;flex-wrap:wrap}
.gp-pick{flex:0 0 260px}
.gp-pick .form-group{margin-bottom:0}
.gp-pick label{font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);font-weight:600;margin-bottom:6px}
.gp-h{flex:1 1 520px;min-width:0}
.gp-h1{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
.gp-sym{margin:0;font-size:22px;font-weight:600;letter-spacing:-.015em;font-family:"IBM Plex Mono",ui-monospace,monospace}
.gp-co{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:13px;color:var(--muted)}
.gp-chips{display:flex;align-items:center;gap:8px;margin-top:12px;flex-wrap:wrap}
.gp-chip{font-size:11.5px;border:1px solid var(--line);color:var(--muted);padding:3px 8px}
.gp-chip.t{border-color:var(--cta);color:var(--cta)}
.gp-act{margin-left:auto;display:flex;gap:8px;align-items:center;padding-top:22px}
.gp-wrap{background:#f4f5f6;border:1px solid var(--line);border-top:0;padding:22px 40px 26px}
.gp-cards{display:flex;gap:20px;align-items:stretch;flex-wrap:wrap}
.gp-c1{flex:1 1 620px;min-width:0}
.gp-c2{flex:1 1 420px;min-width:0}
.gp-tbl{width:100%;border-collapse:collapse;font-size:13px;margin-top:6px}
.gp-tbl th{text-align:center;padding:11px 8px;font-weight:500;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);border-bottom:1px solid var(--line)}
.gp-tbl th:first-child{text-align:left;padding-left:0}
.gp-tbl td{padding:8px;border-bottom:1px solid var(--rule);text-align:center}
.gp-tbl td:first-child{text-align:left;padding-left:0}
.gp-yes{color:#1f7a5c;font-weight:600}
.gp-no{color:var(--faint)}
:root[data-theme="dark"]{--bg:#0e161c;--surface:#15202a;--text:#e8eef2;--muted:#95a5af;--faint:#95a5af;--line:#2b3a45;--rule:#1e2b34;--accent:#5ea0ee;--cta:#ef6f43}
:root[data-theme="dark"] body{background:var(--bg);color:var(--text)}
:root[data-theme="dark"] .navbar,
:root[data-theme="dark"] .dc-hero,
:root[data-theme="dark"] .dc-sec,
:root[data-theme="dark"] .panel,
:root[data-theme="dark"] .dl-main,
:root[data-theme="dark"] .dl-card,
:root[data-theme="dark"] .dl-top,
:root[data-theme="dark"] .doc,
:root[data-theme="dark"] .doc-side,
:root[data-theme="dark"] .doc-body,
:root[data-theme="dark"] .lt-top,
:root[data-theme="dark"] .tc-top,
:root[data-theme="dark"] .tc-card,
:root[data-theme="dark"] .gp-top,
:root[data-theme="dark"] .dc-cell,
:root[data-theme="dark"] .dc-chip,
:root[data-theme="dark"] .dc-btn,
:root[data-theme="dark"] .card{background:var(--surface);color:var(--text)}
:root[data-theme="dark"] .dc-band,
:root[data-theme="dark"] .lt-wrap,
:root[data-theme="dark"] .tc-wrap,
:root[data-theme="dark"] .gp-wrap,
:root[data-theme="dark"] .dl-code,
:root[data-theme="dark"] .rv-bar{background:var(--bg)}
:root[data-theme="dark"] .dc-grid3{background:var(--line)}
:root[data-theme="dark"] .dc-hero-wash{background:radial-gradient(ellipse 58% 70% at 44% 48%,rgba(21,32,42,.96) 0%,rgba(21,32,42,.88) 50%,rgba(21,32,42,.4) 100%)}
:root[data-theme="dark"] .dc-searchwrap{border-color:var(--muted);background:var(--surface)}
:root[data-theme="dark"] .dc-searchwrap .selectize-input{background:var(--surface)!important;color:var(--text)}
:root[data-theme="dark"] .dc-searchbtn{background:var(--accent);color:#0e161c}
:root[data-theme="dark"] .selectize-input,
:root[data-theme="dark"] .selectize-dropdown,
:root[data-theme="dark"] .form-control,
:root[data-theme="dark"] .form-select{background:var(--surface)!important;color:var(--text)!important;border-color:var(--line)!important}
:root[data-theme="dark"] table.dataTable,
:root[data-theme="dark"] table.dataTable td,
:root[data-theme="dark"] table.dataTable th{color:var(--text);border-color:var(--rule)!important}
:root[data-theme="dark"] table.dataTable tbody tr{background:var(--surface)}
:root[data-theme="dark"] table.dataTable tbody tr.odd{background:#182531}
:root[data-theme="dark"] .dataTables_wrapper .dataTables_paginate .paginate_button{color:var(--muted)!important}
:root[data-theme="dark"] .cb-bar{border-left-color:var(--line)}
:root[data-theme="dark"] .lt-ebar{background:#1e2b34}
:root[data-theme="dark"] .tc-svg .tc-bar{fill:#233140;stroke:#2b3a45}
:root[data-theme="dark"] .tc-svg .tc-cen{fill:#0e161c}
:root[data-theme="dark"] .tc-svg .tc-chr{fill:#e8eef2}
:root[data-theme="dark"] .kary rect{fill:#5ea0ee;opacity:.07}
:root[data-theme="dark"] .site-footer{background:var(--surface);color:var(--muted)}
:root[data-theme="dark"] .gp-yes{color:#5fd0a0}
.nav-dark{font-size:11.5px;letter-spacing:.06em;text-transform:uppercase}
:root[data-theme="dark"] .navbar .navbar-brand,
:root[data-theme="dark"] .brand-b,
:root[data-theme="dark"] .navbar .nav-link.active{color:var(--text)}
:root[data-theme="dark"] .navbar .nav-link{color:var(--muted)}
:root[data-theme="dark"] .statline,
:root[data-theme="dark"] .statline *{color:var(--text)}
:root[data-theme="dark"] .sk,
:root[data-theme="dark"] .dc-kicker,
:root[data-theme="dark"] .dc-eyebrow,
:root[data-theme="dark"] .dc-note,
:root[data-theme="dark"] .dc-sub,
:root[data-theme="dark"] .dc-figc{color:var(--muted)}
:root[data-theme="dark"] .bslib-sidebar-layout,
:root[data-theme="dark"] .bslib-sidebar-layout > .sidebar,
:root[data-theme="dark"] .bslib-sidebar-layout > .main,
:root[data-theme="dark"] .sidebar-content,
:root[data-theme="dark"] .tab-content,
:root[data-theme="dark"] .html-fill-container{background:var(--bg);color:var(--text)}
:root[data-theme="dark"] .bslib-sidebar-layout > .sidebar{background:var(--surface);border-color:var(--line)}
:root[data-theme="dark"] .checkbox label,
:root[data-theme="dark"] .radio label,
:root[data-theme="dark"] .control-label,
:root[data-theme="dark"] label{color:var(--muted)}
:root[data-theme="dark"] .selectize-dropdown .active{background:#22313e}
:root[data-theme="dark"] .dataTables_wrapper .dataTables_info,
:root[data-theme="dark"] .dataTables_wrapper .dataTables_length,
:root[data-theme="dark"] .dataTables_wrapper .dataTables_filter{color:var(--muted)}
:root[data-theme="dark"] .doc-body a{border-bottom-color:var(--line)}
:root[data-theme="dark"] .gp-tbl td,
:root[data-theme="dark"] .doc-tbl td,
:root[data-theme="dark"] .dl-tbl td{border-bottom-color:var(--rule)}
:root[data-theme="dark"] .ev-yes{color:#5fd0a0}
:root[data-theme="dark"] .dc-searchbtn{color:#0e161c}

@media (max-width:520px){body{overflow-x:hidden}.dataTables_wrapper{overflow-x:auto}.dc-counts{flex-wrap:wrap}.navbar .navbar-nav{flex-wrap:wrap}}
a.brand{display:inline-flex;align-items:center;gap:10px;text-decoration:none}
a.brand:hover .brand-b{color:var(--cta)}
.tierlegend{margin-top:18px;padding-top:14px;border-top:1px solid var(--rule)}
.tl-grid{display:grid;grid-template-columns:auto auto minmax(230px,1fr) auto auto minmax(230px,1fr);gap:9px 12px;margin-top:11px;align-items:flex-start}
.tl-grid .sw{width:12px;height:12px;display:inline-block;margin-top:3px}
.tl-grid .tk{font-weight:600;font-size:12.5px;color:var(--text)}
.tl-grid .df{font-size:12.5px;line-height:1.45;color:var(--muted)}
@media (max-width:1000px){.tl-grid{grid-template-columns:auto auto minmax(200px,1fr)}}
.cb-bar span{display:flex;align-items:center;justify-content:center;overflow:hidden;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:10.5px;font-weight:500}
.land-tip{position:absolute;pointer-events:none;z-index:20;background:var(--surface);border:1px solid var(--line);padding:9px 11px;font-size:12px;line-height:1.5;min-width:150px}
.land-tip .g{font-family:"IBM Plex Mono",ui-monospace,monospace;font-weight:600;font-size:13px}
.land-tip .m{color:var(--muted)}
.foot-theme{margin-top:12px}
.foot-theme a{font-size:11.5px;border:1px solid var(--line);padding:4px 10px;color:var(--muted)}

.sum2{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:34px;margin-top:4px}
.sum-card .dc-figc{margin-bottom:14px}
.rb{display:grid;grid-template-columns:auto minmax(0,1fr);gap:7px 12px;align-items:center;margin-top:10px}
.rb-l{font-size:12.5px;color:var(--text);white-space:nowrap;font-family:"IBM Plex Mono",ui-monospace,monospace}
.rb-b{display:flex;align-items:center;gap:8px;min-width:0}
.rb-b span{display:flex;align-items:center;justify-content:flex-end;height:18px;padding-right:7px;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:10.5px;font-weight:500;overflow:hidden;white-space:nowrap}
.rb-b .rb-o{background:none;color:var(--muted);padding:0;justify-content:flex-start}
.brand-b{font-family:"Source Serif 4",Georgia,"Times New Roman",serif;font-size:18px;font-weight:600;letter-spacing:-.005em;line-height:1}
.cb-grid{grid-template-columns:168px minmax(160px,1fr) 64px !important;min-width:420px !important}
@media (max-width:820px){
.dc-hero{padding:44px 18px 34px}
.dc-h1{font-size:30px}
.dc-lede{font-size:14px}
.dc-searchwrap{flex-direction:column;align-items:stretch}
.dc-searchbtn{justify-content:center;padding:12px 0}
.dc-counts{gap:26px;flex-wrap:wrap;margin-top:32px}
.dc-count-n{font-size:26px}
.dc-band,.dc-sec,.dl-top,.lt-top,.tc-top,.gp-top,.dl-wrap,.lt-wrap,.tc-wrap,.gp-wrap{padding-left:18px;padding-right:18px}
.dc-grid3{grid-template-columns:1fr}
.sum2{grid-template-columns:1fr;gap:26px}
.dc-sechead{flex-wrap:wrap}
.dc-exp{flex-wrap:wrap}
.doc{flex-direction:column}
.doc-toc{flex:1 1 auto;position:static;border-right:0;border-bottom:1px solid var(--line);padding:18px}
.doc-body{padding:22px 18px 30px}
.doc-tiers,.doc-callout{max-width:100%}
.tl-grid{grid-template-columns:auto auto minmax(0,1fr)}
.cb-wrap{overflow-x:auto}
.dl-wrap,.gp-cards,.tc-bar,.lt-bar,.rv-bar{flex-direction:column;align-items:stretch}
.dl-side,.lt-f,.tc-f,.rv-f,.gp-pick{flex:1 1 auto;width:100%;margin-left:0}
.lt-e,.tc-f3{margin-left:0}
.gp-act{margin-left:0;padding-top:12px}
.gp-top,.tc-top{flex-direction:column}
.dl-tbl td.f{white-space:normal}
.dl-tbl .cl{display:none}
.foot-in{flex-direction:column;gap:18px}
.navbar .navbar-nav{flex-wrap:wrap}
table.dataTable{font-size:12px}
}
@media (max-width:520px){
.dc-h1{font-size:25px}
.dc-try{gap:6px}
.dc-chip{font-size:11px;padding:4px 8px}
.brand-b{font-size:16px}
.brand-logo{height:24px}
}

.brand-logo{position:relative;top:3px}
a.brand{align-items:center}

.cb-grid{grid-template-columns:168px minmax(160px,1fr) !important}
.cb-bar .cb-tot{margin-left:8px;line-height:20px;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;font-weight:600;color:var(--text);white-space:nowrap}
.rb-hov{opacity:0;transition:opacity .12s ease;font-size:12px;color:var(--muted);font-family:"IBM Plex Mono",ui-monospace,monospace;white-space:nowrap}
.rb-b:hover .rb-hov{opacity:1}

.lt-srch{display:grid;grid-template-columns:minmax(260px,1fr) auto auto;gap:14px;align-items:end;margin-bottom:16px}
.lt-srch .form-group{margin-bottom:0}
.lt-sb{display:flex;gap:12px;align-items:center;padding-bottom:7px}
.lt-clear{font-size:13px;color:var(--muted);text-decoration:none}
.lt-sn{font-size:13px;color:var(--muted);padding-bottom:9px}
@media(max-width:720px){.lt-srch{grid-template-columns:1fr}}

.gp-clear{font-size:12.5px;color:var(--muted);text-decoration:none;display:inline-block;margin-top:7px}
.gp-clear:hover{color:var(--cta)}
.ctkey{display:flex;flex-wrap:wrap;gap:7px 20px;margin:0 0 16px;font-size:12.5px;color:var(--muted)}
.ctk-i{display:inline-flex;align-items:center;gap:7px}
.ctk-d{width:9px;height:9px;border-radius:50%;display:inline-block;flex:0 0 auto}
.ctk-h{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);margin-bottom:7px}

.tn-mod,.tl-mod{text-transform:none;letter-spacing:.01em}

.foot-in{display:grid;grid-template-columns:repeat(4,1fr);gap:44px;max-width:1080px;margin:0 auto;align-items:start}
.foot-col{flex:none;min-width:0}
.site-footer{background:#eef1f3;border-top:1px solid var(--line)}
.foot-logos{max-width:1080px;margin:30px auto 0;padding-top:22px;border-top:1px solid var(--line)}
.navbar .navbar-nav{margin-left:auto}
.navbar .navbar-brand,a.brand{margin-right:42px}
.dc-h1{font-family:"Source Serif 4",Georgia,serif;font-size:40px;line-height:1.14;letter-spacing:-.015em}
.dc-counts{flex-wrap:wrap;row-gap:22px}
.dc-lede{max-width:720px;margin:14px auto 0}
@media(max-width:900px){.foot-in{grid-template-columns:repeat(2,1fr);gap:30px}}
@media(max-width:560px){.foot-in{grid-template-columns:1fr}}

.navbar .navbar-collapse{justify-content:flex-end}
.navbar .navbar-nav{margin-left:auto}
.navbar .navbar-brand,a.brand{margin-right:44px}
.tn-srch{margin:4px 0 2px}
.tn-srch .form-group{margin-bottom:6px}

.tn-pick{margin:6px 0 4px}
.tn-pick .form-group{margin-bottom:8px}

.gloss{margin:2px 0 0}
.gloss dt{font-weight:600;font-size:13.5px;margin-top:13px}
.gloss dt:first-child{margin-top:0}
.gloss dd{margin:3px 0 0;color:var(--muted);font-size:13.5px;line-height:1.55;max-width:70ch}

.dc-kick{display:block;font-family:Inter,system-ui,sans-serif;font-size:16px;font-weight:600;letter-spacing:.04em;color:var(--muted);margin-bottom:8px}

.ibub{position:relative;display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;margin-left:6px;border-radius:50%;background:#dbe7f3;color:#1c5c96;font-size:10px;font-weight:700;font-style:normal;cursor:help;vertical-align:middle;line-height:1;font-family:Inter,system-ui,sans-serif;text-transform:none;letter-spacing:0}
.ibub:hover,.ibub:focus{background:#2a78d6;color:#ffffff;outline:none}
.ibub-t{position:absolute;bottom:calc(100% + 9px);left:50%;transform:translateX(-50%);width:262px;background:#233947;color:#ffffff;font-size:12.5px;font-weight:400;line-height:1.5;letter-spacing:0;text-transform:none;text-align:left;white-space:normal;padding:9px 11px;border-radius:6px;opacity:0;visibility:hidden;transition:opacity .12s;z-index:2000;pointer-events:none}
.ibub-t::after{content:"";position:absolute;top:100%;left:50%;transform:translateX(-50%);border:6px solid transparent;border-top-color:#233947}
.ibub:hover .ibub-t,.ibub:focus .ibub-t{opacity:1;visibility:visible}

.bt-wrap{display:grid;grid-template-columns:minmax(0,1fr) 290px;gap:34px;align-items:start;margin-top:20px}
.bt-in textarea{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:13px;line-height:1.6}
.bt-btns{display:flex;align-items:center;gap:18px;margin-top:12px}
.dc-btn-go{background:#2a78d6;color:#ffffff;border-color:#2a78d6}
.dc-btn-go:hover{background:#1c5c96;color:#ffffff}
.bt-side .kv{display:flex;justify-content:space-between;gap:14px;padding:7px 0;border-bottom:1px solid var(--rule);font-size:13.5px}
.bt-side .kv b{font-variant-numeric:tabular-nums}
.bt-side .kv:last-child{border-bottom:0}
@media(max-width:820px){.bt-wrap{grid-template-columns:1fr}}

.ql-sec{margin-top:36px;padding-top:28px;border-top:1px solid var(--line)}
.ql-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:34px;margin-top:26px}
.ql-grid .dc-figc{margin-bottom:10px;min-height:56px}
@media(max-width:900px){.ql-grid{grid-template-columns:1fr}.ql-grid .dc-figc{min-height:0}}

.fl-act{display:flex;flex-wrap:wrap;align-items:center;gap:10px;margin:14px 0 2px;font-size:13px;color:var(--muted)}
.fl-chip{display:inline-flex;align-items:center;gap:6px;background:#eef3f8;border:1px solid #d7e3ef;color:#1c5c96;border-radius:20px;padding:3px 11px;font-size:12.5px;font-weight:500}
.fl-clear{font-size:12.5px;font-weight:600;color:#2a78d6;cursor:pointer;text-decoration:none}
.fl-clear:hover{text-decoration:underline}

a.brand{gap:18px}
.brand-logo{margin-right:4px}
.foot-legal{max-width:1080px;margin:26px auto 0;padding-top:16px;border-top:1px solid var(--line);text-align:right;font-size:12px;color:var(--muted)}

.sb-key{display:flex;flex-wrap:wrap;gap:20px;margin:4px 0 16px}
.sb-ki{display:inline-flex;align-items:center;gap:7px;font-size:12.5px;color:var(--muted)}
.sb-sw{width:11px;height:11px;border-radius:2px;display:inline-block}
.sb-row{display:grid;grid-template-columns:104px 1fr 52px;align-items:center;gap:14px;margin-bottom:8px}
.sb-l{font-size:12px;color:var(--muted);text-align:right}
.sb-track{display:flex;height:17px;gap:2px}
.sb-seg{display:block;height:100%;border-radius:2px;transition:filter .12s}
.sb-seg:hover{filter:brightness(1.12)}
.sb-n{font-size:12.5px;text-align:right;font-variant-numeric:tabular-nums}
@media(max-width:620px){.sb-row{grid-template-columns:78px 1fr 44px;gap:10px}}

.sb-seg{display:flex;align-items:center;justify-content:center;color:#ffffff;font-size:10.5px;font-weight:600;font-variant-numeric:tabular-nums;overflow:hidden;white-space:nowrap}

.t6-strip{display:flex;align-items:baseline;gap:10px;margin-top:16px;padding:12px 14px;background:#f4f6f8;border-radius:7px;max-width:820px}
.t6-sw{width:11px;height:11px;border-radius:2px;background:#10375e;flex:0 0 auto;position:relative;top:2px}
.t6-k{font-weight:700;font-size:12.5px;color:#10375e}
.t6-t{font-size:13px;color:var(--muted);line-height:1.55}

.foot-bottom{max-width:1080px;margin:26px auto 0;padding-top:18px;border-top:1px solid var(--line);display:flex;align-items:center;justify-content:space-between;gap:26px;flex-wrap:wrap}
.foot-bottom .foot-logos{margin:0;display:flex;align-items:center;gap:22px}
.foot-legal{margin:0;padding:0;border:0;text-align:right;font-size:12px;color:var(--muted)}
@media(max-width:560px){.foot-bottom{justify-content:center}.foot-legal{text-align:center}}

.cmp-sec{margin-top:34px;padding-top:26px;border-top:1px solid var(--line)}
.cmp-pick{display:grid;grid-template-columns:1fr 1fr;gap:26px;margin:18px 0 8px;max-width:640px}
.cmp-tbl td:first-child{color:var(--muted);width:190px}
.cmp-tbl th:nth-child(2),.cmp-tbl th:nth-child(3){font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px}
.cmp-same td{opacity:.48}
.gp-rowlink{cursor:pointer}
.gp-rowlink:hover td{background:#eef3f8}
@media(max-width:620px){.cmp-pick{grid-template-columns:1fr;gap:12px}}

.brand-sep{display:inline-block;width:1px;height:19px;background:var(--line);margin:0 4px}
a.brand{gap:14px}
.tc-bar .checkbox-inline,.tc-bar label,.tc-f .checkbox-inline{text-transform:none;letter-spacing:0}
.sb-seg{position:relative}
.sb-bub{position:absolute;bottom:calc(100% + 8px);left:50%;transform:translateX(-50%);width:auto;white-space:nowrap;background:#233947;color:#ffffff;font-size:12px;font-weight:500;line-height:1.4;padding:6px 10px;border-radius:5px;opacity:0;visibility:hidden;transition:opacity .12s;z-index:2000;pointer-events:none}
.sb-bub::after{content:"";position:absolute;top:100%;left:50%;transform:translateX(-50%);border:5px solid transparent;border-top-color:#233947}
.sb-seg:hover .sb-bub{opacity:1;visibility:visible}

.foot-bottom{justify-content:flex-end}
.foot-mark{height:26px;width:auto;display:block}
.foot-h .foot-mark{margin-bottom:2px}
.ti-wrap{display:flex;flex-direction:column;gap:10px;margin-top:8px}
.ti{width:100%;height:auto;max-width:860px}

'))

app_footer <- tags$footer(class = "site-footer",
  div(class = "foot-in",
    div(class = "foot-col",
      div(class = "foot-h", img(src = "logo_statfungen.png", class = "foot-mark",
                            alt = "StatFunGen Lab")),
      div(class = "foot-t",
          "Laboratory of Statistical Functional Genomics, led by Gao Wang, in the ",
          "Gertrude H. Sergievsky Center at Columbia University.")),
    div(class = "foot-col",
      div(class = "foot-h", "Release"),
      div(class = "foot-t", "Data release ", DATA_RELEASE),
      div(class = "foot-t", "Genome build ", GENOME_BUILD)),
    div(class = "foot-col",
      div(class = "foot-h", "Resources"),
      tags$a(href = "https://statfungen.github.io/xqtl-resources/", target = "_blank",
             rel = "noopener", "Dataset catalog"),
      tags$a(href = "https://github.com/StatFunGen/xQTL-AD-loci-explorer", target = "_blank",
             rel = "noopener", "Source code"),
             tags$a(href = "https://www.synapse.org/Synapse:syn68872650", target = "_blank",
                    rel = "noopener", "Data on Synapse"),
             tags$a(href = "https://github.com/StatFunGen/xqtl-protocol", target = "_blank",
                    rel = "noopener", "Analysis protocol")),
    div(class = "foot-col",
      div(class = "foot-h", "Cite"),
      div(class = "foot-t",
          "The Alzheimer’s Disease Functional Genomics (FunGen-AD) Consortium. ",
          tags$em("Broad and deep dissection of Alzheimer’s disease genetics ",
                  "with FunGen-xQTL.")))),
      div(class = "foot-bottom",
        div(class = "foot-legal",
          "© 2026 StatFunGen Lab, Columbia University. All rights reserved.")))

ui <- page_navbar(
  id = "nav",
  title = tags$a(id = "brand_home", href = "#", class = "brand action-button",
              img(src = "logo_fungen.png", class = "brand-logo", alt = "FunGen-AD"),
              span(class = "brand-sep"), span(class = "brand-b", "AD Loci Explorer")),
  theme = app_theme, fillable = FALSE,
  header = tagList(app_css, extra_css, nav_js, a11y_js, vs_css, mob_css, conn_js, cue_js, tour_js, ax_css, ax_js, tbl_js, igv_a11y_js, fig_js, uiOutput("view_strip"), igv_head,
                   tags$link(rel = "stylesheet", href = paste0(
                     "https://fonts.googleapis.com/css2?family=Source+Serif+4:",
                     "opsz,wght@8..60,600;8..60,700&display=swap"))),
  footer = app_footer,

  # ---- Home -----------------------------------------------------------------
  nav_panel(
    "Home", icon = icon("magnifying-glass"),
      div(class = "dc-hero",
        div(class = "dc-hero-wash"),
        div(class = "dc-hero-in",
          h1(class = "dc-h1", span(class = "dc-kick", "ADSP FunGen-xQTL"), "AD Loci Explorer"),
          p(class = "dc-lede",
            "The FunGen-xQTL release links each of ", format(n_loci_total, big.mark = ","), " candidate ",
            "Alzheimer’s disease (AD) loci to its target genes, and shows the brain and immune ",
            "contexts, molecular modalities and evidence tier behind each link. Some loci do not ",
            "reach genome-wide significance on their own; the Loci tab marks them. ",
            "Start with a gene symbol, Ensembl ID, rsID or genomic region."),
          div(class = "note",
            "New in this release: transcriptomic pattern QTL (tpQTL) from three ROSMAP brain ",
            "regions (anterior cingulate, dorsolateral prefrontal and posterior cingulate cortex). ",
            "A pattern describes how read coverage varies along one gene, so its target gene ",
            "is the gene it was built from."),
          div(class = "dc-searchwrap",
            selectizeInput("hero_search", NULL, choices = NULL, width = "100%",
              options = list(placeholder = "Search gene, Ensembl ID, rsID or genomic region",
                maxOptions = 25, create = TRUE, searchField = c("label", "gid"),
                onInitialize = I('function() { this.setValue(""); }'))),
            actionLink("hero_go", "Search", class = "dc-searchbtn")),
          div(class = "dc-try", span("Try"),
            actionLink("ex_bin1", "BIN1", class = "dc-chip dc-mono"),
            actionLink("ex_rs", "rs429358", class = "dc-chip dc-mono"),
            actionLink("ex_region", "chr2:127088929-127118906", class = "dc-chip dc-mono"),
            actionLink("ex_ensg", "ENSG00000130203", class = "dc-chip dc-mono")),
          div(class = "dc-counts",
            div(div(class = "dc-count-n dc-mono", format(n_loci_total, big.mark = ",")),
                div(class = "dc-count-l", "AD loci")),
            div(div(class = "dc-count-n dc-mono", format(n_t1_t6, big.mark = ",")),
                div(class = "dc-count-l", info("genes assigned T1-T6", "Genes with an evidence tier from T1 to T6 at one or more AD loci, as in the gene tier file in the Download tab. T6 genes have TWAS or MR evidence only."))),
            div(div(class = "dc-count-n dc-mono", textOutput("kpi_pct", inline = TRUE)),
                div(class = "dc-count-l", info("loci with evidence", "Share of AD loci with at least one target gene supported by xQTL evidence."))),
            div(div(class = "dc-count-n dc-mono", textOutput("kpi_trans", inline = TRUE)),
                div(class = "dc-count-l", info("distinct trans pairs", "An AD locus and a target gene on another chromosome, linked by trans-xQTL evidence. Each pair counts once here; the Trans tab lists one row per supporting modality, so it shows more rows."))),
            div(div(class = "dc-count-n dc-mono", "263,893"),
                div(class = "dc-count-l", info("cis-xQTL credible sets", "A 95% credible set is the smallest set of variants with at least a 95% probability of containing the causal variant. cis means the variant lies near the gene it affects.")))))),
      div(class = "dc-band",
        div(class = "dc-kicker", "Start here"),
        div(class = "dc-note", style = "margin:-4px 0 14px",
            tags$a(class = "tour-launch", href = "#", "Take the two-minute tour"),
            " if this is your first visit."),
        div(class = "dc-grid3",
          div(class = "dc-cell",
            div(class = "dc-cell-t", "You have a candidate gene"),
            div(class = "dc-cell-b",
              "Enter a gene symbol or Ensembl ID above to open the AD locus where the gene has ",
              "its strongest evidence, with the supporting variants alongside."),
            actionLink("go_genes", "Search a gene →", class = "dc-cell-c")),
          div(class = "dc-cell",
            div(class = "dc-cell-t", "You are following a GWAS signal"),
            div(class = "dc-cell-b",
              "Enter an rsID or a region such as chr2:127088929-127118906 to open the AD locus ",
              "that contains it."),
            actionLink("go_loci", "Search a region →", class = "dc-cell-c")),
          div(class = "dc-cell",
            div(class = "dc-cell-t", "You are browsing"),
            div(class = "dc-cell-b",
              "Open the locus table, then click a gene to see the variants behind it."),
            actionLink("go_explore", "Browse all loci →", class = "dc-cell-c")),
          div(class = "dc-cell",
            div(class = "dc-cell-t", "You have a list of genes"),
            div(class = "dc-cell-b",
                "Paste up to 300 gene symbols, rsIDs or regions to see which are in this release."),
            actionLink("go_batch", "Match a list →", class = "dc-cell-c"))),
        div(class = "dc-note", style = "margin-top:16px",
          "Tiers rank how much evidence links a target gene to an AD locus. They prioritize genes ",
          "for follow-up and do not establish biological causality or clinical relevance.")),
      div(class = "dc-sec",
        div(class = "dc-sechead",
          div(h2(class = "dc-h2", "Release overview"),
              div(class = "dc-sub",
                "Where the AD loci lie in the genome, and how their evidence spreads across contexts."))),
        div(style = "margin-top:26px",
          div(class = "dc-sechead", style = "border-bottom:0;padding-bottom:0",
            div(div(class = "dc-figt", info("AD loci across the genome", "One point per AD locus at its genomic position. Height is the significance of the lead GWAS variant; size is the number of target genes.")),
                div(class = "dc-figc",
                  "Height is the significance of the lead GWAS variant; color is the ",
                  "best evidence tier at the locus.")),
            div(class = "dc-exp", span(class = "dc-kicker", "Export"),
              downloadLink("dl_figA_pdf", "PDF", class = "dc-btn"),
              downloadLink("dl_figA_png", "PNG", class = "dc-btn"),
              downloadLink("dl_figA_csv", "CSV", class = "dc-btn"))),
          div(style = "margin-top:14px",
            div(style = "position:relative",
              withSpinner(plotOutput("p_landscape", height = 340,
                  hover = hoverOpts("land_hover", delay = 60, delayType = "debounce"),
                  click = "land_click"),
                type = 8, color = "#2a78d6", size = 0.6),
              uiOutput("land_tip"))),
          uiOutput("land_kbd"),
          div(class = "dc-note", style = "margin-top:10px;max-width:760px",
            "Point size is the number of target genes. Loci without a variant-level ",
            "p-value are not drawn.")),
        div(style = "margin-top:40px;padding-top:28px;border-top:1px solid var(--line)",
          div(class = "dc-sechead", style = "border-bottom:0;padding-bottom:0",
            div(div(class = "dc-figt", info("From context to modality to evidence tier", "Each ribbon is a group of gene records. Width is the number of records; color is their evidence tier.")),
                div(class = "dc-figc",
                  "Each ribbon follows a group of gene records from the context where they were measured, ",
                  "through the modality, to their evidence tier. The three columns show the same ",
                  format(nrow(.sk_data), big.mark = ","), " records three ways. Hover a ribbon to see ",
                  "its genes and loci, or drag a node to separate crowded ribbons. T6 is left out ",
                  "because T6 records have TWAS or MR evidence only, with no context or modality."),
            div(class = "dc-exp", span(class = "dc-kicker", "Export"),
              downloadLink("dl_figB_csv", "CSV", class = "dc-btn")))),
          div(style = "margin-top:14px",
            plotlyOutput("home_river", height = "600px"), uiOutput("river_kbd")),
          div(class = "t6-strip",
              span(class = "t6-sw"),
              span(class = "t6-k dc-mono", "T6"),
              span(class = "t6-t", uiOutput("t6_n", inline = TRUE),
                   " gene records are T6. They have TWAS or MR evidence only, with no context ",
                   "or modality to show above."))),
        div(class = "tierlegend",
          div(class = "dc-kicker", "Evidence tier"),
          div(class = "tl-grid",
            lapply(TIER_SEQ, function(k)
              tagList(
                span(class = "sw", style = sprintf("background:%s", unname(conf_pal[[k]]))),
                span(class = "tk dc-mono", k),
                span(class = "df", TIER_DEFS[[k]]))),
            span(class = "sw", style = "background:#cfd6db"),
            span(class = "tk dc-mono", "—"),
            span(class = "df", "No tier assigned at this locus.")))),

        div(class = "dc-sec",
          div(class = "dc-figt", info("AD loci with the most trans target genes", "A trans target gene lies on a different chromosome from the AD locus.")),
          div(class = "dc-figc",
            "Each bar counts the distinct trans target genes of one AD locus, with the count ",
            "after the bar. Bar length uses a log scale because the top locus has more than a ",
            "hundred times as many trans target genes as most others. Segments split the count ",
            "by type of trans evidence, so a gene supported by two modalities appears in both."),
          uiOutput("sum_trans")),
  ),

  # ---- Explore -------------------------------------------------------------
    nav_panel(
      "Loci", icon = icon("table-list"),
      div(class = "lt-top",
        h1(class = "doc-h", style = "font-size:22px;margin:0", "All AD loci"),
        div(class = "dc-sub", style = "max-width:860px",
          format(n_loci_total, big.mark = ","),
          " candidate AD loci in this release, each with its strongest xQTL evidence. Sort by any ",
          "column; the evidence bar shows how the target genes at a locus spread across tiers."),
        div(class = "lt-srch",
          div(class = "lt-sf",
            textInput("lt_q", "Find a locus", value = "", width = "100%",
                      placeholder = "Gene symbol, Ensembl ID, rsID or region such as chr2:127088929-127118906")),
          div(class = "lt-sb", actionLink("lt_qgo", "Search", class = "dc-searchbtn"),
              actionLink("lt_qclear", "Clear", class = "lt-clear")),
          div(class = "lt-sn", uiOutput("lt_qnote"))),
        div(class = "lt-bar",
          div(class = "lt-f",
            selectInput("lt_chr", "Chromosome",
              choices = c("All chromosomes" = "", stats::setNames(paste0("chr", CHR_ORD),
                          paste0("chr", CHR_ORD))), width = "100%")),
          div(class = "lt-f",
            selectizeInput("lt_ctx", info("Context", "Brain or immune context of the molecular measurement. Bulk means unsorted tissue rather than a sorted cell type."), choices = NULL, multiple = TRUE,
              width = "100%", options = list(placeholder = "All contexts"))),
          div(class = "lt-f",
            selectizeInput("lt_tier", info("Best tier", "Strongest evidence tier at this locus; T1 is the most stringent. T6 means TWAS or MR evidence without a fine-mapped AD-xQTL variant, so it sits outside the T1 to T5 order."), choices = TIER_SEQ, multiple = TRUE,
              width = "100%", options = list(placeholder = "All tiers"))),
          div(class = "lt-f",
            selectizeInput("lt_mod", info("Modality", "Molecular modality of the evidence. Expression (eQTL), splicing (sQTL), transcriptomic pattern (tpQTL), protein abundance (pQTL), DNA methylation (mQTL), histone acetylation (haQTL), chromatin accessibility (caQTL), glycoprotein abundance (gpQTL) or metabolite abundance (metQTL)."), choices = MOD_ORD, multiple = TRUE,
              width = "100%", options = list(placeholder = "All modalities"))),
              div(class = "lt-f",
                selectizeInput("lt_cis", info("Evidence", "cis means the variant acts on a nearby gene. trans means the target gene lies elsewhere in the genome, often on another chromosome."), choices = c("cis only", "trans only", "cis + trans", "no gene"), multiple = TRUE,
                  width = "100%", options = list(placeholder = "Any evidence"))),
          div(class = "lt-f lt-stat", uiOutput("lt_stat")),
          div(class = "lt-f lt-e",
            div(class = "dc-kicker", "Export"),
            div(class = "dc-exp", style = "margin-top:6px",
              downloadLink("dl_loci_csv", "CSV", class = "dc-btn"))))),
      uiOutput("lt_active"),
      uiOutput("ct_key"),
      div(class = "lt-wrap",
        div(class = "panel",
          withSpinner(DTOutput("loci_tbl"), type = 8, color = "#2a78d6", size = 0.6),
          div(class = "dc-note", style = "margin-top:12px",
            "Genes counts the distinct target genes at the locus; records counts the gene-variant ",
            "rows behind them. Click a row to open the locus."))),
          div(class = "cmp-sec",
            div(class = "dc-sechead",
              div(class = "dc-figt", info("Compare two loci",
                "A side by side read of two loci. Rows that match on both sides are dimmed, so what differs stands out.")),
              div(class = "dc-figc", "Pick any two AD loci to see how they differ.")),
            div(class = "dc-exp", downloadLink("dl_cmp", "CSV", class = "dc-btn")),
            div(class = "cmp-pick",
              selectInput("cmp_a", "First locus", choices = loci, selected = loci[1], width = "100%"),
              selectInput("cmp_b", "Second locus", choices = loci, selected = loci[2], width = "100%")),
            uiOutput("cmp_tbl"))),
  nav_panel(
    "Locus detail", icon = icon("dna"),
    div(class = "sheet",
      div(class = "locus-bar",
        div(class = "locus-pick",
            selectInput("locus_pick", "Select AD locus", choices = loci),
            div(style = "display:flex; gap:6px; margin-top:-4px;",
                actionButton("locus_prev", "← Prev", class = "btn-sm btn-outline-secondary"),
                actionButton("locus_next", "Next →", class = "btn-sm btn-outline-secondary"))),
        uiOutput("locus_meta")),
        uiOutput("locus_interp"),
        uiOutput("cite_box"),
      div(class = "panel",
        div(class = "panel-title", "Genome"),
        div(id = "igv_div", style = "min-height:120px")),
      div(class = "panel",
        div(class = "dc-sechead",
          div(
            div(style = "display:flex;align-items:baseline;gap:12px;flex-wrap:wrap",
              h3(class = "dc-h2", style = "font-size:17px", "Regional evidence"),
              span(class = "dc-mono", style = "font-size:12px;color:var(--muted)",
                   textOutput("rv_win", inline = TRUE))),
            div(class = "dc-figc",
              "Figure 2 · Credible-set variants and the xQTL evidence for each target gene ",
              "in this window.")),
          div(class = "dc-exp", span(class = "dc-kicker", "Export"),
            downloadLink("dl_region_pdf", "PDF", class = "dc-btn"),
            downloadLink("dl_region_png", "PNG", class = "dc-btn"),
            downloadLink("dl_region_csv", "CSV", class = "dc-btn"))),
        div(class = "rv-bar",
          div(class = "rv-f",
            selectizeInput("rv_ctx", info("Context", "Brain or immune context of the molecular measurement."), choices = NULL, multiple = TRUE,
              width = "100%", options = list(placeholder = "All contexts"))),
          div(class = "rv-f",
            textInput("rv_q", "Find", value = "", width = "100%",
                      placeholder = "Gene or rsID")),
          div(class = "rv-f",
            selectizeInput("rv_mod", info("Modality", "Molecular modality of each row, such as expression, splicing, protein abundance, DNA methylation, histone acetylation, chromatin accessibility or glycoprotein abundance."), choices = NULL, multiple = TRUE,
              width = "100%", options = list(placeholder = "All modalities"))),
          div(class = "rv-f rv-s",
            sliderInput("rv_tier", info("Tiers shown", "Show rows up to this tier. Lower numbers mean more stringent evidence."), min = 1, max = 6, value = 6,
              step = 1, ticks = FALSE, width = "100%")),
          div(class = "rv-f rv-b",
            actionLink("rv_reset", "Reset filters", class = "dc-btn"))),
        uiOutput("rv_stat"),
        withSpinner(plotOutput("p_region", height = 620),
                    type = 8, color = "#2a78d6", size = 0.6),
        div(class = "dc-howto",
          div(class = "dc-kicker", "How to read this figure"),
          div(class = "dc-howto-g",
            span(class = "n", "1"),
            span("Find the peak in the upper track. That is where the AD association is strongest."),
            span(class = "n", "2"),
            span("Look straight down. Marks under the peak are xQTL evidence for the gene in that row."),
            span(class = "n", "3"),
            span("The row label gives the context and modality. Dimmed marks are ",
                 "filtered, not deleted.")))),
      div(class = "panel",
        div(class = "dc-sechead",
          div(class = "panel-title", "Genes and evidence at this locus"),
          div(class = "dc-exp", downloadLink("dl_locus_tbl", "CSV", class = "dc-btn"))),
        withSpinner(DTOutput("locus_tbl"), type = 8, color = "#2c7fb8", size = 0.6)),
      div(class = "panel-note", style = "margin-top:-4px",
        "Trans target genes for this locus are drawn on the ",
        actionLink("go_trans", "Trans"), " tab."))
  ),

  # ---- Trans ---------------------------------------------------------------
  nav_panel(
    "Genes", icon = icon("table"),
      div(class = "gp-top",
        div(class = "gp-pick",
          selectizeInput("gp_gene", info("Gene", "Pick a gene to see its AD loci, contexts and modalities."), choices = NULL, width = "100%",
            options = list(placeholder = "Pick a gene", maxOptions = 25)),
          actionLink("gp_clear", "Clear", class = "gp-clear")),
        uiOutput("gp_head")),
      div(class = "gp-wrap",
        div(class = "gp-cards",
          div(class = "panel gp-c1",
            div(class = "dc-sechead",
              div(div(class = "dc-figt", "Evidence by context and modality"),
                  div(class = "dc-figc",
                    "A tick means this release records evidence for the gene in that context ",
                    "and modality; a dash means no record. A dash can mean either tested and null ",
                    "or not tested.")),
              div(class = "dc-exp", downloadLink("dl_gp_matrix", "CSV", class = "dc-btn"))),
            uiOutput("gp_matrix")),
          div(class = "panel gp-c2",
            div(class = "dc-sechead",
              div(div(class = "dc-figt", "AD loci linked to this gene"),
                  div(class = "dc-figc", "Both cis and trans links.")),
              div(class = "dc-exp", downloadLink("dl_gp_loci", "CSV", class = "dc-btn"))),
            uiOutput("gp_loci")))),
    layout_sidebar(
      sidebar = explore_sidebar,
      div(class = "statline",
          span(class = "sk", "Rows shown"), textOutput("e_rows", inline = TRUE),
          span(class = "sk", "Loci"), textOutput("e_loci", inline = TRUE),
          span(class = "sk", "Genes"), textOutput("e_genes", inline = TRUE),
          span(class = "sk", "Trans"), textOutput("e_trans", inline = TRUE)),
        div(class = "viewbar",
            radioButtons("view_mode", NULL, inline = TRUE,
                         choices = c("Genes" = "gene", "Variant-gene rows" = "variant"),
                         selected = "gene")),
        tier_key(),
        ev_key(),
        conditionalPanel("input.view_mode == 'gene'",
          card(full_screen = TRUE,
               card_header(div(style = "display:flex;justify-content:space-between;align-items:center;gap:12px",
                 span("Gene summary  ·  click a gene to see its variant rows"),
                 downloadLink("dl_gene_sum", "CSV", class = "dc-btn"))),
               withSpinner(DTOutput("gene_tbl"), type = 8, color = "#2c7fb8", size = 0.6))),
        conditionalPanel("input.view_mode == 'variant'",
      card(full_screen = TRUE,
           card_header("Variant–gene table  ·  click a row for details"),
           withSpinner(DTOutput("tbl"), type = 8, color = "#2c7fb8", size = 0.6)),
      uiOutput("detail"),
        )
    )
  ),

  # ---- Locus detail ---------------------------------------------------------
    panel_cell_types(),
    nav_panel(
      "Trans", icon = icon("share-nodes"),
      div(class = "tc-top",
        h1(class = "doc-h", style = "font-size:22px;margin:0", "Trans associations"),
        div(class = "dc-sub", style = "max-width:820px",
          "AD loci with trans-xQTL evidence for target genes on other chromosomes. The upper ",
          "chromosome holds the locus and the lower ones hold its trans target genes, drawn to scale."),
        div(class = "tc-bar",
          div(class = "tc-f",
            selectInput("trans_locus", info("Table and export locus", "Sets the locus for the trans pair table and its CSV export, and is shared with the Locus detail tab. The circle figure below has its own selector."), choices = loci, width = "100%"),
            div(class = "scope-note",
                "Sets the trans table and its CSV export. Shared with Locus detail.")),
          div(class = "tc-f tc-f2",
            div(role = "group", `aria-label` = "Trans modality filter",
              checkboxGroupInput("tmod", "Modality",
                choices = c("Trans genes", "Transmap"), inline = TRUE))),
          div(class = "tc-f tc-f3",
            div(class = "dc-kicker", "Export"),
            div(class = "dc-exp", style = "margin-top:6px",
              downloadLink("dl_trans", "CSV", class = "dc-btn"))))),
      div(class = "tc-wrap",
                  div(style = "margin-top:34px;padding-top:26px;border-top:1px solid var(--line)",
            div(class = "dc-figt", info("Trans target genes of one locus", "Each arc joins the locus to one of its trans target genes. Orange marks the locus.")),
            div(class = "dc-figc",
              "Pick a locus to see all its trans target genes. Orange marks the locus; each arc ",
              "runs to a target gene named around the rim and is colored by type of trans evidence."),
            div(class = "tn-pick",
              selectInput("tn_loc", info("Figure locus", "Sets the locus drawn in the circle below. Loci are listed by their number of distinct trans target genes."), choices = trans_loci_choices,
                          selected = if (length(trans_loci_choices)) trans_loci_choices[[1]] else NULL,
                          width = "320px"),
              checkboxInput("tn_sync", "Keep the figure on the table locus", value = TRUE),
              radioButtons("tn_view", info("Figure style", "Circle draws the chromosomes in a ring. Bars shows the same target genes as a bar chart, which is easier to read with a screen reader or in print."),
                           choices = c("Circle", "Bars"), selected = "Circle", inline = TRUE),
              div(class = "scope-note",
                  "This selector changes only the circle figure. Untick the box above to draw ",
                  "a different locus from the table."),
              uiOutput("tn_syncnote")),
            uiOutput("tn_figure")),
        div(class = "panel", style = "margin-top:20px",
          div(class = "dc-sechead",
            div(div(class = "dc-figt", "All trans pairs"),
                div(class = "dc-figc",
                  "One row per locus, target gene and supporting modality."))),
          div(class = "tn-srch",
            textInput("tn_q", info("Find", "Matches the locus, source gene, rsID, target gene or context, together with the scope and modality settings."),
                      value = "", width = "320px",
                      placeholder = "Gene, target, rsID or AD locus"),
            radioButtons("tn_scope", info("Scope", "This locus shows the locus chosen above. All loci shows every trans pair in the release."), inline = TRUE,
                         choices = c(`This locus` = "locus", `All loci` = "all"),
                         selected = "locus")),
          uiOutput("tn_active"),
        withSpinner(DTOutput("trans_tbl"), type = 8, color = "#2a78d6", size = 0.6)))),

  # ---- Documentation --------------------------------------------------
  # ---- Download -------------------------------------------------------------
    panel_batch(),
    nav_panel(
      "Download", icon = icon("download"),
      div(class = "dl-top",
        h1(class = "doc-h", style = "font-size:22px;margin:0", "Downloads"),
        div(class = "dc-sub", style = "max-width:840px",
          "Every table behind this browser as a flat file. Data release 2026-10-02, genome ",
          "build GRCh38, 1-based coordinates. Large source files are linked from the dataset catalog.")),
        div(class = "ql-sec",
          div(class = "dc-sechead",
            div(class = "dc-figt", info("How well resolved is this release",
              "These three views show how precisely fine-mapping resolves the evidence across the release.")),
            div(class = "dc-figc",
              "Check these before relying on the release tables.")),
          uiOutput("ql_stats"),
          div(class = "ql-grid",
            div(div(class = "dc-figt", "Fine-mapping confidence"),
                div(class = "dc-figc",
                  "Highest PIP of any variant in each gene record. Records on the right are ",
                  "resolved to few variants; records on the left are spread over many."),
                plotOutput("ql_pip", height = "195px")),
            div(div(class = "dc-figt", "Variant-to-function score"),
                div(class = "dc-figc",
                  "cV2F scores how likely the variant is to affect function, independently of ",
                  "the xQTL evidence. Records without a score are left out."),
                plotOutput("ql_cv2f", height = "195px")),
            div(div(class = "dc-figt", "GWAS support"),
                div(class = "dc-figc",
                  "GWAS significance of the variant behind each xQTL record. Suggestive and ",
                  "non-significant records are kept."),
                plotOutput("ql_sig", height = "195px"))),
          div(class = "dc-note", style = "margin-top:22px;max-width:860px",
            tags$b("APOE caveat. "),
            "APOE lies in a region of strong, extended linkage disequilibrium (LD), where ",
            "fine-mapping separates neighboring variants poorly. Read a credible set at this ",
            "locus as a region of interest, not a short list of candidate causal variants.")),
      div(class = "dl-wrap",
        div(class = "dl-main",
          div(class = "dc-sechead",
            div(div(class = "dc-figt", "Release files"),
                div(class = "dc-figc",
                  "Served from this browser; sizes are as shipped.")),
            div(class = "dc-exp",
              downloadLink("dl_zip", "Download all files (zip)", class = "dc-btn"))),
          dl_group("Browser tables",
            "what the locus, gene and trans views are built from",
            list(
              dl_row("dlf_browser", "AD locus evidence table",
                "The exact table behind the locus and gene views in this browser.",
                format(nrow(dat), big.mark = ","), "CSV", "2.9 MB",
                "ADlocus, gene, rsid, chr, pos, log10pval, top_confidence, ordered_contexts, trans_genes"),
              dl_row("dlf_genepos", "Gene coordinates",
                "GRCh38 coordinates for every gene referenced in this release.",
                if (is.null(gene_pos)) "—" else format(nrow(gene_pos), big.mark = ","),
                "CSV", "140 KB", "gene, chr, start, end"))),
          dl_group("AD locus tables",
            "the 2026-10-02 unified AD loci release",
            list(
              dl_row("dlf_locsum", "AD locus summary",
                "One row per AD locus: region, lead variant and gene count.",
                "195", "CSV", "19 KB",
                "AD_locus, region, chr, lead_variant, best_tier, top_gene, cell_types, modalities, n_genes, n_records"),
              dl_row("dlf_tier", "Gene tier assignment",
              "The evidence tier assigned to each gene in this release.",
              "495", "CSV", "12 KB",
              "gene_ID, gene_name, tier"),
            dl_row("dlf_varlvl", "AD locus variants",
                "Variant-level unified AD loci table behind the locus definitions.",
                "7,618", "CSV.GZ", "396 KB",
                "chr, pos, ADlocus, variant_ID, GWAS_methods, max_variant_inclusion_probability, gwas_sources, is.cs95, max_zscore"))),
          dl_group("xQTL summary",
            "AD locus by xQTL evidence, as a workbook",
            list(
              dl_row("dlf_xlsx", "AD locus and xQTL summary workbook",
                "Unified AD locus and xQTL summary, one sheet per evidence type.",
                "—", "XLSX", "8.8 MB",
                "see the workbook header row"))),
          dl_group("In the catalog",
            "too large to serve from this browser",
            list(
              dl_ext("Full xQTL overlap with AD loci",
                "Every xQTL method overlapped against the AD loci, unfiltered.",
                "CSV.GZ", "3.2 GB", "https://statfungen.github.io/xqtl-resources/"),
              dl_ext("Trans colocalization summary",
                "Trans-only colocalization summary across the full xQTL set.",
                "BED", "1.9 GB", "https://statfungen.github.io/xqtl-resources/"),
              dl_ext("Trans top loci",
                "Exported trans top loci from the FunGen-xQTL pipeline.",
                "BED.GZ", "159 MB", "https://statfungen.github.io/xqtl-resources/"))),
          div(class = "dl-foot",
            span("Each file carries the release label and genome build in its name."),
            tags$a(href = "https://statfungen.github.io/xqtl-resources/",
                   target = "_blank", rel = "noopener", "Dataset catalog"),
            tags$a(href = "https://github.com/StatFunGen/xQTL-AD-loci-explorer",
                   target = "_blank", rel = "noopener", "Source code"))),
        div(class = "dl-side",
          div(class = "dl-card",
            div(class = "dl-card-t", "Reading a file"),
            div(class = "dl-card-b",
              "The tables are plain delimited text and need no special library."),
            tags$pre(class = "dl-code",
"# R\nd <- read.csv(\"AD_locus_evidence_20261002.csv\")\n\n# python\nimport pandas as pd\nd = pd.read_csv(\"AD_locus_variants_20261002.csv.gz\")")),
          div(class = "dl-card",
            div(class = "dl-card-t", "Evidence tiers"),
            div(class = "dl-tiers",
              lapply(names(TIER_DEFS), function(k)
                tagList(span(class = "k dc-mono", style = sprintf("color:%s", conf_pal[[k]]), k),
                        span(class = "v", TIER_DEFS[[k]]))))),
          div(class = "dl-card",
            div(class = "dl-card-t", "Terms and citation"),
            div(class = "dl-card-b",
              "If you publish with these tables, cite the consortium manuscript and the ",
              "release you downloaded."),
            div(class = "dl-cite",
              "The Alzheimer’s Disease Functional Genomics (FunGen-AD) Consortium. ",
              tags$b("Broad and deep dissection of Alzheimer’s disease genetics with "),
              tags$b("FunGen-xQTL."), " Data release 2026-10-02."))))),
    nav_panel(
      "Documentation", icon = icon("circle-info"),
    div(class = "doc",
      div(class = "doc-toc",
        div(class = "keylab", "Contents"),
        tags$a(href = "#doc-overview", "Overview"),
        tags$a(href = "#doc-study-design-and-data-sources", "Study design and data sources"),
        tags$a(href = "#doc-fine-mapping-molecular-and-disease-signals", "Fine-mapping molecular and disease signals"),
        tags$a(href = "#doc-linking-a-variant-to-a-gene", "Linking a variant to a gene"),
        tags$a(href = "#doc-how-an-ad-locus-is-defined", "How an AD locus is defined"),
        tags$a(href = "#doc-evidence-tiers", "Evidence tiers"),
        tags$a(href = "#doc-reading-the-table", "Reading the table"),
        tags$a(href = "#doc-caveats-worth-knowing", "Caveats worth knowing"),
        tags$a(href = "#doc-column-glossary", "Glossary"),
        tags$a(href = "#doc-sharing-a-view", "Sharing a view"),
        tags$a(href = "#doc-citation-and-data-availability", "Citation and data availability")),
      div(class = "doc-body",

    div(class = "doc-sec", id = "doc-overview", h2(class = "doc-h", "Overview"),
      p("For each Alzheimer's disease (AD) GWAS locus, which gene does the risk variant most ",
        "plausibly act through, and what molecular evidence supports that link? FunGen-xQTL ",
        "answers this by combining AD GWAS with multi-omic molecular quantitative trait loci (xQTL). ",
        "This explorer shows the locus-to-gene results. It summarizes the integration and does not ",
        "redistribute the underlying datasets."),
      p("Both GWAS and xQTL signals are represented as fine-mapped credible sets. A 95% credible set ",
        "is the smallest set of variants with at least 95% posterior probability of containing the ",
        "causal variant. Comparing a GWAS credible set with an xQTL credible set asks whether the two ",
        "signals are consistent with the same variant, a stronger test than positional overlap. ",
        "Compared with marginal testing, fine-mapping reduced candidate variants per gene 7.83-fold ",
        "and candidate eQTL per gene 205-fold, and credible-set xQTL annotations enrich AD GWAS ",
        "heritability 8-fold over marginal QTL."),
      p("Each row in this explorer is one variant and target gene pair at one AD locus, so a gene ",
        "appears once per supporting variant and can appear at several loci. This release covers ",
        format(n_loci_total, big.mark = ","), " candidate AD loci, of ",
        "which ", format(n_loci_qtl, big.mark = ","), " have xQTL evidence, and assigns ",
        format(n_t1_t5, big.mark = ","), " genes to tiers T1 through T5.")),
    div(class = "doc-sec", id = "doc-study-design-and-data-sources", h2(class = "doc-h", "Study design and data sources"),
      p("The Alzheimer's Disease Functional Genomics (FunGen-AD) Consortium produced the molecular ",
        "data, 77 datasets from 2,328 participants with brain and immune-cell profiles. Of these ",
        "participants, 1,769 have multi-omic measurements covering brain, blood and cardiovascular ",
        "tissue. The cohorts are ROSMAP, MSBB, Knight-ADRC, STARNET and MiGA, and MetaBrain serves ",
        "as an external bulk-tissue eQTL replication set."),
      p("The molecular modalities are expression (eQTL), splicing (sQTL), protein abundance (pQTL), ",
        "glycoprotein abundance (gpQTL), metabolite abundance (metQTL), DNA methylation (mQTL), ",
        "histone acetylation (haQTL) and chromatin accessibility (caQTL). Each modality is mapped in ",
        "its own contexts and then represented as credible sets, so modalities can be compared directly."),
      p("The Bellenguez, Jansen, Kunkle and Wightman AD GWAS were fine-mapped with ancestry-matched ",
        "linkage disequilibrium (LD) references, so the GWAS side of each comparison is also a ",
        "credible set rather than a list of marginal associations."),
      p(class = "note", "The full dataset catalog, with per-dataset cohort metadata and ",
        "Synapse download endpoints, is published at ",
        tags$a(href = "https://statfungen.github.io/xqtl-resources/", target = "_blank",
               rel = "noopener", "statfungen.github.io/xqtl-resources"),
        ". The catalog governs data access. Cite the manuscript for the locus-to-gene evidence ",
        "and the catalog for any dataset you download, since the two have separate access terms.")),

    div(class = "doc-sec", id = "doc-fine-mapping-molecular-and-disease-signals", h2(class = "doc-h", "Fine-mapping molecular and disease signals"),
      p("cis and trans QTL are fine-mapped with SuSiE-family models matched to the data type. SuSiE ",
        "is the default gene-level model. fSuSiE handles molecular traits measured as a curve or a ",
        "set of positions. scEEMS is used for low-power single-cell contexts. mvSuSiE jointly ",
        "fine-maps genes in the same TAD-informed region, so neighboring genes are modeled together."),
      p("The AD GWAS were fine-mapped with SuSiE-RSS, using the ADSP-derived LD reference and a model ",
        "robust to LD mismatch. Spurious credible sets arise when the reference panel and the study ",
        "sample differ in correlation structure, which is why the regional view uses this reference ",
        "and not a generic public LD panel."),
          p("The reference is a stochastic LD sketch built from 16,905 European-ancestry ADSP whole ",
            "genomes; multi-ancestry GWAS fine-mapping uses ancestry-specific sketches built the same ",
            "way. For xQTL, SuSiE uses individual-level genotype and phenotype data, SuSiE-RSS uses ",
            "summary statistics, mvSuSiE uses empirical priors for cell-type-specific expression and ",
            "protein QTL, and fSuSiE handles regional epigenomic traits that span an interval. ",
            "Colocalization across traits used ColocBoost, which reports each event as a 95% ",
            "colocalization confidence set."),
      p("A QTL is detectable when its 95% credible set is localized to a high-LD variant set, and ",
        "mappable when that set has three or fewer variants. Unless stated otherwise, xQTL here means ",
        "detectable credible-set QTL. The release contains 263,893 cis-xQTL credible sets targeting ",
        "305,234 molecular events across 18,195 genes.")),

    div(class = "doc-sec", id = "doc-linking-a-variant-to-a-gene", h2(class = "doc-h", "Linking a variant to a gene"),
      p("ColocBoost links molecular and disease signals and reports each colocalization event as a ",
        "95% colocalization confidence set (CoS). Each variant in a CoS has a variant colocalization ",
        "probability, the posterior probability that this variant is the one shared by the GWAS and ",
        "the molecular trait. This is stronger evidence than positional overlap, which two distinct ",
        "causal variants close together can also produce. The integration yields 37,673 ",
        "colocalization events covering 12,623 genes in 20 ROSMAP contexts."),
      p("Variants are also scored with cV2F, an AD consensus variant-to-function model. It extends the ",
        "original cV2F model with FunGen-xQTL features, namely cross-modality PIP, multi-trait ",
        "colocalization and conditional posterior effect estimates. The AUROC is 0.855 with the ",
        "original annotation and 0.918 with the FunGen-xQTL features added.")),

    div(class = "doc-sec", id = "doc-how-an-ad-locus-is-defined", h2(class = "doc-h", "How an AD locus is defined"),
      p("AD loci are 95% credible sets kept across 1,361 predefined European-ancestry LD blocks and ",
        "merged where credible sets overlap. A locus in this table is a credible set or colocalization ",
        "event with p-value 1e-5 or smaller."),
      p("Locus labels such as chr17_57 are numbered within this release and can change between ",
        "releases. The Locus detail tab therefore leads with genomic coordinates and shows the ",
        "label only as a footnote.")),
    div(class = "doc-sec", id = "doc-evidence-tiers", h2(class = "doc-h", "Evidence tiers"),
      p("Target genes are assigned to six evidence tiers. T1 to T4 require 95% posterior support ",
        "from AD-xQTL credible-set overlap or a colocalization confidence set (CoS), and higher tiers ",
        "add independent gene-level evidence. T5 holds lower-stringency localized evidence."),
      tags$table(class = "tier-table",
        tags$thead(tags$tr(tags$th("Tier"), tags$th("Requirement"), tags$th("Genes"))),
        tags$tbody(
          tags$tr(tags$td(tags$span(class="tier-chip t1","T1")),
                  tags$td("95% single-context credible-set overlap between an AD GWAS credible-set and an xQTL credible-set, plus supporting Mendelian randomization, causal TWAS or multi-context causal TWAS"),
                  tags$td(TIERN("T1"))),
          tags$tr(tags$td(tags$span(class="tier-chip t2","T2")),
                  tags$td("95% AD-xQTL CoS from multi-trait colocalization, plus supporting MR or causal TWAS"),
                  tags$td(TIERN("T2"))),
          tags$tr(tags$td(tags$span(class="tier-chip t3","T3")),
                  tags$td("95% credible-set overlap or CoS, plus supporting TWAS"),
                  tags$td(TIERN("T3"))),
          tags$tr(tags$td(tags$span(class="tier-chip t4","T4")),
                  tags$td("95% single-context credible-set overlap as stand-alone localized support"),
                  tags$td(TIERN("T4"))),
          tags$tr(tags$td(tags$span(class="tier-chip t5","T5")),
                  tags$td("Additional localized evidence from 95% CoS, or credible-set overlap computed at 50% or 70% coverage"),
                  tags$td(TIERN("T5"))))),
      p(class = "note",
        strong(sprintf("Genes with TWAS or MR evidence but no localized AD-xQTL support are T6 (%s genes here), counted apart from T1 to T5. ", TIERN("T6"))),
        "Gene-level evidence alone does not resolve the causal variant. A gene missing from the ",
        "tier table has not been ruled out; it lacks localized support.")),

    div(class = "doc-sec", id = "doc-reading-the-table", h2(class = "doc-h", "Reading the table"),
      p("Each row is one variant and target gene pair at one AD locus, so a gene appears once per ",
        "supporting variant. Locus is the AD locus that contains the variant, with a label numbered ",
        "within this release. Tier is the gene's best tier across the whole release, so a row can ",
        "show evidence at one tier while the gene's overall tier differs."),
      p("Sig is the GWAS significance of the variant, genome wide below 5e-8 and suggestive below ",
        "1e-6. Cell types lists the contexts (single-nucleus cell types or bulk tissue) with xQTL ",
        "evidence for the variant and gene. TWAS, MR and cTWAS are gene-level results and do not ",
        "depend on fine-mapping a variant. Trans shows whether the variant has trans-xQTL evidence ",
        "for genes outside the locus.")),
    div(class = "doc-sec", id = "doc-caveats-worth-knowing", h2(class = "doc-h", "Caveats worth knowing"),
      p("T6 marks genes with TWAS or Mendelian randomization evidence but no localized AD-xQTL ",
        "support. These genes stay in a separate gene-level set because such evidence does not ",
        "resolve the causal variant, and a T6 gene is not a negative result."),
      p("The table has two kinds of blank. No gene mapped marks a variant in the locus without a ",
        "target gene. A gene with no tier is linked to the locus but was not assigned a tier during ",
        "prioritization."),
      p("An empty gene-level cell means the annotation is not available for that variant or locus. ",
        "Read it as unknown, not as zero or tested and absent, and use the tier and context columns ",
        "to see what evidence exists.")),
    div(class = "doc-sec", id = "doc-column-glossary", h2(class = "doc-h", "Glossary"),
      accordion(open = FALSE,
        accordion_panel("cV2F score / rank",
          "Consensus variant-to-function score and its rank within the locus. Higher values mean the ",
          "variant is more likely to affect function."),
        accordion_panel("Maximum PIP",
          "Largest fine-mapping posterior inclusion probability (PIP) of the variant across methods ",
          "and contexts, shown with the method that produced it."),
        accordion_panel("Significance",
          "From the minimum GWAS p-value. Genome wide below 5e-8, suggestive below 1e-6, otherwise NS."),
        accordion_panel("Evidence tier (T1-T6)",
          "Best evidence tier of the gene. T1 to T5 have localized AD-xQTL support, strongest first. ",
          "T6 means TWAS or MR evidence without localized AD-xQTL support. No gene mapped marks ",
          "a variant without an xQTL target gene."),
        accordion_panel("TWAS / MR / cTWAS",
          "Whether the gene is significant in transcriptome-wide association, Mendelian ",
          "randomization or causal TWAS analysis. These are gene-level results and do not localize a variant."),
        accordion_panel("Context evidence",
          "A filled dot means the gene has xQTL evidence in that context. Contexts are Brain, ",
          "Excitatory, Inhibitory, Oligodendrocyte, OPC, Astrocyte, Microglia and Bulk Immune."),
        accordion_panel("Trans modalities",
          tags$dl(class = "gloss",
            tags$dt("snRNA"),
            tags$dd("Single-nucleus RNA. Expression is measured one nucleus at a time, so an effect can be assigned to a brain cell type such as neurons, microglia or astrocytes instead of being averaged over bulk tissue."),
            tags$dt("pQTL"),
            tags$dd("Protein abundance QTL. A variant that changes the amount of a protein."),
            tags$dt("gpQTL"),
            tags$dd("Glycoprotein abundance QTL. Sugar chains (glycans) on a protein affect how it folds, moves through the cell and is cleared. A gpQTL changes the glycosylation itself or the amount of glycoprotein. At TMEM106B site N183, the glycosylation change tracks the protein effect and RNA level does not."),
            tags$dt("metQTL"),
            tags$dd("Metabolite abundance QTL. A variant that changes the level of a small molecule (metabolite). In this release, metQTL evidence comes from the transMap colocalization analysis."),
            ),
          p(class = "note", "Each row also gives the number of genes and the contexts of the association.")),
        accordion_panel("Variant resources",
          "Each variant links to gnomAD (allele frequency, constraint, population genetics), ",
          "Open Targets Genetics, dbSNP and the UCSC browser at that position.")),
        accordion_panel("Glossary of statistical terms",
          tags$dl(class = "gloss",
            tags$dt("Fine-mapping"),
            tags$dd("Narrows an association from the many variants inherited together to the few most likely to be causal."),
            tags$dt("Credible set"),
            tags$dd("The short list that fine-mapping produces. A 95% credible set is the smallest set of variants with a 95% probability of containing the causal variant. Smaller sets are more precise."),
            tags$dt("PIP (posterior inclusion probability)"),
            tags$dd("Probability, from 0 to 1, that a given variant is causal. A variant with PIP 0.8 carries most of the probability in its credible set."),
            tags$dt("Colocalization"),
            tags$dd("Tests whether a disease association and a molecular association at the same place come from the same variant, or from two different variants that lie close together."),
            tags$dt("CoS (colocalization confidence set)"),
            tags$dd("The credible-set form of colocalization, the set of variants likely to contain the causal variant shared by both traits."),
            tags$dt("TWAS"),
            tags$dd("Transcriptome-wide association study. It predicts gene expression from genotype, then tests which predicted genes associate with disease risk. It points to genes, not variants."),
            tags$dt("MR (Mendelian randomization)"),
            tags$dd("Uses genetic variants as natural experiments to test whether changing a gene product would change disease risk, beyond simple correlation."),
            tags$dt("LD (linkage disequilibrium)"),
            tags$dd("Nearby variants tend to be inherited together, so they show similar association with disease. LD is why fine-mapping is needed to find the causal variant."),
            tags$dt("cis and trans"),
            tags$dd("A cis effect is a variant acting on a nearby gene. A trans effect is a variant acting on a distant gene, often on another chromosome.")))),

    div(class = "doc-sec", id = "doc-sharing-a-view", h2(class = "doc-h", "Sharing a view"),
      p("Filters on the Loci tab are saved in the web address. Copy the address, or use Copy link ",
        "to this view in the filter panel, and the link opens the same gene, locus, tiers and thresholds."),
      tags$table(class = "tier-table",
        tags$thead(tags$tr(tags$th("Parameter"), tags$th("Meaning"), tags$th("Example"))),
        tags$tbody(
          tags$tr(tags$td("q"), tags$td("gene, rsID or coordinates"), tags$td("q=BIN1")),
            tags$tr(tags$td("locus"), tags$td("AD loci to keep, as coordinates"), tags$td("locus=chr2:127088929-127118906")),
          tags$tr(tags$td("tier"), tags$td("evidence tiers to keep"), tags$td("tier=T1,T2")),
          tags$tr(tags$td("sig"), tags$td("GWAS significance classes"), tags$td("sig=genome%20wide")),
          tags$tr(tags$td("ct"), tags$td("context evidence columns"), tags$td("ct=ct_Mic_xQTL")),
          tags$tr(tags$td("minlp"), tags$td("minimum -log10(p)"), tags$td("minlp=8")),
          tags$tr(tags$td("transonly"), tags$td("trans evidence only"), tags$td("transonly=1")),
            tags$tr(tags$td("gene"), tags$td("locus shown on Locus detail, by gene symbol"), tags$td("gene=BIN1")),
            tags$tr(tags$td("variant"), tags$td("locus shown on Locus detail, by rsID"), tags$td("variant=rs429358")),
            tags$tr(tags$td("region"), tags$td("locus shown on Locus detail, by coordinates"), tags$td("region=chr2:127088929-127118906")),
          tags$tr(tags$td("tab"), tags$td("which tab opens"), tags$td("tab=Loci")))),
      p(class = "note",
        "Links omit default settings, so they stay short. They use genomic coordinates instead of ",
        "locus labels, so a saved view still works after loci are redrawn. Older links with a locus ",
        "label still open and are converted to coordinates.")),


    div(class = "doc-sec", id = "doc-citation-and-data-availability", h2(class = "doc-h", "Citation and data availability"),
      p("The Alzheimer's Disease Functional Genomics (FunGen-AD) Consortium. ",
        em("Broad and deep dissection of Alzheimer's disease genetics with FunGen-xQTL."),
        " Correspondence: Gao Wang and Philip L. De Jager."),
        p("Supported by NIH grant U01AG072572."),
      p(class = "note",
        "The AD loci and xQTL integration pipeline produces the tier assignments (T1 to T6, where ",
        "T6 means TWAS or MR evidence only) and the locus and variant evidence. The project README ",
        "gives the steps to reproduce them.")),
        div(class = "doc-sec", id = "doc-cell-types-and-assays",
          h2(class = "doc-h", "Contexts and modalities"),
          p("Number of gene records in this release for each context and molecular modality."),
          uiOutput("doc_assays"),
          div(class = "dc-exp", style = "margin-top:8px",
              downloadLink("dl_doc_assays", "CSV", class = "dc-btn"))),
      )))

)

# ---- server ----------------------------------------------------------------
