server <- function(input, output, session) {

  # ---- deep links --------------------------------------------------------
  # The filter state is mirrored into the query string and read back on load,
  # so any view can be bookmarked, pasted into a message, or cited in a paper
  # and it reopens exactly as it was. Multi-selects are comma joined; an
  # empty multi-select is written as "none" so "nothing selected" survives a
  # round trip instead of silently reverting to the default.
  qs_multi   <- c("locus", "sig", "tier", "ct", "tmod")
  qs_default <- list(locus = character(0),
                     sig   = c("genome wide", "suggestive", "ns", NO_P),
                     tier = c(paste0("T", 1:6), NO_TIER, NO_GENE),
                     ct    = character(0),
                     tmod  = character(0))
  restored <- reactiveVal(FALSE)

  # Single source of truth for which locus is being inspected. Search, URL,
  # prev/next and the dropdown are all writers; every output is a reader.
  selected_locus <- reactiveVal(loci[1])

  observeEvent(input$locus_pick, {
    if (!identical(input$locus_pick, selected_locus())) selected_locus(input$locus_pick)
  }, ignoreInit = TRUE)

  observeEvent(selected_locus(), {
    if (!identical(input$locus_pick, selected_locus()))
      updateSelectInput(session, "locus_pick", selected = selected_locus())
  }, ignoreNULL = TRUE)


  observe({
    if (!restored()) return()
    parts <- list()
    for (k in qs_multi) {
      v <- input[[k]]; d <- qs_default[[k]]
      if (length(v) == length(d) && setequal(v, d)) next          # at default
      parts[[k]] <- if (!length(v)) "none" else paste(if (k == "locus") loci_to_regions(v) else v, collapse = ",")
    }
    if (!is.null(input$q) && nzchar(input$q))            parts$q <- input$q
    if (!is.null(input$gp_gene) && nzchar(input$gp_gene)) parts$gene <- input$gp_gene
    if (!is.null(input$minlp) && input$minlp > 0)        parts$minlp <- input$minlp
    if (isTRUE(input$transonly))                         parts$transonly <- "1"
    # locus is serialised separately below, as stable keys
    if (!is.null(input$nav) && !identical(input$nav, "Home"))
      parts$tab <- input$nav
    qs <- paste(vapply(names(parts),
                       function(k) paste0(k, "=", URLencode(as.character(parts[[k]]), reserved = TRUE)),
                       character(1)), collapse = "&")
    .lq <- if (!is.null(selected_locus()) && nzchar(selected_locus()))
             sub("^\\?", "", locus_to_query(selected_locus())) else ""
    .all <- c(.lq, qs); .all <- paste(.all[nzchar(.all)], collapse = "&")
    updateQueryString(if (nzchar(.all)) paste0("?", .all) else "?", mode = "replace")
  })

  # ---- Home search ----------------------------------------------------------
  updateSelectizeInput(session, "hero_search",
                       choices  = search_index[, c("value", "label", "gid")],
                       selected = "", server = TRUE)

  observeEvent(list(input$hero_search, input$hero_go), {
    req(nzchar(input$hero_search))
    .tk <- batch_tokens(input$hero_search)
    if (length(.tk) > 1L) {
      updateTextAreaInput(session, "bt_in", value = paste(.tk, collapse = "
"))
      bt_tok(.tk)
      updateSelectizeInput(session, "hero_search", selected = "")
      nav_select("nav", "Batch", session = session)
      return(invisible())
    }
    i <- match(input$hero_search, search_index$value)
    if (!is.na(i)) {
      selected_locus(search_index$ADlocus[i])
      updateTextInput(session, "q", value = search_index$term[i])
      nav_select("nav", "Locus detail", session = session)
      return(invisible())
    }
    .r <- resolve_locus(list(region = input$hero_search))
    if (!is.null(.r)) {
      selected_locus(.r$locus)
      nav_select("nav", "Locus detail", session = session)
    } else {
      showNotification("No gene, variant or region matched that.",
                       type = "warning", duration = 5)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$gp_clear, {
    updateSelectizeInput(session, "gp_gene", choices = genes, selected = "", server = TRUE)
  })

  observeEvent(input$go_explore, nav_select("nav", "Loci", session = session))
  observeEvent(input$go_methods, nav_select("nav", "Documentation", session = session))

  fmt <- function(x) format(x, big.mark = ",")

  # ----- dashboard KPIs & plots (full data) -----
  output$kpi_pct <- renderText({
    n <- nrow(locus_registry)
    k <- length(unique(dat$ADlocus[!is.na(dat$gene) & nzchar(as.character(dat$gene))]))
    sprintf("%.0f%%", 100 * k / n)
  })
  output$kpi_rows  <- renderText(fmt(nrow(dat)))
  output$kpi_loci  <- renderText(fmt(length(loci)))
  output$kpi_genes <- renderText(fmt(length(genes)))
  output$kpi_trans <- renderText(fmt(sum(dat$has_trans, na.rm = TRUE)))
  output$kpi_tier  <- renderText({ m <- median(as.integer(dat$top_confidence[!dat$top_confidence %in% c(NO_GENE, NO_TIER)]), na.rm = TRUE)
                                    if (is.na(m)) "—" else paste0("T", round(m)) })

  # ---- Release overview figures --------------------------------------------
  .locus_summary <- local({
    d <- dat[!is.na(dat$pos) & !is.na(dat$chr), , drop = FALSE]
    d$.chr <- sub("^chr", "", as.character(d$chr))
    ordc <- c(as.character(1:22), "X")
    d <- d[d$.chr %in% ordc, , drop = FALSE]
    tnum <- suppressWarnings(as.integer(sub("^T", "", as.character(d$top_confidence))))
    d$.tr <- ifelse(is.na(tnum), 99L, tnum)
    s <- d %>%
      dplyr::group_by(ADlocus, .chr) %>%
      dplyr::summarise(
        pos = stats::median(pos, na.rm = TRUE),
        lp  = suppressWarnings(max(log10pval, na.rm = TRUE)),
        ng  = dplyr::n_distinct(gene[!is.na(gene) & gene != ""]),
        tr  = min(.tr, na.rm = TRUE),
        .groups = "drop")
    s$lp[!is.finite(s$lp)] <- NA_real_
    s$tier <- ifelse(s$tr <= 6, paste0("T", s$tr), NO_TIER)
    s$.chr <- factor(s$.chr, levels = ordc)
    span <- d %>% dplyr::group_by(.chr) %>%
      dplyr::summarise(m = max(pos, na.rm = TRUE), .groups = "drop")
    span$.chr <- factor(span$.chr, levels = ordc)
    span <- span[order(span$.chr), ]
    span$off <- c(0, cumsum(as.numeric(span$m))[-nrow(span)])
    s$x <- as.numeric(s$pos) + span$off[match(as.character(s$.chr), as.character(span$.chr))]
    list(loci = s, mid = data.frame(chr = as.character(span$.chr),
                                    mid = span$off + as.numeric(span$m) / 2))
  })

  land_df <- reactive({
    o <- chr_offsets
    r <- locus_registry
    lt <- locus_table[match(r$ADlocus, locus_table$ADlocus), , drop = FALSE]
    d <- data.frame(
      ADlocus = as.character(r$ADlocus),
      chr = sub("^chr", "", as.character(r$chr)),
      pos = (as.numeric(r$start) + as.numeric(r$end)) / 2,
      y   = suppressWarnings(as.numeric(
              dat$log10pval[match(as.character(r$ADlocus), dat$ADlocus)])),
      tier = lt$best, ng = lt$n_genes, gene = lt$top_gene,
      stringsAsFactors = FALSE)
    d <- d[d$chr %in% CHR_ORD & !is.na(d$y), , drop = FALSE]
    d$x <- d$pos + o$off[match(d$chr, o$chr)]
    d$tier <- factor(ifelse(is.na(d$tier), "No tier", d$tier),
                     levels = c(TIER_SEQ, "No tier"))
    d
  })

  fig_p_landscape <- reactive({
    o <- chr_offsets
    d <- land_df(); req(nrow(d) > 0)
    pal <- c(conf_pal[TIER_SEQ], "No tier" = "#cfd6db")
    names(pal) <- c(TIER_SEQ, "No tier")
    bands <- data.frame(chr = o$chr, x1 = o$off, x2 = o$off + o$len,
                        mid = o$off + o$len / 2,
                        alt = seq_len(nrow(o)) %% 2 == 0)
    ymax <- max(d$y, na.rm = TRUE) * 1.08
    dk <- FALSE
    inkc <- if (dk) "#e8eef2" else "#233947"
    grdc <- if (dk) "#22313e" else "#edf0f2"
    lab <- d[order(-d$y), , drop = FALSE]
    lab <- lab[!is.na(lab$ADlocus), , drop = FALSE]
    lab <- lab[!duplicated(lab$ADlocus), , drop = FALSE]
    lab <- head(lab, 6)
    ggplot() +
      geom_rect(data = bands[bands$alt, ],
                aes(xmin = x1, xmax = x2, ymin = 0, ymax = ymax),
                fill = "#f4f6f7", colour = NA) +
      geom_hline(yintercept = pretty(c(0, ymax), 4), colour = grdc,
                 linewidth = 0.4) +
      geom_point(data = d, aes(x, y, colour = tier, size = ng), alpha = 0.95) +
      geom_text(data = lab, aes(x, y, label = ADlocus), vjust = -0.9, size = 3,
                colour = inkc, fontface = "bold", check_overlap = TRUE) +
      scale_x_continuous(breaks = bands$mid, labels = bands$chr,
                         expand = c(0.004, 0)) +
      scale_y_continuous(breaks = pretty(c(0, ymax), 4), expand = c(0, 0),
                         limits = c(0, ymax)) +
      scale_colour_manual(values = pal, name = "Best tier", drop = FALSE) +
      scale_size_continuous(range = c(1.6, 6.5), name = "Genes at locus") +
      guides(colour = "none", size = guide_legend(order = 1)) +
      labs(x = NULL, y = expression(-log[10] * P ~ "(Lead variant)")) + g_theme +
      theme(panel.background = element_rect(fill = "#ffffff", colour = NA),
            plot.background = element_rect(fill = "#ffffff", colour = NA),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8.5, colour = "#5c6b75"),
            axis.title.y = element_text(size = 9, colour = "#5c6b75"),
            legend.position = "bottom", legend.box = "horizontal",
            axis.line.x = element_line(colour = "#dde2e6", linewidth = 0.5)) +
      dk_theme(FALSE)
  })
  output$p_landscape <- renderPlot({ fig_p_landscape() })

  fig_p_ctxtier <- reactive({
    g <- dat[!is.na(dat$gene) & dat$gene != "", , drop = FALSE]
    tnum <- suppressWarnings(as.integer(sub("^T", "", as.character(g$top_confidence))))
    g$.tier <- ifelse(is.na(tnum), NA_character_, paste0("T", tnum))
    g <- g[!is.na(g$.tier), , drop = FALSE]
    validate(need(nrow(g) > 0, "No tiered genes to plot."))
    rows <- do.call(rbind, lapply(seq_along(ct_cols), function(i) {
      cc <- ct_cols[i]
      if (!cc %in% names(g)) return(NULL)
      sub <- g[which(g[[cc]] %in% TRUE), , drop = FALSE]
      if (!nrow(sub)) return(NULL)
      tb <- sub %>% dplyr::group_by(.tier) %>%
        dplyr::summarise(n = dplyr::n_distinct(gene), .groups = "drop")
      data.frame(ctx = ct_labels[i], tier = tb$.tier, n = tb$n, stringsAsFactors = FALSE)
    }))
    validate(need(!is.null(rows) && nrow(rows) > 0, "No cell-type evidence to plot."))
    rows$tier <- factor(rows$tier, levels = paste0("T", 1:6))
    ordctx <- rows %>% dplyr::group_by(ctx) %>%
      dplyr::summarise(tot = sum(n), .groups = "drop") %>% dplyr::arrange(tot)
    rows$ctx <- factor(rows$ctx, levels = ordctx$ctx)
    ggplot(rows, aes(x = tier, y = ctx, fill = n)) +
      geom_tile(colour = PAPER, linewidth = 1.4) +
      geom_text(aes(label = n), size = 3.2, colour = INK) +
      scale_fill_gradient(low = "#eaf1fb", high = "#2a78d6", name = "Genes") +
      labs(x = NULL, y = NULL) + g_theme +
      guides(fill = guide_colourbar(barheight = unit(5, "pt"), barwidth = unit(110, "pt"),
                                    title.position = "left", title.vjust = 1)) +
      theme(legend.position = "bottom", panel.grid = element_blank(),
            legend.margin = margin(4, 0, 0, 0),
            axis.text.y = element_text(colour = INK, size = 10),
            axis.text.x = element_text(colour = INK, size = 10))
  })
  output$p_ctxtier <- renderPlot({ fig_p_ctxtier() })

  output$p_tier <- renderPlot({
    d <- dat %>% filter(!is.na(gene), gene != "") %>%
      distinct(gene, top_confidence) %>%
      count(top_confidence, name = "n") %>%
      filter(!is.na(top_confidence))
    rowbar(data.frame(lab = as.character(d$top_confidence), n = d$n),
           pal = conf_pal)
  })
output$p_loci <- renderPlot({
  d <- dat %>% group_by(ADlocus) %>%
    summarise(n = n_distinct(gene[!is.na(gene) & gene != ""]), .groups = "drop") %>%
    arrange(desc(n)) %>% head(5) %>% rename(lab = ADlocus)
  rowbar(d)
})
output$p_ct <- renderPlot({
    .fd <- dat
  d <- tibble(lab = ct_labels,
              n = vapply(ct_cols, function(c) sum(.fd[[c]], na.rm = TRUE), numeric(1))) %>%
    arrange(desc(n)) %>% head(5)
  rowbar(d)
})
output$p_sig <- renderPlot({
  d <- dat %>% count(significance, name = "n") %>%
    filter(!is.na(significance)) %>% arrange(desc(n)) %>%
    rename(lab = significance)
  rowbar(d, pal = SIG_COL)
})
output$p_trans <- renderPlot({
    .fd <- dat
  d <- tibble(lab = names(trans_types),
              n = vapply(trans_types, function(c) sum(.fd[[c]] > 0, na.rm = TRUE), numeric(1))) %>%
    arrange(desc(n)) %>% head(5)
  rowbar(d)
})

  # ----- Explore filtering -----
  filtered <- reactive({
    d <- dat
    if (nzchar(input$q)) {
      qq <- str_trim(input$q)
      # Accept the coordinate forms people actually paste, not just gene/rsID:
      #   chr8:27300000-27400000   a range
      #   chr8:27369273            a single position (matched within 1 kb)
      #   chr8:27369273:A:C        a full variant ID
      # Anything else falls through to the original substring match.
      m_range <- str_match(qq, "^(?:chr)?([0-9XYxy]+)[:[:space:]]+([0-9,]+)[[:space:]]*[-_][[:space:]]*([0-9,]+)$")
      m_pos   <- str_match(qq, "^(?:chr)?([0-9XYxy]+)[:[:space:]]+([0-9,]+)(?::[ACGTacgt]+:[ACGTacgt]+)?$")
      num <- function(x) as.numeric(gsub(",", "", x))
      if (!is.na(m_range[1, 1])) {
        d <- d[na_fill(as.character(d$chr), "") == m_range[1, 2] &
               na_fill(d$pos, -1) >= num(m_range[1, 3]) &
               na_fill(d$pos, -1) <= num(m_range[1, 4]), ]
      } else if (!is.na(m_pos[1, 1])) {
        d <- d[na_fill(as.character(d$chr), "") == m_pos[1, 2] &
               abs(na_fill(d$pos, -1e12) - num(m_pos[1, 3])) <= 1000, ]
      } else {
        pat <- regex(qq, ignore_case = TRUE)
        d <- d[str_detect(na_fill(d$gene, ""), pat) | str_detect(na_fill(d$rsid, ""), pat), ]
      }
    }
    if (length(input$locus)) d <- d[d$ADlocus %in% input$locus, ]
    if (length(input$sig))  d <- d[as.character(d$significance) %in% input$sig, ]
    if (length(input$tier)) d <- d[as.character(d$top_confidence) %in% input$tier, ]
    if (length(input$ct)) for (c in input$ct) d <- d[which(d[[c]]), ]
    d <- d[na_fill(d$log10pval, 0) >= input$minlp, ]
    if (isTRUE(input$transonly)) d <- d[na_fill(d$has_trans, FALSE), ]
    d
  })

  output$e_rows  <- renderText(fmt(nrow(filtered())))
  output$e_loci  <- renderText(fmt(length(unique(filtered()$ADlocus))))
  output$e_genes <- renderText(fmt(length(unique(na.omit(filtered()$gene)))))
  output$e_trans <- renderText(fmt(sum(filtered()$has_trans, na.rm = TRUE)))

  disp <- reactive({
    d <- filtered()
    ## UX: order by tier (T1 strongest -> T6), unmapped rows last,
    ## then by GWAS significance within tier.
    .tord <- suppressWarnings(as.integer(sub('^T', '', as.character(d$top_confidence))))
    .tord[is.na(.tord)] <- 99L
    d <- d[order(.tord, -na_fill(d$log10pval, 0)), , drop = FALSE]
    tibble(
      Locus = d$ADlocus, Gene = d$gene, rsID = d$rsid,
      Pos = ifelse(is.na(d$chr), "—", paste0("chr", d$chr, ":", fmt(d$pos))),
      Sig = sig_badge(d$significance),
      `-log10p` = round(d$log10pval, 2),
      cV2F = round(d$cv2f_score, 3),
      Tier = conf_badge(d$top_confidence),
      `#Ctx` = d$n_contexts,
      `xQTL PIP` = round(d$xqtl_max_inclusion, 3),
      `TWAS z` = round(d$max_twas_z, 2),
      TWAS = ev_mark(d$twas_sig),
      MR   = ev_mark(d$mr_sig),
      cTWAS = ev_mark(d$ctwas_sig),
      Trans = ifelse(d$has_trans, paste0('<span style="color:#0b6e4f">\u25cf</span> ', na_fill(d$trans_n_genes, 0)),
                     '<span style="color:#e2e8f0">\u25cb</span>'),
      `Cell types` = d$ct_dots
    )
  })

  # ---- Gene-level summary (level 1 of the two-level table) -------------------
  gene_summary <- reactive({
    d <- filtered()
    d <- d[!is.na(d$gene) & d$gene != "", , drop = FALSE]
    if (!nrow(d)) return(NULL)
    .tier_num <- suppressWarnings(as.integer(sub("^T", "", as.character(d$top_confidence))))
    d$.tr <- ifelse(is.na(.tier_num), 99L, .tier_num)
    fin <- function(x) if (is.finite(x)) x else NA_real_
    d %>%
      dplyr::group_by(Gene = gene) %>%
      dplyr::summarise(
        .tr       = min(.tr, na.rm = TRUE),
        `Max xQTL PIP` = fin(suppressWarnings(max(xqtl_max_inclusion, na.rm = TRUE))),
        Contexts  = fin(suppressWarnings(max(n_contexts, na.rm = TRUE))),
        Loci      = dplyr::n_distinct(ADlocus),
        Rows      = dplyr::n(),
        TWAS      = dplyr::first(twas_sig),
        MR        = dplyr::first(mr_sig),
        cTWAS     = dplyr::first(ctwas_sig),
        .groups   = "drop") %>%
      dplyr::mutate(Tier = ifelse(.tr <= 6, paste0("T", .tr), NA_character_),
                    `Max xQTL PIP` = round(`Max xQTL PIP`, 3)) %>%
      dplyr::arrange(.tr, dplyr::desc(dplyr::coalesce(`Max xQTL PIP`, -1))) %>%
      as.data.frame()
  })

  output$gene_tbl <- renderDT({
    g <- gene_summary()
    validate(need(!is.null(g) && nrow(g) > 0, "No genes match the current filters."))
    datatable(
      tibble(Gene = g$Gene, Tier = conf_badge(g$Tier), `Max xQTL PIP` = g$`Max xQTL PIP`,
             Contexts = g$Contexts, Loci = g$Loci, `Variant rows` = g$Rows,
             TWAS = ev_mark(g$TWAS), MR = ev_mark(g$MR), cTWAS = ev_mark(g$cTWAS)),
      escape = FALSE, selection = "single", rownames = FALSE,
      options = list(dom = "lrtip", pageLength = 25, scrollX = TRUE, headerCallback = hdr_js,
                     columnDefs = link_defs(0)),
      class = "display compact hover stripe")
  }, server = TRUE)

  observeEvent(input$gene_tbl_rows_selected, {
    i <- input$gene_tbl_rows_selected
    req(length(i) == 1)
    g <- gene_summary()
    req(!is.null(g), i <= nrow(g))
    updateTextInput(session, "q", value = g$Gene[i])
    updateRadioButtons(session, "view_mode", selected = "variant")
  })

  output$tbl <- renderDT({
    datatable(disp(), escape = FALSE, selection = "single", rownames = FALSE,
      options = list(dom = "lrtip", pageLength = 25, scrollX = TRUE, headerCallback = hdr_js,
        columnDefs = c(list(list(orderable = FALSE,
          targets = which(names(disp()) %in% c("Cell types","TWAS","MR","cTWAS","Sig","Trans")) - 1)), link_defs(which(names(disp()) == "Gene") - 1, which(names(disp()) == "rsID") - 1))),
      class = "display compact hover stripe")
  }, server = TRUE)

  output$detail <- renderUI({
    s <- input$tbl_rows_selected
    if (is.null(s)) return(NULL)
    r <- filtered()[s, ]
    tr <- function(lbl, modg, modc) if (!is.na(modg) && nzchar(as.character(modg)))
      tags$li(strong(paste0(lbl, ": ")), as.character(modg),
              if (!is.na(modc) && nzchar(as.character(modc))) tags$em(paste0(" (", modc, ")")))
    card(card_header(sprintf("%s · %s · %s", r$gene, r$ADlocus, r$rsid)),
      layout_columns(col_widths = c(6,6),
        tags$div(
          tags$p(strong("Variant: "), r$variant_ID, " · ", strong("Effect allele: "), r$effect_allele),
          tags$p(strong("Min p: "), signif(r$min_pval,3), " · ", strong("Significance: "), as.character(r$significance)),
          tags$p(strong("cV2F: "), round(r$cv2f_score,3), " (rank ", r$cv2f_rank, ") · ",
                 strong("GWAS PIP: "), round(r$max_inclusion,3), " [", pretty_list(r$max_inclusion_method), "]"),
          tags$p(strong("Tier: "), as.character(r$top_confidence), " · ",
                 strong("# contexts: "), r$n_contexts, " · ",
                 strong("Dist TSS: "), ifelse(is.na(r$dist_tss),"—", paste0(round(r$dist_tss/1000,1)," kb"))),
          tags$p(strong("GWAS: "), pretty_list(r$gwas_assoc)),
          tags$p(strong("Ordered contexts: "), tags$small(pretty_list(ds_relabel(r$ordered_contexts)))),
          ext_links(r$variant_ID, r$rsid),
          gene_links(r$gene, unname(gid_map[r$gene]))
        ),
        tags$div(
          tags$p(strong("Trans evidence:")),
          tags$ul(
            tr("Trans genes", r$trans_genes, r$trans_contexts),
            tr("snRNA",  r$trans_snRNA_genes,  r$trans_snRNA_contexts),
            tr("pQTL",   r$trans_pQTL_genes,   r$trans_pQTL_contexts),
            tr("gpQTL",  r$trans_gpQTL_genes,  r$trans_gpQTL_contexts),
            tr("Hotspot",r$trans_hotspot_programs, r$trans_hotspot_contexts)
          )
        )))
  })

  output$dl <- downloadHandler(
    filename = function() paste0("AD_xQTL_filtered_", Sys.Date(), ".csv"),
    content  = function(f) .prov_csv(filtered() %>% select(-ct_dots), f, "Loci tab, filtered evidence rows"))

  # ----- Locus detail -----
  # Regional view: every variant at the selected locus by position. The app had no
# positional plot at all - only summary bars - so a locus could not be inspected
# the way a genome browser would show it.
# Locus-to-gene dot table, after Figure 5 of the manuscript: rows are the
# genes assigned to this locus, columns are molecular contexts, fill is the
# gene's evidence tier and point size is how many variant rows support it.
# This is also what retires the colour-only cell-type dots column, whose
# eight hues were not separable under deutan simulation.
  locus_win <- reactive({
    loc <- selected_locus(); req(loc)
    i <- match(loc, locus_registry$ADlocus); req(!is.na(i))
    list(chr = as.character(locus_registry$chr[i]),
         lo  = as.numeric(locus_registry$start[i]) - 50000,
         hi  = as.numeric(locus_registry$end[i])   + 50000)
  })

  locus_ev_all <- reactive({
    loc <- selected_locus(); req(loc)
    e <- parse_ctx_tokens(dat[dat$ADlocus == loc, , drop = FALSE])
    if (is.null(e) || !nrow(e)) return(NULL)
    if (!is.null(gene_pos)) {
      j <- match(e$gene, gene_pos$gene)
      e$gx   <- (as.numeric(gene_pos$start[j]) + as.numeric(gene_pos$end[j])) / 2
      e$gchr <- as.character(gene_pos$chr[j])
    } else { e$gx <- NA_real_; e$gchr <- NA_character_ }
    w <- locus_win()
    e$inwin <- !is.na(e$gx) & !is.na(e$gchr) &
      e$gchr == sub("^chr", "", w$chr) & e$gx >= w$lo & e$gx <= w$hi
    e
  })

  locus_ev <- reactive({
    e <- locus_ev_all(); if (is.null(e)) return(NULL)
    if (!is.null(input$rv_ctx) && length(input$rv_ctx)) e <- e[e$ctx %in% input$rv_ctx, , drop = FALSE]
    if (!is.null(input$rv_mod) && length(input$rv_mod)) e <- e[e$mod %in% input$rv_mod, , drop = FALSE]
    tt <- if (is.null(input$rv_tier)) 6L else as.integer(input$rv_tier)
    e[e$tier <= tt, , drop = FALSE]
  })

  observeEvent(locus_ev_all(), {
    e  <- locus_ev_all()
    cx <- if (is.null(e)) character(0) else sort(unique(e$ctx))
    md <- if (is.null(e)) character(0) else sort(unique(e$mod))
    nm <- ifelse(cx %in% names(CTX_FULL), paste0(CTX_FULL[cx], " (", cx, ")"), cx)
    updateSelectizeInput(session, "rv_ctx", choices = stats::setNames(cx, nm),
                         selected = intersect(input$rv_ctx, cx), server = FALSE)
    updateSelectizeInput(session, "rv_mod", choices = md,
                         selected = intersect(input$rv_mod, md), server = FALSE)
  }, ignoreNULL = FALSE)

  observeEvent(input$rv_reset, {
    updateSelectizeInput(session, "rv_ctx", selected = character(0))
    updateSelectizeInput(session, "rv_mod", selected = character(0))
    updateSliderInput(session, "rv_tier", value = 6)
  })

  observeEvent(input$go_trans, nav_select("nav", "Trans", session = session))
  observeEvent(input$go_batch, nav_select("nav", "Batch", session = session))

  output$rv_stat <- renderUI({
    a <- locus_ev_all()
    if (is.null(a)) return(div(class = "rv-stat",
      "No cell-type or modality evidence is recorded for genes at this locus."))
    e   <- locus_ev()
    ev  <- if (is.null(e) || !nrow(e)) a[0, , drop = FALSE] else e[e$inwin, , drop = FALSE]
    tot <- sum(a$inwin); off <- sum(!a$inwin)
    if (tot == 0) return(HTML(paste0('<div class="rv-stat">No evidence marks fall inside the plotted window', if (off > 0) paste0(' &#183; all ', format(off, big.mark = ","), ' marks for genes at this locus sit outside it and are listed in the table below') else '', '</div>')))
    HTML(paste0('<div class="rv-stat"><b>', format(nrow(ev), big.mark = ","), '</b> of ',
      format(tot, big.mark = ","), ' evidence marks shown &#183; ',
      length(unique(ev$gene)), ' genes &#183; ',
      length(unique(ev$ctx)), ' cell types &#183; ',
      length(unique(ev$mod)), ' modalities',
      if (off > 0) paste0(' &#183; ', off,
        ' further marks belong to genes outside this window and stay in the table below')
      else '', '</div>'))
  })

  rv_keep <- function(e) {
    k <- rep(TRUE, nrow(e))
    if (!is.null(input$rv_ctx) && length(input$rv_ctx)) k <- k & e$ctx %in% input$rv_ctx
      qv <- trimws(if (is.null(input$rv_q)) "" else input$rv_q)
      if (nzchar(qv)) {
        pt <- stringr::regex(qv, ignore_case = TRUE)
        hit <- rep(FALSE, length(k))
        if (!is.null(e$gene)) hit <- hit | stringr::str_detect(na_fill(as.character(e$gene), ""), pt)
        if (!is.null(e$rsid)) hit <- hit | stringr::str_detect(na_fill(as.character(e$rsid), ""), pt)
        k <- k & hit
      }
    if (!is.null(input$rv_mod) && length(input$rv_mod)) k <- k & e$mod %in% input$rv_mod
    tt <- if (is.null(input$rv_tier)) 6L else as.integer(input$rv_tier)
    k & e$tier <= tt
  }

  region_fig <- reactive({
    loc <- selected_locus(); req(loc)
    w <- locus_win()
    v <- dat[dat$ADlocus == loc & !is.na(dat$pos) & !is.na(dat$log10pval),
             c("pos", "log10pval", "rsid"), drop = FALSE]
    v <- unique(as.data.frame(v))
    v <- v[v$pos >= w$lo & v$pos <= w$hi, , drop = FALSE]
    vp <- unique(as.numeric(dat$pos[dat$ADlocus == loc & !is.na(dat$pos)]))
    vp <- vp[vp >= w$lo & vp <= w$hi]
    a  <- locus_ev_all()
    ev <- if (is.null(a)) NULL else a[a$inwin, , drop = FALSE]
    if (!is.null(ev) && !nrow(ev)) ev <- NULL
    H <- if (nrow(v)) max(c(v$log10pval, GWS * 1.05), na.rm = TRUE) else 1
    if (!is.finite(H) || H <= 0) H <- 1
    step <- H * 0.075; base <- -H * 0.17
    cats <- character(0)
    if (!is.null(ev)) {
      ev$on   <- rv_keep(ev)
      ev$cell <- ifelse(ev$ctx %in% names(CTX_FULL), CTX_FULL[ev$ctx], ev$ctx)
      ev$cat  <- paste0(ev$cell, "  ", ev$mod)
      o <- order(match(ev$ctx, CTX_ORD), match(ev$mod, MOD_ORD))
      cats <- unique(ev$cat[o])
      ev$y <- base - (match(ev$cat, cats) - 1) * step
      ev$ctx <- factor(ev$ctx, levels = intersect(CTX_ORD, unique(ev$ctx)))
      ev$mod <- factor(ev$mod, levels = intersect(MOD_ORD, unique(ev$mod)))
    }
    k <- length(cats)
    gy <- if (k) base - (k - 1) * step - step * 1.5 else base - step
    ylo <- gy - step * 1.1
    gts <- NULL
    if (!is.null(gene_pos)) {
      gl <- unique(dat$gene[dat$ADlocus == loc & !is.na(dat$gene) & dat$gene != ""])
      m  <- match(gl, gene_pos$gene)
      gts <- data.frame(gene = gl,
                        x = (as.numeric(gene_pos$start[m]) + as.numeric(gene_pos$end[m])) / 2,
                        chr = as.character(gene_pos$chr[m]), stringsAsFactors = FALSE)
      gts <- gts[!is.na(gts$x) & gts$chr == sub("^chr", "", w$chr) &
                 gts$x >= w$lo & gts$x <= w$hi, , drop = FALSE]
      if (!nrow(gts)) gts <- NULL else gts$y <- gy
    }
    dk2 <- FALSE
    ink2 <- if (dk2) "#e8eef2" else "#233947"
    br <- if (nrow(v)) { z <- pretty(c(0, H), 4); z[z >= 0 & z <= H * 1.02] } else numeric(0)
    p <- ggplot() +
      annotate("rect", xmin = w$lo, xmax = w$hi, ymin = ylo, ymax = -H * 0.05,
               fill = "#f6f8fa", colour = NA) +
      geom_hline(yintercept = -H * 0.05, colour = "#dde2e6", linewidth = 0.5) +
      annotate("text", x = w$lo, y = H * 1.05, hjust = 0, size = 3.4, fontface = "bold",
               colour = ink2,
               label = if (nrow(v)) "GWAS credible-set variants" else
                 "Credible-set variant positions (no per-variant p-value recorded)")
    if (nrow(v)) {
      p <- p +
        geom_hline(yintercept = GWS, colour = "#da532c", linetype = "22", linewidth = 0.5) +
        annotate("text", x = w$hi, y = GWS, hjust = 1, vjust = -0.6, size = 2.9,
                 colour = "#da532c", label = "genome-wide significance") +
        geom_point(data = v, aes(pos, log10pval), shape = 21, size = 2.3,
                   fill = "#a8bdd0", colour = "#ffffff", stroke = 0.4)
      ld <- v[which.max(v$log10pval), , drop = FALSE]
      p <- p + geom_point(data = ld, aes(pos, log10pval), shape = 21, size = 3.2,
                          fill = "#233947", colour = "#ffffff", stroke = 0.6) +
        geom_text(data = ld, aes(pos, log10pval,
                  label = ifelse(is.na(rsid) | rsid == "", "lead variant",
                                 paste0("lead variant · ", rsid))),
                  vjust = -1.1, size = 2.9, colour = ink2)
    } else if (length(vp)) {
      p <- p + geom_segment(data = data.frame(x = vp), aes(x = x, xend = x),
                            y = H * 0.18, yend = H * 0.62, colour = "#a8bdd0",
                            linewidth = 1.1)
    }
    if (k) {
      off <- ev[!ev$on, , drop = FALSE]; on <- ev[ev$on, , drop = FALSE]
      if (nrow(off))
        p <- p + geom_point(data = off, aes(gx, y, shape = mod, size = n),
                            fill = "#c9d0d5", colour = "#c9d0d5", alpha = 0.45, stroke = 0.4)
      if (nrow(on))
        p <- p + geom_point(data = on, aes(gx, y, fill = ctx, colour = ctx,
                                           shape = mod, size = n), stroke = 0.6) +
          geom_text(data = on, aes(gx, y, label = gene), vjust = -1.25, size = 2.8,
                    colour = ink2, check_overlap = TRUE)
    }
    if (!is.null(gts))
      p <- p +
        geom_hline(yintercept = gy + step * 0.75, colour = "#dde2e6", linewidth = 0.4) +
        geom_point(data = gts, aes(x, y), shape = 124, size = 3, colour = "#5c6b75") +
        geom_text(data = gts, aes(x, y, label = gene), vjust = 1.9, size = 2.7,
                  colour = "#5c6b75", check_overlap = TRUE)
    p +
      scale_x_continuous(limits = c(w$lo, w$hi), expand = c(0.012, 0),
                         labels = function(z) sprintf("%.2f", z / 1e6)) +
      scale_y_continuous(
        breaks = c(br, if (k) base - (seq_len(k) - 1) * step, if (!is.null(gts)) gy),
        labels = c(as.character(br), cats, if (!is.null(gts)) "genes in window"),
        limits = c(ylo, H * 1.1)) +
      scale_fill_manual(values = CTX_COL, name = "Cell type", drop = TRUE,
                        labels = function(z) unname(CTX_FULL[z])) +
      scale_colour_manual(values = CTX_COL, name = "Cell type", drop = TRUE,
                          labels = function(z) unname(CTX_FULL[z])) +
      scale_shape_manual(values = MOD_SHP, name = "QTL type", drop = TRUE) +
      scale_size_continuous(range = c(2.2, 6.2), name = "Supporting variants") +
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.4), order = 1),
             colour = "none",
             shape = guide_legend(override.aes = list(size = 3.4, fill = "#93a3ad",
                                                      colour = "#93a3ad"), order = 2),
             size = guide_legend(order = 3)) +
      labs(x = paste0("Position on ", w$chr, " (Mb)"), y = NULL) + g_theme +
      theme(legend.position = "bottom", legend.box = "horizontal",
            axis.text.y = element_text(size = 8.5, colour = ink2)) +
      dk_theme(FALSE)
  })

  output$p_region <- renderPlot({ region_fig() })

  output$dl_region_pdf <- downloadHandler(
    filename = function() paste0(selected_locus(), "_regional_evidence.pdf"),
    content  = function(f) ggplot2::ggsave(f, region_fig(), width = 11, height = 7.5,
                                           device = "pdf", useDingbats = FALSE))
  output$dl_region_png <- downloadHandler(
    filename = function() paste0(selected_locus(), "_regional_evidence.png"),
    content  = function(f) ggplot2::ggsave(f, region_fig(), width = 11, height = 7.5,
                                           dpi = 220, device = "png", bg = "white"))
  output$dl_region_csv <- downloadHandler(
    filename = function() paste0(selected_locus(), "_regional_evidence.csv"),
    content  = function(f) {
      e <- locus_ev()
      if (is.null(e) || !nrow(e)) { .prov_csv(data.frame(), f, "locus region evidence rows (empty)"); return(invisible()) }
      w <- locus_win()
      o <- data.frame(ADlocus = selected_locus(), window = paste0(w$chr, ":", round(w$lo), "-", round(w$hi)),
                      gene = e$gene, gene_chr = e$gchr, gene_position = round(e$gx),
                      context = e$ctx, xqtl_type = e$mod, tier = paste0("T", e$tier),
                      supporting_variants = e$n, in_window = e$inwin, stringsAsFactors = FALSE)
      .prov_csv(o[order(o$tier, o$gene), ], f, "locus region evidence rows")
    })

observeEvent(input$locus_prev, {
    i <- match(selected_locus(), loci)
    if (!is.na(i) && i > 1) selected_locus(loci[i - 1])
  })
  observeEvent(input$locus_next, {
    i <- match(selected_locus(), loci)
    if (!is.na(i) && i < length(loci)) selected_locus(loci[i + 1])
  })

  # ---- P0.2 interpretation panel and P0.3 evidence ledger ------------------
  .ev_state <- function(v) {
    v <- as.character(v)
    if (length(v) == 0L || is.na(v) || !nzchar(v)) return("not recorded in this release")
    if (identical(toupper(v), "TRUE")) "yes" else "tested, not significant"
  }
  .gene_tier <- function(d, g) {
    tv <- as.character(d$top_confidence[!is.na(d$gene) & d$gene == g])
    tv <- tv[!is.na(tv) & nzchar(tv)]
    if (!length(tv)) NA_character_ else sort(tv)[1]
  }
  output$locus_interp <- renderUI({
    loc <- selected_locus(); req(loc)
    d <- dat[dat$ADlocus == loc, , drop = FALSE]
    g <- sort(unique(d$gene[!is.na(d$gene) & nzchar(d$gene)]))
    if (!length(g)) return(div(class = "dc-note",
      "No gene is named at this locus in this release."))
    gt <- vapply(g, function(x) .gene_tier(d, x), character(1))
    ok <- !is.na(gt)
    best <- if (any(ok)) sort(gt[ok])[1] else NA_character_
    bg <- g[ok & gt == best]
    i <- match(loc, locus_registry$ADlocus)
    lead <- if (!is.na(i)) locus_registry$lead_rsid[i] else NA_character_
    subj <- if (!is.na(lead)) lead else loc
    s1 <- paste0(subj, " implicates ", length(g), if (length(g) == 1L) " gene at " else " genes at ", loc, ".")
    s2 <- if (is.na(best)) "No gene here carries an assigned tier." else
            paste0(paste(bg, collapse = " and "),
                   if (length(bg) == 1L) " has " else " have ",
                   if (identical(best, "T6")) "T6 support, which rests on TWAS or MR evidence without localized AD-xQTL support."
                   else paste0("localized ", best, " support."))
    t6 <- g[ok & gt == "T6"]
    t6 <- setdiff(t6, bg)
    s3 <- if (!length(t6)) "" else
            paste0(" ", paste(t6, collapse = ", "),
                   if (length(t6) == 1L) " is T6, indicating TWAS or MR evidence without localized AD-xQTL support."
                   else " are T6, indicating TWAS or MR evidence without localized AD-xQTL support.")
    s4 <- if (is.na(best)) "" else paste0(" The locus-level strongest tier is therefore ", best, ".")
    led <- lapply(g[ok], function(x) {
      r <- d[!is.na(d$gene) & d$gene == x, , drop = FALSE][1, , drop = FALSE]
      tt <- .gene_tier(d, x)
      loc_sup <- if (identical(tt, "T6")) "no" else "yes"
      tags$details(class = "led",
        tags$summary(paste0("Why is ", x, " ", tt, "?")),
        tags$ul(
          tags$li(paste0("Localized AD-xQTL support: ", loc_sup)),
          tags$li(paste0("TWAS: ", .ev_state(r$twas_sig))),
          tags$li(paste0("MR: ", .ev_state(r$mr_sig))),
          tags$li(paste0("cTWAS: ", .ev_state(r$ctwas_sig))),
          tags$li(paste0("Tier definition: ", if (tt %in% names(TIER_DEFS)) TIER_DEFS[[tt]] else "not defined"))))
    })
    div(class = "interp",
        div(class = "dc-kicker", "What this locus shows"),
        p(class = "interp-t", paste0(s1, " ", s2, s3, s4)),
        div(class = "interp-led", led),
        p(class = "interp-f",
          "Every statement above is taken from the rows in the table below. ",
          "Tiers order evidence strength; they are not statements of causality."))
  })
  output$locus_meta <- renderUI({
    loc <- selected_locus(); req(loc)
    i   <- match(loc, locus_registry$ADlocus)
    d   <- dat[dat$ADlocus == loc, ]
    reg  <- if (!is.na(i)) locus_registry$region[i]    else NA_character_
    lead <- if (!is.na(i)) locus_registry$lead_rsid[i] else NA_character_
    ng   <- length(unique(stats::na.omit(d$gene[nzchar(d$gene)])))
    nv   <- length(unique(stats::na.omit(d$variant_ID)))
    bt   <- as.character(sort(d$top_confidence)[1])
    bg   <- if (is.na(bt)) character(0) else sort(unique(stats::na.omit(
              d$gene[!is.na(d$top_confidence) & as.character(d$top_confidence) == bt &
                     !is.na(d$gene) & nzchar(d$gene)])))
    bgl  <- if (!length(bg)) NA_character_ else if (length(bg) <= 2L)
              paste(bg, collapse = ", ") else
              paste0(paste(bg[1:2], collapse = ", "), " and ", length(bg) - 2L, " more")
    sg   <- as.character(sort(d$significance)[1])
    kv <- function(k, v) tagList(span(class = "lk", k), span(v))
    kvt <- function(k, v, tip) tagList(span(class = "lk", info(k, tip)), span(v))
    div(class = "locus-record",
      div(class = "locus-title", if (is.na(reg)) loc else pretty_region(reg)),
      if (!is.na(i) && !is.na(locus_registry$start[i]) &&
          locus_registry$start[i] == locus_registry$end[i])
        div(class = "scope-note",
            "One variant position is recorded at this locus, so this is a single coordinate rather than a span."),
      div(class = "locus-sub",
        if (!is.na(lead)) kv("lead variant", tags$code(lead)),
        kvt("genes", ng, "Genes with evidence recorded at this locus. Some of them sit outside the region drawn above: a variant here can act on a gene far away, and that distal evidence is still counted for the locus."), kv("variants", nv),
        if (!is.na(bgl)) kvt("strongest gene here", tags$b(bgl), "The gene with the strongest evidence at this locus. Localized support means a variant here was fine-mapped to the gene. Gene-level support (TWAS, MR) links the gene to AD without pinning a variant. Distal support means the variant acts on a gene outside this locus."),
        if (!is.na(bt)) kvt("its tier", bt, "The best tier reached by any gene at this locus. A tier belongs to a gene, not to the locus, and a gene carries its best tier from anywhere in the release. So a locus can read T4 while a gene in it, such as APOE, reads T6, because that gene has TWAS or MR evidence only and no localized AD-xQTL support."),
        if (!is.na(sg)) kv("GWAS", sg)),
      div(class = "locus-id", "pipeline ID ", tags$code(loc)))
  })
  output$locus_tbl <- renderDT({
    d <- dat[dat$ADlocus == selected_locus(), ]
    datatable(tibble(Gene = d$gene, rsID = d$rsid, Sig = sig_badge(d$significance),
        Tier = conf_badge(d$top_confidence), `xQTL PIP` = round(d$xqtl_max_inclusion,3),
        `TWAS z` = round(d$max_twas_z,2), `Cell types` = d$ct_dots,
        `Ordered contexts` = ds_relabel(d$ordered_contexts)),
      escape = FALSE, rownames = FALSE, options = list(dom = "lrtip", pageLength = 15, scrollX = TRUE, headerCallback = hdr_js, columnDefs = link_defs(0, 1)),
      class = "display compact hover")
  }, server = TRUE)

  # ----- Trans tab -----
  trans_d <- reactive({
    d <- dat[na_fill(dat$has_trans, FALSE), ]
    if (length(input$tmod)) {
      colmap <- c("snRNA"="trans_snRNA_n_genes","pQTL"="trans_pQTL_n_genes",
                  "gpQTL"="trans_gpQTL_n_genes",
                  "Hotspot"="trans_hotspot_n_genes")
      keep <- rep(FALSE, nrow(d))
      for (m in input$tmod) keep <- keep | na_fill(d[[colmap[[m]]]] > 0, FALSE)
      d <- d[keep, ]
    }
    d
  })
  output$t_genes <- renderText(fmt(length(unique(na.omit(trans_d()$gene)))))
  output$t_loci  <- renderText(fmt(length(unique(trans_d()$ADlocus))))
  tn_tbl_d <- reactive({
    d <- trans_all
    if (is.null(d) || !nrow(d)) return(NULL)
    if (identical(input$tn_scope, "locus")) {
      loc <- input$trans_locus
      if (!is.null(loc) && nzchar(loc)) d <- d[d$locus == loc, , drop = FALSE]
    }
    mods <- input$tmod
    if (!is.null(mods) && length(mods)) d <- d[d$modality %in% mods, , drop = FALSE]
    qq <- trimws(if (is.null(input$tn_q)) "" else input$tn_q)
    if (nzchar(qq)) {
      pt  <- stringr::regex(qq, ignore_case = TRUE)
      hay <- paste(d$locus, d$source, d$target, d$rsid, d$context)
      hit <- stringr::str_detect(hay, pt)
      hit[is.na(hit)] <- FALSE
      d <- d[hit, , drop = FALSE]
    }
    d
  })

  output$trans_tbl <- renderDT({
    d <- tn_tbl_d()
    if (is.null(d)) return(NULL)
    cc <- ctx_cols(d$context)
    datatable(tibble(Locus = d$locus, `Source gene` = d$source, rsID = d$rsid,
                     Modality = d$modality, `Distal target` = d$target,
                     `Genes in program` = d$n_genes, Tier = d$tier,
                     `Cell type or region` = cc$ctx, Dataset = cc$dset,
                     `Assay context` = cc$mod),
      escape = FALSE, rownames = FALSE,
      options = list(dom = "lrtip", pageLength = 25, scrollX = TRUE),
      class = "display compact hover")
  }, server = TRUE)
  # ---- Genome view (IGV) ----
  observeEvent(list(selected_locus(), input$nav), {
    req(identical(input$nav, "Locus detail"))
    loc <- selected_locus(); req(loc)
    i <- match(loc, locus_registry$ADlocus); req(!is.na(i))
    s <- max(1, locus_registry$start[i] - 50000)
    e <- locus_registry$end[i] + 50000
    session$sendCustomMessage("igv_locus",
      list(locus = sprintf("%s:%.0f-%.0f", locus_registry$chr[i], s, e)))
  }, ignoreInit = FALSE)

  # ---- Trans network for the selected locus ----
  trans_edges <- reactive({
    loc <- selected_locus(); req(loc)
    d <- dat[dat$ADlocus == loc, ]
    d <- d[!is.na(d$trans_genes) & nzchar(trimws(as.character(d$trans_genes))), ]
    if (!nrow(d)) return(NULL)
    ed <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      tg <- unique(trimws(unlist(strsplit(as.character(d$trans_genes[i]), "[;,|]"))))
      tg <- tg[nzchar(tg)]
      sg <- as.character(d$gene[i])
      if (!length(tg) || is.na(sg) || !nzchar(sg)) return(NULL)
      data.frame(src = sg, tgt = tg, stringsAsFactors = FALSE)
    }))
    if (is.null(ed) || !nrow(ed)) return(NULL)
    unique(ed)
  })

  observeEvent(input$trans_locus, {
    if (!identical(input$trans_locus, selected_locus())) selected_locus(input$trans_locus)
  }, ignoreInit = TRUE)
  observeEvent(selected_locus(), {
    if (!identical(input$trans_locus, selected_locus()))
      updateSelectInput(session, "trans_locus", selected = selected_locus())
  }, ignoreNULL = TRUE)


  output$trans_circos <- renderPlot({
    req(input$tn_loc)
    trans_chord(input$tn_loc)
  }, res = 96)

  .transnet_ui <- function(loc) {
    if (is.null(trans_pairs))
      return(div(class = "note", "Target gene coordinates are not available in this build."))
    d <- trans_pairs[trans_pairs$ADlocus == loc, , drop = FALSE]
    if (!nrow(d))
      return(div(class = "note",
        "No trans evidence with mapped coordinates is recorded for this locus. That is ",
        "not evidence against distal regulation; none was reported here."))
    i <- match(loc, locus_registry$ADlocus)
    srcreg <- if (!is.na(i)) pretty_region(locus_registry$region[i]) else loc
    # one line per distal gene, with every assay that supports it, so the count
    # here matches the gene count used everywhere else
    d$modality <- as.character(d$modality)
    d <- d %>%
      dplyr::group_by(tgt, tgt_chr, tgt_pos) %>%
      dplyr::summarise(
        modality = paste(sort(unique(as.character(modality))), collapse = ", "),
        src      = paste(sort(unique(as.character(src))), collapse = ", "),
        src_chr  = dplyr::first(src_chr),
        src_pos  = dplyr::first(src_pos),
        .groups  = "drop") %>%
      as.data.frame()
    d <- d[order(match(d$tgt_chr, CHR_ORD), d$tgt_pos), , drop = FALSE]
    nc <- length(unique(d$tgt_chr))
    if (nrow(d) <= 6) {
      illus <- lapply(seq_len(nrow(d)), function(k)
        HTML(trans_illu_one(d$src_chr[k], d$src_pos[k], srcreg,
                            d$tgt_chr[k], d$tgt_pos[k], d$tgt[k],
                            d$modality[k], d$src[k])))
      return(tagList(
        div(class = "tl-head",
            div(class = "tl-h-src", srcreg),
            div(class = "tl-h-n",
                sprintf("%d distal gene%s on %d chromosome%s", nrow(d),
                        if (nrow(d) == 1) "" else "s", nc,
                        if (nc == 1) "" else "s"))),
        div(class = "ti-wrap", illus)))
    }
    rows <- lapply(seq_len(nrow(d)), function(k)
      div(class = "tl-row",
        div(class = "tl-coord", sprintf("chr%s:%s", d$tgt_chr[k], mbp(d$tgt_pos[k]))),
        div(class = "tl-gene", d$tgt[k]),
        div(class = "tl-mod", as.character(d$modality[k])),
        div(class = "tl-src", "from ", d$src[k])))
    tagList(
      div(class = "tl-head",
          div(class = "tl-h-src", srcreg),
          div(class = "tl-h-n", sprintf("%d distal target%s on %d chromosome%s",
              nrow(d), if (nrow(d) == 1) "" else "s",
              nc, if (nc == 1) "" else "s"))),
      div(class = "tl-list", rows))
  }

  output$p_transnet_ui <- renderUI({
    loc <- selected_locus(); req(loc)
    .transnet_ui(loc)
  })

  # a circle needs arcs to be worth drawing; below three distal genes the
  # chromosome strands say it more plainly
  output$tn_figure <- renderUI({
    loc <- input$tn_loc
    if (is.null(loc) || !nzchar(loc)) return(div(class = "dc-note", "Pick a locus to draw its distal targets."))
    n <- if (is.null(trans_pairs)) 0L else
      length(unique(trans_pairs$tgt[trans_pairs$ADlocus == loc]))
    if (n >= 3)
      return(tagList(
        div(class = "dc-exp", span(class = "dc-kicker", "Export"),
            downloadLink("dl_trans_png", "PNG", class = "dc-btn"),
            downloadLink("dl_trans_csv", "CSV", class = "dc-btn")),
        if (n > 40) div(class = "dc-note", style = "margin-bottom:10px", paste0("Showing the 40 distal targets with the most assay support, of ", format(n, big.mark = ","), ". The full list is in the table below.")),
        if (identical(input$tn_view, "Bars"))
          withSpinner(plotOutput("trans_bars", height = "680px"),
                      type = 8, color = "#2a78d6", size = 0.6)
        else
          withSpinner(plotOutput("trans_circos", height = "680px"),
                      type = 8, color = "#2a78d6", size = 0.6)))
    tagList(
      div(class = "dc-note", style = "margin-bottom:14px",
          sprintf("%s distal gene%s at this locus, too few for a circle. Shown on the chromosome strands instead.",
                  n, ifelse(n == 1, "", "s"))),
      div(class = "dc-exp", span(class = "dc-kicker", "Export"),
          downloadLink("dl_trans_svg", "SVG", class = "dc-btn"),
          downloadLink("dl_trans_csv", "CSV", class = "dc-btn")),
      .transnet_ui(loc))
  })

  .go_example <- function(term) {
    if (grepl("^ENSG", term) && term %in% gid_map)
      term <- names(gid_map)[match(term, gid_map)][1]
    i <- match(term, search_index$term)
    if (is.na(i)) {
      showNotification(sprintf("%s is not in this release.", term), type = "message", duration = 4)
      return(invisible())
    }
    selected_locus(search_index$ADlocus[i])
    updateTextInput(session, "q", value = term)
    nav_select("nav", "Locus detail", session = session)
  }
  observeEvent(input$ex_bin1,  .go_example("BIN1"))
  observeEvent(input$ex_apoe,  .go_example("APOE"))
  observeEvent(input$ex_trem2, .go_example("TREM2"))
  observeEvent(input$ex_rs,    .go_example("rs429358"))
  observeEvent(input$ex_ensg,   .go_example("ENSG00000130203"))


  output$p_transnet <- renderPlot({
    ed <- trans_edges()
    validate(need(!is.null(ed) && nrow(ed) > 0,
                  "No trans target genes are recorded for this locus in this build."))
    if (nrow(ed) > 250) ed <- ed[seq_len(250), ]
    sN <- sort(unique(ed$src)); tN <- sort(unique(ed$tgt))
    sy <- stats::setNames(seq_along(sN) / (length(sN) + 1), sN)
    ty <- stats::setNames(seq_along(tN) / (length(tN) + 1), tN)
    ed$y1 <- sy[ed$src]; ed$y2 <- ty[ed$tgt]
    ggplot(ed) +
      geom_segment(aes(x = 0, xend = 1, y = y1, yend = y2), colour = LINE, linewidth = 0.4) +
      geom_point(aes(x = 0, y = y1), colour = TEAL, size = 2.6) +
      geom_point(aes(x = 1, y = y2), colour = "#eb6834", size = 2.2) +
      geom_text(data = data.frame(y = as.numeric(sy), lab = names(sy)),
                aes(x = 0, y = y, label = lab), hjust = 1.12, size = 3.4, colour = INK) +
      geom_text(data = data.frame(y = as.numeric(ty), lab = names(ty)),
                aes(x = 1, y = y, label = lab), hjust = -0.12, size = 3.1, colour = MUTE) +
      scale_x_continuous(limits = c(-0.5, 1.5)) +
      labs(x = NULL, y = NULL) + g_theme +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
  })

  observeEvent(input$ex_region, {
    .r <- resolve_locus(list(region = "chr2:127088929-127118906"))
    if (!is.null(.r)) {
      selected_locus(.r$locus)
      nav_select("nav", "Locus detail", session = session)
    }
  })

  output$dl_all <- downloadHandler(
    filename = function() paste0("AD_xQTL_loci_", Sys.Date(), ".csv"),
    content  = function(f) .prov_csv(dat, f, "complete evidence table"))

  output$dl_trans <- downloadHandler(
    filename = function() paste0("AD_xQTL_trans_", Sys.Date(), ".csv"),
    content  = function(f) .prov_csv(tn_tbl_d(), f, "trans pairs matching the table filters"))

  .fig_dl <- function(fn, fig, w = 11, h = 6) list(
    pdf = downloadHandler(
      filename = function() paste0(fn, ".pdf"),
      content  = function(f) ggplot2::ggsave(f, fig(), width = w, height = h,
                                             device = "pdf", useDingbats = FALSE)),
    png = downloadHandler(
      filename = function() paste0(fn, ".png"),
      content  = function(f) ggplot2::ggsave(f, fig(), width = w, height = h,
                                             dpi = 220, device = "png", bg = "white")))
  .dlA <- .fig_dl("AD_locus_landscape", fig_p_landscape, 12, 5)
  output$dl_figA_pdf <- .dlA$pdf
  output$dl_figA_png <- .dlA$png
  .dlB <- .fig_dl("cell_type_by_tier", fig_p_ctxtier, 11, 6)
  output$dl_figB_pdf <- .dlB$pdf
  output$dl_figB_png <- .dlB$png
  output$dl_figA_csv <- downloadHandler(
    filename = function() "AD_locus_landscape.csv",
    content  = function(f) .prov_csv(locus_registry, f, "locus registry"))

  output$rv_win <- renderText({
    w <- locus_win()
    paste0(w$chr, ":", format(round(w$lo), big.mark = ","), "–",
           format(round(w$hi), big.mark = ","),
           " · same window as the browser above")
  })

  output$doc_assays <- renderUI({
    e <- parse_ctx_tokens(dat)
    if (is.null(e) || !nrow(e)) return(div(class = "dc-note", "No assay records."))
    cx <- intersect(CTX_ORD, unique(e$ctx))
    md <- intersect(MOD_ORD, unique(e$mod))
    hdr <- tags$tr(tags$th("Cell type"), lapply(md, tags$th), tags$th("Gene records"))
    rows <- lapply(cx, function(k) {
      sub <- e[e$ctx == k, , drop = FALSE]
      tags$tr(
        tags$td(span(style = sprintf(
          paste0("display:inline-flex;align-items:center;gap:8px")),
          span(style = sprintf(
            "width:9px;height:9px;border-radius:50%%;background:%s;display:inline-block",
            CTX_COL[[k]])),
          if (k %in% names(CTX_FULL)) CTX_FULL[[k]] else k)),
        lapply(md, function(m) {
          n <- sum(sub$mod == m)
          if (n == 0) tags$td(class = "z", HTML("&middot;")) else tags$td(format(n, big.mark = ","))
        }),
        tags$td(tags$b(format(nrow(sub), big.mark = ","))))
    })
    tagList(
      tags$table(class = "doc-tbl", tags$thead(hdr), tags$tbody(rows)),
      div(class = "dc-note", style = "margin-top:10px",
        "Counts are gene records in this release carrying that cell type and assay; a ",
        HTML("&middot;"), " means the combination does not appear."))
  })

  .dl_file <- function(path, name) downloadHandler(
    filename = function() name,
    content  = function(f) file.copy(path, f, overwrite = TRUE))
  output$dlf_browser <- .dl_file("data.csv", "AD_locus_evidence_2026-09.csv")
  output$dlf_genepos <- .dl_file("gene_positions.csv", "gene_coordinates_GRCh38_2026-09.csv")
  output$dlf_locsum  <- .dl_file(
    "downloads/AD_locus_summary_release.csv",
    "AD_locus_summary_2026-09.csv")
  output$dlf_varlvl  <- .dl_file(
    "downloads/AD_locus_variants_release.csv.gz",
    "AD_locus_variants_2026-09.csv.gz")
  output$dlf_tier   <- .dl_file("downloads/gene_tier_assignment_release.csv",
    "AD_gene_tier_assignment_2026-09.csv")
  output$dlf_xlsx    <- .dl_file(
    "downloads/AD_loci_xQTL_summary_release.xlsx",
    "AD_locus_xQTL_summary_2026-09.xlsx")
  output$dl_zip <- downloadHandler(
    filename = function() "AD_loci_explorer_2026-09.zip",
    content  = function(f) {
      fs <- c("data.csv", "gene_positions.csv",
              list.files("downloads", full.names = TRUE))
      fs <- fs[file.exists(fs)]
      if (requireNamespace("zip", quietly = TRUE)) zip::zipr(f, fs)
      else utils::zip(f, fs, flags = "-j9X")
    })

  lt_term <- reactiveVal("")
  observeEvent(input$lt_qgo,   lt_term(trimws(input$lt_q %||% "")))
  observeEvent(input$lt_qclear, { lt_term(""); updateTextInput(session, "lt_q", value = "") })

  lt_hits <- reactive({
    q <- trimws(lt_term())
    if (!nzchar(q)) return(NULL)
    if (grepl("^ENSG", q, ignore.case = TRUE)) {
      g <- names(gid_map)[match(toupper(q), toupper(unname(gid_map)))]
      if (!is.na(g)) q <- g
    }
    i <- match(toupper(q), toupper(search_index$term))
    if (!is.na(i)) return(unique(search_index$ADlocus[i]))
    k <- match(toupper(q), toupper(locus_registry$lead_rsid))
    if (!is.na(k)) return(unique(locus_registry$ADlocus[k]))
    r <- try(resolve_locus(list(region = q)), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r) && !is.null(r$locus)) return(unique(r$locus))
    j <- grep(q, search_index$term, ignore.case = TRUE, fixed = TRUE)
    if (length(j)) return(unique(stats::na.omit(search_index$ADlocus[j])))
    character(0)
  })

  output$lt_qnote <- renderUI({
    h <- lt_hits()
    if (is.null(h)) return(NULL)
    if (!length(h)) return(span(style = "color:var(--cta)",
      sprintf("Nothing in this release matches %s.", lt_term())))
    span(sprintf("%s matches %d locus%s.", lt_term(), length(h),
                 ifelse(length(h) == 1, "", "es")))
  })

  loci_view <- reactive({
    t <- locus_table
    if (!is.null(input$lt_chr) && nzchar(input$lt_chr))
      t <- t[t$chr == input$lt_chr, , drop = FALSE]
    if (!is.null(input$lt_tier) && length(input$lt_tier))
      t <- t[!is.na(t$best) & t$best %in% input$lt_tier, , drop = FALSE]
    if (!is.null(input$lt_ctx) && length(input$lt_ctx))
      t <- t[vapply(strsplit(t$cells, ","), function(v)
               any(input$lt_ctx %in% v), logical(1)), , drop = FALSE]
    if (!is.null(input$lt_mod) && length(input$lt_mod))
      t <- t[vapply(strsplit(t$mods, ","), function(v)
               any(input$lt_mod %in% trimws(v)), logical(1)), , drop = FALSE]
    h <- lt_hits()
    if (!is.null(h)) t <- t[t$ADlocus %in% h, , drop = FALSE]
    t
  })

  observe({
    updateSelectizeInput(session, "lt_ctx",
      choices = stats::setNames(CTX_ORD,
        ifelse(CTX_ORD %in% names(CTX_FULL), CTX_FULL[CTX_ORD], CTX_ORD)),
      selected = isolate(input$lt_ctx), server = FALSE)
  })

  output$ct_key <- renderUI({
    div(
      div(class = "ctk-h", "Cell type colours"),
      div(class = "ctkey",
        lapply(CTX_ORD, function(k)
          span(class = "ctk-i",
            span(class = "ctk-d", style = sprintf("background:%s", unname(CTX_COL[[k]]))),
            unname(CTX_FULL[[k]])))))
  })

  output$lt_stat <- renderUI(HTML(paste0(
    "<b>", format(nrow(loci_view()), big.mark = ","), "</b> of ",
    format(nrow(locus_table), big.mark = ","), " loci shown")))

  output$loci_tbl <- renderDT({
    t <- loci_view()
    if (!nrow(t)) return(datatable(data.frame(), rownames = FALSE))
    chip <- function(k) if (is.na(k)) "<span class=\"lt-na\">&mdash;</span>" else
      sprintf('<span class="tchip" style="background:%s;color:%s">%s</span>',
              unname(conf_pal[[k]]), unname(conf_txt[[k]]), k)
    d <- data.frame(
      Locus = paste0('<div class="lt-k">', t$ADlocus, '</div><div class="lt-co">',
                     vapply(t$region, pretty_region, character(1)), '</div>'),
      "Lead variant" = ifelse(is.na(t$lead) | t$lead == "", "<span class=\"lt-na\">&mdash;</span>",
                              paste0('<span class="dc-mono">', t$lead, '</span>')),
      "Best tier" = vapply(t$best, chip, character(1)),
      "Top gene" = ifelse(is.na(t$top_gene), "<span class=\"lt-na\">&mdash;</span>",
                          paste0('<span class="dc-mono" style="font-weight:600">',
                                 t$top_gene, '</span>')),
      "Cell types" = vapply(t$cells, lt_cells, character(1)),
      Modalities = ifelse(t$mods == "", "<span class=\"lt-na\">&mdash;</span>",
                          paste0('<span class="dc-mono" style="font-size:11.5px">',
                                 t$mods, '</span>')),
      Genes = t$n_genes,
      Records = t$n_rows,
      "Evidence by tier" = vapply(seq_len(nrow(t)), function(i) lt_bar(t[i, ]), character(1)),
      check.names = FALSE, stringsAsFactors = FALSE)
    datatable(d, escape = FALSE, rownames = FALSE, selection = "single",
      options = list(pageLength = 15, dom = "tip", autoWidth = FALSE,
        columnDefs = list(list(className = "dt-right", targets = c(6, 7)),
                          list(orderable = FALSE, targets = c(4, 8)))))
  })

  observeEvent(input$loci_tbl_rows_selected, {
    i <- input$loci_tbl_rows_selected
    req(length(i) == 1)
    selected_locus(loci_view()$ADlocus[i])
    nav_select("nav", "Locus detail", session = session)
  })

  output$dl_loci_csv <- downloadHandler(
    filename = function() "AD_loci_overview_2026-09.csv",
    content  = function(f) {
      t <- loci_view()
      o <- data.frame(ADlocus = t$ADlocus, region = t$region, lead_variant = t$lead,
                      best_tier = t$best, top_gene = t$top_gene, cell_types = t$cells,
                      modalities = t$mods, n_genes = t$n_genes, n_records = t$n_rows,
                      T1 = t$t1, T2 = t$t2, T3 = t$t3, T4 = t$t4, T5 = t$t5, T6 = t$t6,
                      stringsAsFactors = FALSE)
      .prov_csv(o, f, "loci table, current filters")
    })

  output$home_river <- plotly::renderPlotly({
    p <- home_river(FALSE)
    if (is.null(p)) return(plotly::plot_ly(type = "sankey"))
    p
  })

  output$ctx_tier_grid <- renderUI({
    m <- ctx_tier_counts
    if (is.null(m)) return(div(class = "dc-note", "No cell-type records."))
    tot <- rowSums(m); mx <- max(tot)
    hdr <- tagList(div(class = "cb-h"), div(class = "cb-h"),
      NULL)
    rows <- lapply(rownames(m), function(k) {
      v <- m[k, ]; s <- sum(v)
      rtip <- paste(sprintf("%s: %d", TIER_SEQ, as.integer(v)), collapse = "   ")
      segs <- lapply(seq_len(6), function(i) {
        if (v[i] <= 0) return(NULL)
        span(style = sprintf("width:%.3f%%;background:%s;color:%s", 100 * v[i] / mx,
                             unname(conf_pal[[TIER_SEQ[i]]]), unname(conf_txt[[TIER_SEQ[i]]])),
             title = sprintf("%s %s: %d genes", k, TIER_SEQ[i], v[i]),
             if (100 * v[i] / mx > 3.2) as.character(v[i]) else "")
      })
      tagList(
        div(class = "cb-lab",
          span(class = "cb-dot", style = sprintf("background:%s", unname(CTX_COL[[k]]))),
          if (k %in% names(CTX_FULL)) CTX_FULL[[k]] else k),
        div(class = "cb-bar", title = paste0(k, "  ", rtip, "   Total: ", format(s, big.mark = ",")), segs, span(class = "cb-tot", format(s, big.mark = ","))))
    })
    div(class = "cb-wrap", div(class = "cb-grid", hdr, rows))
  })

  output$dl_figB_csv <- downloadHandler(
    filename = function() "cell_type_by_tier_2026-09.csv",
    content  = function(f) {
      m <- ctx_tier_counts
      o <- data.frame(cell_type = rownames(m), m, total = rowSums(m),
                      check.names = FALSE, stringsAsFactors = FALSE)
      .prov_csv(o, f, "cell type by tier counts")
    })

  observe({
    updateSelectizeInput(session, "gp_gene", choices = genes,
                         selected = isolate(input$gp_gene), server = TRUE)
  })

  gp_rows <- reactive({
    g <- input$gp_gene
    if (is.null(g) || !nzchar(g)) return(NULL)
    dat[!is.na(dat$gene) & dat$gene == g, , drop = FALSE]
  })
  gp_ev <- reactive({ d <- gp_rows(); if (is.null(d) || !nrow(d)) NULL else parse_ctx_tokens(d) })

  output$gp_head <- renderUI({
    g <- input$gp_gene
    if (is.null(g) || !nzchar(g))
      return(div(class = "gp-h", div(class = "dc-note",
        "Start with a gene. Type a symbol such as BIN1 or TREM2, or an Ensembl ID, ",
        "in the box above. You will get its tier, the cell types and assays behind ",
        "it, and every AD locus where it is implicated.")))
    d <- gp_rows(); e <- gp_ev()
    k <- match(as.character(d$top_confidence), TIER_SEQ)
    best <- if (any(!is.na(k))) TIER_SEQ[min(k, na.rm = TRUE)] else NA_character_
    co <- ""
    if (!is.null(gene_pos)) {
      i <- match(g, gene_pos$gene)
      if (!is.na(i)) co <- sprintf("chr%s:%s–%s · GRCh38",
        gene_pos$chr[i], format(as.numeric(gene_pos$start[i]), big.mark = ","),
        format(as.numeric(gene_pos$end[i]), big.mark = ","))
    }
    if (g %in% names(gname_map))
      co <- paste0(co, if (nzchar(co)) "  \u00b7  " else "", unname(gname_map[[g]]))
    if (g %in% names(gid_map))
      co <- paste0(co, if (nzchar(co)) "  \u00b7  " else "", unname(gid_map[[g]]))
    nl <- length(unique(d$ADlocus))
    nc <- if (is.null(e)) 0L else length(unique(e$ctx))
    nm <- if (is.null(e)) 0L else length(unique(e$mod))
    tagList(
      div(class = "gp-h",
        div(class = "gp-h1", h1(class = "gp-sym", g), span(class = "gp-co", co)),
        div(class = "gp-chips",
          if (!is.na(best)) span(class = "gp-chip t", best),
          span(class = "gp-chip", paste(nl, if (nl == 1) "locus" else "loci")),
          span(class = "gp-chip", paste(nc, if (nc == 1) "cell type" else "cell types")),
          span(class = "gp-chip", paste(nm, if (nm == 1) "modality" else "modalities")))),
      div(class = "gp-act",
        downloadLink("dl_gp_all", "CSV", class = "dc-btn"),
        actionLink("gp_open", "Open in locus view →", class = "dc-btn")))
  })

output$gp_matrix <- renderUI({
    g <- input$gp_gene
    if (is.null(g) || !nzchar(g)) return(NULL)
    e <- gp_ev()
    if (is.null(e) || !nrow(e))
      return(div(class = "dc-note", style = "margin-top:10px",
                 "No cell-type or assay records for this gene."))
    md <- intersect(MOD_ORD, unique(e$mod))
    cx <- intersect(CTX_ORD, unique(e$ctx))
    tags$table(class = "gp-tbl",
      tags$thead(tags$tr(tags$th("Cell type"), lapply(md, tags$th))),
      tags$tbody(lapply(cx, function(k) {
        s <- e[e$ctx == k, , drop = FALSE]
        tags$tr(
          tags$td(span(style = "display:inline-flex;align-items:center;gap:8px",
            span(style = sprintf(
              "width:9px;height:9px;border-radius:50%%;background:%s;display:inline-block",
              unname(CTX_COL[[k]]))),
            if (k %in% names(CTX_FULL)) CTX_FULL[[k]] else k)),
          lapply(md, function(m) {
            w <- s[s$mod == m, , drop = FALSE]
            if (!nrow(w)) return(tags$td(class = "gp-no", HTML("&mdash;")))
            tn <- paste0("T", min(w$tier, na.rm = TRUE))
            if (!tn %in% names(conf_pal))
              return(tags$td(class = "gp-yes", HTML("&#10003;")))
            tags$td(class = "gp-yes",
              span(class = "tchip",
                   style = sprintf("background:%s;color:%s",
                                   unname(conf_pal[[tn]]), unname(conf_txt[[tn]])),
                   title = sprintf("Best tier reached in %s here", m), tn))
          }))
      })))
  })

  gp_loci_tbl <- reactive({
    g <- input$gp_gene
    if (is.null(g) || !nzchar(g)) return(NULL)
    d <- gp_rows()
    o <- NULL
    if (!is.null(d) && nrow(d)) {
      u <- d[!duplicated(d$ADlocus), , drop = FALSE]
      o <- data.frame(locus = as.character(u$ADlocus),
        region = vapply(as.character(u$ADlocus), function(l)
          pretty_region(locus_registry$region[match(l, locus_registry$ADlocus)]), character(1)),
        tier = as.character(u$top_confidence), assignment = "cis",
        stringsAsFactors = FALSE)
    }
    if (!is.null(trans_pairs)) {
      t <- trans_pairs[trans_pairs$tgt == g, , drop = FALSE]
      if (nrow(t)) {
        t <- t[!duplicated(t$ADlocus), , drop = FALSE]
        o <- rbind(o, data.frame(locus = as.character(t$ADlocus),
          region = vapply(as.character(t$ADlocus), function(l)
            pretty_region(locus_registry$region[match(l, locus_registry$ADlocus)]), character(1)),
          tier = NA_character_, assignment = "trans", stringsAsFactors = FALSE))
      }
    }
    o
  })

  output$gp_loci <- renderUI({
    g <- input$gp_gene
    if (is.null(g) || !nzchar(g))
        return(NULL)
    o <- gp_loci_tbl()
    if (is.null(o) || !nrow(o))
      return(div(class = "dc-note", style = "margin-top:10px", "No loci recorded."))
    tags$table(class = "gp-tbl",
      tags$thead(tags$tr(tags$th("Locus"), tags$th("Region"), tags$th("Tier"),
                         tags$th("Assignment"))),
      tags$tbody(lapply(seq_len(nrow(o)), function(i) tags$tr(
        class = "gp-rowlink",
        onclick = sprintf("Shiny.setInputValue(%s, %s, {priority: %s})",
                          shQuote("gp_row_locus"), shQuote(o$locus[i]), shQuote("event")),
        tags$td(tags$b(o$locus[i])),
        tags$td(class = "dc-mono", style = "font-size:11.5px", o$region[i]),
        tags$td(if (is.na(o$tier[i])) HTML("&mdash;") else
          span(class = "tchip", style = sprintf("background:%s;color:%s",
            unname(conf_pal[[o$tier[i]]]), unname(conf_txt[[o$tier[i]]])), o$tier[i])),
        tags$td(o$assignment[i])))))
  })

  observeEvent(input$gp_open, {
    o <- gp_loci_tbl(); req(!is.null(o), nrow(o) > 0)
    selected_locus(o$locus[1])
    nav_select("nav", "Locus detail", session = session)
  })

  output$dl_gp_matrix <- downloadHandler(
    filename = function() paste0(input$gp_gene, "_evidence_matrix.csv"),
    content  = function(f) .prov_csv(gp_ev(), f, "gene page evidence rows"))
  output$dl_gp_loci <- downloadHandler(
    filename = function() paste0(input$gp_gene, "_loci.csv"),
    content  = function(f) .prov_csv(gp_loci_tbl(), f, "gene page loci"))
  output$dl_gp_all <- downloadHandler(
    filename = function() paste0(input$gp_gene, "_records.csv"),
    content  = function(f) .prov_csv(gp_rows(), f, "gene page variant rows"))

  observeEvent(input$brand_home, nav_select("nav", "Home", session = session))
  observeEvent(list(input$go_search, input$go_search2), {
    nav_select("nav", "Home", session = session)
    session$sendCustomMessage("focus_search", list())
  }, ignoreInit = TRUE)

  output$land_tip <- renderUI({
    h <- input$land_hover
    if (is.null(h)) return(NULL)
    d <- land_df(); req(nrow(d) > 0)
    n <- nearPoints(d, h, xvar = "x", yvar = "y", threshold = 14, maxpoints = 1)
    if (!nrow(n)) return(NULL)
    lft <- h$coords_css$x + 14
    top <- h$coords_css$y - 10
    tier <- as.character(n$tier[1])
    col  <- if (tier %in% names(conf_pal)) unname(conf_pal[[tier]]) else "#cfd6db"
    div(class = "land-tip",
        style = sprintf("left:%.0fpx;top:%.0fpx", lft, top),
        div(class = "g", n$ADlocus[1]),
        div(class = "m", sprintf("chr%s · %s", n$chr[1], mbp(n$pos[1]))),
        div(span(style = sprintf(
              "display:inline-block;width:9px;height:9px;background:%s;margin-right:6px", col)),
            tier, " · ", n$ng[1], if (n$ng[1] == 1) " gene" else " genes"),
        if (!is.na(n$gene[1])) div(class = "m", "top gene ", tags$b(n$gene[1])),
        div(class = "m", sprintf("-log10 p %.1f", n$y[1])),
        div(class = "m", style = "margin-top:5px", "click to open this locus"))
  })

  observeEvent(input$land_click, {
    d <- land_df(); req(nrow(d) > 0)
    n <- nearPoints(d, input$land_click, xvar = "x", yvar = "y",
                    threshold = 16, maxpoints = 1)
    req(nrow(n) == 1)
    selected_locus(as.character(n$ADlocus[1]))
    nav_select("nav", "Locus detail", session = session)
  })

  .kbtn <- function(id, val, label, text) {
    sq <- intToUtf8(39)
    tags$button(type = "button", class = "kbd-pt", `aria-label` = label,
      onclick = paste0("Shiny.setInputValue(", sq, id, sq, ", ", sq, val, sq,
                       ", {priority:", sq, "event", sq, "})"), text)
  }

  output$land_kbd <- renderUI({
    d <- land_df(); req(nrow(d) > 0)
    d <- d[order(-d$y), , drop = FALSE]
    rows <- lapply(seq_len(nrow(d)), function(i) {
      lab <- paste0(d$ADlocus[i], ". Chromosome ", d$chr[i], " at ", mbp(d$pos[i]),
                    ". Tier ", as.character(d$tier[i]), ". ", d$ng[i],
                    if (d$ng[i] == 1L) " gene" else " genes",
                    ". Minus log 10 p ", sprintf("%.1f", d$y[i]),
                    if (!is.na(d$gene[i])) paste0(". Top gene ", d$gene[i]) else "",
                    ". Press Enter to open this locus.")
      tags$tr(
        tags$td(.kbtn("land_pick", d$ADlocus[i], lab, d$ADlocus[i])),
        tags$td(paste0("chr", d$chr[i], " ", mbp(d$pos[i]))),
        tags$td(sprintf("%.1f", d$y[i])),
        tags$td(d$ng[i]),
        tags$td(as.character(d$tier[i])),
        tags$td(if (is.na(d$gene[i])) intToUtf8(8212) else d$gene[i]))
    })
    tags$details(class = "figdata", id = "landscape-data",
      tags$summary(paste0("Read this figure as a table (", nrow(d),
                          " loci, keyboard accessible)")),
      div(class = "figdata-scroll", role = "region", tabindex = "0",
          `aria-label` = "Loci drawn in the landscape figure",
        tags$table(class = "figtbl",
          tags$caption(paste0(nrow(d),
            " loci drawn, strongest lead variant first. Every row opens the same locus page the point opens.")),
          tags$thead(tags$tr(
            tags$th(scope = "col", "Locus"), tags$th(scope = "col", "Position"),
            tags$th(scope = "col", "Minus log10 p"), tags$th(scope = "col", "Genes"),
            tags$th(scope = "col", "Tier"), tags$th(scope = "col", "Top gene"))),
          tags$tbody(rows))))
  })

  observeEvent(input$land_pick, {
    lp <- input$land_pick; req(!is.null(lp), nzchar(lp))
    selected_locus(as.character(lp))
    nav_select("nav", "Locus detail", session = session)
  })

  .river_rows <- local({
    e <- .sk_data
    if (is.null(e) || !nrow(e)) NULL else
      rbind(.sk_links(e, "ctxf", "modg"), .sk_links(e, "modg", "tierl"))
  })

  output$river_kbd <- renderUI({
    r <- .river_rows
    if (is.null(r) || !nrow(r)) return(NULL)
    r <- r[order(-r$Freq), , drop = FALSE]
    rows <- lapply(seq_len(nrow(r)), function(i) {
      key <- paste(r$s[i], r$t[i], r$b[i], sep = "|")
      lab <- paste0(r$s[i], " to ", r$t[i], ", tier ", r$b[i], ". ",
                    format(r$Freq[i], big.mark = ","), " records, ",
                    format(r$ng[i], big.mark = ","), " genes, ",
                    format(r$nl[i], big.mark = ","),
                    " loci. Press Enter to list the supporting genes and loci.")
      tags$tr(
        tags$td(r$s[i]), tags$td(r$t[i]), tags$td(r$b[i]),
        tags$td(format(r$Freq[i], big.mark = ",")),
        tags$td(format(r$ng[i], big.mark = ",")),
        tags$td(format(r$nl[i], big.mark = ",")),
        tags$td(.kbtn("river_pick", key, lab, "View supporting genes and loci")))
    })
    tagList(
      tags$details(class = "figdata", id = "river-data",
        tags$summary(paste0("Read this figure as a table (", nrow(r),
                            " flows, keyboard accessible)")),
        div(class = "figdata-scroll", role = "region", tabindex = "0",
            `aria-label` = "Flows drawn in the cell type to assay to tier figure",
          tags$table(class = "figtbl",
            tags$caption(paste0(nrow(r),
              " flows, widest first. Each row is one ribbon in the figure.")),
            tags$thead(tags$tr(
              tags$th(scope = "col", "From"), tags$th(scope = "col", "To"),
              tags$th(scope = "col", "Tier"), tags$th(scope = "col", "Records"),
              tags$th(scope = "col", "Genes"), tags$th(scope = "col", "Loci"),
              tags$th(scope = "col", "Detail"))),
            tags$tbody(rows)))),
      uiOutput("river_detail"))
  })

  .river_sel <- reactiveVal(NULL)
  observeEvent(input$river_pick, { .river_sel(input$river_pick) })

  output$river_detail <- renderUI({
    k <- .river_sel(); if (is.null(k) || !nzchar(k)) return(NULL)
    p <- strsplit(k, "|", fixed = TRUE)[[1]]
    if (length(p) != 3L) return(NULL)
    e <- .sk_data; if (is.null(e) || !nrow(e)) return(NULL)
    hit <- (e$ctxf == p[1] & e$modg == p[2] & e$tierl == p[3]) |
           (e$modg == p[1] & e$tierl == p[2] & e$tierl == p[3])
    if (p[2] %in% e$tierl) hit <- e$modg == p[1] & e$tierl == p[2]
    hit[is.na(hit)] <- FALSE
    g <- sort(unique(e$gene[hit]))
    l <- if (is.null(.g2l)) character(0) else
         sort(unique(.g2l$loc[.g2l$gene %in% g]))
    div(class = "interp", role = "region", tabindex = "-1",
        `aria-label` = "Supporting genes and loci for the selected flow",
      div(class = "dc-kicker", paste0(p[1], " to ", p[2], ", tier ", p[3])),
      tags$p(class = "interp-t", paste0(
        format(length(g), big.mark = ","), if (length(g) == 1L) " gene" else " genes",
        " across ", format(length(l), big.mark = ","),
        if (length(l) == 1L) " locus." else " loci.")),
      tags$p(class = "interp-f", paste(g, collapse = ", ")),
      if (length(l)) tags$p(class = "interp-f", paste(l, collapse = ", ")))
  })

  observeEvent(input$trans_locus, {
    if (!isTRUE(input$tn_sync)) return(invisible())
    v <- input$trans_locus
    if (is.null(v) || !nzchar(v)) return(invisible())
    if (!v %in% unlist(trans_loci_choices)) return(invisible())
    updateSelectInput(session, "tn_loc", selected = v)
  })

  observeEvent(input$tn_sync, {
    if (!isTRUE(input$tn_sync)) return(invisible())
    v <- input$trans_locus
    if (!is.null(v) && nzchar(v) && v %in% unlist(trans_loci_choices))
      updateSelectInput(session, "tn_loc", selected = v)
  }, ignoreInit = TRUE)

  observeEvent(input$tn_loc, {
    v <- input$tn_loc; t <- input$trans_locus
    if (isTRUE(input$tn_sync) && !is.null(v) && !is.null(t) && !identical(v, t))
      updateCheckboxInput(session, "tn_sync", value = FALSE)
  }, ignoreInit = TRUE)

  output$tn_syncnote <- renderUI({
    tl <- input$trans_locus; fl <- input$tn_loc
    drawable <- !is.null(tl) && nzchar(tl) && tl %in% unlist(trans_loci_choices)
    if (isTRUE(input$tn_sync) && !is.null(tl) && !drawable)
      return(div(class = "ax-status", role = "status",
        paste0("The table locus ", tl, " has no distal targets, so nothing can be drawn for it. ",
               "The figure below shows ", if (is.null(fl)) "another locus" else fl, ".")))
    if (!is.null(tl) && !is.null(fl) && !identical(tl, fl))
      return(div(class = "ax-status", role = "status",
        paste0("Figure locus ", fl, ". Table and export locus ", tl, ".")))
    div(class = "ax-status", role = "status",
        "Figure, table and export are on the same locus. The figure draws up to 40 distal targets.")
  })

  output$trans_bars <- renderPlot({
    loc <- input$tn_loc
    req(!is.null(loc), nzchar(loc), !is.null(trans_pairs))
    d <- trans_pairs[trans_pairs$ADlocus == loc, , drop = FALSE]
    req(nrow(d) > 0)
    d <- unique(d[, c("tgt", "modality", "tgt_chr")])
    n_by <- sort(table(d$tgt), decreasing = TRUE)
    top  <- names(n_by)[seq_len(min(40L, length(n_by)))]
    d    <- d[d$tgt %in% top, , drop = FALSE]
    d$tgt <- factor(d$tgt, levels = rev(top))
    d$modality <- factor(d$modality, levels = names(TMOD_COL)[names(TMOD_COL) %in% d$modality])
    ggplot(d, aes(y = tgt, fill = modality)) +
      geom_bar(width = 0.72) +
      scale_fill_manual(values = TMOD_COL, name = "Assay", drop = TRUE) +
      scale_x_continuous(breaks = function(z) seq(0, ceiling(max(z)), by = 1)) +
      labs(x = "Number of assays supporting the pair", y = NULL) +
      g_theme +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(size = 9))
  })

  .dk <- function(h, f = 0.74) {
    m <- grDevices::col2rgb(h)
    grDevices::rgb(t(round(m * f)), maxColorValue = 255)
  }

  .bars <- function(v, cols = NULL, lab_full = NULL, txtcol = NULL) {
    if (is.null(v) || !length(v)) return(div(class = "dc-note", "No records."))
    mx <- max(v)
    div(class = "rb",
      lapply(seq_along(v), function(i) {
        nm  <- names(v)[i]
        col <- if (!is.null(cols) && nm %in% names(cols)) unname(cols[[nm]]) else "#5b87b5"
        tx  <- if (!is.null(txtcol) && nm %in% names(txtcol)) unname(txtcol[[nm]]) else "#ffffff"
        pc  <- 100 * v[i] / mx
        lab <- if (!is.null(lab_full) && nm %in% names(lab_full)) lab_full[[nm]] else nm
        tagList(
          div(class = "rb-l", lab),
          div(class = "rb-b", title = sprintf("%s: %s", lab, format(v[i], big.mark = ",")),
            span(style = sprintf("width:%.2f%%;background:%s;color:%s", pc, col, tx),
                 if (pc > 11) format(v[i], big.mark = ",") else ""),
            span(class = "rb-hov", format(v[i], big.mark = ","))))
      }))
  }

  output$sum_tier <- renderUI({
    pal <- c(conf_pal[TIER_SEQ], "No tier" = "#cfd6db")
    names(pal) <- c(TIER_SEQ, "No tier")
    txt <- c(conf_txt[TIER_SEQ], "No tier" = "#233947")
    names(txt) <- c(TIER_SEQ, "No tier")
    .bars(locus_tier_counts, pal, NULL, txt)
  })

  output$sum_mod <- renderUI(.bars(mod_counts))

  output$sum_trans <- renderUI({
    d <- trans_all
    if (is.null(d) || !nrow(d))
      return(div(class = "dc-note", "No trans pairs recorded in this release."))
    mods <- c("snRNA", "pQTL", "gpQTL", "Hotspot")
    cols <- c(snRNA = "#2a78d6", pQTL = "#da532c",
              gpQTL = "#17868f", Hotspot = "#a8791b")
    tot <- sort(vapply(split(d$target, d$locus), function(z) length(unique(z)), integer(1)), decreasing = TRUE)
    top <- names(tot)[seq_len(min(12, length(tot)))]
    mx  <- max(as.numeric(tot[top]))
    rows <- lapply(top, function(loc) {
      sub <- d[d$locus == loc, , drop = FALSE]
      n <- vapply(mods, function(m) length(unique(sub$target[sub$modality == m])), numeric(1))
      ttot <- as.numeric(tot[[loc]])
      barw <- 100 * log1p(ttot) / log1p(mx)
      div(class = "sb-row",
          div(class = "sb-l dc-mono", loc),
          div(class = "sb-track",
              lapply(mods, function(m) {
                if (n[[m]] <= 0) return(NULL)
                pc <- barw * n[[m]] / sum(n)
                span(class = "sb-seg",
                     style = sprintf("width:%.3f%%;background:%s", pc, .dk(cols[[m]])),
                     if (pc > 7) format(n[[m]], big.mark = ",") else "",
                     span(class = "sb-bub",
                          sprintf("%s · %s distal target%s at %s", m,
                                  format(n[[m]], big.mark = ","),
                                  ifelse(n[[m]] == 1, "", "s"), loc)))
              })),
          div(class = "sb-n dc-mono", format(ttot, big.mark = ",")))
    })
    tagList(
      div(class = "sb-key",
          lapply(mods, function(m)
            span(class = "sb-ki",
                 span(class = "sb-sw", style = sprintf("background:%s", .dk(cols[[m]]))), m))),
      div(class = "sb-wrap", rows))
  })


  observeEvent(input$go_genes, nav_select("nav", "Genes", session = session))
  observeEvent(input$go_loci,  nav_select("nav", "Loci",  session = session))

  # ---- Cell types ------------------------------------------------------
  cty_d <- reactive({
    d  <- ct_rows(input$cty_ct)
    tr <- input$cty_tier
    if (!is.null(tr) && length(tr))
      d <- d[!is.na(d$top_confidence) & as.character(d$top_confidence) %in% tr, , drop = FALSE]
    qq <- trimws(if (is.null(input$cty_q)) "" else input$cty_q)
    if (nzchar(qq)) {
      pt  <- stringr::regex(qq, ignore_case = TRUE)
      hit <- stringr::str_detect(paste(d$gene, d$ADlocus, d$rsid), pt)
      hit[is.na(hit)] <- FALSE
      d <- d[hit, , drop = FALSE]
    }
    d
  })

  .cnt <- function(n, l) div(div(class = "dc-count-n", format(n, big.mark = ",")),
                             div(class = "dc-count-l", l))

  output$cty_stats <- renderUI({
    d <- cty_d()
    div(class = "dc-counts",
        .cnt(length(unique(d$ADlocus)), "AD loci"),
        .cnt(length(unique(d$gene[!is.na(d$gene) & nzchar(d$gene)])), "genes"),
        .cnt(nrow(d), "gene records"),
        .cnt(length(unique(na.omit(d$context))), "contexts"))
  })

  .tier_bar <- function(d) {
    if (is.null(d) || !nrow(d)) return(NULL)
    tb <- as.data.frame(table(factor(as.character(d$top_confidence), levels = .lv)),
                        stringsAsFactors = FALSE)
    names(tb) <- c("tier", "n")
    tb$tier <- factor(tb$tier, levels = .lv)
    ggplot(tb, aes(x = tier, y = n, fill = tier)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = ifelse(n > 0, format(n, big.mark = ","), "")),
                vjust = -0.55, size = 3.6, colour = "#233947") +
      scale_fill_manual(values = conf_pal, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_blank(),
            plot.background  = element_rect(fill = "#ffffff", colour = NA),
            panel.background = element_rect(fill = "#ffffff", colour = NA))
  }

  output$cty_tier_plot <- renderPlot(.tier_bar(cty_d()), res = 96)

  .rec_tbl <- function(d) {
    cc <- ctx_cols(d$context)
    datatable(
    tibble(Gene = d$gene, Locus = d$ADlocus, rsID = d$rsid,
           Tier = as.character(d$top_confidence),
           `Cell type or region` = cc$ctx, Dataset = cc$dset, Modality = cc$mod,
           `No. contexts` = d$n_contexts,
           `GWAS PIP` = round(d$max_inclusion, 3),
           cV2F = round(d$cv2f_score, 3)),
    escape = FALSE, rownames = FALSE,
    options = list(dom = "lrtip", pageLength = 25, scrollX = TRUE),
    class = "display compact hover")
  }

  output$cty_tbl <- renderDT(.rec_tbl(cty_d()), server = TRUE)

  output$dl_cty <- downloadHandler(
    filename = function() sprintf("AD_loci_%s.csv", gsub("[^A-Za-z0-9]+", "_", input$cty_ct)),
    content  = function(f) .prov_csv(cty_d(), f, "cell type evidence rows"))

  # ---- Batch lookup ----------------------------------------------------
  bt_tok <- reactiveVal(character(0))
  observeEvent(input$bt_go,    bt_tok(batch_tokens(input$bt_in)))
  observeEvent(input$bt_clear, {
    updateTextAreaInput(session, "bt_in", value = "")
    bt_tok(character(0))
  })
  observeEvent(input$bt_in, {
    if (length(bt_tok()) && !identical(batch_tokens(input$bt_in), bt_tok())) bt_tok(character(0))
  }, ignoreInit = TRUE)
  observeEvent(input$bt_demo, {
    updateTextAreaInput(session, "bt_in",
      value = "BIN1\nAPOE\nTREM2\nCLU\nPICALM\nrs429358\nchr2:127088929-127118906\nENSG00000130203")
  })

  bt_m <- reactive({ tk <- bt_tok(); if (!length(tk)) NULL else batch_match(tk) })
  bt_r <- reactive({ tk <- bt_tok(); if (!length(tk)) dat[0, ] else batch_rows(tk) })

  output$bt_summary <- renderUI({
    m <- bt_m()
    if (is.null(m)) return(div(class = "dc-note",
      "Nothing matched yet. Paste a list and press Match list."))
    r <- bt_r()
    div(class = "panel", style = "padding:16px 18px",
      div(class = "kv", span("Entries given"),    tags$b(nrow(m))),
      div(class = "kv", span("Found here"),       tags$b(sum(m$matched))),
      div(class = "kv", span("Not in this release"), tags$b(sum(!m$matched))),
      div(class = "kv", span("AD loci reached"),  tags$b(length(unique(r$ADlocus)))),
      div(class = "kv", span("Gene records"),     tags$b(nrow(r))))
  })

  output$bt_body <- renderUI({
    if (is.null(bt_m())) return(NULL)
    batch_body_ui()
  })

  output$bt_query_tbl <- renderDT({
    m <- bt_m(); req(!is.null(m))
    datatable(
      tibble(Entry = m$query, `Read as` = m$kind,
             Found = ifelse(m$matched, "yes", "no"),
             `AD loci` = m$n_loci, Genes = m$n_genes,
             `Best tier` = m$best_tier, `Loci matched` = m$loci),
      escape = FALSE, rownames = FALSE,
      options = list(dom = "lrtip", pageLength = 15, scrollX = TRUE),
      class = "display compact hover")
  }, server = TRUE)

  output$bt_tier_plot <- renderPlot(.tier_bar(bt_r()), res = 96)
  output$bt_rows_tbl  <- renderDT(.rec_tbl(bt_r()), server = TRUE)

  output$dl_batch <- downloadHandler(
    filename = function() "AD_loci_batch_match.csv",
    content  = function(f) .prov_csv(bt_r(), f, "batch list matches"))


  # ---- Release quality -------------------------------------------------
  .ql_theme <- function() theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill = "#ffffff", colour = NA),
          panel.background = element_rect(fill = "#ffffff", colour = NA),
          axis.title = element_text(size = 9.5, colour = "#5c6b75"))

  output$ql_stats <- renderUI({
    n   <- nrow(dat)
    pip <- dat$max_inclusion
    sharp <- sum(!is.na(pip) & pip >= 0.5)
    gw  <- sum(!is.na(dat$significance) & dat$significance == "genome wide")
    cv  <- sum(!is.na(dat$cv2f_score))
    div(class = "dc-counts", style = "margin-top:20px",
        .cnt(n, "gene records"),
        .cnt(round(100 * sharp / n), "% reach PIP 0.5 or more"),
        .cnt(round(100 * gw / n), "% genome-wide significant"),
        .cnt(round(100 * cv / n), "% carry a cV2F score"))
  })

  output$ql_pip <- renderPlot({
    v <- dat$max_inclusion[!is.na(dat$max_inclusion)]
    ggplot(data.frame(v = v), aes(x = v)) +
      geom_histogram(bins = 26, fill = "#2a78d6", colour = "#ffffff", linewidth = 0.4) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(x = "Highest inclusion probability", y = NULL) +
      .ql_theme()
  }, res = 96)

  output$ql_cv2f <- renderPlot({
    v <- dat$cv2f_score[!is.na(dat$cv2f_score)]
    ggplot(data.frame(v = v), aes(x = v)) +
      geom_histogram(bins = 26, fill = "#17868f", colour = "#ffffff", linewidth = 0.4) +
      labs(x = "cV2F score", y = NULL) +
      .ql_theme()
  }, res = 96)

  output$ql_sig <- renderPlot({
    nice <- c(`genome wide` = "Genome-wide", suggestive = "Suggestive",
              ns = "Not significant", `no p-value` = "No GWAS p-value")
    v <- as.character(dat$significance)
    v[is.na(v) | !nzchar(v)] <- "no p-value"
    keys <- unique(c(names(nice), sort(unique(v))))
    lab_of <- function(k) ifelse(k %in% names(nice), unname(nice[k]), k)
    tb <- as.data.frame(table(factor(v, levels = keys)), stringsAsFactors = FALSE)
    names(tb) <- c("k", "n")
    tb <- tb[tb$n > 0, , drop = FALSE]
    tb$lab <- factor(lab_of(tb$k), levels = lab_of(keys))
    ggplot(tb, aes(x = lab, y = n)) +
      geom_col(width = 0.6, fill = "#7a8b96") +
      geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.55,
                size = 3.3, colour = "#233947") +
      scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(x = NULL, y = NULL) +
      .ql_theme() +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(size = 8.5, lineheight = 0.95))
  }, res = 96)


  # ---- what is currently filtering each tab ----------------------------
  .chips <- function(items, clear_id) {
    if (!length(items)) return(NULL)
    div(class = "fl-act",
        span("Showing:"),
        lapply(items, function(x) span(class = "fl-chip", x)),
        actionLink(clear_id, "Clear all", class = "fl-clear"))
  }

  # ---- persistent current-view strip -------------------------------------
  squo <- intToUtf8(39)
  .cnt3 <- function(d) {
    if (is.null(d) || !is.data.frame(d) || !nrow(d)) return(NULL)
    if (!all(c("ADlocus", "gene") %in% names(d))) return(NULL)
    nl <- length(unique(d$ADlocus[!is.na(d$ADlocus)]))
    ng <- length(unique(d$gene[!is.na(d$gene) & nzchar(d$gene)]))
    paste0(format(nl, big.mark = ","), if (nl == 1L) " locus" else " loci", " · ",
           format(ng, big.mark = ","), if (ng == 1L) " gene" else " genes", " · ",
           format(nrow(d), big.mark = ","), if (nrow(d) == 1L) " record" else " records")
  }
  .tab_counts <- function(tab) {
    g <- function(e) tryCatch(e, error = function(x) NULL)
    if (identical(tab, "Cell types"))   return(.cnt3(g(cty_d())))
    if (identical(tab, "Genes"))        return(.cnt3(g(gp_rows())))
    if (identical(tab, "Batch"))        return(.cnt3(g(bt_r())))
    if (identical(tab, "Loci")) {
      lv <- g(loci_view())
      if (is.null(lv) || !nrow(lv)) return(NULL)
      return(.cnt3(dat[!is.na(dat$ADlocus) & dat$ADlocus %in% lv$ADlocus, , drop = FALSE]))
    }
    if (identical(tab, "Locus detail")) {
      lc <- g(selected_locus())
      if (is.null(lc) || !nzchar(lc)) return(NULL)
      return(.cnt3(dat[!is.na(dat$ADlocus) & dat$ADlocus == lc, , drop = FALSE]))
    }
    if (identical(tab, "Trans")) {
      r <- g(.tn_rows())
      if (is.null(r) || !is.data.frame(r)) return(NULL)
      return(paste0(format(nrow(r), big.mark = ","),
                    if (nrow(r) == 1L) " assay-specific trans record" else " assay-specific trans records"))
    }
    NULL
  }

  # ---- P2.14 / P2.15 provenance ------------------------------------------
  .rel_lab  <- if (exists("RELEASE_LABEL")) as.character(RELEASE_LABEL) else "2026-09"
  .bld_lab  <- if (exists("BUILD_LABEL")) as.character(BUILD_LABEL) else "GRCh38"
  .cite_txt <- "The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium. Broad and deep dissection of Alzheimer disease genetics with FunGen-xQTL."
  .cur_q <- function() {
    v <- tryCatch(session$clientData$url_search, error = function(e) "")
    if (is.null(v) || !nzchar(v)) "none (default view)" else sub("^[?]", "", v)
  }
  .prov_csv <- function(d, f, what) {
    ln <- c(paste0("# AD Loci Explorer export: ", what),
            paste0("# Release: ", .rel_lab),
            paste0("# Genome build: ", .bld_lab),
            paste0("# Generated: ", format(Sys.time(), tz = "America/New_York", usetz = TRUE)),
            paste0("# Query: ", .cur_q()),
            paste0("# Rows: ", nrow(d), "  Columns: ", ncol(d)),
            "# An empty cell means the value is not recorded in this release. It is not a zero",
            "# and not a tested null. TWAS, MR and cTWAS columns do distinguish a tested null.",
            "# Tiers T1 to T6 order evidence strength and are not statements of causality.",
            "# Column definitions are in the Documentation tab of the application.",
            paste0("# Cite: ", .cite_txt),
            "#")
    writeLines(ln, f)
    if (!is.null(d) && ncol(as.data.frame(d)) > 0)
      suppressWarnings(readr::write_csv(d, f, append = TRUE, col_names = TRUE))
  }
  .cite_url <- function() {
    cd <- session$clientData
    pr <- if (is.null(cd$url_protocol)) "https:" else cd$url_protocol
    hn <- if (is.null(cd$url_hostname)) "" else cd$url_hostname
    pt <- if (is.null(cd$url_port)) "" else as.character(cd$url_port)
    pn <- if (is.null(cd$url_pathname)) "/" else cd$url_pathname
    q  <- .cur_q()
    paste0(pr, "//", hn, if (nzchar(pt)) paste0(":", pt) else "", pn,
           if (nzchar(q)) paste0("?", q) else "")
  }

  .man_txt <- paste0(
    "The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium. ",
    "Broad and deep dissection of Alzheimer disease genetics with FunGen-xQTL. ",
    "Manuscript in preparation, ", substr(.rel_lab, 1, 4), ".")

  .rel_id <- function() {
    z <- REL_CITE$zenodo_doi; s <- REL_CITE$synapse_id
    if (nzchar(z)) paste0("https://doi.org/", z)
    else if (nzchar(s)) paste0("https://www.synapse.org/Synapse:", s)
    else ""
  }

  .rel_txt <- function() {
    id <- .rel_id()
    paste0("The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium. ",
           "FunGen-xQTL AD loci evidence release ", .rel_lab,
           " [data set]. ", REL_CITE$publisher, ", ", substr(.rel_lab, 1, 4), ". ",
           if (nzchar(id)) id else paste0("Release identifier not yet minted; cite this version and link: ", .cite_url()))
  }

  .bib_now <- function() {
    loc <- selected_locus(); if (is.null(loc)) loc <- ""
    ver <- gsub("[^A-Za-z0-9]", "", .rel_lab)
    z <- REL_CITE$zenodo_doi; s <- REL_CITE$synapse_id
    paste0(
      "@unpublished{fungen_xqtl_manuscript,\n",
      "  title  = {Broad and deep dissection of Alzheimer disease genetics with FunGen-xQTL},\n",
      "  author = {{The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium}},\n",
      "  year   = {", substr(.rel_lab, 1, 4), "},\n",
      "  note   = {Manuscript in preparation}\n}\n\n",
      "@misc{fungen_xqtl_release_", ver, ",\n",
      "  title     = {FunGen-xQTL AD loci evidence release ", .rel_lab, "},\n",
      "  author    = {{The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium}},\n",
      "  year      = {", substr(.rel_lab, 1, 4), "},\n",
      "  version   = {", .rel_lab, "},\n",
      "  publisher = {", REL_CITE$publisher, "},\n",
      if (nzchar(z)) paste0("  doi       = {", z, "},\n") else "",
      if (nzchar(s)) paste0("  howpublished = {Synapse ", s, "},\n") else "",
      "  url       = {", .cite_url(), "},\n",
      "  note      = {Genome build ", .bld_lab, ". Locus ", loc,
      ". Query: ", .cur_q(), ". Accessed ", format(Sys.Date()),
      ". Supported by NIH grant ", REL_CITE$grant, "}\n}")
  }

  .ris_now <- function() {
    loc <- selected_locus(); if (is.null(loc)) loc <- ""
    yr <- substr(.rel_lab, 1, 4)
    z <- REL_CITE$zenodo_doi; s <- REL_CITE$synapse_id
    man <- c("TY  - UNPB",
             "T1  - Broad and deep dissection of Alzheimer disease genetics with FunGen-xQTL",
             "AU  - The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium",
             paste0("PY  - ", yr),
             "N1  - Manuscript in preparation",
             "ER  - ", "")
    rel <- c("TY  - DATA",
             paste0("T1  - FunGen-xQTL AD loci evidence release ", .rel_lab),
             "AU  - The Alzheimer Disease Functional Genomics (FunGen-AD) Consortium",
             paste0("PY  - ", yr),
             paste0("PB  - ", REL_CITE$publisher),
             paste0("ET  - ", .rel_lab),
             if (nzchar(z)) paste0("DO  - ", z) else NULL,
             if (nzchar(s)) paste0("AN  - ", s) else NULL,
             paste0("UR  - ", .cite_url()),
             paste0("Y2  - ", format(Sys.Date())),
             paste0("N1  - Genome build ", .bld_lab, ". Locus ", loc,
                    ". Query: ", .cur_q(), ". Supported by NIH grant ", REL_CITE$grant),
             "ER  - ")
    paste(c(man, rel), collapse = "\n")
  }

  output$dl_cite_bib <- downloadHandler(
    filename = function() paste0("ad_loci_explorer_", format(Sys.Date()), ".bib"),
    content  = function(f) writeLines(.bib_now(), f))

  output$dl_cite_ris <- downloadHandler(
    filename = function() paste0("ad_loci_explorer_", format(Sys.Date()), ".ris"),
    content  = function(f) writeLines(.ris_now(), f))

  output$cite_box <- renderUI({
    loc <- selected_locus()
    if (is.null(loc) || !nzchar(loc)) return(NULL)
    j   <- match(loc, locus_registry$ADlocus)
    reg <- if (!is.na(j)) pretty_region(locus_registry$region[j]) else NA_character_
    acc <- format(Sys.Date())
    tags$details(class = "cite",
      tags$summary("Cite this result"),
      div(class = "dc-kicker", "Manuscript"),
      tags$p(class = "interp-f", .man_txt),
      div(class = "dc-kicker", "Data release"),
      tags$p(class = "interp-f", .rel_txt()),
      tags$ul(
        tags$li(paste0("Locus ", loc, if (!is.na(reg)) paste0("  ", reg) else "")),
        tags$li(paste0("Release ", .rel_lab, "  Genome build ", .bld_lab)),
        tags$li(paste0("Query: ", .cur_q())),
        tags$li(paste0("Accessed ", acc)),
        tags$li("Stable link for this exact view: ",
                tags$a(href = .cite_url(), .cite_url())),
        tags$li(paste0("Funding: NIH grant ", REL_CITE$grant))),
      div(class = "dc-exp", style = "margin:6px 0 8px",
        span(class = "dc-kicker", "Download citation"),
        downloadLink("dl_cite_bib", "BibTeX", class = "dc-btn"),
        downloadLink("dl_cite_ris", "RIS", class = "dc-btn")),
      tags$pre(class = "cite-bib", .bib_now()))
  })

  .vchip <- function(k, l) {
    span(class = "vs-chip", span(class = "vs-t", l),
         tags$button(class = "vs-x", type = "button",
                     `aria-label` = paste("Remove", l),
                     onclick = paste0("Shiny.setInputValue(", squo,
                                       "vs_drop", squo, ", ", squo, k, squo,
                                       ", {priority:", squo, "event", squo, "})"),
                     HTML("&times;")))
  }

  output$view_strip <- renderUI({
    tab <- input$nav
    G   <- function(e) tryCatch(e, error = function(x) NULL)
    it  <- list()
    add <- function(k, l) it[[length(it) + 1L]] <<- list(k = k, l = l)
    gv  <- function(k) { v <- G(input[[k]]); if (is.null(v)) character(0) else v }
    if (identical(tab, "Loci")) {
      ch <- gv("lt_chr")
      if (length(ch) && nzchar(ch[1])) add(paste0("lt_chr|", ch[1]), paste("chr", sub("^chr", "", ch[1])))
      for (v in gv("lt_ctx"))  add(paste0("lt_ctx|", v), v)
      for (v in gv("lt_tier")) add(paste0("lt_tier|", v), v)
      for (v in gv("lt_mod"))  add(paste0("lt_mod|", v), v)
      q <- G(trimws(lt_term())); if (is.null(q)) q <- ""
      if (nzchar(q)) add("lt_q|", paste0("matching ", q))
    } else if (identical(tab, "Cell types")) {
      for (v in gv("cty_tier")) add(paste0("cty_tier|", v), v)
      q <- gv("cty_q"); q <- if (length(q)) trimws(q[1]) else ""
      if (nzchar(q)) add("cty_q|", paste0("matching ", q))
    } else if (identical(tab, "Trans")) {
      for (v in gv("tmod")) add(paste0("tmod|", v), v)
      q <- gv("tn_q"); q <- if (length(q)) trimws(q[1]) else ""
      if (nzchar(q)) add("tn_q|", paste0("matching ", q))
    } else if (identical(tab, "Locus detail")) {
      for (v in gv("rv_ctx")) add(paste0("rv_ctx|", v), v)
      for (v in gv("rv_mod")) add(paste0("rv_mod|", v), v)
    }
    cnt <- G(.tab_counts(tab))
    if (!length(it) && is.null(cnt)) return(NULL)
    dot <- paste0(" ", intToUtf8(183), " ")
    scope <- if (identical(tab, "Loci") && !length(it))
      paste0("all chromosomes", dot, "all cell types", dot, "all tiers") else NULL
    div(class = "vs-wrap", role = "status",
        span(class = "vs-lab", "Current view"),
        if (!is.null(cnt)) span(class = "vs-n", cnt),
        if (!is.null(scope)) span(class = "vs-n", paste0(dot, scope)),
        lapply(it, function(z) .vchip(z$k, z$l)),
        actionLink("vs_clear", "Clear all", class = "vs-clear"))
  })

  observeEvent(input$vs_drop, {
    p <- strsplit(input$vs_drop, "|", fixed = TRUE)[[1]]
    id <- p[1]; val <- if (length(p) > 1L) p[2] else ""
    multi <- c("lt_ctx", "lt_tier", "lt_mod", "cty_tier", "rv_ctx", "rv_mod")
    if (match(id, multi, nomatch = 0L) > 0L) {
      updateSelectizeInput(session, id, selected = setdiff(input[[id]], val))
    } else if (identical(id, "tmod")) {
      updateCheckboxGroupInput(session, "tmod", selected = setdiff(input$tmod, val))
    } else if (identical(id, "lt_chr")) {
      updateSelectInput(session, "lt_chr", selected = "")
    } else if (identical(id, "lt_q")) {
      updateTextInput(session, "lt_q", value = ""); lt_term("")
    } else if (identical(id, "cty_q")) {
      updateTextInput(session, "cty_q", value = "")
    } else if (identical(id, "tn_q")) {
      updateTextInput(session, "tn_q", value = "")
    }
  })

  observeEvent(input$vs_clear, {
    tab <- input$nav
    if (identical(tab, "Loci")) {
      updateSelectInput(session, "lt_chr", selected = "")
      updateSelectizeInput(session, "lt_ctx",  selected = character(0))
      updateSelectizeInput(session, "lt_tier", selected = character(0))
      updateSelectizeInput(session, "lt_mod",  selected = character(0))
      updateTextInput(session, "lt_q", value = ""); lt_term("")
    } else if (identical(tab, "Cell types")) {
      updateSelectizeInput(session, "cty_tier", selected = character(0))
      updateTextInput(session, "cty_q", value = "")
    } else if (identical(tab, "Trans")) {
      updateCheckboxGroupInput(session, "tmod", selected = character(0))
      updateTextInput(session, "tn_q", value = "")
    } else if (identical(tab, "Locus detail")) {
      updateSelectizeInput(session, "rv_ctx", selected = character(0))
      updateSelectizeInput(session, "rv_mod", selected = character(0))
    }
  })

  output$lt_active <- renderUI({
    it <- character(0)
    if (!is.null(input$lt_chr) && nzchar(input$lt_chr))
      it <- c(it, paste("chromosome", sub("^chr", "", input$lt_chr)))
    if (length(input$lt_ctx))  it <- c(it, paste(input$lt_ctx,  collapse = " or "))
    if (length(input$lt_tier)) it <- c(it, paste(input$lt_tier, collapse = " or "))
    if (length(input$lt_mod))  it <- c(it, paste(input$lt_mod,  collapse = " or "))
    q <- trimws(lt_term())
    if (nzchar(q)) it <- c(it, paste0("matching ", q))
    .chips(it, "lt_reset")
  })
  observeEvent(input$lt_reset, {
    updateSelectInput(session, "lt_chr", selected = "")
    updateSelectizeInput(session, "lt_ctx",  selected = character(0))
    updateSelectizeInput(session, "lt_tier", selected = character(0))
    updateSelectizeInput(session, "lt_mod",  selected = character(0))
    updateTextInput(session, "lt_q", value = "")
    lt_term("")
  })

  output$cty_active <- renderUI({
    it <- character(0)
    if (length(input$cty_tier)) it <- c(it, paste(input$cty_tier, collapse = " or "))
    q <- trimws(if (is.null(input$cty_q)) "" else input$cty_q)
    if (nzchar(q)) it <- c(it, paste0("matching ", q))
    .chips(it, "cty_reset")
  })
  observeEvent(input$cty_reset, {
    updateSelectizeInput(session, "cty_tier", selected = character(0))
    updateTextInput(session, "cty_q", value = "")
  })

  output$tn_active <- renderUI({
    it <- character(0)
    if (identical(input$tn_scope, "locus") && !is.null(input$trans_locus) && nzchar(input$trans_locus))
      it <- c(it, paste("locus", input$trans_locus))
    n_now <- tryCatch(nrow(tn_tbl_d()), error = function(e) NA_integer_)
    if (identical(input$tn_scope, "locus") && !is.na(n_now) && n_now == 0L)
      it <- c(it, "no distal targets for this locus")
    if (length(input$tmod)) it <- c(it, paste(input$tmod, collapse = " or "))
    q <- trimws(if (is.null(input$tn_q)) "" else input$tn_q)
    if (nzchar(q)) it <- c(it, paste0("matching ", q))
    .chips(it, "tn_reset")
  })
  observeEvent(input$tn_reset, {
    updateRadioButtons(session, "tn_scope", selected = "all")
    updateCheckboxGroupInput(session, "tmod", selected = character(0))
    updateTextInput(session, "tn_q", value = "")
  })


  output$t6_n <- renderUI({
    n <- sum(!is.na(dat$top_confidence) & as.character(dat$top_confidence) == "T6")
    tags$b(format(n, big.mark = ","))
  })


  # ---- row clicks, gene deep link, locus compare ------------------------
  observeEvent(input$gp_row_locus, {
    l <- input$gp_row_locus
    req(!is.null(l), nzchar(l))
    selected_locus(l)
    nav_select("nav", "Locus detail", session = session)
  })

  .open_gene <- function(g) {
    if (is.null(g) || !nzchar(g)) return(invisible())
    updateSelectizeInput(session, "gp_gene", choices = genes, selected = g, server = TRUE)
    nav_select("nav", "Genes", session = session)
  }

  observeEvent(input$bt_rows_tbl_rows_selected, {
    i <- input$bt_rows_tbl_rows_selected; req(length(i))
    d <- bt_r(); req(nrow(d) >= i[1])
    .open_gene(as.character(d$gene[i[1]]))
  })
  observeEvent(input$cty_tbl_rows_selected, {
    i <- input$cty_tbl_rows_selected; req(length(i))
    d <- cty_d(); req(nrow(d) >= i[1])
    .open_gene(as.character(d$gene[i[1]]))
  })

  # deep links: gene, variant, region, filters and tab all route on first load
  observeEvent(session$clientData$url_search, {
    q <- shiny::parseQueryString(session$clientData$url_search)
    if (!length(q)) return(invisible())
    gv <- function(k) { v <- q[[k]]; if (is.null(v) || !nzchar(v)) NULL else trimws(v) }
    sp <- function(v) strsplit(v, ",", fixed = TRUE)[[1]]
    if (!is.null(gv("q")))     updateTextInput(session, "q", value = gv("q"))
    if (!is.null(gv("locus"))) updateSelectizeInput(session, "locus", selected = sp(gv("locus")))
    if (!is.null(gv("tier")))  updateSelectizeInput(session, "tier",  selected = sp(gv("tier")))
    if (!is.null(gv("sig")))   updateSelectizeInput(session, "sig",   selected = sp(gv("sig")))
    if (!is.null(gv("ct")))    updateSelectizeInput(session, "ct",    selected = sp(gv("ct")))
    if (!is.null(gv("minlp"))) updateSliderInput(session, "minlp",
          value = suppressWarnings(as.numeric(gv("minlp"))))
    if (!is.null(gv("transonly"))) updateCheckboxInput(session, "transonly",
          value = match(gv("transonly"), c("1", "true", "TRUE"), nomatch = 0L) > 0L)
    r <- tryCatch(resolve_locus(q), error = function(e) NULL)
    if (!is.null(r) && !is.null(r$locus)) selected_locus(r$locus)
    g <- gv("gene")
    if (!is.null(g)) {
      gu  <- toupper(g)
      hit <- genes[toupper(genes) == gu]
      if (!length(hit)) hit <- ens_to_symbol(gu)
      if (length(hit)) .open_gene(hit[1])
    }
    tb <- gv("tab")
    if (!is.null(tb)) nav_select("nav", tb, session = session)
    else if (!is.null(r) && !is.null(r$locus))
      nav_select("nav", "Locus detail", session = session)
  }, once = TRUE)

  # ---- compare two loci ------------------------------------------------
  .loc_facts <- function(loc) {
    d <- dat[!is.na(dat$ADlocus) & dat$ADlocus == loc, , drop = FALSE]
    e <- if (nrow(d)) parse_ctx_tokens(d) else NULL
    k <- match(as.character(d$top_confidence), TIER_SEQ)
    list(
      region = pretty_region(locus_registry$region[match(loc, locus_registry$ADlocus)]),
      best   = if (any(!is.na(k))) TIER_SEQ[min(k, na.rm = TRUE)] else NA_character_,
      genes  = length(unique(d$gene[!is.na(d$gene) & nzchar(d$gene)])),
      recs   = nrow(d),
      ctx    = if (is.null(e)) 0L else length(unique(e$ctx)),
      mods   = if (is.null(e)) 0L else length(unique(e$mod)),
      trans  = if (is.null(trans_all)) 0L else sum(trans_all$locus == loc),
      lead   = if (nrow(d)) as.character(d$rsid[which.max(replace(d$log10pval, is.na(d$log10pval), -Inf))]) else NA_character_,
      top    = if (nrow(d) && any(!is.na(k))) as.character(d$gene[which.min(k)]) else NA_character_)
  }

  output$cmp_tbl <- renderUI({
    a <- input$cmp_a; b <- input$cmp_b
    if (is.null(a) || is.null(b) || !nzchar(a) || !nzchar(b)) return(NULL)
    if (identical(a, b))
      return(div(class = "dc-note", "Pick two different loci."))
    fa <- .loc_facts(a); fb <- .loc_facts(b)
    tchip <- function(t) if (is.na(t)) HTML("&mdash;") else
      span(class = "tchip", style = sprintf("background:%s;color:%s",
           unname(conf_pal[[t]]), unname(conf_txt[[t]])), t)
    txt <- function(v) if (is.na(v) || !nzchar(as.character(v))) HTML("&mdash;") else as.character(v)
    row <- function(lab, va, vb, same_hint = TRUE) {
      hit <- same_hint && identical(as.character(va), as.character(vb))
      tags$tr(class = if (hit) "cmp-same" else NULL,
              tags$td(class = "cmp-l", lab), tags$td(va), tags$td(vb))
    }
    tags$table(class = "gp-tbl cmp-tbl",
      tags$thead(tags$tr(tags$th(""), tags$th(a), tags$th(b))),
      tags$tbody(
        row("Region", txt(fa$region), txt(fb$region), FALSE),
        row("Best tier", tchip(fa$best), tchip(fb$best), FALSE),
        row("Genes implicated", format(fa$genes, big.mark = ","), format(fb$genes, big.mark = ",")),
        row("Gene records", format(fa$recs, big.mark = ","), format(fb$recs, big.mark = ",")),
        row("Cell types", fa$ctx, fb$ctx),
        row("Assays", fa$mods, fb$mods),
        row("Distal links", format(fa$trans, big.mark = ","), format(fb$trans, big.mark = ",")),
        row("Strongest gene", txt(fa$top), txt(fb$top), FALSE),
        row("Lead variant", txt(fa$lead), txt(fb$lead), FALSE)))
  })


  # ---- trans figure exports --------------------------------------------
  .tn_rows <- reactive({
    loc <- input$tn_loc
    if (is.null(loc) || !nzchar(loc) || is.null(trans_all)) return(trans_all[0, ])
    trans_all[trans_all$locus == loc, , drop = FALSE]
  })

  output$dl_trans_csv <- downloadHandler(
    filename = function() sprintf("%s_trans_pairs.csv", input$tn_loc),
    content  = function(f) .prov_csv(.tn_rows(), f, "trans pairs for the drawn locus"))

  output$dl_trans_png <- downloadHandler(
    filename = function() sprintf("%s_trans_circos.png", input$tn_loc),
    content  = function(f) {
      grDevices::png(f, width = 1800, height = 1800, res = 200, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE)
      trans_chord(input$tn_loc)
    })

  output$dl_trans_svg <- downloadHandler(
    filename = function() sprintf("%s_trans_strands.svg", input$tn_loc),
    content  = function(f) {
      loc <- input$tn_loc
      d <- trans_pairs[trans_pairs$ADlocus == loc, , drop = FALSE]
      if (!nrow(d)) { writeLines("", f); return(invisible()) }
      i <- match(loc, locus_registry$ADlocus)
      srcreg <- if (!is.na(i)) pretty_region(locus_registry$region[i]) else loc
      d <- d[!duplicated(d$tgt), , drop = FALSE]
      svgs <- vapply(seq_len(nrow(d)), function(k)
        paste0(trans_illu_one(d$src_chr[k], d$src_pos[k], srcreg,
                              d$tgt_chr[k], d$tgt_pos[k], d$tgt[k],
                              as.character(d$modality[k]), d$src[k])),
        character(1))
      writeLines(svgs, f)
    })

}
