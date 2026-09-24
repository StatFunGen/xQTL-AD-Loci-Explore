# =============================================================================
# FunGen-xQTL · AD Loci Explorer
# Dashboard landing (KPI value boxes + leaderboard plots) + filter-driven
# explorer with definition tooltips & downloads + locus detail + trans + glossary.
# Reads data.csv (sitting beside this file). Deploy with rsconnect::deployApp().
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(plotly)
library(circlize)
library(shinycssloaders)
library(scales)

# ---- module load order ----------------------------------------------------
for (.f in c("10_tokens", "20_data", "30_helpers", "40_registry", "45_views", "50_theme",
             "60_links", "70_css", "80_ui", "90_server")) {
  source(file.path("modules", paste0(.f, ".R")), local = FALSE)
}
rm(.f)

shinyApp(ui, server)
