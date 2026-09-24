suppressMessages({library(data.table); library(readxl)})
d <- fread("shiny_app/data.csv"); b <- function(v) is.na(v)|v==""
sg <- function(v) v %in% c(TRUE,"TRUE","True","Yes",1,"1")
nt <- d[!b(gene) & b(top_confidence)]
g36 <- unique(nt[sg(twas_sig)|sg(mr_sig)|sg(ctwas_sig)]$gene)
message("T6-missing genes: ", length(g36))
f <- "out_20260917_fix/unified_AD_loci_xQTL_summary_20260917_fix.xlsx"
for (s in excel_sheets(f)) {
  x <- suppressMessages(read_excel(f, sheet=s))
  gc <- intersect(c("gene","gene_ID","gene_name","Gene"), names(x))
  if (!length(gc)) { message(s, " : no gene column"); next }
  gv <- unique(as.character(x[[gc[1]]]))
  message(s, " | rows=", nrow(x), " genes=", length(gv),
          " | of the 36 present: ", length(intersect(g36, gv)))
}
