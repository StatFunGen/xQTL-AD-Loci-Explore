suppressMessages({library(data.table); library(readxl)})
d <- fread("shiny_app/data.csv"); b <- function(v) is.na(v)|v==""
sg <- function(v) v %in% c(TRUE,"TRUE","True","Yes",1,"1")
nt <- d[!b(gene) & b(top_confidence)]
g36 <- unique(nt[sg(twas_sig)|sg(mr_sig)|sg(ctwas_sig)]$gene)
f <- "out_20260917_fix/unified_AD_loci_xQTL_summary_20260917_fix.xlsx"
x <- suppressMessages(read_excel(f, sheet=1, skip=1))
gc <- grep("gene", names(x), ignore.case=TRUE, value=TRUE)
message("gene cols: ", paste(gc, collapse=" | "))
if (length(gc)) {
  gv <- unique(as.character(x[[gc[1]]]))
  message("rows=", nrow(x), " genes=", length(gv))
  message("of the 36 present: ", length(intersect(g36, gv)))
  message("absent: ", length(setdiff(g36, gv)))
}
