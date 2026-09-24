suppressMessages(library(data.table))
d <- fread("shiny_app/data.csv"); t <- fread("out_20260917_fix/gene_tier_assignment.csv")
b <- function(v) is.na(v) | v==""
nb <- d[b(gene)]
message("BLANK GENE rows: ", nrow(nb))
message("  with gene_id   : ", nb[!b(gene_id), .N])
message("  with variant_ID: ", nb[!b(variant_ID), .N])
message("  distinct loci  : ", uniqueN(nb$ADlocus))
nt <- d[!b(gene) & b(top_confidence)]
g <- unique(nt$gene)
message("GENE-NO-TIER rows: ", nrow(nt), " genes: ", length(g))
message("  in tier file : ", length(intersect(g, t$gene_name)))
message("  not in it    : ", length(setdiff(g, t$gene_name)))
message("EVIDENCE_TIER set while top_confidence blank: ", d[b(top_confidence) & !b(evidence_tier), .N])
print(table(top=ifelse(b(d$top_confidence),"blank","set"), ev=ifelse(b(d$evidence_tier),"blank","set")))
message("sample untiered genes:")
print(head(sort(g), 18))
