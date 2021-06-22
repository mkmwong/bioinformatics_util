setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
pkg = c("GenomeGraphs")
install_pkg(pkg, "bioconductor")

tmp = fread("data/rb_geneplot_subset.bed") %>% 
  filter(V2 != "gene") %>% filter(V2 != "transcript") %>%
  filter(v2 != "exon")
