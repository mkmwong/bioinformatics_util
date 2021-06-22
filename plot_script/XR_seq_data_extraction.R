setwd("~/Desktop/yr5/summ/bioinformatics_util/data/unbinned/XR-Seq/")
source("../../../draft.R")
pkg <- c("tidyverse")
install_pkg(pkg,"common")
pkg <- c("rtracklayer", "BSgenome.Hsapiens.UCSC.hg19")
install_pkg(pkg,"bioconductor")

file_names <- list.files(".")
new_names <- paste0(tools::file_path_sans_ext(file_names),"_100kb")
#  separate(file_names, c("name"))
gen <- BSgenome.Hsapiens.UCSC.hg19
si.gen <- seqinfo(gen)
for(i in 1:24) {
  print(i)
  dat <- import.bw(con=file_names[i], as="RleList")
  dat[dat==0] <- NA
  si <- si.gen[names(dat)]
  bins <- tileGenome(si,tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)
  binning_GR(bins, dat, TRUE, new_names[i] )
}


#for ( i in seq(1,24,4)) {
#  tab1 <- import.bw("GSM1985845_NHF1CPD_1h_Rep1_MINUS_UNIQUE_NORM_fixedStep_25.bw")
#  tab2 <- import.bw("GSM1985845_NHF1CPD_1h_Rep1_PLUS_UNIQUE_NORM_fixedStep_25.bw")
#  tab3 <- import.bw("GSM1985846_NHF1CPD_1h_Rep2_MINUS_UNIQUE_NORM_fixedStep_25.bw")
#  tab4 <- import.bw("GSM1985846_NHF1CPD_1h_Rep2_PLUS_UNIQUE_NORM_fixedStep_25.bw")
#}
