setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")

pkg <- c("tidyverse","data.table")
install_pkg(pkg,"common")
pkg <- c("rtracklayer")
install_pkg(pkg,"bioconductor")

hic_tab <- fread("hic_compartments_100kb_imr90_2014.txt", select = c(1,2,3,5))
tmp <- hic_tab %>% arrange(chr, domain, start) %>%
  group_by(domain, chr) %>%
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-n()])) %>%
  group_by(domain, indx, chr) %>%
  summarise(start = dplyr::first(start, order_by = indx), end = dplyr::last(end, order_by = indx))
hicGR <- GRanges(seqnames = hic_tab$chr, ranges = IRanges(start = hic_tab$start + 1, end = hic_tab$end + 1 ), domain = hic_tab$domain)
hic_compressed <- reduce(split(hicGR , ~domain))
hic_compressed$closed$domain <- "closed"
hic_compressed$open$domain <- "open"
closed <- as.data.frame(hic_compressed$closed)
open <- as.data.frame(hic_compressed$open)

process_and_bin <- function(fname, bin, outpath) {
  a <- mapply(make_range, bin$seqnames, bin$start, bin$end, SIMPLIFY = FALSE)
  b <- do.call(rbind, a) %>% 
    mutate(id = rep(1:(length(a)), each = 50))
  b <- GRanges(seqnames = b$chr, ranges = IRanges(start = b$st, end = b$en), id = b$id)
  dat <- import.bw(fname, as = "RleList")
  print(head(b))
  print(head(dat))
  dat[dat == 0] <- NA
  dat1 <- remove_chrom(dat,"chrX", column=NULL)
  dat1 <- remove_chrom(dat1,"chrY", column=NULL)
  dat1 <- remove_chrom(dat1,"chrM", column=NULL)
  o <- binning_GR(b, dat1, TRUE, outpath)
}

process_and_bin("unbinned/RbKO_CPD100.fc.signal.bw", closed, "RbKOCPD_closedTAD")
process_and_bin("unbinned/RbKO_CPD100.fc.signal.bw", open, "RbKOCPD_openTAD")
process_and_bin("unbinned/WT_CPD100.fc.signal.bw", closed, "WTCPD_closedTAD")
process_and_bin("unbinned/WT_CPD100.fc.signal.bw", open, "WTCPD_openTAD")
process_and_bin("unbinned/RbKO_H3K27me3.fc.signal.bigwig", closed, "RbKOH3K27me3_closedTAD")
process_and_bin("unbinned/RbKO_H3K27me3.fc.signal.bigwig", open, "RbKOH3K27me3_openTAD")
process_and_bin("unbinned/WT_H3K27me3.fc.signal.bigwig", closed, "WTH3K27me3_closedTAD")
process_and_bin("unbinned/WT_H3K27me3.fc.signal.bigwig", open, "WTH3K27me3_openTAD")
process_and_bin("unbinned/RbKO_H3K9me3.fc.signal.bigwig", closed, "RbKOH3K9me3_closedTAD")
process_and_bin("unbinned/RbKO_H3K9me3.fc.signal.bigwig", open, "RbKOH3K9me3_openTAD")
process_and_bin("unbinned/WT_H3K9me3.fc.signal.bigwig", closed, "WTH3K9me3_closedTAD")
process_and_bin("unbinned/WT_H3K9me3.fc.signal.bigwig", open, "WTH3K9me3_openTAD")




