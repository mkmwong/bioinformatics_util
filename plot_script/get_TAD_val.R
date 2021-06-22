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

process_and_bin <- function(fname, bin, outpath, bin_num) {
  a <- mapply(make_range, bin$seqnames, bin$start, bin$end, bin_num, SIMPLIFY = FALSE)
  b <- do.call(rbind, a) %>% 
    mutate(id = rep(1:(length(a)), each = bin_num))
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

#######
get_coor_GR <- function(tmp, width){
  ret <- as.data.frame(matrix(ncol = 3, nrow = nrow(tmp)*2 )) %>% 
    mutate(V1 = rep(tmp$seqnames, 2),
           V2 = c((tmp$start-width + 1), tmp$end + 1),
           V3 = c(tmp$start, (tmp$end + width)),
           V4 = rep(1:nrow(tmp), each = 2),
           V5 = rep(c("L","R"), nrow(tmp)))
  retGR <- GRanges(seqnames= ret$V1, ranges = IRanges(start = ret$V2, end = ret$V3),
                   which_closed = ret$V4, which_side = ret$V5 )
  return(retGR)
}

closed_side <- as.data.frame(get_coor_GR(closed,20000))
left <- closed_side %>% filter(which_side == "L")
right <- closed_side %>% filter(which_side == "R")
process_and_bin("unbinned/RbKO_CPD100.fc.signal.bw", left, "RbKOCPD_closedTAD_left" ,25)
process_and_bin("unbinned/RbKO_CPD100.fc.signal.bw", right, "RbKOCPD_closedTAD_right")
process_and_bin("unbinned/WT_CPD100.fc.signal.bw", left, "WTCPD_closedTAD_left")
process_and_bin("unbinned/WT_CPD100.fc.signal.bw", right, "WTCPD_closedTAD_right")
process_and_bin("unbinned/RbKO_H3K27me3.fc.signal.bigwig", left, "RbKOH3K27me3_closedTAD_left")
process_and_bin("unbinned/RbKO_H3K27me3.fc.signal.bigwig", right, "RbKOH3K27me3_closedTAD_right")
process_and_bin("unbinned/WT_H3K27me3.fc.signal.bigwig", left, "WTH3K27me3_closedTAD_left")
process_and_bin("unbinned/WT_H3K27me3.fc.signal.bigwig", right, "WTH3K27me3_closedTAD_right")
process_and_bin("unbinned/RbKO_H3K9me3.fc.signal.bigwig", left, "RbKOH3K9me3_closedTAD_left")
process_and_bin("unbinned/RbKO_H3K9me3.fc.signal.bigwig", right, "RbKOH3K9me3_closedTAD_right")
process_and_bin("unbinned/WT_H3K9me3.fc.signal.bigwig", left, "WTH3K9me3_closedTAD_left")
process_and_bin("unbinned/WT_H3K9me3.fc.signal.bigwig", right, "WTH3K9me3_closedTAD_right")


