setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")

pkg <- c("ggplot2", "tidyverse", "data.table")
install_pkg(pkg,"common")
pkg <- c("ggbio","GenomicRanges","BSgenome.Hsapiens.UCSC.hg19")
install_pkg(pkg,"bioconductor")

karyogram_plot <- function(file_1, file_2, outpath, norm){
  tab1 <- fread(file_1, select = c(2,3,4,7))
  tab2 <- fread(file_2, select = c(2,3,4,7))
  summ <- tab1 %>% full_join(tab2, by=c("seqnames","start","end")) %>%
    dplyr::rename(RbKO = average.x,
                  WT = average.y) %>%
    drop_na()
  if(norm == TRUE) {
    summ <- summ %>% 
      mutate(RbKO = RbKO/median(RbKO),
              WT = WT/median(WT))
  }
    summ <- summ %>%
      mutate( FC = RbKO/WT) %>%
      filter(seqnames != "chrY") %>%
      filter(seqnames != "chrM")
  summ$seqnames <- factor(summ$seqnames, levels = c("chr1","chr2","chr3","chr4","chr5","chr6",
                                                    "chr7","chr8","chr9","chr10","chr11","chr12",
                                                    "chr13","chr14","chr15","chr16","chr17","chr18",
                                                    "chr19","chr20","chr21","chr22","chrX"))
  
  genome <- BSgenome.Hsapiens.UCSC.hg19
  genome <- seqlengths(genome)[1:23]
  gr <- GRanges(seqnames = summ$seqnames, range = IRanges(start = summ$start, end = summ$end), sig = summ$FC)
  seqlengths(gr) <- genome
  autoplot(gr,layout="karyogram", aes(color=log2(sig))) + scale_color_gradient2(midpoint=0, low="blue", mid="white",
                                                                                high="red", space ="Lab" )
  ggbio::ggsave(outpath,  dpi = 600, width = 8, height = 8)
}

karyogram_plot("100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv",
               "100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv",
               "../plot/cpdFC_karyogram_plot.pdf", TRUE)

karyogram_plot("10kb/RbKO_H3K9me3.fc.signal.bigwig_binned10kb.csv",
               "10kb/WT_H3K9me3.fc.signal.bigwig_binned10kb.csv",
               "../plot/H3K9me3FC_karyogram_plot.pdf", TRUE)

karyogram_plot("10kb/RbKO_H3K27me3.fc.signal.bigwig_binned10kb.csv",
               "10kb/WT_H3K27me3.fc.signal.bigwig_binned10kb.csv",
               "../plot/H3K27me3FC_karyogram_plot.pdf", TRUE)
