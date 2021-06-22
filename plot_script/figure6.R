setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")

pkg <- c("ggplot2", "tidyverse", "data.table")
install_pkg(pkg,"common")
pkg <- c("ggbio","GenomicRanges","BSgenome.Hsapiens.UCSC.hg19")
install_pkg(pkg,"bioconductor")

build_sig_track <- function(file_1, file_2, norm) {
  tab1 <- fread(file_1, select = c(2,3,4,7))
  tab2 <- fread(file_2, select = c(2,3,4,7))
  summ <- tab1 %>% full_join(tab2, by=c("seqnames","start","end")) %>%
    dplyr::rename(RbKO = average.x,
                  WT = average.y) #%>%
    #drop_na()
  if(norm == TRUE) {
    summ <- summ %>% 
      mutate(RbKO = RbKO/median(na.omit(RbKO)),
             WT = WT/median(na.omit(WT)))
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
  return(gr)
}

karyogram_plot <- function(gr, outpath, h = 8){
  #tab1 <- fread(file_1, select = c(2,3,4,7))
  #tab2 <- fread(file_2, select = c(2,3,4,7))
  #summ <- tab1 %>% full_join(tab2, by=c("seqnames","start","end")) %>%
  #  dplyr::rename(RbKO = average.x,
  #                WT = average.y) %>%
  #  drop_na()
  #if(norm == TRUE) {
  #  summ <- summ %>% 
  #    mutate(RbKO = RbKO/median(RbKO),
  #           WT = WT/median(WT))
  #}
  #  summ <- summ %>%
  #    mutate( FC = RbKO/WT) %>%
  #    filter(seqnames != "chrY") %>%
  #   filter(seqnames != "chrM")
  #summ$seqnames <- factor(summ$seqnames, levels = c("chr1","chr2","chr3","chr4","chr5","chr6",
  #                                                  "chr7","chr8","chr9","chr10","chr11","chr12",
  #                                                  "chr13","chr14","chr15","chr16","chr17","chr18",
  #                                                  "chr19","chr20","chr21","chr22","chrX"))
  
  #genome <- BSgenome.Hsapiens.UCSC.hg19
  #genome <- seqlengths(genome)[1:23]
  #gr <- GRanges(seqnames = summ$seqnames, range = IRanges(start = summ$start, end = summ$end), sig = summ$FC)
  #seqlengths(gr) <- genome
  autoplot(gr,layout="karyogram", aes(color=sig, fill = sig)) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) + 
    theme_bw() + theme(legend.title = element_text(face = "bold"))
  ggbio::ggsave(outpath,  dpi = 600, width = 8, height = h)
}

a <- build_sig_track(file_1 = "100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv",
                     file_2 = "100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", TRUE)
a$sig = log2(a$sig)
a$sig[a$sig >= quantile(a$sig, 0.995, na.rm = TRUE)]= quantile(a$sig, 0.995, na.rm = TRUE)
a$sig[a$sig <= quantile(a$sig, 0.005, na.rm = TRUE)]= quantile(a$sig, 0.005, na.rm = TRUE)
karyogram_plot(a, "../plot/cpdFC_karyogram_plot.tiff", h = 6  )

b <- build_sig_track(file_1 = "10kb/RbKO_H3K9me3.fc.signal.bigwig_binned10kb.csv",
                     file_2 = "10kb/WT_H3K9me3.fc.signal.bigwig_binned10kb.csv", TRUE)
karyogram_plot(b, "../plot/H3K9me3FC_karyogram_plot.pdf")

c <- build_sig_track(file_1 = "10kb/RbKO_H3K27me3.fc.signal.bigwig_binned10kb.csv",
                     file_2 = "10kb/WT_H3K27me3.fc.signal.bigwig_binned10kb.csv", TRUE)
karyogram_plot(c, "../plot/H3K27me3FC_karyogram_plot.pdf" )

b <- build_sig_track(file_1 = "100kb/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv",
                     file_2 = "100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", TRUE)
b$sig = log2(b$sig)
karyogram_plot(b, "../plot/100kb_H3K9me3FC_karyogram_plot.tiff", h= 6 )

c <- build_sig_track(file_1 = "100kb/ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv",
                     file_2 = "100kb/ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", TRUE)
c$sig = log2(c$sig)
karyogram_plot(c, "../plot/100kb_H3K27me3FC_karyogram_plot.tiff",h =6 )

single_chrom_plot <- function(chr, outpath) {
  h = 2.5
  a <- build_sig_track(file_1 = "100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv",
                       file_2 = "100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", TRUE)
  b <- build_sig_track(file_1 = "100kb/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv",
                       file_2 = "100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", TRUE)
  c <- build_sig_track(file_1 = "100kb/ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv",
                       file_2 = "100kb/ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", TRUE)
  p.ideo <- Ideogram(genome = "hg19", xlabel = TRUE) + xlim(GRanges(chr, IRanges(1, 141213431))) 
  #a$sig = log2(a$sig)
  a$sig[a$sig >= quantile(a$sig, 0.995, na.rm = TRUE)]= quantile(a$sig, 0.995, na.rm = TRUE)
  a$sig[a$sig <= quantile(a$sig, 0.005, na.rm = TRUE)]= quantile(a$sig, 0.005, na.rm = TRUE)
  ap <- autoplot(a[which(a@seqnames== chr)], layout="linear", aes(color=log2(sig))) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab",  ) + 
    theme_void() + labs(color="log2(FC)") + theme(legend.direction = "horizontal", legend.position = "right",
                                                  legend.key.size = unit(0.2, "cm"), legend.text=element_text(size=3),
                                                  legend.title = element_text(size = 3))
  bp <- autoplot(b[which(b@seqnames== chr)], layout="linear", aes(color=log2(sig))) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab" ) + 
    theme_void() + labs(color="log2(FC)") + theme(legend.direction = "horizontal", legend.position = "right",
                                                  legend.key.size = unit(0.2, "cm"), legend.text=element_text(size=3),
                                                  legend.title = element_text(size = 3))
  cp <- autoplot(c[which(c@seqnames== chr)], layout="linear", aes(color=log2(sig))) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab" ) + 
    theme_void() + labs(color="log2(FC)") + theme(legend.direction = "horizontal", legend.position = "right",
                                                  legend.key.size = unit(0.2, "cm"), legend.text=element_text(size=3),
                                                  legend.title = element_text(size = 3))
  fixed(p.ideo) <- TRUE
  fixed(ap) <- TRUE
  fixed(bp) <- TRUE
  fixed(cp) <- TRUE
  tks <- tracks(p.ideo, `CPD FC` = ap, `H3K9me3 FC` = bp, `H3K27me3 FC` = cp,label.text.angle = 0, 
                label.width = unit(7, "lines"), heights = c(3,1,1,1), padding = unit(0., "lines"), xlim = IRanges(1, 141213431))
  tks
  ggbio::ggsave(outpath,  dpi = 600, width = 12, height = h)
  
}

single_chrom_plot("chr9", "../plot/chr9_plot.pdf" )
single_chrom_plot("chr1", "../plot/chr1_plot.pdf" )
