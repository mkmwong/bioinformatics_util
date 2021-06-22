setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")

pkg <- c("ggplot2", "data.table", "tidyr", "ggsignif")
install_pkg(pkg, "common")
pkg <- c("rtracklayer")
install_pkg(pkg, "bioconductor")
plot_exp_obs_overlap <- function(peakGR, topGR, bottomGR, outpath, title, ypos) {
  
  topOL = subsetByOverlaps(peakGR, topGR)
  length(topOL)
  bottomOL = subsetByOverlaps(peakGR, bottomGR)
  length(bottomOL)
  df = as.data.frame(matrix(nrow = 6, ncol = 3))
  colnames(df) = c("val","group","oe")
  df$val[1:3] = c(length(topOL), length(bottomOL), (length(peakGR@seqnames) - length(bottomOL) - length(topOL) ))
  df$oe[1:3] = "Observed"
  df$val[4:6] = c(round(length(peakGR@seqnames) * 0.1), round(length(peakGR@seqnames) * 0.1), round(length(peakGR@seqnames) * 0.8))
  df$oe[4:6] = "Expected"
  df$group = factor(rep(c("top10","bottom10","other"),2), levels = c("top10","other","bottom10"))
  
  p = chisq.test(df$val[1:3] ,p = c(0.1,0.8,0.1))$p.val
  print(p)
  ggsave( filename = outpath,
          ggplot(data = df, aes(x = oe, y = val, fill = group, label = val)) + geom_bar(stat = "identity") + 
            xlab("Group") + ylab("Peak Count") + scale_fill_manual(values = c("#CC6677","#888888","#44AA99")) +
            guides(fill=guide_legend(title="Genomic Region")) + geom_text(size = 4, position = position_stack(vjust = 0.5)) +
            ggtitle(title) + geom_signif(y_position=ypos, xmin=1, xmax=2 ,annotation=c("***"), tip_length=0) + 
            theme_bw(), 
          width=4, height=4, dpi=600, units="in", device="pdf")
  # annotate("text", x = 1.5, y = ypos, label = "***")
}

quies <- fread("data/SRR034480.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 >= -log10(0.05))
grow <- read.table("data/SRR034478.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 > -log10(0.05))
growGR <- GRanges(seqnames =  grow$V1, ranges = IRanges(start = grow$V2 + 1, end = grow$V3))
quiesGR <- GRanges(seqnames =  quies$V1, ranges = IRanges(start = quies$V2 + 1, end = quies$V3))

tab1 <- fread("data/100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
tab2 <- fread("data/100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- tab1 %>% full_join(tab2, by=c("seqnames","start","end")) %>%
  dplyr::rename(RbKO = average.x,
                WT = average.y) %>%
  drop_na() %>% 
    mutate(RbKO = RbKO/median(RbKO),
           WT = WT/median(WT)) %>%
  mutate( FC = log2(RbKO/WT)) %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM")

top <- summ %>% 
  filter(FC > quantile(FC,0.9))
bottom <- summ %>% 
  filter(FC < quantile(FC,0.1))
topGR <- GRanges(seqnames = top$seqnames, ranges = IRanges(start = top$start, end = top$end))
bottomGR <- GRanges(seqnames = bottom$seqnames, ranges = IRanges(start = bottom$start, end = bottom$end))

plot_exp_obs_overlap(growGR, topGR, bottomGR, "plot/figure9/growing.pdf", "Growing Rb peaks", 1150)
plot_exp_obs_overlap(quiesGR, topGR, bottomGR, "plot/figure9/quiescent.pdf", "Quiescent Rb peaks",2350)

