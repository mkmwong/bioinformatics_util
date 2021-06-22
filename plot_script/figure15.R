setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("tidyverse","data.table","huxtable","rstatix","ggpubr")
install_pkg(pkg,"common")

signals <- readRDS("100kb/XR-Seq/merged/XRSeq.rds") %>%
  filter(seqnames != "chrM") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrX") 
colnames(signals)[4:9] = c("XR1h", "XR4h", "XR8h", "XR16h", "XR1d", "XR2d") 
signals <- signals %>% 
  mutate(XR2d = (XR2d + XR1d + XR16h + XR8h + XR4h + XR1h),
         XR1d = (XR1d + XR16h + XR8h + XR4h + XR1h),
         XR16h = (XR16h + XR8h + XR4h + XR1h),
         XR8h = (XR8h + XR4h + XR1h),
         XR4h = (XR4h + XR1h))

signalsGR <-  GRanges(seqnames = signals$seqnames, ranges = IRanges(start = signals$start, end = signals$end), 
                      XR1h = signals$XR1h, XR4h = signals$XR4h, XR8h = signals$XR8h, 
                      XR16h = signals$XR16h, XR1d = signals$XR1d, XR2d = signals$XR2d)

gap <- fread("hg19gap.txt")
centro_coor <- get_gap_GR(gap, "centromere", 1000000)
telo_coor <- get_gap_GR(gap, "telomere", 100000)

gap_signal <- as.data.frame(subsetByOverlaps(signalsGR,centro_coor)) %>%
  dplyr::mutate(group = "Pericentric") %>%
  dplyr::select(-width, -strand)

rand_loc <- sample(1:nrow(signals),nrow(gap_signal),replace=F)
rand_region <- signals[rand_loc,] %>% 
  mutate(group = "Random") 

for (i in 1:999) {
  tmp_loc <- sample(1:nrow(signals),nrow(gap_signal),replace=F)
  tmp_region <- signals[tmp_loc,] %>% 
    mutate(group = "Random") 
  rand_region <- rbind(rand_region, tmp_region)
}

summ <- na.omit(rbind(gap_signal, rand_region)) %>% 
  gather(key = "Treatment", value = "Signal", -c(seqnames,start, end, group))
summ$Treatment = factor(summ$Treatment, levels = c("XR1h", "XR4h", "XR8h", "XR16h" ,"XR1d", "XR2d" ))
stat_test <- summ %>%
  group_by(Treatment) %>%
  wilcox_test(Signal ~ group, exact = FALSE, paired = FALSE, 
              p.adjust.method = "BH") %>%
  add_significance() %>% 
  mutate(group1 = group1, group2 = group2) %>% 
  add_xy_position(x = "group")
bold_text <- element_text(face = "bold")
ggplot(summ, aes(x = group, y = Signal, fill = group)) + geom_boxplot(outlier.colour = NA) + 
   xlab("Group") + ylab("") + facet_wrap(~ Treatment, nrow=1) + ylim(0,100) + 
  theme_bw() + stat_pvalue_manual(stat_test, inherit.aes = FALSE , y.position = 95) +
  theme(axis.text.x = element_text(angle = 30, hjust =1)) + 
  scale_fill_manual(values = c("darkorchid","white"))
ggsave("../fig15v2_pericentric.pdf", w = 5, h = 4)

gap_signal <- as.data.frame(subsetByOverlaps(signalsGR,telo_coor)) %>%
  dplyr::mutate(group = "Subtelomeric") %>%
  dplyr::select(-width, -strand)
rand_loc <- sample(1:nrow(signals),nrow(gap_signal),replace=F)
rand_region <- signals[rand_loc,] %>% 
  mutate(group = "Random") 

for (i in 1:999) {
  tmp_loc <- sample(1:nrow(signals),nrow(gap_signal),replace=F)
  tmp_region <- signals[tmp_loc,] %>% 
    mutate(group = "Random") 
  rand_region <- rbind(rand_region, tmp_region)
}

summ <- na.omit(rbind(gap_signal, rand_region)) %>% 
  gather(key = "Treatment", value = "Signal", -c(seqnames,start, end, group))
summ$Treatment = factor(summ$Treatment, levels = c("XR1h", "XR4h", "XR8h", "XR16h" ,"XR1d", "XR2d" ))
summ$group = factor(summ$group, levels = c("Subtelomeric","Random"))
stat_test <- summ %>%
  group_by(Treatment) %>%
  wilcox_test(Signal ~ group, exact = FALSE, paired = FALSE, 
              p.adjust.method = "BH") %>%
  add_significance() %>% 
  mutate(group1 = group1, group2 = group2) %>% 
  add_xy_position(x = "group")
bold_text <- element_text(face = "bold")
ggplot(summ, aes(x = group, y = Signal, fill = group)) + geom_boxplot(outlier.colour = NA) + 
  xlab("Group") + ylab("") + facet_wrap(~ Treatment, nrow=1) + ylim(0,100) + 
  theme_bw() + stat_pvalue_manual(stat_test, inherit.aes = FALSE , y.position = 95) +
  theme(axis.text.x = element_text(angle = 30, hjust =1)) + 
  scale_fill_manual(values = c("plum2","white"))
ggsave("../fig15v2_subtelomeric.pdf", w = 5, h = 4)
