setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("tidyverse","data.table","huxtable")
install_pkg(pkg,"common")

rbko <- fread("100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  mutate(WT = WT/median(na.omit(WT)), 
         RbKO = RbKO/median(na.omit(RbKO)),
         FC = RbKO/WT) %>% 
  dplyr::select(seqnames, start, end, FC) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(FC), log2) 

rbko_hist <- fread("100kb/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt_hist <- fread("100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ_hist <- rbko_hist %>% full_join(wt_hist, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  mutate(WT = WT/median(na.omit(WT)), 
         RbKO = RbKO/median(na.omit(RbKO)),
         FC = RbKO/WT) %>% 
  dplyr::select(seqnames, start, end, FC) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(FC), log2) 
summ_plot <- summ %>% full_join(summ_hist, by=c("seqnames", "start","end")) %>%
  drop_na()

ylabel = bquote(atop(bold("DNA Lesions"),log2(RbKO/WT)))
xlabel = bquote(atop(bold("H3K9me3"),log2(RbKO/WT)))

corr_plot(summ_plot$FC.x, summ_plot$FC.y, xlabel, ylabel, 
          "../plot/figure4_H3K9me3", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=0.5, x2 = 0.5, y1 = -0.6, y2 = -0.8,
          xmin = -0.5, xmax = 0.75, ymin = -1, ymax = 0.75)

rbko_hist <- fread("100kb/ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt_hist <- fread("100kb/ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ_hist <- rbko_hist %>% full_join(wt_hist, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  mutate(WT = WT/median(na.omit(WT)), 
         RbKO = RbKO/median(na.omit(RbKO)),
         FC = RbKO/WT) %>% 
  dplyr::select(seqnames, start, end, FC) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(FC), log2) 

summ_plot <- summ %>% full_join(summ_hist, by=c("seqnames", "start","end")) %>%
  drop_na() 

ylabel = bquote(atop(bold("DNA Lesions"),log2(RbKO/WT)))
xlabel = bquote(atop(bold("H3K27me3"),log2(RbKO/WT)))

corr_plot(summ_plot$FC.x, summ_plot$FC.y, ylabel, 
          xlabel, "../plot/figure4_H3K27me3", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=0.5, x2 = 0.5, y1 = 0.8, y2 = 0.5,
          xmin = -0.5, xmax = 0.75, ymin = -1, ymax = 1)


  
