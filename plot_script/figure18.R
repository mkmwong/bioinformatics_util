setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("tidyverse","data.table","huxtable")
install_pkg(pkg,"common")

rbko <- fread("100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(WT, RbKO), log2)  %>%
  drop_na()

ylabel = bquote(atop(bold("DNA Lesions"),log2(RbKO)))
xlabel = bquote(atop(bold("H3K9me3"),log2(RbKO)))

corr_plot(summ$RbKO, summ$WT, ylabel, 
          xlabel,  "../plot/RbKO_CPDvsH3K9me3", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=-1.5, x2 = -1.5, y1 = 0.5, y2 = 0.2,
          xmin = -2, xmax = 0.5, ymin = -2, ymax = 0.8)

rbko <- fread("100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()

ylabel = bquote(atop(bold("DNA Lesions"),log2(WT)))
xlabel = bquote(atop(bold("H3K9me3"),log2(WT)))

corr_plot(summ$RbKO, summ$WT, ylabel, 
          xlabel,  "../plot/WT_CPDvsH3K9me3", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=-1.5, x2 = -1.5, y1 = 1, y2 = 0.6,
          xmin = -2, xmax = 0.5, ymin = -2, ymax = 1.2)

rbko <- fread("100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()

ylabel = bquote(atop(bold("DNA Lesions"),log2(WT)))
xlabel = bquote(atop(bold("H3K27me3"),log2(WT)))

corr_plot(summ$RbKO, summ$WT, ylabel, 
          xlabel,  "../plot/WT_CPDvsH3K27me3", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=-1.5, x2 = -1.5, y1 = 1.8, y2 = 1.2,
          xmin = -2, xmax = 0.5, ymin = -2, ymax = 2.5)

rbko <- fread("100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()

ylabel = bquote(atop(bold("DNA Lesions"),log2(RbKO)))
xlabel = bquote(atop(bold("H3K27me3"),log2(RbKO)))

corr_plot(summ$RbKO, summ$WT, ylabel, 
          xlabel,  "../plot/RbKO_CPDvsH3K27me3", marg = 0, 
          w = 3 , h = 2.5, un = "in", d = 600, de = "tiff",
          addLab=FALSE, x1=-1.5, x2 = -1.5, y1 = 1.8, y2 = 1.2,
          xmin = -2, xmax = 0.5, ymin = -2.2, ymax = 2)

##### comparison purpose with old data ####
rbko <- fread("100kb/E017-H3K9me3.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(WT, RbKO), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "H3K9me3(Roadmap)", 
          "H3K9me3(This paper)",  "../plot/H3K9me3_roadmapvsourpaper", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=2, x2 = 2, y1 = -1, y2 = -1.4)

rbko <- fread("100kb/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(WT, RbKO), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "H3K9me3(RbKO)", 
          "H3K9me3(WT)",  "../plot/H3K9me3_wt_v_rbko", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=0.6, x2 = 0.6, y1 = -1, y2 = -1.4)


rbko <- fread("100kb/E017-H3K27me3.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(WT, RbKO), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "H3K27me3(Roadmap)", 
          "H3K27me3(This paper)",  "../plot/H3K27me3_roadmapvsourpaper", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=2, x2 = 2, y1 = -1, y2 = -1.4)

rbko <- fread("100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/E017-H3K9me3.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "CPD(WT)", 
          "H3K9me3(WT)",  "../plot/WT_CPDvsH3K9me3_roadmap", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=1, x2 = 1, y1 = 2, y2 = 2.5)

rbko <- fread("1mb/WT_CPD_fc.signal.bw_binned1mb.csv", select = c(2,3,4,7))
wt <- fread("1mb/E017-H3K9me3.fc.signal.bigwig_binned1mb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "CPD(WT)", 
          "H3K9me3(WT)",  "../plot/WT_CPDvsH3K9me3_roadmap_1mb", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=0.6, x2 = 0.6, y1 = 1.5, y2 = 2)

rbko <- fread("100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("100kb/E017-H3K27me3.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "CPD(WT)", 
          "H3K27me3(Roadmap)",  "../plot/WT_CPDvsH3K27me3_roadmap", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=0.6, x2 = 0.6, y1 = -1, y2 = -1.4)

rbko <- fread("1mb/WT_CPD100_fc_signal.bw_binned1mb.csv", select = c(2,3,4,7))
wt <- fread("1mb/E017-H3K27me3.fc.signal.bigwig_binned1mb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO = average.x, WT = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO, WT), log2)  %>%
  drop_na()
corr_plot(summ$RbKO, summ$WT, "CPD(WT)", 
          "H3K27me3(Roadmap)",  "../plot/WT_CPDvsH3K27me3_roadmap_1mb", marg = 0, 
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=0.6, x2 = 0.6, y1 = -1, y2 = -1.4)

##### comparison purpose with old data ####