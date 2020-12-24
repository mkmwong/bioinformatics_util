#### Corresponding to fig.1
source("~/Desktop/work_tmp/bioinformatics_util/draft.R")
setwd("~/Desktop/work_tmp/bioinformatics_util/100kb/")

pkg <- c("dplyr","tidyr")
install_pkg(pkg,"common")

rbko <- fread("ML-Rb_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  rename(RbKO_Sig = average.x, WT_Sig = average.y) %>%
  mutate(WT_Sig = WT_Sig/median(na.omit(WT_Sig)), 
         RbKO_Sig = RbKO_Sig/median(na.omit(RbKO_Sig))) %>% 
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  mutate_at(vars(RbKO_Sig, WT_Sig), log2) %>% 
  drop_na() 

ylabel = bquote(atop(bold("DNA Lesions (RbKO)"),log2(FC)))
xlabel = bquote(atop(bold("DNA Lesions (WT)"),log2(FC)))

### Use without label to see where to place the pval/r label
corr_plot(summ$WT_Sig, summ$RbKO_Sig, xlabel, ylabel, "~/Desktop/scatter", marg = 0.1,
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=FALSE, x1=NULL, x2 = NULL, y1 = NULL, y2 = NULL)

### Now place the label onto desired location
corr_plot(summ$WT_Sig, summ$RbKO_Sig, xlabel, ylabel, "~/Desktop/scatter_2", marg = 0.1,
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=-1.5, x2 = -0.5, y1 = 2, y2 = 1.5)

