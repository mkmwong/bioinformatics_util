source("~/Desktop/yr5/summ/bioinformatics_util/draft.R")
setwd("~/Desktop/yr5/summ/bioinformatics_util/data/100kb/")

pkg <- c("dplyr","tidyr","data.table")
install_pkg(pkg,"common")

rbko <- fread("ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO_Sig = average.x, WT_Sig = average.y) %>%
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
corr_plot(summ$WT_Sig, summ$RbKO_Sig, xlabel, ylabel, "~/Desktop/yr5/summ/bioinformatics_util/plot/figure1", marg = 0.1,
          w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=-1.5, x2 = -1.5, y1 = 2, y2 = 1.5)


#### coloring of top and bottom

rbko <- fread("ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
  dplyr::rename(RbKO_Sig = average.x, WT_Sig = average.y) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  drop_na() 


summ <- summ %>% mutate(FC = RbKO_Sig/WT_Sig) 
summ$Group = "Other"
summ$Group[which(summ$FC >= quantile(summ$FC, 0.9))] = "Top 10%"
summ$Group[which(summ$FC <= quantile(summ$FC, 0.1))] = "Bottom 10%"
summ$Group <- factor(summ$Group, levels = c("Top 10%", "Bottom 10%", "Other"))

corre = paste("r =", round(cor.test(summ$RbKO_Sig, summ$WT_Sig, method = "spearman", exact = FALSE)$estimate,2))
corre.pval = pval_lab(cor.test(summ$RbKO_Sig, summ$WT_Sig)$p.val)

ggplot(data = summ, aes(x = log2(WT_Sig), y = log2(RbKO_Sig), color = Group)) + 
  geom_point(alpha = 0.2, size = 0.2) + xlab(xlabel) + ylab(ylabel) +
  theme_bw() + scale_color_manual(values = c("#CC6677","#44AA99","#888888")) + 
  annotate(geom="text",x = -1.5, y = 0.8, label=corre,size=4, fontface="bold") + 
  annotate(geom="text",x = -1.5, y = 0.5, label=corre.pval, size=4) +
  geom_smooth(color="blue",method = "lm", size = 0.5) + xlim(-2.2,1) + ylim(-2.2,1)
ggsave("../../plot/figure1v2.pdf", w =4.5, h = 3)


