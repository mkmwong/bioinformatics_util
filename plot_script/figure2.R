##figure2
setwd("~/Desktop/yr5/summ/bioinformatics_util/data/100kb/")
source("../../draft.R")
install_pkg(c("ggpubr","ggplot2", "tidyr", "data.table","rstatix"), "common")
set.seed(123)

get_df <- function(rbko, wt, gap_type, gap_width, group_name){
  summ <- rbko %>% full_join(wt, by=c("seqnames", "start","end")) %>%
    dplyr::rename(RbKO = average.x, WT = average.y) %>%
    mutate(WT = WT/median(na.omit(WT)), 
           RbKO = RbKO/median(na.omit(RbKO)),
           FC = RbKO/WT) %>% 
    filter(seqnames != "chrX") %>%
    filter(seqnames != "chrY") %>%
    filter(seqnames != "chrM") %>%
    mutate_at(vars(RbKO, WT, FC), log2) %>% 
    drop_na() 
  
  gap <- fread("../hg19gap.txt", stringsAsFactors = FALSE)
  centroGR <- get_gap_GR(gap, gap_type, gap_width)
  summGR <- GRanges(seqnames = summ$seqnames, ranges = IRanges(start = summ$start, end = summ$end), 
                            RbKO = summ$RbKO, WT = summ$WT, FC = summ$FC)
  centro_regionGR <- subsetByOverlaps(summGR,centroGR)
  centro_region <- as.data.frame(centro_regionGR) %>%
    mutate(group = group_name)
  centro_region$seqnames <- as.character(centro_region$seqnames)
  rand_loc <- sample(1:nrow(summ),nrow(centro_region),replace=F)
  rand_region <- summ[rand_loc,] %>% 
    mutate(group = "Random")
  for (i in 1:999) {
    tmp_loc <- sample(1:nrow(summ),nrow(centro_region),replace=F)
    tmp_region <- summ[tmp_loc,] %>% 
      mutate(group = "Random")
    rand_region <- rbind(rand_region, tmp_region)
  }
  
  all_df <- rbind(centro_region[,c(1:3,6:9)], rand_region) %>%
    mutate(group = as.factor(group)) %>% 
    gather(Sig_type, Sig, -c(seqnames, start, end, group ))
  all_df$Sig_type <- factor(all_df$Sig_type, levels = c("RbKO","WT","FC"))
  return(all_df)
}

plot_boxplot <- function(all_df, outpath, y.position=NULL, ylab = NULL,rotate = FALSE, colors = NULL,
                        de = "pdf", w = 5, h = 3.5, dp = 600) {
  stat_test <- all_df %>%
    group_by(Sig_type) %>%
    wilcox_test(Sig ~ group, exact = FALSE, paired = FALSE, 
                p.adjust.method = "BH") %>%
    add_significance() %>% 
    mutate(group1 = group1, group2 = group2) %>% 
    add_xy_position(x = "group")
  bold_text <- element_text(face = "bold")
  all_df <- all_df %>% 
    group_by(Sig_type, group) %>% 
    filter(Sig < (quantile(Sig,0.75)-quantile(Sig,0.25))*1.5) %>% 
    filter(Sig> (quantile(Sig,0.75)-quantile(Sig,0.25))*-1.5)
  a = ggplot(data = all_df, aes(x = group, y = Sig, fill = group)) + 
    geom_boxplot( outlier.color = NA) +
    facet_wrap(~Sig_type, scale = "free", nrow = 1) + theme_bw()  + 
    xlab("Location") + ylab(ylab) + scale_fill_discrete(name = "Group") + 
    theme(axis.title = bold_text, legend.title = bold_text, strip.text = bold_text ) + 
    stat_pvalue_manual(stat_test, inherit.aes = FALSE, y.position = y.position )  
  if (rotate == TRUE) {
    a = a + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  if (!is.null(colors)) {
    a = a + scale_fill_manual(values = colors)
  }
  ggsave(outpath, device = de, width = w, height = h, dpi = dp)
}

rbko <- fread("ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "centromere", 1000000, "Pericentric")
ylabel <- "CPD log2(FC)"
plot_boxplot(all_df, "../../plot/CPD_Centromere.pdf", c(0.8,0.8,0.45), ylabel, TRUE, colors = c("darkorchid", "white"))

rbko <- fread("ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "centromere", 1000000, "Pericentric")
ylabel <- "H3K9me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K9me3_Centromere_new.pdf", c(1.4,1.4,0.4), ylabel, TRUE, colors = c("darkorchid", "white"))

rbko <- fread("ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "centromere", 1000000, "Pericentric")
ylabel <- "H3K27me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K27me3_Centromere_new.pdf", c(1.9,1.9,0.5), ylabel, TRUE, colors = c("darkorchid", "white"))

rbko <- fread("ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "telomere", 100000, "Subtelomeric") 
all_df$group = factor(all_df$group, levels = c("Subtelomeric","Random"))
ylabel <- "CPD log2(FC)"
plot_boxplot(all_df, "../../plot/CPD_Telomere.pdf", c(0.7,0.7, 0.45), ylabel, TRUE, colors = c("plum2", "white"))

rbko <- fread("ML-RbKO_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "telomere", 100000, "Subtelomeric")
all_df$group = factor(all_df$group, levels = c("Subtelomeric","Random"))
ylabel <- "H3K9me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K9me3_Telomere_new.pdf", c(1.6,1.6,0.4), ylabel, TRUE , colors = c("plum2", "white"))

rbko <- fread("ML-RbKO_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K27me3_new.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "telomere", 100000, "Subtelomeric")
all_df$group = factor(all_df$group, levels = c("Subtelomeric","Random"))
ylabel <- "H3K27me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K27me3_Telomere_new.pdf", c(1.9,1.9,0.5), ylabel, TRUE , colors = c("plum2", "white"))

#### Excluded on final version
rbko <- fread("ML-RbKO_H3K9me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K9me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "centromere", 1000000, "Pericentric")
ylabel <- "H3K9me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K9me3_Centromere_old.pdf", c(1.4,1.2,0.5), ylabel, TRUE)

rbko <- fread("ML-RbKO_H3K27me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K27me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "centromere", 1000000, "Pericentric")
ylabel <- "H3K27me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K27me3_Centromere_old.pdf", c(1.5,1.2,0.5), ylabel, TRUE)

rbko <- fread("ML-RbKO_H3K9me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K9me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "telomere", 100000, "Subtelomeric")
all_df$group = factor(all_df$group, levels = c("Subtelomeric","Random"))
ylabel <- "H3K9me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K9me3_Telomere_old.pdf", c(1.4,1.2,0.45), ylabel, TRUE)

rbko <- fread("ML-RbKO_H3K27me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
wt <- fread("ML-WT_H3K27me3_old.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4,7))
all_df <- get_df(rbko, wt, "telomere", 100000, "Subtelomeric")
all_df$group = factor(all_df$group, levels = c("Subtelomeric","Random"))
ylabel <- "H3K27me3 log2(FC)"
plot_boxplot(all_df, "../../plot/H3K27me3_Telomere_old.pdf", c(1.7,1.8,0.5), ylabel, TRUE)

