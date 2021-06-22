setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")

pkg <- c("ggplot2","data.table","tidyverse")
install_pkg(pkg,"common")
pkg <- c("DESeq2","IHW","biomaRt","vsn")
install_pkg(pkg,"bioconductor")

read_binned_repeats <- function(name, col) {
  tmp <- fread(name) %>%
    group_by_(col) %>% 
    drop_na() %>%
    summarize(mean = mean(average))
  return(tmp)
}

plot_7b <- function(col, outpath, w) {
  if(col == "V12") {
    order_sel = "repeat_class"
  } else {
    order_sel = "repeat_family"
  }
  H3K9ac <- read_binned_repeats("binned_repeats/E017-H3K9ac.fc.signal.bigwig_binnedrepeat.csv",col)
  H3K9me3 <- read_binned_repeats("binned_repeats/E017-H3K9me3.fc.signal.bigwig_binnedrepeat.csv",col)
  WT_H3K9me3 <- read_binned_repeats("binned_repeats/ML-WT_H3K9me3_new.fc.signal.bigwig_binnedrepeat.csv",col)
  RbKO_H3K9me3 <- read_binned_repeats("binned_repeats/ML-RbKO_H3K9me3_new.fc.signal.bigwig_binnedrepeat.csv",col)
  WT_H3K27me3 <- read_binned_repeats("binned_repeats/ML-WT_H3K27me3_new.fc.signal.bigwig_binnedrepeat.csv",col)
  RbKO_H3K27me3 <- read_binned_repeats("binned_repeats/ML-RbKO_H3K27me3_new.fc.signal.bigwig_binnedrepeat.csv",col)
  
  ordering <- summ %>% group_by_(order_sel) %>% summarize(medianfc = median(log2FoldChange)) 
  print(ordering)
  colnames(ordering)[1] = "repeat_t"
  summ2 <- WT_H3K9me3 %>%
    full_join(RbKO_H3K9me3,col) %>%
    full_join(WT_H3K27me3,col) %>%
    full_join(RbKO_H3K27me3,col) %>%
    full_join(H3K9ac, col) %>%
    full_join(H3K9me3, col) %>%
    dplyr::rename(WT_H3K9me3 = mean.x,
                  RbKO_H3K9me3 = mean.y,
                  WT_H3K27me3 = mean.x.x,
                  RbKO_H3K27me3 = mean.y.y,
                  Roadmap_H3K9ac = mean.x.x.x,
                  Roadmap_H3K9me3 = mean.y.y.y,
                  repeat_t = col) %>%
    mutate(FC_H3K9me3 = RbKO_H3K9me3/WT_H3K9me3,
           FC_H3K27me3 = RbKO_H3K27me3/WT_H3K27me3) %>%
    #### here is to normalize by genome median??(100kb bins)
    mutate(Roadmap_H3K9ac = Roadmap_H3K9ac/0.9884352,
           Roadmap_H3K9me3 = Roadmap_H3K9me3/1.042893,
           WT_H3K9me3 = WT_H3K9me3/0.7524252,
           WT_H3K27me3 = WT_H3K27me3/0.7545076, 
           RbKO_H3K9me3 = RbKO_H3K9me3/0.7408611,
           RbKO_H3K27me3 = RbKO_H3K27me3/0.6983734) %>%
    full_join(ordering, "repeat_t") %>%
    drop_na() %>%
    gather(type, Sig, -medianfc, -repeat_t) %>% 
    mutate(Sig = log2(Sig)/median(Sig))
  
  summ2$type <- factor(summ2$type, levels=c("Roadmap_H3K9ac", "Roadmap_H3K9me3",
                                           "RbKO_H3K9me3", "WT_H3K9me3", "FC_H3K9me3",
                                           "RbKO_H3K27me3","WT_H3K27me3", "FC_H3K27me3"))
  summ2 <- summ2 %>% 
    filter(repeat_t %in% selected) %>%
    mutate(repeat_t = replace(repeat_t, repeat_t=="Other","SVA"),
           repeat_t = replace(repeat_t, repeat_t=="centr","Centromere"),
           repeat_t = replace(repeat_t, repeat_t=="acro","Acromere"),
           repeat_t = replace(repeat_t, repeat_t=="telo","Telomere"))
  
  ggplot(summ2, aes(reorder(repeat_t,medianfc), type)) + geom_tile(aes(fill = Sig)) +
    scale_fill_gradient2( low = "blue", high = "red", mid = "white", midpoint = 0) + 
    xlab("Types of repeats") + ylab("Histone Modification") +labs (fill="Signal") + 
    theme_bw() + theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))  
  ggsave(outpath,width=w, height=3, dpi=600, units="in", device="pdf")
}

#### plot 7a ####
summ <- fread("repenrich_out/lane2_TCAG_L002_R1_ALL/lane2_TCAG_L002_R1_ALL_fraction_counts.txt",  stringsAsFactors = FALSE)
tmp <- fread("repenrich_out/lane3_GCTT_L003_R1_ALL/lane3_GCTT_L003_R1_ALL_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp1 <- fread("repenrich_out/SCGPM_DAK-CPD-01_H7H52_L1_ACAGTG_R1/SCGPM_DAK-CPD-01_H7H52_L1_ACAGTG_R1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp2 <- fread("repenrich_out/SCGPM_DAK-CPD-01_H7H52_L1_GCCAAT_R1/SCGPM_DAK-CPD-01_H7H52_L1_GCCAAT_R1_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)

summ <- summ %>% 
  full_join(tmp,by="V1") %>%
  full_join(tmp1,by="V1") %>%
  full_join(tmp2, by="V1") %>% 
  tibble::column_to_rownames("V1") %>%
  dplyr::rename(RB1_1 = V4.x,
         RB1_2 = V4.y,
         IMR_1 = V4.x.x,
         IMR_2 = V4.y.y)

coldata <- as.data.frame(matrix(ncol = 1, nrow = 4))
colnames(coldata) <- c("type")
coldata$type <- c("RbKO","RbKO","WT","WT")
rownames(coldata) <- c("RB1_1","RB1_2","IMR_1","IMR_2")
tmp1 <- run_deseq(summ, coldata, "type", 10, 'WT', "../plot/fig7pca.pdf")

tmp1 <- as.data.frame(tmp1)
tmp1$repeats <- rownames(tmp1)
tmp <- fread("repenrich_out/lane3_GCTT_L003_R1_ALL/lane3_GCTT_L003_R1_ALL_fraction_counts.txt", select=(c(1,2,3))) %>%
  dplyr::rename(repeats = V1,
         repeat_class = V2,
         repeat_family = V3) %>% 
  full_join(tmp1) %>% 
  drop_na()
summ <- tmp
summ$dir = "unchanged"
summ$dir[which(summ$pvalue <= 0.05 & summ$log2FoldChange >= 0)] = "more"
summ$dir[which(summ$pvalue <= 0.05 & summ$log2FoldChange <= 0)] = "less"
write.csv(summ,file="../fig7_deseq2testout.csv")

#### removing classes with ? after its name
summ <- summ %>% filter(!grepl("\\?", repeat_family)) 
ggplot(summ,aes(x= reorder(repeat_class,log2FoldChange, median), y = log2FoldChange, fill = repeat_class )) + 
  geom_boxplot() + xlab("Repeat class") + ylab("log2FC(RbKO/WT CPDSeq)") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none", 
        axis.title = element_text(face = "bold"))
ggsave("../plot/cpd_repenrich_classFC.pdf", width = 5, height = 5, dpi = 600)

ggplot(summ,aes(x= reorder(repeat_family,log2FoldChange, median), y = log2FoldChange, fill = repeat_class )) + 
  geom_boxplot() + xlab("Repeat family") + ylab("log2FC(RbKO/WT CPDSeq)") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(face = "bold"))
ggsave("../plot/cpd_repenrich_familyFC.pdf", width = 15, height = 5, dpi = 600)

#### plot 7b ####
col <- "V13"
col <- "V12"

plot_7b("V13", "../plot/fig7b_family.pdf", 12)
plot_7b("V12", "../plot/fig7b_class.pdf", 5)


#### pick some important ones to plot
selected <- c("L1", "PiggyBac","TcMar", "LTR", "ERV1", "Alu","Other","ERVK","L2","Satellite","centr","telo")
summ1 <- summ %>% 
  filter(repeat_family %in% selected) %>%
  mutate(repeat_family = replace(repeat_family, repeat_family=="Other","SVA"),
         repeat_family = replace(repeat_family, repeat_family=="centr","Centromere"),
         repeat_family = replace(repeat_family, repeat_family=="acro","Acromere"),
         repeat_family = replace(repeat_family, repeat_family=="telo","Telomere"))
ggplot(summ1,aes(x= reorder(repeat_family,log2FoldChange, median), y = log2FoldChange, fill = repeat_class )) + 
  geom_boxplot() + xlab("Repeat family") + ylab("log2FC(RbKO/WT CPDSeq)") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) + labs(fill = "Repeat class") 
ggsave("../plot/cpd_repenrich_familyFC_subset.pdf", width = 6, height = 5, dpi = 600)
plot_7b("V13", "../plot/fig7b_family_subset.pdf", 6)

#### plot 7c #### 
a <- fread("100kb/E017-H3K9me3.fc.signal.bigwig_binned100kb.csv", select=c(2,3,4,7))
b <- fread("100kb/ML-WT_H3K9me3_new.fc.signal.bigwig_binned100kb.csv", select=c(2,3,4,7))
summ <- a %>%
  full_join(b, by=c("seqnames","start","end")) %>%
  dplyr::rename(Encode_H3K9me3 = average.x,
                Morr_H3K9me3 = average.y) %>%
  drop_na()

corr_plot(summ$Encode_H3K9me3, summ$Morr_H3K9me3, "Encode_H3K9me3", "MorrisonLab_H3K9me3", 
          "../plot/fig7c", marg = 0, w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
          addLab=TRUE, x1=1, x2 = 1, y1 = 2.8, y2 = 2.5)

