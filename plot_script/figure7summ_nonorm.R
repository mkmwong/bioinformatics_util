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

#### plot 7supp ####
summ <- fread("repenrich_out/RbKOH3K9me3_1_out/RbKOH3K9me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp <- fread("repenrich_out/RbKOH3K9me3_2_out/RbKOH3K9me3_2_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp1 <- fread("repenrich_out/WTH3K9me3_1_out/WTH3K9me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp2 <- fread("repenrich_out/WTH3K9me3_2_out/WTH3K9me3_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)

summ <- fread("repenrich_out/RbKOH3K27me3_1_out/RbKOH3K27me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp <- fread("repenrich_out/RbKOH3K27me3_2_out/RbKOH3K27me3_2_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp1 <- fread("repenrich_out/WTH3K27me3_1_out/WTH3K27me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp2 <- fread("repenrich_out/WTH3K27me3_2_out/WTH3K27me3_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)

plot_7supp <- function(in1, in2, in3, in4, out1, out2, out3) {
  summ <- fread(in1, select=c(1,4), stringsAsFactors = FALSE)
  tmp <- fread(in2, select=c(1,4), stringsAsFactors = FALSE)
  tmp1 <- fread(in3, select=c(1,4), stringsAsFactors = FALSE)
  tmp2 <- fread(in4,  select=c(1,4), stringsAsFactors = FALSE)
  
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
  write.csv(summ,file=out1)
  summ <- summ %>% filter(!grepl("\\?", repeat_family)) 
  ggplot(summ,aes(x= reorder(repeat_family,log2FoldChange, median), y = log2FoldChange, fill = repeat_class )) + 
    geom_boxplot() + xlab("Repeat family") + ylab("log2FC(RbKO/WT H3K9me3)") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none",
          axis.title = element_text(face = "bold"))
  ggsave(out2, width = 15, height = 5, dpi = 600)
  
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
  ggsave(out3, width = 6, height = 5, dpi = 600)
}
plot_7supp("repenrich_out/RbKOH3K27me3_1_out/RbKOH3K27me3_1_fraction_counts.txt", 
           "repenrich_out/RbKOH3K27me3_2_out/RbKOH3K27me3_2_fraction_counts.txt", 
           "repenrich_out/WTH3K27me3_1_out/WTH3K27me3_1_fraction_counts.txt", 
           "repenrich_out/WTH3K27me3_2_out/WTH3K27me3_2_fraction_counts.txt", 
           "../fig7k27me3_deseq2testout.csv", 
           "../plot/k27_repenrich_familyFC.pdf", 
           "../plot/k27_repenrich_familyFC_subset.pdf")
plot_7supp("repenrich_out/RbKOH3K9me3_1_out/RbKOH3K9me3_1_fraction_counts.txt", 
           "repenrich_out/RbKOH3K9me3_2_out/RbKOH3K9me3_2_fraction_counts.txt", 
           "repenrich_out/WTH3K9me3_1_out/WTH3K9me3_1_fraction_counts.txt", 
           "repenrich_out/WTH3K9me3_2_out/WTH3K9me3_2_fraction_counts.txt", 
           "../fig7k9me3_deseq2testout.csv", 
           "../plot/k9_repenrich_familyFC.pdf", 
           "../plot/k9_repenrich_familyFC_subset.pdf")
