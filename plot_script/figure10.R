setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")

pkg <- c("ggplot2", "data.table", "ggsignif")
install_pkg(pkg, "common")

plot_10b <- function( outpath, w, h) {
  WT_H3K27me3_Genome_med <- fread("data/binned_chrom_state/WT_H3K27me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8)) %>%
    mutate(stat = "GenomeMed")
  WT_H3K27me3<- fread("data/binned_chrom_state/WT_H3K27me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))
  WT_H3K27me3 <- rbind(WT_H3K27me3_Genome_med, WT_H3K27me3) %>% 
    group_by(stat) %>% 
    drop_na() %>% 
    mutate(average = average/0.8605049) %>%
    summarize(mean = mean(average))
  WT_H3K9me3_Genome_med <- fread("data/binned_chrom_state/WT_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8)) %>%
    mutate(stat = "GenomeMed")
  WT_H3K9me3<- fread("data/binned_chrom_state/WT_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))
  WT_H3K9me3 <- rbind(WT_H3K9me3_Genome_med, WT_H3K9me3) %>% 
    group_by(stat) %>% 
    drop_na() %>% 
    mutate(average = average/0.7406521) %>%
    summarize(mean = mean(average))
  
  ordering <- summ %>% group_by(stat) %>% summarize(medianfc = median(FC)) 
  print(ordering)
  colnames(ordering)[1] = "stat"
  summ2 <- WT_H3K9me3 %>%
    full_join(WT_H3K27me3, by=c("stat")) %>%
    dplyr::rename(WT_H3K9me3 =mean.x,
                  WT_H3K27me3 = mean.y) %>%
    select(WT_H3K9me3, WT_H3K27me3, stat) %>%
    full_join(ordering, "stat") %>%
    drop_na() %>%
    gather(type, Sig, -medianfc, -stat) %>% 
    mutate(Sig = log2(Sig)/median(Sig))
  
  summ2$type <- factor(summ2$type, levels=c("WT_H3K9me3", "WT_H3K27me3"))
  
  ggplot(summ2, aes(reorder(stat,medianfc, median), type)) + geom_tile(aes(fill = Sig)) +
    scale_fill_gradient2( low = "blue", high = "red", mid = "white", midpoint = 0) + 
    xlab("") + ylab("") + labs(fill="Signal") + theme_bw() + 
    theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1),
          axis.ticks.length.x = unit(0, "cm"), legend.key.size = unit(0.3, "cm"), 
          legend.text=element_text(size=8), legend.title = element_text(size = 8))  
  ggsave(outpath,width=w, height=h, dpi=600, units="in")
}
#### H3K9me3 
x_fill = c("#00FF56","#114662","#6FFFB9","#99C238","#F64844",
           "#FF00AB","#F07202","#00983A","#17D56D","#004A1E",
           "#FB81F1","#40A0CB","#CFE6F6","#999999","#F80000","#FEFEFE")

tmp <- fread("data/binned_chrom_state/RbKO_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))
tmp1 <- fread("data/binned_chrom_state/WT_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))

summ <- tmp %>% full_join(tmp1, by=c("seqnames", "start","end","stat")) %>%
  dplyr::rename(Rb_sig = average.x, WT_sig = average.y) %>%
  filter(seqnames != "chrM") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrX") %>% 
  drop_na() %>%
  mutate(Rb_sig = Rb_sig/median(Rb_sig)) %>%
  mutate(WT_sig = WT_sig/median(WT_sig)) %>% 
  mutate(FC = log2(Rb_sig/WT_sig))
summ2 <- summ
summ2$stat <- "GenomeMed"
summ <- rbind(summ,summ2)

ylabel = bquote(atop(bold("H3K9me3"),log2(FC)))
x_lab_col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","firebrick",
              "chartreuse4","chartreuse4","chartreuse4","firebrick","chartreuse4",
              "firebrick","chartreuse4", "Black","chartreuse4","chartreuse4",
              "firebrick")
x_lab = c("Bivalent TSS","Flanking Bivalent TSS/Enhancer","Bivalent Enhancer",
          "Zinc Finger","Repressed PolyComb","Transcription","Genic Enhancer",
          "Transcription 5' and 3'","Weak Repressed PolyComb","Active TSS","Heterochromatin",
          "Flanking TSS","Genome Median","Weak Transcription", "Enhancer","Quiescent")

all_state <- unique(summ$stat)
stat_sig <- as.data.frame(matrix(nrow = length(all_state), ncol = 2))
stat_sig$V1 <- all_state
for (i in 1:length(all_state)) {
  subsett = summ %>% 
    filter(stat == all_state[i])
  stat_sig$V2[i] <- wilcox.test(subsett$FC, summ2$FC)$p.val
}
stat_sig$V2 <- round(stat_sig$V2, 5)

ggplot(data = summ, aes(x = reorder(stat, FC, median), y = FC, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,position = "top") + scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), 
        axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 120, 30, 0)) + 
  geom_signif(y_position=c(2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4), 
              xmin=c(1,2,3,4,5,6,7,8,9,14,15,16), xmax=c(1,2,3,4,5,6,7,8,9,14,15,16),
              annotation=c("*"), tip_length=0, linetype = "blank")
ggsave("plot/H3K9me3_chrom_with_sig.pdf" , width=5.5, height=4, dpi=600, units="in")
plot_10b("plot/10b_H3K9me3.pdf", w = 5.5, h = 1.5)

### H3K27me3
## For H3K27me3
#order: 10,11,12,4,8,
#       6,3,13,5,MED
#       1,2,4,7,15,9
#labels = c("Active TSS","Flanking TSS","Transcription 5' and 3'",
#"Transcription","Weak Transcription","Genic Enhancer",
#"Enhancer","Zinc Finger","Heterochromatin","Bivalent TSS",
#"Flanking Bivalent TSS/Enhancer","Bivalent Enhancer","Repressed PolyComb",
#"Weak Repressed PolyComb","Quiescent/","Genome Median")

tmp <- fread("data/binned_chrom_state/RbKO_H3K27me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))
tmp1 <- fread("data/binned_chrom_state/WT_H3K27me3.fc.signal.bigwig_binnedchromatinState.csv", select = c(2,3,4,7,8))

summ <- tmp %>% full_join(tmp1, by=c("seqnames", "start","end","stat")) %>%
  dplyr::rename(Rb_sig = average.x, WT_sig = average.y) %>%
  filter(seqnames != "chrM") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrX") %>% 
  drop_na() %>%
  mutate(Rb_sig = Rb_sig/median(Rb_sig)) %>%
  mutate(WT_sig = WT_sig/median(WT_sig)) %>% 
  mutate(FC = log2(Rb_sig/WT_sig))
summ2 <- summ
summ2$stat <- "GenomeMed"
summ <- rbind(summ,summ2)

ylabel = bquote(atop(bold("H3K27me3"),log2(FC)))

x_lab_col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","chartreuse4",
              "chartreuse4","firebrick","chartreuse4","chartreuse4","Black", 
              "chartreuse4", "chartreuse4","chartreuse4","firebrick","firebrick",
              "firebrick")

x_lab = c("Bivalent TSS", "Flanking Bivalent TSS/Enhancer", "Bivalent Enhancer",
          "Strong transcription", "Zinc Finger", "Genic enhancers",
          "Repressed PolyComb",  "Weak transcription", "Transcription 5' and 3'",
          "Genome Median", "Active TSS", "Flanking TSS", 
          "Enhancer", "Heterochromatin", "Quiescent", "Weak Repressed PolyComb")

all_state <- unique(summ$stat)
stat_sig <- as.data.frame(matrix(nrow = length(all_state), ncol = 2))
stat_sig$V1 <- all_state
for (i in 1:length(all_state)) {
  subsett = summ %>% 
    filter(stat == all_state[i])
  stat_sig$V2[i] <- wilcox.test(subsett$FC, summ2$FC)$p.val
}
stat_sig$V2 <- round(stat_sig$V2, 5)

ggplot(data = summ, aes(x = reorder(stat, FC, median), y = FC, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,
                   position = "top") +
  scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 100, 30, 10)) + 
  geom_signif(y_position=c(2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4), 
              xmin=c(1,2,3,4,5,6,7,8,9,12,13,14,15,16), xmax=c(1,2,3,4,5,6,7,8,9,12,13,14,15,16),
              annotation=c("*"), tip_length=0, linetype = "blank")
ggsave("plot/H3K27me3_chrom_with_sig.pdf" , width=5.5, height=4, dpi=600, units="in")
plot_10b("plot/10b_H3K27me3.pdf", w = 5.5, h = 1.5)

tmp <- fread("data/binned_chrom_state/RbKO_CPD100.fc.signal.bw_binnedchromatinState.csv", select = c(2,3,4,7,8))
tmp1 <- fread("data/binned_chrom_state/WT_CPD100.fc.signal.bw_binnedchromatinState.csv", select = c(2,3,4,7,8))

summ <- tmp %>% full_join(tmp1, by=c("seqnames", "start","end","stat")) %>%
  dplyr::rename(Rb_sig = average.x, WT_sig = average.y) %>%
  filter(seqnames != "chrM") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrX") %>% 
  drop_na() %>%
  mutate(Rb_sig = Rb_sig/median(Rb_sig)) %>%
  mutate(WT_sig = WT_sig/median(WT_sig)) %>% 
  mutate(FC = log2(Rb_sig/WT_sig))
summ2 <- summ
summ2$stat <- "GenomeMed"
summ <- rbind(summ,summ2)

ylabel = bquote(atop(bold("DNA Lesion"),log2(FC)))

x_lab_col = c("chartreuse4","chartreuse4","firebrick","chartreuse4","Black",
              "firebrick","chartreuse4","firebrick","chartreuse4","chartreuse4",
              "firebrick","chartreuse4","chartreuse4","chartreuse4","chartreuse4",
              "chartreuse4",)

x_lab = c("Weak Transcription","Transcription","Weak Repressed PolyComb",
          "Genic Enhancer", "Genome Median", "Quiescent",
          "Enhancer", "Repressed PolyComb", "Zinc Finger",
          "Transcription 5' and 3'", "Heterochromatin", "Active TSS",
          "Bivalent Enhancer", "Flanking TSS", "Bivalent TSS", "Flanking Bivalent TSS/Enhancer" )

all_state <- unique(summ$stat)
stat_sig <- as.data.frame(matrix(nrow = length(all_state), ncol = 2))
stat_sig$V1 <- all_state
for (i in 1:length(all_state)) {
  subsett = summ %>% 
    filter(stat == all_state[i])
  stat_sig$V2[i] <- wilcox.test(subsett$FC, summ2$FC)$p.val
}
stat_sig$V2 <- round(stat_sig$V2, 5)

ggplot(data = summ, aes(x = reorder(stat, FC), y = FC, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,
                   position = "top") +
  scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 100, 30, 10)) + 
  geom_signif(y_position=c(2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4), 
              xmin=c(1,2,3,6,7,8,9,10,11,12,13,14,15,16), xmax=c(1,2,3,6,7,8,9,10,11,12,13,14,15,16),
              annotation=c("*"), tip_length=0, linetype = "blank")
ggsave("plot/CPD_chrom_with_sig.pdf" , width=5.5, height=4, dpi=600, units="in")
plot_10b("plot/10b_CPD.pdf", w = 5.5, h = 1)
