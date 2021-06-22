setwd("~/Desktop/final_figures/")
library(ggplot2)
rbko = read.csv("data/binned_bigwig/chromatinState/RbKO_CPD100.fc.signal.bw_binnedchromatinState.csv")
wt = read.csv("data/binned_bigwig/chromatinState/WT_CPD100.fc.signal.bw_binnedchromatinState.csv")
rbko = read.csv("data/binned_bigwig/chromatinState/RbKO_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv")
wt = read.csv("data/binned_bigwig/chromatinState/WT_H3K9me3.fc.signal.bigwig_binnedchromatinState.csv")
sum = merge(rbko, wt[,c(2,3,4,8)], by=c("seqnames","start","end"))
sum = sum[which(sum$seqnames!="chrX"),]
sum = sum[which(sum$seqnames!="chrY"),]
colnames(sum)[c(8,9)] = c("RbKO","WT")
sum$RbKO = sum$RbKO/median(na.omit(sum$RbKO))
sum$WT = sum$WT/median(na.omit(sum$WT))
sum$ratio = sum$RbKO/sum$WT
sum[,c(8:10)] = log2(sum[,c(8:10)])
sum2 = sum
sum2$stat = "GenomeMed"
sum = rbind(sum,sum2)
sum = na.omit(sum)
summary(rbko$width)
hist(rbko$width, breaks = 1000000, xlim=c(0,10000))
#fill is common
x_fill = c("#00FF56","#114662","#6FFFB9","#99C238","#F64844",
           "#FF00AB","#F07202","#00983A","#17D56D","#004A1E",
           "#FB81F1","#40A0CB","#CFE6F6","#999999","#F80000","#FEFEFE")

#For CPD
#order: 5,4,14,6,MED,
#       15,7,13,8,3,
#       9,1,12,2,10,11
#labels = c("Active TSS","Flanking TSS","Transcription 5' and 3'",
#"Transcription","Weak Transcription","Genic Enhancer",
#"Enhancer","Zinc Finger","Heterochromatin","Bivalent TSS",
#"Flanking Bivalent TSS/Enhancer","Bivalent Enhancer","Repressed PolyComb",
#"Weak Repressed PolyComb","Quiescent","Genome Median")
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


ggsave( filename = "figure11_13_20/CPD_ratio.pdf",
ggplot(data = sum, aes(x = reorder(stat, ratio), y = ratio, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,
                   position = "top") +
  scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 110, 30, 0)),
width=5.5, height=4, dpi=600, units="in", device="pdf")


#For H3K9me3
#order: 10,12,11,8,13,
#       4,6,14,3,Med
#       9,5,1,2,7,15
#labels = c("Active TSS","Flanking TSS","Transcription 5' and 3'",
#"Transcription","Weak Transcription","Genic Enhancer",
#"Enhancer","Zinc Finger","Heterochromatin","Bivalent TSS",
#"Flanking Bivalent TSS/Enhancer","Bivalent Enhancer","Repressed PolyComb",
#"Weak Repressed PolyComb","Quiescent","Genome Median")
ylabel = bquote(atop(bold("H3K9me3"),log2(FC)))

x_lab_col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","firebrick",
              "chartreuse4","chartreuse4","firebrick","Black","chartreuse4",
              "firebrick","chartreuse4","chartreuse4","chartreuse4","chartreuse4",
              "firebrick")

x_lab = c("Bivalent TSS","Bivalent Enhancer","Flanking Bivalent TSS/Enhancer",
          "Zinc Finger","Repressed PolyComb","Transcription","Genic Enhancer",
          "Weak Repressed PolyComb","Transcription 5' and 3'","Genome Median",
          "Heterochromatin","Weak Transcription","Active TSS","Flanking TSS",
          "Enhancer","Quiescent")
ggsave( filename = "figure11_13_20/H3K9me3_ratio.pdf",
ggplot(data = sum, aes(x = reorder(stat, ratio), y = ratio, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,
                   position = "top") +
  scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 100, 30, 10)),
width=5.5, height=4, dpi=600, units="in", device="pdf")


## For H3K27me3
#order: 10,11,4,6,8,
#       12,3,5,1,2,
#       MED,7,13,15,9,14
#labels = c("Active TSS","Flanking TSS","Transcription 5' and 3'",
#"Transcription","Weak Transcription","Genic Enhancer",
#"Enhancer","Zinc Finger","Heterochromatin","Bivalent TSS",
#"Flanking Bivalent TSS/Enhancer","Bivalent Enhancer","Repressed PolyComb",
#"Weak Repressed PolyComb","Quiescent/","Genome Median")
ylabel = bquote(atop(bold("H3K27me3"),log2(FC)))

x_lab_col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","chartreuse4",
              "chartreuse4","chartreuse4","chartreuse4","chartreuse4","chartreuse4",
              "Black","chartreuse4","firebrick","firebrick","firebrick",
              "firebrick")

x_lab = c("Bivalent TSS", "Flanking Bivalent TSS/Enhancer", "Strong transcription",
          "Genic enhancers", "Zinc Finger","Bivalent Enhancer",
          "Transcription 5' and 3'", "Weak transcription","Active TSS",
          "Flanking TSS", "Genome Median", "Enhancer",
          "Repressed PolyComb","Quiescent", "Heterochromatin", "Weak Repressed PolyComb")
ggsave( filename = "figure11_13_20/H3K27me3_ratio.pdf",
ggplot(data = sum, aes(x = reorder(stat, ratio), y = ratio, fill = stat)) + geom_boxplot(outlier.color = NA) +
  ylim(-2.5,2.5) + theme_bw() + ylab(ylabel) + xlab("") +
  scale_x_discrete(labels = x_lab ,
                   position = "top") +
  scale_fill_manual(values= x_fill) +
  theme(axis.text.x= element_text(angle = 30, vjust=0, hjust=0, color = x_lab_col ), axis.title.x = element_text(face = "bold"),
        legend.position = "none",  plot.margin = margin(10, 100, 30, 10)),
width=5.5, height=4, dpi=600, units="in", device="pdf")
