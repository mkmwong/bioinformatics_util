setwd("~/Desktop/final_figures/")
library(scales)
library(reshape2)
library(plyr)
H3K9me3 = read.csv("data/binned_bigwig/repeats/E017-H3K9me3.fc.signal.bigwig_binnedrepeat.csv")
H3K9me3 = na.omit(H3K9me3)
#H3K9me3$average = H3K9me3$average/median(H3K9me3$average)
summe3 = ddply(H3K9me3, .(V13), summarize,  sig = mean(average))

H3K9ac = read.csv("data/binned_bigwig/repeats/E017-H3K9ac.fc.signal.bigwig_binnedrepeat.csv")
H3K9ac = na.omit(H3K9ac)
#H3K9ac$average = H3K9ac$average/median(H3K9ac$average)
sumac = ddply(H3K9ac, .(V13), summarize,  sig = mean(average))

summ = merge(summe3, sumac, by=c("V13"))
colnames(summ) = c("Category","H3K9me3","H3K9ac")
summ[,c(2:3)] = log2(summ[,c(2:3)])
summ = summ[-c(7,14,16,21,23,25,29,34,36,45,52,55,56),]
summm = melt(summ)

ggsave(filename = "figure4-in_progress/fig4b.pdf",
ggplot(summm, aes(reorder(Category, rat$value), variable)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient2( low = "blue", high = "red", mid = "white", midpoint = 0) + 
  xlab("Types of repeats") + ylab("Histone Modification") +labs (fill="Signal") + 
  theme_bw() + theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1)),
width=8, height=2, dpi=600, units="in", device="pdf")
