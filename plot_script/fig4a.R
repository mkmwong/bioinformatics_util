setwd("~/Desktop/final_figures/")
library(scales)
library(cowplot)
library(reshape2)

pt = list()
ratio = read.csv("data/repeat_class_summary.csv")
ratio = ratio[-c(7,14,16,21,23,25,29,34,36,45,52,55,56),]
ratio$Category = factor(ratio$Category, levels = levels(ratio$Category))
for(i in 1:length(ratio$Category)) {
  print(as.character(ratio$Category[i]))
  ratio$top_over_pval[i] = phyper(ratio$Top[i]-1,ratio$All[i], 2.7*10^9-ratio$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE) #* nrow(ratio)
  ratio$bot_over_pval[i] = phyper(ratio$Bottom[i]-1,ratio$All[i], 2.7*10^9-ratio$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE ) #* nrow(ratio)
  ratio$top_dept_pval[i] = phyper(ratio$Top[i],ratio$All[i], 2.7*10^9-ratio$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(ratio)
  ratio$bot_dept_pval[i] = phyper(ratio$Bottom[i],ratio$All[i], 2.7*10^9-ratio$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(ratio)
}
ratio[,c(6:9)] = -(ratio[,c(6:9)]/log(10))

ratio$top_dept_pval[ratio$top_dept_pval < -log10(0.99)] = 0
ratio$top_over_pval[ratio$top_over_pval < -log10(0.99)] = 0
ratio$top_dept_pval[ratio$top_dept_pval > 0] = -ratio$top_dept_pval[ratio$top_dept_pval > 0]

weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log(abs(x)),
                           inverse=function(x) sign(x)*exp(abs(x)))

rat = melt(ratio[,c(2,6,8)])
ylabel = bquote(atop(bold(Depletion %<->% Enrichment), "Bonferroni p-value"))
a = ggplot(data = rat, aes(x = reorder(Category, value), y = value, fill = variable)) + geom_bar(stat = "identity") + 
  scale_y_continuous( trans = weird, breaks = c(-1e8, -1e6, -1e4, -1e2, 0, 1e2, 1e4, 1e6, 1e8),
                     labels = c("1e+08", "1e+06", "1e+04", "1e+02", "0", "1e+02", "1e+04", "1e+06", "1e+08")) +
  ylab(ylabel) + xlab("Types of repeats") + geom_hline(yintercept=1.3, linetype="dashed", color = "black") +
  geom_hline(yintercept=-1.3, linetype="dashed", color = "black") + ggtitle("RbKO/WT Top 10% 1MB bins") +
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

wt = read.csv("data/repeat_class_summary_wt.csv")
wt = wt[-c(7,14,16,21,23,25,29,34,36,45,52,55,56),]
wt$Category = factor(wt$Category, levels = levels(wt$Category))
for(i in 1:length(wt$Category)) {
  print(as.character(wt$Category[i]))
  wt$top_over_pval[i] = phyper(wt$Top[i]-1,wt$All[i], 2.7*10^9-wt$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE) #* nrow(wt)
  wt$bot_over_pval[i] = phyper(wt$Bottom[i]-1,wt$All[i], 2.7*10^9-wt$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE ) #* nrow(wt)
  wt$top_dept_pval[i] = phyper(wt$Top[i],wt$All[i], 2.7*10^9-wt$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(wt)
  wt$bot_dept_pval[i] = phyper(wt$Bottom[i],wt$All[i], 2.7*10^9-wt$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(wt)
}
wt[,c(6:9)] = -(wt[,c(6:9)]/log(10))
wt$top_dept_pval[wt$top_dept_pval < -log10(0.99)] = 0
wt$top_over_pval[wt$top_over_pval < -log10(0.99)] = 0
wt$top_dept_pval[wt$top_dept_pval > 0] = -wt$top_dept_pval[wt$top_dept_pval > 0]

wtt = melt(wt[,c(2,6,8)])

b = ggplot(data = wtt, aes(x = reorder(Category, rat$value), y = value, fill = variable)) + geom_bar(stat = "identity") + 
  scale_y_continuous(trans = weird, breaks = c(-1e8, -1e6, -1e4, -1e2, 0, 1e2, 1e4, 1e6, 1e8),
                     labels = c("1e+08", "1e+06", "1e+04", "1e+02", "0", "1e+02", "1e+04", "1e+06", "1e+08")) +
  ylab(ylabel) + xlab("Types of repeats") + geom_hline(yintercept=1.3, linetype="dashed", color = "black") +
  geom_hline(yintercept=-1.3, linetype="dashed", color = "black") + ggtitle("WT Top 10% 1MB bins") +
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

rbko = read.csv("data/repeat_class_summary_rbko.csv")
rbko = rbko[-c(7,14,16,21,23,25,29,34,36,45,52,55,56),]
rbko$Category = factor(rbko$Category, levels = levels(rbko$Category))
for(i in 1:length(rbko$Category)) {
  print(as.character(rbko$Category[i]))
  rbko$top_over_pval[i] = phyper(rbko$Top[i]-1,rbko$All[i], 2.7*10^9-rbko$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE) #* nrow(rbko)
  rbko$bot_over_pval[i] = phyper(rbko$Bottom[i]-1,rbko$All[i], 2.7*10^9-rbko$All[i], 2.7*10^8, lower.tail = FALSE, log.p = TRUE ) #* nrow(rbko)
  rbko$top_dept_pval[i] = phyper(rbko$Top[i],rbko$All[i], 2.7*10^9-rbko$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(rbko)
  rbko$bot_dept_pval[i] = phyper(rbko$Bottom[i],rbko$All[i], 2.7*10^9-rbko$All[i], 2.7*10^8, lower.tail = TRUE, log.p = TRUE) #* nrow(rbko)
}
rbko[,c(6:9)] = -(rbko[,c(6:9)]/log(10))
rbko$top_dept_pval[rbko$top_dept_pval < -log10(0.99)] = 0
rbko$top_over_pval[rbko$top_over_pval < -log10(0.99)] = 0
rbko$top_dept_pval[rbko$top_dept_pval > 0] = -rbko$top_dept_pval[rbko$top_dept_pval > 0]

rbkot = melt(rbko[,c(2,6,8)])

c = ggplot(data = rbkot, aes(x = reorder(Category, rat$value), y = value, fill = variable)) + geom_bar(stat = "identity") + 
  scale_y_continuous(trans = weird, breaks = c(-1e8, -1e6, -1e4, -1e2, 0, 1e2, 1e4, 1e6, 1e8),
                     labels = c("1e+08", "1e+06", "1e+04", "1e+02", "0", "1e+02", "1e+04", "1e+06", "1e+08")) +
  ylab(ylabel) + xlab("Types of repeats") + geom_hline(yintercept=1.3, linetype="dashed", color = "black") +
  geom_hline(yintercept=-1.3, linetype="dashed", color = "black") + ggtitle("RbKO Top 10% 1MB bins") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10), legend.position = "none")

plot_grid(a,b,c, align = "v", nrow = 3, rel_heights = c(0.3,0.3,0.4))
