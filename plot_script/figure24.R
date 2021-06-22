setwd("~/Desktop/yr5/summ/bioinformatics_util/data/telomer_data/")
source("../../draft.R")
pkg <- c("ggplot2", "data.table")
install_pkg(pkg,"common")
dat <- fread("qPCR_quant.csv")

dat$Passage = factor(dat$Passage, levels = c("Mid-1","Mid-2","Mid-3","Late-1","Late-2"))
ggplot(dat, aes(x=Passage, y = Ratio, fill = Passage)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(x = Passage, group = Passage, ymin = LowerBound, ymax = UpperBound),
                width = 0.2,  position = position_dodge(width = 1))  + theme_bw() + 
  geom_hline(yintercept=1, linetype="dashed") + ylab("Telomere content (Fold change)") + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("../../plot/fig24.pdf", width = 4, height = 3)
