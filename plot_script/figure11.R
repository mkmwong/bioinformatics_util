setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")

pkg <- c("ggplot2", "data.table", "ggsignif")
install_pkg(pkg, "common")

tab <- fread("data/qPCR_ashby/forPlot.csv") %>% 
  drop_na()
tab$Sample = factor(tab$Sample, levels = c("RbKO", "WT"))

ggplot(tab, aes(x = Sample, y=Per_inhibi, fill = Sample)) + geom_bar(stat = "identity") +
  geom_errorbar(aes(x = Sample, ymin = LowerBound, ymax = UpperBound, width = 0.2)) + 
  theme_bw() + ylab("TERT Percentage inhibition") + 
  theme(legend.position = "None", axis.title = element_text(face = "bold"))
ggsave("plot/TERTinhibition_v2.pdf", width = 2, height = 3, dpi = 600)
