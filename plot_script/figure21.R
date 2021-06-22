setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")

pkg <- c("ggplot2", "data.table", "ggsignif")
install_pkg(pkg, "common")

tab <- fread("data/rafael_data/forplot.csv") %>% 
  drop_na()
#tab$Sample = factor(tab$Sample, levels = c("WT", "RbKO"))

ggplot(tab, aes(x = Sample, y=Relative_amplification, fill = Sample)) + geom_bar(stat = "identity") +
  geom_errorbar(aes(x = Sample, ymin = Lower_bound, ymax = Upper_bound, width = 0.2)) + 
  theme_bw() + ylab("MECOM Relative Amplification (%)") + xlab("UV dosage") +
  theme(legend.position = "None", axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plot/MECOMinhibition.pdf", width = 2, height = 4, dpi = 600)
