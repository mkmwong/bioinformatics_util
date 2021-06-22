setwd("~/Desktop/yr5/summ/bioinformatics_util/")

dat <- fread("data/telomer_data/teltale_summary.csv")
ggplot(data = dat, aes(x = V1, y = `Portion of telomeric read` , fill = V1)) + 
  geom_bar(stat="identity") + theme_bw() + ylab("Normalized telomeric reads") + xlab("") + 
  geom_errorbar(aes(x = V1, ymin = LowerBound, ymax = UpperBound, width = 0.2)) + 
  theme(legend.position = "None", axis.title = element_text(face = "bold")) 
ggsave("plot/teltale.pdf", w=2, h=3.5)

dat <- fread("data/telomer_data/telomere_hunter_summary.csv")
dat1 <- dat %>% filter(V1 == "Total")
ggplot(data = dat1, aes(x = Group, y = `Normalized count` , fill = Group)) + 
  geom_bar(stat="identity") + theme_bw() + ylab("Normalized telomere content") + xlab("") + 
  geom_errorbar(aes(x = Group, ymin = LowerBound, ymax = UpperBound, width = 0.2)) + 
  theme(legend.position = "None", axis.title = element_text(face = "bold")) 
ggsave("plot/telomere_hunter.pdf", w=2, h=3.5)

dat2 <- dat %>% filter(V1 != "Total")
dat2$V1 = factor(dat2$V1, levels = c("TTAGGG","TGAGGG","TCAGGG","GTTGGG","Others"))
ggplot(data = dat2, aes(x=V1, y = `Normalized count`, fill = Group, group = Group)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.9) + theme_bw() + 
  geom_errorbar(aes(x = V1, group = Group, ymin = LowerBound, ymax = UpperBound),
                width = 0.2,  position = position_dodge(width = 1)) + xlab("") + 
  theme(legend.position = "None", axis.title = element_text(face = "bold")) + 
  ylab("Normalized telomere content")
ggsave("plot/telomere_hunter_2.pdf", w=4, h=3.5)

ggplot(b, aes(x = variable, y = value, fill = Pattern)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("firebrick3","chartreuse3","dodgerblue1","goldenrod1","black")) +
  ylab("Telomere Content (Input)") + scale_x_discrete(labels=c("RbKO", "WT")) + xlab("") + 
  theme_bw() + theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold"))
