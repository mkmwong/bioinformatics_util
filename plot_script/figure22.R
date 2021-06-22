setwd("~/Desktop/tmp/MorrisonLab/Computational/06102019-TERT/")
library(ggplot2)
TERT = read.table("TERT_mut.tsv", sep="\t")
TERT = TERT[,c(1,2,7,8,9,10,11,13,14,16,17,29)]
TERT = unique(TERT)

rb_patient = fread("~/Desktop/rblist.csv")
not_rb_patient = fread("~/Desktop/notrblist.csv")

TERT_rb <- TERT %>% filter(V2 %in% rb_patient$x)
TERT_notrb <- TERT %>% filter(V2 %in% not_rb_patient$x)

TERT_rb1 <- TERT_rb %>% filter((V10 == "1295250" & V16 == "G" & V17 == "A") | (V10 == "1295228" & V16 =="G" & V17 == "A") | (V10 == "1295242" & V11 == "1295243"))
length(unique(TERT_rb1$V2))
TERT_notrb1 <- TERT_notrb %>% filter((V10 == "1295250" & V16 == "G" & V17 == "A") | (V10 == "1295228" & V16 =="G" & V17 == "A") | (V10 == "1295242" & V11 == "1295243"))
length(unique(TERT_notrb1$V2))

dat <- fread("~/Desktop/yr5/summ/bioinformatics_util/data/fig5e.csv")
ggplot(dat, aes(x = Type, y = `Porportion of patients`,  fill = `Rb-pathway status`)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1)) + xlab("Mutation")
ggsave("~/Desktop/yr5/summ/bioinformatics_util/fig22.pdf", width = 3.5, height = 3.5)
