setwd("~/Desktop/yr5/summ/bioinformatics_util/")
a <- fread("data/MELA-AU_individual_gene_data/Book2.csv")
a$fisherp <- 1
for (i in 1:length(a$V1)) {
  tmp <- rbind(c(a$rb_mutated[i], a$rb_not_mutated[i]), c(37-a$rb_mutated[i], 103-a$rb_not_mutated[i]))
  out <- fisher.test(tmp)$p.value
  a$fisherp[i] <- out
}
a$V1 = factor(a$V1, levels= c("Upstream","Downstream","Exon","Intron","Missense"))

tert <- a %>% filter(gene == "tert") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
tert <- tert %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(tert, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.92,0.75,0.4,0.95,0.15), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("***","*","NS","***","NS"), tip_length=0.02) +
  labs(title = "TERT") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/tert.pdf", dpi = 600, width =6, height = 4)

kmt2c <- a %>% filter(gene == "kmt2c") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
kmt2c <- kmt2c %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(kmt2c, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.6, 0.6, 0.34, 0.98,0.26), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("NS","**","NS","NS","NS"), tip_length=0.02) +
  labs(title = "kmt2c") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/kmt2c.pdf", dpi = 600, width =6, height = 4)

brd9 <- a %>% filter(gene == "brd9") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
brd9 <- brd9 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(brd9, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.6, 0.7, 0.3, 0.81,0.09), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("*","**","NS","***","NS"), tip_length=0.02) +
  labs(title = "brd9") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/brd9.pdf", dpi = 600, width =6, height = 4)

prdm16 <- a %>% filter(gene == "prdm16") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
prdm16 <- prdm16 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(prdm16, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.82,0.75,0.65,1.07,0.29), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("***","*","*","NS","NS"), tip_length=0.02) +
  labs(title = "prdm16") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/prdm16.pdf", dpi = 600, width =6, height = 4)

dnmt3l <- a %>% filter(gene == "dnmt3l") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
dnmt3l <- dnmt3l %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(dnmt3l, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.45, 0.55, 0.05, 0.67, 0.12), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("*","*","NS","***","NS"), tip_length=0.02) +
  labs(title = "dnmt3l") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/dnmt3l.pdf", dpi = 600, width =6, height = 4)

col5a1 <- a %>% filter(gene == "col5a1") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
col5a1 <- col5a1 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(col5a1, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.96,0.95,0.38,1.08,0.5), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("***","*","NS","*","NS"), tip_length=0.02) +
  labs(title = "col5a1") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/col5a1.pdf", dpi = 600, width =6, height = 4)

bcl2 <- a %>% filter(gene == "bcl2") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
bcl2 <- bcl2 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(bcl2, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.3, 0.44, 0.08, 0.9,0.08), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("*", "*", "NS","NS","NS"), tip_length=0.02) +
  labs(title = "bcl2") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/bcl2.pdf", dpi = 600, width =6, height = 4)

fgfr2 <- a %>% filter(gene == "fgfr2") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
fgfr2 <- fgfr2 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(fgfr2, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.74,0.77, 0.21, 1.03, 0.2), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("NS", "NS", "NS","*","NS"), tip_length=0.02) +
  labs(title = "fgfr2") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/fgfr2.pdf", dpi = 600, width =6, height = 4)

tpte <- a %>% filter(gene == "tpte") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
tpte <- tpte %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(tpte, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.9, 0.98, 0.7, 1.07, 0.6), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("*","*","*","NS","**"), tip_length=0.02) +
  labs(title = "tpte") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/tpte.pdf", dpi = 600, width =6, height = 4)

ankrd30a <- a %>% filter(gene == "ankrd30a") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
ankrd30a <- ankrd30a %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(ankrd30a, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.9, 0.98, 0.7, 1.07, 0.6), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("NS", "**", "NS","NS","NS"), tip_length=0.02) +
  labs(title = "ankrd30a") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/ankrd30a.pdf", dpi = 600, width =6, height = 4)

tpte2 <- a %>% filter(gene == "tpte2") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
tpte2 <- tpte2 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(tpte2, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.67, 0.5, 0.2, 1.07, 0.17), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("NS","NS","NS","*","NS"), tip_length=0.02) +
  labs(title = "tpte2") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/tpte2.pdf", dpi = 600, width =6, height = 4)

col4a2 <- a %>% filter(gene == "col4a2") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
col4a2 <- col4a2 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(col4a2, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.71, 0.74, 0.28, 1.07, 0.25), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("**","*","NS","*","NS"), tip_length=0.02) +
  labs(title = "col4a2") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/col4a2.pdf", dpi = 600, width =6, height = 4)

herc2 <- a %>% filter(gene == "herc2") %>% 
  mutate(rb_mutated = rb_mutated/37,
         rb_not_mutated = rb_not_mutated/103)
herc2 <- herc2 %>% gather("group", "val", -c("gene", "fisherp", "V1"))
ggplot(herc2, aes(x = V1, y = val, group = group, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Rb-pathway status", labels = c("Mutated", "Not Mutated")) +
  geom_signif(y_position=c(0.5, 0.69, 0.2, 0.96, 0.26), xmin=c(0.75,1.75,2.75,3.75,4.75), 
              xmax=c(1.25,2.25,3.25,4.25,5.25),annotation=c("NS","*", "NS", "*", "NS"), tip_length=0.02) +
  labs(title = "herc2") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("plot/figure19/herc2.pdf", dpi = 600, width =6, height = 4)












