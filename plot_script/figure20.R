setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
pkg <- c("ggplot2","dplyr","tidyr","data.table","ggsignif","stringr")
install_pkg(pkg,"common")
#### getting just the file names
dat <- fread("/Users/kamanwong/Desktop/tmp/MorrisonLab/Computational/tool/annovar/annnotated_stjude_tert.variant_function")
dat1 <- dat %>% count(V1, V8) %>% 
  separate(V8, into =c("name"), sep="\\.") %>% 
  filter(!grepl("M", name)) %>%
  filter(!grepl("X", name)) %>% 
  mutate(group=ifelse(grepl('D',name), 'Tumor', 'Germline')) 

dat2 <- dat1 %>% select(name, group)
dat2 <- unique(dat2)
dat2 <- dat2 %>%
  separate(name, into =c("first", "second"))
name_counts <- as.data.frame(table(dat2$first))
name_counts <- name_counts$Var1[which(name_counts$Freq==2)]

dat3 <- as.data.frame(matrix(ncol = 7, nrow = 0))
for (i in 1:length(name_counts)) {
  print(name_counts[i])
  normal <- dat %>% filter(grepl(paste0(name_counts[i],"_G"), V8)) 
  tumor <- dat %>% filter(grepl(paste0(name_counts[i],"_D"), V8))
  normal_name <- unique(normal$V8)
  tumor_name <- unique(tumor$V8)
  normal <- normal %>% select(-V8)
  tumor <- tumor %>% select(-V8)
  print(head(normal))
  print(head(tumor))
  normal_unique <- setdiff(normal, tumor) %>% 
    mutate(name = normal_name)
  tumor_unique <- setdiff( tumor, normal) %>%
    mutate(name = tumor_name)
  dat3 <- rbind(dat3, normal_unique)
  dat3 <- rbind(dat3, tumor_unique)
}
dat3 <- dat3 %>% 
  mutate(group=ifelse(grepl('D',name), 'Tumor', 'Germline')) 
dat4 <- dat3 %>% select(V1, name, group)
dat4 <- unique(dat4) %>% mutate(V1 = str_to_title(V1))
dat4$V1 <- factor(dat4$V1 , levels = c("Upstream","Downstream","Intronic"))
plotting <- as.data.frame(table(dat4$V1, dat4$group)) 
plotting$fisherp <- 1
for (i in 1:(length(plotting$Var1)/2)) {
  tmp <- rbind(c(plotting$Freq[i], plotting$Freq[i+(length(plotting$Var1)/2)]), 
               c(34-plotting$Freq[i], 34-plotting$Freq[i+(length(plotting$Var1)/2)]))
  print(tmp)
  out <- fisher.test(tmp)$p.value
  print(out)
  plotting$fisherp[i] <- out
  plotting$fisherp[i+(length(plotting$Var1)/2)] <- out
}

plotting$Freq[which(plotting$Var2=="Tumor")] <- plotting$Freq[which(plotting$Var2=="Tumor")]/34
plotting$Freq[which(plotting$Var2=="Germline")] <- plotting$Freq[which(plotting$Var2=="Germline")]/34
#plotting <- plotting %>% mutate(Var1 = str_to_title(Var1))
#p#lotting$Var1 <- factor(plotting$Var1 , levels = c("Upstream","Downstream","Intronic"))

ggplot(plotting, aes(x = Var1, y = Freq, group = Var2, fill = Var2)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  xlab("Type") + ylab("Percentage of patients") + ylim(0,1.1) +
  scale_fill_discrete(name = "Sample", labels = c("Tumor", "Matched-normal")) +
  geom_signif(y_position=c(0.2,0.45,1.08), xmin=c(0.75,1.75,2.75), xmax=c(1.25,2.25,3.25),
              annotation=c("NS","NS","NS"), tip_length=0.02) +
  labs(title = "TERT") + theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold")) 
ggsave("plot/figure20.pdf", width = 4.5, height = 4)



