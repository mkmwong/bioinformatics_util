setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")

tab <- fread("cooccurrence_v3.tsv", select = c(1,2,7,9)) %>% 
  mutate(`Log2 Odds Ratio`=replace(`Log2 Odds Ratio`, `Log2 Odds Ratio`==">3", 3),
         `Log2 Odds Ratio`=replace(`Log2 Odds Ratio`, `Log2 Odds Ratio`=="<-3", -3),
         `q-Value` = replace(`q-Value`, `q-Value`=="<0.001", 0.0009)) %>%
  mutate(`Log2 Odds Ratio` = as.numeric(`Log2 Odds Ratio`)) %>%
  mutate(`q-Value` = as.numeric(`q-Value`))
tab <- tab %>% mutate(label = case_when(`q-Value` > 0.05 ~ " ",
                                        `q-Value` <= 0.05 & `q-Value` > 0.01 ~ "*",
                                        `q-Value` <= 0.01 & `q-Value` > 0.001 ~ "**",
                                        `q-Value` <= 0.001 ~ "***"))
tab2 <- tab %>% mutate(tmp = A, 
                        A = B,
                        B = tmp) %>%
  dplyr::select (-tmp)
tab3 <- as.data.frame(matrix(nrow = 9, ncol = 5)) 
colnames(tab3) <- colnames(tab)
tab3 <- tab3 %>% mutate(A = unique(c(tab$A, tab$B)),
                B = unique(c(tab$A, tab$B)),
                `Log2 Odds Ratio` = NA,
                `q-Value` = NA,
                label = NA)
summ <- rbind(tab, tab2, tab3) 

ggplot(summ, aes(x = A, y =B, fill = `Log2 Odds Ratio`)) + geom_tile() +
  scale_fill_gradient2( low = "blue", high = "red", mid = "white", midpoint = 0) + 
  geom_text(aes(label=(label)), size = 5) + theme_bw() + 
  theme(axis.title = element_blank(), legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("../plot/figure14_v3.pdf", width = 5, height = 3.2)
