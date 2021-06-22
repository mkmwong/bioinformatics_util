setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")

pkg <- c("ggplot2","data.table","tidyverse")
install_pkg(pkg,"common")
pkg <- c("DESeq2","IHW","biomaRt","vsn")
install_pkg(pkg,"bioconductor")

read_binned_repeats <- function(name, col) {
  tmp <- fread(name) %>%
    group_by_(col) %>% 
    drop_na() %>%
    summarize(mean = mean(average))
  return(tmp)
}

summ <- fread("repenrich_out/RbKOH3K9me3_1_out/RbKOH3K9me3_1_fraction_counts.txt", stringsAsFactors = FALSE)
tmp <- fread("repenrich_out/RbKOH3K9me3_2_out/RbKOH3K9me3_2_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp1 <- fread("repenrich_out/WTH3K9me3_1_out/WTH3K9me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp2 <- fread("repenrich_out/WTH3K9me3_2_out/WTH3K9me3_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)
tmp3 <- fread("repenrich_out/RbKOInput_1_out/RbKOInput_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp4 <- fread("repenrich_out/RbKOInput_2_out/RbKOInput_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)
tmp5 <- fread("repenrich_out/WTInput_1_out/WTInput_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp6 <- fread("repenrich_out/WTInput_2_out/WTInput_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)
tmp7 <- fread("repenrich_out/RbKOH3K27me3_1_out/RbKOH3K27me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp8 <- fread("repenrich_out/RbKOH3K27me3_2_out/RbKOH3K27me3_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)
tmp9 <- fread("repenrich_out/WTH3K27me3_1_out/WTH3K27me3_1_fraction_counts.txt", select=c(1,4), stringsAsFactors = FALSE)
tmp10 <- fread("repenrich_out/WTH3K27me3_2_out/WTH3K27me3_2_fraction_counts.txt",  select=c(1,4), stringsAsFactors = FALSE)

summ <- summ %>% 
  full_join(tmp,by="V1") %>%
  full_join(tmp1,by="V1") %>%
  full_join(tmp2, by="V1") %>% 
  full_join(tmp3,by="V1") %>%
  full_join(tmp4,by="V1") %>%
  full_join(tmp5, by="V1") %>% 
  full_join(tmp6, by="V1") %>% 
  full_join(tmp7,by="V1") %>%
  full_join(tmp8,by="V1") %>%
  full_join(tmp9, by="V1") %>% 
  full_join(tmp10, by="V1") %>% 
  tibble::column_to_rownames("V1") %>%
  dplyr::rename(RBKOH3K9me3_1 = V4.x,
                RBKOH3K9me3_2  = V4.y,
                WTH3K9me3_1  = V4.x.x,
                WTH3K9me3_2  = V4.y.y,
                RBKOInput_1  = V4.x.x.x,
                RBKOInput_2 = V4.y.y.y,
                WTInput_1 = V4.x.x.x.x,
                WTInput_2 = V4.y.y.y.y,
                RBKOH3K27me3_1 = V4.x.x.x.x.x,
                RBKOH3K27me3_2  = V4.y.y.y.y.y,
                WTH3K27me3_1 = V4.x.x.x.x.x.x,
                WTH3K27me3_2 = V4.y.y.y.y.y.y)
                
summ1 <- summ %>% 
  mutate(RbK9FC = (RBKOH3K9me3_1 + RBKOH3K9me3_2)/(RBKOInput_1 + RBKOInput_2)) %>%
  mutate(WTK9FC = (WTH3K9me3_1 + WTH3K9me3_2)/(WTInput_1 + WTInput_2)) %>% 
  mutate(RbK27FC = (RBKOH3K27me3_1 + RBKOH3K27me3_2)/(RBKOInput_1 + RBKOInput_2)) %>%
  mutate(WTK27FC = (WTH3K27me3_1 + WTH3K27me3_2)/(WTInput_1 + WTInput_2)) %>% 
  mutate(K9FC = log2(RbK9FC/WTK9FC)) %>%
  mutate(K27FC = log2(RbK27FC/WTK27FC)) %>%
  dplyr::select(c(V2, V3, K9FC, K27FC))

summ1 <- summ1 %>% filter(!grepl("\\?", V3)) 
ggplot(summ1,aes(x= reorder(V3, K9FC, median), y = K9FC, fill = V2 )) + 
  geom_boxplot() + xlab("Repeat family") + ylab("log2FC(RbKO/WT H3K9me3)") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none",
        axis.title = element_text(face = "bold")) + ylim(-1.5, 1.5)
ggsave("../plot/k9_repenrich_familyFC.pdf", width = 15, height = 5, dpi = 600)

ggplot(summ1,aes(x= reorder(V3, K27FC, median), y = K27FC, fill = V2 )) + 
  geom_boxplot() + xlab("Repeat family") + ylab("log2FC(RbKO/WT H3K27me3)") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none",
        axis.title = element_text(face = "bold")) + ylim(-1.5, 1.5)
ggsave("../plot/k27_repenrich_familyFC.pdf", width = 15, height = 5, dpi = 600)


