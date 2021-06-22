setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
install_pkg(c("Biostrings","BSgenome.Hsapiens.UCSC.hg19"), "bioconductor")
install_pkg(c("ggpubr","ggplot2", "tidyr", "data.table","rstatix"), "common")

rbko <- fread("100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select = c(2,3,4)) %>% 
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM")

genome <- BSgenome.Hsapiens.UCSC.hg19
gap <- fread("hg19gap.txt")
centGR <- get_gap_GR(gap, "centromere", 1000000)
rbkoGR <- GRanges(seqnames = rbko$seqnames, ranges = IRanges(start = rbko$start, end = rbko$end))
centro_regionGR <- subsetByOverlaps(rbkoGR,centGR)
cent_count <- get_dipy_counts(centro_regionGR, genome)
cent_count <- cent_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>%
  gather() %>%
  mutate(group = "Pericentric") %>% 
  drop_na()
  
rand_loc <- sample(1:nrow(rbko),length(centro_regionGR),replace=F)
rand_region <- rbko[rand_loc,]
randGR <- GRanges(seqnames = rand_region$seqnames, ranges = IRanges(start = rand_region$start, end = rand_region$end))
rand_count <- get_dipy_counts(randGR, genome)
rand_count <- rand_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>%
  mutate(group = "Random") %>% 
  drop_na()

for (i in 1:999) {
  rand_loc <- sample(1:nrow(rbko),length(centro_regionGR),replace=F)
  rand_region <- rbko[rand_loc,]
  randGR <- GRanges(seqnames = rand_region$seqnames, ranges = IRanges(start = rand_region$start, end = rand_region$end))
  tmp_count <- get_dipy_counts(randGR, genome)
  tmp_count <- tmp_count %>%
    dplyr::select(TpT, TpC, CpT, CpC) %>%
    mutate(group = "Random") %>% 
    drop_na()
  rand_count <- rbind(rand_count, tmp_count)
}

rand_count <- rand_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>% 
  gather() %>%
  mutate(group = "Random")

summ <- rbind(cent_count, rand_count)
stat_test <- summ %>%
  group_by(key) %>%
  wilcox_test(value ~ group, exact = FALSE, paired = FALSE, p.adjust.method = "BH") %>%
  add_significance() %>% 
  mutate(group1 = group1, group2 = group2) %>% 
  add_xy_position(x = "group")

bold_text = element_text(face = "bold")
ggplot(summ, aes(x = group, y=value, fill = group)) + geom_boxplot(outlier.size = 0.2) + 
  facet_wrap(~key, nrow = 1) + theme_bw() + ylab("Porportion") + xlab("") + 
  stat_pvalue_manual(stat_test, inherit.aes = FALSE) + labs(fill = "Group") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1), 
        axis.title = bold_text, legend.title = bold_text, strip.text = bold_text) +
  scale_fill_manual(values = c("darkorchid", "white"))
saveRDS(summ, "figure2_supp_centromere.rds")
ggsave("../plot/figure2_percentric_dipy_wbootstrap.pdf", width = 4.5, height  = 3.5)

### making the same plot for subtelomeric
teloGR <- get_gap_GR(gap, "telomere", 100000)
telo_regionGR <- subsetByOverlaps(rbkoGR,teloGR)
telo_count <- get_dipy_counts(telo_regionGR, genome)
telo_count <- telo_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>%
  gather() %>%
  mutate(group = "Subtelomeric") %>% 
  drop_na()

rand_loc <- sample(1:nrow(rbko),length(telo_regionGR),replace=F)
rand_region <- rbko[rand_loc,]
randGR <- GRanges(seqnames = rand_region$seqnames, ranges = IRanges(start = rand_region$start, end = rand_region$end))
rand_count <- get_dipy_counts(randGR, genome)
rand_count <- rand_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>%
  gather() %>%
  mutate(group = "Random") %>% 
  drop_na()

rand_loc <- sample(1:nrow(rbko),length(telo_regionGR),replace=F)
rand_region <- rbko[rand_loc,]
randGR <- GRanges(seqnames = rand_region$seqnames, ranges = IRanges(start = rand_region$start, end = rand_region$end))
rand_count <- get_dipy_counts(randGR, genome)
rand_count <- rand_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>%
  mutate(group = "Random") %>% 
  drop_na()

for (i in 1:999) {
  rand_loc <- sample(1:nrow(rbko),length(telo_regionGR),replace=F)
  rand_region <- rbko[rand_loc,]
  randGR <- GRanges(seqnames = rand_region$seqnames, ranges = IRanges(start = rand_region$start, end = rand_region$end))
  tmp_count <- get_dipy_counts(randGR, genome)
  tmp_count <- tmp_count %>%
    dplyr::select(TpT, TpC, CpT, CpC) %>%
    mutate(group = "Random") %>% 
    drop_na()
  rand_count <- rbind(rand_count, tmp_count)
}

rand_count <- rand_count %>%
  dplyr::select(TpT, TpC, CpT, CpC) %>% 
  gather() %>%
  mutate(group = "Random")

summ <- rbind(telo_count, rand_count)
summ$group = factor(summ$group, levels = c("Subtelomeric", "Random"))
stat_test <- summ %>%
  group_by(key) %>%
  wilcox_test(value ~ group, exact = FALSE, paired = FALSE, p.adjust.method = "BH") %>%
  add_significance() %>% 
  mutate(group1 = group1, group2 = group2) %>% 
  add_xy_position(x = "group")

saveRDS(summ, "figure2_supp_telomere.rds")
bold_text = element_text(face = "bold")
ggplot(summ, aes(x = group, y=value, fill = group)) + geom_boxplot(outlier.size = 0.2) + 
  facet_wrap(~key, nrow = 1) + theme_bw() + ylab("Porportion") + xlab("") + 
  stat_pvalue_manual(stat_test, inherit.aes = FALSE) + labs(fill = "Group") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1), 
        axis.title = bold_text, legend.title = bold_text, strip.text = bold_text) + 
  scale_fill_manual(values = c("plum2", "white"))
ggsave("../plot/figure2_subtelo_dipy.pdf", width = 4.5, height  = 3.5)
