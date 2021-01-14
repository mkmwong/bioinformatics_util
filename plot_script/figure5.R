#### Corresponding to fig.9
source("~/Desktop/yr5/summ/bioinformatics_util/draft.R")

pkg <- c("dplyr","tidyr")
install_pkg(pkg,"common")

a<-merge_df("~/Desktop/yr5/summ/bioinformatics_util/data/100kb/", c("seqnames","start","end"), 
                      4, TRUE, "~/Desktop/yr5/summ/bioinformatics_util/data/100kb_summary", c(2,3,4,7))
#### Alternative could be reaad in after running merge_df.
#a <- readRDS(file = "~/Desktop/yr5/summ/bioinformatics_util/data/100kb_summary.rds")
a <- a %>% 
  dplyr::rename( H2AZ = H2A, 
         Rep_Time = wgEncodeFsuRepliChipImr90WaveSignalRep1,
         RbKO_H3K27me3 = RbKO_H3K27me3_new,
         WT_H3K27me3 = WT_H3K27me3_new,
         RbKO_H3K9me3 = RbKO_H3K9me3_new,
         WT_H3K9me3 = WT_H3K9me3_new) %>%
  mutate(WT_CPD = WT_CPD/median(na.omit(WT_CPD)), 
         RbKO_CPD = RbKO_CPD/median(na.omit(RbKO_CPD)),
         FC_CPD = RbKO_CPD/WT_CPD) %>% 
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  dplyr::select(-c(seqnames, start, end, RbKO_H3K27me3_old,
            WT_H3K27me3_old, RbKO_H3K9me3_old, WT_H3K9me3_old)) %>% 
  mutate_at(vars(-Rep_Time), log2) %>% 
  drop_na()
a <- data.matrix(a)

col <- c("lightsalmon1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "skyblue1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "skyblue1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "skyblue1","skyblue1","skyblue1","skyblue1","skyblue1","lightsalmon1",
        "skyblue1","skyblue1","skyblue1", "skyblue1", "skyblue1","lightsalmon1","lightsalmon1","lightsalmon1",
        "lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1")
all_out <- pca_plot (a, "~/Desktop/yr5/summ/bioinformatics_util/PCA", col, 0.2)
