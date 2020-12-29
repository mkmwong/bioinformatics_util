#### Corresponding to fig.9
source("~/Desktop/work_tmp/bioinformatics_util/draft.R")

pkg <- c("dplyr","tidyr")
install_pkg(pkg,"common")

a<-merge_df("~/Desktop/work_tmp/bioinformatics_util/100kb/", c("seqnames","start","end"), 
                      4, TRUE, "~/Desktop/summary", c(2,3,4,7))
#### Alternative could be reaad in after running merge_df.
#a <- readRDS(file = "~/Desktop/summary.rds")
a <- a %>% 
  rename(RbKO_CPD = Rb_CPD, H2AZ = H2A, 
         Rep_Time = wgEncodeFsuRepliChipImr90WaveSignalRep1 ) %>%
  mutate(WT_CPD = WT_CPD/median(na.omit(WT_CPD)), 
         RbKO_CPD = RbKO_CPD/median(na.omit(RbKO_CPD)),
         FC_CPD = RbKO_CPD/WT_CPD) %>% 
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  select(-c(seqnames, start, end)) %>% 
  mutate_at(vars(-Rep_Time), log2) %>% 
  drop_na()
a <- data.matrix(a)

col <- c("skyblue1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "skyblue1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "lightsalmon1","skyblue1","skyblue1","skyblue1","skyblue1","skyblue1",
        "skyblue1","skyblue1","skyblue1","skyblue1","lightsalmon1","skyblue1",
        "skyblue1","skyblue1","skyblue1", "skyblue1", "lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1",
        "lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1","lightsalmon1")
all_out <- pca_plot (a, "~/Desktop/work_tmp/bioinformatics_util/PCA", col, 0.2)