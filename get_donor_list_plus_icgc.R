setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
pkg <- c("data.table","dplyr","tidyr")
install_pkg(pkg)
##extracting importing data from icgc mela-au
'%ni%' <- Negate('%in%')

exon <- fread("/Users/kamanwong/Desktop/tmp/MorrisonLab/Computational/tool/annovar/annnotated_MelaAU.exonic_variant_function",  header = FALSE)
exon <- exon %>% filter(V2 != "synonymous SNV")
exon_sep <- separate_rows(exon,V3) %>% 
  filter(grepl("ENSG", V3))
exon_sep = unique(exon_sep)
exon_tmp <- unique(exon_sep %>% dplyr::select(V3, V9))

####
#all <- fread("/Users/kamanwong/Desktop/tmp/MorrisonLab/Computational/tool/annovar/annnotated_MelaAU.variant_function", header = FALSE)
#all <- all %>% filter(V2 != "synonymous SNV")
#all_sep <- separate_rows(all,V2) %>% 
#  filter(grepl("ENSG", V2))
#saveRDS(all_sep, "data/annotated_mutation.rds")
####
# ENSG00000139687 = rb1; ENSG00000147889 = cdkn2a; ENSG00000124762 = cdkn1a;
# ENSG00000110092 = ccnd1; ENSG00000105173 = ccne1; ENSG00000147883 = cdkn2b;
# ENSG00000111276 = cdkn1b
donor <- fread("data/donorid.tsv", header = FALSE)
rbpath <- exon_tmp %>% 
  filter(V3 == "ENSG00000105173" | V3 == "ENSG00000110092" | V3 == "ENSG00000124762" |
           V3 == "ENSG00000147889" | V3 == "ENSG00000139687" | V3 == "ENSG00000147883" | 
           V3 == "ENSG00000135446" | V3 == "ENSG00000111276")
mut = fread("~/Desktop/tmp/MorrisonLab/Computational/05092019-Rb_related_genes/donor_w_rb_mut.csv", header = FALSE)
length(unique(rbpath$V9))
(unique(rbpath$V9) %in% mut$V1)
write.csv(unique(rbpath$V9),"~/Desktop/rblist.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
donor$V1[which(donor$V1 %ni% rbpath$V9)]
write.csv(donor$V1[which(donor$V1 %ni% rbpath$V9)],"~/Desktop/notrblist.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
#tert <- exon_sep %>% filter(V3 == "ENSG00000164362")

all <- fread("/Users/kamanwong/Desktop/tmp/MorrisonLab/Computational/tool/annovar/annnotated_MelaAU.variant_function", header = FALSE)
tert <- all %>% filter(V3 == "5") 
tert <- separate_rows(tert, V2, sep=",") #%>% 
  #filter(grepl("ENSG", V2))
tert <- tert %>% filter(grepl("ENSG00000164362",V2))
tert_rb <- tert %>% filter(V8 %in% unique(rbpath$V9))
tert_notrb <- tert %>% filter(V8 %ni% unique(rbpath$V9))

not_rb <- table(tert_notrb$V1, tert_notrb$V8)
not_rb[not_rb>0] = 1
rowSums(not_rb)/103

rb <- table(tert_rb$V1, tert_rb$V8)
rb[rb>0] = 1
rowSums(rb)/37

tert <- all %>% filter(V3 == "5") 
tert <- separate_rows(tert, V2) %>% 
  filter(grepl("ENSG", V2))
tert <- tert %>% filter(V2 == "ENSG00000164362")
tert_rb <- tert %>% filter(V8 %in% unique(rbpath$V9))
tert_notrb <- tert %>% filter(V8 %ni% unique(rbpath$V9))

not_rb <- table(tert_notrb$V1, tert_notrb$V8)
not_rb[not_rb>0] = 1
rowSums(not_rb)/107

rb <- table(tert_rb$V1, tert_rb$V8)
rb[rb>0] = 1
rowSums(rb)/33

C250T <- tert %>% filter(V4 == 1295250 & V6 == "G" & V7 == "A")
length(which((C250T$V8 %in% rbpath$V9 )== FALSE))/103
length(which((C250T$V8 %in% rbpath$V9 )== TRUE))/37
C228T <- tert %>% filter(V4 == 1295228 & V6 == "G" & V7 == "A")
length(which((C228T$V8 %in% rbpath$V9 )== FALSE))/103
length(which((C228T$V8 %in% rbpath$V9 )== TRUE))/37
###KMT2C
tert <- all %>% filter(V3 == "7") 
tert <- separate_rows(tert, V2) %>% 
  filter(grepl("ENSG", V2))
tert <- tert %>% filter(V2 == "ENSG00000055609")
tert_rb <- tert %>% filter(V8 %in% unique(rbpath$V9))
tert_notrb <- tert %>% filter(V8 %ni% unique(rbpath$V9))

not_rb <- table(tert_notrb$V1, tert_notrb$V8)
not_rb[not_rb>0] = 1
rowSums(not_rb)/107

rb <- table(tert_rb$V1, tert_rb$V8)
rb[rb>0] = 1
rowSums(rb)/33

###MAPK11
mapk <- all %>% filter(V3 == "22") 
mapk <- separate_rows(mapk, V2) %>% 
  filter(grepl("ENSG", V2))
mapk <- mapk %>% filter(V2 == "ENSG00000185386")
mapk_rb <- mapk %>% filter(V8 %in% unique(rbpath$V9))
mapk_notrb <- mapk %>% filter(V8 %ni% unique(rbpath$V9))

not_rb <- table(mapk_notrb$V1, mapk_notrb$V8)
not_rb[not_rb>0] = 1
rowSums(not_rb)/107

rb <- table(mapk_rb$V1, mapk_rb$V8)
rb[rb>0] = 1
rowSums(rb)/33






