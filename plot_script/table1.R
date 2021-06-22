setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
source("figure12.R")
install_pkg(c("dplyr","tidyr","data.table"),"common")
install_pkg(c("biomaRt","GenomicRanges"),"bioconductor")

## getting protein coding gene information from biomart
hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
pcg <- getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"), 
             filters = "biotype", values="protein_coding", mart=hg19)
pcg <- pcg %>% 
  filter(!grepl('H', chromosome_name)) %>%
  filter(!grepl('G', chromosome_name)) %>%
  filter(!grepl('MT', chromosome_name)) %>%
  filter(!grepl('X', chromosome_name)) %>%
  filter(!grepl('Y', chromosome_name)) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name))
pcgGR <- GRanges(seqnames = pcg$chromosome_name, range = IRanges(start = pcg$start_position, end = pcg$end_position), 
                 hgnc_symbol = pcg$hgnc_symbol, ensembl = pcg$ensembl_gene_id)

## cal percent C>T mutation for each gene
donor <- fread("data/donorid.tsv", header = FALSE)
tmp <- fread("data/donorid_mutationid_geneaffected.tsv") 
tmp <- unique(tmp) %>% 
  filter(mutated_from_allele == 'C' | mutated_from_allele == 'T' | mutated_from_allele == 'A' | mutated_from_allele == 'G') %>% 
  filter(mutated_to_allele == 'C' | mutated_to_allele == 'T' | mutated_to_allele == 'A' | mutated_to_allele == 'G') %>% 
  mutate(loc = paste(chromosome, chromosome_start, chromosome_end))
tmp <- tmp %>% filter(icgc_donor_id %in% donor$V1)
saveRDS(tmp, "data/filtered_unique_mutations.rds")

tmp <- tmp %>% 
  dplyr::select(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele) %>%
  mutate(chromosome = paste0("chr", chromosome))
#tmp <- unique(tmp)
tmp <- tmp %>% count(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele)
tmp$CtoT = FALSE
tmp$CtoT[which(tmp$mutated_from_allele == "C" & tmp$mutated_to_allele == "T")] = TRUE
tmp$CtoT[which(tmp$mutated_from_allele == "G" & tmp$mutated_to_allele == "A")] = TRUE

tmp_false = tmp %>% filter(CtoT == FALSE)
tmp_fGR = GRanges(seqnames = tmp_false$chromosome, range = IRanges(start = tmp_false$chromosome_start, end = tmp_false$chromosome_end), n = tmp_false$n)
#pcgGR$not_CtoT_count = countOverlaps(pcgGR, tmp_fGR)
fol <- as.data.frame(findOverlaps(pcgGR, tmp_fGR))
fol$n <- tmp_fGR$n[fol$subjectHits]
fol <- aggregate(n~ queryHits, fol, sum)
pcgGR$not_CtoT_count = 0
pcgGR$not_CtoT_count[fol$queryHits] = fol$n
tmp_true = tmp %>% filter(CtoT == TRUE)
tmp_tGR = GRanges(seqnames = tmp_true$chromosome, range = IRanges(start = tmp_true$chromosome_start, end = tmp_true$chromosome_end), n = tmp_true$n)
#pcgGR$CtoT_count = countOverlaps(pcgGR, tmp_tGR)
tol <- as.data.frame(findOverlaps(pcgGR, tmp_tGR))
tol$n <- tmp_tGR$n[tol$subjectHits]
tol <- aggregate(n~ queryHits, tol, sum)
pcgGR$CtoT_count = 0
pcgGR$CtoT_count[tol$queryHits] = tol$n
pcgGR$percent_CtoT = pcgGR$CtoT_count/(pcgGR$not_CtoT_count + pcgGR$CtoT_count)
saveRDS(pcgGR, "data/pcg_ctot.rds")

### export tmp(unique mutations for each donor) for annovar 
### col = chrom, start, end, from, to , donor
tmp = readRDS("data/filtered_unique_mutations.rds")
tmp = tmp %>% dplyr::select(chromosome, chromosome_start, chromosome_end, 
                            mutated_from_allele, mutated_to_allele, icgc_donor_id)
write.table(tmp, "data/MELA-AU_Input_for_annovar.tsv", quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE  )

### count percent mutated melanoma for each gene
exon <- fread("/Users/kamanwong/Desktop/tmp/MorrisonLab/Computational/tool/annovar/annnotated_MelaAU.exonic_variant_function",  header = FALSE)
exon <- exon %>% filter(V2 != "synonymous SNV")
exon_sep <- separate_rows(exon,V3) %>% 
  filter(grepl("ENSG", V3))
exon_sep = unique(exon_sep)
exon_tmp <- unique(exon_sep %>% dplyr::select(V3, V9))
mutation_rate_per_gene <- as.data.frame(table(exon_tmp$V3))

#tmp = readRDS("data/filtered_unique_mutations.rds")
##tmpGR = GRanges(seqnames = paste0("chr",tmp$chromosome), range = IRanges(start = tmp$chromosome_start, end = tmp$chromosome_end), 
#                donor = tmp$icgc_donor_id)
#test = findOverlaps(tmpGR, pcgGR)
#tmpGR$gene = NA
#tmpGR[test@from]$gene = pcgGR[test@to]$hgnc_symbol
#tmp1 <- as.data.frame(tmpGR) %>%
#  drop_na() %>% 
#  dplyr::select(donor, gene)
#tmp1 <- unique(tmp1)
#tmp2 <- table(tmp1$gene)
#View(tmp2)

cdg1 <- fread("data/cancer_driver_gene/Bailey2018.csv", select = 1, col.names = "gene")
cdg2 <- fread("data/cancer_driver_gene/COSMIC_CGS.csv", select = 1, col.names = "gene")
cdg3 <- fread("data/cancer_driver_gene/IntOGen-DriverGenesv2020-02-01.tsv", select = 1, col.names = "gene")
cdg4 <- fread("data/cancer_driver_gene/Dietlein_2020.csv")
summ <- unique(unlist(cdg1 %>%  union(cdg2) %>% union(cdg3) %>% union(cdg4)))

pcg$dist_to_LAD <- mapply(dist_to_LAD,pcg$chromosome_name, pcg$start_position, pcg$end_position, SIMPLIFY = TRUE)
pcg$dist_to_centro_telo <- mapply(dist_to_centrotelo,pcg$chromosome_name, pcg$start_position, pcg$end_position, SIMPLIFY = TRUE)
pcg$is_cancer_driver <- FALSE
pcg$is_cancer_driver[which(pcg$hgnc_symbol %in% summ)] <- TRUE

##### merging the data together into one genomic ranges
CToT <- readRDS( "data/pcg_ctot.rds")
CToT <- as.data.frame(CToT) 
colnames(CToT)[7] <- c("ensembl_gene_id")
colnames(mutation_rate_per_gene) <- c("ensembl_gene_id", "melanoma_mutation_rate")
tmp2 <- pcg %>% left_join(mutation_rate_per_gene, by ="ensembl_gene_id") %>%
  mutate(melanoma_mutation_rate = melanoma_mutation_rate/140) %>%
  full_join(CToT, by=c("ensembl_gene_id")) %>%
  dplyr::select(seqnames, start, end, ensembl_gene_id, hgnc_symbol.x, is_cancer_driver, 
                dist_to_LAD, dist_to_centro_telo, percent_CtoT, melanoma_mutation_rate)
write.csv(tmp2, "plot/table1_gene_overlap_with_bins/table1_v2.csv")
#pcgGR <- GRanges(seqnames = pcg$chromosome_name, ranges = IRanges(start = pcg$start_position, end = pcg$end_position), 
#                 ensembl = pcg$ensembl_gene_id, symbol = pcg$hgnc_symbol, dist_to_LAD = pcg$dist_to_LAD, 
#                 dist_to_centro_telo = pcg$dist_to_centro_telo)
pcgGR <- GRanges(seqnames = tmp2$seqnames, ranges = IRanges(tmp2$start, tmp2$end),
                hgnc_symbol = tmp2$hgnc_symbol.x, is_cancer_driver = tmp2$is_cancer_driver,
                ensembl_gene_id = tmp2$ensembl_gene_id, dist_to_LAD = tmp2$dist_to_LAD,
                dist_to_centro_telo = tmp2$dist_to_centro_telo, percentage_CtoT = tmp2$percent_CtoT,
                melanoma_mutation_rate = tmp2$melanoma_mutation_rate)
tmp <- fread("data/100kb/ML-RbKO_CPD.fc.signal.bigwig_binned100kb.csv", select=c(2,3,4,7))
tmp1 <- fread("data/100kb/ML-WT_CPD.fc.signal.bigwig_binned100kb.csv", select=c(2,3,4,7))
summ <- tmp %>% 
  full_join(tmp1, by=c("seqnames","start","end")) %>%
  filter(seqnames != "chrX") %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrM") %>%
  drop_na() %>%
  mutate(average.x = average.x/median(average.x)) %>%
  mutate(average.y = average.y/median(average.y)) %>%
  mutate(FC = log2(average.x/average.y))

top <- summ %>%
  filter(FC >= quantile(FC, 0.9))
bottom <- summ %>%
  filter(FC <= quantile(FC, 0.1))
topGR <- GRanges(seqnames = top$seqnames, ranges = IRanges(top$start, top$end), sig = top$FC)
bottomGR <- GRanges(seqnames = bottom$seqnames, ranges = IRanges(bottom$start, bottom$end), sig = bottom$FC)

hits <- findOverlaps(topGR, pcgGR)
overlaps <- pintersect(topGR[queryHits(hits)], pcgGR[subjectHits(hits)])
pcgGR$topOverlap <- 0
pcgGR$FC <- 0 
pcgGR[subjectHits(hits)]$topOverlap <- overlaps@ranges@width / pcgGR[subjectHits(hits)]@ranges@width
pcgGR[subjectHits(hits)]$FC <- topGR[queryHits(hits)]$sig
top_gene <- as.data.frame(pcgGR) %>%
  dplyr::filter(topOverlap >= 0.1)
write.csv(top_gene, "plot/table1_gene_overlap_with_bins/top10per_overlap10per_v2.csv")
quantile(summ$FC, 0.95)
top_gene <- top_gene %>% filter(FC >= quantile(summ$FC, 0.95))
write.csv(top_gene, "plot/table1_gene_overlap_with_bins/top5per_overlap10per_v2.csv")

hits <- findOverlaps(bottomGR, pcgGR)
overlaps <- pintersect(bottomGR[queryHits(hits)], pcgGR[subjectHits(hits)])
pcgGR$topOverlap <- NULL
pcgGR$bottomOverlap <- 0
pcgGR$FC <- 0 
pcgGR[subjectHits(hits)]$bottomOverlap <- overlaps@ranges@width / pcgGR[subjectHits(hits)]@ranges@width
pcgGR[subjectHits(hits)]$FC <- bottomGR[queryHits(hits)]$sig
bottom_gene <- as.data.frame(pcgGR) %>%
  dplyr::filter(bottomOverlap >= 0.1)
write.csv(bottom_gene, "plot/table1_gene_overlap_with_bins/bottom10per_overlap10per_v2.csv")
bottom_gene <- bottom_gene %>% filter(FC <= quantile(summ$FC, 0.05))
write.csv(bottom_gene, "plot/table1_gene_overlap_with_bins/bottom5per_overlap10per_v2.csv")





