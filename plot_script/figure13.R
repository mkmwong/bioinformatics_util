setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("tidyverse","data.table","huxtable","ggsignif")
install_pkg(pkg,"common")
pkg <- c("biomaRt","GenomicRanges","BSgenome.Hsapiens.UCSC.hg19","Biostrings")
install_pkg(pkg,"bioconductor")

'%ni%' <- Negate('%in%')

tmp <- fread("donorid_mutationid_geneaffected.tsv") 
tmp <- unique(tmp) %>% 
  filter((mutated_from_allele == "C" & mutated_to_allele == "T") | (mutated_from_allele == "G" & mutated_to_allele == "A")) %>% 
  dplyr::select(chromosome, chromosome_start, chromosome_end) %>% 
  mutate(loc = paste(chromosome, chromosome_start, chromosome_end))
tmp <- table(tmp$loc)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% separate(Var1, c("chromosome","start","end"))
tmpGR <- GRanges(seqnames = tmp$chromosome, ranges = IRanges(start = as.numeric(tmp$start), end = as.numeric(tmp$end)), freq = tmp$Freq)

hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
pcg <- getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"), 
             filters = "biotype", values="protein_coding", mart=hg19)
pcg <- pcg %>% 
  filter(!grepl('H', chromosome_name)) %>%
  filter(!grepl('G', chromosome_name)) %>%
  filter(!grepl('MT', chromosome_name)) %>%
  filter(!grepl('X', chromosome_name)) 
pcgGR = GRanges(seqnames = pcg$chromosome_name, ranges = IRanges(start = pcg$start_position, end = pcg$end_position), ensembl_gene_id = pcg$ensembl_gene_id)

tmp2 <- as.data.frame(findOverlaps(pcgGR, tmpGR))
pcgGR[tmp2@from]$ensembl_gene_id
tmpGR[tmp2@to]$freq

summ <- as.data.frame(matrix(nrow = nrow(tmp2), ncol = 2))
colnames(summ) <- c("gene","freq")
summ$gene <- pcgGR[tmp2$queryHits]$ensembl_gene_id
summ$freq <- tmpGR[tmp2$subjectHits]$freq

b <- summ %>% group_by(gene) %>%
  summarise(freq = sum(freq))
write.csv(b, "Num_ctotmut_per_gene.csv")

b <- fread("Num_ctotmut_per_gene.csv")

gen <- BSgenome.Hsapiens.UCSC.hg19
pcg <- pcg %>% mutate(chromosome_name = paste0("chr",chromosome_name))
pcgGR = GRanges(seqnames = pcg$chromosome_name, ranges = IRanges(start = pcg$start_position, end = pcg$end_position), ensembl_gene_id = pcg$ensembl_gene_id)
gap_a <- getSeq(gen, pcgGR)
gap_tt <- as.data.frame(letterFrequency(gap_a, letters = c("C","G")))
gap_tt$gene <- pcgGR$ensembl_gene_id
gap_tt$gene_width <- pcgGR@ranges@width
gap_tt$chrom <- pcg$chromosome_name
gap_tt$start <- pcgGR@ranges@start
gap_tt$end <- gap_tt$start + gap_tt$gene_width - 1

summ <- b %>% full_join(gap_tt, by="gene") %>%
  mutate(mutation_freq = freq/(C+G))

plot_mutation_frequency <- function(table1, table2, table3, name1, name2,name3, ypos1, ypos2, ypos3, ylimit, outpath, w, h) {
  tab1 <- fread(table1, select = 9)
  tab2 <- fread(table2, select = 9)
  tab3 <- fread(table3, select = 5) #%>% 
  
  top_gene <- summ %>% filter(gene %in% tab1$ensembl_gene_id) %>% drop_na() %>%
    mutate(group = name1)
  bottom_gene <- summ %>% filter(gene %in% tab2$ensembl_gene_id) %>% drop_na() %>%
    mutate(group = name2)
  middle_gene <- summ %>% filter(gene %in% tab3$ensembl_gene_id) %>% drop_na() %>%
    mutate(group = name3)
  
  plotting <- rbind(top_gene, bottom_gene, middle_gene)
  plotting$group <- factor(plotting$group, levels=c(name1, name2, name3))
  ggplot(plotting, aes(x = group, y= mutation_freq, fill = group)) + geom_boxplot(outlier.colour = NA) + 
    geom_signif(comparisons = list(c(name1, name2)), map_signif_level=TRUE, y_position = ypos1) +
    geom_signif(comparisons = list(c(name3, name2)), map_signif_level=TRUE, y_position = ypos2) +
    geom_signif(comparisons = list(c(name1, name3)), map_signif_level=TRUE, y_position = ypos3) +
    xlab("Susceptibility") + ylab("Melanoma mutation frequency") + ylim(ylimit) + 
    scale_fill_manual(values = c("#CC6677","#44AA99","#888888")) +
    theme_bw() + theme(legend.position = "None", axis.title = element_text(face = "bold"),
                       axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(outpath, width = w, height = h,  dpi=600)
}

plot_mutation_frequency("../plot/table1_gene_overlap_with_bins/top10per_overlap10per_v2.csv",
                        "../plot/table1_gene_overlap_with_bins/bottom10per_overlap10per_v2.csv",
                        "../plot/table1_gene_overlap_with_bins/table1_v2.csv",
                        "Top 10%", "Bottom 10%", "Whole genome", 0.037,0.025,0.045, ylimit = c(0,0.05),
                        "../plot/Mutation_frequency_10percent.pdf", 2,4)

plot_mutation_frequency("../plot/table1_gene_overlap_with_bins/top5per_overlap10per_v2.csv",
                        "../plot/table1_gene_overlap_with_bins/bottom5per_overlap10per_v2.csv",
                        "../plot/table1_gene_overlap_with_bins/table1_v2.csv",
                        "Top 5%", "Bottom 5%", "Whole genome", 0.037,0.025,0.045, ylimit = c(0,0.05),                   
                        "../plot/Mutation_frequency_5percent.pdf", 2,4)
