setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
install_pkg(c("data.table", "dplyr", "ggplot2", "ggrepel", "ggbeeswarm"),"common")
install_pkg(c("biomaRt"), "bioconductor")

'%ni%' <- Negate('%in%')

## getting list of CDGs and their coordinates
cdg1 <- fread("data/cancer_driver_gene/Bailey2018.csv", select = 1, col.names = "gene")
cdg2 <- fread("data/cancer_driver_gene/COSMIC_CGS.csv", select = 1, col.names = "gene")
cdg3 <- fread("data/cancer_driver_gene/IntOGen-DriverGenesv2020-02-01.tsv", select = 1, col.names = "gene")
summ <- unique(unlist(cdg1 %>%  union(cdg2) %>% union(cdg3)))
summ <- c(summ,"TPTE")
hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
allgene <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = "biotype", values = "protein_coding", mart = hg19) %>%
  dplyr::filter(!grepl("H", chromosome_name)) %>% 
  dplyr::filter(!grepl("G", chromosome_name)) %>%
  dplyr::filter(!grepl("M", chromosome_name)) %>% 
  mutate(chromosome_name = paste0("chr", chromosome_name))

cdg_summ <- allgene %>%
  dplyr::filter(hgnc_symbol %in% summ)
reg_summ <- allgene %>%
  dplyr::filter(hgnc_symbol %ni% summ)

## select centromeric/telomeric region
cent <- fread("data/hg19gap.txt") %>% 
  dplyr::filter(type == "centromere")
telo <- fread("data/hg19gap.txt") %>% 
  dplyr::filter(type == "telomere")
ct_summ <- rbind(cent, telo) %>% 
  dplyr::select(chrom, chromStart, chromEnd)
setorder(ct_summ, chrom, chromStart)


dist_to_centrotelo <- function(sn, st, en) {
  # case1
  selected <- ct_summ %>% 
    dplyr::filter(chrom== sn) 
  subselected <- selected %>%
    dplyr::filter((chromStart <= st & chromEnd >= en) | (chromStart >= st & chromStart <= en) | (chromEnd >= st & chromEnd <= en)) 
  if(nrow(subselected) > 0) {
    return(0)
  }
  sub1 <- selected %>% dplyr::filter(chromEnd <= st)
  sub2 <- selected %>% dplyr::filter(chromStart >= en)
  if(nrow(sub1) == 0) {
    return(sub2$chromStart[1] - en )
  } else if (nrow(sub2) == 0 ) {
    return(st - sub1$chromEnd[nrow(sub1)])
  } else {
    return(min((sub2$chromStart[1] - en), (st - sub1$chromEnd[nrow(sub1)]) ))
  }
}

cdg_summ$dist <- mapply(dist_to_centrotelo, cdg_summ$chromosome_name, cdg_summ$start_position, cdg_summ$end_position)
reg_summ$dist <- mapply(dist_to_centrotelo, reg_summ$chromosome_name, reg_summ$start_position, reg_summ$end_position)
cdg_summ$group = "CDG"
reg_summ$group = "All others"
plotting <- rbind(cdg_summ, reg_summ)
ggplot(plotting, aes(x = group, y = dist) ) + geom_boxplot()

### individual investigation
allgene$dist <- mapply(dist_to_centrotelo, allgene$chromosome_name, allgene$start_position, allgene$end_position)
gene_care <- c("ALK", "APC","ARID1A", "ARID2", "ATM", "ATRX", "BCOR", "BRAF", "BRCA1", "BRCA2",
               "CCND1", "CDK4", "CDKN2A", "EZH2", "IDH1", "IDH2", "JAK1", "KRAS", "MECOM", "MTOR",
               "MYC", "NF1", "NOTCH2", "NRAS", "PTEN","RB1", "SMAD4", "TP53", "TERT", "HRAS", "AXIN1", "TPTE")

plot_cdg <- function(df, xlabel, outpath ) {
  df <- df %>% mutate(plotting = "") 
  df$plotting[which(df$hgnc_symbol %in% gene_care)] = df$hgnc_symbol[which(df$hgnc_symbol %in% gene_care)]
  ggplot(df, aes(x=1, y=dist)) + geom_quasirandom(size = 0.2) + 
    xlim(0.5,1.5) + geom_hline(yintercept=median(allgene$dist),linetype="dashed") +
    geom_label_repel(df[which(df$plotting!=""),], mapping = aes(label=plotting), size = 2, max.overlaps = Inf) + 
    theme_bw() + xlab(xlabel) + ylab("Distance from Centromer or Telomere") + 
    annotate("text", label = "Genome Median", x = 0.62, y = 18000000, size = 2.5) +
    theme( axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
           axis.title = element_text(face = "bold"))# + 
    #scale_color_manual(values = c("#CC6677","#888888","#44AA99"))
  ggsave(outpath, width = 4.5, height = 6, dpi = 600)
}
cdg1 <- fread("data/cancer_driver_gene/Bailey2018.csv") 
cdg1plot <-  allgene %>% dplyr::filter(hgnc_symbol %in% cdg1$Gene) 
plot_cdg(cdg1plot, "Cancer Driver Genes (Bailey 2018)", "plot/bailey2018.pdf")

cdg2 <- fread("data/cancer_driver_gene/COSMIC_CGS.csv") 
cdg2plot <-  allgene %>% filter(hgnc_symbol %in% cdg2$Gene) 
plot_cdg(cdg2plot, "Cancer Driver Genes (COSMIC)", "plot/COSMIC.pdf")

cdg3 <- fread("data/cancer_driver_gene/IntOGen-DriverGenesv2020-02-01.tsv")
cdg3plot <-  allgene %>% filter(hgnc_symbol %in% cdg3$Symbol) 
plot_cdg(cdg2plot, "Cancer Driver Genes (Intogen)", "plot/Intogen.pdf")

cdgall <- allgene %>% dplyr::filter(hgnc_symbol %in% summ) 
top <- fread("plot/table1_gene_overlap_with_bins/top10per_overlap10per_v2.csv")
bottom <- fread("plot/table1_gene_overlap_with_bins/bottom10per_overlap10per_v2.csv")
cdgall$Group <- "Others"
cdgall$Group[which(cdgall$hgnc_symbol %in% top$hgnc_symbol )] <- "Top 10%"
cdgall$Group[which(cdgall$hgnc_symbol %in% bottom$hgnc_symbol )] <- "Bottom 10%"

plot_cdg(cdgall, "Cancer Driver Genes (Compiled)", "plot/Compiled.pdf")
0
