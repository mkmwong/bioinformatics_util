setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("ChIPseeker","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
         "ReactomePA","org.Hs.eg.db")
install_pkg(pkg,"bioconductor")
install_pkg("data.table","common")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
quies <- fread("SRR034480.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 >= -log10(0.05))
grow <- read.table("SRR034478.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 > -log10(0.05))

growGR <- GRanges(seqnames =  grow$V1, ranges = IRanges(start = grow$V2 + 1, end = grow$V3))
quiesGR <- GRanges(seqnames =  quies$V1, ranges = IRanges(start = quies$V2 + 1, end = quies$V3))

export_annopie(growGR, txdb, "../plot/grow_RbPeak", TRUE)
export_annopie(quiesGR, txdb, "../plot/quies_RbPeak", TRUE)

enrich_path_dotplot(growGR, "../plot/grow_RbPeak_promoter_pathway", TRUE,"Promoter")
enrich_path_dotplot(quiesGR, "../plot/quies_RbPeak_promoter_pathway", TRUE, "Promoter")
