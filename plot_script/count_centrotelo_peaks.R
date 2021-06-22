### get number of RB peaks that binds to pericentric/subtelomeric regions.

gap <- fread("../hg19gap.txt", stringsAsFactors = FALSE)
centroGR <- get_gap_GR(gap, "centromere", 1000000)
teloGR <- get_gap_GR(gap, "centromere", 100000)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
quies <- fread("../SRR034480.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 >= -log10(0.05))
grow <- read.table("../SRR034478.merged.nodup.pooled_x_SRR034492.merged.nodup.pooled.pval0.01.500K.bfilt.narrowPeak") %>%
  filter(V9 > -log10(0.05))

growGR <- GRanges(seqnames =  grow$V1, ranges = IRanges(start = grow$V2 + 1, end = grow$V3))
quiesGR <- GRanges(seqnames =  quies$V1, ranges = IRanges(start = quies$V2 + 1, end = quies$V3))
peakAnnoGrow <- annotatePeak(growGR, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
#peakAnnoGrow <- peakAnnoGrow@anno[peakAnnoGrow@anno$annotation == "Distal Intergenic"]
sum(countOverlaps(peakAnnoGrow@anno,centroGR))
sum(countOverlaps(peakAnnoGrow@anno,teloGR))
peakAnnoQuies <- annotatePeak(quiesGR, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
sum(countOverlaps(peakAnnoQuies@anno,centroGR))
sum(countOverlaps(peakAnnoQuies@anno,teloGR))
