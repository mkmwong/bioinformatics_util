setwd("~/Desktop/yr5/summ/bioinformatics_util/")
source("draft.R")
pkg = c("Sushi","biomaRt")
install_pkg(pkg, "bioconductor")
data(Sushi_genes.bed)

chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000
chrom_biomart    = 15
plotGenes(Sushi_genes.bed,chrom_biomart,chromstart,chromend ,types=Sushi_genes.bed$type,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=FALSE,col="blue",
          labeloffset=1,fontsize=1.2)

tmp = fread("data/rb_geneplot_subset.bed") %>% 
  filter(V2 != "gene") %>% filter(V2 != "transcript") %>%
  filter(V2 != "exon")
tmp$gene = "Rb1"
tmp = tmp[,c(1,3,4,7,5,6,2)]
colnames(tmp) = c("chrom","start","stop","gene","score","strand","type")
tmp$start = as.numeric(tmp$start)
tmp$stop = as.numeric(tmp$stop)
#tmp$chrom = substring(tmp$chrom,4)
chrom = "chr13"
chromstart = 48877000
chromend = 49057000
chrom_biomart = 13
plotGenes(tmp,chrom = chrom_biomart,chromstart = chromstart,chromend = chromend ,types=tmp$type,
          maxrows=2,height=0.05,plotgenetype="box",bentline=FALSE,col="brown",
          labeloffset=0.2,fontsize=1, labeltext = TRUE)
labelgenome( chrom, chromstart,chromend,n=3,scale="Mb")
