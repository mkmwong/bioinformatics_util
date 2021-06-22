setwd("~/Desktop/yr5/summ/bioinformatics_util/")
pathi = "data/stjude_retinoblastoma_snv/tert/"
a<- list.files(path=pathi)
starter <- fread(paste0(pathi,a[1]), sep="\t") %>% 
  separate(V5, into = c("mutated_to", "type"), sep=",") %>%
  separate(V1, into = c("not_used","chrom"), sep="r") %>% 
  mutate(filename = a[1]) %>% 
  separate(V10, into = c("genotype"), sep=":")
for (i in 2:length(a)) {
  print(i)
  tmp <- fread(paste0(pathi,a[i]), sep="\t") %>% 
    separate(V5, into = c("mutated_to", "type"), sep=",") %>%
    separate(V1, into = c("not_used","chrom"), sep="r") %>% 
    mutate(filename = a[i])  %>% 
    separate(V10, into = c("genotype"), sep=":")
  starter <- rbind(starter, tmp)
}
starter <- starter %>% mutate(end = V2 + nchar(V4)-1) %>%
  select(chrom, V2, end,  V4, mutated_to, filename, genotype)
write.table(starter, "data/stjude_retinoblastoma_snv/tert_for_annovar.tsv", quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE  )


View(starter %>% filter(nchar(V4) <=2))
