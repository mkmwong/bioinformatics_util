setwd("~/Desktop/yr5/summ/")
source("bioinformatics_util/draft.R")

plot_DAVID_enrichment <- function(df, top_path, bottom_path) { 
  df$Category <- factor(df$Category, levels = c("GO","KEGG","IPR"))
  up_df <- df %>% filter(`signed enrichment score` > 0 )
  down_df <- df %>% filter(`signed enrichment score` < 0 )
  #down_df <- down_df %>% mutate(`signed enrichment score` = -`signed enrichment score`)
  ggplot(up_df, aes(x = reorder(term, `signed enrichment score`) , y = `signed enrichment score`, fill = Category)) + geom_bar(stat = "identity") + 
    coord_flip() + theme_bw() +  xlab("Gene Ontology Term") + ylab("Signed Enrichment Score") + 
    theme(legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) + 
    geom_hline(yintercept = 1.3, linetype="dashed") + scale_x_discrete(position = "top") + scale_fill_discrete(drop=FALSE) + 
    ylim(-10,18)
  ggsave(top_path, w=7, h= 2, dpi = 600)
  ggplot(down_df, aes(x = reorder(term, `signed enrichment score`) , y = `signed enrichment score`, fill = Category)) + geom_bar(stat = "identity") + 
    coord_flip() + theme_bw() +  xlab("Gene Ontology Term") + ylab("Signed Enrichment Score") + 
    theme(legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) + 
    geom_hline(yintercept = -1.3, linetype="dashed") + scale_x_discrete(position = "top") + 
    scale_y_continuous(trans = "reverse") + scale_fill_discrete(drop=FALSE) + ylim(-10,18)
  ggsave(bottom_path, w=7, h= 2, dpi = 600)
}

plot_DAVID_enrichment_oneimage <- function(df, out_path, breaks ) { 
  df$Category <- factor(df$Category, levels = c("GO","KEGG","GAD","KW"))
  #up_df <- df %>% filter(`signed enrichment score` > 0 )
  #down_df <- df %>% filter(`signed enrichment score` < 0 )
  #down_df <- down_df %>% mutate(`signed enrichment score` = -`signed enrichment score`)
  df$dir = "Up"
  df$dir[which(df$`signed enrichment score` <= 0)] = "Down"
  ggplot(df, aes(x = reorder(term, `signed enrichment score`) , y = `signed enrichment score`, fill = Category)) + geom_bar(stat = "identity") + 
    coord_flip() + theme_bw() +  xlab("Gene Ontology Term") + ylab("Signed Enrichment Score") + 
    theme(legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) + 
    geom_hline(yintercept = 1.3, linetype="dashed") +  geom_hline(yintercept = -1.3, linetype="dashed") + 
    scale_x_discrete(position = "top") + scale_fill_discrete(drop=FALSE) +
    scale_y_continuous(name="Stopping distance", breaks=breaks) + facet_wrap(~dir)
  ggsave(out_path, w=7, h= 3, dpi = 600)
  #ggplot(down_df, aes(x = reorder(term, `signed enrichment score`) , y = `signed enrichment score`, fill = Category)) + geom_bar(stat = "identity") + 
  #  coord_flip() + theme_bw() +  xlab("Gene Ontology Term") + ylab("Signed Enrichment Score") + 
  #  theme(legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) + 
  #  geom_hline(yintercept = 1.3, linetype="dashed") +  scale_y_continuous(trans = "reverse") + scale_fill_discrete(drop=FALSE)
  #ggsave(bottom_path, w=7, h= 2, dpi = 600)
}

sus_gene <- fread("bioinformatics_util/data/ashbys_analysis/Susceptibility_gene.csv") %>%
  separate(V1, sep=":", into = c("Category", "term"))

plot_DAVID_enrichment(sus_gene,
                      "bioinformatics_util/plot/top_susceptible.pdf",
                      "bioinformatics_util/plot/bottom_susceptible.pdf")

#plot_DAVID_enrichment_oneimage(sus_gene, "bioinformatics_util/plot/susceptible.pdf", breaks = c(-2,0,2,4,6,8))

rnaseq <- fread("bioinformatics_util/data/ashbys_analysis/RNASeq_plot.csv") %>%
  separate(V1, sep=":", into = c("Category", "term"))

plot_DAVID_enrichment(rnaseq,
                      "bioinformatics_util/plot/top_rnaseq.pdf",
                      "bioinformatics_util/plot/bottom_rnaseq.pdf")
