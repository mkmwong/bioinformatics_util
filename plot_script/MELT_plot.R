### Changing directory and installing/loading all required packages
setwd("~/Desktop/work_tmp/retinoblastoma/")
source("~/Desktop/work_tmp/bioinformatics_util/draft.R")
pkg <- c("ggplot2","ggpubr","tidyverse","dplyr","data.table","rstatix")
install_pkg(pkg, "common")

### Importing metadata from the sequencing data
tmp <- fread("SAMPLE_INFO.txt") %>% 
  filter(sequencing_type == "WGS") %>%
  filter( !grepl("bai",file_path)) %>%
  filter( file_type == "BAM") %>%
  filter( sj_datasets != "Clinical Pilot")
all_f <- list.dirs("samp_out/",recursive = FALSE, full.names= FALSE)
r_type <- c("LINE1","ALU","SVA")

### Function to return count of each type of insertion per sample and to 
### parse the INFO column of .vcf file. 
### Filtering recommended by MELT in MELT google group was used.
### Specifying out to return either count for each type of retrotransposon 
### (out == "COUNT") or to return count for each region in a specific type 
### of retrotransposon (out == "REGION_COUNT")
MELT_count <- function(type, dir, out) {
  ret = NULL
  vcf <- fread(paste0(dir,"/", type,".final_comp.vcf")) %>%
    filter(FILTER == "PASS") 
  if(type == "LINE1") {
    vcf <- vcf %>%
      extract(INFO, into =c("TSD","ASSESS","INTERNAL", "SVTYPE", "SVLEN","MEINFO","DIFF", "LP", "RP","RA","ISTP","PRIOR","SR"),
              regex="=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*)")
  }
  else {
    vcf <- vcf %>% 
      extract(INFO, into =c("TSD","ASSESS","INTERNAL", "SVTYPE", "SVLEN","MEINFO","DIFF", "LP", "RP","RA","PRIOR","SR"),
                       regex="=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*);[:alpha:]+=([^;]*)")
  }
  vcf <- vcf %>%  
    mutate(ASSESS = as.numeric(ASSESS), SR = as.numeric(SR),
           SVLEN = as.numeric(SVLEN)) %>%
    filter(ASSESS >= 3) %>% 
    filter(SR >= 3) %>% 
    filter(PRIOR == "false") %>%
    filter(SVLEN != 0 )
  if (out == "COUNT") {
    ret = length(vcf$POS)
  }
  else {
    vcf <- vcf %>% 
      separate(INTERNAL, sep=",", into=c("GENE","LOC"))
    ret = as.list(table(vcf$LOC))
  }
  return(ret)
}

### Function to make boxplot or density plot with summary of .vcf file 
### that is output from MELT. 
MELT_plot <- function(summary, plot_type, outpath, group1, group2, y.position=NULL, 
                      xlab=NULL, ylab = NULL, de = "pdf", w = 6, h = 3.5, dp = 600) {
  stat_test <- summary %>%
    group_by(repeat_type) %>% 
    wilcox_test(count ~ group, exact = FALSE, paired = TRUE, 
                p.adjust.method = "BH") %>%
    add_significance() %>% 
    mutate(group1 = group1, group2 = group2) %>% 
    add_xy_position(x = "group")
  if(length(y.position) > 0){
    stat_test$y.position = y.position
  }
  bold_text <- element_text(face = "bold")
  
  if(plot_type == "density") {
    ggplot() + geom_density(data = summary, aes(x=count, group = group, fill = group ),alpha = 0.5) +
      theme_bw() + xlab (xlab) + ylab(ylab) + scale_fill_discrete(name = "Group") + 
      theme(axis.title = bold_text, legend.title = bold_text, strip.text = bold_text ) + 
      facet_wrap(~repeat_type, scale = "free", nrow = 1)
    ggsave(paste0(outpath,".",de), device = de , width = w, height = h, dpi = dp )
  }
  ### boxplot
  else {
    print(head(stat_test))
    ggplot(summary, aes(x = group, y = count, fill = group, group = group )) + geom_violin() +
      geom_boxplot(width = 0.1, outlier.shape = NA) + facet_wrap(~repeat_type, scale = "free", nrow = 1) + 
      theme_bw() + ylab(ylab) + xlab(xlab) + scale_fill_discrete(name = "Group") + 
        theme(axis.title = bold_text, legend.title = bold_text, strip.text = bold_text ) + 
        stat_pvalue_manual(stat_test, inherit.aes = FALSE )
    ggsave(paste0(outpath,".",de), device = de , width = w, height = h, dpi = dp)
    
  }
}

### Getting count for each type of retrotransposon insertion per sample
count_df <- as.data.frame(sapply(r_type, function(i) sapply(all_f, function(j) MELT_count(i,paste0("samp_out/",j),"COUNT"))))
count_df <- count_df %>% 
  tibble::rownames_to_column(var = "sample_name") %>% 
  inner_join(tmp, by = "sample_name") %>%
  mutate(LINE1 = as.numeric(LINE1), ALU = as.numeric(ALU), SVA = as.numeric(SVA)) 

### Separating tumor sample from control sample
diagnosis <- count_df %>% 
  filter(sample_type == "Diagnosis") %>% 
  select(subject_name, LINE1, ALU, SVA) %>% 
  gather( key = "repeat_type", value = "count", -subject_name) %>%
  mutate(group = "Tumor")
germline <- count_df %>% 
  filter(sample_type=="Germline")
### Only keeping subject that has both tumor and control sample for paired comparison 
germline <- germline[which(germline$subject_name %in% diagnosis$subject_name ),] 
germline <- germline %>% 
  select(subject_name, LINE1, ALU, SVA) %>% 
  gather( key = "repeat_type", value = "count", -subject_name) %>%
  mutate(group = "Control")

### Merging sample and creating plots
summary <- rbind(germline, diagnosis) 
summary$group <- factor(summary$group, levels = c("Control","Tumor"))
MELT_plot(summary,"density","density_plot","Tumor","Control", NULL, 
          "Count of repeat types", "Density", "pdf", 7, 3, 600 )
MELT_plot(summary,"boxplot","box_plot","Tumor","Control", c(1000,125,40), 
          "", "Count of repeat types", "pdf", 6, 3.5, 600 )

### Since only ALU has a significant difference between tumor and control,
### we move forward with only ALU. However, the following analysis can be 
### applied to SVA or LINE1 by changing the type parameter on the next line.
ALU <- sapply(all_f, function(j) MELT_count("ALU",paste0("samp_out/",j), "REGION_COUNT"))
tmp_name <- names(ALU)
ALU <- bind_rows(ALU) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  mutate(across(everything(), ~replace_na(.x,0))) %>%
  mutate(EXON = rowSums(select(., starts_with("EXON_")), na.rm = TRUE)) %>%
  mutate(sample_name = tmp_name) %>% 
  inner_join(tmp, by = "sample_name") %>% 
  select(-starts_with("EXON_"))

### Separating tumor sample from control sample
diagnosis <- ALU %>% 
  filter(sample_type == "Diagnosis") %>% 
  select("3_UTR", INTRONIC, null, PROMOTER, TERMINATOR, "5_UTR", EXON, subject_name) %>%
  gather( key = "repeat_type", value = "count", -subject_name) %>%
  mutate(group = "Tumor")
germline <- ALU %>% 
  filter(sample_type=="Germline")
germline <- germline[which(germline$subject_name %in% diagnosis$subject_name ),] 
germline <- germline %>% 
  select("3_UTR", INTRONIC, null, PROMOTER, TERMINATOR, "5_UTR", EXON, subject_name) %>%
  gather( key = "repeat_type", value = "count", -subject_name) %>%
  mutate(group = "Control")

### Merging sample and creating plots
summary <- rbind(germline, diagnosis) 
summary$group <- factor(summary$group, levels = c("Control","Tumor"))
MELT_plot(summary,"boxplot","box_plot_bytype","Tumor","Control", c(10,5,3,450,520,50,40), 
          "", "Count of repeat types", "pdf", 10, 3.5, 600 )






