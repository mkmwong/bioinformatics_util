#### install_pkg                                                           ####
#### This is a supporting function to check for other functions if all     ####
#### their dependencies has been installed. If not, install them.          ####
#### parameters: pkg; list of packages required for the function           ####
####            pkg_type; common for CRAN packages, bioconductor for       ####
####                      bioconductor packages.                           ####
#### return: NULL.                                                         ####
install_pkg <- function(pkg, pkg_type) {
  uninstalled <- which( !pkg %in% rownames(installed.packages()) )
  if( length(uninstalled) ) {
    print("One or several required package(s) are not installed. Installing the package(s)...")
    if (pkg_type == "common") {
      install.packages(pkg[uninstalled])
    }
    else {
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        print("BiocManager is not installed. Installing BiocManager...")
        install.packages("BiocManager")
      }
      BiocManager::install(pkg[uninstalled])
    }
  }
  else {
    print("All required packages are installed.")
  }
  lapply(pkg, library, character.only = TRUE)
  return(NULL)
}

#### remove_chrom                                                          ####
#### This is a pre-processing step prior to running binning_GR.            ####
#### This function should remove all ranges that are in specific           ####
#### chromosome.                                                           ####
#### parameters:input; could be either a data.frame object or a            ####
####                   SimpleRleList object.                               ####
####            value: the chromosome to be removed.                       ####
####            column(optional): the column where the chromosome info is  ####
####                              stored.                                  ####
#### return: input with ranges of the specific chromosome removed.         ####
remove_chrom <- function(input, value, column=NULL) {
  pkg <- c("dplyr")
  install_pkg(pkg,"common")
  if(class(input) == "data.frame") {
    args <- as.list(match.call())
    tmp <- input %>% filter((!!sym(column)) != value)
  }
  if(class(input) == "SimpleRleList") {
    tmp <- input
    tmp@listData[value] = NULL
  }
  return(tmp)
}

#### binning_GR                                                            ####
#### This is a function to bin any bigwig file to the size of supplied     ####
#### bins. The output could be saved as rds (and csv if sepcified).        ####
#### parameters: binsGR; GRanges object defining the bins to bin the data  ####
####                     with.                                             ####
####             datrle; Rlelist object that needs to be binned by binsGR. ####
####             csv: Boolean variable indicating whether or not .csv file ####
####                  should be exported in addition to .rds file.         ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            #### 
#### return: binned GRanges object.                                        ####
# example usage:
# bins <- read.table("~/Desktop/tmp/Old_Stuffs/final_figures/data/hg19gap.txt", sep="\t")
# bins <- remove_chrom(bins, "chrM", "V2")
# binsGR <- GRanges(seqnames = bins$V1, ranges = IRanges(bins$V2+1, bins$V3), stat = bins$V4)
# dat <- import.bw(con="~/Desktop/tmp/Old_Stuffs/final_figures/data/unbinned_bigwig/wgEncodeFsuRepliChipImr90WaveSignalRep1.bigWig", as='RleList')
# dat[dat == 0] <- NA
# dat <- remove_chrom(dat, "chrM")
# a <- binning_GR(binsGR, dat, TRUE,"~/Desktop/something.csv")

binning_GR <- function(binsGR, datrle, csv, outpath) {

  bio_pkg <- c("rtracklayer","GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19")
  install_pkg(bio_pkg,"bioconductor")
  gen <- BSgenome.Hsapiens.UCSC.hg19
  si.gen <- seqinfo(gen)
  si <- si.gen[names(datrle)]
  seqlevels(binsGR) = names(datrle)
  ba <- binnedAverage(binsGR,numvar=datrle,varname='average', na.rm = TRUE)
  save_files(as.data.frame(ba), csv, outpath)
  return(ba)
}


#### merge_df                                                              ####
#### This is a function to merge all files in a directory (some delimited  ####
#### files ) into one single dataframe. The output would be saved as a     #### 
#### variable, and would be saved as rds (and csv if specified).           ####
#### parameters: directory: A string of the path for the directory that    ####
####                        contains all the files that need to be merged. ####
####             by_col: list of column names to use for merging.          ####
####             sig_col: numeric value of the column which the column     ####
####                      name is used as the name of column after merging.####
####             csv: Boolean variable indicating whether or not .csv file ####
####                  should be exported in addition to .rds file.         ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            #### 
####             select_col(optional): Columns in each file to be included ####
#### returned: A dataframe object contains all files from the specified    ####
####           directory merged together.                                  ####
merge_df <- function(directory, by_col, sig_col, csv, outpath, select_col=NULL) {
  pkg <- c("data.table", "stringr")
  install_pkg(pkg, "common")
  all_files <- list.files(path = directory)
  path_files <- paste0(directory, all_files)
  names <- sapply(lapply(all_files, str_match,"-\\s*(.*?)\\s*\\."), "[[", 2)
  all_df <- lapply(path_files, fread, select=select_col)
  all_df <- mapply( function(x,y) {colnames(x)[sig_col] <- y; x}, all_df, names, SIMPLIFY = FALSE)
  merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE, by=by_col), all_df)
  save_files(merge_df, csv, outpath)
  return(merge_df)
}
# example usage
# a=merge_df("~/Desktop/work_tmp/bioinformatics_util/100kb/", c("seqnames","start","end"), 
#           4, TRUE, "~/Desktop/summary", c(2,3,4,7))

#### save_files                                                            ####
#### This is a function to save a dataframe as .rds file, with the option  ####
#### to save the datafraame as a .csv file.                                ####
#### parameters: df: A dataframe object to be saved.                       ####
####             csv: Boolean variable indicating whether or not .csv file ####
####                  should be exported in addition to .rds file.         ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            ####
#### return: NULL.                                                         ####
save_files <- function(df, csv, outpath) {
  if (csv == TRUE ) {
    write.csv(df,paste0(outpath,".csv"))
  }
  saveRDS(df, paste0(outpath,".rds"))
}

#### pca_plot                                                              ####
#### This is a function to perform PCA on a named matrix, and returning    ####
#### principal component loading vectors for graphing.                     ####
#### parameters: df: A dataframe object to preform principal component     ####
####                 analysis on.                                          ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            ####
####             col: Color for each variable in plot                      ####
####             marg(optional): Multiplier to set the xmin, xmax, ymin,   ####
####                             ymax for plotting                         ####
####             w,h,d,un,de(all optional): width, height, dpi, unit and   ####
####                                        device for saving the plot.    ####
####                                                                       ####
#### return: A prcomp object from PCA.                                     ####
pca_plot <- function(df,  outpath, col, marg=0.2, w = 5, h = 5, d = 600, un = "in", de = "pdf" ) {
  pkg <- c("ggrepel","ggplot2")
  install_pkg(pkg, "common")
  all_pca <- prcomp(df, center = TRUE, scale. = TRUE)
  all_out <- as.data.frame(all_pca$rotation)
  all_out$feature <- row.names(all_out)
  var <- all_pca$sdev^2
  var_ex <- var/sum(var)
  x_int <- max(all_out$PC1) - min(all_out$PC1)
  y_int <- max(all_out$PC2) - min(all_out$PC2)
  x_max <- max(all_out$PC1) + x_int * marg
  x_min <- min(all_out$PC1) - x_int * marg
  y_max <- max(all_out$PC2) + y_int * marg
  y_min <-  min(all_out$PC2) - x_int * marg
  ggplot(all_out,aes(x=PC1,y=PC2,label=feature )) + geom_point(color = col)  +
    xlab(paste("PC1 (var. explained =", round(var_ex[1],2), ")")) + 
    ylab(paste("PC2 (var. explained =", round(var_ex[2],2), ")")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + theme_bw() + 
    theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
          axis.title.y = element_text(color = "black", size = 14, face = "bold")) + geom_text_repel(size=3) +
    xlim(x_min,x_max) + ylim(y_min, y_max)
  ggsave(paste0(outpath,".",de), width=w, height=h, dpi=d, units=un, device=de)
  save_files(all_pca,FALSE, outpath)
  return(all_pca)
}

#### pval_lab                                                              ####
#### This is a function to return corresponding p value labels for         ####
#### plotting purpose.                                                     ####
#### parameters: pval: p-value to convert to label.                        ####
#### return: A string, the label.                                          ####
pval_lab <- function(pval) {
  lab <- NULL
  if ( pval < 0.05 & pval > 0.01) {
    lab <- "p < 0.05"
  } else if ( pval < 0.01 & pval > 0.001) {
      lab <- "p < 0.01"
  }
  else if (pval < 0.001) {
    lab <- "p < 0.001"
  }
  else {
    lab <- "p > 0.05"
  }
  return(lab)
}

#### corr_plot                                                             ####
#### This function makes a scatter plot from the two supplied lists of     ####
#### data. Correlation coefficient and p-value can be plotted optionally.  ####                                    ####
#### parameters: col1, col2: Two lists of data for plotting scatter plot.  ####
####                         The two list must have same length.           ####
####             xlab, ylab: String; X-axis label and y-axis label for the ####
####                         plot.                                         ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            ####
####             marg(optional): Multiplier to set the xmin, xmax, ymin,   ####
####                             ymax for plotting                         ####
####             w,h,d,un,de(all optional): width, height, dpi, unit and   ####
####                                        device for saving the plot.    ####
####                                                                       ####
####             addLab(optional): Boolean, whether or not to add R and    ####
####                               p-value to the plot. Default is FALSE.  ####                                   ####
####             x1,y1,x2,y2: x and y coordinates for R and p-value        ####
####                          respectively.                                ####
#### return: NULL.                                                         ####
corr_plot <- function(col1, col2, xlab, ylab, outpath, marg = 0.2, 
                      w = 3, h = 2.5, un = "in", d = 600, de = "pdf",
                      addLab=FALSE, x1=NULL, x2 = NULL, y1 = NULL, y2 = NULL) {
  if (length(col1) != length(col2) ) {
    stop("Length of two supplied list must be equal!")
  }
  pkg <- c("ggplot2")
  install_pkg(pkg, "common")
  x_int <- max(col1) - min(col1)
  y_int <- max(col2) - min(col2)
  x_max <- max(col1) + x_int * marg
  x_min <- min(col1) - x_int * marg
  y_max <- max(col2) + y_int * marg
  y_min <-  min(col2) - x_int * marg
  a = ggplot(data = as.data.frame(cbind(col1, col2)), 
             aes(x = col1, y = col2)) + geom_point(alpha = 0.1) +
    xlim(x_min, x_max) + ylim(y_min, y_max) + geom_smooth(method = "lm") + theme_bw() +
    xlab(xlab) + ylab(ylab)
  if( addLab) {
    corre = paste("r =", round(cor.test(col1, col2, method = "spearman", exact = FALSE)$estimate,2))
    corre.pval = pval_lab(cor.test(col1, col2)$p.val)
    a = a + annotate(geom="text",x = x1, y = y1, label=corre,size=4, fontface="bold") + 
    annotate(geom="text",x = x2, y = y2, label=corre.pval, size=4)
  }
  ggsave(paste0(outpath,".",de) ,width = w, height = h, unit = un, dpi = d, device = de)
  return(NULL)
}

## function to return coordinate of genomic range of given feature of gap file 
## gap file optained from ucsc table browser
get_gap_GR <- function(gap, type_gap, width){
  tmp <- gap %>% filter(type == type_gap)
  ret <- as.data.frame(matrix(ncol = 3, nrow = nrow(tmp)*2 )) %>% 
    mutate(V1 = rep(tmp$chrom, 2),
           V2 = c((tmp$chromStart-width + 1), tmp$chromEnd + 1),
           V3 = c(tmp$chromStart, (tmp$chromEnd + width)))
  retGR <- GRanges(seqnames= ret$V1, ranges = IRanges(start = ret$V2, end = ret$V3))
  return(retGR)
}

### function to export anno pie plot directly from genomic ranges
export_annopie <- function(gr,txdb,outpath, writefile, h= 4, w = 8) {
  peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  pdf(file = paste0(outpath,".pdf"), height = h, width = w)
  print(plotAnnoPie(peakAnno))
  dev.off()
  if (writefile == TRUE) {
    write.csv(as.data.frame(peakAnno),paste0(outpath,".csv"))
  }
}

#### function to export enriche pathway from peak anno output
enrich_path_dotplot <- function(gr, outpath, writefile, type = NULL, h=4, w=12) { 
  df <- as.data.frame(annotatePeak(gr, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db"))
  if (!is.null(type)) {
    df <- df %>%
      filter(grepl(type, annotation))
  }
  pathway <- enrichPathway(df$geneId)
  print(as.data.frame(pathway))
  pdf(file = paste0(outpath,".pdf"), height = h, width = w)
  print(dotplot(pathway))
  dev.off()
  if (writefile == TRUE) {
    write.csv(as.data.frame(pathway@result),paste0(outpath,".csv"))
  }
}

### function to run deseq###
run_deseq <- function(ctd, cd, des, cutoff, ref, outpath) {
  dds <- DESeqDataSetFromMatrix(countData = ctd,
                                colData = cd,
                                design = formula(paste0("~", des)))
  keep <- rowSums(counts(dds)) >= cutoff
  dds <- dds[keep,]
  dds$type <- relevel(dds$type, ref = ref)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind=FALSE)
  plotPCA(vsd, intgroup=c(des))
  ggsave(outpath)
  res <- results(dds)
  return(res)
}

#### making same n bins from ranges
make_range <- function(chr, start, end) {
  range <- end - start
  jump <- (end - start)/50
  st <- round(c(start, seq(start+jump+1, end, jump)))
  en <- round(seq(start+jump-1, end, jump)+1)
  data.frame(chr, st, en)
}
