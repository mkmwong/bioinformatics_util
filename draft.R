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
  si <- si.gen[names(dat)]
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

### pca_plot
#### This is a function to perform PCA on a named matrix, and returning    ####
#### principal component loading vectors for graphing.                     ####
#### parameters: df: A dataframe object to preform principal component     ####
####                 analysis on.                                          ####
####             outpath; A string of the name of the output directory,    ####
####                      including the prefix of the new file.            ####
#### return: A dataframe containing PC loading vectors.                    ####
pca_plot <- function(df, outpath ) {
  all_pca <- prcomp(df, center = TRUE, scale. = TRUE)
  all_out <- as.data.frame(all.pca$rotation)
  all_out$feature <- row.names(all_out)
  save_files(all_out,csv, outpath)
  return(all_out)
}
#b <- na.omit(a)
#b <- data.matrix(b[,4:40])
#c <- prcomp(b, center=TRUE, scale.=TRUE)
