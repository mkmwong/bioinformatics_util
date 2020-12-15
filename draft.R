#### install_pkg                                                           ####
#### This is a supporting function to check for other functions if all     ####
#### their dependencies has been installed. If not, install them.          ####
#### parameter: pkg; list of packages required for the function            ####
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
#### parameter: input; could be either a data.frame object or a            ####
####                   SimpleRleList object.                               ####
####            value: the chromosome to be removed.                       ####
####            colum: the column where the chromosome info is stored.     ####
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
  return( tmp )
}

### Function to bin bigwig file into a 
binning_GR <- function(binsGR, datrle, outpath) {

  bio_pkg <- c("rtracklayer","GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19")
  install_pkg(bio_pkg,"bioconductor")
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  si.gen <- seqinfo(gen)
  si <- si.gen[names(dat)]
  ba <- binnedAverage(binsGR,numvar=datrle,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba),outpath,sep=""))
  return(ba)
}

### testing usage ###
bins <- read.table("~/Desktop/tmp/Old_Stuffs/final_figures/data/hg19gap.txt", sep="\t")
bins <- remove_chrom(bins, "chrM", "V2")
binsGR = GRanges(seqnames = bins$V1, ranges = IRanges(bins$V2+1, bins$V3), stat = bins$V4)
dat <- import.bw(con="~/Desktop/tmp/Old_Stuffs/final_figures/data/unbinned_bigwig/wgEncodeFsuRepliChipImr90WaveSignalRep1.bigWig", as='RleList')
dat[dat == 0] <- NA
dat <- remove_chrom(dat, "chrM")
a <- binning_GR(binsGR, datrle, "~/Desktop/something.csv")
