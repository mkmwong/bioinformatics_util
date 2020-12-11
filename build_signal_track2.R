#### This is partial code from Devin King's pending TRAD-Seq computational ####
#### package. Devin devking@stanford.edu #
build_signal_track2 <- function(bam_in, bed_out, tmp_out, mapQ=0, ranges_operation=function(x) promoters(x,upstream=2,downstream=0), include_qname=TRUE, include_seqs='reverse_complement', genome=Hsapiens, filter_by_seqs=c('TT','TC','CT','CC'), overwrite=FALSE, verbose=TRUE) {
  dat <- mapply(.build_signal_track2, bam_in, bed_out, tmp_out, MoreArgs=list(mapQ=mapQ, ranges_operation=ranges_operation, include_seqs=include_seqs, genome=genome, filter_by_seqs=filter_by_seqs, overwrite=overwrite, verbose=verbose))
}  

.build_signal_track2 <- function(input, output, tmp_out, mapQ=0, ranges_operation=function(x) promoters(x,upstream=2,downstream=0), include_seqs='reverse_complement', genome=NULL, filter_by_seqs=c('TT','TC','CT','CC'), include_qname=TRUE, overwrite=FALSE, verbose=TRUE) {
  require('R.utils')
  require('GenomicAlignments')
  require('rtracklayer')  
  library(BSgenome.Hsapiens.UCSC.hg19) 
#  if(include_qname) stop('currently not supported')
  ### Check arguments
  if(!overwrite) if(any(file.exists(output))) stop('One or more output files already exist and will not be overwritten.')
  if(overwrite) if(file.exists(output)) unlink(output)
  if(isTRUE(include_seqs) | as.character(include_seqs) == 'reverse_complement') if(is.null(genome)) stop('A "genome" must be provided if using "include_seqs"')
  
  ### Define import parameters
  sbf <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  if(include_qname) sbp <- ScanBamParam( flag=sbf, what='qname') else {
    sbp <- ScanBamParam(flag=sbf)
  }
  ### Open connection to input BAM file
  #IN <- BamFile(input, yieldSize=chunk_size, asMates=TRUE)
  #open(IN)
  #on.exit(.close_connection(IN))
  OUT <- tempfile(tmpdir=dirname(output))
  #OUT <- tempfile(tmpdir=getwd())
  
  ### Process
  total_reads <- 0
  before_filt <- table(NULL)
  after_filt <- table(NULL)
  
  if(verbose) message('Building signal track for BAM file ',input,' ...')
  ga <- readGAlignmentPairs(input, param=sbp); gc()
  total_reads <- total_reads + length(ga)
  print(head(ga))
  if(verbose) message('Imported ',total_reads,' reads')
  # perform ranges operation
  if(verbose) message('Performing ranges operation ...')
  if(include_qname) {
      mcols(ga) <- mcols(GenomicAlignments::first(ga))
      ga <- granges(ga,on.discordant.seqnames='drop',use.mcols=TRUE)
    } else {
    ga <- granges(ga,on.discordant.seqnames='drop'); gc()
  }
  ga <- ga[strand(ga)!='*']; gc()
  ga <- ranges_operation(ga); gc()
  # Note, if the ranges_operation creates ranges that do not reside on the genome, those ranges must be identified and removed
  ga <- GenomicRanges::trim(ga); gc()
  rmd <- .mode(width(ga))
  ga <- ga[width(ga)==rmd]; gc()
  # get seqs
  if(verbose) message('Annotating ranges with seqs ...')
  ga <- .annotate_ranges_with_seqs(ga, include_seqs=include_seqs, genome=genome, filter_by_seqs=filter_by_seqs); gc()
  stats <- ga$tab
  ga <- ga$x;gc()
  # write to file
  #ga <- GenomicRanges::sort(ga)
  print(mcols(ga))
  if(include_qname) names(ga) <- mcols(ga)$qname
  print(head(as.data.frame(ga)))
  write.csv(as.data.frame(ga), tmp_out, sep='\t',quote=FALSE)
  export.bed(object=ga, con=OUT)
  # update stats
  before_filt <- .add_tables(before_filt,stats$before_filt)
  after_filt <- .add_tables(after_filt,stats$after_filt)
  
  ### gzip output
  if(.is_compressed(output)) {
    if(verbose) message('Compressing output file ...')
    gzip(OUT)
    file.rename(paste0(OUT,'.gz'),output)
  } else {
    file.rename(OUT,output)
  }
  
  ### Output stats
  total_reads <- setNames(total_reads,'total')
  stats <- .ba(total_reads=total_reads,imported_reads=before_filt,final_reads=after_filt)
  stats['imported_reads','total'] <- sum(stats['imported_reads',])
  stats['final_reads','total'] <- sum(stats['final_reads',])
  stats_file <- append_to_filename(append_to_filename(output,insert='_stats'),'txt',as_ext = T)
  if(verbose) message('Metrics are available in ',stats_file)
  write.table(stats,file=stats_file,sep='\t',quote=FALSE)
  
}
library(BSgenome.Hsapiens.UCSC.hg19)
source("utilities.R")

#build_signal_track2("../test.bam","../test.bed",tmp_out="../test.tab")
#build_signal_track2("../TCATGGCA_rmdup.bam","../TCATGGCA.bed",tmp_out="../TCATGGCA.tab")
#build_signal_track2("../CACAGTGT_rmdup.bam","../CACAGTGT.bed",tmp_out="../CACAGTGT.tab")
#build_signal_track2("../TCCAGACT_rmdup.bam","../TCCAGACT.bed",tmp_out="../TCCAGACT.tab")
#build_signal_track2("../TTCTGTGG_rmdup.bam","../TTCTGTGG.bed",tmp_out="../TTCTGTGG.tab")



