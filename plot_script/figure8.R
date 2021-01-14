setwd("~/Desktop/yr5/summ/bioinformatics_util/data/")
source("../draft.R")
pkg <- c("boot", "ggpubr")
install_pkg(pkg, "common")

merge_TADs <- function(file1, file2, norm = FALSE) {
  tmp <- fread(file1) %>%
    mutate(order = rep(1:50, length(id)/50)) %>% 
    select(c(order, id), average) %>% 
    spread(order, average)
  
  tmp_c <- fread(file2) %>%
    mutate(order = rep(1:50, length(id)/50)) %>% 
    select(c(order, id), average) %>% 
    spread(order, average)
  
  diff_len <- length(tmp_c$id) - length(tmp$id)
  
  if (diff_len > 0) {
    supp <- as.data.frame(matrix(nrow = diff_len, ncol = ncol(tmp)))
    tmp <- rbind(tmp, supp, use.names = FALSE)
  }
  else {
    supp <- as.data.frame(matrix(nrow = -diff_len, ncol = ncol(tmp)))
    tmp <- rbind(tmp_c, supp, use.names = FALSE)
  }
  
  summ <- cbind(tmp[,27:51], tmp_c[,2:51], tmp[,2:26])
  if(norm == TRUE){
    summ <- summ/median(na.omit(as.matrix(summ)))
  }
  return(summ)
}

sampleMean <- function(x,d) {
  return(mean(na.omit(x[d])))
}

bootstrap_for_plot <- function(summ) {
  int <- apply(summ,2,function(y){ 
    b<-boot(y,sampleMean,R=1000)
    c(mean(b$t),boot.ci(b,type="perc", conf=0.95)$percent[4:5])
  })
  intt <- as.data.frame(t(int))
  intt$pos <- 1:100
  return(intt)
}

save_plot <- function(f1, f2, f3, f4, norm, outpath, ylabel="", ylim1=c(-1,1), ylim2=c(-1,1)) {
  
  summ <- merge_TADs(f1, f2, norm)
  intt <- bootstrap_for_plot(summ)
  a <- ggplot(intt,aes(y=log2(V1),x=pos,group=1)) + 
    annotate("rect", xmin = 25,xmax = 75,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) + 
    #geom_rect(aes(xmin = 26,xmax = 76,ymin = -Inf, ymax = Inf),fill = "indianred", alpha = .35) +
    geom_ribbon(data=intt,aes(ymin=log2(V2),ymax=log2(V3),x=pos),alpha=0.3) +
    geom_line() + theme_bw() + xlab("")+ ylab(paste(ylabel,"(log2(RbKO))")) + 
    guides(fill=FALSE) + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.title.y = element_text(face="bold"))
  
  summ1 <- merge_TADs(f3,f4,norm)
  intt1 <- bootstrap_for_plot(summ1)
  b <- ggplot(intt1,aes(y=log2(V1),x=pos,group=1)) + 
    annotate("rect", xmin = 26,xmax = 76,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) + 
    geom_ribbon(data=intt1,aes(ymin=log2(V2),ymax=log2(V3),x=pos),alpha=0.3) +
    geom_line() + theme_bw() + xlab("")+ ylab(paste(ylabel,"(log2(WT))")) +
    guides(fill=FALSE) + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.title.y = element_text(face="bold"))
  
  summ2 <- summ/summ1
  intt2 <- bootstrap_for_plot(summ2)
  c <- ggplot(intt2,aes(y=log2(V1),x=pos,group=1)) + 
    annotate("rect", xmin = 26,xmax = 76,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) + 
    geom_ribbon(data=intt2,aes(ymin=log2(V2),ymax=log2(V3),x=pos),alpha=0.3) +
    geom_line() + theme_bw() + xlab("")+ ylab(paste(ylabel,"(log2(RbKO/WT))")) + 
    guides(fill=FALSE) + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.title.y = element_text(face="bold"))
  
  a = a+ylim(ylim1)
  b = b+ylim(ylim1)
  c = c+ylim(ylim2)
  ggsave(outpath,ggarrange(c,a,b, ncol=3, nrow=1, align="hv"), device = "pdf", height = 3, width = 6)
  
  
}

save_plot("TAD_data/RbKOCPD_openTAD.csv", "TAD_data/RbKOCPD_closedTAD.csv",
          "TAD_data/WTCPD_openTAD.csv", "TAD_data/WTCPD_closedTAD.csv", TRUE, 
          "../plot/figure8_CPD_normalized.pdf",  ylab="CPD", ylim1=c(-0.1,0.2), ylim2=c(-0.05,0.1))

save_plot("TAD_data/RbKOH3K9me3_openTAD.csv", "TAD_data/RbKOH3K9me3_closedTAD.csv",
          "TAD_data/WTH3K9me3_openTAD.csv", "TAD_data/WTH3K9me3_closedTAD.csv", TRUE, 
          "../plot/figure8_H3K9me3_normalized.pdf",  ylab="H3K9me3", ylim1=c(-0.25,0.35), ylim2=c(-0.1,0.1))

save_plot("TAD_data/RbKOH3K27me3_openTAD.csv", "TAD_data/RbKOH3K27me3_closedTAD.csv",
          "TAD_data/WTH3K27me3_openTAD.csv", "TAD_data/WTH3K27me3_closedTAD.csv", TRUE, 
          "../plot/figure8_H3K27me3_normalized.pdf",  ylab="H3K27me3", ylim1=c(-0.25,0.8), ylim2=c(-.3,.2))

save_plot("TAD_data/RbKOCPD_openTAD.csv", "TAD_data/RbKOCPD_closedTAD.csv",
          "TAD_data/WTCPD_openTAD.csv", "TAD_data/WTCPD_closedTAD.csv", FALSE, 
          "../plot/figure8_CPD.pdf",  ylab="CPD", ylim1=c(-0.45,-0.1), ylim2=c(-0.05,0.1))

save_plot("TAD_data/RbKOH3K9me3_openTAD.csv", "TAD_data/RbKOH3K9me3_closedTAD.csv",
          "TAD_data/WTH3K9me3_openTAD.csv", "TAD_data/WTH3K9me3_closedTAD.csv", FALSE, 
          "../plot/figure8_H3K9me3.pdf",  ylab="H3K9me3", ylim1=c(-0.7,-0.05), ylim2=c(-0.1,0.1))

save_plot("TAD_data/RbKOH3K27me3_openTAD.csv", "TAD_data/RbKOH3K27me3_closedTAD.csv",
          "TAD_data/WTH3K27me3_openTAD.csv", "TAD_data/WTH3K27me3_closedTAD.csv", FALSE, 
          "../plot/figure8_H3K27me3.pdf",  ylab="H3K27me3", ylim1=c(-1,0.2), ylim2=c(-.3,.2))

