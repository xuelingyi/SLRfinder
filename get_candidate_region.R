#### this is the modified version to avoid error when no significant candidate is found #####
library(parallel)
library(cowplot)
library(ggplot2)
library(data.table)
library(SNPRelate)

eucl <- function(i, j){ sqrt((i[,1]-j[1])^2 + (i[,2]-j[2])^2 ) }

get_candidate_regions <- function(data_cls, GT, map, pop, ranks=c("Dext_max_rank","r2_rank","nSNPs_rank","chi2_rank"),nPerm=10000, cores=10, alpha=0.05){

  cat("Generating gds file \n\n")

  name <- paste0("file.gds")
  snpgdsCreateGeno(name, genmat = t(GT),sample.id = 1:nrow(GT), snp.id = map$SNP, snpfirstdim=TRUE)
  file_gds <- snpgdsOpen(name)

  cl <- 1
  cat("Processing data \n\n")
  data_out <-  rbindlist(mclapply(1:nrow(data_cls), function(cl){
    #print(cl)
    cl_info <- data_cls[cl,]
    SNPs <- cl_info$SNPs[[1]]
    pca <- snpgdsPCA(file_gds,snp.id = cl_info$SNPs[[1]],verbose = FALSE)
    gt <- GT[,which(map$SNP %in% SNPs)]
    PC1 <- as.matrix(pca$eigenvect[,1])
    PVE <- pca$eigenval[1]/sum(na.omit(pca$eigenval))
    rm(pca)
    het <- apply(gt, 1, function(x) length(which(x==1))/length(na.omit(x))/2)

    ## polarize so correlation always positive
    if(cor(het,PC1,use = "pair")[1,1]<0) PC1 <- -PC1
    data <- data.table(PC1=as.vector(PC1),Het=het)
    r2 <- summary(lm(PC1~Het,data))$adj.r.squared
    cl_info[,r2:=r2]
    cl_info[,PVE:=PVE]

    ## distance from the extremes
    data[,PC_scaled:=(PC1-min(PC1))/max((PC1-min(PC1)))]
    c1 <- c(PC1=0,het=0) ## bottom left corner
    c2 <- c(PC1=1,het=0.5) ## top right corner

    # get all distances from the corners
    Dists <-  cbind(d_c1=eucl(data[,.(PC_scaled,Het)],c1),d_c2=eucl(data[,.(PC_scaled,Het)],c2))
    # divide the minimum with the maximum possible so values are scaled betweeen 0 and 1
    Dext <- apply(Dists,1,min)/((eucl(matrix(c1,nrow=1),c2))/2)
    #data[,plot(PC_scaled,Het)]

    cl_info[,Dext_mean:=mean(Dext)] # mean
    cl_info[,Dext_max:=max(Dext)] # max

    #data[,plot(PC_scaled,Het)]
    data[,Pop:=pop]
    tbl <- data[,table(Pop,PC_scaled<0.5)]
    #tbl[,] <- abs(as.integer(rnorm(tbl)))*10

    cl_info[,chi2:=sum(apply(tbl,1,function(x){
      #print(x)
      chisq.test(x)$statistic
    }))] # max

    return(cl_info[,.(chr,nSNPs,mean_LD,nE,c,r2,PVE,Dext_mean,Dext_max,chi2,SNPs,data=list(data))])

  },mc.cores=cores))

  data_out[,cluster:=1:nrow(data_out)]
  data_out[order(nSNPs,decreasing = TRUE),nSNPs_rank := 1:nrow(data_out)]
  data_out[order(r2,decreasing = TRUE),r2_rank := 1:nrow(data_out)]
  data_out[order(Dext_max,decreasing = FALSE),Dext_max_rank := 1:nrow(data_out)]
  data_out[order(chi2,decreasing = FALSE),chi2_rank := 1:nrow(data_out)]
  rank=rowSums(data_out[,..ranks])
  data_out[,rank:=rank]
  setorder(data_out,rank)

  #data_out[1:10,]
  #data_out[,cor(r2_rank,Dext_max_rank)]

  cat("Estimating p-values (rank permutation) \n\n")
  exp <- do.call(cbind, mclapply(1:nPerm, function(x){
    sort(rowSums(apply(data_out[,..ranks],2,sample)),decreasing = FALSE)
  },mc.cores = cores))

  data_out[,rank_exp:=as.integer(rowMeans(exp))]

  obs <- data_out[,rank]

  null <- unlist(exp)
  p <- Pval(unlist(exp),obs,alternative = "less")

  ## check and corret for p-value inflation
  cat("Checking and correcting p-value inflation \n\n")
  qq_data <- data.table(exp=-log10(ppoints(p)),obs=-log10(sort(p)))
  lambda <- lm(obs~exp+0,qq_data)$coefficients

  p_gc <- 1/10^(-log10(p)/lambda)

  #table(p.adjust(p_gc,"fdr")<alpha)

  data_out[,p:=p]
  data_out[,p_gc:=p_gc]
  data_out[,p_gc_adj:=p.adjust(p_gc,"fdr")]

  col=rep("steelblue",nrow(data_out))
  col[data_out[,which(p_gc_adj<alpha)]] <- "indianred"


  plot(qq_data, pch=20,col=col,main=paste("QQ-plot | lambda=",round(lambda,2)),xlab="Expected -log10(P)",ylab="Observed -log10(P)")
  abline(0,1)
  abline(lm(obs~exp+0,qq_data),col="salmon")
  qq_plot <- recordPlot()
  #data_out[,plot(r2_rank,Dext_max_rank)]

  cat("Finding candidates \n\n")
  candidates <- data_out[which(p_gc_adj<alpha),]
  if(nrow(candidates) == 0){
    print("no candidates found by the current alpha!")
    PCA_het_data = NULL
    p1 = NULL
  } else {
    PCA_het_data <- rbindlist(apply(candidates,1,function(x){

      region=paste(x$chr, paste(range(as.numeric(do.call(rbind, strsplit(x$SNPs,"_",fixed=TRUE))[,2])),collapse =  "-"), sep=":")

      label = paste(region, paste0("nSNPs=",x$nSNPs),sep=" | ")
      label = paste(label,paste0("p=",signif(x$p_gc_adj,3)),sep=" | ")

      data.table(x$data,rank=x$rank,p=x$p_gc_adj,nSNPs=x$nSNPs,region=region,label=label)

    }))


    p1 <- ggplot(PCA_het_data, aes(PC_scaled,Het)) +
      geom_point() +
      geom_smooth(method = "lm",se=FALSE,col="black")+
      scale_color_viridis_c()+
      theme_bw() +
      facet_wrap(label~.,ncol = 2) +
      theme(aspect.ratio = 1,
            strip.background = element_blank(),
            legend.background = element_blank(),
            strip.text.x = element_text(hjust = 0),
            legend.key = element_blank()) +
      xlab("PC1 (scaled)") +
      ylab("Proportion heterozygous loci")


    p1

    candidates[,region:=apply(candidates,1,function(x)region=paste(x$chr, paste(range(as.numeric(do.call(rbind,strsplit(x$SNPs,"_",fixed=TRUE))[,2])),collapse =  ":"),sep=""))]

    }

  cat("Closing gds file and returning data \n\n")
  snpgdsClose(file_gds)
  system("rm file.gds")
  return(list(data=data_out,candidates=candidates,qq_data=qq_data,plot=p1,qq_plot=qq_plot,plot_data=PCA_het_data))
}
