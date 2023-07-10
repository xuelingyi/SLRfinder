library(parallel)
library(cowplot)
library(ggplot2)
library(data.table)
library(SNPRelate)

eucl <- function(i, j){ sqrt((i[,1]-j[1])^2 + (i[,2]-j[2])^2 ) }

get_candidate_regions <- function(data_cls, GT, map, pop, ranks=c("Dext_max_rank","R2_rank","nSNPs_rank","chi2_rank"), nPerm=10000, cores=1, alpha=0.05){

  cat("Generating gds file \n")
  name <- paste0("file.gds")
  snpgdsCreateGeno(name, genmat = t(GT),sample.id = 1:nrow(GT), snp.id = map$SNP, snpfirstdim=TRUE)
  file_gds <- snpgdsOpen(name)

  cat("Processing data \n")
  data_out <-  rbindlist(mclapply(1:nrow(data_cls), function(cl){
    cl_info <- data_cls[cl,]
    SNPs <- cl_info$SNPs[[1]]
    pca <- snpgdsPCA(file_gds,snp.id = cl_info$SNPs[[1]],verbose = FALSE)
    gt <- GT[,which(map$SNP %in% SNPs)]
    PC1 <- as.matrix(pca$eigenvect[,1])
    PVE <- pca$eigenval[1]/sum(na.omit(pca$eigenval))
    PC2 <- as.matrix(pca$eigenvect[,2])
    PVE2 <- pca$eigenval[2]/sum(na.omit(pca$eigenval))
    rm(pca)
    het <- apply(gt, 1, function(x) length(which(x==1))/length(na.omit(x)))

    ## polarize so correlation always positive
    if(cor(het,PC1,use = "pair")[1,1]<0) PC1 <- -PC1
    data <- data.table(PC1=as.vector(PC1),Het=het)
    data$PC2 = PC2
    R2 <- summary(lm(PC1~Het,data))$adj.r.squared
    cl_info[,R2:=R2]
    cl_info[,PVE:=PVE]

    ## re-scale PC1 wihtin [0,1]
    data[,PC_scaled:=(PC1-min(PC1))/max((PC1-min(PC1)))]

    # get all distances from the corners
    c1 <- c(PC1=0,het=0) ## bottom left corner
    c2 <- c(PC1=1,het=1) ## top right corner
    Dists <-  cbind(d_c1=eucl(data[,.(PC_scaled,Het)],c1),d_c2=eucl(data[,.(PC_scaled,Het)],c2))
    # divide the minimum with the maximum possible so values are scaled betweeen 0 and 1
    Dext <- apply(Dists,1,min)/((eucl(matrix(c1,nrow=1),c2))/2)

    cl_info[,Dext_mean:=mean(Dext)] # mean
    cl_info[,Dext_max:=max(Dext)] # max

    data[,Ind:=ind]
    data[,Pop:=pop]
    tbl <- data[,table(Pop,PC_scaled<0.5)]
    cl_info[,chi2:=sum(apply(tbl,1,function(x){
    chisq.test(x)$statistic
     }))] # max
    ## NOTE: this step might give the warning "In chisq.test(x) : Chi-squared approximation may be incorrect" if some of the populations have too few individuals (e.g., below 10 in NOR-TYR). However, this warning should not impact the final result because 1) when the individual distribution is uneven the chi-square is still significant, and 2) the final chi2 value is the sum of the chi-square statistic across all populations and thus should not be biased by one or a few tests. 

    return(cl_info[, .(chr, nSNPs, mean_LD, nE, c, R2, PVE, PVE2, Dext_mean, Dext_max, chi2, SNPs, data=list(data))])

  },mc.cores=cores))

  data_out[,cluster:=1:nrow(data_out)]
  data_out[order(nSNPs,decreasing = TRUE),nSNPs_rank := 1:nrow(data_out)]
  data_out[order(R2,decreasing = TRUE),R2_rank := 1:nrow(data_out)]
  data_out[order(Dext_max,decreasing = FALSE),Dext_max_rank := 1:nrow(data_out)]
  data_out[order(chi2,decreasing = FALSE),chi2_rank := 1:nrow(data_out)]
  rank=rowSums(data_out[,..ranks])
  data_out[,rank:=rank]
  setorder(data_out,rank)

  cat("Estimating p-values (rank permutation) \n\n")
  exp <- do.call(cbind, mclapply(1:nPerm, function(x){
    sort(rowSums(apply(data_out[,..ranks],2,sample)),decreasing = FALSE)
  },mc.cores = cores))

  data_out[,rank_exp:=as.integer(rowMeans(exp))]
  obs <- data_out[,rank]
  null <- unlist(exp)
  p <- unlist(mclapply(obs,function(x) (sum(null<x)+1)/(length(null)+1),mc.cores = cores))

  ## check and corret for p-value inflation
  cat("Checking and correcting p-value inflation \n\n")
  qq_data <- data.table(exp=-log10(ppoints(p)),obs=-log10(sort(p)))
  lambda <- lm(obs~exp+0,qq_data)$coefficients
  p_gc <- 1/10^(-log10(p)/lambda)

  data_out[,p:=p]
  data_out[,p_gc:=p_gc]
  data_out[,p_gc_adj:=p.adjust(p_gc,"fdr")]

  cat("Finding candidates \n\n")
  candidates <- data_out[which(p_gc_adj<alpha),]
  if(nrow(candidates) == 0){
    print("No candidates found! Report the top 5 ranked LD clusters.")
    PCA_het_data = data_out[1:5, ]
  } else {
    PCA_het_data <- rbindlist(apply(candidates,1,function(x){
      region=paste(x$chr, paste(range(as.numeric(do.call(rbind, strsplit(x$SNPs,"_",fixed=TRUE))[,2])),collapse =  "-"), sep=":")
      label = paste(region, paste0("nSNPs=",x$nSNPs),sep=" | ")
      label = paste(label,paste0("p=",signif(x$p_gc_adj,3)),sep=" | ")
      data.table(x$data,rank=x$rank,p=x$p_gc_adj,nSNPs=x$nSNPs,region=region,label=label)
    }))

  candidates[,region:=apply(candidates,1,function(x)region=paste(x$chr, paste(range(as.numeric(do.call(rbind,strsplit(x$SNPs,"_",fixed=TRUE))[,2])),collapse=":"),sep=":"))]
  }

  cat("Closing gds file and returning data \n\n")
  snpgdsClose(file_gds)
  system("rm file.gds")
  return(list(data_out=data_out, candidates=candidates, qq_data=qq_data, PCA_het_data=PCA_het_data))
}
