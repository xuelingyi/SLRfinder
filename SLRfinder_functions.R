library(igraph)
library(data.table)

library(parallel)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(data.table)
library(SNPRelate)

get_single_LD_cluster <- function(geno.LD, min_LD = 0.85, min.cl.size = 20){
  
  chr = unique(geno.LD$CHR) 
  #Remove edges lower than min_LD
  white.list <- subset(geno.LD, as.numeric(r2) > min_LD)
  
  #Parse the edge data to create a graph object
  g <- graph_from_edgelist(apply(white.list[,c("from", "to")], 2, function(o) as.character(o)), directed = FALSE)
  #Decompose the graph object into a list of components (i.e., LD clusters) and only keep the components containing at least min.cl.size number of vertices (i.e., SNPs)
  d_g <- decompose(g, min.vertices = min.cl.size)
  #For each component (i.e., LD cluster), create a vertex sequence (vs) containing all graph vertices (i.e., SNPs named by positions on this chr)
  d_g <- lapply(d_g, function(cl) V(cl)$name)
  
  #Summerize data for each LD cluster
  out <- rbindlist(lapply(d_g, function(dg){
    #number of SNPs (nodes)
    nSNPs <- length(dg)
    #mean LD among SNPs in the LD cluster (edge weight)
    mean_LD <- mean(as.numeric(white.list[white.list$from %in% dg & white.list$to %in% dg, "r2"]))
    #number of SNP pairs having high LD in this cluseter (retained edges)
    nE <- nrow(white.list[white.list$from %in% dg & white.list$to %in% dg, ])
    #retained number of edges/nodes as an approximate of clustering coefficiency (higher c means tighter clustering)
    c <- nE/nSNPs
    data.table(chr, nSNPs, mean_LD, nE, c, SNPs=list(dg))
  }))
  
  return(out)  
}


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
    my.cor = cor(het,PC1,use = "pair")[1,1]
    if(!(is.na(my.cor)) & my.cor<0) PC1 <- -PC1   
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
  
  qq_data$col=rep("steelblue", nrow(data_out))
  qq_data$col[which(data_out$p_gc_adj<alpha)] <- "indianred"
  
  cat("Finding candidates \n\n")
  candidates <- data_out[which(p_gc_adj<alpha),]
  if(nrow(candidates) == 0){
    print("No candidates found! Report the top 5 ranked LD clusters.")
    candidates <- data_out[1:5,]    
  } 
  
  myindex = length(grep("_", candidates[1, "chr"])) + 2
  
  PCA_het_data <- rbindlist(apply(candidates,1,function(x){
    region=paste(x$chr, paste(range(as.numeric(do.call(rbind, strsplit(x$SNPs,"_",fixed=TRUE))[,myindex])),collapse =  "-"), sep=":")
    label = paste(region, paste0("nSNPs=",x$nSNPs),sep=" | ")
    label = paste(label,paste0("p=",signif(x$p_gc_adj,3)),sep=" | ")
    data.table(x$data,rank=x$rank,p=x$p_gc_adj,nSNPs=x$nSNPs,region=region,label=label)
  }))
  
  candidates[,region:=apply(candidates,1,function(x)region=paste(x$chr, paste(range(as.numeric(do.call(rbind,strsplit(x$SNPs,"_",fixed=TRUE))[,myindex])),collapse="-"),sep=":"))]
  
  cat("Closing gds file and returning data \n\n")
  snpgdsClose(file_gds)
  system("rm file.gds")
  return(list(data_out=data_out, candidates=candidates, qq_data=qq_data, PCA_het_data=PCA_het_data, lambda=lambda))
}
