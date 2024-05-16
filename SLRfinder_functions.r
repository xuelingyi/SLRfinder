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

  out = NULL
  if(nrow(white.list) > 1){
  #Parse the edge data to create a graph object
  g <- graph_from_edgelist(apply(white.list[,c("from", "to")], 2, function(o) as.character(o)), directed = FALSE)
  #Decompose the graph object into a list of components (i.e., LD clusters) and only keep the components containing at least min.cl.size number of vertices (i.e., SNPs)
  d_g <- decompose(g, min.vertices = min.cl.size)
  #For each component (i.e., LD cluster), create a vertex sequence (vs) containing all graph vertices (i.e., SNPs named by positions on this chr)
  d_g <- lapply(d_g, function(cl) V(cl)$name)
  
  #Summerize data for each LD cluster
  out <- rbindlist(mclapply(d_g, function(dg){
    #number of SNPs (nodes)
    nSNPs <- length(dg)
    #mean LD among SNPs in the LD cluster (edge weight)
    mean_LD <- mean(as.numeric(white.list[white.list$from %in% dg & white.list$to %in% dg, "r2"]))
    #number of SNP pairs having high LD in this cluseter (retained edges)
    nE <- nrow(white.list[white.list$from %in% dg & white.list$to %in% dg, ])
    #retained number of edges/nodes as an approximate of clustering coefficiency (higher c means tighter clustering)
    c <- nE/nSNPs
    data.table(chr, nSNPs, mean_LD, nE, c, SNPs=list(dg))
  },mc.cores=ncores))
  }
  return(out)  
}


eucl <- function(i, j){ sqrt((i[,1]-j[1])^2 + (i[,2]-j[2])^2 ) }


## process data
get_data_output = function(data_cls, GT, map, pop, sex_info=T, heterog_homog = c(0.5, 0.5), cores=1){
  
  cat("Generating gds file \n")
  name <- paste0("file.gds")
  snpgdsCreateGeno(name, genmat = t(GT),sample.id = 1:nrow(GT), snp.id = map$SNP, snpfirstdim=TRUE)
  file_gds <- snpgdsOpen(name)
  
  cat("Processing data \n")
  data_out <- rbindlist(mclapply(1:nrow(data_cls), function(cl){
    cl_info <- data_cls[cl,]
    SNPs <- cl_info$SNPs[[1]]
    gt <- GT[,which(map$SNP %in% SNPs)]
    
    #if(nrow(gt) > 0){}
    pca <- snpgdsPCA(file_gds, snp.id = cl_info$SNPs[[1]], verbose = FALSE)
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
    
    ## Note: individuals having only missing data would have het == NaN, which would be discarded in R2 estimations
    R2 <- summary(lm(PC1~Het,data))$adj.r.squared
    cl_info[,R2:=R2]
    cl_info[,PVE:=PVE]
    
    ## re-scale PC1 wihtin [0,1]
    data[,PC_scaled:=(PC1-min(PC1))/max((PC1-min(PC1)))]
    
    # get all distances from the corners
    ## bottom left corner
    c1 <- c(min(data[data$Het == min(na.omit(data$Het)), "PC_scaled"]), min(na.omit(data$Het))) 
    ## top right corner
    c2 <- c(max(data[data$Het == max(na.omit(data$Het)), "PC_scaled"]), max(na.omit(data$Het)))  
    d_c1=eucl(data[,.(PC_scaled,Het)],c1)
    d_c2=eucl(data[,.(PC_scaled,Het)],c2)
    Dists <-  cbind(d_c1, d_c2)
    #standarize distances for between LD comparisons 
    Dext <- apply(Dists,1,min)/((eucl(matrix(c1,nrow=1),c2))/2)
    cl_info[,Dext_mean:=mean(Dext[!is.na(Dext)])] # mean
    cl_info[,Dext_max:=max(Dext[!is.na(Dext)])] # max
    cl_info[,Dext_var:=var(Dext[!is.na(Dext)])]
    
    data[,Ind:=ind]
    data[,Pop:=pop]
    tbl <- data[,table(Pop,PC_scaled<0.5)]
    ## F, T
    if(all(heterog_homog == 0.5)){
      cl_info[,chi2:=sum(apply(tbl,1,function(x){
        chisq.test(x)$statistic
      }))] # max
    } else {
      cl_info[,chi2:=sum(apply(tbl,1,function(x){
        chisq.test(x, p=heterog_homog)$statistic
      }))]
    }
    
    cl_info[,Sex_g:=1]
    if(sex_info){
      data = merge(data, sif[, c("SampleID", "sex")], by.x="Ind", by.y="SampleID", sort=F)
      cl_info[,Sex_g:=(
        (min(sum(data[which(d_c1 < d_c2), "sex"][[1]] %in% c("female", "Female", "F")), sum(data[which(d_c1 < d_c2), "sex"][[1]] %in% c("male", "Male", "M"))) + 
          min(sum(data[which(d_c1 >= d_c2), "sex"][[1]] %in% c("female", "Female", "F")), sum(data[which(d_c1 >= d_c2), "sex"][[1]] %in% c("male", "Male", "M"))))/nrow(data[data$sex %in% c("female", "Female", "F", "male", "Male", "M"),])
        )]}
    
    return(as.data.frame(cl_info[, .(chr, nSNPs, mean_LD, nE, c, R2, PVE, PVE2, Dext_mean, Dext_max, Dext_var, Sex_g, chi2, SNPs, data=list(data))]))
  },mc.cores=cores))
  
  cat("Closing gds file and returning data \n\n")
  snpgdsClose(file_gds)
  system("rm file.gds")
  
  # cluster ID
  data_out[,cluster:=1:nrow(data_out)]
  
  # rank clusters: use the same rank when having the same values
  data_out$Dext_mean_rank=rank(data_out$Dext_mean, ties.method="min")
  data_out$Dext_max_rank=rank(data_out$Dext_max, ties.method="min")
  data_out$Dext_var_rank=rank(data_out$Dext_var, ties.method="min")
  data_out$chi2_rank=rank(data_out$chi2, ties.method="min")
  # rank by decreasing order!
  data_out$nSNPs_rank=rank(-data_out$nSNPs, ties.method="min")
  data_out$R2_rank=rank(-data_out$R2, ties.method="min")
  return(data_out)
}

## get candidate regions
get_candidate_regions <- function(data_out, ranks=c("Dext_var_rank", "R2_rank","nSNPs_rank","chi2_rank"), nPerm=10000, cores=1, alpha=0.05){
  
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
    if(nrow(data_out) < 5){ candidates <- data_out  }
    if(nrow(data_out) >= 5){ candidates <- data_out[1:5,] } 
  } 
  
  myindex = length(grep("_", candidates[1, "chr"])) + 2
  
  PCA_het_data <- rbindlist(apply(candidates,1,function(x){
    region=paste(x$chr, paste(range(as.numeric(do.call(rbind, strsplit(x$SNPs,"_",fixed=TRUE))[,myindex])),collapse =  "-"), sep=":")
    label = paste(region, paste0("nSNPs=",x$nSNPs),sep=" | ")
    label = paste(label,paste0("p=",signif(x$p_gc_adj,3)),sep=" | ")
    data.table(x$data,rank=x$rank,p=x$p_gc_adj,nSNPs=x$nSNPs,region=region,label=label)
  }))
  
  candidates[,region:=apply(candidates,1,function(x)region=paste(x$chr, paste(range(as.numeric(do.call(rbind,strsplit(x$SNPs,"_",fixed=TRUE))[,myindex])),collapse="-"),sep=":"))]
  
  print(paste0(nrow(candidates), " candidate(s)"))
  candidates = merge(candidates, LG, by="chr")
  
  write.csv(candidates[, c("chr", "lg", "region", "Sex_g", "rank", "p_gc_adj",
                           "nSNPs", "R2", "chi2", "Dext_var", "Dext_mean", "Dext_max", 
                           "nSNPs_rank", "R2_rank", "chi2_rank", "Dext_var_rank", "Dext_mean_rank", "Dext_max_rank",
                           "mean_LD")], 
            "candidates.csv", row.names = F)
  
  return(list(data_out=data_out, candidates=candidates, qq_data=qq_data, PCA_het_data=PCA_het_data, lambda=lambda))
}

