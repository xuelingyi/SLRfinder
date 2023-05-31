## This function takes the edge list of LD estiamtes on SNP pairs and outputs a white list of SNPs that are in LD above 0.9 (min_LD = 0.9) and in LD clusters larger than 10 SNPs (min.cl.size=10). 

# NOTE: LD.data is the full name (including path) of the *.geno.ld file putput from vcftools 
# NOTE: the function works on each chr separately 

get_single_LD_cluster <- function(LD.data, min_LD = 0.9, min.cl.size=10){
  
  geno.LD<-read.table(LD.data, header = T)
  names(geno.LD) = c("CHR",	"from",	"to",	"N_INDV",	"r2")
  chr = unique(geno.LD$CHR)
  
  cat("Removing edges data lower than min_LD \n")
  white.list <- subset(geno.LD, as.numeric(r2)>min_LD)
  
  library("igraph")
  cat("Parsing to graph object \n")
  g <- graph_from_edgelist(apply(white.list[,c("from", "to")], 2, function(o) as.character(o)), directed = FALSE)
  E(g)$weight <- as.numeric(white.list[,"r2"])
  cat("Decomposing graph \n")
  d_g <- decompose.graph(g)
  d_g <- d_g[sapply(d_g, length)>min.cl.size]
  d_g <- lapply(d_g, function(cl) V(cl)$name)
  #length(d_g)
  
  cat("Summerizing data \n")
  library(data.table)
  out <- rbindlist(lapply(d_g, function(g){
    nSNPs <- length(g)
    mean_LD <- mean(as.numeric(white.list[white.list$from %in% g & white.list$to %in% g,"r2"]))
    nE <- nrow(white.list[white.list$from %in% g & white.list$to %in% g, ])
    nE/nSNPs
    data.table(chr, nSNPs, mean_LD,nE,c=nE/nSNPs,SNPs=list(g))
  }))
  cat("Done \n")
  
  return(out)  
}
