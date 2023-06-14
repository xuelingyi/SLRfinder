### get the top candidates for each parameter combiantion of each dataset
source("get_candidate_region.R")
source("Pval.R")

mydata="BEL_MAL"



setwd(mydata)
cores=1

## run inside each dataset folder
#

for (min_LD in c(8, 8.5, 9)) {
  for (min.cl.size in c(10, 20)) {
    # move into the parameter folder
    setwd(paste0("LD", min_LD, "cl", min.cl.size))
    
    # load data
    data_cls <- readRDS(file=paste0("data_cls.rds"))
    pos_files <- paste0(mydata, "_v7LG", 1:21, "_a15m75_LD", min_LD, "cl", min.cl.size, ".012.pos")
    GT_files <- paste0(mydata, "_v7LG", 1:21, "_a15m75_LD", min_LD, "cl", min.cl.size, ".012")
    indv <- fread(paste0(mydata, "_v7LG1_a15m75_LD", min_LD, "cl", min.cl.size, ".012.indv"), header = F)
    pop_info <- fread(paste0("../", mydata, ".csv"))
    if(all(indv$V1 == pop_info$SampleID)) {pop <- pop_info$Population 
    }else{
    print("indv and pop do not match!")}

    map <- rbindlist(lapply(1:length(pos_files),function(i){
      pos <- fread(pos_files[i])
      colnames(pos) <- c("Chr","Pos")
      pos$SNP <- paste0(pos$Chr, "_", pos$Pos)
      return(pos)
    }))
    GT <- do.call(cbind, lapply(1:length(pos_files),function(i){ as.matrix(fread(GT_files[i])[,-1]) }))
    GT[GT==-1] <- NA
    save(data_cls, GT, map, pop, file=paste0(mydata, "_LD", min_LD, "cl", min.cl.size, ".RData"))
    
    # load(paste0(mydata, "_LD", min_LD, "cl", min.cl.size, ".RData"))
    ## identify candidate regions
    ## the ranks used for defining candidate regions
    ranks = c("Dext_max_rank","r2_rank","nSNPs_rank","chi2_rank")
    
    cand_regions <- get_candidate_regions(data_cls, GT, map, pop, ranks=ranks, nPerm=10000, cores=cores, alpha=0.05)
    cand_regions$candidates
    cand_regions$qq_plot
    cand_regions$plot
    
    saveRDS(cand_regions, paste0("cand_regions.rds"))

    pairs(cand_regions$data[,.(r2,Dext_max,nSNPs,chi2)])
    
    setwd("../")
  }
}




