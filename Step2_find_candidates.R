### get the top candidates for each parameter combiantion of each dataset
mydata="BEL_MAL"

source("get_candidate_region.R")
source("Pval.R")

## run inside each dataset folder
setwd(mydata)
cores=1
## the ranks used for defining candidate regions
ranks = c("Dext_max_rank","r2_rank","nSNPs_rank","chi2_rank")

# run below within each parameter folder
for (min_LD in c(8, 8.5, 9)) {
  for (min.cl.size in c(10, 20)) {
    setwd(paste0("LD", min_LD, "cl", min.cl.size))
    
    ### alternatively, save this script as candidate.R and run it in a unix loop:
    #cd LD${ld}cl${cl}
    #echo "min_LD="${ld} >> param
    #echo "min.cl.size="${cl} >> param
    #cat param ../candidate.R > candidate_LD${ld}cl${cl}.R
    #srun singularity_wrapper exec Rscript --no-save candidate_LD${ld}cl${cl}.R
    #cd ../
    
    load(paste0("LD", min_LD, "cl", min.cl.size, ".RData"))
   
    cand_regions <- get_candidate_regions(data_cls, GT, map, pop, ranks=ranks, nPerm=10000, cores=cores, alpha=0.05)
    pdf(paste0("../", mydata, "_LD", min_LD, "cl", min.cl.size, "_Rplots.pdf"))
    print(cand_regions$qq_plot)
    print(cand_regions$plot)
    pairs(cand_regions$data[,.(r2,Dext_max,nSNPs,chi2)])
    dev.off()

    saveRDS(cand_regions, paste0(mydata, "_LD", min_LD, "cl", min.cl.size, "_cand_regions.rds"))
    setwd("../")
  }
}
