### run this script within each parameter folder of the dataset folder
source("../../candidate_region_functions.R")

### load data
#LD clusters
data_cls <- readRDS(file=paste0("data_cls.rds"))
#SNPs
files <- paste0("./file012/", list.files("file012"))
indv_files <- files[grep(".indv",files)]
pos_files <- files[grep(".pos",files)]
GT_files <- files[!grepl(".log",files) & !grepl(".indv",files) & !grepl(".pos",files)]

map <- rbindlist(lapply(1:length(pos_files),function(i){
  pos <- fread(pos_files[i])
  colnames(pos) <- c("Chr","Pos")
  pos$SNP <- paste0(pos$Chr, "_", pos$Pos)
  return(pos)
}))
GT <- do.call(cbind, lapply(1:length(pos_files),function(i){ as.matrix(fread(GT_files[i])[,-1]) }))
#recode missing data
GT[GT==-1] <- NA

#pop and sif info
sif = list.files("../")
pop_info <- fread(paste0("../", sif[grep(".csv", sif)]))
indv <- fread(indv_files[1], header=F)
pop_info = pop_info[factor(pop_info$SampleID, levels = ind),]
if(all(indv$V1 == pop_info$SampleID)) {
  pop <- pop_info$Population
  ind <- pop_info$SampleID
}else{
  print("indv and pop do not match!")}

save(data_cls, GT, map, ind, pop, file="GT.RData")

### get the top candidate LD clusters
## the ranks used for defining candidate regions
ranks = c("Dext_max_rank", "R2_rank", "nSNPs_rank", "chi2_rank")
cand_regions <- get_candidate_regions(data_cls, GT, map, pop, ranks=ranks, nPerm=10000, cores=1, alpha=0.05)
saveRDS(cand_regions, "cand_regions.rds")

