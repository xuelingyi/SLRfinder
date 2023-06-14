library(data.table)
### run this R script within the LDcl folder within each dataset
    # load data
    data_cls <- readRDS(file=paste0("data_cls.rds"))

    files <- paste0("./file012/", list.files("file012"))
    indv_files <- files[grep(".indv",files)]
    pos_files <- files[grep(".pos",files)]
    GT_files <- files[!grepl(".indv",files) & !grepl(".pos",files)]

    ## get and check pop info
    indv <- fread(indv_files[1], header=F)
    pop_info <- fread("sif.csv")
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
    save(data_cls, GT, map, pop, file=paste0("LD", min_LD, "cl", min.cl.size, ".RData"))
