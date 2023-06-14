## run this R script within the dataset folder
source("../get_single_LD_cluster.R")

## test different parameter combinations
for (min_LD in c(0.8, 0.85, 0.9)) {
  for (min.cl.size in c(10, 20)) {

    ## create the paramter data folder and change directory inside
    system(paste0("mkdir ", "LD", min_LD*10, "cl", min.cl.size))
    setwd(paste0("LD", min_LD*10, "cl", min.cl.size))

    system(paste0("mkdir ", "whitelist"))

    ## generate data using the current parameters
    data_cls <- NULL
    for(i in c(1:21)){
      LD.data = paste0("../GenoLD.snp100/v7LG", i, "_a15m75.geno.ld")
      chr=paste0("LG", i)

      sink("./whitelist/mylog", append=T)
      out = get_single_LD_cluster(LD.data, min_LD = min_LD, min.cl.size=min.cl.size)
      print(paste0(out$chr, ": ", out$nSNPs))
      sink()

      position = as.data.frame(unlist(out$SNPs))
      position = cbind(rep(chr, sum(out$nSNPs)), position)
      write.table(position, paste0("./whitelist/position.", chr, ".list"), sep="\t", quote = F, row.names = F)

      data_cls <- rbind(data_cls, out)
    }

    data_cls$SNPs = apply(data_cls, 1, function(cl){ paste0(cl$chr, "_", cl$SNPs)})
    saveRDS(data_cls, file="data_cls.rds")
    setwd("../")
  }
}
