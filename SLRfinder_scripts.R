## run this script in the dataset folder

mydata = readline(prompt="The dataset name: ")
min_LD = as.numeric(readline(prompt="min_LD (default 0.85): "))
min.cl.size = as.numeric(readline(prompt="min.cl.size (default 20): "))
myranks = readline(prompt="use default (Dext_var_rank, R2_rank, nSNPs_rank, chi2_rank)?")
sex_info = readline(prompt="Sex available? (default F)")
sex_filter = as.numeric(readline(prompt="filtering threshold of misplaced sexes (default 0.1): "))
my_sex_ratio = readline(prompt="sex ratio (space separated; default 0.5 0.5): ")
ncores = as.numeric(readline(prompt="cores to use (default 1): "))

if(myranks) { myranks = c("Dext_var_rank", "R2_rank", "nSNPs_rank", "chi2_rank")} else {
  myranks = readline(prompt="specify myranks (space separated): ")
  myranks = unlist(strsplit(myranks, split = " "))
}
my_sex_ratio = as.numeric(unlist(strsplit(my_sex_ratio, split = " ")))


#mydata = "mydata"
#min_LD=0.85
#min.cl.size=20
#ncores=1
#myranks=c("Dext_var_rank", "R2_rank", "nSNPs_rank", "chi2_rank")
#sex_info=F
#sex_filter = 0.1

########### read data information ###########
print("read data information")
## dataset name
#mydata = "mydata"

sif = read.csv(paste0(mydata, ".csv"))
# get genome information (contig names in column1, chromosome names in column2)
LG = read.table("reference.list", header = F)
names(LG) = c("chr", "lg")

# default parameters for whole-genome sequencing data
# min_LD=0.85
# min.cl.size=20 
# use loser thresholds for RADseq data which are much sparser
# min_LD=0.2
# min.cl.size=5

source("SLRfinder_functions.r")
########## step 1. get LD clusters ##########
print("step 1. get LD clusters")
## if all chr combined 
# geno.LD <- read.table(paste0(mydata, "_a15m75.geno.ld"), header = T)
# names(geno.LD) = c("CHR", "from", "to", "N_INDV", "r2")

system(paste0("mkdir ", "LD", min_LD*10, "cl", min.cl.size))
setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
system(paste0("mkdir ", "whitelist"))

data_cls <- NULL
for (i in 1:nrow(LG)) {
  chr = LG[i, "chr"]
  lg = LG[i, "lg"]
  
  data = read.table(paste0("../GenoLD.snp100/", mydata, "_", lg, "_a15m75.geno.ld"), header = T)
  names(data) = c("CHR", "from", "to", "N_INDV", "r2")
  
  # if all chr combined 
  # data = geno.LD[geno.LD$CHR == LG[i, "chr"], ]
  
  out = get_single_LD_cluster(data, min_LD = min_LD, min.cl.size=min.cl.size)
  
  if(!is.null(out)){
    if(nrow(out) != 0){
    position = as.data.frame(unlist(out$SNPs))
    position = cbind(rep(chr, sum(out$nSNPs)), position)
    write.table(position, paste0("./whitelist/position.", lg, ".list"), sep="\t", quote = F, row.names = F)
    data_cls <- rbind(data_cls, out)
  }}
}

#range(data_cls$nSNPs)
data_cls$SNPs = apply(data_cls, 1, function(cl){ paste0(cl$chr, "_", cl$SNPs)})
print(paste0("total number of LD clusters: ", nrow(data_cls)))
saveRDS(data_cls, file="data_cls.rds")
#setwd("../")

system("mkdir file012")
## note: this requires that vcftools has been installed in the local environment
for (i in 1:nrow(LG)){
  lg = LG[i, "lg"]
  system(paste0("vcftools --vcf ../a15m75/", mydata, "_", lg, "_a15m75.recode.vcf --positions ./whitelist/position.", lg, ".list", " --012 --out ./file012/", mydata, "_", lg, "_a15m75_LD", min_LD, "cl", min.cl.size))
}

####### step 2 process LD clusters #######
print("step 2 process LD clusters")
## if starting R from new: the data information needs to be read in again; run the script below within the dataset folder. 
# setwd(paste0("LD", min_LD*10, "cl", min.cl.size))

data_cls = readRDS("data_cls.rds")

files <- paste0("./file012/", list.files("file012"))
indv_files <- files[grep(".indv",files)]
pos_files <- files[grep(".pos",files)]
GT_files <- files[!grepl(".log",files) & !grepl(".indv",files) & !grepl(".pos",files)]

map <- rbindlist(lapply(1:length(pos_files),function(i){
  pos <- fread(pos_files[i], sep="\t")
  if(ncol(pos) > 1 ) {
    colnames(pos) <- c("Chr","Pos")
    pos$SNP <- paste0(pos$Chr, "_", pos$Pos)}
  if(ncol(pos) ==0 ) {pos = NULL}
  return(pos)
}))

GT <- do.call(cbind, lapply(1:length(pos_files),function(i){ 
  gt.matrix = as.matrix(fread(GT_files[i])[,-1]) 
  if (nrow(gt.matrix) == 0){return(NULL)}
  if (nrow(gt.matrix) > 0){return(gt.matrix)}
}))
#recode missing data
GT[GT==-1] <- NA

# check pop and sif info
indv <- fread(indv_files[1], header=F)
pop_info = sif[order(factor(sif$SampleID, levels = indv$V1)),]
if(all(indv$V1 == pop_info$SampleID)) {
  pop <- pop_info$Population
  ind <- pop_info$SampleID
}else{
  print("indv and pop do not match!")}

save(data_cls, GT, map, ind, pop, file="GT.RData")

data_all = get_data_output(data_cls, GT, map, pop, sex_info, heterog_homog = c(0.5, 0.5), cores=ncores)
saveRDS(data_all, "data_all.rds")
## if expecting skewed sex ratios, the expected probability of sampling the heterogametic sex and the homogametic sex can be specified using the heterog_homog parameter. This will impact the chi-square tests of goodness of fit.


####### step 3 identify SLR candidates #######
print("step 3 identify SLR candidates")
## if starting R from new: the data information needs to be read in again; run the script below within the dataset folder. 
# setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
data_all = readRDS("data_all.rds")

## print sex-related regions if know sex_info
if(sex_info){
  print(paste0("filter LD clusters by sex (", sex_filter*100, "% misgrouped)"))
  data_sex = data_all[data_all$Sex_g <= sex_filter,]
  myindex = length(grep("_", data_sex[1, "chr"])) + 2
  
  if(nrow(data_sex) > 0){
    if(length(unique(data_sex$chr)) > 10){
      print(paste0("Sex filtering retained more than 10 chr: ", paste(unique(data_sex$chr), collapse=", ")))
    } else {
      print(paste0(nrow(data_sex), " cluster(s) remained on ", paste0(data_sex$chr, collapse = ", ")))
      data_sex[,region:=apply(data_sex,1,function(x)region=paste(x$chr, paste(range(as.numeric(do.call(rbind,strsplit(x$SNPs,"_",fixed=TRUE))[,myindex])),collapse="-"),sep=":"))]
      
      pdf(paste0(mydata, "_sexg.pdf"))
      for (i in 1:nrow(data_sex)) {
        data = as.data.frame(data_sex$data[i])
        region=data_sex$region[i]
        title = paste0(region, "\nnSNPs=",data_sex$nSNPs[i])
        
        print(ggplot(data, aes(x=PC_scaled, y=Het)) + geom_point(aes(color=sex), alpha=0.6, size=2.5) + geom_smooth(method = "lm",se=FALSE, col="black") + theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme(title = element_text(size=10)))
      }
      dev.off()
      
      write.csv(data_sex[, c("chr", "region", "Sex_g", 
                             "nSNPs", "Dext_mean", "R2", "chi2", 
                             "nSNPs_rank", "Dext_mean_rank", "R2_rank", "chi2_rank", 
                             "Dext_max_rank", "mean_LD", "Dext_max")], 
                "sex_filter.csv", row.names = F)
    }
  } else {
    print(paste0("No cluster remained after filtering sex_g <= ", sex_filter, "!"))
  }
} else {
  print("Sex info unknown. Identify candidates by ranks only.")
}

#data_all = readRDS("data_all.rds")
print(paste0("identify LD clusters by ranks: ", paste0(myranks, collapse = ", ")))
cand_regions <- get_candidate_regions(data_all, ranks=myranks, nPerm=10000, cores=ncores, alpha=0.05)
saveRDS(cand_regions, "cand_regions.rds")

## plot results
#cand_regions = readRDS("cand_regions.rds")
list2env(cand_regions, globalenv())
alpha=0.05
lambda <- lm(obs~exp+0, cand_regions$qq_data)$coefficients
qq_data$col=rep("grey40", nrow(data_out))
qq_data$col[which(data_out$p_gc_adj<alpha)] <- "#ff9727"

## the script below will print all results (the QQ-plot and all significant candidates or the five
## top-ranked LD clusters) into one PDF file, with one plot per page
pdf(paste0(mydata, "_LD", min_LD, "cl", min.cl.size, ".pdf"), width = 6, height = 4)
print(ggplot(qq_data, aes(x=exp, y=obs)) +
        geom_point(col=qq_data$col) + theme_bw() + theme(legend.position = "none") +
        labs(title=paste0("min_LD=", min_LD, ", min.cl.size=", min.cl.size, "\n",
                          "lambda=", round(lambda,2)), x="Expected -log10(P)", y="Observed -log10(P)") +
        geom_abline(slope = 1, intercept = 0, linewidth=0.5) + geom_smooth(method = "lm", aes(x=exp, y=obs+0),
                                                                           col="#ff9727", linewidth=0.5))

for(r in unique(PCA_het_data$region)) {
  pca = PCA_het_data[PCA_het_data$region == r, ]
  
  label = unlist(strsplit(unique(pca$label), split=" | ", fixed = T))
  chr = unlist(strsplit(label[1],  ":"))[1]
  lg = LG[LG$chr == chr, "lg"]
  
  ## replace the original contig name with the more informative chromosome ID (LGx)
  title = paste0(sub(chr, lg, label[1]), "\n", label[2], " ", label[3])
  
  ## color individuals by population; can color by sex if known
  if(sex_info){
    print (ggplot(pca, aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") + geom_point(aes(x=PC_scaled, y=Het, color=sex), alpha=0.6, size=2.5) + theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme(title = element_text(size=10)))
  } else {
    print (ggplot(pca, aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") + geom_point(aes(x=PC_scaled, y=Het, color=Pop), alpha=0.6, size=2.5) + theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme(title = element_text(size=10)))
  }
}
dev.off()

setwd("../../")
