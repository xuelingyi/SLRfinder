# SLRfinder

This method is aimed to identify candidate sex-linked regions (SLRs) based on linkage and heterozygosity using SNP genotypes. Individual sexes can be used to further validate the candidate regions but are not required for this method to identify candidates. The method is mostly written in R scripts and can be readily applied to any vcf datasets.  

Required R packages: igraph, data.table, SNPRelate, ggplot2, ggpubr, cowplot, parallel
<br/> </br>

**Step0: generate the input vcf dataset and the LD edge list**

Create a directory for each dataset using the dataset name. The dataset folder should contain:
1. the SNP genotypes in the vcf format (mydata.vcf or mydata.vcf.gz) 
2. the sample information file (dataset.csv), including at least two columns named "SampleID" (the same as the sample names in the sequencing data) and "Population"
3. the genome information (reference.list), including two columns of the contig/scaffold ID in the reference genome and the corresponding chromosome names that are more informative (e.g., LGx). 
<br/> </br>

If using large datasets (e.g., whole-genome resequencing), it will be faster to process data in parallel by chromosome (if using chromosome-level reference genomes; the unassembled contigs may not be necessary to be included) or by contig/scaffold (if using low-quality genomes). See below for an example of filtering and LD estimation.  
```
## create folders to save the output of processed data of each chromosome
mkdir a15m75 GenoLD.snp100

## run below in an array job
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $1}')
lg=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $2}')
ind=mydata

## filtering 
bcftools view --regions ${chr} ${ind}.vcf.gz \
| vcftools --vcf - --minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 --recode --recode-INFO-all -c \
| bcftools view -Oz -o ./a15m75/${ind}_${lg}_a15m75.vcf.gz

## LD calculation
vcftools --gzvcf ./a15m75/${ind}_${lg}_a15m75.vcf.gz \
--geno-r2 --ld-window 100 \
--out ./GenoLD.snp100/${ind}_${lg}_a15m75
```
<br/> </br>

If the dataset is small (e.g., RADseq data), all chromosomes can be processed together. Filtered can be done popoulations in Stacks, e.g.:
```
populations -P ./ -M popmap --min-maf 0.15 -R 0.75 --ordered-export --vcf
vcftools --vcf populations.snps.vcf --geno-r2 --ld-window 100 --out mydata_a15m75
```
<br/> </br>

The output LD edge list (mydata_LGx_a15m75.geno.ld, or mydata_a15m75.geno.ld) will be input in the R codes below. 
<br/> </br>

**Step1: get LD clusters**

Run the R codes below in the dataset folder. If using large datasets (such as whole-genome resequencing), the script can be made faster by running on each chr in parallel (rather than in loops or combined as in the script below). 
```
####### input data information #######
## dataset name
mydata = "mydata"

# sample information 
sif = read.csv(paste0(mydata, ".csv"))
# genome information (contig names in column1, chromosome names in column2)
LG = read.table("reference.list", header = F)
names(LG) = c("chr", "lg")

## set parameters
# default for whole-genome sequencing data
min_LD=0.85
min.cl.size=20 
# use loser thresholds for RADseq data which are much sparser
# min_LD=0.2
# min.cl.size=5

source("SLRfinder_functions.R")

## if all chr combined 
# geno.LD <- read.table(paste0(mydata, "_a15m75.geno.ld"), header = T)
# names(geno.LD) = c("CHR", "from", "to", "N_INDV", "r2")

####### step 1. get LD clusters #######
system(paste0("mkdir ", "LD", min_LD*10, "cl", min.cl.size))
setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
system(paste0("mkdir ", "whitelist"))

data_cls <- NULL
for (i in 1:nrow(LG)) {
  lg = LG[i, "lg"]

  data = read.table(paste0("./GenoLD.snp100/", mydata, "_", lg, "_a15m75.geno.ld"), header = T)
  names(data) = c("CHR", "from", "to", "N_INDV", "r2")

  # if all chr combined 
  # data = geno.LD[geno.LD$CHR == LG[i, "chr"], ]

  out = get_single_LD_cluster(data, min_LD = min_LD, min.cl.size=min.cl.size)
  position = as.data.frame(unlist(out$SNPs))
  position = cbind(rep(lg, sum(out$nSNPs)), position)
  write.table(position, paste0("./whitelist/position.", lg, ".list"), sep="\t", quote = F, row.names = F)
  data_cls <- rbind(data_cls, out)
}

#range(data_cls$nSNPs)
data_cls$SNPs = apply(data_cls, 1, function(cl){ paste0(cl$chr, "_", cl$SNPs)})
saveRDS(data_cls, file="data_cls.rds")
setwd("../")
```
This script generates a parameter-named folder (e.g., LD8.5cl20) which contains a file named data_cls.rds (the LD clusters identified by this parameter), and a subfolder named whitelist which contains positions (by chromosome) of the SNPs included in these LD clusters. These SNPs were then extracted from the filtered vcf files and output in the 012 format using vcftools: 
```
# This unix script can be modified into an R script if VCFtools is installed on the local computer
# run this script within the dataset folder

module load vcftools
param=LD8.5cl20
ind=mydata

cd ${param}
mkdir file012

## loop by chromosome (here total 21 chromosomes)
for chr in {1..21}
do
vcftools --gzvcf ../a15m75/${ind}_LG${chr}_a15m75.vcf.gz \
--positions ./whitelist/position.LG${chr}.list \
--012 \
--out ./file012/${ind}_LG${chr}_a15m75_LD${LD}cl${cl}
done

cd ../
```
This script will generate in the parameter-named folder another subfolder named file012 containing the vcf files (by chromosome) of the SNPs included in the LD clusters.  
<br/> </br>

**Step2: identify candidate SLRs and plot the results**

Run the R script below within the dataset folder. 
```
##
####### skip this section if all data in step 1 are still active in the R environment #######
## input data information (same in step1)
## dataset name
mydata = "mydata"
# sample information 
sif = read.csv(paste0(mydata, ".csv"))
# genome information (contig names in column1, chromosome names in column2)
LG = read.table("reference", header = F)
names(LG) = c("chr", "CHR")

## parameters
min_LD=0.85
min.cl.size=20 
setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
## load the data of LD clusters generated in step 1
data_cls <- readRDS("data_cls.rds")

source("../SLRfinder_functions.R")

####### continue after the R script in step 1 #######
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

# check pop and sif info
indv <- fread(indv_files[1], header=F)
pop_info = sif[order(factor(sif$SampleID, levels = indv$V1)),]
if(all(indv$V1 == pop_info$SampleID)) {
  pop <- pop_info$Population
  ind <- pop_info$SampleID
}else{
  print("indv and pop do not match!")}

save(data_cls, GT, map, ind, pop, file="GT.RData")

### get the candidate LD clusters: output all significant candidates; if no significant result is found, output the five LD clusters that have the lowest adjusted p-values (the top-ranked LD clusters)
## the ranks used for defining candidate regions
ranks = c("Dext_max_rank", "R2_rank", "nSNPs_rank", "chi2_rank")
cand_regions <- get_candidate_regions(data_cls, GT, map, pop, ranks=ranks, nPerm=10000, cores=1, alpha=0.05)
saveRDS(cand_regions, "cand_regions.rds")
#print(cand_regions$candidates)

#cand_regions = readRDS("cand_regions.rds")
### plot results ###
alpha=0.05
list2env(cand_regions, globalenv())
lambda <- lm(obs~exp+0, cand_regions$qq_data)$coefficients
qq_data$col=rep("steelblue", nrow(data_out))
qq_data$col[which(data_out$p_gc_adj<alpha)] <- "indianred"
PCA_het_data = merge(PCA_het_data, sif, by.x="Ind", by.y="Run", sort=F)

candidates = merge(candidates, LG, by="chr")
write.csv(candidates[, c("chr", "CHR", "region", "p_gc_adj", "nSNPs", "rank",  "nSNPs_rank", "R2_rank", "Dext_max_rank", "chi2_rank")], paste0(mydata, "_can.csv"), row.names = F)
```

