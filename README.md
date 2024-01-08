# SLRfinder

This method is aimed at identifying candidate sex-linked regions (SLRs) based on heterozygosity and linkage disequilibrium (LD) using SNP genotypes in a vcf file. If sexes are known for some or all individual samples, then the sex information can be used to filter the LD clusters to identify candidate sex chromosomes. If sexes are unknown or the filtering step fails (e.g., no LD cluster remained or too many clusters remained after filtering), ranks of several other parameters can be used to identify candidate sex-linked LD clusters. 

Required R packages: igraph, data.table, SNPRelate, ggplot2, ggpubr, cowplot, parallel

## Input files
Create a directory for each dataset using the dataset name (e.g., mydata). This dataset folder should contain:
1. **mydata.csv**: the sample information file named after the dataset. This should include at least two columns named "SampleID" (same as names in the sequencing data) and "Population". Sex information (e.g., male, female, unknown) can be included in a column named "sex".
2. **reference.list**: the genome information including two space-delimited columns: column1 is the contig/scaffold ID in the reference genome, column2 is the more informative chromosome names (e.g., LGx). The two columns can be identical if the contigs have already been renamed into the human-informative version in the vcf file. If using the whole-genome resequencing (WGS) data, the unassembled contigs may not be necessary to be included in the analyses and do not need to be listed in this file. 
3. **SLRfinder_functions.r**: the R script for running SLRfinder, available on Github
4. **./a15m75**: a folder containing the filtered vcf files (can be generated using the step0_script below)
5. **./GenoLD.snp100**: a folder containing the LD edge lists estimated using the filtered vcf files (can be generated using the step0_script below)
<br>

**Step0: prepare the input vcf datasets and the LD edge lists**

This is an example **unix** script (can be adapted to R scripts) for vcf filtering and LD estimation. If using large datasets (e.g., WGS), it will be faster to process data in parallel by chromosome (if using chromosome-level reference genomes) or by contig/scaffold (if using low-quality genomes). Run the following scripts in the dataset folder that should also contain the unfiltered SNP genotypes (**mydata.vcf**). This script will generate a folder **a15m75** to save the filtered vcf files, and a folder **GenoLD.snp100** to save the LD edge lists (mydata_LGx_a15m75.geno.ld, or mydata_a15m75.geno.ld).
```
## create folders to save the output of processed data of each chromosome
mkdir a15m75 GenoLD.snp100

## below is an array job 
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $1}')
lg=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $2}')
ind=mydata

## filtering 
vcftools --gzvcf ${ind}.vcf.gz \
--chr ${chr} --minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 \
--recode --recode-INFO-all --out a15m75/${ind}_${lg}_a15m75

## LD calculation
vcftools --vcf ./a15m75/${ind}_${lg}_a15m75.recode.vcf \
--geno-r2 --ld-window 100 \
--out ./GenoLD.snp100/${ind}_${lg}_a15m75
```
If the dataset is small (e.g., RADseq data), all chromosomes can be processed together. Filtered can be done popoulations in Stacks, e.g.:
```
populations -P ./ -M popmap --min-maf 0.15 -R 0.75 --ordered-export --vcf
vcftools --vcf populations.snps.vcf --geno-r2 --ld-window 100 --out mydata_a15m75
```
<br>

## Quick Start
SLRfinder can be readily applied to the prepared data by saving the R script **SLRfinder_scripts.R** in the dataset folder and run the R code below:
```
## change the working directory to the dataset folder
## This script has a readline prompt asking for the inputs (dataset name, min_LD, min.cl.size, ncores).
## If you don't like the prompt, you can modify the first few lines to specify your inputs.
source("SLRfinder_scripts.R")
```

<br>

## Detailed tutorial
Below is a more detailed tutorial that can be followed step by step. Basically, the script SLRfinder_scripts.R is divided into sections with more detailed explanations on each section. 

**read data information**

The basic information for SLRfinder to process the data. 
```
## dataset name
mydata = "mydata"

## Run the R scripts in the dataset folder
# setwd(mydata)
sif = read.csv(paste0(mydata, ".csv"))
# get genome information (contig names in column1, chromosome names in column2)
LG = read.table("reference.list", header = F)
names(LG) = c("chr", "lg")

# default parameters for whole-genome sequencing data
min_LD=0.85
min.cl.size=20
ncores=1
# use more cores can speed up the analyses
# use loser thresholds for RADseq data which are much sparser
# min_LD=0.2
# min.cl.size=5

source("SLRfinder_functions.r")
```

**Step1: get LD clusters**

Use the LD edge list and the parameters set above to identify LD clusters. This script processes the input data by chromosome (so that it can be faster when processing WGS data), but it can also be easily modified to process all data combined (see the notes in the codes below). This script will generate: 
1. a parameter-named folder (e.g., **LD8.5cl20**) that contains information of the identified LD clusters (**data_cls.rds**)
2. a subfolder (**LD8.5cl20/whitelist**) that contains positions (by chromosome) of the SNPs included in the identified LD clusters
3. a subfolder (**LD8.5cl20/file012**) that contains the .012 genotypes of the SNPs (by chromosome) included in the LD clusters

This step can take a while so it might be good to use more cores.
```
########## step 1. get LD clusters ##########
## if all chr combined 
# geno.LD <- read.table(paste0(mydata, "_a15m75.geno.ld"), header = T)
# names(geno.LD) = c("CHR", "from", "to", "N_INDV", "r2")

## if using scaffolds
# LG = LG[LG$chr %in% unique(geno.LD$CHR), ]

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
  data_cls <- rbind(data_cls, out)}}
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
  system(paste0("vcftools --gzvcf ../a15m75/", mydata, "_", lg, "_a15m75.recode.vcf --positions ./whitelist/position.", lg, ".list", " --012 --out ./file012/", mydata, "_", lg, "_a15m75_LD", min_LD, "cl", min.cl.size))
}

# The above vcftools script can also be done in unix using the code below
# module load vcftools
# ind=mydata
# param=LD8.5cl20
# cd ${param}
# mkdir file012
## loop by chromosome (here total 21 chromosomes)
# for chr in {1..21}
# do
# vcftools --gzvcf ../a15m75/${ind}_LG${chr}_a15m75.vcf.gz --positions ./whitelist/position.LG${chr}.list \
# --012 --out ./file012/${ind}_LG${chr}_a15m75_LD${LD}cl${cl}
# done
# cd ../
```

**Step2.1: process the LD clusters**

This step will This script will get the 012 genotypes of all LD clusters and generate two files in the dataset/parameter folder: a genotype file (**GT.RData**) and a data file (**data_all.rds**) that includes the estimated criteria parameters. 
```
####### step 2.1 process LD clusters #######
## if starting R from new: the data information needs to be read in again; run the script below within the dataset folder. 
# setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
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

data_all = get_data_output(data_cls, GT, map, pop, sex_info, cores=ncores)
saveRDS(data_all, "data_all.rds")
```

**Step2.2 identify candidate SLRs and plot the results**

If sex information is known, the LD clusters will be filtered and only clusters having <10% misplaced sexes among all samples of known sexes will be retained and plotted (**mydata_sexg.pdf**).

Candidate LD clusters will also be identified based on rankings of default four criteria parameters:

1. nSNPs: the number of SNPs in this cluster. The true SLR cluster should have more SNPs. 
2. R2: the correlation coefficient between heterozygosity and PC1 scores. The true SLR cluster should have a higher R2. 
3. Dext_var: on the heterozygosity~PC1 plot, the Euclidean distance between each dot and its closer corner (bottom-left or top-right) divided by half of the distance between the corners. The true SLR cluster should have a small variance of this distance (i.e., dense clustering). 
4. chi2: the chi-square value in the goodness of fit test on the even split of individuals into two clusters in each population. The true SLR cluster should have an even split of individuals into the two clusters representing males and females, thus a smaller chi2 (not different from the null expectation). 

The script will write out data of the candidate regions (**cand_regions.rds**, **mydata_can.csv**) in the dataset/parameter folder and plotted the results (**mydata_LD0.85cl20.pdf**).

```
####### step 2.2 identify SLR candidates #######
print("step 2.2 identify SLR candidates")
## if starting R from new: the data information needs to be read in again; run the script below within the dataset folder. 
# setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
## print sex-related regions if know sex_info
myranks=c("Dext_var_rank", "R2_rank","nSNPs_rank", "chi2_rank")
sex_info=T

if(sex_info){
  print("filter LD clusters by sex (10% misgrouped)")
  sex_filter = 0.1
  data_sex = data_all[data_all$Sex_g <= sex_filter,]
  myindex = length(grep("_", data_sex[1, "chr"])) + 2
  
  if(nrow(data_sex) > 0){
    print(paste0(nrow(data_sex), " cluster(s) remained on ", 
                 paste0(data_sex$chr, collapse = ", ")))
    
    
    pdf(paste0(mydata, "_sexg.pdf"))
    
    for (i in 1:nrow(data_sex)) {
      data = as.data.frame(data_sex$data[i])
      
      region=paste(LG[LG$chr == data_sex$chr[i], "lg"], 
                   paste(range(as.numeric(do.call(rbind, strsplit(data_sex$SNPs[i][[1]],"_",fixed=TRUE))[,myindex])),collapse =  "-"), sep=":")
      
      title = paste0(region, "\nnSNPs=",data_sex$nSNPs[i])
      
      print(ggplot(data, aes(x=PC_scaled, y=Het)) + geom_point(aes(color=sex), alpha=0.6, size=2.5) +
              geom_smooth(method = "lm",se=FALSE, col="black") +
              theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci",
                                title=title) + theme(title = element_text(size=10)))
    }
    dev.off()
  } else {
    print(paste0("No cluster remained after filtering sex_g <= ", sex_filter, "!"))
  }
} else {
  print("Sex info unknown. Identify candidates by ranks only.")
}

print(paste0("identify LD clusters by ranks: ", paste0(myranks, collapse = ", ")))
#data_all = readRDS("data_all.rds")
cand_regions <- get_candidate_regions(data_all, ranks=myranks, nPerm=10000, cores=ncores, alpha=0.05)
saveRDS(cand_regions, "cand_regions.rds")

## plot results
#cand_regions = readRDS("cand_regions.rds")
list2env(cand_regions, globalenv())
alpha=0.05
lambda <- lm(obs~exp+0, cand_regions$qq_data)$coefficients
qq_data$col=rep("grey40", nrow(data_out))
qq_data$col[which(data_out$p_gc_adj<alpha)] <- "#ff9727"
PCA_het_data = merge(PCA_het_data, sif, by.x="Ind", by.y="SampleID", sort=F)

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
  print (ggplot(pca, aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") +
           geom_point(aes(x=PC_scaled, y=Het, color=sex.y), alpha=0.6, size=2.5) + theme_bw() +
           labs(x="PC1 (scaled)", y="Proportion heterozygous loci",
                title=title) + theme(title = element_text(size=10)))
}
dev.off()
```

