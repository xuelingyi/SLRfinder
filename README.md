# SLRfinder

This method is aimed at identifying candidate sex-linked regions (SLRs) based on heterozygosity and linkage disequilibrium (LD) using SNP genotypes of populations (in VCF format). SLR candidates are identified among LD clusters based on ranks of parameters (see details below). If sexes are known for some or all individual samples, then the sex information can also be incorporated to filter LD clusters. Please refer to the bioRxiv preprint for the citation of this method.   

Required R packages: igraph, data.table, SNPRelate, ggplot2, ggpubr, cowplot, parallel

## Prepare the input files
Create a directory for each dataset using the dataset name (e.g., mydata). This dataset folder should contain:
1. **mydata.csv**: the sample information file named after the dataset. This should include at least two columns named "SampleID" (same as names in the sequencing data) and "Population". Sex information (e.g., male, female, unknown) can be included in a column named "sex".
2. **reference.list**: the genome information including two space-delimited columns: column1 is the contig/scaffold ID in the reference genome, column2 is the more informative chromosome names (e.g., LGx). The two columns can be identical if the contigs have already been renamed into the human-informative version in the vcf file. If using the whole-genome resequencing (WGS) data, the unassembled contigs may not be necessary to be included in the analyses and do not need to be listed in this file. 
3. **SLRfinder_functions.r**: the R script for running SLRfinder, available on Github
4. **./a15m75**: a folder containing the filtered vcf files (can be generated using the script below)
5. **./GenoLD.snp100**: a folder containing the LD edge lists estimated using the filtered vcf files (can be generated using the script below)

The file datasets_sif.zip includes sample information and reference lists of the datasets used in our paper, which can be example inputs. 
<br>

**prepare the input vcf datasets and the LD edge lists**

This is an example **unix** script for vcf filtering and LD estimation. If using large datasets (e.g., WGS), it will be faster to process data in parallel by chromosome (if using chromosome-level reference genomes) or by contig/scaffold (if using low-quality genomes). Run the following scripts in the dataset folder that should also contain the unfiltered SNP genotypes (**mydata.vcf**). This script will generate a folder **a15m75** to save the filtered vcf files, and a folder **GenoLD.snp100** to save the LD edge lists (mydata_LGx_a15m75.geno.ld, or mydata_a15m75.geno.ld).
```
## create folders to save the output of processed data of each chromosome
mkdir a15m75 GenoLD.snp100

## below is an array job 
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $1}')
lg=$(sed -n ${SLURM_ARRAY_TASK_ID}p reference.list | awk '{print $2}')
ind=mydata

## SNP filtering (e.g., input data can be the outputs from GATK GenotypeGVCFs)
bcftools view -m2 -M2 -v snps --min-ac=1 ../vcf/${ind}_${lg}.vcf.gz \
| vcftools --vcf - --minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 --recode --recode-INFO-all --out a15m75/${ind}_${lg}_a15m75


## LD edge list
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
## NB! This script requires that vcftools is installed and can be directly called in your R environment using system
source("SLRfinder_scripts.R")
```

<br>

<img width="840" alt="Screenshot 2024-01-16 at 1 21 03 PM" src="https://github.com/xuelingyi/SLRfinder/assets/76267685/83c45375-0c4a-44cc-aace-7acf16563ddc">

<br>

## Detailed tutorial
Below is a more detailed tutorial that can be followed step by step. Basically, the script SLRfinder_scripts.R is divided into sections with more detailed explanations on each section. 

**read data information**

The basic information for SLRfinder to process the data. 
```
## dataset name
mydata = "mydata"

## Run the R scripts in the dataset folder
setwd(mydata)
sif = read.csv(paste0(mydata, ".csv"))
# get genome information (contig names in column1, chromosome names in column2)
LG = read.table("reference.list", header = F)
names(LG) = c("chr", "lg")

sex_info=F
myranks=c("Dext_var_rank", "R2_rank","nSNPs_rank", "chi2_rank")

# default parameters for whole-genome sequencing data
min_LD=0.85
min.cl.size=20
ncores=1 # use more cores can speed up the analyses
# use loser thresholds (e.g., min_LD=0.2, min.cl.size=5) for RADseq data which are much sparser

source("SLRfinder_functions.r")
```

**Step1: identify LD clusters**

Use the LD edge list and the parameters set above to identify LD clusters. This script processes the input data by chromosome (so that it can be faster when processing WGS data), but it can also be easily modified to process all data combined (see the notes in the codes below). This script will generate: 
1. a parameter-named folder (e.g., **LD8.5cl20**) that contains information of the identified LD clusters (**data_cls.rds**)
2. a subfolder (**LD8.5cl20/whitelist**) that contains positions (by chromosome) of the SNPs included in the identified LD clusters
3. a subfolder (**LD8.5cl20/file012**) that contains the .012 genotypes of the SNPs (by chromosome) included in the LD clusters

This step can take a while so it might be good to use more cores.
```
########## step 1. identify LD clusters ##########
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
```

**Step2: process LD clusters**

This script will get the 012 genotypes of all LD clusters and generate two files in the dataset/parameter folder: a genotype file (**GT.RData**) and a file of the processed data including estimated criteria parameters (**data_all.rds**). 
```
####### step 2 process LD clusters #######
## if starting R from new: read data information again (see scripts above) 
## run the script below within the dataset folder. 
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
```

**Step3. Identify candidate SLRs and plot the results**

If sex information is known, the LD clusters will be filtered and only clusters having <10% misplaced sexes among all samples of known sexes will be retained and plotted (**mydata_sexg.pdf**).

Candidate LD clusters will also be identified based on rankings of the default four criteria parameters:

1. nSNPs: the number of SNPs in this cluster. The true SLR cluster should have more SNPs. 
2. R2: the adjusted R-squared values of the fitted regression of heterozygosity to PC1 scores. The true SLR cluster should have a higher R2. 
3. chi2: the chi-square value in the goodness of fit test on the even split of individuals into two clusters in each population. The true SLR cluster should have an even split of individuals into the two clusters representing males and females, thus a smaller chi2 (not different from the null expectation).
4. Dext_var: on the heterozygosity~PC1 plot, the Euclidean distance between each dot and its closer corner (bottom-left or top-right) divided by half of the distance between the corners. The true SLR cluster should have a small variance of this distance (i.e., dense clustering). 

This step can be rerun on the processed data file (**data_all.rds**) using different ranks that might get better results than default. 

The script will write out data of the candidate regions (**cand_regions.rds**, **mydata_can.csv**) in the dataset/parameter folder and plotted the results (**mydata_LD0.85cl20.pdf**).

```
####### step 3. identify SLR candidates #######
print("step 3 identify SLR candidates")
## if starting R from new: the data information needs to be read in again; run the script below within the dataset folder. 
# setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
data_all = readRDS("data_all.rds")

## print sex-related regions if know sex_info
if(sex_info){
  print("filter LD clusters by sex (10% misgrouped)")
  sex_filter = 0.1
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

## identify candidate LD clusters based on ranks (choose among: nSNPs_rank, R2_rank, chi2_rank, Dext_var_rank, Dext_max_rank, Dext_mean_rank).
# this step can be rerun using different rank combinations
# data_all = readRDS("data_all.rds")
print(paste0("identify LD clusters by ranks: ", paste0(myranks, collapse = ", ")))
cand_regions <- get_candidate_regions(data_all, ranks=myranks, nPerm=10000, cores=ncores, alpha=0.05)
saveRDS(cand_regions, "cand_regions.rds")

## plot results: the QQ-plot and all significant candidates or the five top-ranked LD clusters will be plotted in one PDF file
# cand_regions = readRDS("cand_regions.rds")
list2env(cand_regions, globalenv())
alpha=0.05
lambda <- lm(obs~exp+0, cand_regions$qq_data)$coefficients
qq_data$col=rep("grey40", nrow(data_out))
qq_data$col[which(data_out$p_gc_adj<alpha)] <- "#ff9727"

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
  
  ## color individuals by population; can color by phenotypic sex if known
  if(sex_info){
    print (ggplot(pca, aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") + geom_point(aes(x=PC_scaled, y=Het, color=sex), alpha=0.6, size=2.5) + theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme(title = element_text(size=10)))
  } else {
    print (ggplot(pca, aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") + geom_point(aes(x=PC_scaled, y=Het, color=Pop), alpha=0.6, size=2.5) + theme_bw() + labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme(title = element_text(size=10)))
  }

### extract identified genetic sex
#pca$SLRfinder_sex = "unknown"
#pca[pca$PC_scaled < 0.5, "SLRfinder_sex"] = "homogametic"
#pca[pca$PC_scaled > 0.5, "SLRfinder_sex"] = "heterogametic"
#sif = merge(sif, pca[, c("Ind", "SLRfinder_sex")], by.x="SampleID", by.y="Ind", all.x=T, sort=F)
#ggplot(pca[order(pca$SLRfinder_sex, decreasing = T),], aes(PC_scaled,Het)) + geom_smooth(method = "lm",se=FALSE, col="black") + geom_point(aes(x=PC_scaled, y=Het, color=SLRfinder_sex), alpha=0.6, size=2.5) +
#labs(x="PC1 (scaled)", y="Proportion heterozygous loci", title=title) + theme_bw() + theme(title = element_text(size=10)))
}
dev.off()

```


