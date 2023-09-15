# SLRfinder

This method is aimed to identify candidate sex-linked regions (SLRs) based on linkage and heterozygosity using SNP genotypes. Individual sexes can be used to further validate the candidate regions but are not required for this method to identify candidates. The method is written in R scripts and can be readily applied to any vcf datasets.  

Required R packages: igraph, data.table, SNPRelate, ggplot2, ggpubr, cowplot, parallel

**Step0: generate the input vcf dataset and the LD edge list**

Create a directory for each dataset using the dataset name. The dataset folder should contain two files: 

1. the list of samples to keep (named as dataset.list)

2. the sample information table (named as dataset.csv), including at least two columns named "SampleID" (basically the same as the dataset.list) and "Population".

The vcf dataset can be filtered using VCFtools, e.g.:

```
$ vcftools --vcf myinput.vcf \
--minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 \
--recode --recode-INFO-all --out myoutput
```

or using popoulations in Stacks, e.g.:

```
$ populations -P ./ -M popmap --min-maf 0.15 -R 0.75 --ordered-export --vcf
```

Then LD is calculated using the retained SNPs in VCFtools. The filtering and LD estimation can be down on all SNPs combined or separately on each chromosome if the dataset is too large.
```
$ vcftools --vcf mydata.vcf --geno-r2 --ld-window 100 --out mydata
```
The output LD edge list (mydata.geno.ld) will be used in the R codes below. 

**Step1: get LD clusters**

Run the R codes below in the dataset folder

```
## input data
# sample information 
sif = read.csv("mydata.csv")
# genome information (if contig names in column1 differ from the chromosome names in column2)
LG = read.table("reference", header = F)
names(LG) = c("chr", "CHR")
# LD edge list
geno.LD <- read.table("mydata.geno.ld", header = T)
names(geno.LD) = c("CHR", "from", "to", "N_INDV", "r2")

## set parameters
# default for whole-genome sequencing data
min_LD=0.85
min.cl.size=20 
# use loser thresholds for RADseq data which are much sparser
min_LD=0.2
min.cl.size=5

source("SLRfinder_functions.R")
system(paste0("mkdir ", "LD", min_LD*10, "cl", min.cl.size))
setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
system(paste0("mkdir ", "whitelist"))
data_cls <- NULL
for (i in 1:nrow(LG)) {
  chr = LG[i, "chr"]
  data = geno.LD[geno.LD$CHR == chr, ]
  out = get_single_LD_cluster(data, min_LD = min_LD, min.cl.size=min.cl.size)
  position = as.data.frame(unlist(out$SNPs))
  position = cbind(rep(chr, sum(out$nSNPs)), position)
  write.table(position, paste0("./whitelist/position.LG", i, ".list"), sep="\t", quote = F, row.names = F)
  data_cls <- rbind(data_cls, out)
}
#range(data_cls$nSNPs)
data_cls$SNPs = apply(data_cls, 1, function(cl){ paste0(cl$chr, "_", cl$SNPs)})
saveRDS(data_cls, file="data_cls.rds")
#setwd("../")
```


ls | grep cl | grep LD > subfolders
module load vcftools

while read subfolder
do
cd ${subfolder}
mkdir file012

for chr in {1..37}
do
vcftools --vcf ../../stacks2/populations.snps.vcf \
--positions ./whitelist/position.LG${chr}.list \
--012 \
--out ./file012/McK2020_LG${chr}_a15m75_LD${LD}cl${cl}
done

cd ../
done < subfolders



![image](https://github.com/xuelingyi/SLRfinder/assets/76267685/1e3c41c5-b9c5-4e34-82dd-5254f1d18fc7)

**Step2: identify candidate SLRs**






