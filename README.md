# SLRfinder

This method is aimed to identify candidate sex-determining regions (SDRs) based on linkage and heterozygosity using SNP genotypes. Individual sexes can be used to further validate the candidate regions but are not required in this method to identify candidates. The method is written in R scripts and can be readily applied to other datasets.  

Required R packages: igraph, data.table, SNPRelate, ggplot2, cowplot, parallel

Software used to process SNP data: samtools (bcftools), vcftools

Data were obtained from the published whole-genome sequencing of 887 wild individuals of Pungitius pungitius in Feng et al. (Feng, X., Merilä, J., & Löytynoja, A. (2022). Complex population history affects admixture analyses in nine‐spined sticklebacks. Molecular Ecology, 31(20), 5386-5401. https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16651). 

Datasets are generated based on the characterized sex-determining regions (SDRs) in Yi et al. (in prep). 

**Step0: generate the input vcf dataset and estimate LD**
The vcf dataset can be filtered using VCFtools, e.g.:
$ vcftools --vcf myinput.vcf --minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 --recode --recode-INFO-all --out myoutput

or using Stacks - popoulations, e.g.:
$ populations -P ./ -M popmap --min-maf 0.15 -R 0.75 --ordered-export --vcf

Then LD is calculated using VCFtools:
$ vcftools --vcf mydata.vcf --geno-r2 --ld-window 100 --out mydata

**Step1: get LD clusters **




