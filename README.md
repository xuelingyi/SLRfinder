# SLRfinder

This method is aimed to identify candidate sex-linked regions (SLRs) based on linkage and heterozygosity using SNP genotypes. Individual sexes can be used to further validate the candidate regions but are not required for this method to identify candidates. The method is written in R scripts and can be readily applied to any vcf datasets.  

Required R packages: igraph, data.table, SNPRelate, ggplot2, ggpubr, cowplot, parallel


**Step0: generate the input vcf dataset and estimate LD**

create a directory for each dataset using the dataset name. The dataset folder should contain two files: 
1. the list of samples to keep (named as dataset.list)
2. the sample information (named as dataset.csv) including at least two columns named "SampleID" (basically the same as the dataset.list) and "Population".

Filtering of the used datasets is archived in the file "Input_datasets". Basically, the vcf dataset can be filtered using VCFtools, e.g.:

$ vcftools --vcf myinput.vcf --minGQ 20 --minQ 30 --maf 0.15 --max-missing 0.75 --recode --recode-INFO-all --out myoutput


or using popoulations in Stacks, e.g.:

$ populations -P ./ -M popmap --min-maf 0.15 -R 0.75 --ordered-export --vcf


Then LD is calculated using VCFtools:

$ vcftools --vcf mydata.vcf --geno-r2 --ld-window 100 --out mydata


**Step1: get LD clusters**


default min.cl.size=20 for whole-genome sequencing SNPs; use larger sizes for stack loci (e.g., min.cl.size=100) so that stacks (should be physically linked) are not further broken.

**Step2: identify candidate SLRs**






