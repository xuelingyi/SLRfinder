# SDR_LD

Data were obtained from the published whole-genome sequencing of 887 wild individuals of Pungitius pungitius in Feng et al. (Feng, X., Merilä, J., & Löytynoja, A. (2022). Complex population history affects admixture analyses in nine‐spined sticklebacks. Molecular Ecology, 31(20), 5386-5401. https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16651). 

Datasets are generated based on the characterized sex-determining regions (SDRs) in Yi et al. (in prep). 

Total 13 datasets: 
1. ELpop27: 
#length(unique(ELpop27$Population))
# pop=27, n=558

## WL (LG3)
WLpop7 = sif[!(sif$Population %in% unclear) & sif$sex_region == "LG3", ]
#length(unique(WLpop7$Population))
# pop=7, n=146

## 11 by pops
for(pop in unclear) {
  assign(gsub("-", "_", pop), sif[sif$Population == pop, ])
}
unclear = c("RUS-LEN", "JAP-BIW","USA-HLA", "CAN-FLO", "CAN-TEM",  #non_EU
            "SCO-HAR", "GBR-GRO", "FRA-VEY", "SWE-NAV", #unknown
            "FIN-KRK", #female-only
            "POL-GDY" #2SDRs
            )
