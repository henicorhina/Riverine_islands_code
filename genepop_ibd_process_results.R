#!/usr/bin/env Rscript

library(dplyr)
library(stringr)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/genepop_final_all_snps/")

# insufficient data: "elaenia_pelzelni", 
species <- c("conirostrum_bicolor", "conirostrum_margaritae", "cranioleuca_vulpecula", "dendroplex_kienerii", 
             "furnarius_minor", "knipolegus_orenocensis", "leucippus_chlorocercus", "mazaria_propinqua", "myrmoborus_lugubris", 
             "myrmochanes_hemileucus", "myrmotherula_assimilis", "myrmotherula_klagesi", "ochthornis_littoralis", "serpophaga_hypoleuca", 
             "stigmatura_napensis", "thamnophilus", "Campephilus_melanoleucos", "Campephilus_rubricollis", "Cantorchilus_leucotis", 
             "Celeus_flavus", "Celeus_grammicus", "Crypturellus_undulatus", "Crypturellus_variegatus", "Formicarius_analis", 
             "Formicarius_colma", "Glaucidium_brasilianum", "Glaucidium_hardyi", "Hylophylax_naevia", "Hylophylax_punctulata", 
             "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus", "Monasa_nigrifrons", "Myrmeciza_fortis", "Myrmeciza_hyperythra", 
             "Myrmoborus_leucophrys", "Myrmoborus_myotherinus", "Phaethornis_bourcieri", "Phaethornis_hispidus", "Pheugopedius_coraya", 
             "Piaya_cayana", "Piaya_melanogaster", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens", "Saltator_grossus", 
             "Schiffornis_major", "Schiffornis_turdina", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus", 
             "Tachyphonus_luctuosus", "Trogon_collaris", "Trogon_rufus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus")

df.res <- data.frame(species=character(),
                       e_statistic=double(), 
                       e_slope=double(), 
                       a_statistic=double(), 
                       a_slope=double(), 
                       stringsAsFactors=FALSE) 

#testing
# spec <- "Cantorchilus_leucotis"

for (spec in species) {
  filename.e <- paste0("1_ibd_results/", spec, "_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.ISO.e")
  
    df.e <- read.table(filename.e, sep = '\t', header = FALSE)
    res.e <- df.e %>% dplyr::filter(str_detect(V1, "a = "))
    vals.e.1 <- unlist(strsplit(res.e[1,1],","))
    vals.e.2 <- as.data.frame(unlist(strsplit(vals.e.1, " ")))
    colnames(vals.e.2) <- "V1"
    vals.e.3 <- vals.e.2 %>% dplyr::filter(str_detect(V1, "0"))
    # t(vals.e.3)
  
    filename.a <- paste0("1_ibd_results/", spec, "_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.ISO.a")
    df.a <- read.table(filename.a, sep = '\t', header = FALSE)
    res.a <- df.a %>% dplyr::filter(str_detect(V1, "a = "))
    vals.a.1 <- unlist(strsplit(res.a[1,1],","))
    vals.a.2 <- as.data.frame(unlist(strsplit(vals.a.1, " ")))
    colnames(vals.a.2) <- "V1"
    vals.a.3 <- vals.a.2 %>% dplyr::filter(str_detect(V1, "0"))
    # t(vals.a.3)
  
  res.sp <- bind_cols(str_to_title(spec), 
                      as.double(vals.e.3[1,1]), as.double(vals.e.3[2,1]), 
                      as.double(vals.a.3[1,1]), as.double(vals.a.3[2,1]))

  # clear results
  vals.a.3 <- vals.e.3 <- NA
  
  colnames(res.sp) <- colnames(df.res)
  df.res <- bind_rows(df.res, res.sp)
  
}
  



fname.out <- "/Volumes/Brumfield_Lab_Drive/River_islands/3_results/genepop_ibd.results.v1.csv"
write.csv(df.res, file = fname.out, row.names = FALSE)


  