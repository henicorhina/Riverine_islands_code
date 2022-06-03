#!/usr/bin/env Rscript

library("genepop")

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/genepop_final_all_snps/")

#--------------------------------------------------------------------------------------------
# # testing. these worked
# locinfile <- 'Campephilus_rubricollis_SNPs_phased_rmIndels_75_QC_DP_random_GenePop_all.IBD.txt'
# outfile <- 'Campephilus_rubricollis_SNPs_phased_rmIndels_75_QC_DP_random_GenePop_all.ISO'
# out <- ibd(locinfile,
#            outputFile = outfile, 
#            statistic = 'a')
# # run this, with sizes = true (allelic Fst) and = false (overall Fst)
# locinfile <- 'conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_random_GenePop_all.txt'
# outfile <- 'conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_random_GenePop_all.fst'
# Fst(locinfile,
#     outputFile = outfile,
#     sizes = TRUE,
#     pairs = TRUE)



#--------------------------------------------------------------------------------------------

# all species
species <- c("conirostrum_bicolor", "conirostrum_margaritae", "cranioleuca_vulpecula", "dendroplex_kienerii",
    "furnarius_minor", "knipolegus_orenocensis", "leucippus_chlorocercus", "mazaria_propinqua", "myrmoborus_lugubris", "myrmochanes_hemileucus",
    "myrmotherula_assimilis", "myrmotherula_klagesi", "ochthornis_littoralis", "serpophaga_hypoleuca", "stigmatura_napensis", "thamnophilus",
    "Campephilus_melanoleucos", "Campephilus_rubricollis", "Cantorchilus_leucotis", "Celeus_flavus", "Celeus_grammicus", "Crypturellus_undulatus",
    "Crypturellus_variegatus", "Formicarius_analis", "Formicarius_colma", "Glaucidium_brasilianum", "Glaucidium_hardyi", "Hylophylax_naevia",
    "Hylophylax_punctulata", "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus", "Monasa_nigrifrons", "Myrmeciza_fortis",
    "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_myotherinus", "Phaethornis_bourcieri", "Phaethornis_hispidus",
    "Pheugopedius_coraya", "Piaya_cayana", "Piaya_melanogaster", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens",
    "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus",
    "Tachyphonus_luctuosus", "Trogon_collaris", "Trogon_rufus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus")

# working list
# species <- c("myrmotherula_klagesi", "ochthornis_littoralis", "serpophaga_hypoleuca", "stigmatura_napensis", "thamnophilus", 
#              "Celeus_flavus", "Celeus_grammicus", "Crypturellus_undulatus", "Crypturellus_variegatus", "Formicarius_analis", 
#              "Formicarius_colma", "Glaucidium_brasilianum", "Glaucidium_hardyi", "Hylophylax_naevia", "Hylophylax_punctulata", 
#              "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus", "Monasa_nigrifrons", "Myrmeciza_fortis", 
#              "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_myotherinus", "Phaethornis_bourcieri", "Phaethornis_hispidus", 
#              "Pheugopedius_coraya", "Piaya_cayana", "Piaya_melanogaster", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens", 
#              "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "Synallaxis_gujanensis", "Synallaxis_rutilans", 
#              "Tachyphonus_cristatus", "Tachyphonus_luctuosus", "Trogon_collaris", "Trogon_rufus", "Xiphorhynchus_elegans", 
#              "Xiphorhynchus_obsoletus")

# done
# species <- c("myrmotherula_assimilis", "conirostrum_margaritae", "cranioleuca_vulpecula", "dendroplex_kienerii", "elaenia_pelzelni", 
#              "furnarius_minor", "knipolegus_orenocensis", "leucippus_chlorocercus", "mazaria_propinqua", "myrmoborus_lugubris", 
#              "myrmochanes_hemileucus", "conirostrum_bicolor", "Campephilus_melanoleucos", "Campephilus_rubricollis", "Cantorchilus_leucotis")

# insufficient data: "elaenia_pelzelni"

# testing
# i <- "myrmotherula_assimilis"

for (i in species){
    print(i)
    
    #IBD
    locinfile <- paste0(i, '_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.IBD.txt')
    outfile <- paste0('1_ibd_results/', i, '_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.ISO.a')
    ibd(locinfile,
        outputFile = outfile, 
        statistic = 'a')
    
    outfile <- paste0('1_ibd_results/', i, '_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.ISO.e')
    ibd(locinfile,
        outputFile = outfile, 
        statistic = 'e')
    
    #overall Fst
    outfile <- paste0('1_ibd_results/', i, '_SNPs_phased_rmIndels_75_QC_DP_GenePop_all.overall.fst')
    Fst(locinfile,
        outputFile = outfile,
        sizes = FALSE,
        pairs = TRUE)
    
}
