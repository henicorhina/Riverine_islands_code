
setwd("/Users/home/Documents/Research/River_Islands/results/")

all = read.csv("results.csv", header = TRUE, row.names = 1)
head(all)

# compute mass corrected substitution rate per method of Nabholz et al 2016
# this is for 3rd codon position using the slope and intercept of their 'analysis 2'
# and the min and max 95% confidence interval slope and intercepts
slope_min95 <- -.141
slope <- -.145
slope_max95 <- -.143
intercept_min95 <- .367
intercept <- .459
intercept_max95 <- .559

#test mean mass of zosterops, result should be 0.0206 for mean slope and intercept
mean_z <- 10

# mean mass of your organism
mean <- 30.0
(10^(slope*log10(mean)+(intercept)))/100


substitution_rate_min95 <- (10^(slope_min95*log10(mean)+(intercept_min95)))/100
substitution_rate <- (10^(slope*log10(mean)+(intercept)))/100
substitution_rate_max95 <- (10^(slope_max95*log10(mean)+(intercept_max95)))/100

print(substitution_rate_min95)
print(substitution_rate)
print(substitution_rate_max95)






species <- c("Campephilus_melanoleucos", "Campephilus_rubricollis", "Cantorchilus_leucotis", "Celeus_flavus", "Celeus_grammicus", "conirostrum_bicolor_full", "conirostrum_bicolor", "conirostrum_margaritae", "cranioleuca_vulpecula", "Crypturellus_undulatus", "Crypturellus_variegatus", "dendroplex_kienerii", "elaenia_pelzelni", "Formicarius_analis", "Formicarius_colma", "furnarius_minor_full", "furnarius_minor", "Glaucidium_brasilianum", "Glaucidium_hardyi", "Hylophylax_naevia", "Hylophylax_punctulata", "knipolegus_orenocensis", "leucippus_chlorocercus", "mazaria_propinqua", "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus", "Monasa_nigrifrons", "Myrmeciza_fortis", "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "myrmoborus_lugubris", "Myrmoborus_myotherinus", "myrmochanes_hemileucus", "myrmotherula_assimilis", "myrmotherula_klagesi", "ochthornis_littoralis", "Phaethornis_bourcieri", "Phaethornis_hispidus", "Pheugopedius_coraya", "Piaya_cayana", "Piaya_melanogaster", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens", "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "serpophaga_hypoleuca", "stigmatura_full", "stigmatura_napensis", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus", "Tachyphonus_luctuosus", "thamnophilus", "Trogon_collaris", "Trogon_rufus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus")
for (s in species){
  mean <-
  substitution_rate <- (10^(slope*log10(mean)+(intercept)))/100
  
  
}


m <- subset(all, Av_mass == "Campephilus_melanoleucos")
c('true_read_depth')
c <- all[ , all$species]
m <- all[all$ses == 1, ]
  