library("phytools", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

setwd('/Volumes/Brumfield_Lab_Drive/data/1_analysis/UCE_gene_trees/')
filename="/Volumes/Brumfield_Lab_Drive/data/1_analysis/UCE_gene_trees/gene_tree_output.txt"
cat("species", "av_UCE_gene_tree_depth", "\n", file=filename, append=FALSE, sep='\t')
species <- c("conirostrum_bicolor", "conirostrum_margaritae", "cranioleuca_vulpecula", "dendroplex_kienerii", "elaenia_pelzelni", "furnarius_minor", "furnarius_minor_full", "knipolegus_orenocensis", "leucippus_chlorocercus", "mazaria_propinqua", "myrmoborus_lugubris", "myrmochanes_hemileucus", "myrmotherula_assimilis", "myrmotherula_klagesi", "ochthornis_littoralis", "serpophaga_hypoleuca", "stigmatura_napensis", "stigmatura_napensis_full", "thamnophilus", "Campephilus_melanoleucos", "Campephilus_rubricollis", "Cantorchilus_leucotis", "Celeus_flavus", "Celeus_grammicus", "Crypturellus_undulatus", "Crypturellus_variegatus", "Formicarius_analis", "Formicarius_colma", "Glaucidium_brasilianum", "Glaucidium_hardyi", "Hylophylax_naevia", "Hylophylax_punctulata", "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus", "Monasa_nigrifrons", "Myrmeciza_fortis", "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_myotherinus", "Phaethornis_bourcieri", "Phaethornis_hispidus", "Pheugopedius_coraya", "Piaya_cayana", "Piaya_melanogaster", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens", "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus", "Tachyphonus_luctuosus", "Trogon_collaris", "Trogon_rufus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus")

for (s in species){
  files <- list.files(path=s, pattern="*.tre", full.names=TRUE, recursive=FALSE)
  l<-c()
  for (i in files){
    t <- read.tree(i)
    temp <- sum(nodeHeights(t))/length(t$tip.label)
    l <- c(l, temp)
    }
  cat(s, mean(l), "\n", file=filename, append=TRUE, sep='\t')
  print(s)
  length(l)
  mean(l)
  }

#files <- list.files(path="/Volumes/Brumfield_Lab_Drive/data/1_analysis/UCE_gene_trees/conirostrum_bicolor", pattern="*.tre", full.names=TRUE, recursive=FALSE)
#l <- vector("list", length(files))
#l<-c()
#for (i in files){
#  t <- read.tree(i)
#  temp <- sum(nodeHeights(t))/length(t$tip.label)
#  l <- c(l, temp)
#}
#length(l)
#mean(l)

