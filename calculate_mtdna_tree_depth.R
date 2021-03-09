library("phytools", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

setwd("/Volumes/Backup_Plus/data/mitogenomes/Final_trees/")
filename="/Volumes/Backup_Plus/data/mitogenomes/mtdna_tree_length_output.txt"
cat("species", "av_mtDNA_tree_depth", "\n", file=filename, append=FALSE, sep='\t')
#done: "RAxML_bipartitions.Crypturellus_undulatus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Campephilus_melanoleucos_mtdna_final_trimmed.tre", "RAxML_bipartitions.Campephilus_rubricollis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Cantorchilus_leucotis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Celeus_flavus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Celeus_grammicus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Conirostrum_bicolor_mtdna_final_trimmed.tre", "RAxML_bipartitions.Conirostrum_margaritae_mtdna_final_trimmed.tre", "RAxML_bipartitions.Cranioleuca_vulpecula_mtdna_final_trimmed.tre", "RAxML_bipartitions.Dendroplex_kienerii_mtdna_final_trimmed.tre", 
species <- c("RAxML_bipartitions.Formicarius_analis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Formicarius_colma_mtdna_final_trimmed.tre", "RAxML_bipartitions.Furnarius_mtdna_final_trimmed.tre", "RAxML_bipartitions.Glaucidium_brasilianum_mtdna_final_trimmed.tre", "RAxML_bipartitions.Glaucidium_hardyi_mtdna_final_trimmed.tre", "RAxML_bipartitions.Hylophylax_naevia_mtdna_final_trimmed.tre", "RAxML_bipartitions.Hylophylax_punctulata_mtdna_final_trimmed.tre", "RAxML_bipartitions.Knipolegus_full_mtdna_final_trimmed.tre", "RAxML_bipartitions.Leucippus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Mazaria_mtdna_final_trimmed.tre", "RAxML_bipartitions.Megascops_choliba_mtdna_final_trimmed.tre", "RAxML_bipartitions.Megascops_watsonii_mtdna_final_trimmed.tre", "RAxML_bipartitions.Monasa_morphoeus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Monasa_nigrifrons_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmeciza_fortis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmeciza_hyperythra_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmoborus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmoborus_leucophrys_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmoborus_myotherinus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmochanes_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmotherula_assimilis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Myrmotherula_klagesi_mtdna_final_trimmed.tre", "RAxML_bipartitions.Ochthornis_full_mtdna_final_trimmed.tre", "RAxML_bipartitions.Phaethornis_bourcieri_mtdna_final_trimmed.tre", "RAxML_bipartitions.Phaethornis_hispidus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Pheugopedius_coraya_mtdna_final_trimmed.tre", "RAxML_bipartitions.Piaya_cayana_mtdna_final_trimmed.tre", "RAxML_bipartitions.Piaya_melanogaster_mtdna_final_trimmed.tre", "RAxML_bipartitions.Pipra_filicauda_mtdna_final_trimmed.tre", "RAxML_bipartitions.Pipra_rubrocapilla_mtdna_final_trimmed.tre", "RAxML_bipartitions.Saltator_coerulescens_mtdna_final_trimmed.tre", "RAxML_bipartitions.Saltator_grossus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Schiffornis_major_mtdna_final_trimmed.tre", "RAxML_bipartitions.Schiffornis_turdina_mtdna_final_trimmed.tre", "RAxML_bipartitions.Serpophaga_mtdna_final_trimmed.tre", "RAxML_bipartitions.Stigmatura_full_mtdna_final_trimmed.tre", "RAxML_bipartitions.Synallaxis_gujanensis_mtdna_final_trimmed.tre", "RAxML_bipartitions.Synallaxis_rutilans_mtdna_final_trimmed.tre", "RAxML_bipartitions.Tachyphonus_cristatus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Tachyphonus_luctuosus_mtdna_final_trimmed.tre", "RAxML_bipartitions.thamnophilus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Trogon_collaris_mtdna_final_trimmed.tre", "RAxML_bipartitions.Trogon_rufus_mtdna_final_trimmed.tre", "RAxML_bipartitions.Xiphorhynchus_elegans_mtdna_final_trimmed.tre", "RAxML_bipartitions.Xiphorhynchus_obsoletus_mtdna_final_trimmed.tre")
species <- c("RAxML_bipartitions.Crypturellus_variegatus_mtdna_final_trimmed.tre")

for (s in species){
  t <- read.tree(s)
  temp <- sum(nodeHeights(t))/length(t$tip.label)
  cat(s, temp, "\n", file=filename, append=TRUE, sep='\t')
  print(s)
  print(temp)
  #length(temp)
  #mean(temp)
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

