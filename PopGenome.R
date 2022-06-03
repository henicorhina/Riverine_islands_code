#!/usr/bin/env Rscript


library(PopGenome)
setwd("/Volumes/Brumfield_Lab_Drive/River_islands")

# to hold results
res.df <- data.frame(species=character(),
                     haplotype.F_ST=double(), 
                     nucleotide.F_ST=double(), 
                     Nei.G_ST=double(), 
                     Hudson.G_ST=double(), 
                     Hudson.H_ST=double(), 
                     Hudson.K_ST=double(), 
                     Dxy=double(), 
                     stringsAsFactors=FALSE) 

# testing
#spec <- "Campephilus_melanoleucos"

# for species that have '-AmazonOnly' filenames
for (spec in c("Campephilus_melanoleucos", "Cantorchilus_leucotis", "Formicarius_analis", "Formicarius_colma", "Glaucidium_brasilianum", 
  "Megascops_choliba", "Monasa_nigrifrons", "Pheugopedius_coraya", "Piaya_cayana", "Pipra_erythrocephala", "Saltator_coerulescens", 
  "Saltator_grossus", "Schiffornis_turdina", "Tachyphonus_cristatus", "Tachyphonus_luctuosus", "Trogon_collaris", "Trogon_rufus")) {

  print(spec)
  fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
  GENOME.class <- readData(fname, format = 'nexus')
  sp <- get.individuals(GENOME.class)[[1]]
  
  # population assignments based on DAPC clusters
  # GENOME.class <- set.populations(GENOME.class,list(
  #    c("Campephilus_melanoleucos_LSUMNS2802_0","Campephilus_melanoleucos_LSUMNS2802-1"),
  #    c("Campephilus_melanoleucos_LSUMNS26102_0", "Campephilus_melanoleucos_LSUMNS26102_1", "Campephilus_melanoleucos_MPEG16102_0",  
  #      "Campephilus_melanoleucos_MPEG16102_1",   "Campephilus_melanoleucos_AMNH12719_0",   "Campephilus_melanoleucos_AMNH12719_1",  
  #      "Campephilus_melanoleucos_FMNH433276_0",  "Campephilus_melanoleucos_FMNH433276_1",  "Campephilus_melanoleucos_FMNH456682_0", "Campephilus_melanoleucos_FMNH456682_1",
  #      "Campephilus_melanoleucos_LSUMNS33338_0", "Campephilus_melanoleucos_LSUMNS33338_1", "Campephilus_melanoleucos_MPEG11746_0",  
  #      "Campephilus_melanoleucos_MPEG11746_1",   "Campephilus_melanoleucos_MPEG16666_0",   "Campephilus_melanoleucos_MPEG16666_1",  
  #      "Campephilus_melanoleucos_MPEG7571_0",    "Campephilus_melanoleucos_MPEG7571_1",    "Campephilus_melanoleucos_USNMB22230_0", 
  #      "Campephilus_melanoleucos_USNMB22230_1" 
  #    )
  #    ))
  
  # set each individual as its own population
  GENOME.class <- set.populations(GENOME.class,unname(
    split(sp, rep(1:(length(sp)/2), each = 2))
    ))
  
  # set all individuals to the same population
  # GENOME.class <- set.populations(GENOME.class,list(sp))
  
  # Fst / Gst
  GENOME.class <- F_ST.stats(GENOME.class)
  # get.F_ST(GENOME.class, pairwise=TRUE)
  # head(get.F_ST(GENOME.class))
  res <- as.data.frame(get.F_ST(GENOME.class))
  Fst <- colMeans(res, na.rm = TRUE)
  
  # Dxy
  GENOME.class <- diversity.stats.between(GENOME.class)
  dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)
  
  res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
  colnames(res.sp)[1] <- "species"
  colnames(res.sp)[8] <- "Dxy"
  res.df <- bind_rows(res.df, res.sp)
  
}


# for species that DON'T have '-AmazonOnly' filenames
for (spec in c("Campephilus_rubricollis", "Celeus_flavus", "Celeus_grammicus", 
               "Crypturellus_undulatus", "Crypturellus_variegatus", "Glaucidium_hardyi", "Hylophylax_naevia", "Hylophylax_punctulata", 
               "Megascops_watsonii", "Monasa_morphoeus", "Myrmeciza_fortis", "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_myotherinus", 
               "Phaethornis_bourcieri", "Phaethornis_hispidus", "Piaya_melanogaster", "Pipra_filicauda", "Schiffornis_major", "Synallaxis_gujanensis", 
               "Synallaxis_rutilans", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus")) {
  
  print(spec)
  fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
  GENOME.class <- readData(fname, format = 'nexus')
  sp <- get.individuals(GENOME.class)[[1]]
  
  # set each individual as its own population
  GENOME.class <- set.populations(GENOME.class,unname(
    split(sp, rep(1:(length(sp)/2), each = 2))
  ))
  
  # set all individuals to the same population
  # GENOME.class <- set.populations(GENOME.class,list(sp))
  
  # Fst / Gst
  GENOME.class <- F_ST.stats(GENOME.class)
  # get.F_ST(GENOME.class, pairwise=TRUE)
  # head(get.F_ST(GENOME.class))
  res <- as.data.frame(get.F_ST(GENOME.class))
  Fst <- colMeans(res, na.rm = TRUE)
  
  # Dxy
  GENOME.class <- diversity.stats.between(GENOME.class)
  dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)
  
  res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
  colnames(res.sp)[1] <- "species"
  colnames(res.sp)[8] <- "Dxy"
  res.df <- bind_rows(res.df, res.sp)
  
}


# river islands birds (different folder names)
for (spec in c("Conirostrum_bicolor", "Conirostrum_margaritae", "Cranioleuca_vulpecula", "Dendroplex_kienerii", "Elaenia_pelzelni", "Furnarius_minor", 
               "Knipolegus_orenocensis", "Talaphorus_chlorocercus", "Mazaria_propinqua", "Myrmoborus_lugubris", "Myrmochanes_hemileucus", 
               "Myrmotherula_assimilis", "Myrmotherula_klagesi", "Ochthornis_littoralis", "Serpophaga_hypoleuca", "Stigmatura_napensis", 
               "Thamnophilus")) {
  
  print(spec)
  fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
  GENOME.class <- readData(fname, format = 'nexus')
  sp <- get.individuals(GENOME.class)[[20]]
  
  # set each individual as its own population
  GENOME.class <- set.populations(GENOME.class,unname(
    split(sp, rep(1:(length(sp)/2), each = 2))
  ))
  
  # set all individuals to the same population
  # GENOME.class <- set.populations(GENOME.class,list(sp))
  
  # Fst / Gst
  GENOME.class <- F_ST.stats(GENOME.class)
  # get.F_ST(GENOME.class, pairwise=TRUE)
  # head(get.F_ST(GENOME.class))
  res <- as.data.frame(get.F_ST(GENOME.class))
  Fst <- colMeans(res, na.rm = TRUE)
  
  # Dxy
  GENOME.class <- diversity.stats.between(GENOME.class)
  dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)
  
  res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
  colnames(res.sp)[1] <- "species"
  colnames(res.sp)[8] <- "Dxy"
  res.df <- bind_rows(res.df, res.sp)
  
}




fname.out <- paste0("3_results/PopGenome_Fst.results.v1.csv")
write.csv(res.df, file = fname.out, row.names = FALSE)


