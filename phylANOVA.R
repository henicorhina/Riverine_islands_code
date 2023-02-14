#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(ape)
library(phytools)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(forcats)
library(caper)
library(qpcR)
library(geiger)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands")

theta <- read.csv("3_results/dendropy/theta.csv")
df <- read.csv("3_results/results_for_pgls_take4.formatted.csv")
new.fst <- read.csv("3_results/PopGenome_Fst.results.v1.DAPCassignments.csv")
island.Tree <- read.tree(file = "all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.newick.phy")
# island.Tree <- read.tree(file = "/Volumes/Brumfield_Lab_Drive/River_islands/all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.phy")
# island.Tree <- read.tree(file = "/Volumes/Brumfield_Lab_Drive/River_islands/all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.tre")
island.Tree <- drop.tip(island.Tree, "Thamnophilus_cryptoleucus2")
# avonet <- read.csv("/Volumes/Backup_Plus/AVONET/Supplementary_dataset_1_raw_data.csv") # doesn't include trophic data / diet
# avonet <- read.csv("/Volumes/Backup_Plus/AVONET/Supplementary_dataset_1.csv") # preliminary df
avonet <- read.csv("/Users/oscar/Documents/Projects/AVONET/ELEData/TraitData/AVONET_Raw_Data.csv") # final df, raw data
avonet.bl <- read.csv("/Users/oscar/Documents/Projects/AVONET/ELEData/TraitData/AVONET1_BirdLife.csv") # final df, birdlife
all_species_key <- read.csv("All_species_key.csv")
PopGen <- read.csv("3_results/PopGenome_Fst.results.v1.csv")
genepop <- read.csv("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/genepop_ibd.results.v1.csv")

#------------------------------------------------------------------------

phy_anova_island.char <- function(df.t, char) {
  # df.t = dataframe
  # char = column number of the character in question
  
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(habitat, species, char)
  temp <- temp %>% dplyr::filter(!is.na(char))
  
  habitats <- as.vector(temp$habitat)
  param <- as.vector(temp$char)
  names(habitats) <- temp$species
  names(param) <- temp$species
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  # res.island.aov <- phylANOVA(island.Tree.t, habitats, log(param), nsim=1000, posthoc=TRUE, p.adj="holm") #log transformed
  res.island.aov <- phylANOVA(island.Tree.t, habitats, param, nsim=1000, posthoc=TRUE, p.adj="holm")
  
  res_list <- list("anova" = res.island.aov, 
                   "p" = res.island.aov$Pf, 
                   "F" = res.island.aov$F,
                   "island.floodplain.P" = res.island.aov$Pt["island","floodplain"], 
                   "island.upland.P" = res.island.aov$Pt["island","upland"], 
                   "floodplain.upland.P" = res.island.aov$Pt["floodplain","upland"],
                   "island.floodplain.T" = abs(res.island.aov$T["island","floodplain"]), 
                   "island.upland.T" = abs(res.island.aov$T["island","upland"]), 
                   "floodplain.upland.T" = abs(res.island.aov$T["floodplain","upland"])
  )
  
  return(res_list)
}


phy_anova_island.islandsplit <- function(df.t, char) {
  # df.t = dataframe
  # char = column number of the character in question
  
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(habitat, species, char)
  temp <- temp %>% dplyr::filter(!is.na(char))
  
  habitats <- as.vector(temp$habitat)
  param <- as.vector(temp$char)
  names(habitats) <- temp$species
  names(param) <- temp$species
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  # res.island.aov <- phylANOVA(island.Tree.t, habitats, log(param), nsim=1000, posthoc=TRUE, p.adj="holm") #log transformed
  res.island.aov <- phylANOVA(island.Tree.t, habitats, param, nsim=1000, posthoc=TRUE, p.adj="holm")
  
  res_list <- list("anova" = res.island.aov, 
                   "p" = res.island.aov$Pf, 
                   "F" = res.island.aov$F,
                   "island.early.island.late.P" = res.island.aov$Pt["island.early","island.late"], 
                   "island.early.floodplain.P" = res.island.aov$Pt["island.early","floodplain"], 
                   "island.early.upland.P" = res.island.aov$Pt["island.early","upland"], 
                   "island.late.floodplain.P" = res.island.aov$Pt["island.late","floodplain"], 
                   "island.late.upland.P" = res.island.aov$Pt["island.late","upland"], 
                   "floodplain.upland.P" = res.island.aov$Pt["floodplain","upland"],
                   "island.early.island.late.T" = abs(res.island.aov$T["island.early","island.late"]), 
                   "island.early.floodplain.T" = abs(res.island.aov$T["island.early","floodplain"]), 
                   "island.early.upland.T" = abs(res.island.aov$T["island.early","upland"]),
                   "island.late.floodplain.T" = abs(res.island.aov$T["island.late","floodplain"]), 
                   "island.late.upland.T" = abs(res.island.aov$T["island.late","upland"]), 
                   "floodplain.upland.T" = abs(res.island.aov$T["floodplain","upland"])
  )
  
  return(res_list)
}

phy_anova_island.model <- function(df.t, char) {
  # don't use this
  # df.t = dataframe
  # char = column number of the character in question
  # df.t <- df.temp
  # char = colnum
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(habitat, species, char)
  temp <- temp %>% dplyr::filter(!is.na(char))
  
  habitats <- as.factor(temp$habitat)
  param <- as.vector(temp$char)
  names(habitats) <- temp$species
  names(param) <- temp$species
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  # res.island.aov <- phylANOVA(island.Tree.t, habitats, param, nsim=1000, posthoc=TRUE, p.adj="holm")
  res.island.aov <- aov.phylo(param~habitats, island.Tree.t, nsim=50)
  
  res.island.aov
}

pgls.char <- function(df.t, char) {
  # df.t = dataframe
  # char = column number of the character in question
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(habitat, species, char) %>% 
    dplyr::filter(!is.na(char)) %>%
    droplevels()
  
  temp$habitat <- gsub("island", "0", temp$habitat)
  temp$habitat <- gsub("floodplain", "1", temp$habitat)
  temp$habitat <- gsub("upland", "2", temp$habitat)
  temp$habitat <- as.numeric(temp$habitat)

  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  hab.dat <- comparative.data(island.Tree.t, temp, species)

  # data have been rescaled, so no need to use log()
  res.pgls <- pgls(habitat ~ char, hab.dat)
  
  # if (char == 25) { # absolute value for log-transform Tajima's D. all are negative
  #   res.pgls <- pgls(log(abs(char)) ~ habitat, hab.dat)
  # } else {
  #   # res.pgls <- pgls(log(char) ~ habitat, hab.dat)
  #   res.pgls <- pgls(habitat ~ log(char), hab.dat)
  # }
  # summary(res.pgls)
  # coef(res.pgls)
  
  
  # plot(log(temp$char) ~ temp$habitat)
  # abline(a = coef(res.pgls)[1], b = coef(res.pgls)[2])
  
  
  res_list <- list("pgls" = res.pgls, 
                   "estimate" = summary(res.pgls)[5]$coefficients[2,1], 
                   "p" = summary(res.pgls)[5]$coefficients[2,4], 
                   "t" = summary(res.pgls)[5]$coefficients[2,3], 
                   "aicc" = res.pgls$aicc,
                   "slope" = coef(res.pgls)[2],
                   "int" = coef(res.pgls)[1],
                   "sterr" = res.pgls$sterr)
  
  
  return(res_list)
}

pgls.hwi <- function(df.t, char) {
  # df.t = dataframe
  # char = column number of the character in question
  
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(Hand.wing.Index, species, char) %>% 
    dplyr::filter(!is.na(Hand.wing.Index)) %>%
    dplyr::filter(!is.na(char)) %>%
    droplevels()
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  hab.dat <- comparative.data(island.Tree.t, temp, species)
  
  if (char == 25) { # absolute value for log-transformed Tajima's D, because all are negative
    res.pgls <- pgls(log(abs(char)) ~ Hand.wing.Index, hab.dat)
  } else {
    res.pgls <- pgls(log(char) ~ Hand.wing.Index, hab.dat)
  }
  # summary(res.pgls)
  # coef(res.pgls)
  
  
  # plot(log(temp$char) ~ temp$habitat)
  # abline(a = coef(res.pgls)[1], b = coef(res.pgls)[2])
  
  
  res_list <- list("pgls" = res.pgls, 
                   "estimate" = summary(res.pgls)[5]$coefficients[2,1], 
                   "p" = summary(res.pgls)[5]$coefficients[2,4], 
                   "t" = summary(res.pgls)[5]$coefficients[2,3], 
                   "aicc" = res.pgls$aicc,
                   "slope" = coef(res.pgls)[2],
                   "int" = coef(res.pgls)[1],
                   "sterr" = res.pgls$sterr)
  
  
  return(res_list)
}


#------------------------------------------------------------------------
# format df to include the avonet HWI data

avonet$Species2_eBird<-gsub(" ", "_", avonet$Species2_eBird)
avonet <- left_join(avonet, all_species_key, by = c("Species2_eBird" = "Species_avonet"))
avonet <- avonet %>% dplyr::filter(!is.na(Species_group))
avonet <- avonet %>% dplyr::select(Species2_eBird, Beak.Length_Culmen:Tail.Length, Species_group)

avonet <- avonet %>% 
  dplyr::group_by(Species_group) %>% 
  dplyr::summarise(dplyr::across(Beak.Length_Culmen:Tail.Length, list(~ mean(., na.rm = TRUE))))


# add avonet.birdlife, which includes diet data
avonet.bl$Species1<-gsub(" ", "_", avonet.bl$Species1)
all_species_key[62,1] <- "Islerothraupis_cristata"
all_species_key[63,1] <- "Islerothraupis_luctuosa"
all_species_key[64,1] <- "Leucippus_chlorocercus"
avonet.bl <- left_join(avonet.bl, all_species_key, by = c("Species1" = "Species_avonet"))
avonet.bl <- avonet.bl %>% dplyr::filter(!is.na(Species_group))
avonet.bl <- avonet.bl %>% 
  dplyr::filter(Species1 != "Ceratopipra_mentalis") %>% 
  dplyr::filter(Species1 != "Schiffornis_aenea") %>% 
  dplyr::filter(Species1 != "Schiffornis_olivacea") %>% 
  dplyr::filter(Species1 != "Schiffornis_stenorhyncha") %>% 
  dplyr::filter(Species1 != "Schiffornis_veraepacis")
avonet.bl <- avonet.bl %>% 
  dplyr::select(Species_group, Habitat:Primary.Lifestyle, Range.Size, Mass)

avonet.bl.2 <- avonet.bl %>% 
  dplyr::group_by(Species_group) %>% 
  dplyr::summarise(Mass.m = mean(Mass, na.rm = TRUE))
avonet.bl.3 <- avonet.bl %>% 
  dplyr::group_by(Species_group) %>% 
  dplyr::summarise(Range.Size.m = sum(Range.Size))
avonet.bl <- avonet.bl %>% 
  dplyr::distinct(Species_group, .keep_all= TRUE) %>% 
  dplyr::select(Species_group:Primary.Lifestyle)
avonet.bl <- left_join(avonet.bl, avonet.bl.2, by = "Species_group")
avonet.bl <- left_join(avonet.bl, avonet.bl.3, by = "Species_group")
avonet.bl <- avonet.bl %>%
  dplyr::rename(Mass = Mass.m) %>%
  dplyr::rename(Range.Size = Range.Size.m)

df <- left_join(df, avonet, by = c("species" = "Species_group"))
df <- left_join(df, avonet.bl, by = c("species" = "Species_group"))

df <- df %>% dplyr::select(-IBDa_genepop, -IBDe_genepop, -Fst_genepop, -reanalyzed)
df <- df %>% 
  dplyr::rename(Hand.wing.Index = Hand.wing.Index_1) %>%
  dplyr::rename(Beak.Length_Nares = Beak.Length_Nares_1)


# and PopGenome / genepop data

PopGen[48,1] <- "Leucippus_chlorocercus"
PopGen[57,1] <- "Thamnophilus_cryptoleucus"
df <- left_join(df, PopGen, by = "species")

genepop[16,1] <- "Thamnophilus_cryptoleucus"
df <- left_join(df, genepop, by = "species")

write.csv(df, file = "3_results/pgls.anova.database.formatted.csv", row.names = FALSE)
# df <- read.csv("3_results/pgls.anova.database.formatted.csv")

df["Trophic.Level"][df["Trophic.Level"] == "Carnivore"] <- 1
df["Trophic.Level"][df["Trophic.Level"] == "Herbivore"] <- 2
df["Trophic.Level"][df["Trophic.Level"] == "Omnivore"] <- 3
df$Trophic.Level <- as.integer(df$Trophic.Level)

df["Trophic.Niche"][df["Trophic.Niche"] == "Invertivore"] <- 1
df["Trophic.Niche"][df["Trophic.Niche"] == "Nectarivore"] <- 2
df["Trophic.Niche"][df["Trophic.Niche"] == "Omnivore"] <- 3
df["Trophic.Niche"][df["Trophic.Niche"] == "Frugivore"] <- 4
df$Trophic.Niche <- as.integer(df$Trophic.Niche)

df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Aerial"] <- 1
df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Terrestrial"] <- 2
df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Insessorial"] <- 3
df$Primary.Lifestyle <- as.integer(df$Primary.Lifestyle)

df["Habitat"][df["Habitat"] == "Forest"] <- 1
df["Habitat"][df["Habitat"] == "Woodland"] <- 2
df["Habitat"][df["Habitat"] == "Shrubland"] <- 3
df["Habitat"][df["Habitat"] == "Human Modified"] <- 4
df["Habitat"][df["Habitat"] == "Riverine"] <- 5
df$Habitat <- as.integer(df$Habitat)

# add in new Fst data, with individuals assigned to populations from DAPC

new.fst[48,1] <- "Leucippus_chlorocercus"
new.fst[57,1] <- "Thamnophilus_cryptoleucus"
new.fst <- new.fst %>% 
  dplyr::rename(Nei.G_ST.pop = Nei.G_ST) %>%
  dplyr::rename(nucleotide.F_ST.pop = nucleotide.F_ST) %>%
  dplyr::rename(Dxy.pop = Dxy)

new.fst <- new.fst[c("species", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop")]

df <- left_join(df, new.fst, by = "species")

# final trait dataset for upload
write.csv(df, file = "3_results/trait.database.formatted.final.csv", row.names = FALSE)
# df <- read.csv("3_results/trait.database.formatted.final.csv")

sister_clade <- read.csv("4_River_Island_tables_figures/Supplemental_data_table.csv")
sister_clade <- sister_clade %>% dplyr::select(species, Sister_clade_richness)
# sister_clade <- sister_clade %>% 
#   dplyr::mutate(Sister_clade_richness.bin = if_else(Sister_clade_richness == 1, 0, 1))
df <- left_join(df, sister_clade, by = "species")

het.all <- read.csv("1_analysis/heterozygosity_stats/allspecies_heterozygosity.formatted.csv")
df <- left_join(df, het.all, by = "species")


#------------------------------------------------------------------------
# Make new trait; kipp's/bill depth
# and theta/bp

df <- df %>% dplyr::mutate(Kipps.BillDepth = Kipps.Distance_1 / Beak.Depth_1)
df <- df %>% dplyr::mutate(theta.bp = theta / total_bp)
df <- df %>% dplyr::mutate(heterozygosity.bp = heterozygosity.nsites / total_bp)
write.csv(df, file = "3_results/trait.database.formatted.v2.final.csv", row.names = FALSE)
df <- read.csv("3_results/trait.database.formatted.v2.final.csv")
rownames(df) <- df$species

traits <- c("Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
           "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
           "D", "theta", "seg_sites", "pairwise_diffs", "nuc_div",
           "seg_sites_per_bp", "subtending_branch", "stem_length",
           "mtDNA.branch.length", "Beak.Length_Nares", "Hand.wing.Index",
           "Trophic.Niche", "Mass", "Range.Size", "nucleotide.F_ST",
           "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
           "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop", 
           "Kipps.BillDepth", "Kipps.Distance_1", "Beak.Depth_1", "Beak.Width_1",
           "Secondary1_1", "Wing.Length_1", "theta.bp", "heterozygosity.bp",
           "heterozygosity", "inbreeding", "N_SITES", "Sister_clade_richness")
#, "Sister_clade_richness.bin"

#rescale for effect sizes
df.temp <- df[,traits]
rownames(df.temp) <- df$species
df.temp.2 <- as.matrix(df.temp)
df.temp.2 <- scale(df.temp.2)
df.rescaled <- as.data.frame(df.temp.2)
df.rescaled$species <- df$species
df.rescaled$habitat <- df$habitat

df.anova <- data.frame(trait=character(),
                      F=double(), 
                      P=double(), 
                      island.floodplain.T=double(), 
                      island.floodplain.P=double(), 
                      island.upland.T=double(), 
                      island.upland.P=double(), 
                      floodplain.upland.T=double(), 
                      floodplain.upland.P=double(), 
                      stringsAsFactors=FALSE) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  
  colnum <- which(colnames(df.rescaled) == i)
  res <- phy_anova_island.char(df.rescaled, char=colnum)
  
  # colnum <- which(colnames(df.temp) == i)
  # res <- phy_anova_island.char(df.temp, char=colnum)
  
  anova.res <- as.data.frame(round(as.numeric(res$F), digits = 3))
  colnames(anova.res) <- "F"
  anova.res$P <- round(as.numeric(res$p), digits = 4)
  anova.res$island.floodplain.T <- round(as.numeric(res$island.floodplain.T), digits = 3) 
  anova.res$island.floodplain.P <- round(as.numeric(res$island.floodplain.P), digits = 3) 
  anova.res$island.upland.T <- round(as.numeric(res$island.upland.T), digits = 3) 
  anova.res$island.upland.P <- round(as.numeric(res$island.upland.P), digits = 3) 
  anova.res$floodplain.upland.T <- round(as.numeric(res$floodplain.upland.T), digits = 3) 
  anova.res$floodplain.upland.P <- round(as.numeric(res$floodplain.upland.P), digits = 3) 
  anova.res$trait <- i
  anova.res <- anova.res[,c(9,1:8)]
  df.anova <- bind_rows(df.anova, anova.res)
  
}


fname <- paste0("3_results/phylANOVA.results.v9.csv")
write.csv(df.anova, file = fname, row.names = FALSE)
df.anova <- read.csv(fname)

# fname <- paste0("3_results/phylANOVA.results.v5.AllTraits.csv")
# write.csv(df.anova, file = fname, row.names = FALSE)

#------------------------------------------------------------------------
# run with Ochthornis and Dendroplex reassigned

df.temp <- df.rescaled
df.temp["Ochthornis_littoralis","habitat"] <- "floodplain"
df.temp["Dendroplex_kienerii","habitat"] <- "floodplain"

df.anova <- data.frame(trait=character(),
                       F=double(), 
                       P=double(), 
                       island.floodplain.T=double(), 
                       island.floodplain.P=double(), 
                       island.upland.T=double(), 
                       island.upland.P=double(), 
                       floodplain.upland.T=double(), 
                       floodplain.upland.P=double(), 
                       stringsAsFactors=FALSE) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.temp) == i)
  
  res <- phy_anova_island.char(df.temp, char=colnum)
  anova.res <- as.data.frame(round(as.numeric(res$F), digits = 3))
  colnames(anova.res) <- "F"
  anova.res$P <- round(as.numeric(res$p), digits = 4)
  anova.res$island.floodplain.T <- round(as.numeric(res$island.floodplain.T), digits = 3) 
  anova.res$island.floodplain.P <- round(as.numeric(res$island.floodplain.P), digits = 3) 
  anova.res$island.upland.T <- round(as.numeric(res$island.upland.T), digits = 3) 
  anova.res$island.upland.P <- round(as.numeric(res$island.upland.P), digits = 3) 
  anova.res$floodplain.upland.T <- round(as.numeric(res$floodplain.upland.T), digits = 3) 
  anova.res$floodplain.upland.P <- round(as.numeric(res$floodplain.upland.P), digits = 3) 
  anova.res$trait <- i
  anova.res <- anova.res[,c(9,1:8)]
  df.anova <- bind_rows(df.anova, anova.res)
  
}


fname <- paste0("3_results/phylANOVA.results.v6.Ochthornis.Dendroplex.csv")
write.csv(df.anova, file = fname, row.names = FALSE)
df.anova <- read.csv(fname)


#------------------------------------------------------------------------
# run anova with island split into early and late successional stages
# note that a different function is necessary here: phy_anova_island.islandsplit


df.temp <- df.rescaled
df.temp["habitat"][df.temp["habitat"] == "island"] <- "island.early"

df.temp["Dendroplex_kienerii","habitat"] <- "island.late"
df.temp["Myrmotherula_assimilis","habitat"] <- "island.late"
df.temp["Myrmotherula_klagesi","habitat"] <- "island.late"
df.temp["Thamnophilus_cryptoleucus","habitat"] <- "island.late"
df.temp["Myrmoborus_lugubris","habitat"] <- "island.late"
df.temp["Elaenia_pelzelni","habitat"] <- "island.late"
df.temp["Conirostrum_bicolor","habitat"] <- "island.late"
df.temp["Conirostrum_margaritae","habitat"] <- "island.late"

df.anova <- data.frame(trait=character(),
                       F=double(), 
                       P=double(), 
                       island.early.island.late.T=double(), 
                       island.early.island.late.P=double(), 
                       island.early.floodplain.T=double(), 
                       island.early.floodplain.P=double(), 
                       island.early.upland.T=double(), 
                       island.early.upland.P=double(), 
                       island.late.floodplain.T=double(), 
                       island.late.floodplain.P=double(), 
                       island.late.upland.T=double(), 
                       island.late.upland.P=double(), 
                       floodplain.upland.T=double(), 
                       floodplain.upland.P=double(), 
                       stringsAsFactors=FALSE) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.temp) == i)
  
  res <- phy_anova_island.islandsplit(df.temp, char=colnum)
  anova.res <- as.data.frame(round(as.numeric(res$F), digits = 3))
  colnames(anova.res) <- "F"
  anova.res$P <- round(as.numeric(res$p), digits = 4)
  anova.res$island.early.island.late.T <- round(as.numeric(res$island.early.island.late.T), digits = 3) 
  anova.res$island.early.island.late.P <- round(as.numeric(res$island.early.island.late.P), digits = 3) 
  anova.res$island.early.floodplain.T <- round(as.numeric(res$island.early.floodplain.T), digits = 3) 
  anova.res$island.early.floodplain.P <- round(as.numeric(res$island.early.floodplain.P), digits = 3) 
  anova.res$island.early.upland.T <- round(as.numeric(res$island.early.upland.T), digits = 3) 
  anova.res$island.early.upland.P <- round(as.numeric(res$island.early.upland.P), digits = 3) 
  anova.res$island.late.floodplain.T <- round(as.numeric(res$island.late.floodplain.T), digits = 3) 
  anova.res$island.late.floodplain.P <- round(as.numeric(res$island.late.floodplain.P), digits = 3) 
  anova.res$island.late.upland.T <- round(as.numeric(res$island.late.upland.T), digits = 3) 
  anova.res$island.late.upland.P <- round(as.numeric(res$island.late.upland.P), digits = 3) 
  anova.res$floodplain.upland.T <- round(as.numeric(res$floodplain.upland.T), digits = 3) 
  anova.res$floodplain.upland.P <- round(as.numeric(res$floodplain.upland.P), digits = 3) 
  anova.res$trait <- i
  anova.res <- anova.res[,c(15,1:14)]
  df.anova <- bind_rows(df.anova, anova.res)
  
}


fname <- paste0("3_results/phylANOVA.results.v6.islandsplit.csv")
write.csv(df.anova, file = fname, row.names = FALSE)
df.anova <- read.csv(fname)


#------------------------------------------------------------------------

# AICC Weights

df.temp <- df.temp %>% 
  dplyr::filter(species != "Leucippus_chlorocercus") %>% 
  dplyr::filter(species != "Elaenia_pelzelni")

traits.gen <- traits[c(1:2,4,7:14,17,23:24,26)]
anova.l <- list() 

for (x in 1:length(traits.gen)) {
  i <- traits.gen[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.temp) == i)
  
  res <- phy_anova_island.model(df.temp, char=colnum)
  res.l <- list(res)
  anova.l <- c(anova.l, res.l)
  
  # res.l <- list(res.island.aov)
  # anova.l <- c(anova.l, (res.island.aov))
  
}

res.iffy <- aictab(cand.set = anova.l, modnames = traits.gen)
res.iffy %>% arrange(desc(Cum.Wt))

#------------------------------------------------------------------------
# run HWI separately for species with HWI < & > 30

df.l <- df.rescaled %>% dplyr::filter(HWI <= 2) #units in effect size
df.h <- df.rescaled %>% dplyr::filter(HWI > 2)
i <- "Hand.wing.Index"
colnum <- which(colnames(df.l) == i)
res <- phy_anova_island.char(df.l, char=colnum)
res$anova

colnum <- which(colnames(df.h) == i)
res <- phy_anova_island.char(df.h, char=colnum)
res$anova


# with passerines vs non-passerines
df.l <- df.rescaled[c("Campephilus_melanoleucos", "Campephilus_rubricollis", "Celeus_flavus", "Celeus_grammicus",
                      "Crypturellus_undulatus", "Crypturellus_variegatus", "Glaucidium_brasilianum", "Glaucidium_hardyi",
                      "Leucippus_chlorocercus", "Megascops_choliba", "Megascops_watsonii", "Monasa_morphoeus",
                      "Monasa_nigrifrons", "Phaethornis_bourcieri", "Phaethornis_hispidus", "Piaya_cayana", "Piaya_melanogaster",
                      "Trogon_collaris", "Trogon_rufus"),]
                    
df.h <- df.rescaled[c("Cantorchilus_leucotis", "Conirostrum_bicolor", "Conirostrum_margaritae", "Cranioleuca_vulpecula",
                      "Dendroplex_kienerii", "Elaenia_pelzelni", "Formicarius_analis", "Formicarius_colma", "Furnarius_minor",
                      "Hylophylax_naevia", "Hylophylax_punctulata", "Knipolegus_orenocensis", "Mazaria_propinqua",
                      "Myrmeciza_fortis", "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_lugubris", "Myrmoborus_myotherinus",
                      "Myrmochanes_hemileucus", "Myrmotherula_assimilis", "Myrmotherula_klagesi", "Ochthornis_littoralis",
                      "Pheugopedius_coraya", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens",
                      "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "Serpophaga_hypoleuca",
                      "Stigmatura_napensis", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus",
                      "Tachyphonus_luctuosus", "Thamnophilus_cryptoleucus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus"),]

i <- "Hand.wing.Index"
colnum <- which(colnames(df.l) == i)
res <- phy_anova_island.char(df.l, char=colnum)
res$anova

colnum <- which(colnames(df.h) == i)
res <- phy_anova_island.char(df.h, char=colnum)
res$anova

# no significance across habitats


#------------------------------------------------------------------------

Av_groups <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Av_groups)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nNumber of groups (average)") 

Dxy <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Dxy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nDxy") 

Fst <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Nei.G_ST)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") 

SNPs <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nTotal SNPs") 

SNPs_per_locus <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs_per_locus)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSNPs per locus") 

SNPs_per_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs_per_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSNPs per bp") 

loci <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = loci)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nNumber of loci") 

av_contig_length <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = av_contig_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage contig length") 

total_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = total_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nTotal bp") 

Av_gene_tree <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Av_UCE_gene_tree_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nAverage UCE gene tree length") 

Taj_D <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = D)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nTajima's D") 

theta <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = theta)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nTheta") 


theta.bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = theta.bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nTheta/bp") 

seg_sites <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = seg_sites)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSegregating sites") 

seg_sites_per_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = seg_sites_per_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSegregating sites per bp") 

pairwise_diffs <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = pairwise_diffs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage pairwise differences") 

nuc_div <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = nuc_div)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nNucleotide diversity") 

Av_mtdna_tree <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = mtDNA.branch.length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage mtDNA tree length") 


IBD <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = e_slope)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nIBD (slope)") 



heterozygosity <- df %>% dplyr::filter(species != "Elaenia_pelzelni") %>% 
  # dplyr::filter(heterozygosity < 0.3) %>% # removing three species with low sample size, which seems to strongly affect heterozygosity estimates 
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = heterozygosity)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHeterozygosity") 

ggsave("3_results/1_plots/heterozygosity_boxplot.pdf", 
       width = 6, height = 6, units = "in")


heterozygosity.bp <- df %>% dplyr::filter(species != "Elaenia_pelzelni") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = heterozygosity.bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHeterozygosity / bp") 

ggsave("3_results/1_plots/heterozygosity.bp_boxplot.pdf", 
       width = 6, height = 6, units = "in")



inbreeding <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = inbreeding)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nInbreeding coefficient") 




# subtending branch

stem <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = stem_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nStem branch length") 

subtending <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = subtending_branch)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSubtending branch length") 


Sister_clade_richness <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log10(Sister_clade_richness))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSister clade richness (log)") 

ggsave("3_results/1_plots/Sister_clade_richness_boxplot.pdf", 
       width = 6, height = 6, units = "in")



# there's a decreasing trend in the outliers, but this is kind of reaching for a pattern
# removing the extreme outlier (Leucippus; 40) dampens the pattern
df %>% dplyr::filter(Sister_clade_richness >= 3) %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log10(Sister_clade_richness))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSister clade richness (log)") 


# traits 

hwi <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand-wing Index") 

Range.Size <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log(Range.Size))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nRange size (log)") 

mass <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log(Mass))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nMass (log)") 

Beak.Length_Nares <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Length_Nares)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak length at nares") 


Trophic.Niche <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = scale(Trophic.Niche))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nDiet categories") 


hwi.30 <- df %>% filter(Hand.wing.Index < 30) %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand-wing Index (HWI < 30)") 



rownames(df) <- df$species
df.passerine <- df[c("Cantorchilus_leucotis", "Conirostrum_bicolor", "Conirostrum_margaritae", "Cranioleuca_vulpecula",
                     "Dendroplex_kienerii", "Elaenia_pelzelni", "Formicarius_analis", "Formicarius_colma", "Furnarius_minor",
                     "Hylophylax_naevia", "Hylophylax_punctulata", "Knipolegus_orenocensis", "Mazaria_propinqua",
                     "Myrmeciza_fortis", "Myrmeciza_hyperythra", "Myrmoborus_leucophrys", "Myrmoborus_lugubris", "Myrmoborus_myotherinus",
                     "Myrmochanes_hemileucus", "Myrmotherula_assimilis", "Myrmotherula_klagesi", "Ochthornis_littoralis",
                     "Pheugopedius_coraya", "Pipra_erythrocephala", "Pipra_filicauda", "Saltator_coerulescens",
                     "Saltator_grossus", "Schiffornis_major", "Schiffornis_turdina", "Serpophaga_hypoleuca",
                     "Stigmatura_napensis", "Synallaxis_gujanensis", "Synallaxis_rutilans", "Tachyphonus_cristatus",
                     "Tachyphonus_luctuosus", "Thamnophilus_cryptoleucus", "Xiphorhynchus_elegans", "Xiphorhynchus_obsoletus"),]

hwi.passerine <- df.passerine %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand-wing Index (passerines)") 


# plots for a few morphometric traits, not used in paper


df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Kipps.BillDepth)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nKipps.BillDepth") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Kipps.Distance_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nKipps.Distance") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Depth_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak.Depth") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Width_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak.Width") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Secondary1_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSecondary1") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length.Mass)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length / Mass") 


df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length.Mass.residuals)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length / Mass, residuals") 

df %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index.Mass)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand.wing.Index / Mass") 

df %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index.Mass.residuals)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand.wing.Index / Mass, residuals") 



#------------------------------------------------------------------------

# ggarrange(Fst, Dxy, nuc_div, Taj_D,
#           log.theta, pairwise_diffs, seg_sites, seg_sites_per_bp,
#           loci, total_bp, av_contig_length, Av_groups, Av_gene_tree,
#           Av_gene_tree, SNPs, SNPs_per_locus, SNPs_per_bp, 
#           ncol = 4, nrow = 4                                   
# ) 
# text.p <- ggparagraph(text = "Habitat", size = 11, color = "black")
# 
# 
# grid.arrange(Fst, Dxy, nuc_div, Taj_D,
#              log.theta, pairwise_diffs, seg_sites, seg_sites_per_bp,
#              loci, total_bp, av_contig_length, Av_groups, Av_gene_tree,
#              SNPs, SNPs_per_locus, SNPs_per_bp,
#              ncol = 4)   

gt <- arrangeGrob(heterozygosity, SNPs_per_locus, Av_mtdna_tree, Fst, 
                  Av_gene_tree, IBD, Dxy, inbreeding,
                  SNPs_per_bp, nuc_div, seg_sites, theta.bp, 
                  pairwise_diffs, Av_groups, Taj_D, Range.Size,
                  nrow = 4, ncol = 4)

as_ggplot(gt) +
  draw_plot_label(label = c("A*", "B*", "C*", "D*",
                            "E*", "F*", "G*", "H*",
                            "I*", "J*", "K*", "L*",
                            "M", "N", "O", "P"), 
                  size = 15,
                  x = c(0.001, 0.25, 0.5, 0.75,
                        0.001, 0.25, 0.5, 0.75,
                        0.001, 0.25, 0.5, 0.75,
                        0.001, 0.25, 0.5, 0.75), 
                  y = c(1,    1,    1,    1,
                        0.75, 0.75, 0.75, 0.75,
                        0.5, 0.5, 0.5, 0.5,
                        0.25, 0.25, 0.25, 0.25))

ggsave("3_results/1_plots/multipage_boxplot.v8.pdf", 
       width = 12, height = 12, units = "in")


gt.sub <- arrangeGrob(subtending, stem, Sister_clade_richness,
                      nrow = 1, ncol = 3)
as_ggplot(gt.sub) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 15,
                  x = c(0.001, 0.33, 0.66), 
                  y = c(1,  1,  1))
ggsave("3_results/1_plots/multipage_boxplot.v7.subtending.pdf", 
       width = 9, height = 3, units = "in")

# traits
gt.trait <- arrangeGrob(hwi, hwi.30, hwi.passerine, 
                        mass, Beak.Length_Nares, Trophic.Niche,
                        nrow = 2, ncol = 3)
as_ggplot(gt.trait) +
  draw_plot_label(label = c("A", "B", "C",
                            "D", "E", "F"), 
                  size = 15,
                  x = c(0.001, 0.33, 0.65,
                        0.001, 0.33, 0.65), 
                  y = c(1,     1,    1,
                        0.5,   0.5,  0.5))
ggsave("3_results/1_plots/multipage_boxplot.v7.traits.png", 
       width = 9, height = 6, units = "in")

# descriptive stats
gt.descriptive <- arrangeGrob(total_bp, av_contig_length, 
                              loci,
                        nrow = 2, ncol = 2)
as_ggplot(gt.descriptive) +
  draw_plot_label(label = c("A*", "B*",
                            "C", ""), 
                  size = 15,
                  x = c(0.001, 0.5,
                        0.001, 0.5), 
                  y = c(1,     1,
                        0.5,   0.5))
ggsave("3_results/1_plots/multipage_boxplot.v6.descriptive.png", 
       width = 6, height = 6, units = "in")

#------------------------------------------------------------------------
# pgls

traits <- c("Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
            "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
            "D", "theta", "seg_sites", "pairwise_diffs", "nuc_div",
            "seg_sites_per_bp", "subtending_branch", "stem_length",
            "mtDNA.branch.length", "Beak.Length_Nares", "Hand.wing.Index",
            "Trophic.Niche", "Mass", "Range.Size", "nucleotide.F_ST",
            "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
            "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop", 
            "Kipps.BillDepth", "Kipps.Distance_1", "Beak.Depth_1", "Beak.Width_1",
            "Secondary1_1", "Wing.Length_1", "theta.bp", "heterozygosity.bp",
             "inbreeding", "heterozygosity", "N_SITES", "Sister_clade_richness")

df.o <- df.rescaled 
# df.o <- df

# traits <- colnames(df)[c(7, 11:25, 29, 31, 38, 44, 46:47, 49:50, 54:61)]

df.pgls <- data.frame(trait=character(),
                      estimate=double(), 
                      p=double(), 
                      t=double(), 
                      aicc=double(), 
                      slope=double(), 
                      int=double(), 
                      sterr.int=double(), 
                      sterr.slope=double() 
                      ) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.o) == i)
  
  res <- pgls.char(df.o, char=colnum)
  pgls.res <- as.data.frame(round(as.numeric(res$estimate), digits = 3))
  colnames(pgls.res) <- "estimate"
  pgls.res$p <- round(as.numeric(res$p), digits = 4)
  pgls.res$t <- round(as.numeric(res$t), digits = 3) 
  pgls.res$aicc <- round(as.numeric(res$aicc[1,1]), digits = 3) 
  pgls.res$slope <- round(as.numeric(res$slope), digits = 3) 
  pgls.res$int <- round(as.numeric(res$int), digits = 3) 
  pgls.res$sterr.int <- round(as.numeric(res$sterr[1]), digits = 3) 
  pgls.res$sterr.slope <- round(as.numeric(res$sterr[2]), digits = 3) 
  pgls.res$trait <- i
  pgls.res <- pgls.res[,c(9,1:8)]
  df.pgls <- bind_rows(df.pgls, pgls.res)
  
}

fname <- paste0("3_results/pgls.results.v5.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)
df.pgls <- read.csv(fname)

df.pgls.1 <- df.pgls %>%
  dplyr::filter(trait != "log.theta") %>%
  dplyr::filter(trait != "log.total_bp") %>%
  dplyr::filter(trait != "loci") %>%
  dplyr::filter(trait != "av_contig_length") %>%
  dplyr::filter(trait != "total_bp") %>%
  dplyr::filter(trait != "HWI")

df.pgls.1$deltaAIC <- akaike.weights(df.pgls.1$aicc)$deltaAIC
df.pgls.1$rel.LL <- akaike.weights(df.pgls.1$aicc)$rel.LL
df.pgls.1$weights <- akaike.weights(df.pgls.1$aicc)$weights

fname <- paste0("3_results/pgls.results.aicw.v3.csv")
write.csv(df.pgls.1, file = fname, row.names = FALSE)

# Remove more traits so that AIC weights are for traits with the same sample sizes
# and separate for each group of traits

# genetics
df.pgls.2 <- df.pgls.1 %>%
  dplyr::filter(trait != "subtending_branch") %>%
  dplyr::filter(trait != "stem_length") %>%
  dplyr::filter(trait != "Beak.Length_Nares") %>%
  dplyr::filter(trait != "Hand.wing.Index") %>%
  dplyr::filter(trait != "Mass") %>%
  dplyr::filter(trait != "Range.Size") %>%
  dplyr::filter(trait != "nucleotide.F_ST") %>%
  dplyr::filter(trait != "e_statistic") %>%
  dplyr::filter(trait != "a_statistic") %>%
  dplyr::filter(trait != "a_slope")


df.pgls.2$deltaAIC <- akaike.weights(df.pgls.2$aicc)$deltaAIC
df.pgls.2$rel.LL <- akaike.weights(df.pgls.2$aicc)$rel.LL
df.pgls.2$weights <- akaike.weights(df.pgls.2$aicc)$weights

fname <- paste0("3_results/pgls.results.aicw.v3.genetics.csv")
write.csv(df.pgls.2, file = fname, row.names = FALSE)


# traits

df.pgls.3 <- df.pgls.1 %>%
  dplyr::filter(trait == "Beak.Length_Nares" | trait == "Hand.wing.Index" | trait == "Mass" | trait == "Range.Size" | trait == "Trophic.Niche")

df.pgls.3$deltaAIC <- akaike.weights(df.pgls.3$aicc)$deltaAIC
df.pgls.3$rel.LL <- akaike.weights(df.pgls.3$aicc)$rel.LL
df.pgls.3$weights <- akaike.weights(df.pgls.3$aicc)$weights

fname <- paste0("3_results/pgls.results.aicw.v3.traits.csv")
write.csv(df.pgls.3, file = fname, row.names = FALSE)


# subtending branch length

df.pgls.4 <- df.pgls.1 %>%
  dplyr::filter(trait == "subtending_branch" | trait == "stem_length")

df.pgls.4$deltaAIC <- akaike.weights(df.pgls.4$aicc)$deltaAIC
df.pgls.4$rel.LL <- akaike.weights(df.pgls.4$aicc)$rel.LL
df.pgls.4$weights <- akaike.weights(df.pgls.4$aicc)$weights

fname <- paste0("3_results/pgls.results.aicw.v3.subtending_branch.csv")
write.csv(df.pgls.4, file = fname, row.names = FALSE)

#------------------------------------------------------------------------
# run with Ochthornis and Dendroplex reassigned

df.o["Ochthornis_littoralis","habitat"] <- "floodplain"
df.o["Dendroplex_kienerii","habitat"] <- "floodplain"


df.pgls <- data.frame(trait=character(),
                      estimate=double(), 
                      p=double(), 
                      t=double(), 
                      aicc=double(), 
                      slope=double(), 
                      int=double(), 
                      sterr.int=double(), 
                      sterr.slope=double() 
) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.o) == i)
  
  res <- pgls.char(df.o, char=colnum)
  pgls.res <- as.data.frame(round(as.numeric(res$estimate), digits = 3))
  colnames(pgls.res) <- "estimate"
  pgls.res$p <- round(as.numeric(res$p), digits = 4)
  pgls.res$t <- round(as.numeric(res$t), digits = 3) 
  pgls.res$aicc <- round(as.numeric(res$aicc[1,1]), digits = 3) 
  pgls.res$slope <- round(as.numeric(res$slope), digits = 3) 
  pgls.res$int <- round(as.numeric(res$int), digits = 3) 
  pgls.res$sterr.int <- round(as.numeric(res$sterr[1]), digits = 3) 
  pgls.res$sterr.slope <- round(as.numeric(res$sterr[2]), digits = 3) 
  pgls.res$trait <- i
  pgls.res <- pgls.res[,c(9,1:8)]
  df.pgls <- bind_rows(df.pgls, pgls.res)
  
}

fname <- paste0("3_results/pgls.results.v4.Ochthornis.Dendroplex.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)



#------------------------------------------------------------------------
# plot Fst and HWI separately
fname <- paste0("3_results/pgls.results.v4.csv")
df.pgls <- read.csv(fname)
rownames(df.pgls) <- df.pgls$trait

label.Fst <- paste0("p < 0.0001, t = ", prettyNum(df.pgls["Nei.G_ST","t"], digits = 3), 
                       "\nestimate = ", prettyNum(df.pgls["Nei.G_ST","estimate"], digits = 3))


df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Nei.G_ST)) + 
  geom_abline(slope = ,
              intercept = df.pgls[3,7]) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") +
  annotate("text", x=1.2, y=0.55, label = label.Fst, size = 3.6)

ggsave("3_results/1_plots/Fst.v4.pdf", 
       width = 6,
       height = 6,
       units = c("in"))


label.HWI <- paste0("p = ", prettyNum(df.pgls["Hand.wing.Index","p"], digits = 2),
                    ", t = ", prettyNum(df.pgls["Hand.wing.Index","t"], digits = 3), 
                    "\nestimate = ", prettyNum(df.pgls["Hand.wing.Index","estimate"], digits = 3))

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log10(Hand.wing.Index))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  scale_y_continuous(limits = c(0.75, 2.25),
                     breaks = c(0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25)) +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nHand-wing index (log)") #+
  # annotate("text", x=1, y=4.4, label = label.HWI, size = 3.6)

ggsave("3_results/1_plots/HWI.v5.png", width = 6,
       height = 6,
       units = c("in"))


# removing hummingbirds
df %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log10(Hand.wing.Index))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nHand-wing index (log)") 

ggsave("3_results/1_plots/HWI.v4.pdf", width = 6,
       height = 6,
       units = c("in"))

#------------------------------------------------------------------------
# plot Fst and HWI separately
# but reverse order of habitats on x axis for defense slides

fname <- paste0("3_results/pgls.results.v3.csv")
df.pgls <- read.csv(fname)

label.Fst <- paste0("PGLS\np < 0.001, t = ", prettyNum(df.pgls[24,4], digits = 3), 
                    "\nestimate = -", prettyNum(df.pgls[24,2], digits = 3))


df %>%
  mutate(habitat = fct_relevel(habitat, "upland", "floodplain", "island")) %>%
  ggplot( aes(x = habitat, y = Nei.G_ST)) + 
  geom_abline(slope = ,
              intercept = df.pgls[3,7]) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") +
  annotate("text", x=3, y=0.55, label = label.Fst, size = 3.6)

ggsave("3_results/1_plots/Fst.for_talk.v1.pdf")
  

label.HWI <- paste0("PGLS\np = ", prettyNum(df.pgls[19,3], digits = 2),
                    ", t = ", prettyNum(df.pgls[19,4], digits = 3), 
                    "\nestimate = -", prettyNum(df.pgls[19,2], digits = 3))

df %>%
  mutate(habitat = fct_relevel(habitat, "upland", "floodplain", "island")) %>%
  ggplot( aes(x = habitat, y = log(Hand.wing.Index))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  scale_y_continuous(limits=c(1.8, 4.6), 
                     breaks = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nHand-wing index (log)") +
  annotate("text", x=3, y=4.5, label = label.HWI, size = 3.6) 

ggsave("3_results/1_plots/HWI.for_talk.v1.pdf")

  
  # df.o %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  #   ggplot( aes(x = Hand.wing.Index, y = nucleotide.F_ST.pop)) + 
  #   geom_point(aes(colour = factor(habitat))) 
    
    
#------------------------------------------------------------------------
# pgls
# HWI vs genetic traits


df.o <- df

traits <- c("Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
              "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
              "theta", "seg_sites", "pairwise_diffs", "nuc_div",
              "seg_sites_per_bp", "subtending_branch", "stem_length",
              "mtDNA.branch.length", "Range.Size", "nucleotide.F_ST",
              "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
              "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop")
  
df.pgls <- data.frame(trait=character(),
                      estimate=double(), 
                      p=double(), 
                      t=double(), 
                      aicc=double(), 
                      slope=double(), 
                      int=double(), 
                      sterr.int=double(), 
                      sterr.slope=double() 
) 

for (x in 1:length(traits)) {
  i <- traits[x]
  
  # call trait function for each trait
  print(paste0("working on trait: ", i, collapse = ", "))
  cat("\n")
  colnum <- which(colnames(df.o) == i)
  
  res <- pgls.hwi(df.o, char=colnum)
  pgls.res <- as.data.frame(round(as.numeric(res$estimate), digits = 3))
  colnames(pgls.res) <- "estimate"
  pgls.res$p <- round(as.numeric(res$p), digits = 4)
  pgls.res$t <- round(as.numeric(res$t), digits = 3) 
  pgls.res$aicc <- round(as.numeric(res$aicc[1,1]), digits = 3) 
  pgls.res$slope <- round(as.numeric(res$slope), digits = 3) 
  pgls.res$int <- round(as.numeric(res$int), digits = 3) 
  pgls.res$sterr.int <- round(as.numeric(res$sterr[1]), digits = 3) 
  pgls.res$sterr.slope <- round(as.numeric(res$sterr[2]), digits = 3) 
  pgls.res$trait <- i
  pgls.res <- pgls.res[,c(9,1:8)]
  df.pgls <- bind_rows(df.pgls, pgls.res)
  
}

fname <- paste0("3_results/pgls.results.HWIvsGenetics.v1.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)


df.pgls <- df.pgls %>%
  dplyr::filter(trait != "log.theta") %>%
  dplyr::filter(trait != "log.total_bp") %>%
  dplyr::filter(trait != "loci") %>%
  dplyr::filter(trait != "av_contig_length") %>%
  dplyr::filter(trait != "total_bp") %>%
  dplyr::filter(trait != "HWI")

df.pgls$deltaAIC <- akaike.weights(df.pgls$aicc)$deltaAIC
df.pgls$rel.LL <- akaike.weights(df.pgls$aicc)$rel.LL
df.pgls.1$weights <- akaike.weights(df.pgls.1$aicc)$weights

fname <- paste0("3_results/pgls.results.HWIvsGenetics.v2.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)

#------------------------------------------------------------------------
# sister clade richness sign test (per Robb's request)
df.sister_clade <- df
df.sister_clade <- df.sister_clade %>% dplyr::select(species, habitat, Sister_clade_richness)
df.sister_clade <- df.sister_clade %>% 
  dplyr::mutate(Sister_clade_richness.bin = if_else(Sister_clade_richness == 1, 0, 1))
df.sister_clade <- df.sister_clade %>% filter(!is.na(Sister_clade_richness))



####
df.sister_clade.r <- df.sister_clade %>% filter(habitat == "island")
df.sister_clade.f <- df.sister_clade %>% filter(habitat == "floodplain")
df.sister_clade.u <- df.sister_clade %>% filter(habitat == "upland")

df.sister_clade.r$ID <- c(1:nrow(df.sister_clade.r))
df.sister_clade.f$ID <- c(1:nrow(df.sister_clade.f))
df.sister_clade.u$ID <- c(1:nrow(df.sister_clade.u))

df.sister_clade.rf <- left_join(df.sister_clade.r[,c(3,5)], df.sister_clade.f[,c(3,5)], by = "ID")
df.sister_clade.rf <- df.sister_clade.rf %>% dplyr::select(Sister_clade_richness.x, Sister_clade_richness.y)

df.sister_clade.rf$Sister_clade_richness.x <- sample(df.sister_clade.rf$Sister_clade_richness.x)
df.sister_clade.rf$Sister_clade_richness.y <- sample(df.sister_clade.rf$Sister_clade_richness.y)

# not sure that this is an appropriate test
diversity.contrast.test(df.sister_clade.rf)

####

# count N of 1's and N total
df.sc.t <- df.sister_clade %>% dplyr::group_by(habitat) %>% count(Sister_clade_richness.bin)
df.sc.n <- df.sister_clade %>% dplyr::group_by(habitat) %>% count()

binom.test(df.sc.t[4,3][[1]], df.sc.n[2,2][[1]]) # island
binom.test(df.sc.t[2,3][[1]], df.sc.n[1,2][[1]]) # floodplain
binom.test(df.sc.t[6,3][[1]], df.sc.n[3,2][[1]]) # upland


#------------------------------------------------------------------------
quit()

