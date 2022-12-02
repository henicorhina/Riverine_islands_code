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
library(car)
library(here)


# old file paths, when using external drive
# setwd("/Volumes/Brumfield_Lab_Drive/River_islands")
# theta <- read.csv("3_results/dendropy/theta.csv")
# island.Tree <- read.tree(file = "all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.newick.phy")
# island.Tree <- read.tree(file = "/Volumes/Brumfield_Lab_Drive/River_islands/all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.phy")
# island.Tree <- read.tree(file = "/Volumes/Brumfield_Lab_Drive/River_islands/all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.tre")
# avonet <- read.csv("/Volumes/Backup_Plus/AVONET/Supplementary_dataset_1_raw_data.csv") # doesn't include trophic data / diet
# avonet <- read.csv("/Volumes/Backup_Plus/AVONET/Supplementary_dataset_1.csv") # preliminary df
# avonet <- read.csv("/Users/oscar/Documents/Projects/AVONET/ELEData/TraitData/AVONET_Raw_Data.csv") # final df, raw data
# avonet.bl <- read.csv("/Users/oscar/Documents/Projects/AVONET/ELEData/TraitData/AVONET1_BirdLife.csv") # final df, birdlife
# genepop <- read.csv("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/genepop_ibd.results.v1.csv")

# relative file paths, with here()
theta <- read.csv("2_data/theta.csv")
df <- read.csv("2_data/results_for_pgls_take4.formatted.csv")
island.Tree <- read.tree(file = "2_data/RAxML_bipartitions.all_species_95percent_final.newick.phy")
island.Tree <- drop.tip(island.Tree, "Thamnophilus_cryptoleucus2") # duplicate tip
avonet <- read.csv("2_data/AVONET_Raw_Data.csv") # final df, raw data
avonet.bl <- read.csv("2_data/AVONET1_BirdLife.csv") # final df, birdlife
all_species_key <- read.csv("2_data/All_species_key.csv")
PopGen <- read.csv("2_data/PopGenome_Fst.results.v1.csv")
new.fst <- read.csv("2_data/PopGenome_Fst.results.v1.DAPCassignments.csv")
genepop <- read.csv("2_data/genepop_ibd.results.v1.csv")
het <- read.csv("2_data/Supplemental_data_table.csv")

#------------------------------------------------------------------------
# functions

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
    dplyr::select(HWI, species, char) %>% 
    dplyr::filter(!is.na(HWI)) %>%
    dplyr::filter(!is.na(char)) %>%
    droplevels()
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  hab.dat <- comparative.data(island.Tree.t, temp, species)
  
  if (char == 25) { # absolute value for log-transform Tajima's D. all are negative
    res.pgls <- pgls(log(abs(char)) ~ HWI, hab.dat)
  } else {
    res.pgls <- pgls(log(char) ~ HWI, hab.dat)
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

# more carnivorous
df["Trophic.Level"][df["Trophic.Level"] == "Herbivore"] <- 1
df["Trophic.Level"][df["Trophic.Level"] == "Omnivore"] <- 2
df["Trophic.Level"][df["Trophic.Level"] == "Carnivore"] <- 3
df$Trophic.Level <- as.integer(df$Trophic.Level)

# more carnivorous (invertivorous?)
df["Trophic.Niche"][df["Trophic.Niche"] == "Frugivore"] <- 1
df["Trophic.Niche"][df["Trophic.Niche"] == "Nectarivore"] <- 2
df["Trophic.Niche"][df["Trophic.Niche"] == "Omnivore"] <- 3
df["Trophic.Niche"][df["Trophic.Niche"] == "Invertivore"] <- 4
df$Trophic.Niche <- as.integer(df$Trophic.Niche)

# more terrestrial
df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Aerial"] <- 1
df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Terrestrial"] <- 2
df["Primary.Lifestyle"][df["Primary.Lifestyle"] == "Insessorial"] <- 3
df$Primary.Lifestyle <- as.integer(df$Primary.Lifestyle)

# more riverine? I don't think we used this anyways
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
df <- left_join(df, het[,c(1,25,59:60)], by = "species")

# some new variables
df <- df %>% dplyr::mutate(Kipps.BillDepth = Kipps.Distance_1 / Beak.Depth_1)
df <- df %>% dplyr::mutate(theta.bp = theta / total_bp)
df <- df %>% dplyr::mutate(heterozygosity.bp = heterozygosity / total_bp)

# final trait dataset for upload
write.csv(df, file = "3_results/trait.database.formatted.final.csv", row.names = FALSE)
# df <- read.csv("3_results/trait.database.formatted.final.csv")

#------------------------------------------------------------------------
# calculate variance inflation factor (GVIF) 
# discard values greater than 5

traits <- c("species", "Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
            "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
            "D", "theta", "theta.bp", "seg_sites", "pairwise_diffs", "nuc_div",
            "seg_sites_per_bp", "subtending_branch", "stem_length",
            "mtDNA.branch.length", "Beak.Length_Nares", "Hand.wing.Index",
            "Trophic.Niche", "Mass", "Range.Size", "nucleotide.F_ST",
            "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
            "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop",
            "Sister_clade_richness", "heterozygosity", "inbreeding",  "Kipps.BillDepth",
            "DAPC", "STRUCTURE", "BAPS")

df.temp <- df[,traits]
rownames(df.temp) <- df$species
df.temp <- df.temp %>% dplyr::filter(species != "Elaenia_pelzelni" | species != "Leucippus_chlorocercus")

# Data attributes:
model <- lm(av_contig_length ~ loci + total_bp, data = df.temp)
vif(model)
# these are all below 5

# Genetic structure
model <- lm(Dxy ~ Av_groups + SNPs_per_bp + SNPs_per_locus + SNPs + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs, data = df.temp)
vif(model)
# Three are higher than 5: SNPs_per_bp, SNPs_per_locus, SNPs. Remove latter two for next:
model <- lm(Dxy ~ Av_groups + SNPs_per_bp + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs, data = df.temp)
vif(model)
# these are all below 5

# gene flow metrics
model <- lm(Nei.G_ST ~ e_statistic + Dxy, data = df.temp)
vif(model)
# these are all below 5

# Genetic diversity and population size
model <- lm(theta ~ nuc_div + seg_sites + D + heterozygosity + Range.Size, data = df.temp)
vif(model)
# nuc_div and seg_sites are above 5
model <- lm(theta ~ D + heterozygosity + Range.Size, data = df.temp)
vif(model)
# these are all below 5


# Species traits
model <- lm(Mass ~ Beak.Length_Nares + Hand.wing.Index + Trophic.Niche, data = df.temp)
vif(model)
# these are all below 5

# Speciation dynamics
model <- lm(stem_length ~ subtending_branch + Sister_clade_richness, data = df.temp)
vif(model)
# these are all below 5



# all together now (tons above 5. Weird.) Try trimming down to just genetics metrics to get a better set
model <- lm(av_contig_length ~ loci + total_bp + Dxy + Av_groups + SNPs_per_bp + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs + Nei.G_ST + e_statistic + theta + D + heterozygosity + Range.Size + Mass + Beak.Length_Nares + Hand.wing.Index + Trophic.Niche + stem_length + subtending_branch + Sister_clade_richness, data = df.temp)
# remove contig length descriptive stats
model <- lm(Dxy ~ Av_groups + SNPs_per_bp + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs + Nei.G_ST + e_statistic + theta + D + heterozygosity + Range.Size + Mass + Beak.Length_Nares + Hand.wing.Index + Trophic.Niche + stem_length + subtending_branch + Sister_clade_richness, data = df.temp)
# remove traits
model <- lm(Dxy ~ Av_groups + SNPs_per_bp + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs + Nei.G_ST + e_statistic + theta + D + heterozygosity + stem_length + subtending_branch + Sister_clade_richness, data = df.temp)
# remove subtending branch stuff, so just pop gen stats now
model <- lm(Dxy ~ Av_groups + SNPs_per_bp + mtDNA.branch.length + Av_UCE_gene_tree_length + pairwise_diffs + Nei.G_ST + e_statistic + theta + D + heterozygosity, data = df.temp)
# remove pairwise_diffs and SNPs_per_bp
model <- lm(Dxy ~ Av_groups + mtDNA.branch.length + Av_UCE_gene_tree_length + Nei.G_ST + e_statistic + theta + D + heterozygosity, data = df.temp)
model <- lm(mtDNA.branch.length ~ Av_groups + Dxy + Av_UCE_gene_tree_length + Nei.G_ST + e_statistic + theta + D + heterozygosity, data = df.temp)
model <- lm(Nei.G_ST ~ Dxy + mtDNA.branch.length + Av_UCE_gene_tree_length + Av_groups + e_statistic + theta + D + heterozygosity, data = df.temp)
vif(model)

# this last one looks like a good set. Plot
vif_values <- vif(model)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
abline(v = 5, lwd = 3, lty = 2)
vif_values[vif_values > 5]

# correlation matrix for the final set
traits <- c("Dxy", "Av_groups", "mtDNA.branch.length", "Av_UCE_gene_tree_length", "Nei.G_ST", "e_statistic", "theta", "D", "heterozygosity")
df.temp <- df %>% dplyr::filter(species != "Elaenia_pelzelni" & species != "Leucippus_chlorocercus")
df.temp <- df.temp[,traits]
cor(df.temp)

# traits to remove from final analysis
remove <- c("pairwise_diffs", "SNPs_per_bp", "SNPs_per_locus", "SNPs", "nuc_div", "seg_sites")

#------------------------------------------------------------------------
# traits <- colnames(df)[c(13, 16:30, 32, 36:37)] # old version
# traits <- colnames(df)[c(7, 11:25, 29, 31, 38, 44, 46:47, 49:50, 54:61)]
# traits <- c("Av_groups", "BAPS", "DAPC", "STRUCTURE")

traits <- c("Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
            "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
            "D", "theta", "theta.bp", "seg_sites", "pairwise_diffs", "nuc_div",
            "seg_sites_per_bp", "subtending_branch", "stem_length",
            "mtDNA.branch.length", "Beak.Length_Nares", "Hand.wing.Index",
            "Trophic.Niche", "Mass", "Range.Size", "nucleotide.F_ST",
            "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
            "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop",
            "Sister_clade_richness", "heterozygosity", "inbreeding",  "Kipps.BillDepth",
            "DAPC", "STRUCTURE", "BAPS")

# just functional traits
# traits <- colnames(df)[30:47] 

#rescale for effect sizes
df.temp <- df[,traits]
rownames(df.temp) <- df$species
df.temp <- as.matrix(df.temp)
df.temp <- scale(df.temp)
df.rescaled <- as.data.frame(df.temp)
df.rescaled$species <- df$species
df.rescaled$habitat <- df$habitat
# df.temp$habitat <- df$habitat
# df.temp$species <- df$species

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


fname <- paste0("3_results/phylANOVA.results.v6.csv")
write.csv(df.anova, file = fname, row.names = FALSE)

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


fname <- paste0("3_results/phylANOVA.results.v6.Ochthornis.Dendroplex.reassigned.csv")
write.csv(df.anova, file = fname, row.names = FALSE)
df.anova <- read.csv(fname)

#------------------------------------------------------------------------
# run with Ochthornis and Dendroplex removed

#rescale for effect sizes
df.temp <- df[,traits]
rownames(df.temp) <- df$species
df.temp <- as.matrix(df.temp)
df.temp <- scale(df.temp)
df.rescaled <- as.data.frame(df.temp)
df.rescaled$species <- df$species
df.rescaled$habitat <- df$habitat

df.temp <- df.rescaled %>% 
  dplyr::filter(species != "Ochthornis_littoralis" & species != "Dendroplex_kienerii") 

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


fname <- paste0("3_results/phylANOVA.results.v6.Ochthornis.Dendroplex.removed.csv")
write.csv(df.anova, file = fname, row.names = FALSE)
df.anova <- read.csv(fname)



#------------------------------------------------------------------------
# run anova with island split into early and late successional stages
# note that a different function is necessary here: phy_anova_island.islandsplit

#rescale for effect sizes
df.temp <- df[,traits]
rownames(df.temp) <- df$species
df.temp <- as.matrix(df.temp)
df.temp <- scale(df.temp)
df.rescaled <- as.data.frame(df.temp)
df.rescaled$species <- df$species
df.rescaled$habitat <- df$habitat

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
# pgls


df.o <- df.rescaled 

traits <- c("Av_groups", "SNPs", "loci", "SNPs_per_locus", "av_contig_length",
            "total_bp", "SNPs_per_bp", "Av_UCE_gene_tree_length",
            "D", "theta", "theta.bp", "seg_sites", "pairwise_diffs", "nuc_div",
            "seg_sites_per_bp", "subtending_branch", "stem_length",
            "mtDNA.branch.length", "Beak.Length_Nares", "Hand.wing.Index",
            "Trophic.Niche", "Mass", "Range.Size", "nucleotide.F_ST",
            "Nei.G_ST", "Dxy", "e_statistic", "e_slope", "a_statistic",
            "a_slope", "Nei.G_ST.pop", "nucleotide.F_ST.pop", "Dxy.pop",
            "Sister_clade_richness", "heterozygosity", "inbreeding",  "Kipps.BillDepth",
            "DAPC", "STRUCTURE", "BAPS")

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

fname <- paste0("3_results/pgls.results.v4.csv")
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

fname <- paste0("3_results/pgls.results.v4.Ochthornis.Dendroplex.reassigned.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)


# run with Ochthornis and Dendroplex removed


df.o <- df.o %>% 
  dplyr::filter(species != "Ochthornis_littoralis" & species != "Dendroplex_kienerii") 

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

fname <- paste0("3_results/pgls.results.v4.Ochthornis.Dendroplex.removed.csv")
write.csv(df.pgls, file = fname, row.names = FALSE)




#------------------------------------------------------------------------
#plotting

Av_groups <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Av_groups)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nNumber of groups (average)") 

Dxy <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Dxy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nDxy") 

Fst <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Nei.G_ST)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") 

SNPs <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nTotal SNPs") 

SNPs_per_locus <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs_per_locus)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSNPs per locus") 

SNPs_per_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = SNPs_per_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSNPs per bp") 

loci <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = loci)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nNumber of loci") 

av_contig_length <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = av_contig_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage contig length") 

total_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = total_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nTotal bp") 

Av_gene_tree <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Av_UCE_gene_tree_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nAverage UCE gene tree length") 

Taj_D <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = D)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nTajima's D") 

theta <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = theta)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nTheta") 

theta.bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = theta.bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nTheta/bp") 

seg_sites <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = seg_sites)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSegregating sites") 

seg_sites_per_bp <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = seg_sites_per_bp)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSegregating sites per bp") 

pairwise_diffs <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = pairwise_diffs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage pairwise differences") 

nuc_div <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = nuc_div)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\nNucleotide diversity") 

Av_mtdna_tree <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = mtDNA.branch.length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nAverage mtDNA tree length") 


IBD <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = e_slope)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nIBD (slope)") 

heterozygosity <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = heterozygosity)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nheterozygosity") 

inbreeding <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = inbreeding)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\ninbreeding") 


# subtending branch

stem <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = stem_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nStem branch length") 

subtending <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = subtending_branch)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSubtending branch length") 

Sister_clade_richness <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Sister_clade_richness)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSister_clade_richness") 



# traits 

hwi <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand-wing Index") 

Range.Size <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log(Range.Size))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nRange size (log)") 

mass <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log(Mass))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nMass (log)") 

Beak.Length_Nares <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Length_Nares)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak length at nares") 


Trophic.Niche <- df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Trophic.Niche)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nDiet categories") 

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

# for reviewer
df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = STRUCTURE)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSTRUCTURE groups") 
ggsave("4_plots/STRUCTURE.png")

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = DAPC)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nDAPC groups") 
ggsave("4_plots/DAPC.png")

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = BAPS)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBAPS groups") 
ggsave("4_plots/BAPS.png")



# plots for a few morphometric traits, not used in paper

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Kipps.BillDepth)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nKipps.BillDepth") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Kipps.Distance_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nKipps.Distance") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Depth_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak.Depth") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Beak.Width_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nBeak.Width") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Secondary1_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nSecondary1") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length_1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length") 

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length.Mass)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length / Mass") 


df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Wing.Length.Mass.residuals)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nWing.Length / Mass, residuals") 

df %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index.Mass)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand.wing.Index / Mass") 

df %>% dplyr::filter(species != "Leucippus_chlorocercus" & species != "Phaethornis_bourcieri" & species != "Phaethornis_hispidus") %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = Hand.wing.Index.Mass.residuals)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\nHand.wing.Index / Mass, residuals") 



#------------------------------------------------------------------------
# plot as grid

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

ggsave("4_plots/multipage_boxplot.pdf", 
       width = 12, height = 12, units = "in")


gt.sub <- arrangeGrob(subtending, stem, Sister_clade_richness,
                      nrow = 1, ncol = 3)
as_ggplot(gt.sub) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 15,
                  x = c(0.001, 0.33, 0.66), 
                  y = c(1,  1,  1))
ggsave("4_plots/multipage_boxplot.subtending.pdf", 
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
ggsave("4_plots/multipage_boxplot.traits.png", 
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
ggsave("4_plots/multipage_boxplot.descriptive.png", 
       width = 6, height = 6, units = "in")

#------------------------------------------------------------------------
# plot Fst and HWI separately for figures
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
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") +
  annotate("text", x=1.2, y=0.55, label = label.Fst, size = 3.6)

ggsave("4_plots/Fst.pdf", 
       width = 6,
       height = 6,
       units = c("in"))


label.HWI <- paste0("p = ", prettyNum(df.pgls["Hand.wing.Index","p"], digits = 2),
                    ", t = ", prettyNum(df.pgls["Hand.wing.Index","t"], digits = 3), 
                    "\nestimate = ", prettyNum(df.pgls["Hand.wing.Index","estimate"], digits = 3))

df %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log(Hand.wing.Index))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nHand-wing index (log)") +
  annotate("text", x=1, y=4.4, label = label.HWI, size = 3.6)

ggsave("4_plots/HWI.pdf", width = 6,
       height = 6,
       units = c("in"))

#------------------------------------------------------------------------
# plot Fst and HWI separately
# but reverse order of habitats on x axis for presentation slides

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
  geom_jitter(position=position_jitter(height = 0.1, width = 0.35)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nFst") +
  annotate("text", x=3, y=0.55, label = label.Fst, size = 3.6)

  ggsave("4_plots/Fst.for_talk.v1.pdf")
  

label.HWI <- paste0("PGLS\np = ", prettyNum(df.pgls[19,3], digits = 2),
                    ", t = ", prettyNum(df.pgls[19,4], digits = 3), 
                    "\nestimate = -", prettyNum(df.pgls[19,2], digits = 3))

df %>%
  mutate(habitat = fct_relevel(habitat, "upland", "floodplain", "island")) %>%
  ggplot( aes(x = habitat, y = log(Hand.wing.Index))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(height = 0.1, width = 0.35)) +
  scale_y_continuous(limits=c(1.8, 4.6), 
                     breaks = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  ggtitle("\n") +
  xlab("\n") +
  ylab("\n\nHand-wing index (log)") +
  annotate("text", x=3, y=4.5, label = label.HWI, size = 3.6) 

  ggsave("4_plots/HWI.for_talk.v1.pdf")


#------------------------------------------------------------------------
# pgls
# HWI vs genetic traits
# just an idea

df.o <- df

traits <- colnames(df)[c(13, 16:30, 33, 35:36)]

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

save.image()

quit()

