library(dplyr)
library("adegenet")
library("ade4")
library("factoextra")
library("hierfstat")
library(nlme)
library(pegas)
# library(multiNe)
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_all_snps_for_adegenet/")


df.pop.het <- data.frame(species=character(),
                         population=character(),
                         Hobs=double(), 
                         # theta.pop=double(),
                         stringsAsFactors=FALSE) 

#------------------------------------------------------------------------

# Campephilus_melanoleucos 
Campephilus_melanoleucos <- read.structure("Campephilus_melanoleucos_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                           n.ind=11, n.loc=5657, onerowperind=FALSE, 
                                           col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Campephilus_melanoleucos.grp <- find.clusters(Campephilus_melanoleucos, n.pca = 100, 
                            method = "kmeans", choose.n.clust = FALSE, 
                            criterion = "min")
Campephilus_melanoleucos$pop <- Campephilus_melanoleucos.grp$grp

Campephilus_melanoleucos.div <- summary(Campephilus_melanoleucos)
Campephilus_melanoleucos.heteroz <- mean(Campephilus_melanoleucos.div$Hobs)

# within pop theta
# Campephilus_melanoleucos.loci <- as.loci(Campephilus_melanoleucos)
# Campephilus_melanoleucos.theta <- sapply(Campephilus_melanoleucos.loci, function(x) theta.h(x))
# Campephilus_melanoleucos.theta <- Campephilus_melanoleucos.theta[-1]
# Campephilus_melanoleucos.theta.mean <- mean(Campephilus_melanoleucos.theta, na.rm = TRUE)

Campephilus_melanoleucos.df <- as.data.frame(list("Campephilus_melanoleucos", "all", Campephilus_melanoleucos.heteroz)) #, Campephilus_melanoleucos.theta.mean))
colnames(Campephilus_melanoleucos.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Campephilus_melanoleucos.df)

Campephilus_melanoleucos.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Campephilus_melanoleucos.df)

#------------------------------------------------------------------------

# Campephilus_rubricollis
Campephilus_rubricollis <- read.structure("Campephilus_rubricollis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                          n.ind=11, n.loc=6286, onerowperind=FALSE, 
                                          col.lab=1, col.pop=0, col.others=0, row.marknames=1)

Campephilus_rubricollis.grp <- find.clusters(Campephilus_rubricollis, n.pca = 100, 
                            method = "kmeans", choose.n.clust = FALSE, 
                            criterion = "min")
Campephilus_rubricollis$pop <- Campephilus_rubricollis.grp$grp

Campephilus_rubricollis.div <- summary(Campephilus_rubricollis)
Campephilus_rubricollis.heteroz <- mean(Campephilus_rubricollis.div$Hobs)

Campephilus_rubricollis.div.1 <- summary(Campephilus_rubricollis[pop=1])
Campephilus_rubricollis.heteroz.1 <- mean(Campephilus_rubricollis.div.1$Hobs, na.rm = TRUE)

Campephilus_rubricollis.div.2 <- summary(Campephilus_rubricollis[pop=2])
Campephilus_rubricollis.heteroz.2 <- mean(Campephilus_rubricollis.div.2$Hobs, na.rm = TRUE)


# within pop theta
# Campephilus_rubricollis.loci <- as.loci(Campephilus_rubricollis)
# Campephilus_rubricollis.theta <- sapply(Campephilus_rubricollis.loci, function(x) theta.h(x))
# Campephilus_rubricollis.theta.mean <- mean(Campephilus_rubricollis.theta[-1], na.rm = TRUE)
# 
# Campephilus_rubricollis.loci.1 <- as.loci(Campephilus_rubricollis[pop=1])
# Campephilus_rubricollis.theta.1 <- sapply(Campephilus_rubricollis.loci.1, function(x) theta.h(x))
# Campephilus_rubricollis.theta.mean.1 <- mean(Campephilus_rubricollis.theta.1[-1], na.rm = TRUE)
# 
# Campephilus_rubricollis.loci.2 <- as.loci(Campephilus_rubricollis[pop=2])
# Campephilus_rubricollis.theta.2 <- sapply(Campephilus_rubricollis.loci.2, function(x) theta.h(x))
# Campephilus_rubricollis.theta.mean.2 <- mean(Campephilus_rubricollis.theta.2[-1], na.rm = TRUE)

Campephilus_rubricollis.df <- as.data.frame(list("Campephilus_rubricollis", "all", Campephilus_rubricollis.heteroz)) #, Campephilus_rubricollis.theta.mean))
Campephilus_rubricollis.df.1 <- as.data.frame(list("Campephilus_rubricollis", "one", Campephilus_rubricollis.heteroz.1)) #, Campephilus_rubricollis.theta.mean.1))
Campephilus_rubricollis.df.2 <- as.data.frame(list("Campephilus_rubricollis", "two", Campephilus_rubricollis.heteroz.2)) #, Campephilus_rubricollis.theta.mean.2))

colnames(Campephilus_rubricollis.df) <- colnames(Campephilus_rubricollis.df.1) <- colnames(Campephilus_rubricollis.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Campephilus_rubricollis.df, Campephilus_rubricollis.df.1, Campephilus_rubricollis.df.2)



# n <- length(Campephilus_rubricollis.loci)-1
# 
# s <- length(seg.sites(Campephilus_rubricollis$tab))
# theta.s(s, n)
# 
# s <- length(seg.sites(Campephilus_rubricollis[pop=1]$tab))
# theta.s(s, n)
# 
# s <- length(seg.sites(Campephilus_rubricollis[pop=2]$tab))
# theta.s(s, n)


#------------------------------------------------------------------------

# Cantorchilus_leucotis
Cantorchilus_leucotis <- read.structure("Cantorchilus_leucotis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                        n.ind=10, n.loc=9635, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Cantorchilus_leucotis.grp <- find.clusters(Cantorchilus_leucotis, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Cantorchilus_leucotis$pop <- Cantorchilus_leucotis.grp$grp

Cantorchilus_leucotis.div <- summary(Cantorchilus_leucotis)
Cantorchilus_leucotis.heteroz <- mean(Cantorchilus_leucotis.div$Hobs)


# within pop theta
# Cantorchilus_leucotis.loci <- as.loci(Cantorchilus_leucotis)
# Cantorchilus_leucotis.theta <- sapply(Cantorchilus_leucotis.loci, function(x) theta.h(x))
# Cantorchilus_leucotis.theta <- Cantorchilus_leucotis.theta[-1]
# Cantorchilus_leucotis.theta.mean <- mean(Cantorchilus_leucotis.theta, na.rm = TRUE)

Cantorchilus_leucotis.df <- as.data.frame(list("Cantorchilus_leucotis", "all", Cantorchilus_leucotis.heteroz)) #, Cantorchilus_leucotis.theta.mean))

colnames(Cantorchilus_leucotis.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Cantorchilus_leucotis.df)

Cantorchilus_leucotis.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Cantorchilus_leucotis.df)

#------------------------------------------------------------------------

# Celeus_flavus
Celeus_flavus <- read.structure("Celeus_flavus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                n.ind=9, n.loc=3507, onerowperind=FALSE, 
                                col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Celeus_flavus.grp <- find.clusters(Celeus_flavus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Celeus_flavus$pop <- Celeus_flavus.grp$grp

Celeus_flavus.div <- summary(Celeus_flavus)
Celeus_flavus.heteroz <- mean(Celeus_flavus.div$Hobs)


# within pop theta
Celeus_flavus.loci <- as.loci(Celeus_flavus)
Celeus_flavus.theta <- sapply(Celeus_flavus.loci, function(x) theta.h(x))
Celeus_flavus.theta <- Celeus_flavus.theta[-1]
Celeus_flavus.theta.mean <- mean(Celeus_flavus.theta, na.rm = TRUE)

Celeus_flavus.df <- as.data.frame(list("Celeus_flavus", "all", Celeus_flavus.heteroz)) #, Celeus_flavus.theta.mean))
colnames(Celeus_flavus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Celeus_flavus.df)
Celeus_flavus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Celeus_flavus.df)

#------------------------------------------------------------------------
# Celeus_grammicus
Celeus_grammicus <- read.structure("Celeus_grammicus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                   n.ind=11, n.loc=6690, onerowperind=FALSE, 
                                   col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Celeus_grammicus.grp <- find.clusters(Celeus_grammicus, n.pca = 100, 
                                   method = "kmeans", choose.n.clust = FALSE, 
                                   criterion = "min")
Celeus_grammicus$pop <- Celeus_grammicus.grp$grp

Celeus_grammicus.div <- summary(Celeus_grammicus)
Celeus_grammicus.heteroz <- mean(Celeus_grammicus.div$Hobs)

# within pop theta
# Celeus_grammicus.loci <- as.loci(Celeus_grammicus)
# Celeus_grammicus.theta <- sapply(Celeus_grammicus.loci, function(x) theta.h(x))
# Celeus_grammicus.theta <- Celeus_grammicus.theta[-1]
# Celeus_grammicus.theta.mean <- mean(Celeus_grammicus.theta, na.rm = TRUE)

Celeus_grammicus.df <- as.data.frame(list("Celeus_grammicus", "all", Celeus_grammicus.heteroz)) #, Celeus_grammicus.theta.mean))
colnames(Celeus_grammicus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Celeus_grammicus.df)

Celeus_grammicus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Celeus_grammicus.df)

#------------------------------------------------------------------------
# Conirostrum_bicolor
Conirostrum_bicolor <- read.structure("conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=7, n.loc=7991, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Conirostrum_bicolor.grp <- find.clusters(Conirostrum_bicolor, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Conirostrum_bicolor$pop <- Conirostrum_bicolor.grp$grp

Conirostrum_bicolor.div <- summary(Conirostrum_bicolor)
Conirostrum_bicolor.heteroz <- mean(Conirostrum_bicolor.div$Hobs)

Conirostrum_bicolor.div.1 <- summary(Conirostrum_bicolor[pop=1])
Conirostrum_bicolor.heteroz.1 <- mean(Conirostrum_bicolor.div.1$Hobs, na.rm = TRUE)

Conirostrum_bicolor.div.2 <- summary(Conirostrum_bicolor[pop=2])
Conirostrum_bicolor.heteroz.2 <- mean(Conirostrum_bicolor.div.2$Hobs, na.rm = TRUE)

# within pop theta
# Conirostrum_bicolor.loci <- as.loci(Conirostrum_bicolor)
# Conirostrum_bicolor.theta <- sapply(Conirostrum_bicolor.loci, function(x) theta.h(x))
# Conirostrum_bicolor.theta.mean <- mean(Conirostrum_bicolor.theta[-1], na.rm = TRUE)
# 
# Conirostrum_bicolor.loci.1 <- as.loci(Conirostrum_bicolor[pop=1])
# Conirostrum_bicolor.theta.1 <- sapply(Conirostrum_bicolor.loci.1, function(x) theta.h(x))
# Conirostrum_bicolor.theta.mean.1 <- mean(Conirostrum_bicolor.theta.1[-1], na.rm = TRUE)
# 
# Conirostrum_bicolor.loci.2 <- as.loci(Conirostrum_bicolor[pop=2])
# Conirostrum_bicolor.theta.2 <- sapply(Conirostrum_bicolor.loci.2, function(x) theta.h(x))
# Conirostrum_bicolor.theta.mean.2 <- mean(Conirostrum_bicolor.theta.2[-1], na.rm = TRUE)

# n <- length(Conirostrum_bicolor.loci)-1
# n <- df["Conirostrum_bicolor", "total_bp"]
# s <- length(seg.sites(Conirostrum_bicolor$tab))
# Conirostrum_bicolor.theta.mean <- theta.s(s, n) 
# s <- length(seg.sites(Conirostrum_bicolor[pop=1]$tab))
# Conirostrum_bicolor.theta.mean.1 <- theta.s(s, n) 
# s <- length(seg.sites(Conirostrum_bicolor[pop=2]$tab))
# Conirostrum_bicolor.theta.mean.2 <- theta.s(s, n) 

Conirostrum_bicolor.df <- as.data.frame(list("Conirostrum_bicolor", "all", Conirostrum_bicolor.heteroz)) #, Conirostrum_bicolor.theta.mean))
Conirostrum_bicolor.df.1 <- as.data.frame(list("Conirostrum_bicolor", "one", Conirostrum_bicolor.heteroz.1)) #, Conirostrum_bicolor.theta.mean.1))
Conirostrum_bicolor.df.2 <- as.data.frame(list("Conirostrum_bicolor", "two", Conirostrum_bicolor.heteroz.2)) #, Conirostrum_bicolor.theta.mean.2))

colnames(Conirostrum_bicolor.df) <- colnames(Conirostrum_bicolor.df.1) <- colnames(Conirostrum_bicolor.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Conirostrum_bicolor.df, Conirostrum_bicolor.df.1, Conirostrum_bicolor.df.2)




#------------------------------------------------------------------------
# Conirostrum_margaritae
Conirostrum_margaritae <- read.structure("conirostrum_margaritae_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=4, n.loc=1430, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Conirostrum_margaritae.grp <- find.clusters(Conirostrum_margaritae, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Conirostrum_margaritae$pop <- Conirostrum_margaritae.grp$grp

Conirostrum_margaritae.div <- summary(Conirostrum_margaritae)
Conirostrum_margaritae.heteroz <- mean(Conirostrum_margaritae.div$Hobs)

# Conirostrum_margaritae.loci <- as.loci(Conirostrum_margaritae)
# Conirostrum_margaritae.theta <- sapply(Conirostrum_margaritae.loci, function(x) theta.h(x))
# Conirostrum_margaritae.theta <- Conirostrum_margaritae.theta[-1]
# Conirostrum_margaritae.theta.mean <- mean(Conirostrum_margaritae.theta, na.rm = TRUE)

Conirostrum_margaritae.df <- as.data.frame(list("Conirostrum_margaritae", "all", Conirostrum_margaritae.heteroz)) #, Conirostrum_margaritae.theta.mean))
colnames(Conirostrum_margaritae.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Conirostrum_margaritae.df)

Conirostrum_margaritae.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Conirostrum_margaritae.df)


#------------------------------------------------------------------------
# Cranioleuca_vulpecula
Cranioleuca_vulpecula <- read.structure("cranioleuca_vulpecula_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=7, n.loc=1636, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Cranioleuca_vulpecula.grp <- find.clusters(Cranioleuca_vulpecula, n.pca = 100, 
                            method = "kmeans", choose.n.clust = FALSE, 
                            criterion = "min")
Cranioleuca_vulpecula$pop <- Cranioleuca_vulpecula.grp$grp

Cranioleuca_vulpecula.div <- summary(Cranioleuca_vulpecula)
Cranioleuca_vulpecula.heteroz <- mean(Cranioleuca_vulpecula.div$Hobs)

# within pop theta
# Cranioleuca_vulpecula.loci <- as.loci(Cranioleuca_vulpecula)
# Cranioleuca_vulpecula.theta <- sapply(Cranioleuca_vulpecula.loci, function(x) theta.h(x))
# Cranioleuca_vulpecula.theta <- Cranioleuca_vulpecula.theta[-1]
# Cranioleuca_vulpecula.theta.mean <- mean(Cranioleuca_vulpecula.theta, na.rm = TRUE)

Cranioleuca_vulpecula.df <- as.data.frame(list("Cranioleuca_vulpecula", "all", Cranioleuca_vulpecula.heteroz)) #, Cranioleuca_vulpecula.theta.mean))
colnames(Cranioleuca_vulpecula.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Cranioleuca_vulpecula.df)

Cranioleuca_vulpecula.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Cranioleuca_vulpecula.df)

#------------------------------------------------------------------------
# Crypturellus_undulatus
Crypturellus_undulatus <- read.structure("Crypturellus_undulatus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                         n.ind=9, n.loc=6513, onerowperind=FALSE, 
                                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Crypturellus_undulatus.grp <- find.clusters(Crypturellus_undulatus, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Crypturellus_undulatus$pop <- Crypturellus_undulatus.grp$grp

Crypturellus_undulatus.div <- summary(Crypturellus_undulatus)
Crypturellus_undulatus.heteroz <- mean(Crypturellus_undulatus.div$Hobs)

Crypturellus_undulatus.div.1 <- summary(Crypturellus_undulatus[pop=1])
Crypturellus_undulatus.heteroz.1 <- mean(Crypturellus_undulatus.div.1$Hobs, na.rm = TRUE)

Crypturellus_undulatus.div.2 <- summary(Crypturellus_undulatus[pop=2])
Crypturellus_undulatus.heteroz.2 <- mean(Crypturellus_undulatus.div.2$Hobs, na.rm = TRUE)

# within pop theta
# Crypturellus_undulatus.loci <- as.loci(Crypturellus_undulatus)
# Crypturellus_undulatus.theta <- sapply(Crypturellus_undulatus.loci, function(x) theta.h(x))
# Crypturellus_undulatus.theta.mean <- mean(Crypturellus_undulatus.theta[-1], na.rm = TRUE)
# 
# Crypturellus_undulatus.loci.1 <- as.loci(Crypturellus_undulatus[pop=1])
# Crypturellus_undulatus.theta.1 <- sapply(Crypturellus_undulatus.loci.1, function(x) theta.h(x))
# Crypturellus_undulatus.theta.mean.1 <- mean(Crypturellus_undulatus.theta.1[-1], na.rm = TRUE)
# 
# Crypturellus_undulatus.loci.2 <- as.loci(Crypturellus_undulatus[pop=2])
# Crypturellus_undulatus.theta.2 <- sapply(Crypturellus_undulatus.loci.2, function(x) theta.h(x))
# Crypturellus_undulatus.theta.mean.2 <- mean(Crypturellus_undulatus.theta.2[-1], na.rm = TRUE)

Crypturellus_undulatus.df <- as.data.frame(list("Crypturellus_undulatus", "all", Crypturellus_undulatus.heteroz)) #, Crypturellus_undulatus.theta.mean))
Crypturellus_undulatus.df.1 <- as.data.frame(list("Crypturellus_undulatus", "one", Crypturellus_undulatus.heteroz.1)) #, Crypturellus_undulatus.theta.mean.1))
Crypturellus_undulatus.df.2 <- as.data.frame(list("Crypturellus_undulatus", "two", Crypturellus_undulatus.heteroz.2)) #, NA))

colnames(Crypturellus_undulatus.df) <- colnames(Crypturellus_undulatus.df.1) <- colnames(Crypturellus_undulatus.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Crypturellus_undulatus.df, Crypturellus_undulatus.df.1, Crypturellus_undulatus.df.2)

#------------------------------------------------------------------------

# Crypturellus_variegatus
Crypturellus_variegatus <- read.structure("Crypturellus_variegatus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                          n.ind=10, n.loc=11149, onerowperind=FALSE, 
                                          col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Crypturellus_variegatus.grp <- find.clusters(Crypturellus_variegatus, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Crypturellus_variegatus$pop <- Crypturellus_variegatus.grp$grp

Crypturellus_variegatus.div <- summary(Crypturellus_variegatus)
Crypturellus_variegatus.heteroz <- mean(Crypturellus_variegatus.div$Hobs)

Crypturellus_variegatus.div.1 <- summary(Crypturellus_variegatus[pop=1])
Crypturellus_variegatus.heteroz.1 <- mean(Crypturellus_variegatus.div.1$Hobs, na.rm = TRUE)

Crypturellus_variegatus.div.2 <- summary(Crypturellus_variegatus[pop=2])
Crypturellus_variegatus.heteroz.2 <- mean(Crypturellus_variegatus.div.2$Hobs, na.rm = TRUE)

# within pop theta
# Crypturellus_variegatus.loci <- as.loci(Crypturellus_variegatus)
# Crypturellus_variegatus.theta <- sapply(Crypturellus_variegatus.loci, function(x) theta.h(x))
# Crypturellus_variegatus.theta.mean <- mean(Crypturellus_variegatus.theta[-1], na.rm = TRUE)
# 
# Crypturellus_variegatus.loci.1 <- as.loci(Crypturellus_variegatus[pop=1])
# Crypturellus_variegatus.theta.1 <- sapply(Crypturellus_variegatus.loci.1, function(x) theta.h(x))
# Crypturellus_variegatus.theta.mean.1 <- mean(Crypturellus_variegatus.theta.1[-1], na.rm = TRUE)
# 
# Crypturellus_variegatus.loci.2 <- as.loci(Crypturellus_variegatus[pop=2])
# Crypturellus_variegatus.theta.2 <- sapply(Crypturellus_variegatus.loci.2, function(x) theta.h(x))
# Crypturellus_variegatus.theta.mean.2 <- mean(Crypturellus_variegatus.theta.2[-1], na.rm = TRUE)

Crypturellus_variegatus.df <- as.data.frame(list("Crypturellus_variegatus", "all", Crypturellus_variegatus.heteroz)) #, Crypturellus_variegatus.theta.mean))
Crypturellus_variegatus.df.1 <- as.data.frame(list("Crypturellus_variegatus", "one", Crypturellus_variegatus.heteroz.1)) #, Crypturellus_variegatus.theta.mean.1))
Crypturellus_variegatus.df.2 <- as.data.frame(list("Crypturellus_variegatus", "two", Crypturellus_variegatus.heteroz.2)) #, Crypturellus_variegatus.theta.mean.2))

colnames(Crypturellus_variegatus.df) <- colnames(Crypturellus_variegatus.df.1) <- colnames(Crypturellus_variegatus.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Crypturellus_variegatus.df, Crypturellus_variegatus.df.1, Crypturellus_variegatus.df.2)

#------------------------------------------------------------------------

# Dendroplex_kienerii
Dendroplex_kienerii <- read.structure("dendroplex_kienerii_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=7, n.loc=2905, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Dendroplex_kienerii.grp <- find.clusters(Dendroplex_kienerii, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Dendroplex_kienerii$pop <- Dendroplex_kienerii.grp$grp

Dendroplex_kienerii.div <- summary(Dendroplex_kienerii)
Dendroplex_kienerii.heteroz <- mean(Dendroplex_kienerii.div$Hobs)

# within pop theta
# Dendroplex_kienerii.loci <- as.loci(Dendroplex_kienerii)
# Dendroplex_kienerii.theta <- sapply(Dendroplex_kienerii.loci, function(x) theta.h(x))
# Dendroplex_kienerii.theta <- Dendroplex_kienerii.theta[-1]
# Dendroplex_kienerii.theta.mean <- mean(Dendroplex_kienerii.theta, na.rm = TRUE)

Dendroplex_kienerii.df <- as.data.frame(list("Dendroplex_kienerii", "all", Dendroplex_kienerii.heteroz)) #, Dendroplex_kienerii.theta.mean))
colnames(Dendroplex_kienerii.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Dendroplex_kienerii.df)

Dendroplex_kienerii.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Dendroplex_kienerii.df)


#------------------------------------------------------------------------

# Elaenia_pelzelni
# can't find clusters
Elaenia_pelzelni <- read.structure("elaenia_pelzelni_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                   n.ind=2, n.loc=1914, onerowperind=FALSE, 
                                   col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Elaenia_pelzelni.grp <- find.clusters(Elaenia_pelzelni, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Elaenia_pelzelni$pop <- Elaenia_pelzelni.grp$grp

Elaenia_pelzelni.div <- summary(Elaenia_pelzelni)
Elaenia_pelzelni.heteroz <- mean(Elaenia_pelzelni.div$Hobs)
Elaenia_pelzelni.df <- as.data.frame(list("Elaenia_pelzelni", "all", Elaenia_pelzelni.heteroz))
colnames(Elaenia_pelzelni.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Elaenia_pelzelni.df)

Elaenia_pelzelni.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Elaenia_pelzelni.df)

#------------------------------------------------------------------------
# Formicarius_analis
Formicarius_analis <- read.structure("Formicarius_analis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                     n.ind=11, n.loc=11601, onerowperind=FALSE, 
                                     col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Formicarius_analis.grp <- find.clusters(Formicarius_analis, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Formicarius_analis$pop <- Formicarius_analis.grp$grp

Formicarius_analis.div <- summary(Formicarius_analis)
Formicarius_analis.heteroz <- mean(Formicarius_analis.div$Hobs)

Formicarius_analis.div.1 <- summary(Formicarius_analis[pop=1])
Formicarius_analis.heteroz.1 <- mean(Formicarius_analis.div.1$Hobs, na.rm = TRUE)

Formicarius_analis.div.2 <- summary(Formicarius_analis[pop=2])
Formicarius_analis.heteroz.2 <- mean(Formicarius_analis.div.2$Hobs, na.rm = TRUE)

Formicarius_analis.df <- as.data.frame(list("Formicarius_analis", "all", Formicarius_analis.heteroz))
Formicarius_analis.df.1 <- as.data.frame(list("Formicarius_analis", "one", Formicarius_analis.heteroz.1))
Formicarius_analis.df.2 <- as.data.frame(list("Formicarius_analis", "two", Formicarius_analis.heteroz.2))

colnames(Formicarius_analis.df) <- colnames(Formicarius_analis.df.1) <- colnames(Formicarius_analis.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Formicarius_analis.df, Formicarius_analis.df.1, Formicarius_analis.df.2)



#------------------------------------------------------------------------
# Formicarius_colma
Formicarius_colma <- read.structure("Formicarius_colma_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                    n.ind=11, n.loc=7051, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Formicarius_colma.grp <- find.clusters(Formicarius_colma, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Formicarius_colma$pop <- Formicarius_colma.grp$grp

Formicarius_colma.div <- summary(Formicarius_colma)
Formicarius_colma.heteroz <- mean(Formicarius_colma.div$Hobs)

Formicarius_colma.div.1 <- summary(Formicarius_colma[pop=1])
Formicarius_colma.heteroz.1 <- mean(Formicarius_colma.div.1$Hobs, na.rm = TRUE)

Formicarius_colma.div.2 <- summary(Formicarius_colma[pop=2])
Formicarius_colma.heteroz.2 <- mean(Formicarius_colma.div.2$Hobs, na.rm = TRUE)

Formicarius_colma.df <- as.data.frame(list("Formicarius_colma", "all", Formicarius_colma.heteroz))
Formicarius_colma.df.1 <- as.data.frame(list("Formicarius_colma", "one", Formicarius_colma.heteroz.1))
Formicarius_colma.df.2 <- as.data.frame(list("Formicarius_colma", "two", Formicarius_colma.heteroz.2))

colnames(Formicarius_colma.df) <- colnames(Formicarius_colma.df.1) <- colnames(Formicarius_colma.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Formicarius_colma.df, Formicarius_colma.df.1, Formicarius_colma.df.2)


#------------------------------------------------------------------------
# Furnarius_minor
Furnarius_minor <- read.structure("furnarius_minor_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                  n.ind=5, n.loc=1792, onerowperind=FALSE, 
                                  col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Furnarius_minor.grp <- find.clusters(Furnarius_minor, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Furnarius_minor$pop <- Furnarius_minor.grp$grp

Furnarius_minor.div <- summary(Furnarius_minor)
Furnarius_minor.heteroz <- mean(Furnarius_minor.div$Hobs)
Furnarius_minor.df <- as.data.frame(list("Furnarius_minor", "all", Furnarius_minor.heteroz))
colnames(Furnarius_minor.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Furnarius_minor.df)

Furnarius_minor.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Furnarius_minor.df)



#------------------------------------------------------------------------
# Glaucidium_brasilianum
Glaucidium_brasilianum <- read.structure("Glaucidium_brasilianum_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                         n.ind=8, n.loc=2167, onerowperind=FALSE, 
                                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Glaucidium_brasilianum.grp <- find.clusters(Glaucidium_brasilianum, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Glaucidium_brasilianum$pop <- Glaucidium_brasilianum.grp$grp

Glaucidium_brasilianum.div <- summary(Glaucidium_brasilianum)
Glaucidium_brasilianum.heteroz <- mean(Glaucidium_brasilianum.div$Hobs)

Glaucidium_brasilianum.div.1 <- summary(Glaucidium_brasilianum[pop=1])
Glaucidium_brasilianum.heteroz.1 <- mean(Glaucidium_brasilianum.div.1$Hobs, na.rm = TRUE)

Glaucidium_brasilianum.df <- as.data.frame(list("Glaucidium_brasilianum", "all", Glaucidium_brasilianum.heteroz))
Glaucidium_brasilianum.df.1 <- as.data.frame(list("Glaucidium_brasilianum", "one", Glaucidium_brasilianum.heteroz.1))

colnames(Glaucidium_brasilianum.df) <- colnames(Glaucidium_brasilianum.df.1) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Glaucidium_brasilianum.df, Glaucidium_brasilianum.df.1)


#------------------------------------------------------------------------

# Glaucidium_hardyi
Glaucidium_hardyi <- read.structure("Glaucidium_hardyi_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                    n.ind=10, n.loc=3711, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Glaucidium_hardyi.grp <- find.clusters(Glaucidium_hardyi, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Glaucidium_hardyi$pop <- Glaucidium_hardyi.grp$grp

Glaucidium_hardyi.div <- summary(Glaucidium_hardyi)
Glaucidium_hardyi.heteroz <- mean(Glaucidium_hardyi.div$Hobs)
Glaucidium_hardyi.df <- as.data.frame(list("Glaucidium_hardyi", "all", Glaucidium_hardyi.heteroz))
colnames(Glaucidium_hardyi.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Glaucidium_hardyi.df)

Glaucidium_hardyi.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Glaucidium_hardyi.df)


#------------------------------------------------------------------------

# Hylophylax_naevia
Hylophylax_naevia <- read.structure("Hylophylax_naevia_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                    n.ind=9, n.loc=11906, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Hylophylax_naevia.grp <- find.clusters(Hylophylax_naevia, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Hylophylax_naevia$pop <- Hylophylax_naevia.grp$grp

Hylophylax_naevia.div <- summary(Hylophylax_naevia)
Hylophylax_naevia.heteroz <- mean(Hylophylax_naevia.div$Hobs)

Hylophylax_naevia.div.1 <- summary(Hylophylax_naevia[pop=1])
Hylophylax_naevia.heteroz.1 <- mean(Hylophylax_naevia.div.1$Hobs, na.rm = TRUE)

Hylophylax_naevia.div.2 <- summary(Hylophylax_naevia[pop=2])
Hylophylax_naevia.heteroz.2 <- mean(Hylophylax_naevia.div.2$Hobs, na.rm = TRUE)

Hylophylax_naevia.df <- as.data.frame(list("Hylophylax_naevia", "all", Hylophylax_naevia.heteroz))
Hylophylax_naevia.df.1 <- as.data.frame(list("Hylophylax_naevia", "one", Hylophylax_naevia.heteroz.1))
Hylophylax_naevia.df.2 <- as.data.frame(list("Hylophylax_naevia", "two", Hylophylax_naevia.heteroz.2))

colnames(Hylophylax_naevia.df) <- colnames(Hylophylax_naevia.df.1) <- colnames(Hylophylax_naevia.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Hylophylax_naevia.df, Hylophylax_naevia.df.1, Hylophylax_naevia.df.2)

#------------------------------------------------------------------------

# Hylophylax_punctulata
Hylophylax_punctulata <- read.structure("Hylophylax_punctulata_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=9, n.loc=11994, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Hylophylax_punctulata.grp <- find.clusters(Hylophylax_punctulata, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Hylophylax_punctulata$pop <- Hylophylax_punctulata.grp$grp

Hylophylax_punctulata.div <- summary(Hylophylax_punctulata)
Hylophylax_punctulata.heteroz <- mean(Hylophylax_punctulata.div$Hobs)

Hylophylax_punctulata.div.1 <- summary(Hylophylax_punctulata[pop=1])
Hylophylax_punctulata.heteroz.1 <- mean(Hylophylax_punctulata.div.1$Hobs, na.rm = TRUE)

Hylophylax_punctulata.div.2 <- summary(Hylophylax_punctulata[pop=2])
Hylophylax_punctulata.heteroz.2 <- mean(Hylophylax_punctulata.div.2$Hobs, na.rm = TRUE)

Hylophylax_punctulata.df <- as.data.frame(list("Hylophylax_punctulata", "all", Hylophylax_punctulata.heteroz))
Hylophylax_punctulata.df.1 <- as.data.frame(list("Hylophylax_punctulata", "one", Hylophylax_punctulata.heteroz.1))
Hylophylax_punctulata.df.2 <- as.data.frame(list("Hylophylax_punctulata", "two", Hylophylax_punctulata.heteroz.2))

colnames(Hylophylax_punctulata.df) <- colnames(Hylophylax_punctulata.df.1) <- colnames(Hylophylax_punctulata.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Hylophylax_punctulata.df, Hylophylax_punctulata.df.1, Hylophylax_punctulata.df.2)

#------------------------------------------------------------------------
# Knipolegus_orenocensis
Knipolegus_orenocensis <- read.structure("knipolegus_orenocensis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=12, n.loc=5874, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Knipolegus_orenocensis.grp <- find.clusters(Knipolegus_orenocensis, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Knipolegus_orenocensis$pop <- Knipolegus_orenocensis.grp$grp

Knipolegus_orenocensis.div <- summary(Knipolegus_orenocensis)
Knipolegus_orenocensis.heteroz <- mean(Knipolegus_orenocensis.div$Hobs)

Knipolegus_orenocensis.div.1 <- summary(Knipolegus_orenocensis[pop=1])
Knipolegus_orenocensis.heteroz.1 <- mean(Knipolegus_orenocensis.div.1$Hobs, na.rm = TRUE)

Knipolegus_orenocensis.div.2 <- summary(Knipolegus_orenocensis[pop=2])
Knipolegus_orenocensis.heteroz.2 <- mean(Knipolegus_orenocensis.div.2$Hobs, na.rm = TRUE)

Knipolegus_orenocensis.df <- as.data.frame(list("Knipolegus_orenocensis", "all", Knipolegus_orenocensis.heteroz))
Knipolegus_orenocensis.df.1 <- as.data.frame(list("Knipolegus_orenocensis", "one", Knipolegus_orenocensis.heteroz.1))
Knipolegus_orenocensis.df.2 <- as.data.frame(list("Knipolegus_orenocensis", "two", Knipolegus_orenocensis.heteroz.2))

colnames(Knipolegus_orenocensis.df) <- colnames(Knipolegus_orenocensis.df.1) <- colnames(Knipolegus_orenocensis.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Knipolegus_orenocensis.df, Knipolegus_orenocensis.df.1, Knipolegus_orenocensis.df.2)


#------------------------------------------------------------------------

# Leucippus_chlorocercus
Leucippus_chlorocercus <- read.structure("leucippus_chlorocercus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                         n.ind=4, n.loc=2326, onerowperind=FALSE, 
                                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Leucippus_chlorocercus.grp <- find.clusters(Leucippus_chlorocercus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Leucippus_chlorocercus$pop <- Leucippus_chlorocercus.grp$grp

Leucippus_chlorocercus.div <- summary(Leucippus_chlorocercus)
Leucippus_chlorocercus.heteroz <- mean(Leucippus_chlorocercus.div$Hobs)
Leucippus_chlorocercus.df <- as.data.frame(list("Leucippus_chlorocercus", "all", Leucippus_chlorocercus.heteroz))
colnames(Leucippus_chlorocercus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Leucippus_chlorocercus.df)

Leucippus_chlorocercus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Leucippus_chlorocercus.df)


#------------------------------------------------------------------------

# Mazaria_propinqua
Mazaria_propinqua <- read.structure("mazaria_propinqua_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                    n.ind=7, n.loc=3094, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Mazaria_propinqua.grp <- find.clusters(Mazaria_propinqua, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Mazaria_propinqua$pop <- Mazaria_propinqua.grp$grp

Mazaria_propinqua.div <- summary(Mazaria_propinqua)
Mazaria_propinqua.heteroz <- mean(Mazaria_propinqua.div$Hobs)
Mazaria_propinqua.df <- as.data.frame(list("Mazaria_propinqua", "all", Mazaria_propinqua.heteroz))
colnames(Mazaria_propinqua.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Mazaria_propinqua.df)

Mazaria_propinqua.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Mazaria_propinqua.df)


#------------------------------------------------------------------------

# Megascops_choliba
Megascops_choliba <- read.structure("Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                    n.ind=11, n.loc=3370, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Megascops_choliba.grp <- find.clusters(Megascops_choliba, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Megascops_choliba$pop <- Megascops_choliba.grp$grp

Megascops_choliba.div <- summary(Megascops_choliba)
Megascops_choliba.heteroz <- mean(Megascops_choliba.div$Hobs)
Megascops_choliba.df <- as.data.frame(list("Megascops_choliba", "all", Megascops_choliba.heteroz))
colnames(Megascops_choliba.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Megascops_choliba.df)

Megascops_choliba.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Megascops_choliba.df)


#------------------------------------------------------------------------
# Megascops_watsonii
Megascops_watsonii <- read.structure("Megascops_watsonii_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                     n.ind=11, n.loc=6389, onerowperind=FALSE, 
                                     col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Megascops_watsonii.grp <- find.clusters(Megascops_watsonii, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Megascops_watsonii$pop <- Megascops_watsonii.grp$grp

Megascops_watsonii.div <- summary(Megascops_watsonii)
Megascops_watsonii.heteroz <- mean(Megascops_watsonii.div$Hobs)
Megascops_watsonii.df <- as.data.frame(list("Megascops_watsonii", "all", Megascops_watsonii.heteroz))
colnames(Megascops_watsonii.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Megascops_watsonii.df)

Megascops_watsonii.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Megascops_watsonii.df)



#------------------------------------------------------------------------
# Monasa_morphoeus
Monasa_morphoeus <- read.structure("Monasa_morphoeus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                   n.ind=10, n.loc=6331, onerowperind=FALSE, 
                                   col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Monasa_morphoeus.grp <- find.clusters(Monasa_morphoeus, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Monasa_morphoeus$pop <- Monasa_morphoeus.grp$grp

Monasa_morphoeus.div <- summary(Monasa_morphoeus)
Monasa_morphoeus.heteroz <- mean(Monasa_morphoeus.div$Hobs)

Monasa_morphoeus.div.1 <- summary(Monasa_morphoeus[pop=1])
Monasa_morphoeus.heteroz.1 <- mean(Monasa_morphoeus.div.1$Hobs, na.rm = TRUE)

Monasa_morphoeus.div.2 <- summary(Monasa_morphoeus[pop=2])
Monasa_morphoeus.heteroz.2 <- mean(Monasa_morphoeus.div.2$Hobs, na.rm = TRUE)

Monasa_morphoeus.df <- as.data.frame(list("Monasa_morphoeus", "all", Monasa_morphoeus.heteroz))
Monasa_morphoeus.df.1 <- as.data.frame(list("Monasa_morphoeus", "one", Monasa_morphoeus.heteroz.1))
Monasa_morphoeus.df.2 <- as.data.frame(list("Monasa_morphoeus", "two", Monasa_morphoeus.heteroz.2))

colnames(Monasa_morphoeus.df) <- colnames(Monasa_morphoeus.df.1) <- colnames(Monasa_morphoeus.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Monasa_morphoeus.df, Monasa_morphoeus.df.1, Monasa_morphoeus.df.2)


#------------------------------------------------------------------------

# Monasa_nigrifrons
Monasa_nigrifrons <- read.structure("Monasa_nigrifrons_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                    n.ind=11, n.loc=6873, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Monasa_nigrifrons.grp <- find.clusters(Monasa_nigrifrons, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Monasa_nigrifrons$pop <- Monasa_nigrifrons.grp$grp

Monasa_nigrifrons.div <- summary(Monasa_nigrifrons)
Monasa_nigrifrons.heteroz <- mean(Monasa_nigrifrons.div$Hobs)

Monasa_nigrifrons.div.1 <- summary(Monasa_nigrifrons[pop=1])
Monasa_nigrifrons.heteroz.1 <- mean(Monasa_nigrifrons.div.1$Hobs, na.rm = TRUE)

Monasa_nigrifrons.div.2 <- summary(Monasa_nigrifrons[pop=2])
Monasa_nigrifrons.heteroz.2 <- mean(Monasa_nigrifrons.div.2$Hobs, na.rm = TRUE)

Monasa_nigrifrons.df <- as.data.frame(list("Monasa_nigrifrons", "all", Monasa_nigrifrons.heteroz))
Monasa_nigrifrons.df.1 <- as.data.frame(list("Monasa_nigrifrons", "one", Monasa_nigrifrons.heteroz.1))
Monasa_nigrifrons.df.2 <- as.data.frame(list("Monasa_nigrifrons", "two", Monasa_nigrifrons.heteroz.2))

colnames(Monasa_nigrifrons.df) <- colnames(Monasa_nigrifrons.df.1) <- colnames(Monasa_nigrifrons.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Monasa_nigrifrons.df, Monasa_nigrifrons.df.1, Monasa_nigrifrons.df.2)

#------------------------------------------------------------------------

# Myrmeciza_fortis
Myrmeciza_fortis <- read.structure("Myrmeciza_fortis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                   n.ind=7, n.loc=6339, onerowperind=FALSE, 
                                   col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Myrmeciza_fortis.grp <- find.clusters(Myrmeciza_fortis, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Myrmeciza_fortis$pop <- Myrmeciza_fortis.grp$grp

Myrmeciza_fortis.div <- summary(Myrmeciza_fortis)
Myrmeciza_fortis.heteroz <- mean(Myrmeciza_fortis.div$Hobs)

Myrmeciza_fortis.div.1 <- summary(Myrmeciza_fortis[pop=1])
Myrmeciza_fortis.heteroz.1 <- mean(Myrmeciza_fortis.div.1$Hobs, na.rm = TRUE)

Myrmeciza_fortis.div.2 <- summary(Myrmeciza_fortis[pop=2])
Myrmeciza_fortis.heteroz.2 <- mean(Myrmeciza_fortis.div.2$Hobs, na.rm = TRUE)

Myrmeciza_fortis.df <- as.data.frame(list("Myrmeciza_fortis", "all", Myrmeciza_fortis.heteroz))
Myrmeciza_fortis.df.1 <- as.data.frame(list("Myrmeciza_fortis", "one", Myrmeciza_fortis.heteroz.1))
Myrmeciza_fortis.df.2 <- as.data.frame(list("Myrmeciza_fortis", "two", Myrmeciza_fortis.heteroz.2))

colnames(Myrmeciza_fortis.df) <- colnames(Myrmeciza_fortis.df.1) <- colnames(Myrmeciza_fortis.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmeciza_fortis.df, Myrmeciza_fortis.df.1, Myrmeciza_fortis.df.2)


#------------------------------------------------------------------------

# Myrmeciza_hyperythra
Myrmeciza_hyperythra <- read.structure("Myrmeciza_hyperythra_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                       n.ind=7, n.loc=4381, onerowperind=FALSE, 
                                       col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Myrmeciza_hyperythra.grp <- find.clusters(Myrmeciza_hyperythra, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Myrmeciza_hyperythra$pop <- Myrmeciza_hyperythra.grp$grp

Myrmeciza_hyperythra.div <- summary(Myrmeciza_hyperythra)
Myrmeciza_hyperythra.heteroz <- mean(Myrmeciza_hyperythra.div$Hobs)
Myrmeciza_hyperythra.df <- as.data.frame(list("Myrmeciza_hyperythra", "all", Myrmeciza_hyperythra.heteroz))
colnames(Myrmeciza_hyperythra.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmeciza_hyperythra.df)

Myrmeciza_hyperythra.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Myrmeciza_hyperythra.df)


#------------------------------------------------------------------------
# Myrmoborus_leucophrys
Myrmoborus_leucophrys <- read.structure("Myrmoborus_leucophrys_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=7, n.loc=9508, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Myrmoborus_leucophrys.grp <- find.clusters(Myrmoborus_leucophrys, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Myrmoborus_leucophrys$pop <- Myrmoborus_leucophrys.grp$grp

Myrmoborus_leucophrys.div <- summary(Myrmoborus_leucophrys)
Myrmoborus_leucophrys.heteroz <- mean(Myrmoborus_leucophrys.div$Hobs)

Myrmoborus_leucophrys.div.1 <- summary(Myrmoborus_leucophrys[pop=1])
Myrmoborus_leucophrys.heteroz.1 <- mean(Myrmoborus_leucophrys.div.1$Hobs, na.rm = TRUE)

Myrmoborus_leucophrys.div.2 <- summary(Myrmoborus_leucophrys[pop=2])
Myrmoborus_leucophrys.heteroz.2 <- mean(Myrmoborus_leucophrys.div.2$Hobs, na.rm = TRUE)

Myrmoborus_leucophrys.df <- as.data.frame(list("Myrmoborus_leucophrys", "all", Myrmoborus_leucophrys.heteroz))
Myrmoborus_leucophrys.df.1 <- as.data.frame(list("Myrmoborus_leucophrys", "one", Myrmoborus_leucophrys.heteroz.1))
Myrmoborus_leucophrys.df.2 <- as.data.frame(list("Myrmoborus_leucophrys", "two", Myrmoborus_leucophrys.heteroz.2))

colnames(Myrmoborus_leucophrys.df) <- colnames(Myrmoborus_leucophrys.df.1) <- colnames(Myrmoborus_leucophrys.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmoborus_leucophrys.df, Myrmoborus_leucophrys.df.1, Myrmoborus_leucophrys.df.2)


#------------------------------------------------------------------------
# Myrmoborus_lugubris
Myrmoborus_lugubris <- read.structure("myrmoborus_lugubris_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                      n.ind=8, n.loc=3280, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Myrmoborus_lugubris.grp <- find.clusters(Myrmoborus_lugubris, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Myrmoborus_lugubris$pop <- Myrmoborus_lugubris.grp$grp

Myrmoborus_lugubris.div <- summary(Myrmoborus_lugubris)
Myrmoborus_lugubris.heteroz <- mean(Myrmoborus_lugubris.div$Hobs)

Myrmoborus_lugubris.div.1 <- summary(Myrmoborus_lugubris[pop=1])
Myrmoborus_lugubris.heteroz.1 <- mean(Myrmoborus_lugubris.div.1$Hobs, na.rm = TRUE)

Myrmoborus_lugubris.div.2 <- summary(Myrmoborus_lugubris[pop=2])
Myrmoborus_lugubris.heteroz.2 <- mean(Myrmoborus_lugubris.div.2$Hobs, na.rm = TRUE)

Myrmoborus_lugubris.df <- as.data.frame(list("Myrmoborus_lugubris", "all", Myrmoborus_lugubris.heteroz))
Myrmoborus_lugubris.df.1 <- as.data.frame(list("Myrmoborus_lugubris", "one", Myrmoborus_lugubris.heteroz.1))
Myrmoborus_lugubris.df.2 <- as.data.frame(list("Myrmoborus_lugubris", "two", Myrmoborus_lugubris.heteroz.2))

colnames(Myrmoborus_lugubris.df) <- colnames(Myrmoborus_lugubris.df.1) <- colnames(Myrmoborus_lugubris.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmoborus_lugubris.df, Myrmoborus_lugubris.df.1, Myrmoborus_lugubris.df.2)

#------------------------------------------------------------------------
# Myrmoborus_myotherinus
Myrmoborus_myotherinus <- read.structure("Myrmoborus_myotherinus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                         n.ind=7, n.loc=10230, onerowperind=FALSE, 
                                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Myrmoborus_myotherinus.grp <- find.clusters(Myrmoborus_myotherinus, n.pca = 100, 
                                      method = "kmeans", choose.n.clust = FALSE, 
                                      criterion = "min")
Myrmoborus_myotherinus$pop <- Myrmoborus_myotherinus.grp$grp

Myrmoborus_myotherinus.div <- summary(Myrmoborus_myotherinus)
Myrmoborus_myotherinus.heteroz <- mean(Myrmoborus_myotherinus.div$Hobs)

Myrmoborus_myotherinus.div.1 <- summary(Myrmoborus_myotherinus[pop=1])
Myrmoborus_myotherinus.heteroz.1 <- mean(Myrmoborus_myotherinus.div.1$Hobs, na.rm = TRUE)

Myrmoborus_myotherinus.div.2 <- summary(Myrmoborus_myotherinus[pop=2])
Myrmoborus_myotherinus.heteroz.2 <- mean(Myrmoborus_myotherinus.div.2$Hobs, na.rm = TRUE)

Myrmoborus_myotherinus.df <- as.data.frame(list("Myrmoborus_myotherinus", "all", Myrmoborus_myotherinus.heteroz))
Myrmoborus_myotherinus.df.1 <- as.data.frame(list("Myrmoborus_myotherinus", "one", Myrmoborus_myotherinus.heteroz.1))
Myrmoborus_myotherinus.df.2 <- as.data.frame(list("Myrmoborus_myotherinus", "two", Myrmoborus_myotherinus.heteroz.2))

colnames(Myrmoborus_myotherinus.df) <- colnames(Myrmoborus_myotherinus.df.1) <- colnames(Myrmoborus_myotherinus.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmoborus_myotherinus.df, Myrmoborus_myotherinus.df.1, Myrmoborus_myotherinus.df.2)

#------------------------------------------------------------------------
# Myrmochanes_hemileucus
Myrmochanes_hemileucus <- read.structure("myrmochanes_hemileucus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                         n.ind=6, n.loc=1463, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Myrmochanes_hemileucus.grp <- find.clusters(Myrmochanes_hemileucus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Myrmochanes_hemileucus$pop <- Myrmochanes_hemileucus.grp$grp

Myrmochanes_hemileucus.div <- summary(Myrmochanes_hemileucus)
Myrmochanes_hemileucus.heteroz <- mean(Myrmochanes_hemileucus.div$Hobs)
Myrmochanes_hemileucus.df <- as.data.frame(list("Myrmochanes_hemileucus", "all", Myrmochanes_hemileucus.heteroz))
colnames(Myrmochanes_hemileucus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmochanes_hemileucus.df)

Myrmochanes_hemileucus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Myrmochanes_hemileucus.df)

#------------------------------------------------------------------------
# Myrmotherula_assimilis
Myrmotherula_assimilis <- read.structure("myrmotherula_assimilis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                         n.ind=7, n.loc=3433, onerowperind=FALSE, 
                                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Myrmotherula_assimilis.grp <- find.clusters(Myrmotherula_assimilis, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Myrmotherula_assimilis$pop <- Myrmotherula_assimilis.grp$grp

Myrmotherula_assimilis.div <- summary(Myrmotherula_assimilis)
Myrmotherula_assimilis.heteroz <- mean(Myrmotherula_assimilis.div$Hobs)
Myrmotherula_assimilis.df <- as.data.frame(list("Myrmotherula_assimilis", "all", Myrmotherula_assimilis.heteroz))
colnames(Myrmotherula_assimilis.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmotherula_assimilis.df)

Myrmotherula_assimilis.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Myrmotherula_assimilis.df)

#------------------------------------------------------------------------
# Myrmotherula_klagesi
Myrmotherula_klagesi <- read.structure("myrmotherula_klagesi_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                       n.ind=7, n.loc=3703, onerowperind=FALSE, 
                                       col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Myrmotherula_klagesi.grp <- find.clusters(Myrmotherula_klagesi, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Myrmotherula_klagesi$pop <- Myrmotherula_klagesi.grp$grp

Myrmotherula_klagesi.div <- summary(Myrmotherula_klagesi)
Myrmotherula_klagesi.heteroz <- mean(Myrmotherula_klagesi.div$Hobs)
Myrmotherula_klagesi.df <- as.data.frame(list("Myrmotherula_klagesi", "all", Myrmotherula_klagesi.heteroz))
colnames(Myrmotherula_klagesi.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Myrmotherula_klagesi.df)

Myrmotherula_klagesi.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Myrmotherula_klagesi.df)


#------------------------------------------------------------------------
# Ochthornis_littoralis
Ochthornis_littoralis <- read.structure("ochthornis_littoralis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=11, n.loc=3335, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Ochthornis_littoralis.grp <- find.clusters(Ochthornis_littoralis, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Ochthornis_littoralis$pop <- Ochthornis_littoralis.grp$grp

Ochthornis_littoralis.div <- summary(Ochthornis_littoralis)
Ochthornis_littoralis.heteroz <- mean(Ochthornis_littoralis.div$Hobs)
Ochthornis_littoralis.df <- as.data.frame(list("Ochthornis_littoralis", "all", Ochthornis_littoralis.heteroz))
colnames(Ochthornis_littoralis.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Ochthornis_littoralis.df)

Ochthornis_littoralis.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Ochthornis_littoralis.df)


#------------------------------------------------------------------------

# Phaethornis_bourcieri
Phaethornis_bourcieri <- read.structure("Phaethornis_bourcieri_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=7, n.loc=9184, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Phaethornis_bourcieri.grp <- find.clusters(Phaethornis_bourcieri, n.pca = 100, 
                                           method = "kmeans", choose.n.clust = FALSE, 
                                           criterion = "min")
Phaethornis_bourcieri$pop <- Phaethornis_bourcieri.grp$grp

Phaethornis_bourcieri.div <- summary(Phaethornis_bourcieri)
Phaethornis_bourcieri.heteroz <- mean(Phaethornis_bourcieri.div$Hobs, na.rm = TRUE)

Phaethornis_bourcieri.div.1 <- summary(Phaethornis_bourcieri[pop=1])
Phaethornis_bourcieri.heteroz.1 <- mean(Phaethornis_bourcieri.div.1$Hobs, na.rm = TRUE)

Phaethornis_bourcieri.div.2 <- summary(Phaethornis_bourcieri[pop=2])
Phaethornis_bourcieri.heteroz.2 <- mean(Phaethornis_bourcieri.div.2$Hobs, na.rm = TRUE)

#Hs(Phaethornis_bourcieri)

Phaethornis_bourcieri.df <- as.data.frame(list("Phaethornis_bourcieri", "all", Phaethornis_bourcieri.heteroz))
Phaethornis_bourcieri.df.1 <- as.data.frame(list("Phaethornis_bourcieri", "one", Phaethornis_bourcieri.heteroz.1))
Phaethornis_bourcieri.df.2 <- as.data.frame(list("Phaethornis_bourcieri", "two", Phaethornis_bourcieri.heteroz.2))

colnames(Phaethornis_bourcieri.df) <- colnames(Phaethornis_bourcieri.df.1) <- colnames(Phaethornis_bourcieri.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Phaethornis_bourcieri.df, Phaethornis_bourcieri.df.1, Phaethornis_bourcieri.df.2)

#------------------------------------------------------------------------

# Phaethornis_hispidus
Phaethornis_hispidus <- read.structure("Phaethornis_hispidus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                       n.ind=7, n.loc=7099, onerowperind=FALSE, 
                                       col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Phaethornis_hispidus.grp <- find.clusters(Phaethornis_hispidus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Phaethornis_hispidus$pop <- Phaethornis_hispidus.grp$grp

Phaethornis_hispidus.div <- summary(Phaethornis_hispidus)
Phaethornis_hispidus.heteroz <- mean(Phaethornis_hispidus.div$Hobs)
Phaethornis_hispidus.df <- as.data.frame(list("Phaethornis_hispidus", "all", Phaethornis_hispidus.heteroz))
colnames(Phaethornis_hispidus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Phaethornis_hispidus.df)

Phaethornis_hispidus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Phaethornis_hispidus.df)



#------------------------------------------------------------------------
# Pheugopedius_coraya
Pheugopedius_coraya <- read.structure("Pheugopedius_coraya_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                      n.ind=11, n.loc=14067, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Pheugopedius_coraya.grp <- find.clusters(Pheugopedius_coraya, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Pheugopedius_coraya$pop <- Pheugopedius_coraya.grp$grp

Pheugopedius_coraya.div <- summary(Pheugopedius_coraya)
Pheugopedius_coraya.heteroz <- mean(Pheugopedius_coraya.div$Hobs)

Pheugopedius_coraya.div.1 <- summary(Pheugopedius_coraya[pop=1])
Pheugopedius_coraya.heteroz.1 <- mean(Pheugopedius_coraya.div.1$Hobs, na.rm = TRUE)

Pheugopedius_coraya.div.2 <- summary(Pheugopedius_coraya[pop=2])
Pheugopedius_coraya.heteroz.2 <- mean(Pheugopedius_coraya.div.2$Hobs, na.rm = TRUE)

Pheugopedius_coraya.df <- as.data.frame(list("Pheugopedius_coraya", "all", Pheugopedius_coraya.heteroz))
Pheugopedius_coraya.df.1 <- as.data.frame(list("Pheugopedius_coraya", "one", Pheugopedius_coraya.heteroz.1))
Pheugopedius_coraya.df.2 <- as.data.frame(list("Pheugopedius_coraya", "two", Pheugopedius_coraya.heteroz.2))

colnames(Pheugopedius_coraya.df) <- colnames(Pheugopedius_coraya.df.1) <- colnames(Pheugopedius_coraya.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Pheugopedius_coraya.df, Pheugopedius_coraya.df.1, Pheugopedius_coraya.df.2)


#------------------------------------------------------------------------
# Piaya_cayana
Piaya_cayana <- read.structure("Piaya_cayana_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                               n.ind=10, n.loc=11862, onerowperind=FALSE, 
                               col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Piaya_cayana.grp <- find.clusters(Piaya_cayana, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Piaya_cayana$pop <- Piaya_cayana.grp$grp

Piaya_cayana.div <- summary(Piaya_cayana)
Piaya_cayana.heteroz <- mean(Piaya_cayana.div$Hobs)
Piaya_cayana.df <- as.data.frame(list("Piaya_cayana", "all", Piaya_cayana.heteroz))
colnames(Piaya_cayana.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Piaya_cayana.df)

Piaya_cayana.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Piaya_cayana.df)


#------------------------------------------------------------------------
# Piaya_melanogaster
Piaya_melanogaster <- read.structure("Piaya_melanogaster_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                     n.ind=7, n.loc=6550, onerowperind=FALSE, 
                                     col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Piaya_melanogaster.grp <- find.clusters(Piaya_melanogaster, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Piaya_melanogaster$pop <- Piaya_melanogaster.grp$grp

Piaya_melanogaster.div <- summary(Piaya_melanogaster)
Piaya_melanogaster.heteroz <- mean(Piaya_melanogaster.div$Hobs)

Piaya_melanogaster.div.1 <- summary(Piaya_melanogaster[pop=1])
Piaya_melanogaster.heteroz.1 <- mean(Piaya_melanogaster.div.1$Hobs, na.rm = TRUE)

Piaya_melanogaster.div.2 <- summary(Piaya_melanogaster[pop=2])
Piaya_melanogaster.heteroz.2 <- mean(Piaya_melanogaster.div.2$Hobs, na.rm = TRUE)

Piaya_melanogaster.df <- as.data.frame(list("Piaya_melanogaster", "all", Piaya_melanogaster.heteroz))
Piaya_melanogaster.df.1 <- as.data.frame(list("Piaya_melanogaster", "one", Piaya_melanogaster.heteroz.1))
Piaya_melanogaster.df.2 <- as.data.frame(list("Piaya_melanogaster", "two", Piaya_melanogaster.heteroz.2))

colnames(Piaya_melanogaster.df) <- colnames(Piaya_melanogaster.df.1) <- colnames(Piaya_melanogaster.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Piaya_melanogaster.df, Piaya_melanogaster.df.1, Piaya_melanogaster.df.2)


#------------------------------------------------------------------------
# Pipra_erythrocephala
Pipra_erythrocephala <- read.structure("Pipra_erythrocephala_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                       n.ind=11, n.loc=10659, onerowperind=FALSE, 
                                       col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Pipra_erythrocephala.grp <- find.clusters(Pipra_erythrocephala, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Pipra_erythrocephala$pop <- Pipra_erythrocephala.grp$grp

Pipra_erythrocephala.div <- summary(Pipra_erythrocephala)
Pipra_erythrocephala.heteroz <- mean(Pipra_erythrocephala.div$Hobs)

Pipra_erythrocephala.div.1 <- summary(Pipra_erythrocephala[pop=1])
Pipra_erythrocephala.heteroz.1 <- mean(Pipra_erythrocephala.div.1$Hobs, na.rm = TRUE)

Pipra_erythrocephala.div.2 <- summary(Pipra_erythrocephala[pop=2])
Pipra_erythrocephala.heteroz.2 <- mean(Pipra_erythrocephala.div.2$Hobs, na.rm = TRUE)

Pipra_erythrocephala.df <- as.data.frame(list("Pipra_erythrocephala", "all", Pipra_erythrocephala.heteroz))
Pipra_erythrocephala.df.1 <- as.data.frame(list("Pipra_erythrocephala", "one", Pipra_erythrocephala.heteroz.1))
Pipra_erythrocephala.df.2 <- as.data.frame(list("Pipra_erythrocephala", "two", Pipra_erythrocephala.heteroz.2))

colnames(Pipra_erythrocephala.df) <- colnames(Pipra_erythrocephala.df.1) <- colnames(Pipra_erythrocephala.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Pipra_erythrocephala.df, Pipra_erythrocephala.df.1, Pipra_erythrocephala.df.2)



#------------------------------------------------------------------------
# Pipra_filicauda
Pipra_filicauda <- read.structure("Pipra_filicauda_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                  n.ind=7, n.loc=6416, onerowperind=FALSE, 
                                  col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Pipra_filicauda.grp <- find.clusters(Pipra_filicauda, n.pca = 100, 
                                      method = "kmeans", choose.n.clust = FALSE, 
                                      criterion = "min")
Pipra_filicauda$pop <- Pipra_filicauda.grp$grp

Pipra_filicauda.div <- summary(Pipra_filicauda)
Pipra_filicauda.heteroz <- mean(Pipra_filicauda.div$Hobs)

Pipra_filicauda.div.1 <- summary(Pipra_filicauda[pop=1])
Pipra_filicauda.heteroz.1 <- mean(Pipra_filicauda.div.1$Hobs, na.rm = TRUE)

Pipra_filicauda.div.2 <- summary(Pipra_filicauda[pop=2])
Pipra_filicauda.heteroz.2 <- mean(Pipra_filicauda.div.2$Hobs, na.rm = TRUE)

Pipra_filicauda.df <- as.data.frame(list("Pipra_filicauda", "all", Pipra_filicauda.heteroz))
Pipra_filicauda.df.1 <- as.data.frame(list("Pipra_filicauda", "one", Pipra_filicauda.heteroz.1))
Pipra_filicauda.df.2 <- as.data.frame(list("Pipra_filicauda", "two", Pipra_filicauda.heteroz.2))

colnames(Pipra_filicauda.df) <- colnames(Pipra_filicauda.df.1) <- colnames(Pipra_filicauda.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Pipra_filicauda.df, Pipra_filicauda.df.1, Pipra_filicauda.df.2)



#------------------------------------------------------------------------
# Saltator_coerulescens
Saltator_coerulescens <- read.structure("Saltator_coerulescens_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                        n.ind=11, n.loc=13197, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Saltator_coerulescens.grp <- find.clusters(Saltator_coerulescens, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Saltator_coerulescens$pop <- Saltator_coerulescens.grp$grp

Saltator_coerulescens.div <- summary(Saltator_coerulescens)
Saltator_coerulescens.heteroz <- mean(Saltator_coerulescens.div$Hobs)

Saltator_coerulescens.div.1 <- summary(Saltator_coerulescens[pop=1])
Saltator_coerulescens.heteroz.1 <- mean(Saltator_coerulescens.div.1$Hobs, na.rm = TRUE)

Saltator_coerulescens.div.2 <- summary(Saltator_coerulescens[pop=2])
Saltator_coerulescens.heteroz.2 <- mean(Saltator_coerulescens.div.2$Hobs, na.rm = TRUE)

Saltator_coerulescens.df <- as.data.frame(list("Saltator_coerulescens", "all", Saltator_coerulescens.heteroz))
Saltator_coerulescens.df.1 <- as.data.frame(list("Saltator_coerulescens", "one", Saltator_coerulescens.heteroz.1))
Saltator_coerulescens.df.2 <- as.data.frame(list("Saltator_coerulescens", "two", Saltator_coerulescens.heteroz.2))

colnames(Saltator_coerulescens.df) <- colnames(Saltator_coerulescens.df.1) <- colnames(Saltator_coerulescens.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Saltator_coerulescens.df, Saltator_coerulescens.df.1, Saltator_coerulescens.df.2)



#------------------------------------------------------------------------
# Saltator_grossus
Saltator_grossus <- read.structure("Saltator_grossus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                   n.ind=11, n.loc=10453, onerowperind=FALSE, 
                                   col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Saltator_grossus.grp <- find.clusters(Saltator_grossus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Saltator_grossus$pop <- Saltator_grossus.grp$grp

Saltator_grossus.div <- summary(Saltator_grossus)
Saltator_grossus.heteroz <- mean(Saltator_grossus.div$Hobs)
Saltator_grossus.df <- as.data.frame(list("Saltator_grossus", "all", Saltator_grossus.heteroz))
colnames(Saltator_grossus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Saltator_grossus.df)

Saltator_grossus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Saltator_grossus.df)



#------------------------------------------------------------------------
# Schiffornis_major
Schiffornis_major <- read.structure("Schiffornis_major_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                    n.ind=7, n.loc=9188, onerowperind=FALSE, 
                                    col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Schiffornis_major.grp <- find.clusters(Schiffornis_major, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Schiffornis_major$pop <- Schiffornis_major.grp$grp

Schiffornis_major.div <- summary(Schiffornis_major)
Schiffornis_major.heteroz <- mean(Schiffornis_major.div$Hobs)
Schiffornis_major.df <- as.data.frame(list("Schiffornis_major", "all", Schiffornis_major.heteroz))
colnames(Schiffornis_major.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Schiffornis_major.df)

Schiffornis_major.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Schiffornis_major.df)



#------------------------------------------------------------------------

# Schiffornis_turdina
Schiffornis_turdina <- read.structure("Schiffornis_turdina_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                      n.ind=11, n.loc=11697, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Schiffornis_turdina.grp <- find.clusters(Schiffornis_turdina, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Schiffornis_turdina$pop <- Schiffornis_turdina.grp$grp

Schiffornis_turdina.div <- summary(Schiffornis_turdina)
Schiffornis_turdina.heteroz <- mean(Schiffornis_turdina.div$Hobs)

Schiffornis_turdina.div.1 <- summary(Schiffornis_turdina[pop=1])
Schiffornis_turdina.heteroz.1 <- mean(Schiffornis_turdina.div.1$Hobs, na.rm = TRUE)

Schiffornis_turdina.div.2 <- summary(Schiffornis_turdina[pop=2])
Schiffornis_turdina.heteroz.2 <- mean(Schiffornis_turdina.div.2$Hobs, na.rm = TRUE)

Schiffornis_turdina.df <- as.data.frame(list("Schiffornis_turdina", "all", Schiffornis_turdina.heteroz))
Schiffornis_turdina.df.1 <- as.data.frame(list("Schiffornis_turdina", "one", Schiffornis_turdina.heteroz.1))
Schiffornis_turdina.df.2 <- as.data.frame(list("Schiffornis_turdina", "two", Schiffornis_turdina.heteroz.2))

colnames(Schiffornis_turdina.df) <- colnames(Schiffornis_turdina.df.1) <- colnames(Schiffornis_turdina.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Schiffornis_turdina.df, Schiffornis_turdina.df.1, Schiffornis_turdina.df.2)


#------------------------------------------------------------------------

# Serpophaga_hypoleuca
Serpophaga_hypoleuca <- read.structure("serpophaga_hypoleuca_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                       n.ind=9, n.loc=2259, onerowperind=FALSE, 
                                       col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Serpophaga_hypoleuca.grp <- find.clusters(Serpophaga_hypoleuca, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Serpophaga_hypoleuca$pop <- Serpophaga_hypoleuca.grp$grp

Serpophaga_hypoleuca.div <- summary(Serpophaga_hypoleuca)
Serpophaga_hypoleuca.heteroz <- mean(Serpophaga_hypoleuca.div$Hobs)
Serpophaga_hypoleuca.df <- as.data.frame(list("Serpophaga_hypoleuca", "all", Serpophaga_hypoleuca.heteroz))
colnames(Serpophaga_hypoleuca.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Serpophaga_hypoleuca.df)

Serpophaga_hypoleuca.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Serpophaga_hypoleuca.df)


#------------------------------------------------------------------------
# Stigmatura_napensis
Stigmatura_napensis <- read.structure("stigmatura_napensis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                      n.ind=6, n.loc=2198, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Stigmatura_napensis.grp <- find.clusters(Stigmatura_napensis, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Stigmatura_napensis$pop <- Stigmatura_napensis.grp$grp

Stigmatura_napensis.div <- summary(Stigmatura_napensis)
Stigmatura_napensis.heteroz <- mean(Stigmatura_napensis.div$Hobs)
Stigmatura_napensis.df <- as.data.frame(list("Stigmatura_napensis", "all", Stigmatura_napensis.heteroz))
colnames(Stigmatura_napensis.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Stigmatura_napensis.df)

Stigmatura_napensis.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Stigmatura_napensis.df)



#------------------------------------------------------------------------
# Synallaxis_gujanensis
Synallaxis_gujanensis <- read.structure("Synallaxis_gujanensis_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=7, n.loc=8192, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Synallaxis_gujanensis.grp <- find.clusters(Synallaxis_gujanensis, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Synallaxis_gujanensis$pop <- Synallaxis_gujanensis.grp$grp

Synallaxis_gujanensis.div <- summary(Synallaxis_gujanensis)
Synallaxis_gujanensis.heteroz <- mean(Synallaxis_gujanensis.div$Hobs)

Synallaxis_gujanensis.div.1 <- summary(Synallaxis_gujanensis[pop=1])
Synallaxis_gujanensis.heteroz.1 <- mean(Synallaxis_gujanensis.div.1$Hobs, na.rm = TRUE)

Synallaxis_gujanensis.div.2 <- summary(Synallaxis_gujanensis[pop=2])
Synallaxis_gujanensis.heteroz.2 <- mean(Synallaxis_gujanensis.div.2$Hobs, na.rm = TRUE)

Synallaxis_gujanensis.df <- as.data.frame(list("Synallaxis_gujanensis", "all", Synallaxis_gujanensis.heteroz))
Synallaxis_gujanensis.df.1 <- as.data.frame(list("Synallaxis_gujanensis", "one", Synallaxis_gujanensis.heteroz.1))
Synallaxis_gujanensis.df.2 <- as.data.frame(list("Synallaxis_gujanensis", "two", Synallaxis_gujanensis.heteroz.2))

colnames(Synallaxis_gujanensis.df) <- colnames(Synallaxis_gujanensis.df.1) <- colnames(Synallaxis_gujanensis.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Synallaxis_gujanensis.df, Synallaxis_gujanensis.df.1, Synallaxis_gujanensis.df.2)



#------------------------------------------------------------------------
# Synallaxis_rutilans
Synallaxis_rutilans <- read.structure("Synallaxis_rutilans_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                      n.ind=7, n.loc=6304, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Synallaxis_rutilans.grp <- find.clusters(Synallaxis_rutilans, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Synallaxis_rutilans$pop <- Synallaxis_rutilans.grp$grp

Synallaxis_rutilans.div <- summary(Synallaxis_rutilans)
Synallaxis_rutilans.heteroz <- mean(Synallaxis_rutilans.div$Hobs)

Synallaxis_rutilans.div.1 <- summary(Synallaxis_rutilans[pop=1])
Synallaxis_rutilans.heteroz.1 <- mean(Synallaxis_rutilans.div.1$Hobs, na.rm = TRUE)

Synallaxis_rutilans.div.2 <- summary(Synallaxis_rutilans[pop=2])
Synallaxis_rutilans.heteroz.2 <- mean(Synallaxis_rutilans.div.2$Hobs, na.rm = TRUE)

Synallaxis_rutilans.df <- as.data.frame(list("Synallaxis_rutilans", "all", Synallaxis_rutilans.heteroz))
Synallaxis_rutilans.df.1 <- as.data.frame(list("Synallaxis_rutilans", "one", Synallaxis_rutilans.heteroz.1))
Synallaxis_rutilans.df.2 <- as.data.frame(list("Synallaxis_rutilans", "two", Synallaxis_rutilans.heteroz.2))

colnames(Synallaxis_rutilans.df) <- colnames(Synallaxis_rutilans.df.1) <- colnames(Synallaxis_rutilans.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Synallaxis_rutilans.df, Synallaxis_rutilans.df.1, Synallaxis_rutilans.df.2)



#------------------------------------------------------------------------
# Tachyphonus_cristatus
Tachyphonus_cristatus <- read.structure("Tachyphonus_cristatus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                        n.ind=8, n.loc=9958, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Tachyphonus_cristatus.grp <- find.clusters(Tachyphonus_cristatus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Tachyphonus_cristatus$pop <- Tachyphonus_cristatus.grp$grp

Tachyphonus_cristatus.div <- summary(Tachyphonus_cristatus)
Tachyphonus_cristatus.heteroz <- mean(Tachyphonus_cristatus.div$Hobs)
Tachyphonus_cristatus.df <- as.data.frame(list("Tachyphonus_cristatus", "all", Tachyphonus_cristatus.heteroz))
colnames(Tachyphonus_cristatus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Tachyphonus_cristatus.df)

Tachyphonus_cristatus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Tachyphonus_cristatus.df)



#------------------------------------------------------------------------

# Tachyphonus_luctuosus
Tachyphonus_luctuosus <- read.structure("Tachyphonus_luctuosus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                        n.ind=11, n.loc=10680, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)


Tachyphonus_luctuosus.grp <- find.clusters(Tachyphonus_luctuosus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Tachyphonus_luctuosus$pop <- Tachyphonus_luctuosus.grp$grp

Tachyphonus_luctuosus.div <- summary(Tachyphonus_luctuosus)
Tachyphonus_luctuosus.heteroz <- mean(Tachyphonus_luctuosus.div$Hobs)
Tachyphonus_luctuosus.df <- as.data.frame(list("Tachyphonus_luctuosus", "all", Tachyphonus_luctuosus.heteroz))
colnames(Tachyphonus_luctuosus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Tachyphonus_luctuosus.df)

Tachyphonus_luctuosus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Tachyphonus_luctuosus.df)

#------------------------------------------------------------------------
# Thamnophilus_cryptoleucus
Thamnophilus_cryptoleucus <- read.structure("thamnophilus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                               n.ind=18, n.loc=7015, onerowperind=FALSE, 
                               col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Thamnophilus_cryptoleucus.grp <- find.clusters(Thamnophilus_cryptoleucus, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Thamnophilus_cryptoleucus$pop <- Thamnophilus_cryptoleucus.grp$grp

Thamnophilus_cryptoleucus.div <- summary(Thamnophilus_cryptoleucus)
Thamnophilus_cryptoleucus.heteroz <- mean(Thamnophilus_cryptoleucus.div$Hobs)

Thamnophilus_cryptoleucus.div.1 <- summary(Thamnophilus_cryptoleucus[pop=1])
Thamnophilus_cryptoleucus.heteroz.1 <- mean(Thamnophilus_cryptoleucus.div.1$Hobs, na.rm = TRUE)

Thamnophilus_cryptoleucus.div.2 <- summary(Thamnophilus_cryptoleucus[pop=2])
Thamnophilus_cryptoleucus.heteroz.2 <- mean(Thamnophilus_cryptoleucus.div.2$Hobs, na.rm = TRUE)

Thamnophilus_cryptoleucus.df <- as.data.frame(list("Thamnophilus_cryptoleucus", "all", Thamnophilus_cryptoleucus.heteroz))
Thamnophilus_cryptoleucus.df.1 <- as.data.frame(list("Thamnophilus_cryptoleucus", "one", Thamnophilus_cryptoleucus.heteroz.1))
Thamnophilus_cryptoleucus.df.2 <- as.data.frame(list("Thamnophilus_cryptoleucus", "two", Thamnophilus_cryptoleucus.heteroz.2))

colnames(Thamnophilus_cryptoleucus.df) <- colnames(Thamnophilus_cryptoleucus.df.1) <- colnames(Thamnophilus_cryptoleucus.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Thamnophilus_cryptoleucus.df, Thamnophilus_cryptoleucus.df.1, Thamnophilus_cryptoleucus.df.2)


#------------------------------------------------------------------------
# Trogon_collaris
Trogon_collaris <- read.structure("Trogon_collaris_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                                  n.ind=10, n.loc=10788, onerowperind=FALSE, 
                                  col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Trogon_collaris.grp <- find.clusters(Trogon_collaris, n.pca = 100, 
                                             method = "kmeans", choose.n.clust = FALSE, 
                                             criterion = "min")
Trogon_collaris$pop <- Trogon_collaris.grp$grp

Trogon_collaris.div <- summary(Trogon_collaris)
Trogon_collaris.heteroz <- mean(Trogon_collaris.div$Hobs)

Trogon_collaris.div.1 <- summary(Trogon_collaris[pop=1])
Trogon_collaris.heteroz.1 <- mean(Trogon_collaris.div.1$Hobs, na.rm = TRUE)

Trogon_collaris.div.2 <- summary(Trogon_collaris[pop=2])
Trogon_collaris.heteroz.2 <- mean(Trogon_collaris.div.2$Hobs, na.rm = TRUE)

Trogon_collaris.df <- as.data.frame(list("Trogon_collaris", "all", Trogon_collaris.heteroz))
Trogon_collaris.df.1 <- as.data.frame(list("Trogon_collaris", "one", Trogon_collaris.heteroz.1))
Trogon_collaris.df.2 <- as.data.frame(list("Trogon_collaris", "two", Trogon_collaris.heteroz.2))

colnames(Trogon_collaris.df) <- colnames(Trogon_collaris.df.1) <- colnames(Trogon_collaris.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Trogon_collaris.df, Trogon_collaris.df.1, Trogon_collaris.df.2)

#------------------------------------------------------------------------
# Trogon_rufus
Trogon_rufus <- read.structure("Trogon_rufus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.AmazonOnly.edited.str", 
                               n.ind=11, n.loc=13083, onerowperind=FALSE, 
                               col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Trogon_rufus.grp <- find.clusters(Trogon_rufus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Trogon_rufus$pop <- Trogon_rufus.grp$grp

Trogon_rufus.div <- summary(Trogon_rufus)
Trogon_rufus.heteroz <- mean(Trogon_rufus.div$Hobs)
Trogon_rufus.df <- as.data.frame(list("Trogon_rufus", "all", Trogon_rufus.heteroz))
colnames(Trogon_rufus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Trogon_rufus.df)

Trogon_rufus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Trogon_rufus.df)


#------------------------------------------------------------------------
# Xiphorhynchus_elegans
Xiphorhynchus_elegans <- read.structure("Xiphorhynchus_elegans_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                        n.ind=7, n.loc=7514, onerowperind=FALSE, 
                                        col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Xiphorhynchus_elegans.grp <- find.clusters(Xiphorhynchus_elegans, n.pca = 100, 
                                            method = "kmeans", choose.n.clust = FALSE, 
                                            criterion = "min")
Xiphorhynchus_elegans$pop <- Xiphorhynchus_elegans.grp$grp

Xiphorhynchus_elegans.div <- summary(Xiphorhynchus_elegans)
Xiphorhynchus_elegans.heteroz <- mean(Xiphorhynchus_elegans.div$Hobs)

Xiphorhynchus_elegans.div.1 <- summary(Xiphorhynchus_elegans[pop=1])
Xiphorhynchus_elegans.heteroz.1 <- mean(Xiphorhynchus_elegans.div.1$Hobs, na.rm = TRUE)

Xiphorhynchus_elegans.div.2 <- summary(Xiphorhynchus_elegans[pop=2])
Xiphorhynchus_elegans.heteroz.2 <- mean(Xiphorhynchus_elegans.div.2$Hobs, na.rm = TRUE)

Xiphorhynchus_elegans.df <- as.data.frame(list("Xiphorhynchus_elegans", "all", Xiphorhynchus_elegans.heteroz))
Xiphorhynchus_elegans.df.1 <- as.data.frame(list("Xiphorhynchus_elegans", "one", Xiphorhynchus_elegans.heteroz.1))
Xiphorhynchus_elegans.df.2 <- as.data.frame(list("Xiphorhynchus_elegans", "two", Xiphorhynchus_elegans.heteroz.2))

colnames(Xiphorhynchus_elegans.df) <- colnames(Xiphorhynchus_elegans.df.1) <- colnames(Xiphorhynchus_elegans.df.2) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Xiphorhynchus_elegans.df, Xiphorhynchus_elegans.df.1, Xiphorhynchus_elegans.df.2)



#------------------------------------------------------------------------
# Xiphorhynchus_obsoletus
Xiphorhynchus_obsoletus <- read.structure("Xiphorhynchus_obsoletus_SNPs_phased_rmIndels_75_QC_DP_seqcap_pop_structure.edited.str", 
                                          n.ind=7, n.loc=5307, onerowperind=FALSE, 
                                          col.lab=1, col.pop=0, col.others=0, row.marknames=0)

Xiphorhynchus_obsoletus.grp <- find.clusters(Xiphorhynchus_obsoletus, n.pca = 100, 
                                              method = "kmeans", choose.n.clust = FALSE, 
                                              criterion = "min")
Xiphorhynchus_obsoletus$pop <- Xiphorhynchus_obsoletus.grp$grp

Xiphorhynchus_obsoletus.div <- summary(Xiphorhynchus_obsoletus)
Xiphorhynchus_obsoletus.heteroz <- mean(Xiphorhynchus_obsoletus.div$Hobs)
Xiphorhynchus_obsoletus.df <- as.data.frame(list("Xiphorhynchus_obsoletus", "all", Xiphorhynchus_obsoletus.heteroz))
colnames(Xiphorhynchus_obsoletus.df) <- colnames(df.pop.het)

df.pop.het <- bind_rows(df.pop.het, Xiphorhynchus_obsoletus.df)

Xiphorhynchus_obsoletus.df[1,2] <- "one"
df.pop.het <- bind_rows(df.pop.het, Xiphorhynchus_obsoletus.df)



#------------------------------------------------------------------------

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/")
write.csv(df.pop.het, file = "3_results/per_population_heterozygosity.allSNPs.csv", row.names = FALSE)

# some plotting

temp <- df[,c(1:2,15)]
df.pop.het.t <- left_join(df.pop.het, temp, by = "species")

# do a regular gls, as we can't match species to tree tips
# because of multiple values per species
df.pop.het.tt <- df.pop.het %>% 
  filter(population != "all") %>% # restrict to per-population values
  filter(Hobs < 0.4) # remove two outliers with super low per-population sample size
df.pop.het.tt$habitat <- gsub("island", "0", df.pop.het.tt$habitat)
df.pop.het.tt$habitat <- gsub("floodplain", "1", df.pop.het.tt$habitat)
df.pop.het.tt$habitat <- gsub("upland", "2", df.pop.het.tt$habitat)
df.pop.het.tt$habitat <- as.numeric(df.pop.het.tt$habitat)

het.gls <- gls(habitat ~ Hobs, df.pop.het.tt)
summary(het.gls)

pval <- summary(het.gls)[18]$tTable[2,4]

label.het <- paste0("GLS\np = ", prettyNum(summary(het.gls)[18]$tTable[2,4], digits = 3),
                    ", t = ", prettyNum(summary(het.gls)[18]$tTable[2,3], digits = 3), 
                    "\nestimate = ", prettyNum(summary(het.gls)[4]$coefficients[[2]], digits = 3))

df.pop.het.t %>% 
  filter(population != "all") %>% # restrict to per-population values
  filter(Hobs < 0.4) %>% # remove two outliers with super low per-population sample size
  ggplot( aes(x=habitat, y=Hobs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.3)) +
  theme_classic() +
  scale_x_discrete(limits=c("island", "floodplain", "upland")) + 
  ggtitle("Per-population Heterozygosity\n") +
  labs(y = "Observed heterozygosity\n", x = "\nHabitat") +
  annotate("text", x=2.6, y=0.35, label = label.het, size = 3.6) 

ggsave("3_results/1_plots/per-population_heterozygosity.allSNPs.png")

df.pop.het.t %>% 
  filter(population == "all") %>% 
  filter(Hobs < 0.5) %>% # remove Elaenia pelzelni
  ggplot( aes(x=habitat, y=Hobs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.3)) +
  theme_classic() +
  scale_x_discrete(limits=c("island", "floodplain", "upland")) + 
  ggtitle("Per-species Heterozygosity\n") +
  labs(y = "Observed heterozygosity\n", x = "\nHabitat") 

ggsave("3_results/1_plots/per-species_heterozygosity.allSNPs.png")


