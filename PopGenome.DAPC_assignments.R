#!/usr/bin/env Rscript

library(PopGenome)
setwd("/Volumes/Brumfield_Lab_Drive/River_islands")

# each line break below separates species
# can't loop through each one because samples need to be assigned to 
# populations from DAPC analysis

# data frame to hold results
res.df <- data.frame(species=character(),
                     haplotype.F_ST=double(), 
                     nucleotide.F_ST=double(), 
                     Nei.G_ST=double(), 
                     Hudson.G_ST=double(), 
                     Hudson.H_ST=double(), 
                     Hudson.K_ST=double(), 
                     Dxy=double(), 
                     stringsAsFactors=FALSE) 

#------------------------------------------------------------------------
# Campephilus_melanoleucos

spec <- "Campephilus_melanoleucos"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
   c("Campephilus_melanoleucos_LSUMNS2802_0","Campephilus_melanoleucos_LSUMNS2802-1"),
   c("Campephilus_melanoleucos_LSUMNS26102_0", "Campephilus_melanoleucos_LSUMNS26102_1", "Campephilus_melanoleucos_MPEG16102_0",
     "Campephilus_melanoleucos_MPEG16102_1",   "Campephilus_melanoleucos_AMNH12719_0",   "Campephilus_melanoleucos_AMNH12719_1",
     "Campephilus_melanoleucos_FMNH433276_0",  "Campephilus_melanoleucos_FMNH433276_1",  "Campephilus_melanoleucos_FMNH456682_0", "Campephilus_melanoleucos_FMNH456682_1",
     "Campephilus_melanoleucos_LSUMNS33338_0", "Campephilus_melanoleucos_LSUMNS33338_1", "Campephilus_melanoleucos_MPEG11746_0",
     "Campephilus_melanoleucos_MPEG11746_1",   "Campephilus_melanoleucos_MPEG16666_0",   "Campephilus_melanoleucos_MPEG16666_1",
     "Campephilus_melanoleucos_MPEG7571_0",    "Campephilus_melanoleucos_MPEG7571_1",    "Campephilus_melanoleucos_USNMB22230_0",
     "Campephilus_melanoleucos_USNMB22230_1"
   )
   ))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Cantorchilus_leucotis

spec <- "Cantorchilus_leucotis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Cantorchilus_leucotis_FMNH389261_0", "Cantorchilus_leucotis_FMNH391580_0", "Cantorchilus_leucotis_FMNH392699_0", "Cantorchilus_leucotis_INPA2275_0", 
    "Cantorchilus_leucotis_LSUMNS35385_0", "Cantorchilus_leucotis_USNMB06896_0", "Cantorchilus_leucotis_FMNH389261_1", "Cantorchilus_leucotis_FMNH391580_1", 
    "Cantorchilus_leucotis_FMNH392699_1", "Cantorchilus_leucotis_INPA2275_1", "Cantorchilus_leucotis_LSUMNS35385_1", "Cantorchilus_leucotis_USNMB06896_1"
  ),
  c("Cantorchilus_leucotis_LSUMNS7343_0", "Cantorchilus_leucotis_LSUMNS9529_0", "Cantorchilus_leucotis_LSUMNS10899_0", "Cantorchilus_leucotis_LSUMNS46117_0", 
    "Cantorchilus_leucotis_LSUMNS7343_1", "Cantorchilus_leucotis_LSUMNS9529_1", "Cantorchilus_leucotis_LSUMNS10899_1", "Cantorchilus_leucotis_LSUMNS46117_1"
  )
))



# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Formicarius_analis

spec <- "Formicarius_analis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Formicarius_analis_AMNH11914_1", "Formicarius_analis_MPEG6767_1", "Formicarius_analis_INPA1349_1", "Formicarius_analis_AMNH11914_0", 
    "Formicarius_analis_MPEG6767_0", "Formicarius_analis_INPA1349_0"
  ),
  c("Formicarius_analis_FMNH433399_0", "Formicarius_analis_LSUMNS36683_0", "Formicarius_analis_LSUMNS4204_0", "Formicarius_analis_LSUMNS42862_0",
  "Formicarius_analis_LSUMNS81110_0", "Formicarius_analis_MPEG12249_0", "Formicarius_analis_MPEG12559_0", "Formicarius_analis_MPEG2999_0", 
  "Formicarius_analis_FMNH433399_1", "Formicarius_analis_LSUMNS36683_1", "Formicarius_analis_LSUMNS4204_1", "Formicarius_analis_LSUMNS42862_1", 
  "Formicarius_analis_LSUMNS81110_1", "Formicarius_analis_MPEG12249_1", "Formicarius_analis_MPEG12559_1", "Formicarius_analis_MPEG2999"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Formicarius_colma

spec <- "Formicarius_colma"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Formicarius_colma_MPEG9519_0", "Formicarius_colma_MPEG12536_0", "Formicarius_colma_LSUMNS13068_0", "Formicarius_colma_MPEG9519_1", 
    "Formicarius_colma_MPEG12536_1", "Formicarius_colma_LSUMNS13068"
  ),
  c("Formicarius_colma_AMNH11860_0", "Formicarius_colma_AMNH12722_0", "Formicarius_colma_INPA702_0", "Formicarius_colma_INPA7991_0", 
    "Formicarius_colma_LSUMNS42332_0", "Formicarius_colma_MPEG12382_0", "Formicarius_colma_MPEG328_0", "Formicarius_colma_MPEG6705_0",
    "Formicarius_colma_AMNH11860_1", "Formicarius_colma_AMNH12722_1", "Formicarius_colma_INPA702_1", "Formicarius_colma_INPA7991_1", 
    "Formicarius_colma_LSUMNS42332_1", "Formicarius_colma_MPEG12382_1", "Formicarius_colma_MPEG328_1", "Formicarius_colma_MPEG6705_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Glaucidium_brasilianum

spec <- "Glaucidium_brasilianum"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Glaucidium_brasilianum_AMNH2937_0", "Glaucidium_brasilianum_AMNH2937_1"
  ),
  c("Glaucidium_brasilianum_LSUMNS15361_0", "Glaucidium_brasilianum_LSUMNS15361_1"
  ),
  c("Glaucidium_brasilianum_FMNH456484_1", "Glaucidium_brasilianum_LSUMNS75626_1", "Glaucidium_brasilianum_LSUMNS980_1", 
    "Glaucidium_brasilianum_MPEG13997_1", "Glaucidium_brasilianum_MPEG14295_1", "Glaucidium_brasilianum_MPEG16187_1", 
    "Glaucidium_brasilianum_FMNH456484_0", "Glaucidium_brasilianum_LSUMNS75626_0", "Glaucidium_brasilianum_LSUMNS980_0", 
    "Glaucidium_brasilianum_MPEG13997_0", "Glaucidium_brasilianum_MPEG14295_0", "Glaucidium_brasilianum_MPEG16187_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Megascops_choliba

spec <- "Megascops_choliba"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Megascops_choliba_MPEG17491_0", "Megascops_choliba_MPEG16711_0", "Megascops_choliba_LSUMNS15318_0", "Megascops_choliba_MPEG16509_0", 
    "Megascops_choliba_LSUMNS9593_0", "Megascops_choliba_LSUMNS38217_0", "Megascops_choliba_MPEG17491_1", "Megascops_choliba_MPEG16711_1", 
    "Megascops_choliba_LSUMNS15318_1", "Megascops_choliba_MPEG16509_1", "Megascops_choliba_LSUMNS9593_1", "Megascops_choliba_LSUMNS38217_1"
  ),
  c("Megascops_choliba_AMNH4811_0", "Megascops_choliba_FMNH392673_0", "Megascops_choliba_LSUMNS42284_0", "Megascops_choliba_LSUMNS7420_0", 
    "Megascops_choliba_MPEG13970_0", "Megascops_choliba_AMNH4811_1", "Megascops_choliba_FMNH392673_1", "Megascops_choliba_LSUMNS42284_1", 
    "Megascops_choliba_LSUMNS7420_1", "Megascops_choliba_MPEG13970_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Monasa_nigrifrons

spec <- "Monasa_nigrifrons"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Monasa_nigrifrons_MPEG10345_1", "Monasa_nigrifrons_MPEG11481_1", "Monasa_nigrifrons_MPEG2312_1", "Monasa_nigrifrons_MPEG17327_1", 
    "Monasa_nigrifrons_USNMB07074_1", "Monasa_nigrifrons_MPEG10345_0", "Monasa_nigrifrons_MPEG11481_0", "Monasa_nigrifrons_MPEG2312_0", 
    "Monasa_nigrifrons_MPEG17327_0", "Monasa_nigrifrons_USNMB07074_0"
  ),
  c("Monasa_nigrifrons_FMNH321039_1", "Monasa_nigrifrons_FMNH456622_1", "Monasa_nigrifrons_INPA5607_1", "Monasa_nigrifrons_LSUMNS42834_1", 
    "Monasa_nigrifrons_LSUMNS5042_1", "Monasa_nigrifrons_MPEG16028_1", "Monasa_nigrifrons_FMNH321039_0", "Monasa_nigrifrons_FMNH456622_0", 
    "Monasa_nigrifrons_INPA5607_0", "Monasa_nigrifrons_LSUMNS42834_0", "Monasa_nigrifrons_LSUMNS5042_0", "Monasa_nigrifrons_MPEG16028_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Pheugopedius_coraya

spec <- "Pheugopedius_coraya"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Pheugopedius_coraya_LSUMNS4133_0", "Pheugopedius_coraya_AMNH4240_0", "Pheugopedius_coraya_MPEG9534_0", "Pheugopedius_coraya_KU16858_0", 
    "Pheugopedius_coraya_MPEG12256_0", "Pheugopedius_coraya_LSUMNS4133_1", "Pheugopedius_coraya_AMNH4240_1", "Pheugopedius_coraya_MPEG9534_1", 
    "Pheugopedius_coraya_KU16858_1", "Pheugopedius_coraya_MPEG12256_1"
  ),
  c("Pheugopedius_coraya_LSUMNS46172_0", "Pheugopedius_genibarbis_INPA1333_0", "Pheugopedius_genibarbis_LSUMNS58336_0", 
    "Pheugopedius_genibarbis_MPEG10477_0", "Pheugopedius_genibarbis_MPEG12881_0", "Pheugopedius_genibarbis_MPEG9500_0", 
    "Pheugopedius_coraya_LSUMNS46172_1", "Pheugopedius_genibarbis_INPA1333_1", "Pheugopedius_genibarbis_LSUMNS58336_1", 
    "Pheugopedius_genibarbis_MPEG10477_1", "Pheugopedius_genibarbis_MPEG12881_1", "Pheugopedius_genibarbis_MPEG9500_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Piaya_cayana

spec <- "Piaya_cayana"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Piaya_cayana_USNMB13940_0", "Piaya_cayana_USNMB13940_1"
  ),
  c("Piaya_cayana_FMNH392671_0", "Piaya_cayana_KU16728_0", "Piaya_cayana_LSUMNS31482_0", "Piaya_cayana_LSUMNS35624_0", 
    "Piaya_cayana_LSUMNS44417_0", "Piaya_cayana_LSUMNS75766_0", "Piaya_cayana_MPEG1237_0", "Piaya_cayana_MPEG2161_0", 
    "Piaya_cayana_MPEG3662_0", "Piaya_cayana_FMNH392671_1", "Piaya_cayana_KU16728_1", "Piaya_cayana_LSUMNS31482_1", 
    "Piaya_cayana_LSUMNS35624_1", "Piaya_cayana_LSUMNS44417_1", "Piaya_cayana_LSUMNS75766_1", "Piaya_cayana_MPEG1237_1", 
    "Piaya_cayana_MPEG2161_1", "Piaya_cayana_MPEG3662_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Pipra_erythrocephala

spec <- "Pipra_erythrocephala"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Pipra_chloromeros_KU18533_1", "Pipra_chloromeros_LSUMNS106768_1", "Pipra_chloromeros_KU18533_0", "Pipra_chloromeros_LSUMNS106768_0"
  ),
  c("Pipra_erythrocephala_AMNH14356_0", "Pipra_erythrocephala_LSUMNS3092_0", "Pipra_erythrocephala_MPEG6699_0", 
    "Pipra_erythrocephala_MVZ165267_0", "Pipra_rubrocapilla_FMNH457363_0", "Pipra_rubrocapilla_MPEG14204_0", 
    "Pipra_rubrocapilla_MPEG14850_0", "Pipra_rubrocapilla_MPEG7551_0", "Pipra_rubrocapilla_MPEG854_0", "Pipra_erythrocephala_AMNH14356_1", 
    "Pipra_erythrocephala_LSUMNS3092_1", "Pipra_erythrocephala_MPEG6699_1", "Pipra_erythrocephala_MVZ165267_1", 
    "Pipra_rubrocapilla_FMNH457363_1", "Pipra_rubrocapilla_MPEG14204_1", "Pipra_rubrocapilla_MPEG14850_1", "Pipra_rubrocapilla_MPEG7551_1", 
    "Pipra_rubrocapilla_MPEG854_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Saltator_coerulescens

spec <- "Saltator_coerulescens"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Saltator_coerulescens_USNMB11140_1", "Saltator_coerulescens_INPA2123_1", "Saltator_coerulescens_USNMB11140_0", "Saltator_coerulescens_INPA2123"
  ),
  c("Saltator_coerulescens_FMNH391613_0", "Saltator_coerulescens_INPA6946_0", "Saltator_coerulescens_LSUMNS6816_0", 
    "Saltator_coerulescens_LSUMNS7298_0", "Saltator_coerulescens_LSUMNS75970_0", "Saltator_coerulescens_MPEG15265_0", 
    "Saltator_coerulescens_MPEG16382_0", "Saltator_coerulescens_MPEG3969_0", "Saltator_coerulescens_MPEG5844_0", 
    "Saltator_coerulescens_FMNH391613_1", "Saltator_coerulescens_INPA6946_1", "Saltator_coerulescens_LSUMNS6816_1", 
    "Saltator_coerulescens_LSUMNS7298_1", "Saltator_coerulescens_LSUMNS75970_1", "Saltator_coerulescens_MPEG15265_1", 
    "Saltator_coerulescens_MPEG16382_1", "Saltator_coerulescens_MPEG3969_1", "Saltator_coerulescens_MPEG5844_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Saltator_grossus

spec <- "Saltator_grossus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Saltator_grossus_USNMB09368_0", "Saltator_grossus_FMNH433813_0", "Saltator_grossus_USNMB09368_1", "Saltator_grossus_FMNH433813"
  ),
  c("Saltator_grossus_AMNH12692_1", "Saltator_grossus_FMNH457551_1", "Saltator_grossus_LSUMNS76026_1", "Saltator_grossus_MPEG10988_1", 
    "Saltator_grossus_MPEG11309_1", "Saltator_grossus_MPEG2492_1", "Saltator_grossus_MPEG5696_1", "Saltator_grossus_MPEG6513_1", 
    "Saltator_grossus_MPEG7434_1", "Saltator_grossus_AMNH12692_0", "Saltator_grossus_FMNH457551_0", "Saltator_grossus_LSUMNS76026_0", 
    "Saltator_grossus_MPEG10988_0", "Saltator_grossus_MPEG11309_0", "Saltator_grossus_MPEG2492_0", "Saltator_grossus_MPEG5696_0", 
    "Saltator_grossus_MPEG6513_0", "Saltator_grossus_MPEG7434_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Schiffornis_turdina

spec <- "Schiffornis_turdina"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Schiffornis_turdina_LSUMNS27903_1", "Schiffornis_turdina_LSUMNS27903_0"
  ),
  c("Schiffornis_turdina_AMNH4312_0", "Schiffornis_turdina_INPA725_0", "Schiffornis_turdina_KU18445_0", "Schiffornis_turdina_MPEG10904_0", 
    "Schiffornis_turdina_MPEG11026_0", "Schiffornis_turdina_MPEG12581_0", "Schiffornis_turdina_MPEG2157_0", "Schiffornis_turdina_MPEG2281_0", 
    "Schiffornis_turdina_MPEG3417_0", "Schiffornis_turdina_MPEG6768_0", "Schiffornis_turdina_AMNH4312_1", "Schiffornis_turdina_INPA725_1", 
    "Schiffornis_turdina_KU18445_1", "Schiffornis_turdina_MPEG10904_1", "Schiffornis_turdina_MPEG11026_1", "Schiffornis_turdina_MPEG12581_1", 
    "Schiffornis_turdina_MPEG2157_1", "Schiffornis_turdina_MPEG2281_1", "Schiffornis_turdina_MPEG3417_1", "Schiffornis_turdina_MPEG6768_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Tachyphonus_cristatus

spec <- "Tachyphonus_cristatus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Tachyphonus_cristatus_FMNH390058_0", "Tachyphonus_cristatus_INPA1490_0", "Tachyphonus_cristatus_LSUMNS35475_0", 
    "Tachyphonus_cristatus_LSUMNS9548_0", "Tachyphonus_cristatus_MPEG2482_0", "Tachyphonus_cristatus_FMNH390058_1", 
    "Tachyphonus_cristatus_INPA1490_1", "Tachyphonus_cristatus_LSUMNS35475_1", "Tachyphonus_cristatus_LSUMNS9548_1", 
    "Tachyphonus_cristatus_MPEG2482_1"
  ),
  c("Tachyphonus_cristatus_LSUMNS2693_1", "Tachyphonus_cristatus_INPA1636_1", "Tachyphonus_cristatus_AMNH2985_1", 
    "Tachyphonus_cristatus_LSUMNS2693_0", "Tachyphonus_cristatus_INPA1636_0", "Tachyphonus_cristatus_AMNH2985_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Tachyphonus_luctuosus

spec <- "Tachyphonus_luctuosus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Tachyphonus_luctuosus_LSUMNS18277_1", "Tachyphonus_luctuosus_LSUMNS27928_1", "Tachyphonus_luctuosus_MPEG10610_1", 
    "Tachyphonus_luctuosus_MPEG2517_1", "Tachyphonus_luctuosus_MPEG4059_1", "Tachyphonus_luctuosus_MSUZP95518_1", 
    "Tachyphonus_luctuosus_MZUSP87733_1", "Tachyphonus_luctuosus_USNMB09569_1", "Tachyphonus_luctuosus_USNMB10690_1", 
    "Tachyphonus_luctuosus_LSUMNS18277_0", "Tachyphonus_luctuosus_LSUMNS27928_0", "Tachyphonus_luctuosus_MPEG10610_0", 
    "Tachyphonus_luctuosus_MPEG2517_0", "Tachyphonus_luctuosus_MPEG4059_0", "Tachyphonus_luctuosus_MSUZP95518_0", 
    "Tachyphonus_luctuosus_MZUSP87733_0", "Tachyphonus_luctuosus_USNMB09569_0", "Tachyphonus_luctuosus_USNMB10690_0"
  ),
  c("Tachyphonus_luctuosus_MVZ169546_1", "Tachyphonus_luctuosus_MPEG4450_1", "Tachyphonus_luctuosus_MVZ169546_0", "Tachyphonus_luctuosus_MPEG4450_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Trogon_collaris

spec <- "Trogon_collaris"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Trogon_collaris_MPEG8909_0", "Trogon_collaris_USNMB22128_0", "Trogon_collaris_MPEG8909_1", "Trogon_collaris_USNMB22128_1"
  ),
  c("Trogon_collaris_FMNH456552_0", "Trogon_collaris_INPA6982_0", "Trogon_collaris_LSUMNS10657_0", "Trogon_collaris_LSUMNS18342_0", 
    "Trogon_collaris_LSUMNS22827_0", "Trogon_collaris_LSUMNS27866_0", "Trogon_collaris_MPEG17320_0", "Trogon_collaris_MPEG3657_0", 
    "Trogon_collaris_FMNH456552_1", "Trogon_collaris_INPA6982_1", "Trogon_collaris_LSUMNS10657_1", "Trogon_collaris_LSUMNS18342_1", 
    "Trogon_collaris_LSUMNS22827_1", "Trogon_collaris_LSUMNS27866_1", "Trogon_collaris_MPEG17320_1", "Trogon_collaris_MPEG3657_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Trogon_rufus

spec <- "Trogon_rufus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-AmazonOnly")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Trogon_rufus_FMNH389730_1", "Trogon_rufus_FMNH456559_1", "Trogon_rufus_INPA5284_1", "Trogon_rufus_MPEG16761_1", 
    "Trogon_rufus_MPEG16986_1", "Trogon_rufus_MPEG1878_1", "Trogon_rufus_FMNH389730_0", "Trogon_rufus_FMNH456559_0", 
    "Trogon_rufus_INPA5284_0", "Trogon_rufus_MPEG16761_0", "Trogon_rufus_MPEG16986_0", "Trogon_rufus_MPEG1878_0"
  ),
  c("Trogon_rufus_LSUMNS4256_1", "Trogon_rufus_INPA1170_1", "Trogon_rufus_INPA1668_1", "Trogon_rufus_MPEG6808_1", 
    "Trogon_rufus_LSUMNS27570_1", "Trogon_rufus_LSUMNS4256_0", "Trogon_rufus_INPA1170_0", "Trogon_rufus_INPA1668_0", 
    "Trogon_rufus_MPEG6808_0", "Trogon_rufus_LSUMNS27570_0"

  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Campephilus_rubricollis

spec <- "Campephilus_rubricollis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Campephilus_rubricollis_FMNH456685_0", "Campephilus_rubricollis_LSUMNS15139_0", "Campephilus_rubricollis_LSUMNS39164_0", 
    "Campephilus_rubricollis_LSUMNS42486_0", "Campephilus_rubricollis_LSUMNS5075_0", "Campephilus_rubricollis_LSUMNS75519_0", 
    "Campephilus_rubricollis_LSUMNS80981_0", "Campephilus_rubricollis_MPEG3577_0", "Campephilus_rubricollis_FMNH456685_1", 
    "Campephilus_rubricollis_LSUMNS15139_1", "Campephilus_rubricollis_LSUMNS39164_1", "Campephilus_rubricollis_LSUMNS42486_1", 
    "Campephilus_rubricollis_LSUMNS5075_1", "Campephilus_rubricollis_LSUMNS75519_1", "Campephilus_rubricollis_LSUMNS80981_1", 
    "Campephilus_rubricollis_MPEG3577_1"
  ),
  c("Campephilus_rubricollis_MPEG7906_0", "Campephilus_rubricollis_LSUMNS48392_0", "Campephilus_rubricollis_AMNH3870_0", 
    "Campephilus_rubricollis_MPEG7906_1", "Campephilus_rubricollis_LSUMNS48392_1", "Campephilus_rubricollis_AMNH3870_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Celeus_flavus

spec <- "Celeus_flavus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Celeus_flavus_MPEG9460_1", "Celeus_flavus_USNMB06880_1", "Celeus_flavus_LSUMNS77905_1", "Celeus_flavus_KU5738_1", 
    "Celeus_flavus_MPEG9460_0", "Celeus_flavus_USNMB06880_0", "Celeus_flavus_LSUMNS77905_0", "Celeus_flavus_KU5738_0"
  ),
  c("Celeus_flavus_KU675_1", "Celeus_flavus_LSUMNS34966_1", "Celeus_flavus_LSUMNS4385_1", "Celeus_flavus_LSUMNS75724_1", 
    "Celeus_flavus_MPEG3876_1", "Celeus_flavus_KU675_0", "Celeus_flavus_LSUMNS34966_0", "Celeus_flavus_LSUMNS4385_0", 
    "Celeus_flavus_LSUMNS75724_0", "Celeus_flavus_MPEG3876_0"

  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Celeus_grammicus

spec <- "Celeus_grammicus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Celeus_undatus_USNMB05167_0", "Celeus_undatus_USNMB11945_0", "Celeus_undatus_USNMB05167_1", "Celeus_undatus_USNMB11945_1"
  ),
  c("Celeus_grammicus_AMNH3861_0", "Celeus_grammicus_FMNH389781_0", "Celeus_grammicus_LSUMNS23798_0", "Celeus_grammicus_LSUMNS42309_0", 
    "Celeus_grammicus_MPEG10819_0", "Celeus_grammicus_MPEG10847_0", "Celeus_grammicus_MPEG13886_0", "Celeus_grammicus_MPEG14455_0", 
    "Celeus_grammicus_MPEG2850_0", "Celeus_grammicus_AMNH3861_1", "Celeus_grammicus_FMNH389781_1", "Celeus_grammicus_LSUMNS23798_1", 
    "Celeus_grammicus_LSUMNS42309_1", "Celeus_grammicus_MPEG10819_1", "Celeus_grammicus_MPEG10847_1", "Celeus_grammicus_MPEG13886_1", 
    "Celeus_grammicus_MPEG14455_1", "Celeus_grammicus_MPEG2850_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Crypturellus_undulatus

spec <- "Crypturellus_undulatus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Crypturellus_undulatus_MZUSP79467_1", "Crypturellus_undulatus_MZUSP94561_1", "Crypturellus_undulatus_MZUSP79467_0", "Crypturellus_undulatus_MZUSP94561_0"
  ),
  c("Crypturellus_undulatus_MPEG7863_1", "Crypturellus_undulatus_MPEG14258_1", "Crypturellus_undulatus_MPEG7863_0", "Crypturellus_undulatus_MPEG14258_0"
  ),
  c("Crypturellus_undulatus_AMNH2312_1", "Crypturellus_undulatus_FMNH395852_1", "Crypturellus_undulatus_MPEG15255_1", 
    "Crypturellus_undulatus_MPEG7713_1", "Crypturellus_undulatus_YPM136994_1", "Crypturellus_undulatus_AMNH2312_0", 
    "Crypturellus_undulatus_FMNH395852_0", "Crypturellus_undulatus_MPEG15255_0", "Crypturellus_undulatus_MPEG7713_0", 
    "Crypturellus_undulatus_YPM136994_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)


#------------------------------------------------------------------------
# Crypturellus_variegatus

spec <- "Crypturellus_variegatus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Crypturellus_variegatus_LSUMNS36774_0", "Crypturellus_variegatus_LSUMNS40657_0", "Crypturellus_variegatus_LSUMNS76609_0", 
    "Crypturellus_variegatus_LSUMNS9460_0", "Crypturellus_variegatus_MPEG14165_0", "Crypturellus_variegatus_USNMB07094_0", 
    "Crypturellus_variegatus_LSUMNS36774_1", "Crypturellus_variegatus_LSUMNS40657_1", "Crypturellus_variegatus_LSUMNS76609_1", 
    "Crypturellus_variegatus_LSUMNS9460_1", "Crypturellus_variegatus_MPEG14165_1", "Crypturellus_variegatus_USNMB07094_1"
  ),
  c("Crypturellus_variegatus_MPEG8335_0", "Crypturellus_variegatus_LSUMNS20328_0", "Crypturellus_variegatus_USNMB09220_0", 
    "Crypturellus_variegatus_INPA131_0", "Crypturellus_variegatus_MPEG8335_1", "Crypturellus_variegatus_LSUMNS20328_1", 
    "Crypturellus_variegatus_USNMB09220_1", "Crypturellus_variegatus_INPA131_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Glaucidium_hardyi

spec <- "Glaucidium_hardyi"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Glaucidium_hardyi_FMNH389707_1", "Glaucidium_hardyi_KU476_1", "Glaucidium_hardyi_LSUMNS31382_1", "Glaucidium_hardyi_LSUMNS40662_1", 
    "Glaucidium_hardyi_MPEG11657_1", "Glaucidium_hardyi_MPEG1762_1", "Glaucidium_hardyi_MPEG694_1", "Glaucidium_hardyi_FMNH389707_0", 
    "Glaucidium_hardyi_KU476_0", "Glaucidium_hardyi_LSUMNS31382_0", "Glaucidium_hardyi_LSUMNS40662_0", "Glaucidium_hardyi_MPEG11657_0", 
    "Glaucidium_hardyi_MPEG1762_0", "Glaucidium_hardyi_MPEG694_0"
  ),
  c("Glaucidium_hardyi_MPEG14844_1", "Glaucidium_hardyi_USNMB12705_1", "Glaucidium_hardyi_LSUMNS20184_1", "Glaucidium_hardyi_MPEG14844_0", 
    "Glaucidium_hardyi_USNMB12705_0", "Glaucidium_hardyi_LSUMNS20184_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Hylophylax_naevia

spec <- "Hylophylax_naevia"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Hylophylax_naevia_AMNH12315_0", "Hylophylax_naevia_INPA1301_0", "Hylophylax_naevia_LSUMNS5972_0", "Hylophylax_naevia_USNMB22064_0", 
    "Hylophylax_naevia_AMNH12315_1", "Hylophylax_naevia_INPA1301_1", "Hylophylax_naevia_LSUMNS5972_1", "Hylophylax_naevia_USNMB22064_1"
  ),
  c("Hylophylax_naevia_FMNH457166_0", "Hylophylax_naevia_MPEG12164_0", "Hylophylax_naevia_MPEG10842_0", "Hylophylax_naevia_FMNH457166_1", 
    "Hylophylax_naevia_MPEG12164_1", "Hylophylax_naevia_MPEG10842_1"
  ),
  c("Hylophylax_naevia_MPEG828_0", "Hylophylax_naevia_LSUMNS106776_0", "Hylophylax_naevia_MPEG828_1", "Hylophylax_naevia_LSUMNS106776_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Hylophylax_punctulata

spec <- "Hylophylax_punctulata"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class,list(
  c("Hylophylax_punctulata_FMNH457204_1", "Hylophylax_punctulata_INPA5623_1", "Hylophylax_punctulata_LSUMNS42951_1", 
    "Hylophylax_punctulata_LSUMNS4399_1", "Hylophylax_punctulata_MPEG3422_1", "Hylophylax_punctulata_FMNH457204_0", 
    "Hylophylax_punctulata_INPA5623_0", "Hylophylax_punctulata_LSUMNS42951_0", "Hylophylax_punctulata_LSUMNS4399_0", 
    "Hylophylax_punctulata_MPEG3422_0"
  ),
  c("Hylophylax_punctulata_LSUMNS12357_1", "Hylophylax_punctulata_LSUMNS81100_1", "Hylophylax_punctulata_USNMB06906_1", 
    "Hylophylax_punctulata_MPEG12824_1", "Hylophylax_punctulata_LSUMNS12357_0", "Hylophylax_punctulata_LSUMNS81100_0", 
    "Hylophylax_punctulata_USNMB06906_0", "Hylophylax_punctulata_MPEG12824_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Megascops_watsonii

spec <- "Megascops_watsonii"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Megascops_watsonii_AMNH3868_0", "Megascops_watsonii_LSUMNS42986_0", "Megascops_watsonii_LSUMNS4521_0", "Megascops_watsonii_LSUMNS75828_0", 
    "Megascops_watsonii_LSUMNS80997_0", "Megascops_watsonii_LSUMNS947_0", "Megascops_watsonii_MPEG5114_0", "Megascops_watsonii_MPEG632_0", 
    "Megascops_watsonii_USNMB06991_0", "Megascops_watsonii_AMNH3868_1", "Megascops_watsonii_LSUMNS42986_1", "Megascops_watsonii_LSUMNS4521_1", 
    "Megascops_watsonii_LSUMNS75828_1", "Megascops_watsonii_LSUMNS80997_1", "Megascops_watsonii_LSUMNS947_1", "Megascops_watsonii_MPEG5114_1", 
    "Megascops_watsonii_MPEG632_1", "Megascops_watsonii_USNMB06991_1"
  ),
  c("Megascops_watsonii_LSUMNS20185_0", "Megascops_watsonii_MPEG9198_0", "Megascops_watsonii_LSUMNS20185_1", "Megascops_watsonii_MPEG9198_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Monasa_morphoeus

spec <- "Monasa_morphoeus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Monasa_morphoeus_LSUMNS35587_0", "Monasa_morphoeus_LSUMNS35587_1"
  ),
  c("Monasa_morphoeus_MPEG4939_1", "Monasa_morphoeus_INPA3002_1", "Monasa_morphoeus_MPEG4939_0", "Monasa_morphoeus_INPA3002_0"
  ),
  c("Monasa_morphoeus_AMNH12750_1", "Monasa_morphoeus_LSUMNS2686_1", "Monasa_morphoeus_LSUMNS28006_1", "Monasa_morphoeus_LSUMNS42358_1", 
    "Monasa_morphoeus_AMNH12750_0", "Monasa_morphoeus_LSUMNS2686_0", "Monasa_morphoeus_LSUMNS28006_0", "Monasa_morphoeus_LSUMNS42358_0"
  ),
  c("Monasa_morphoeus_MPEG1432_1", "Monasa_morphoeus_MPEG2245_1", "Monasa_morphoeus_MPEG2423_1", "Monasa_morphoeus_MPEG1432_0", 
    "Monasa_morphoeus_MPEG2245_0", "Monasa_morphoeus_MPEG2423_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmeciza_fortis

spec <- "Myrmeciza_fortis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Myrmeciza_fortis_LSUMNS27685_0", "Myrmeciza_fortis_LSUMNS43063_0", "Myrmeciza_fortis_LSUMNS4299_0", 
    "Myrmeciza_fortis_FMNH457129_0", "Myrmeciza_fortis_LSUMNS27685_1", "Myrmeciza_fortis_LSUMNS43063_1", 
    "Myrmeciza_fortis_LSUMNS4299_1", "Myrmeciza_fortis_FMNH457129_1"
  ),
  c("Myrmeciza_fortis_INPA431_0", "Myrmeciza_fortis_LSUMNS11047_0", "Myrmeciza_fortis_LSUMNS78791_0", 
    "Myrmeciza_fortis_LSUMNS9094_0", "Myrmeciza_fortis_MPEG3804_0", "Myrmeciza_fortis_MPEG3868_0", 
    "Myrmeciza_fortis_MPEG5024_0", "Myrmeciza_fortis_INPA431_1", "Myrmeciza_fortis_LSUMNS11047_1", 
    "Myrmeciza_fortis_LSUMNS78791_1", "Myrmeciza_fortis_LSUMNS9094_1", "Myrmeciza_fortis_MPEG3804_1", 
    "Myrmeciza_fortis_MPEG3868_1", "Myrmeciza_fortis_MPEG5024_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmeciza_hyperythra

spec <- "Myrmeciza_hyperythra"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Myrmeciza_hyperythra_FMNH433423_0", "Myrmeciza_hyperythra_INPA31_0", "Myrmeciza_hyperythra_KU439_0", 
    "Myrmeciza_hyperythra_LSUMNS35684_0", "Myrmeciza_hyperythra_LSUMNS75668_0", "Myrmeciza_hyperythra_MPEG13975_0", 
    "Myrmeciza_hyperythra_MPEG3395_0", "Myrmeciza_hyperythra_MPEG4370_0", "Myrmeciza_hyperythra_FMNH433423_1", 
    "Myrmeciza_hyperythra_INPA31_1", "Myrmeciza_hyperythra_KU439_1", "Myrmeciza_hyperythra_LSUMNS35684_1", 
    "Myrmeciza_hyperythra_LSUMNS75668_1", "Myrmeciza_hyperythra_MPEG13975_1", "Myrmeciza_hyperythra_MPEG3395_1", 
    "Myrmeciza_hyperythra_MPEG4370_1"
  ),
  c("Myrmeciza_hyperythra_LSUMNS42497_1", "Myrmeciza_hyperythra_LSUMNS7342_1", "Myrmeciza_hyperythra_LSUMNS27370_1", 
    "Myrmeciza_hyperythra_LSUMNS42497_0", "Myrmeciza_hyperythra_LSUMNS7342_0", "Myrmeciza_hyperythra_LSUMNS27370_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmoborus_leucophrys

spec <- "Myrmoborus_leucophrys"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Myrmoborus_leucophrys_MPEG11444_0", "Myrmoborus_leucophrys_MPEG1868_0", "Myrmoborus_leucophrys_MPEG16732_0", 
    "Myrmoborus_leucophrys_MPEG11444_1", "Myrmoborus_leucophrys_MPEG1868_1", "Myrmoborus_leucophrys_MPEG16732_1"
  ),
  c("Myrmoborus_leucophrys_LSUMNS10883_0", "Myrmoborus_leucophrys_LSUMNS46168_0", "Myrmoborus_leucophrys_MPEG16066_0", 
    "Myrmoborus_leucophrys_LSUMNS1005_0", "Myrmoborus_leucophrys_INPA32_0", "Myrmoborus_leucophrys_KU18514_0", 
    "Myrmoborus_leucophrys_LSUMNS10883_1", "Myrmoborus_leucophrys_LSUMNS46168_1", "Myrmoborus_leucophrys_MPEG16066_1", 
    "Myrmoborus_leucophrys_LSUMNS1005_1", "Myrmoborus_leucophrys_INPA32_1", "Myrmoborus_leucophrys_KU18514_1"
  ),
  c("Myrmoborus_leucophrys_MPEG6752_0", "Myrmoborus_leucophrys_USNMB12301_0", 
    "Myrmoborus_leucophrys_MPEG6752_1", "Myrmoborus_leucophrys_USNMB12301_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmoborus_myotherinus

spec <- "Myrmoborus_myotherinus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Myrmoborus_myotherinus_MPEG10893_1", "Myrmoborus_myotherinus_MPEG9173_1", "Myrmoborus_myotherinus_MPEG9511_1", 
    "Myrmoborus_myotherinus_MPEG4187_1", "Myrmoborus_myotherinus_MPEG10893_0", "Myrmoborus_myotherinus_MPEG9173_0", 
    "Myrmoborus_myotherinus_MPEG9511_0", "Myrmoborus_myotherinus_MPEG4187_0"
  ),
  c("Myrmoborus_myotherinus_LSUMNS1039_1", "Myrmoborus_myotherinus_LSUMNS42235_1", "Myrmoborus_myotherinus_MPEG3884_1", 
    "Myrmoborus_myotherinus_LSUMNS1039_0", "Myrmoborus_myotherinus_LSUMNS42235_0", "Myrmoborus_myotherinus_MPEG3884_0"
  ),
  c("Myrmoborus_myotherinus_MPEG390_1", "Myrmoborus_myotherinus_MPEG390_0"
  ),
  c("Myrmoborus_myotherinus_LSUMNS4345_1", "Myrmoborus_myotherinus_AMNH8777_1", "Myrmoborus_myotherinus_AMNH8837_1", 
    "Myrmoborus_myotherinus_LSUMNS4345_0", "Myrmoborus_myotherinus_AMNH8777_0", "Myrmoborus_myotherinus_AMNH8837_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Phaethornis_bourcieri

spec <- "Phaethornis_bourcieri"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Phaethornis_bourcieri_AMNH11853_0", "Phaethornis_bourcieri_AMNH11853_1"
  ),
  c("Phaethornis_bourcieri_AMNH14368_0", "Phaethornis_bourcieri_LSUMNS27643_0", "Phaethornis_bourcieri_LSUMNS4314_0", 
    "Phaethornis_bourcieri_LSUMNS75469_0", "Phaethornis_bourcieri_AMNH14368_1", "Phaethornis_bourcieri_LSUMNS27643_1", 
    "Phaethornis_bourcieri_LSUMNS4314_1", "Phaethornis_bourcieri_LSUMNS75469_1"
  ),
  c("Phaethornis_philippii_LSUMNS9442_0", "Phaethornis_philippii_MPEG3369_0", "Phaethornis_philippii_MPEG386_0", 
    "Phaethornis_philippii_LSUMNS9442_1", "Phaethornis_philippii_MPEG3369_1", "Phaethornis_philippii_MPEG386_1"
  ),
  c("Phaethornis_bourcieri_MPEG10930_0", "Phaethornis_bourcieri_MPEG1805_0", 
    "Phaethornis_bourcieri_MPEG10930_1", "Phaethornis_bourcieri_MPEG1805_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Phaethornis_hispidus

spec <- "Phaethornis_hispidus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Phaethornis_hispidus_MPEG13651_1", "Phaethornis_hispidus_MPEG10447_1", "Phaethornis_hispidus_MPEG15892_1", 
    "Phaethornis_hispidus_MPEG16385_1", "Phaethornis_hispidus_MPEG13651_0", "Phaethornis_hispidus_MPEG10447_0", 
    "Phaethornis_hispidus_MPEG15892_0", "Phaethornis_hispidus_MPEG16385_0"
  ),
  c("Phaethornis_hispidus_FMNH433168_1", "Phaethornis_hispidus_LSUMNS103606_1", "Phaethornis_hispidus_LSUMNS10677_1", 
    "Phaethornis_hispidus_LSUMNS42786_1", "Phaethornis_hispidus_LSUMNS6802_1", "Phaethornis_hispidus_MPEG13967_1", 
    "Phaethornis_hispidus_MPEG3535_1", "Phaethornis_hispidus_FMNH433168_0", "Phaethornis_hispidus_LSUMNS103606_0", 
    "Phaethornis_hispidus_LSUMNS10677_0", "Phaethornis_hispidus_LSUMNS42786_0", "Phaethornis_hispidus_LSUMNS6802_0", 
    "Phaethornis_hispidus_MPEG13967_0", "Phaethornis_hispidus_MPEG3535_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Piaya_melanogaster

spec <- "Piaya_melanogaster"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Piaya_melanogaster_INPA1392_0", "Piaya_melanogaster_LSUMNS36684_0", "Piaya_melanogaster_LSUMNS42540_0", 
    "Piaya_melanogaster_LSUMNS76234_0", "Piaya_melanogaster_LSUMNS78489_0", "Piaya_melanogaster_MPEG2461_0", 
    "Piaya_melanogaster_MPEG3831_0", "Piaya_melanogaster_MPEG5177_0", "Piaya_melanogaster_USNMB14512_0", 
    "Piaya_melanogaster_INPA1392_1", "Piaya_melanogaster_LSUMNS36684_1", "Piaya_melanogaster_LSUMNS42540_1", 
    "Piaya_melanogaster_LSUMNS76234_1", "Piaya_melanogaster_LSUMNS78489_1", "Piaya_melanogaster_MPEG2461_1", 
    "Piaya_melanogaster_MPEG3831_1", "Piaya_melanogaster_MPEG5177_1", "Piaya_melanogaster_USNMB14512_1"
  ),
  c("Piaya_melanogaster_MPEG9244_0", "Piaya_melanogaster_KU3921_0", 
    "Piaya_melanogaster_MPEG9244_1", "Piaya_melanogaster_KU3921_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Pipra_filicauda

spec <- "Pipra_filicauda"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Pipra_aureola_MPEG11971_1", "Pipra_aureola_MPEG1199_1", "Pipra_aureola_MPEG11971_0", "Pipra_aureola_MPEG1199_0"
  ),
  c("Pipra_aureola_MPEG16053_1", "Pipra_aureola_MPEG16053_0"
  ),
  c("Pipra_fasciicauda_FMNH433698_1", "Pipra_fasciicauda_MPEG10317_1", "Pipra_fasciicauda_MPEG17411_1", 
    "Pipra_fasciicauda_FMNH433698_0", "Pipra_fasciicauda_MPEG10317_0", "Pipra_fasciicauda_MPEG17411_0"
  ),
  c("Pipra_filicauda_AMNH4246_1", "Pipra_filicauda_INPA2095_1", "Pipra_filicauda_LSUMNS42835_1", 
    "Pipra_filicauda_MPEG3278_1", "Pipra_filicauda_AMNH4246_0", "Pipra_filicauda_INPA2095_0", 
    "Pipra_filicauda_LSUMNS42835_0", "Pipra_filicauda_MPEG3278_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Schiffornis_major

spec <- "Schiffornis_major"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Schiffornis_major_FMNH457274_0", "Schiffornis_major_LSUMNS103592_0", "Schiffornis_major_LSUMNS42991_0", 
    "Schiffornis_major_LSUMNS4462_0", "Schiffornis_major_LSUMNS7503_0", "Schiffornis_major_MPEG15841_0", 
    "Schiffornis_major_MPEG15853_0", "Schiffornis_major_MPEG16040_0", "Schiffornis_major_MPEG6239_0", 
    "Schiffornis_major_FMNH457274_1", "Schiffornis_major_LSUMNS103592_1", "Schiffornis_major_LSUMNS42991_1", 
    "Schiffornis_major_LSUMNS4462_1", "Schiffornis_major_LSUMNS7503_1", "Schiffornis_major_MPEG15841_1", 
    "Schiffornis_major_MPEG15853_1", "Schiffornis_major_MPEG16040_1", "Schiffornis_major_MPEG6239_1"
  ),
  c("Schiffornis_major_KU1426_0", "Schiffornis_major_KU1426_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Synallaxis_gujanensis

spec <- "Synallaxis_gujanensis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Synallaxis_gujanensis_AMNH11985_1", "Synallaxis_gujanensis_INPA2220_1", "Synallaxis_gujanensis_INPA939_1", 
    "Synallaxis_gujanensis_LSUMNS7276_1", "Synallaxis_gujanensis_MPEG1885_1", "Synallaxis_gujanensis_MPEG2288_1", 
    "Synallaxis_gujanensis_MPEG7070_1", "Synallaxis_gujanensis_USNMB06925_1", "Synallaxis_gujanensis_AMNH11985_0", 
    "Synallaxis_gujanensis_INPA2220_0", "Synallaxis_gujanensis_INPA939_0", "Synallaxis_gujanensis_LSUMNS7276_0", 
    "Synallaxis_gujanensis_MPEG1885_0", "Synallaxis_gujanensis_MPEG2288_0", "Synallaxis_gujanensis_MPEG7070_0", 
    "Synallaxis_gujanensis_USNMB06925_0"
  ),
  c("Synallaxis_gujanensis_LSUMNS46034_1", "Synallaxis_gujanensis_FMNH321495_1", "Synallaxis_gujanensis_AMNH6193_1", 
    "Synallaxis_gujanensis_LSUMNS46034_0", "Synallaxis_gujanensis_FMNH321495_0", "Synallaxis_gujanensis_AMNH6193_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Synallaxis_rutilans

spec <- "Synallaxis_rutilans"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Synallaxis_rutilans_FMNH391101_0", "Synallaxis_rutilans_LSUMNS31362_0", "Synallaxis_rutilans_MPEG12491_0", 
    "Synallaxis_rutilans_MPEG1345_0", "Synallaxis_rutilans_MPEG16900_0", "Synallaxis_rutilans_MPEG440_0", 
    "Synallaxis_rutilans_USNMB06831_0", "Synallaxis_rutilans_FMNH391101_1", "Synallaxis_rutilans_LSUMNS31362_1", 
    "Synallaxis_rutilans_MPEG12491_1", "Synallaxis_rutilans_MPEG1345_1", "Synallaxis_rutilans_MPEG16900_1", 
    "Synallaxis_rutilans_MPEG440_1", "Synallaxis_rutilans_USNMB06831_1"
  ),
  c("Synallaxis_rutilans_LSUMNS42924_0", "Synallaxis_rutilans_MPEG5105_0", "Synallaxis_rutilans_INPA1617_0", 
    "Synallaxis_rutilans_AMNH11928_0", "Synallaxis_rutilans_LSUMNS42924_1", "Synallaxis_rutilans_MPEG5105_1", 
    "Synallaxis_rutilans_INPA1617_1", "Synallaxis_rutilans_AMNH11928_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Xiphorhynchus_elegans

spec <- "Xiphorhynchus_elegans"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Xiphorhynchus_elegans_INPA8144_1", "Xiphorhynchus_elegans_LSUMNS36629_1", "Xiphorhynchus_elegans_LSUMNS18545_1", 
    "Xiphorhynchus_elegans_INPA8144_0", "Xiphorhynchus_elegans_LSUMNS36629_0", "Xiphorhynchus_elegans_LSUMNS18545_0"
  ),
  c("Xiphorhynchus_spixii_FMNH456868_1", "Xiphorhynchus_spixii_MPEG10312_1", "Xiphorhynchus_spixii_MPEG12556_1", 
    "Xiphorhynchus_spixii_MPEG3585_1", "Xiphorhynchus_spixii_FMNH456868_0", "Xiphorhynchus_spixii_MPEG10312_0", 
    "Xiphorhynchus_spixii_MPEG12556_0", "Xiphorhynchus_spixii_MPEG3585_0"
  ),
  c("Xiphorhynchus_elegans_LSUMNS35681_1", "Xiphorhynchus_elegans_LSUMNS42949_1", "Xiphorhynchus_elegans_MPEG7671_1", 
    "Xiphorhynchus_elegans_LSUMNS35681_0", "Xiphorhynchus_elegans_LSUMNS42949_0", "Xiphorhynchus_elegans_MPEG7671_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Xiphorhynchus_obsoletus

spec <- "Xiphorhynchus_obsoletus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("Xiphorhynchus_obsoletus_MPEG11024_0", "Xiphorhynchus_obsoletus_MPEG17364_0", "Xiphorhynchus_obsoletus_MPEG2456_0", 
    "Xiphorhynchus_obsoletus_LSUMNS35670_0", "Xiphorhynchus_obsoletus_MPEG11024_1", "Xiphorhynchus_obsoletus_MPEG17364_1", 
    "Xiphorhynchus_obsoletus_MPEG2456_1", "Xiphorhynchus_obsoletus_LSUMNS35670_1"
  ),
  c("Xiphorhynchus_obsoletus_AMNH12343_0", "Xiphorhynchus_obsoletus_FMNH391089_0", "Xiphorhynchus_obsoletus_FMNH456853_0", 
    "Xiphorhynchus_obsoletus_LSUMNS12752_0", "Xiphorhynchus_obsoletus_LSUMNS35642_0", "Xiphorhynchus_obsoletus_LSUMNS35734_0", 
    "Xiphorhynchus_obsoletus_LSUMNS42526_0", "Xiphorhynchus_obsoletus_AMNH12343_1", "Xiphorhynchus_obsoletus_FMNH391089_1", 
    "Xiphorhynchus_obsoletus_FMNH456853_1", "Xiphorhynchus_obsoletus_LSUMNS12752_1", "Xiphorhynchus_obsoletus_LSUMNS35642_1", 
    "Xiphorhynchus_obsoletus_LSUMNS35734_1", "Xiphorhynchus_obsoletus_LSUMNS42526_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Conirostrum_bicolor

spec <- "Conirostrum_bicolor"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[1]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("conirostrum_bicolor_66_UFPE_T1161_1", "conirostrum_bicolor_MPEG_T611_1", 
    "conirostrum_bicolor_66_UFPE_T1161_0", "conirostrum_bicolor_MPEG_T611_0"
  ),
  c("conirostrum_bicolor_7_Goeldi22702_1", "conirostrum_bicolor_9_LSU7282_1", "conirostrum_bicolor_10_LSU43023_1", 
    "conirostrum_bicolor_11_LSU93328_1", "conirostrum_bicolor_12_INPA1341_1", "conirostrum_bicolor_7_Goeldi22702_0", 
    "conirostrum_bicolor_9_LSU7282_0", "conirostrum_bicolor_10_LSU43023_0", "conirostrum_bicolor_11_LSU93328_0", 
    "conirostrum_bicolor_12_INPA1341_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Conirostrum_margaritae

spec <- "Conirostrum_margaritae"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("conirostrum_margaritae_13_LSU25428_0", "conirostrum_margaritae_14_LSU7293_0", 
    "conirostrum_margaritae_13_LSU25428_1", "conirostrum_margaritae_14_LSU7293_1"
  ),
  c("conirostrum_margaritae_15_LSU93298_0", "conirostrum_margaritae_16_Goeldi15984_0", 
    "conirostrum_margaritae_15_LSU93298_1", "conirostrum_margaritae_16_Goeldi15984_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Cranioleuca_vulpecula

spec <- "Cranioleuca_vulpecula"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[21]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("cranioleuca_vulpecula_22_LSU79790_1", "cranioleuca_vulpecula_25_INPA_A123_1", 
    "cranioleuca_vulpecula_22_LSU79790_0", "cranioleuca_vulpecula_25_INPA_A123_0"
  ),
  c("cranioleuca_vulpecula_19_LSU25424_1", "cranioleuca_vulpecula_20_LSU3637_1", "cranioleuca_vulpecula_21_LSU74799_1", 
    "cranioleuca_vulpecula_23_LSU3181_1", "cranioleuca_vulpecula_24_LSU93331_1", "cranioleuca_vulpecula_19_LSU25424_0", 
    "cranioleuca_vulpecula_20_LSU3637_0", "cranioleuca_vulpecula_21_LSU74799_0", "cranioleuca_vulpecula_23_LSU3181_0", 
    "cranioleuca_vulpecula_24_LSU93331_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Dendroplex_kienerii

spec <- "Dendroplex_kienerii"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("dendroplex_kienerii_26_LSU80430_0", "dendroplex_kienerii_28_LSU35627_0", "dendroplex_kienerii_30_LSU20237_0", 
    "dendroplex_kienerii_31_LSU25413_0", "dendroplex_kienerii_26_LSU80430_1", "dendroplex_kienerii_28_LSU35627_1", 
    "dendroplex_kienerii_30_LSU20237_1", "dendroplex_kienerii_31_LSU25413_1"
  ),
  c("dendroplex_kienerii_27_LSU29023_0", "dendroplex_kienerii_29_LSU35692_0", "dendroplex_kienerii_32_LSU93459_0", 
    "dendroplex_kienerii_27_LSU29023_1", "dendroplex_kienerii_29_LSU35692_1", "dendroplex_kienerii_32_LSU93459_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Elaenia_pelzelni

spec <- "Elaenia_pelzelni"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("elaenia_pelzelni_36_Goeldi22693_0", "elaenia_pelzelni_36_Goeldi22693_1"
  ),
  c("elaenia_pelzelni_35_LSU7249_0", "elaenia_pelzelni_35_LSU7249_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Furnarius_minor

spec <- "Furnarius_minor"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[21]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("furnarius_minor_160_INPA1339_0", "furnarius_minor_38_LSU3641_0", "furnarius_minor_40_LSU74803_0", 
    "furnarius_minor_42_LSU7265_0", "furnarius_minor_160_INPA1339_1", "furnarius_minor_38_LSU3641_1", 
    "furnarius_minor_40_LSU74803_1", "furnarius_minor_42_LSU7265_1"
  ),
  c("furnarius_minor_39_LSU45888_0", "furnarius_minor_39_LSU45888_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Knipolegus_orenocensis

spec <- "Knipolegus_orenocensis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[21]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("knipolegus_orenocensis_46_COP882_0", "knipolegus_orenocensis_47_COP629_0", 
    "knipolegus_orenocensis_46_COP882_1", "knipolegus_orenocensis_47_COP629_1"
  ),
  c("knipolegus_orenocensis_48_MPEG_T15658_0", "knipolegus_orenocensis_49_USNM303_0", 
    "knipolegus_orenocensis_48_MPEG_T15658_1", "knipolegus_orenocensis_49_USNM303_1"
  ),
  c("knipolegus_orenocensis_50_LSU43080_0", "knipolegus_orenocensis_51_MPEG_T16338_0", "knipolegus_orenocensis_52_INPA_A11140_0", 
    "knipolegus_orenocensis_53_INPA_A18068_0", "knipolegus_orenocensis_54_LSU75980_0", "knipolegus_orenocensis_55_LSU3178_0", 
    "knipolegus_orenocensis_56_LSU3647_0", "knipolegus_orenocensis_57_LSU93474_0", "knipolegus_orenocensis_50_LSU43080_1", 
    "knipolegus_orenocensis_51_MPEG_T16338_1", "knipolegus_orenocensis_52_INPA_A11140_1", "knipolegus_orenocensis_53_INPA_A18068_1", 
    "knipolegus_orenocensis_54_LSU75980_1", "knipolegus_orenocensis_55_LSU3178_1", "knipolegus_orenocensis_56_LSU3647_1", 
    "knipolegus_orenocensis_57_LSU93474_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Talaphorus_chlorocercus

spec <- "Talaphorus_chlorocercus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[100]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("leucippus_chlorocercus_62_LSU74801_1", "leucippus_chlorocercus_63_LSU43036_1", "leucippus_chlorocercus_64_LSU93457_1", 
    "leucippus_chlorocercus_62_LSU74801_0", "leucippus_chlorocercus_63_LSU43036_0", "leucippus_chlorocercus_64_LSU93457_0"
  ),
  c("leucippus_chlorocercus_58_LSU7264_1", "leucippus_chlorocercus_58_LSU7264_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Mazaria_propinqua

spec <- "Mazaria_propinqua"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[100]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("mazaria_propinqua_132_LSU7289_0", "mazaria_propinqua_133_LSU75947_0", "mazaria_propinqua_134_LSU79789_0", 
    "mazaria_propinqua_135_LSU43083_0", "mazaria_propinqua_136_LSU93326_0", "mazaria_propinqua_137_MPEG_T16312_0", 
    "mazaria_propinqua_138_UFPE_T1027_0", "mazaria_propinqua_132_LSU7289_1", "mazaria_propinqua_133_LSU75947_1", 
    "mazaria_propinqua_134_LSU79789_1", "mazaria_propinqua_135_LSU43083_1", "mazaria_propinqua_136_LSU93326_1", 
    "mazaria_propinqua_137_MPEG_T16312_1", "mazaria_propinqua_138_UFPE_T1027_1"
  ),
  c("mazaria_propinqua_241_UFPE5254_0", "mazaria_propinqua_241_UFPE5254_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmoborus_lugubris

spec <- "Myrmoborus_lugubris"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("myrmoborus_lugubris_76_LSU93417_1", "myrmoborus_lugubris_242_MPEG_T22612_1", 
    "myrmoborus_lugubris_76_LSU93417_0", "myrmoborus_lugubris_242_MPEG_T22612_0"
  ),
  c("myrmoborus_lugubris_70_INPA_A2193_1", "myrmoborus_lugubris_74_LSU25513_1", 
    "myrmoborus_lugubris_70_INPA_A2193_0", "myrmoborus_lugubris_74_LSU25513_0"
  ),
  c("myrmoborus_lugubris_169_MPEG73846_1", "myrmoborus_lugubris_71_MPEG_T22926_1", "myrmoborus_lugubris_72_MPEG_T22751_1", 
    "myrmoborus_lugubris_73_LSU80385_1", "myrmoborus_lugubris_169_MPEG73846_0", "myrmoborus_lugubris_71_MPEG_T22926_0", 
    "myrmoborus_lugubris_72_MPEG_T22751_0", "myrmoborus_lugubris_73_LSU80385_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmochanes_hemileucus

spec <- "Myrmochanes_hemileucus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("myrmochanes_hemileucus_77_LSU74774_0", "myrmochanes_hemileucus_78_LSU43093_0", "myrmochanes_hemileucus_79_LSU3649_0", 
    "myrmochanes_hemileucus_77_LSU74774_1", "myrmochanes_hemileucus_78_LSU43093_1", "myrmochanes_hemileucus_79_LSU3649_1"
  ),
  c("myrmochanes_hemileucus_80_MPEG_T16307_1", "myrmochanes_hemileucus_81_MPEG_T22712_1", "myrmochanes_hemileucus_83_INPA_A1340_1", 
    "myrmochanes_hemileucus_84_INPA_A11135_1", "myrmochanes_hemileucus_80_MPEG_T16307_0", "myrmochanes_hemileucus_81_MPEG_T22712_0", 
    "myrmochanes_hemileucus_83_INPA_A1340_0", "myrmochanes_hemileucus_84_INPA_A11135_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmotherula_assimilis

spec <- "Myrmotherula_assimilis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("myrmotherula_assimilis_87_LSU93460_0", "myrmotherula_assimilis_88_MPEG_T14289_0", 
    "myrmotherula_assimilis_87_LSU93460_1", "myrmotherula_assimilis_88_MPEG_T14289_1"
  ),
  c("myrmotherula_assimilis_243_MPEG_T22595_0", "myrmotherula_assimilis_85_LSU23703_0", "myrmotherula_assimilis_89_MPEG_T23587_0", 
    "myrmotherula_assimilis_93_INPA_A8418_0", "myrmotherula_assimilis_94_MPEG_T73294_0", "myrmotherula_assimilis_243_MPEG_T22595_1", 
    "myrmotherula_assimilis_85_LSU23703_1", "myrmotherula_assimilis_89_MPEG_T23587_1", "myrmotherula_assimilis_93_INPA_A8418_1", 
    "myrmotherula_assimilis_94_MPEG_T73294_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Myrmotherula_klagesi

spec <- "Myrmotherula_klagesi"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("myrmotherula_klagesi_96_LSU25562_1", "myrmotherula_klagesi_97_LSU81357_1", "myrmotherula_klagesi_101_INPA_A048_1", 
    "myrmotherula_klagesi_96_LSU25562_0", "myrmotherula_klagesi_97_LSU81357_0", "myrmotherula_klagesi_101_INPA_A048_0"
  ),
  c("myrmotherula_klagesi_100_INPA_A15925_1", "myrmotherula_klagesi_95_LSU20250_1", "myrmotherula_klagesi_98_LSU25511_1", 
    "myrmotherula_klagesi_99_INPA_A8276_1", "myrmotherula_klagesi_100_INPA_A15925_0", "myrmotherula_klagesi_95_LSU20250_0", 
    "myrmotherula_klagesi_98_LSU25511_0", "myrmotherula_klagesi_99_INPA_A8276_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Ochthornis_littoralis

spec <- "Ochthornis_littoralis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("ochthornis_littoralis_103_LSU968_0", "ochthornis_littoralis_107_LSU93193_0", "ochthornis_littoralis_110_INPA_A4762_0", 
    "ochthornis_littoralis_111_INPA_A8301_0", "ochthornis_littoralis_113_MPEG_T963_0", "ochthornis_littoralis_114_MPEG_T22557_0", 
    "ochthornis_littoralis_115_MPEG_T1998_0", "ochthornis_littoralis_116_MPEG_T23235_0", "ochthornis_littoralis_170_LSU40723_0", 
    "ochthornis_littoralis_103_LSU968_1", "ochthornis_littoralis_107_LSU93193_1", "ochthornis_littoralis_110_INPA_A4762_1", 
    "ochthornis_littoralis_111_INPA_A8301_1", "ochthornis_littoralis_113_MPEG_T963_1", "ochthornis_littoralis_114_MPEG_T22557_1", 
    "ochthornis_littoralis_115_MPEG_T1998_1", "ochthornis_littoralis_116_MPEG_T23235_1", "ochthornis_littoralis_170_LSU40723_1"
  ),
  c("ochthornis_littoralis_102_LSU81242_0", "ochthornis_littoralis_105_LSU35432_0", "ochthornis_littoralis_112_MPEG_T19792_0", 
    "ochthornis_littoralis_102_LSU81242_1", "ochthornis_littoralis_105_LSU35432_1", "ochthornis_littoralis_112_MPEG_T19792_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Serpophaga_hypoleuca

spec <- "Serpophaga_hypoleuca"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("serpophaga_hypoleuca_118_LSU43033_0", "serpophaga_hypoleuca_119_LSU74789_0", "serpophaga_hypoleuca_120_LSU93475_0", 
    "serpophaga_hypoleuca_121_LSU79793_0", "serpophaga_hypoleuca_123_INPA_A15959_0", "serpophaga_hypoleuca_125_COP632_0", 
    "serpophaga_hypoleuca_167_INPA_A1342_0", "serpophaga_hypoleuca_118_LSU43033_1", "serpophaga_hypoleuca_119_LSU74789_1", 
    "serpophaga_hypoleuca_120_LSU93475_1", "serpophaga_hypoleuca_121_LSU79793_1", "serpophaga_hypoleuca_123_INPA_A15959_1", 
    "serpophaga_hypoleuca_125_COP632_1", "serpophaga_hypoleuca_167_INPA_A1342_1"
  ),
  c("serpophaga_hypoleuca_T122_1", "serpophaga_hypoleuca_122_UFPE_T1166_1", 
    "serpophaga_hypoleuca_T122_0", "serpophaga_hypoleuca_122_UFPE_T1166_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Stigmatura_napensis

spec <- "Stigmatura_napensis"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("stigmatura_napensis_199_UFPE_5249_0", "stigmatura_napensis_199_UFPE_5249_1"
  ),
  c("stigmatura_napensis_126_LSU3175_1", "stigmatura_napensis_duplicateLoc_172_LSU7240_1", "stigmatura_napensis_159_LSU89217_1", 
    "stigmatura_napensis_127_LSU43079_1", "stigmatura_napensis_130_INPA_A15960_1", "stigmatura_napensis_128_MPEG_T16303_1", 
    "stigmatura_napensis_126_LSU3175_0", "stigmatura_napensis_duplicateLoc_172_LSU7240_0", "stigmatura_napensis_159_LSU89217_0", 
    "stigmatura_napensis_127_LSU43079_0", "stigmatura_napensis_130_INPA_A15960_0", "stigmatura_napensis_128_MPEG_T16303_0"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------
# Thamnophilus

spec <- "Thamnophilus"
fname <- paste0("2_phasing/2_incomplete-taxon-set/", spec, "-phased-mafft-nexus-edge-trimmed-clean-nooutgroup")
GENOME.class <- readData(fname, format = 'nexus')
sp <- get.individuals(GENOME.class)[[20]]

# population assignments based on DAPC clusters
GENOME.class <- set.populations(GENOME.class, list(
  c("thamnophilus_cryptoleucus_139_LSU25431_0", "thamnophilus_cryptoleucus_141_LSU74103_0", "thamnophilus_cryptoleucus_142_LSU7285_0", 
    "thamnophilus_cryptoleucus_144_LSU93318_0", "thamnophilus_cryptoleucus_145_MPEG_T22677_0", "thamnophilus_cryptoleucus_147_MPEG_T22714_0", 
    "thamnophilus_cryptoleucus_244_MPEG_T22585_0", "thamnophilus_cryptoleucus_139_LSU25431_1", "thamnophilus_cryptoleucus_141_LSU74103_1",
    "thamnophilus_cryptoleucus_142_LSU7285_1", "thamnophilus_cryptoleucus_144_LSU93318_1", "thamnophilus_cryptoleucus_145_MPEG_T22677_1", 
    "thamnophilus_cryptoleucus_147_MPEG_T22714_1", "thamnophilus_cryptoleucus_244_MPEG_T22585_1"
  ),
  c("thamnophilus_cryptoleucus_146_MPEG_T23595_0", "thamnophilus_nigrocinereus_151_INPA_A10843_0", 
    "thamnophilus_nigrocinereus_175_MZUSP92841_0", "thamnophilus_nigrocinereus_176_MZUSP93302_0", 
    "thamnophilus_nigrocinereus_158_INPA_A295_0", "thamnophilus_nigrocinereus_177_MZUSP_OI001_0", 
    "thamnophilus_cryptoleucus_146_MPEG_T23595_1", "thamnophilus_nigrocinereus_151_INPA_A10843_1", 
    "thamnophilus_nigrocinereus_175_MZUSP92841_1", "thamnophilus_nigrocinereus_176_MZUSP93302_1", 
    "thamnophilus_nigrocinereus_158_INPA_A295_1", "thamnophilus_nigrocinereus_177_MZUSP_OI001_1"
  ),
  c("thamnophilus_nigrocinereus_154_INPA_A2165_0", "thamnophilus_nigrocinereus_155_INPA_A3093_0", 
    "thamnophilus_nigrocinereus_165_MPEG_T77261_0", "thamnophilus_nigrocinereus_173_ICN39301_0", 
    "thamnophilus_nigrocinereus_174_LSU20234_0", "thamnophilus_nigrocinereus_154_INPA_A2165_1", 
    "thamnophilus_nigrocinereus_155_INPA_A3093_1", "thamnophilus_nigrocinereus_165_MPEG_T77261_1", 
    "thamnophilus_nigrocinereus_173_ICN39301_1", "thamnophilus_nigrocinereus_174_LSU20234_1"
  )
))

# Fst / Gst
GENOME.class <- F_ST.stats(GENOME.class)
res <- as.data.frame(get.F_ST(GENOME.class))
Fst <- colMeans(res, na.rm = TRUE)

# Dxy
GENOME.class <- diversity.stats.between(GENOME.class)
dxy <- mean(GENOME.class@nuc.diversity.between / GENOME.class@n.sites)

res.sp <- bind_cols(spec, t(as.data.frame(Fst)), dxy)
colnames(res.sp)[1] <- "species"
colnames(res.sp)[8] <- "Dxy"
res.df <- bind_rows(res.df, res.sp)

#------------------------------------------------------------------------

fname.out <- paste0("3_results/PopGenome_Fst.results.v1.DAPCassignments.csv")
write.csv(res.df, file = fname.out, row.names = FALSE)

