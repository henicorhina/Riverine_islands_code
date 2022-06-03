library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# stigmatura_napensis
# Stigmatura napensis
# stigmatura_napensis
# 1208 with your number of snps
# xlim=c(-12, 8), ylim=c(-8, 12), 

stigmatura_napensis <- read.structure("stigmatura_napensis_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1208, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/stigmatura_napensis")
stigmatura_napensis.X <- tab(stigmatura_napensis, freq=TRUE, NA.method = "mean")
stigmatura_napensis.X <- scaleGen(stigmatura_napensis, NA.method = "mean")
sum(is.na(stigmatura_napensis$tab)) #amount of missing data
stigmatura_napensis.X[1:2,1:2]
stigmatura_napensis.pca <- dudi.pca(stigmatura_napensis.X, scale=FALSE, scannf = FALSE, nf = 4)
stigmatura_napensis.pca$eig[1]/sum(stigmatura_napensis.pca$eig) * 100
stigmatura_napensis.pca$eig[2]/sum(stigmatura_napensis.pca$eig) * 100

stigmatura_napensis.pca$li

stigmatura_napensis.grp <- find.clusters(stigmatura_napensis, n.pca = 100)

xval <- xvalDapc(stigmatura_napensis.X, stigmatura_napensis.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(stigmatura_napensis.X, stigmatura_napensis.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Plot.png", width = 6, height = 6, units = 'in', 
    res = 800)
bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
       type="tifflzw", res=800)
# run the desired format above then the plot, then this dev.off command
dev.off()
par(mfrow = c(1,1))


# pca with colored dots
tiff("stigmatura_napensis_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(stigmatura_napensis.pca$li, stigmatura_napensis.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 8), ylim=c(-8, 12), 
          xlab="PC 1 (20.5%)", ylab="PC 2 (18.2%)",
          main = "PCA of Stigmatura napensis 1208 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("stigmatura_napensis_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis.pca,
      	     xlim=c(-12, 8), ylim=c(-8, 12), 
             repel = TRUE,
             title = "PCA of Stigmatura napensis 1208 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("stigmatura_napensis_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis.pca,
             xlim=c(-12, 8), ylim=c(-8, 12),
             geom = "point",
             title = "PCA of Stigmatura napensis 1208 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("stigmatura_napensis_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis.pca,
             geom = "point",
             xlim=c(-12, 8), ylim=c(-8, 12), 
             title = "PCA of Stigmatura napensis 1208 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(stigmatura_napensis)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr)
contrib <- loadingplot(dapc1$var.contr, threshold = 0.002, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign







# STINAP FULL

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# stigmatura_napensis_full
# Stigmatura napensis full
# stigmatura_napensis_full
# 1208 with your number of snps
# xlim=c(-10, 10), ylim=c(-10, 10), 

stigmatura_napensis_full <- read.structure("stigmatura_napensis_full_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                                      n.ind=11, n.loc=1208, onerowperind=FALSE, 
                                      col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/stigmatura_napensis_full")
stigmatura_napensis_full.X <- tab(stigmatura_napensis_full, freq=TRUE, NA.method = "mean")
stigmatura_napensis_full.X <- scaleGen(stigmatura_napensis_full, NA.method = "mean")
sum(is.na(stigmatura_napensis_full$tab)) #amount of missing data
stigmatura_napensis_full.X[1:2,1:2]
stigmatura_napensis_full.pca <- dudi.pca(stigmatura_napensis_full.X, scale=FALSE, scannf = FALSE, nf = 4)
stigmatura_napensis_full.pca$eig[1]/sum(stigmatura_napensis_full.pca$eig) * 100
stigmatura_napensis_full.pca$eig[2]/sum(stigmatura_napensis_full.pca$eig) * 100

stigmatura_napensis_full.pca$li

stigmatura_napensis_full.grp <- find.clusters(stigmatura_napensis_full, n.pca = 100)

xval <- xvalDapc(stigmatura_napensis_full.X, stigmatura_napensis_full.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(stigmatura_napensis_full.X, stigmatura_napensis_full.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Plot.png", width = 6, height = 6, units = 'in', 
    res = 800)
bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
       type="tifflzw", res=800)
# run the desired format above then the plot, then this dev.off command
dev.off()
par(mfrow = c(1,1))


# pca with colored dots
tiff("stigmatura_napensis_full_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(stigmatura_napensis_full.pca$li, stigmatura_napensis_full.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (14.3%)", ylab="PC 2 (12.3%)",
          main = "PCA of Stigmatura napensis full 1208 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("stigmatura_napensis_full_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis_full.pca,
             xlim=c(-10, 10), ylim=c(-10, 10), 
             repel = TRUE,
             title = "PCA of Stigmatura napensis full 1208 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("stigmatura_napensis_full_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis_full.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             geom = "point",
             title = "PCA of Stigmatura napensis full 1208 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("stigmatura_napensis_full_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(stigmatura_napensis_full.pca,
             geom = "point",
             xlim=c(-10, 10), ylim=c(-10, 10), 
             title = "PCA of Stigmatura napensis full 1208 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(stigmatura_napensis_full)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr)
contrib <- loadingplot(dapc1$var.contr, threshold = 0.003, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign




