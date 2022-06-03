library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# furnarius_minor
# Furnarius minor
# furnarius_minor
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

furnarius_minor <- read.structure("furnarius_minor_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=5, n.loc=1080, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/furnarius_minor")
furnarius_minor.X <- tab(furnarius_minor, freq=TRUE, NA.method = "mean")
furnarius_minor.X <- scaleGen(furnarius_minor, NA.method = "mean")
sum(is.na(furnarius_minor$tab)) #amount of missing data
furnarius_minor.X[1:2,1:2]
furnarius_minor.pca <- dudi.pca(furnarius_minor.X, scale=FALSE, scannf = FALSE, nf = 4)
furnarius_minor.pca$eig[1]/sum(furnarius_minor.pca$eig) * 100
furnarius_minor.pca$eig[2]/sum(furnarius_minor.pca$eig) * 100

furnarius_minor.pca$li

furnarius_minor.grp <- find.clusters(furnarius_minor, n.pca = 100)

xval <- xvalDapc(furnarius_minor.X, furnarius_minor.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc1 <- dapc(furnarius_minor.X, furnarius_minor.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp1$tab[1:2,1:2]
temp1$pop.score


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
tiff("furnarius_minor_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(furnarius_minor.pca$li, furnarius_minor.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (34.8%)", ylab="PC 2 (23.8%)",
          main = "PCA of Furnarius minor 1080 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("furnarius_minor_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Furnarius minor 1080 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("furnarius_minor_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Furnarius minor 1080 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("furnarius_minor_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Furnarius minor 1080 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(furnarius_minor)
xval$DAPC$assign
xval$DAPC$n.pca
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







# FURMIN FULL


furnarius_minor_full <- read.structure("furnarius_minor_full_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                                  n.ind=7, n.loc=1080, onerowperind=FALSE, 
                                  col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/furnarius_minor_full")
furnarius_minor_full.X <- tab(furnarius_minor_full, freq=TRUE, NA.method = "mean")
furnarius_minor_full.X <- scaleGen(furnarius_minor_full, NA.method = "mean")
sum(is.na(furnarius_minor_full$tab)) #amount of missing data
furnarius_minor_full.X[1:2,1:2]
furnarius_minor_full.pca <- dudi.pca(furnarius_minor_full.X, scale=FALSE, scannf = FALSE, nf = 4)
furnarius_minor_full.pca$eig[1]/sum(furnarius_minor_full.pca$eig) * 100
furnarius_minor_full.pca$eig[2]/sum(furnarius_minor_full.pca$eig) * 100

furnarius_minor_full.pca$li

furnarius_minor_full.grp <- find.clusters(furnarius_minor_full, n.pca = 100)

xval <- xvalDapc(furnarius_minor_full.X, furnarius_minor_full.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc1 <- dapc(furnarius_minor_full.X, furnarius_minor_full.grp$grp)

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
tiff("furnarius_minor_full_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(furnarius_minor_full.pca$li, furnarius_minor_full.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (25.6%)", ylab="PC 2 (19.4%)",
          main = "PCA of Furnarius minor full 1080 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("furnarius_minor_full_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor_full.pca,
             xlim=c(-12, 12), ylim=c(-12, 12), 
             repel = TRUE,
             title = "PCA of Furnarius minor full 1080 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("furnarius_minor_full_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor_full.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Furnarius minor full 1080 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("furnarius_minor_full_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(furnarius_minor_full.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12), 
             title = "PCA of Furnarius minor full 1080 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(furnarius_minor_full)
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

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
