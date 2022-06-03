library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Cantorchilus_leucotis
# Cantorchilus leucotis
# Cantorchilus_leucotis
# 1985 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 20), 

Cantorchilus_leucotis <- read.structure("Cantorchilus_leucotis_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=10, n.loc=1985, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Cantorchilus_leucotis_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Cantorchilus_leucotis_AmazonOnly")
Cantorchilus_leucotis.X <- tab(Cantorchilus_leucotis, freq=TRUE, NA.method = "mean")
Cantorchilus_leucotis.X <- scaleGen(Cantorchilus_leucotis, NA.method = "mean")
sum(is.na(Cantorchilus_leucotis$tab)) #amount of missing data
Cantorchilus_leucotis.X[1:2,1:2]
Cantorchilus_leucotis.pca <- dudi.pca(Cantorchilus_leucotis.X, scale=FALSE, scannf = FALSE, nf = 4)
Cantorchilus_leucotis.pca$eig[1]/sum(Cantorchilus_leucotis.pca$eig) * 100
Cantorchilus_leucotis.pca$eig[2]/sum(Cantorchilus_leucotis.pca$eig) * 100

Cantorchilus_leucotis.pca$li

Cantorchilus_leucotis.grp <- find.clusters(Cantorchilus_leucotis, n.pca = 100)

xval <- xvalDapc(Cantorchilus_leucotis.X, Cantorchilus_leucotis.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Cantorchilus_leucotis.X, Cantorchilus_leucotis.grp$grp)

temp <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
# tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
#      res = 800, compression = 'none')
# png("Plot.png", width = 6, height = 6, units = 'in', 
#     res = 800)
# bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
#        type="tifflzw", res=800)
# # run the desired format above then the plot, then this dev.off command
# dev.off()
# par(mfrow = c(1,1))


# pca with colored dots
tiff("Cantorchilus_leucotis_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Cantorchilus_leucotis_colors.png", width = 6, height = 6, units = 'in',
    res = 800)
colorplot(Cantorchilus_leucotis.pca$li, Cantorchilus_leucotis.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-45, 45), ylim=c(-40, 50), 
          xlab="PC 1 (21.5%)", ylab="PC 2 (13.7%)",
          main = "PCA of Cantorchilus leucotis 1985 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Cantorchilus_leucotis_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Cantorchilus_leucotis_labeled.png", width = 12, height = 12, units = 'in',
    res = 800)
fviz_pca_ind(Cantorchilus_leucotis.pca,
             xlim=c(-45, 45), ylim=c(-40, 50), 
             repel = TRUE,
             title = "PCA of Cantorchilus leucotis 1985 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Cantorchilus_leucotis_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Cantorchilus_leucotis_groups.png", width = 6, height = 6, units = 'in',
    res = 800)
fviz_pca_ind(Cantorchilus_leucotis.pca,
             xlim=c(-45, 45), ylim=c(-40, 50), 
             geom = "point",
             title = "PCA of Cantorchilus leucotis 1985 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Cantorchilus_leucotis_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Cantorchilus_leucotis_black.png", width = 6, height = 6, units = 'in',
    res = 800)
fviz_pca_ind(Cantorchilus_leucotis.pca,
             geom = "point",
             xlim=c(-45, 45), ylim=c(-40, 50), 
             title = "PCA of Cantorchilus leucotis 1985 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Cantorchilus_leucotis)
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
