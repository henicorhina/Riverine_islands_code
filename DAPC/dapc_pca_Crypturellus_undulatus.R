library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Crypturellus_undulatus
# Crypturellus undulatus
# Crypturellus_undulatus
# 1818 with your number of snps
# xlim=c(-20, 20), ylim=c(-20, 20), 

Crypturellus_undulatus <- read.structure("Crypturellus_undulatus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=9, n.loc=1818, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Crypturellus_undulatus")
Crypturellus_undulatus.X <- tab(Crypturellus_undulatus, freq=TRUE, NA.method = "mean")
Crypturellus_undulatus.X <- scaleGen(Crypturellus_undulatus, NA.method = "mean")
sum(is.na(Crypturellus_undulatus$tab)) #amount of missing data
Crypturellus_undulatus.X[1:2,1:2]
Crypturellus_undulatus.pca <- dudi.pca(Crypturellus_undulatus.X, scale=FALSE, scannf = FALSE, nf = 4)
Crypturellus_undulatus.pca$eig[1]/sum(Crypturellus_undulatus.pca$eig) * 100
Crypturellus_undulatus.pca$eig[2]/sum(Crypturellus_undulatus.pca$eig) * 100

Crypturellus_undulatus.pca$li

Crypturellus_undulatus.grp <- find.clusters(Crypturellus_undulatus, n.pca = 100)

xval <- xvalDapc(Crypturellus_undulatus.X, Crypturellus_undulatus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Crypturellus_undulatus.X, Crypturellus_undulatus.grp$grp)

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
tiff("Crypturellus_undulatus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Crypturellus_undulatus.pca$li, Crypturellus_undulatus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-20, 20), ylim=c(-20, 20), 
          xlab="PC 1 (38.7%)", ylab="PC 2 (18.3%)",
          main = "PCA of Crypturellus undulatus 1818 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Crypturellus_undulatus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_undulatus.pca,
      	     xlim=c(-20, 20), ylim=c(-20, 20), 
             repel = TRUE,
             title = "PCA of Crypturellus undulatus 1818 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Crypturellus_undulatus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_undulatus.pca,
             xlim=c(-20, 20), ylim=c(-20, 20),
             geom = "point",
             title = "PCA of Crypturellus undulatus 1818 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Crypturellus_undulatus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_undulatus.pca,
             geom = "point",
             xlim=c(-20, 20), ylim=c(-20, 20), 
             title = "PCA of Crypturellus undulatus 1818 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Crypturellus_undulatus)
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
contrib <- loadingplot(dapc1$var.contr, threshold = 0.0015, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
