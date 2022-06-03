library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# serpophaga_hypoleuca
# Serpophaga hypoleuca
# serpophaga_hypoleuca
# 1156 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12), 

serpophaga_hypoleuca <- read.structure("serpophaga_hypoleuca_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=9, n.loc=1156, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/serpophaga_hypoleuca")
serpophaga_hypoleuca.X <- tab(serpophaga_hypoleuca, freq=TRUE, NA.method = "mean")
serpophaga_hypoleuca.X <- scaleGen(serpophaga_hypoleuca, NA.method = "mean")
sum(is.na(serpophaga_hypoleuca$tab)) #amount of missing data
serpophaga_hypoleuca.X[1:2,1:2]
serpophaga_hypoleuca.pca <- dudi.pca(serpophaga_hypoleuca.X, scale=FALSE, scannf = FALSE, nf = 4)
serpophaga_hypoleuca.pca$eig[1]/sum(serpophaga_hypoleuca.pca$eig) * 100
serpophaga_hypoleuca.pca$eig[2]/sum(serpophaga_hypoleuca.pca$eig) * 100

serpophaga_hypoleuca.pca$li

serpophaga_hypoleuca.grp <- find.clusters(serpophaga_hypoleuca, n.pca = 100)

xval <- xvalDapc(serpophaga_hypoleuca.X, serpophaga_hypoleuca.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(serpophaga_hypoleuca.X, serpophaga_hypoleuca.grp$grp)

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
tiff("serpophaga_hypoleuca_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(serpophaga_hypoleuca.pca$li, serpophaga_hypoleuca.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (17.0%)", ylab="PC 2 (16.6%)",
          main = "PCA of Serpophaga hypoleuca 1156 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("serpophaga_hypoleuca_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(serpophaga_hypoleuca.pca,
      	     xlim=c(-12, 12), ylim=c(-12, 12), 
             repel = TRUE,
             title = "PCA of Serpophaga hypoleuca 1156 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("serpophaga_hypoleuca_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(serpophaga_hypoleuca.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Serpophaga hypoleuca 1156 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("serpophaga_hypoleuca_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(serpophaga_hypoleuca.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12), 
             title = "PCA of Serpophaga hypoleuca 1156 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(serpophaga_hypoleuca)
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
