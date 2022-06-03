library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# myrmotherula_assimilis
# Myrmotherula assimilis
# myrmotherula_assimilis
# 1508 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

myrmotherula_assimilis <- read.structure("myrmotherula_assimilis_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1508, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/myrmotherula_assimilis")
myrmotherula_assimilis.X <- tab(myrmotherula_assimilis, freq=TRUE, NA.method = "mean")
myrmotherula_assimilis.X <- scaleGen(myrmotherula_assimilis, NA.method = "mean")
sum(is.na(myrmotherula_assimilis$tab)) #amount of missing data
myrmotherula_assimilis.X[1:2,1:2]
myrmotherula_assimilis.pca <- dudi.pca(myrmotherula_assimilis.X, scale=FALSE, scannf = FALSE, nf = 4)
myrmotherula_assimilis.pca$eig[1]/sum(myrmotherula_assimilis.pca$eig) * 100
myrmotherula_assimilis.pca$eig[2]/sum(myrmotherula_assimilis.pca$eig) * 100

myrmotherula_assimilis.pca$li

myrmotherula_assimilis.grp <- find.clusters(myrmotherula_assimilis, n.pca = 100)

xval <- xvalDapc(myrmotherula_assimilis.X, myrmotherula_assimilis.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(myrmotherula_assimilis.X, myrmotherula_assimilis.grp$grp)

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
tiff("myrmotherula_assimilis_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(myrmotherula_assimilis.pca$li, myrmotherula_assimilis.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (23.8%)", ylab="PC 2 (17.9%)",
          main = "PCA of Myrmotherula assimilis 1508 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("myrmotherula_assimilis_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmotherula_assimilis.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Myrmotherula assimilis 1508 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("myrmotherula_assimilis_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmotherula_assimilis.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Myrmotherula assimilis 1508 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("myrmotherula_assimilis_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmotherula_assimilis.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Myrmotherula assimilis 1508 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(myrmotherula_assimilis)
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
