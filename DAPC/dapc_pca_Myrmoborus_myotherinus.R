library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Myrmoborus_myotherinus
# Myrmoborus myotherinus
# Myrmoborus_myotherinus
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Myrmoborus_myotherinus <- read.structure("Myrmoborus_myotherinus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1018, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Myrmoborus_myotherinus")
Myrmoborus_myotherinus.X <- tab(Myrmoborus_myotherinus, freq=TRUE, NA.method = "mean")
Myrmoborus_myotherinus.X <- scaleGen(Myrmoborus_myotherinus, NA.method = "mean")
sum(is.na(Myrmoborus_myotherinus$tab)) #amount of missing data
Myrmoborus_myotherinus.X[1:2,1:2]
Myrmoborus_myotherinus.pca <- dudi.pca(Myrmoborus_myotherinus.X, scale=FALSE, scannf = FALSE, nf = 4)
Myrmoborus_myotherinus.pca$eig[1]/sum(Myrmoborus_myotherinus.pca$eig) * 100
Myrmoborus_myotherinus.pca$eig[2]/sum(Myrmoborus_myotherinus.pca$eig) * 100

Myrmoborus_myotherinus.pca$li

Myrmoborus_myotherinus.grp <- find.clusters(Myrmoborus_myotherinus, n.pca = 100)

xval <- xvalDapc(Myrmoborus_myotherinus.X, Myrmoborus_myotherinus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Myrmoborus_myotherinus.X, Myrmoborus_myotherinus.grp$grp)

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
tiff("Myrmoborus_myotherinus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Myrmoborus_myotherinus.pca$li, Myrmoborus_myotherinus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (35.3%)", ylab="PC 2 (33.5%)",
          main = "PCA of Myrmoborus myotherinus 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Myrmoborus_myotherinus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Myrmoborus_myotherinus.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Myrmoborus myotherinus 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Myrmoborus_myotherinus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Myrmoborus_myotherinus.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Myrmoborus myotherinus 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Myrmoborus_myotherinus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Myrmoborus_myotherinus.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Myrmoborus myotherinus 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Myrmoborus_myotherinus)
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
