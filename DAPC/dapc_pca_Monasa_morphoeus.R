library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Monasa_morphoeus
# Monasa morphoeus
# Monasa_morphoeus
# 1692 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12), 

Monasa_morphoeus <- read.structure("Monasa_morphoeus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=10, n.loc=1692, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Monasa_morphoeus")
Monasa_morphoeus.X <- tab(Monasa_morphoeus, freq=TRUE, NA.method = "mean")
Monasa_morphoeus.X <- scaleGen(Monasa_morphoeus, NA.method = "mean")
sum(is.na(Monasa_morphoeus$tab)) #amount of missing data
Monasa_morphoeus.X[1:2,1:2]
Monasa_morphoeus.pca <- dudi.pca(Monasa_morphoeus.X, scale=FALSE, scannf = FALSE, nf = 4)
Monasa_morphoeus.pca$eig[1]/sum(Monasa_morphoeus.pca$eig) * 100
Monasa_morphoeus.pca$eig[2]/sum(Monasa_morphoeus.pca$eig) * 100

Monasa_morphoeus.pca$li

Monasa_morphoeus.grp <- find.clusters(Monasa_morphoeus, n.pca = 100)

xval <- xvalDapc(Monasa_morphoeus.X, Monasa_morphoeus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Monasa_morphoeus.X, Monasa_morphoeus.grp$grp)

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
tiff("Monasa_morphoeus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Monasa_morphoeus.pca$li, Monasa_morphoeus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (28.5%)", ylab="PC 2 (18.7%)",
          main = "PCA of Monasa morphoeus 1692 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Monasa_morphoeus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Monasa_morphoeus.pca,
      	     xlim=c(-12, 12), ylim=c(-12, 12), 
             repel = TRUE,
             title = "PCA of Monasa morphoeus 1692 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Monasa_morphoeus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Monasa_morphoeus.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Monasa morphoeus 1692 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Monasa_morphoeus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Monasa_morphoeus.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12), 
             title = "PCA of Monasa morphoeus 1692 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Monasa_morphoeus)
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
