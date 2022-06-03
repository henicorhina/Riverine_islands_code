library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Crypturellus_variegatus
# Crypturellus variegatus
# Crypturellus_variegatus
# 1984 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Crypturellus_variegatus <- read.structure("Crypturellus_variegatus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=10, n.loc=1984, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Crypturellus_variegatus")
Crypturellus_variegatus.X <- tab(Crypturellus_variegatus, freq=TRUE, NA.method = "mean")
Crypturellus_variegatus.X <- scaleGen(Crypturellus_variegatus, NA.method = "mean")
sum(is.na(Crypturellus_variegatus$tab)) #amount of missing data
Crypturellus_variegatus.X[1:2,1:2]
Crypturellus_variegatus.pca <- dudi.pca(Crypturellus_variegatus.X, scale=FALSE, scannf = FALSE, nf = 4)
Crypturellus_variegatus.pca$eig[1]/sum(Crypturellus_variegatus.pca$eig) * 100
Crypturellus_variegatus.pca$eig[2]/sum(Crypturellus_variegatus.pca$eig) * 100

Crypturellus_variegatus.pca$li

Crypturellus_variegatus.grp <- find.clusters(Crypturellus_variegatus, n.pca = 100)

xval <- xvalDapc(Crypturellus_variegatus.X, Crypturellus_variegatus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Crypturellus_variegatus.X, Crypturellus_variegatus.grp$grp)

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
tiff("Crypturellus_variegatus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Crypturellus_variegatus.pca$li, Crypturellus_variegatus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (33.1%)", ylab="PC 2 (12.5%)",
          main = "PCA of Crypturellus variegatus 1984 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Crypturellus_variegatus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_variegatus.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Crypturellus variegatus 1984 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Crypturellus_variegatus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_variegatus.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Crypturellus variegatus 1984 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Crypturellus_variegatus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Crypturellus_variegatus.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Crypturellus variegatus 1984 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Crypturellus_variegatus)
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
