library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# elaenia_pelzelni
# Elaenia pelzelni
# elaenia_pelzelni
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

elaenia_pelzelni <- read.structure("elaenia_pelzelni_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=2, n.loc=1036, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/elaenia_pelzelni")
elaenia_pelzelni.X <- tab(elaenia_pelzelni, freq=TRUE, NA.method = "mean")
elaenia_pelzelni.X <- scaleGen(elaenia_pelzelni, NA.method = "mean")
sum(is.na(elaenia_pelzelni$tab)) #amount of missing data
elaenia_pelzelni.X[1:2,1:2]
elaenia_pelzelni.pca <- dudi.pca(elaenia_pelzelni.X, scale=FALSE, scannf = FALSE, nf = 4)
elaenia_pelzelni.pca$eig[1]/sum(elaenia_pelzelni.pca$eig) * 100
elaenia_pelzelni.pca$eig[2]/sum(elaenia_pelzelni.pca$eig) * 100

elaenia_pelzelni.pca$li

elaenia_pelzelni.grp <- find.clusters(elaenia_pelzelni, n.pca = 100)

xval <- xvalDapc(elaenia_pelzelni.X, elaenia_pelzelni.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc1 <- dapc(elaenia_pelzelni.X, elaenia_pelzelni.grp$grp)

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
tiff("elaenia_pelzelni_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(elaenia_pelzelni.pca$li, elaenia_pelzelni.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (100%)", ylab="PC 2 (0%)",
          main = "PCA of Elaenia pelzelni 1036 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("elaenia_pelzelni_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(elaenia_pelzelni.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Elaenia pelzelni 1036 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("elaenia_pelzelni_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(elaenia_pelzelni.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Elaenia pelzelni 1036 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("elaenia_pelzelni_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(elaenia_pelzelni.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Elaenia pelzelni 1036 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(elaenia_pelzelni)
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
