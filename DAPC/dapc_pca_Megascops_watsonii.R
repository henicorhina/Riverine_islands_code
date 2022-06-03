library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Megascops_watsonii
# Megascops watsonii
# Megascops_watsonii
# 1872 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12), 

Megascops_watsonii <- read.structure("Megascops_watsonii_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=11, n.loc=1872, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Megascops_watsonii")
Megascops_watsonii.X <- tab(Megascops_watsonii, freq=TRUE, NA.method = "mean")
Megascops_watsonii.X <- scaleGen(Megascops_watsonii, NA.method = "mean")
sum(is.na(Megascops_watsonii$tab)) #amount of missing data
Megascops_watsonii.X[1:2,1:2]
Megascops_watsonii.pca <- dudi.pca(Megascops_watsonii.X, scale=FALSE, scannf = FALSE, nf = 4)
Megascops_watsonii.pca$eig[1]/sum(Megascops_watsonii.pca$eig) * 100
Megascops_watsonii.pca$eig[2]/sum(Megascops_watsonii.pca$eig) * 100

Megascops_watsonii.pca$li

Megascops_watsonii.grp <- find.clusters(Megascops_watsonii, n.pca = 100)

xval <- xvalDapc(Megascops_watsonii.X, Megascops_watsonii.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Megascops_watsonii.X, Megascops_watsonii.grp$grp)

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
tiff("Megascops_watsonii_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Megascops_watsonii.pca$li, Megascops_watsonii.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (15.2%)", ylab="PC 2 (10.8%)",
          main = "PCA of Megascops watsonii 1872 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Megascops_watsonii_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Megascops_watsonii.pca,
      	     xlim=c(-12, 12), ylim=c(-12, 12), 
             repel = TRUE,
             title = "PCA of Megascops watsonii 1872 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Megascops_watsonii_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Megascops_watsonii.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Megascops watsonii 1872 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Megascops_watsonii_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Megascops_watsonii.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12), 
             title = "PCA of Megascops watsonii 1872 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Megascops_watsonii)
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
