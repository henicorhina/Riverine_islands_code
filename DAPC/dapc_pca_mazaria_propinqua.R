library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# mazaria_propinqua
# Mazaria propinqua
# mazaria_propinqua
# 1413 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

mazaria_propinqua <- read.structure("mazaria_propinqua_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=8, n.loc=1413, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/mazaria_propinqua")
mazaria_propinqua.X <- tab(mazaria_propinqua, freq=TRUE, NA.method = "mean")
mazaria_propinqua.X <- scaleGen(mazaria_propinqua, NA.method = "mean")
sum(is.na(mazaria_propinqua$tab)) #amount of missing data
mazaria_propinqua.X[1:2,1:2]
mazaria_propinqua.pca <- dudi.pca(mazaria_propinqua.X, scale=FALSE, scannf = FALSE, nf = 4)
mazaria_propinqua.pca$eig[1]/sum(mazaria_propinqua.pca$eig) * 100
mazaria_propinqua.pca$eig[2]/sum(mazaria_propinqua.pca$eig) * 100

mazaria_propinqua.pca$li

mazaria_propinqua.grp <- find.clusters(mazaria_propinqua, n.pca = 100)

xval <- xvalDapc(mazaria_propinqua.X, mazaria_propinqua.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca
dapc1 <- dapc(mazaria_propinqua.X, mazaria_propinqua.grp$grp)

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
tiff("mazaria_propinqua_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(mazaria_propinqua.pca$li, mazaria_propinqua.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (24.2%)", ylab="PC 2 (16.7%)",
          main = "PCA of Mazaria propinqua 1413 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("mazaria_propinqua_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(mazaria_propinqua.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Mazaria propinqua 1413 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("mazaria_propinqua_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(mazaria_propinqua.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Mazaria propinqua 1413 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("mazaria_propinqua_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(mazaria_propinqua.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Mazaria propinqua 1413 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(mazaria_propinqua)
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
