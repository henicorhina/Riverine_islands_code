library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# thamnophilus_cryptoleucus
# Thamnophilus_cryptoleucus
# thamnophilus_cryptoleucus
# 1864 with your number of snps
# xlim=c(-5, 12), ylim=c(-5, 12), 

thamnophilus_cryptoleucus <- read.structure("thamnophilus_cryptoleucus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1864, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/thamnophilus_cryptoleucus")
thamnophilus_cryptoleucus.X <- tab(thamnophilus_cryptoleucus, freq=TRUE, NA.method = "mean")
thamnophilus_cryptoleucus.X <- scaleGen(thamnophilus_cryptoleucus, NA.method = "mean")
sum(is.na(thamnophilus_cryptoleucus$tab)) #amount of missing data
thamnophilus_cryptoleucus.X[1:2,1:2]
thamnophilus_cryptoleucus.pca <- dudi.pca(thamnophilus_cryptoleucus.X, scale=FALSE, scannf = FALSE, nf = 4)
thamnophilus_cryptoleucus.pca$eig[1]/sum(thamnophilus_cryptoleucus.pca$eig) * 100
thamnophilus_cryptoleucus.pca$eig[2]/sum(thamnophilus_cryptoleucus.pca$eig) * 100

thamnophilus_cryptoleucus.pca$li

thamnophilus_cryptoleucus.grp <- find.clusters(thamnophilus_cryptoleucus, n.pca = 100)

xval <- xvalDapc(thamnophilus_cryptoleucus.X, thamnophilus_cryptoleucus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(thamnophilus_cryptoleucus.X, thamnophilus_cryptoleucus.grp$grp)

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
tiff("thamnophilus_cryptoleucus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(thamnophilus_cryptoleucus.pca$li, thamnophilus_cryptoleucus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-5, 12), ylim=c(-5, 12), 
          xlab="PC 1 (28.6%)", ylab="PC 2 (19.8%)",
          main = "PCA of Thamnophilus cryptoleucus 1864 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("thamnophilus_cryptoleucus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_cryptoleucus.pca,
      	     xlim=c(-5, 12), ylim=c(-5, 12), 
             repel = TRUE,
             title = "PCA of Thamnophilus cryptoleucus 1864 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("thamnophilus_cryptoleucus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_cryptoleucus.pca,
             xlim=c(-5, 12), ylim=c(-5, 12),
             geom = "point",
             title = "PCA of Thamnophilus cryptoleucus 1864 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("thamnophilus_cryptoleucus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_cryptoleucus.pca,
             geom = "point",
             xlim=c(-5, 12), ylim=c(-5, 12), 
             title = "PCA of Thamnophilus cryptoleucus 1864 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(thamnophilus_cryptoleucus)
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
