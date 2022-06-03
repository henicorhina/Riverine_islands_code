library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# thamnophilus_nigrocinereus
# Thamnophilus nigrocinereus
# thamnophilus_nigrocinereus
# 1864 with your number of snps
# xlim=c(-10, 10), ylim=c(-10, 10), 

thamnophilus_nigrocinereus <- read.structure("thamnophilus_nigrocinereus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=11, n.loc=1864, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/thamnophilus_nigrocinereus")
thamnophilus_nigrocinereus.X <- tab(thamnophilus_nigrocinereus, freq=TRUE, NA.method = "mean")
thamnophilus_nigrocinereus.X <- scaleGen(thamnophilus_nigrocinereus, NA.method = "mean")
sum(is.na(thamnophilus_nigrocinereus$tab)) #amount of missing data
thamnophilus_nigrocinereus.X[1:2,1:2]
thamnophilus_nigrocinereus.pca <- dudi.pca(thamnophilus_nigrocinereus.X, scale=FALSE, scannf = FALSE, nf = 4)
thamnophilus_nigrocinereus.pca$eig[1]/sum(thamnophilus_nigrocinereus.pca$eig) * 100
thamnophilus_nigrocinereus.pca$eig[2]/sum(thamnophilus_nigrocinereus.pca$eig) * 100

thamnophilus_nigrocinereus.pca$li

thamnophilus_nigrocinereus.grp <- find.clusters(thamnophilus_nigrocinereus, n.pca = 100)

xval <- xvalDapc(thamnophilus_nigrocinereus.X, thamnophilus_nigrocinereus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(thamnophilus_nigrocinereus.X, thamnophilus_nigrocinereus.grp$grp)

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
tiff("thamnophilus_nigrocinereus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(thamnophilus_nigrocinereus.pca$li, thamnophilus_nigrocinereus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (23.6%)", ylab="PC 2 (15.3%)",
          main = "PCA of Thamnophilus nigrocinereus 1864 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("thamnophilus_nigrocinereus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_nigrocinereus.pca,
      	     xlim=c(-10, 10), ylim=c(-10, 10), 
             repel = TRUE,
             title = "PCA of Thamnophilus nigrocinereus 1864 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("thamnophilus_nigrocinereus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_nigrocinereus.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             geom = "point",
             title = "PCA of Thamnophilus nigrocinereus 1864 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("thamnophilus_nigrocinereus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(thamnophilus_nigrocinereus.pca,
             geom = "point",
             xlim=c(-10, 10), ylim=c(-10, 10), 
             title = "PCA of Thamnophilus nigrocinereus 1864 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(thamnophilus_nigrocinereus)
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
contrib <- loadingplot(dapc1$var.contr, threshold = 0.0035, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
