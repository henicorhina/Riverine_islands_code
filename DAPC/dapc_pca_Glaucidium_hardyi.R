library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Glaucidium_hardyi
# Glaucidium hardyi
# Glaucidium_hardyi
# 1561 with your number of snps
# xlim=c(-10, 10), ylim=c(-10, 10), 

Glaucidium_hardyi <- read.structure("Glaucidium_hardyi_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=10, n.loc=1561, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Glaucidium_hardyi")
Glaucidium_hardyi.X <- tab(Glaucidium_hardyi, freq=TRUE, NA.method = "mean")
Glaucidium_hardyi.X <- scaleGen(Glaucidium_hardyi, NA.method = "mean")
sum(is.na(Glaucidium_hardyi$tab)) #amount of missing data
Glaucidium_hardyi.X[1:2,1:2]
Glaucidium_hardyi.pca <- dudi.pca(Glaucidium_hardyi.X, scale=FALSE, scannf = FALSE, nf = 4)
Glaucidium_hardyi.pca$eig[1]/sum(Glaucidium_hardyi.pca$eig) * 100
Glaucidium_hardyi.pca$eig[2]/sum(Glaucidium_hardyi.pca$eig) * 100

Glaucidium_hardyi.pca$li

Glaucidium_hardyi.grp <- find.clusters(Glaucidium_hardyi, n.pca = 100)

xval <- xvalDapc(Glaucidium_hardyi.X, Glaucidium_hardyi.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Glaucidium_hardyi.X, Glaucidium_hardyi.grp$grp)

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
tiff("Glaucidium_hardyi_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Glaucidium_hardyi.pca$li, Glaucidium_hardyi.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (22.3%)", ylab="PC 2 (12.7%)",
          main = "PCA of Glaucidium hardyi 1561 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Glaucidium_hardyi_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Glaucidium_hardyi.pca,
      	     xlim=c(-10, 10), ylim=c(-10, 10), 
             repel = TRUE,
             title = "PCA of Glaucidium hardyi 1561 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Glaucidium_hardyi_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Glaucidium_hardyi.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             geom = "point",
             title = "PCA of Glaucidium hardyi 1561 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Glaucidium_hardyi_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Glaucidium_hardyi.pca,
             geom = "point",
             xlim=c(-10, 10), ylim=c(-10, 10), 
             title = "PCA of Glaucidium hardyi 1561 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Glaucidium_hardyi)
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
