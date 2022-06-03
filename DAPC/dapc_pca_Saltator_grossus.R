library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Saltator_grossus
# Saltator grossus
# Saltator_grossus
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Saltator_grossus <- read.structure("Saltator_grossus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2021, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Saltator_grossus_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Saltator_grossus_AmazonOnly")
Saltator_grossus.X <- tab(Saltator_grossus, freq=TRUE, NA.method = "mean")
Saltator_grossus.X <- scaleGen(Saltator_grossus, NA.method = "mean")
sum(is.na(Saltator_grossus$tab)) #amount of missing data
Saltator_grossus.X[1:2,1:2]
Saltator_grossus.pca <- dudi.pca(Saltator_grossus.X, scale=FALSE, scannf = FALSE, nf = 4)
Saltator_grossus.pca$eig[1]/sum(Saltator_grossus.pca$eig) * 100
Saltator_grossus.pca$eig[2]/sum(Saltator_grossus.pca$eig) * 100

Saltator_grossus.pca$li

Saltator_grossus.grp <- find.clusters(Saltator_grossus, n.pca = 100)

xval <- xvalDapc(Saltator_grossus.X, Saltator_grossus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Saltator_grossus.X, Saltator_grossus.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Plot .png", width = 6, height = 6, units = 'in', 
    res = 800)
bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
       type="tifflzw", res=800)
# run the desired format above then the plot, then this dev.off command
dev.off()
par(mfrow = c(1,1))


# pca with colored dots
tiff("Saltator_grossus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_grossus_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Saltator_grossus.pca$li, Saltator_grossus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-45, 25), ylim=c(-30, 45), 
          xlab="PC 1 (14.0%)", ylab="PC 2 (12.0%)",
          main = "PCA of Saltator grossus 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Saltator_grossus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_grossus_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_grossus.pca,
             xlim=c(-45, 25), ylim=c(-30, 45), 
             repel = TRUE,
             title = "PCA of Saltator grossus 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Saltator_grossus_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_grossus_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_grossus.pca,
             xlim=c(-45, 25), ylim=c(-30, 45), 
             geom = "point",
             title = "PCA of Saltator grossus 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Saltator_grossus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_grossus_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_grossus.pca,
             geom = "point",
             xlim=c(-45, 25), ylim=c(-30, 45), 
             title = "PCA of Saltator grossus 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Saltator_grossus)
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
