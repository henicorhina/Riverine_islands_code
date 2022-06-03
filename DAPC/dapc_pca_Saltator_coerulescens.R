library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Saltator_coerulescens
# Saltator coerulescens
# Saltator_coerulescens
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Saltator_coerulescens <- read.structure("Saltator_coerulescens_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2065, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Saltator_coerulescens_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Saltator_coerulescens_AmazonOnly")
Saltator_coerulescens.X <- tab(Saltator_coerulescens, freq=TRUE, NA.method = "mean")
Saltator_coerulescens.X <- scaleGen(Saltator_coerulescens, NA.method = "mean")
sum(is.na(Saltator_coerulescens$tab)) #amount of missing data
Saltator_coerulescens.X[1:2,1:2]
Saltator_coerulescens.pca <- dudi.pca(Saltator_coerulescens.X, scale=FALSE, scannf = FALSE, nf = 4)
Saltator_coerulescens.pca$eig[1]/sum(Saltator_coerulescens.pca$eig) * 100
Saltator_coerulescens.pca$eig[2]/sum(Saltator_coerulescens.pca$eig) * 100

Saltator_coerulescens.pca$li

Saltator_coerulescens.grp <- find.clusters(Saltator_coerulescens, n.pca = 100)

xval <- xvalDapc(Saltator_coerulescens.X, Saltator_coerulescens.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Saltator_coerulescens.X, Saltator_coerulescens.grp$grp)

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
tiff("Saltator_coerulescens_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_coerulescens_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Saltator_coerulescens.pca$li, Saltator_coerulescens.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-20, 65), ylim=c(-25, 40), 
          xlab="PC 1 (23.0%)", ylab="PC 2 (9.7%)",
          main = "PCA of Saltator coerulescens 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Saltator_coerulescens_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_coerulescens_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_coerulescens.pca,
             xlim=c(-20, 65), ylim=c(-25, 40), 
             repel = TRUE,
             title = "PCA of Saltator coerulescens 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Saltator_coerulescens_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_coerulescens_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_coerulescens.pca,
             xlim=c(-20, 65), ylim=c(-25, 40), 
             geom = "point",
             title = "PCA of Saltator coerulescens 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Saltator_coerulescens_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Saltator_coerulescens_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Saltator_coerulescens.pca,
             geom = "point",
             xlim=c(-20, 65), ylim=c(-25, 40), 
             title = "PCA of Saltator coerulescens 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Saltator_coerulescens)
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
