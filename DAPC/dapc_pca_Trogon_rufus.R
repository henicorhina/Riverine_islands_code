library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Trogon_rufus
# Trogon rufus
# Trogon_rufus
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Trogon_rufus <- read.structure("Trogon_rufus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2074, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Trogon_rufus_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Trogon_rufus_AmazonOnly")
Trogon_rufus.X <- tab(Trogon_rufus, freq=TRUE, NA.method = "mean")
Trogon_rufus.X <- scaleGen(Trogon_rufus, NA.method = "mean")
sum(is.na(Trogon_rufus$tab)) #amount of missing data
Trogon_rufus.X[1:2,1:2]
Trogon_rufus.pca <- dudi.pca(Trogon_rufus.X, scale=FALSE, scannf = FALSE, nf = 4)
Trogon_rufus.pca$eig[1]/sum(Trogon_rufus.pca$eig) * 100
Trogon_rufus.pca$eig[2]/sum(Trogon_rufus.pca$eig) * 100

Trogon_rufus.pca$li

Trogon_rufus.grp <- find.clusters(Trogon_rufus, n.pca = 100)

xval <- xvalDapc(Trogon_rufus.X, Trogon_rufus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Trogon_rufus.X, Trogon_rufus.grp$grp)

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
tiff("Trogon_rufus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Trogon_rufus_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Trogon_rufus.pca$li, Trogon_rufus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-30, 50), ylim=c(-30, 50), 
          xlab="PC 1 (12.7%)", ylab="PC 2 (11.4%)",
          main = "PCA of Trogon rufus 2074 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Trogon_rufus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Trogon_rufus_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Trogon_rufus.pca,
             xlim=c(-30, 50), ylim=c(-30, 50), 
             repel = TRUE,
             title = "PCA of Trogon rufus 2074 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Trogon_rufus_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Trogon_rufus_groups_forced_to_2.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Trogon_rufus.pca,
             xlim=c(-30, 50), ylim=c(-30, 50), 
             geom = "point",
             title = "PCA of Trogon rufus 2074 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Trogon_rufus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Trogon_rufus_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Trogon_rufus.pca,
             geom = "point",
             xlim=c(-30, 50), ylim=c(-30, 50), 
             title = "PCA of Trogon rufus 2074 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Trogon_rufus)
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
