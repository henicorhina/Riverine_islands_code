library("adegenet") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Tachyphonus_luctuosus
# Tachyphonus luctuosus
# Tachyphonus_luctuosus
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Tachyphonus_luctuosus <- read.structure("Tachyphonus_luctuosus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2035, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Tachyphonus_luctuosus_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Tachyphonus_luctuosus_AmazonOnly")
Tachyphonus_luctuosus.X <- tab(Tachyphonus_luctuosus, freq=TRUE, NA.method = "mean")
Tachyphonus_luctuosus.X <- scaleGen(Tachyphonus_luctuosus, NA.method = "mean")
sum(is.na(Tachyphonus_luctuosus$tab)) #amount of missing data
Tachyphonus_luctuosus.X[1:2,1:2]
Tachyphonus_luctuosus.pca <- dudi.pca(Tachyphonus_luctuosus.X, scale=FALSE, scannf = FALSE, nf = 4)
Tachyphonus_luctuosus.pca$eig[1]/sum(Tachyphonus_luctuosus.pca$eig) * 100
Tachyphonus_luctuosus.pca$eig[2]/sum(Tachyphonus_luctuosus.pca$eig) * 100

Tachyphonus_luctuosus.pca$li

Tachyphonus_luctuosus.grp <- find.clusters(Tachyphonus_luctuosus, n.pca = 100)

xval <- xvalDapc(Tachyphonus_luctuosus.X, Tachyphonus_luctuosus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Tachyphonus_luctuosus.X, Tachyphonus_luctuosus.grp$grp)

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
tiff("Tachyphonus_luctuosus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_luctuosus_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Tachyphonus_luctuosus.pca$li, Tachyphonus_luctuosus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-30, 50), ylim=c(-50, 40), 
          xlab="PC 1 (12.8%)", ylab="PC 2 (12.1%)",
          main = "PCA of Tachyphonus luctuosus 2035 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Tachyphonus_luctuosus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_luctuosus_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_luctuosus.pca,
             xlim=c(-30, 50), ylim=c(-50, 40), 
             repel = TRUE,
             title = "PCA of Tachyphonus luctuosus 2035 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Tachyphonus_luctuosus_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_luctuosus_groups_forced_to_2.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_luctuosus.pca,
             xlim=c(-30, 50), ylim=c(-50, 40), 
             geom = "point",
             title = "PCA of Tachyphonus luctuosus 2035 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Tachyphonus_luctuosus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_luctuosus_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_luctuosus.pca,
             geom = "point",
             xlim=c(-30, 50), ylim=c(-50, 40), 
             title = "PCA of Tachyphonus luctuosus 2035 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Tachyphonus_luctuosus)
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
