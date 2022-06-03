library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Tachyphonus_cristatus
# Tachyphonus cristatus
# Tachyphonus_cristatus
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Tachyphonus_cristatus <- read.structure("Tachyphonus_cristatus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=8, n.loc=1997, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Tachyphonus_cristatus_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Tachyphonus_cristatus_AmazonOnly")
Tachyphonus_cristatus.X <- tab(Tachyphonus_cristatus, freq=TRUE, NA.method = "mean")
Tachyphonus_cristatus.X <- scaleGen(Tachyphonus_cristatus, NA.method = "mean")
sum(is.na(Tachyphonus_cristatus$tab)) #amount of missing data
Tachyphonus_cristatus.X[1:2,1:2]
Tachyphonus_cristatus.pca <- dudi.pca(Tachyphonus_cristatus.X, scale=FALSE, scannf = FALSE, nf = 4)
Tachyphonus_cristatus.pca$eig[1]/sum(Tachyphonus_cristatus.pca$eig) * 100
Tachyphonus_cristatus.pca$eig[2]/sum(Tachyphonus_cristatus.pca$eig) * 100

Tachyphonus_cristatus.pca$li

Tachyphonus_cristatus.grp <- find.clusters(Tachyphonus_cristatus, n.pca = 100)

xval <- xvalDapc(Tachyphonus_cristatus.X, Tachyphonus_cristatus.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Tachyphonus_cristatus.X, Tachyphonus_cristatus.grp$grp)

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
tiff("Tachyphonus_cristatus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_cristatus_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Tachyphonus_cristatus.pca$li, Tachyphonus_cristatus.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-30, 55), ylim=c(-60, 35), 
          xlab="PC 1 (17.4%)", ylab="PC 2 (15.4%)",
          main = "PCA of Tachyphonus cristatus 1997 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Tachyphonus_cristatus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_cristatus_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_cristatus.pca,
             xlim=c(-30, 55), ylim=c(-60, 35), 
             repel = TRUE,
             title = "PCA of Tachyphonus cristatus 1997 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Tachyphonus_cristatus_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_cristatus_groups_forced_to_2.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_cristatus.pca,
             xlim=c(-30, 55), ylim=c(-60, 35), 
             geom = "point",
             title = "PCA of Tachyphonus cristatus 1997 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Tachyphonus_cristatus_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Tachyphonus_cristatus_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Tachyphonus_cristatus.pca,
             geom = "point",
             xlim=c(-30, 55), ylim=c(-60, 35), 
             title = "PCA of Tachyphonus cristatus 1997 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Tachyphonus_cristatus)
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
