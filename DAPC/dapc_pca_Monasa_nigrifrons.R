library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Monasa_nigrifrons
# Monasa nigrifrons
# Monasa_nigrifrons
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Monasa_nigrifrons <- read.structure("Monasa_nigrifrons_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=1852, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Monasa_nigrifrons_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Monasa_nigrifrons_AmazonOnly")
Monasa_nigrifrons.X <- tab(Monasa_nigrifrons, freq=TRUE, NA.method = "mean")
Monasa_nigrifrons.X <- scaleGen(Monasa_nigrifrons, NA.method = "mean")
sum(is.na(Monasa_nigrifrons$tab)) #amount of missing data
Monasa_nigrifrons.X[1:2,1:2]
Monasa_nigrifrons.pca <- dudi.pca(Monasa_nigrifrons.X, scale=FALSE, scannf = FALSE, nf = 4)
Monasa_nigrifrons.pca$eig[1]/sum(Monasa_nigrifrons.pca$eig) * 100
Monasa_nigrifrons.pca$eig[2]/sum(Monasa_nigrifrons.pca$eig) * 100

Monasa_nigrifrons.pca$li

Monasa_nigrifrons.grp <- find.clusters(Monasa_nigrifrons, n.pca = 100)

xval <- xvalDapc(Monasa_nigrifrons.X, Monasa_nigrifrons.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Monasa_nigrifrons.X, Monasa_nigrifrons.grp$grp)

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
tiff("Monasa_nigrifrons_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Monasa_nigrifrons_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Monasa_nigrifrons.pca$li, Monasa_nigrifrons.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-30, 35), ylim=c(-20, 60), 
          xlab="PC 1 (16.1%)", ylab="PC 2 (11.4%)",
          main = "PCA of Monasa nigrifrons 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Monasa_nigrifrons_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Monasa_nigrifrons_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Monasa_nigrifrons.pca,
             xlim=c(-30, 35), ylim=c(-20, 60), 
             repel = TRUE,
             title = "PCA of Monasa nigrifrons 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Monasa_nigrifrons_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Monasa_nigrifrons_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Monasa_nigrifrons.pca,
             xlim=c(-30, 35), ylim=c(-20, 60), 
             geom = "point",
             title = "PCA of Monasa nigrifrons 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Monasa_nigrifrons_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Monasa_nigrifrons_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Monasa_nigrifrons.pca,
             geom = "point",
             xlim=c(-30, 35), ylim=c(-20, 60), 
             title = "PCA of Monasa nigrifrons 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Monasa_nigrifrons)
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
