library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Megascops_choliba
# Megascops choliba
# Megascops_choliba
# 1554 with your number of snps
# xlim=c(-10, 10), ylim=c(-10, 10), 

Megascops_choliba <- read.structure("Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=1554, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Megascops_choliba_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Megascops_choliba_AmazonOnly")
Megascops_choliba.X <- tab(Megascops_choliba, freq=TRUE, NA.method = "mean")
Megascops_choliba.X <- scaleGen(Megascops_choliba, NA.method = "mean")
sum(is.na(Megascops_choliba$tab)) #amount of missing data
Megascops_choliba.X[1:2,1:2]
Megascops_choliba.pca <- dudi.pca(Megascops_choliba.X, scale=FALSE, scannf = FALSE, nf = 4)
Megascops_choliba.pca$eig[1]/sum(Megascops_choliba.pca$eig) * 100
Megascops_choliba.pca$eig[2]/sum(Megascops_choliba.pca$eig) * 100

Megascops_choliba.pca$li

Megascops_choliba.grp <- find.clusters(Megascops_choliba, n.pca = 100)

xval <- xvalDapc(Megascops_choliba.X, Megascops_choliba.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Megascops_choliba.X, Megascops_choliba.grp$grp)

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
tiff("Megascops_choliba_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Megascops_choliba_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Megascops_choliba.pca$li, Megascops_choliba.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-55, 25), ylim=c(-25, 55), 
          xlab="PC 1 (12.2%)", ylab="PC 2 (11.9%)",
          main = "PCA of Megascops choliba 1554 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Megascops_choliba_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Megascops_choliba_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Megascops_choliba.pca,
             xlim=c(-55, 25), ylim=c(-25, 55), 
             repel = TRUE,
             title = "PCA of Megascops choliba 1554 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Megascops_choliba_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Megascops_choliba_groups_forced_to_2.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Megascops_choliba.pca,
             xlim=c(-55, 25), ylim=c(-25, 55), 
             geom = "point",
             title = "PCA of Megascops choliba 1554 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Megascops_choliba_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Megascops_choliba_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Megascops_choliba.pca,
             geom = "point",
             xlim=c(-55, 25), ylim=c(-25, 55), 
             title = "PCA of Megascops choliba 1554 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Megascops_choliba)
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
