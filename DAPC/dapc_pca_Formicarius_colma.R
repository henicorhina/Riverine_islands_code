library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Formicarius_colma
# Formicarius colma
# Formicarius_colma
# 1901 with your number of snps
# xlim=c(-18, 18), ylim=c(-18, 18), 

Formicarius_colma <- read.structure("Formicarius_colma_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=1901, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Formicarius_colma_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Formicarius_colma_AmazonOnly")
Formicarius_colma.X <- tab(Formicarius_colma, freq=TRUE, NA.method = "mean")
Formicarius_colma.X <- scaleGen(Formicarius_colma, NA.method = "mean")
sum(is.na(Formicarius_colma$tab)) #amount of missing data
Formicarius_colma.X[1:2,1:2]
Formicarius_colma.pca <- dudi.pca(Formicarius_colma.X, scale=FALSE, scannf = FALSE, nf = 4)
Formicarius_colma.pca$eig[1]/sum(Formicarius_colma.pca$eig) * 100
Formicarius_colma.pca$eig[2]/sum(Formicarius_colma.pca$eig) * 100

Formicarius_colma.pca$li

Formicarius_colma.grp <- find.clusters(Formicarius_colma, n.pca = 100)

xval <- xvalDapc(Formicarius_colma.X, Formicarius_colma.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Formicarius_colma.X, Formicarius_colma.grp$grp)

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
tiff("Formicarius_colma_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Formicarius_colma_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Formicarius_colma.pca$li, Formicarius_colma.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-25, 60), ylim=c(-35, 55), 
          xlab="PC 1 (21.8%)", ylab="PC 2 (12.3%)",
          main = "PCA of Formicarius colma 1901 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Formicarius_colma_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Formicarius_colma_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Formicarius_colma.pca,
             xlim=c(-25, 60), ylim=c(-35, 55), 
             repel = TRUE,
             title = "PCA of Formicarius colma 1901 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Formicarius_colma_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Formicarius_colma_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Formicarius_colma.pca,
             xlim=c(-25, 60), ylim=c(-35, 55), 
             geom = "point",
             title = "PCA of Formicarius colma 1901 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Formicarius_colma_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Formicarius_colma_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Formicarius_colma.pca,
             geom = "point",
             xlim=c(-25, 60), ylim=c(-35, 55), 
             title = "PCA of Formicarius colma 1901 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Formicarius_colma)
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
