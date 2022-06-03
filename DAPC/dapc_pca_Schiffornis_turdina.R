library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Schiffornis_turdina
# Schiffornis turdina
# Schiffornis_turdina
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Schiffornis_turdina <- read.structure("Schiffornis_turdina_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2024, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Schiffornis_turdina_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Schiffornis_turdina_AmazonOnly")
Schiffornis_turdina.X <- tab(Schiffornis_turdina, freq=TRUE, NA.method = "mean")
Schiffornis_turdina.X <- scaleGen(Schiffornis_turdina, NA.method = "mean")
sum(is.na(Schiffornis_turdina$tab)) #amount of missing data
Schiffornis_turdina.X[1:2,1:2]
Schiffornis_turdina.pca <- dudi.pca(Schiffornis_turdina.X, scale=FALSE, scannf = FALSE, nf = 4)
Schiffornis_turdina.pca$eig[1]/sum(Schiffornis_turdina.pca$eig) * 100
Schiffornis_turdina.pca$eig[2]/sum(Schiffornis_turdina.pca$eig) * 100

Schiffornis_turdina.pca$li

Schiffornis_turdina.grp <- find.clusters(Schiffornis_turdina, n.pca = 100)

xval <- xvalDapc(Schiffornis_turdina.X, Schiffornis_turdina.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Schiffornis_turdina.X, Schiffornis_turdina.grp$grp)

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
tiff("Schiffornis_turdina_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Schiffornis_turdina_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Schiffornis_turdina.pca$li, Schiffornis_turdina.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-70, 15), ylim=c(-45, 25),
          xlab="PC 1 (14.9%)", ylab="PC 2 (11.7%)",
          main = "PCA of Schiffornis turdina 2024 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Schiffornis_turdina_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Schiffornis_turdina_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Schiffornis_turdina.pca,
             xlim=c(-70, 15), ylim=c(-45, 25),
             repel = TRUE,
             title = "PCA of Schiffornis turdina 2024 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Schiffornis_turdina_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Schiffornis_turdina_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Schiffornis_turdina.pca,
             xlim=c(-70, 15), ylim=c(-45, 25),
             geom = "point",
             title = "PCA of Schiffornis turdina 2024 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Schiffornis_turdina_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Schiffornis_turdina_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Schiffornis_turdina.pca,
             geom = "point",
             xlim=c(-70, 15), ylim=c(-45, 25),
             title = "PCA of Schiffornis turdina 2024 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Schiffornis_turdina)
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
