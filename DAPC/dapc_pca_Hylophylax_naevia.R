library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Hylophylax_naevia
# Hylophylax naevia
# Hylophylax_naevia
# 2027 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Hylophylax_naevia <- read.structure("Hylophylax_naevia_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=9, n.loc=2027, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Hylophylax_naevia")
Hylophylax_naevia.X <- tab(Hylophylax_naevia, freq=TRUE, NA.method = "mean")
Hylophylax_naevia.X <- scaleGen(Hylophylax_naevia, NA.method = "mean")
sum(is.na(Hylophylax_naevia$tab)) #amount of missing data
Hylophylax_naevia.X[1:2,1:2]
Hylophylax_naevia.pca <- dudi.pca(Hylophylax_naevia.X, scale=FALSE, scannf = FALSE, nf = 4)
Hylophylax_naevia.pca$eig[1]/sum(Hylophylax_naevia.pca$eig) * 100
Hylophylax_naevia.pca$eig[2]/sum(Hylophylax_naevia.pca$eig) * 100

Hylophylax_naevia.pca$li

Hylophylax_naevia.grp <- find.clusters(Hylophylax_naevia, n.pca = 100)

xval <- xvalDapc(Hylophylax_naevia.X, Hylophylax_naevia.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Hylophylax_naevia.X, Hylophylax_naevia.grp$grp)

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
tiff("Hylophylax_naevia_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Hylophylax_naevia.pca$li, Hylophylax_naevia.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-16, 15), 
          xlab="PC 1 (48.4%)", ylab="PC 2 (16.0%)",
          main = "PCA of Hylophylax naevia 2027 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Hylophylax_naevia_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Hylophylax_naevia.pca,
      	     xlim=c(-15, 15), ylim=c(-16, 15), 
             repel = TRUE,
             title = "PCA of Hylophylax naevia 2027 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Hylophylax_naevia_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Hylophylax_naevia.pca,
             xlim=c(-15, 15), ylim=c(-16, 15),
             geom = "point",
             title = "PCA of Hylophylax naevia 2027 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Hylophylax_naevia_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Hylophylax_naevia.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-16, 15), 
             title = "PCA of Hylophylax naevia 2027 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Hylophylax_naevia)
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
contrib <- loadingplot(dapc1$var.contr, threshold = 0.002, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
