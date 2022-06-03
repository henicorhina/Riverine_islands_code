library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace cravul with your species
# replace Cranioleuca vulpecula
# replace cranioleuca_vulpecula
# replace 1018 with your number of snps

cravul <- read.structure("cranioleuca_vulpecula_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1018, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/cranioleuca_vulpecula")
cravul.X <- tab(cravul, freq=TRUE, NA.method = "mean")
cravul.X <- scaleGen(cravul, NA.method = "mean")
sum(is.na(cravul$tab)) #amount of missing data
cravul.X[1:2,1:2]
cravul.pca <- dudi.pca(cravul.X, scale=FALSE, scannf = FALSE, nf = 4)
cravul.pca$eig[1]/sum(cravul.pca$eig) * 100
cravul.pca$eig[2]/sum(cravul.pca$eig) * 100

cravul.pca$li

cravul.grp <- find.clusters(cravul, n.pca = 100)

xval <- xvalDapc(cravul.X, cravul.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca
dapc1 <- dapc(cravul.X, cravul.grp$grp)

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
tiff("cranioleuca_vulpecula_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(cravul.pca$li, cravul.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (22.8%)", ylab="PC 2 (17.8%)",
          main = "PCA of Cranioleuca vulpecula 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("cranioleuca_vulpecula_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(cravul.pca,
             xlim=c(-15, 15), ylim=c(-12, 12),
             repel = TRUE,
             title = "PCA of Cranioleuca vulpecula 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("cranioleuca_vulpecula_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(cravul.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Cranioleuca vulpecula 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("cranioleuca_vulpecula_black.tiff", width = 10, height = 10, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(cravul.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-12, 12),
             title = "PCA of Cranioleuca vulpecula 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(cravul)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.003, 
                       thres=.07, lab.jitter=1)
contrib$var.values
compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
