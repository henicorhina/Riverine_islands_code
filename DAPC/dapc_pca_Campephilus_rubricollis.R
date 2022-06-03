library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Campephilus_rubricollis
# Campephilus rubricollis
# Campephilus_rubricollis
# 1760 with your number of snps
# xlim=c(-10, 10), ylim=c(-10, 10), 

Campephilus_rubricollis <- read.structure("Campephilus_rubricollis_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=11, n.loc=1760, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/Campephilus_rubricollis")
Campephilus_rubricollis.X <- tab(Campephilus_rubricollis, freq=TRUE, NA.method = "mean")
Campephilus_rubricollis.X <- scaleGen(Campephilus_rubricollis, NA.method = "mean")
sum(is.na(Campephilus_rubricollis$tab)) #amount of missing data
Campephilus_rubricollis.X[1:2,1:2]
Campephilus_rubricollis.pca <- dudi.pca(Campephilus_rubricollis.X, scale=FALSE, scannf = FALSE, nf = 4)
Campephilus_rubricollis.pca$eig[1]/sum(Campephilus_rubricollis.pca$eig) * 100
Campephilus_rubricollis.pca$eig[2]/sum(Campephilus_rubricollis.pca$eig) * 100

Campephilus_rubricollis.pca$li

Campephilus_rubricollis.grp <- find.clusters(Campephilus_rubricollis, n.pca = 100)

xval <- xvalDapc(Campephilus_rubricollis.X, Campephilus_rubricollis.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Campephilus_rubricollis.X, Campephilus_rubricollis.grp$grp)

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
tiff("Campephilus_rubricollis_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(Campephilus_rubricollis.pca$li, Campephilus_rubricollis.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (22.2%)", ylab="PC 2 (13.8%)",
          main = "PCA of Campephilus rubricollis 1760 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Campephilus_rubricollis_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Campephilus_rubricollis.pca,
      	     xlim=c(-10, 10), ylim=c(-10, 10), 
             repel = TRUE,
             title = "PCA of Campephilus rubricollis 1760 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Campephilus_rubricollis_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Campephilus_rubricollis.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             geom = "point",
             title = "PCA of Campephilus rubricollis 1760 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Campephilus_rubricollis_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(Campephilus_rubricollis.pca,
             geom = "point",
             xlim=c(-10, 10), ylim=c(-10, 10), 
             title = "PCA of Campephilus rubricollis 1760 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Campephilus_rubricollis)
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
