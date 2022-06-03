library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

# replace conmar with your species
# replace Conirostrum margaritae
# replace conirostrum_margaritae
# replace 859 with your number of snps

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")
conmar <- read.structure("conirostrum_margaritae_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=4, n.loc=859, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/conirostrum_margaritae")

conmar.X <- tab(conmar, freq=TRUE, NA.method = "mean")
#conmar.X <- scaleGen(conmar, NA.method = "mean")
sum(is.na(conmar$tab)) #amount of missing data
conmar.X[1:2,1:2]
conmar.pca <- dudi.pca(conmar.X, scale=FALSE, scannf = FALSE, nf = 3)
conmar.pca$eig[1]/sum(conmar.pca$eig) * 100
conmar.pca$eig[2]/sum(conmar.pca$eig) * 100

conmar.pca$li

conmar.grp <- find.clusters(conmar, n.pca = 100)

xval <- xvalDapc(conmar.X, conmar.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc1 <- dapc(conmar.X, conmar.grp$grp)

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
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(conmar.pca$li, conmar.pca$li, add.plot=FALSE, cex=3, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (35.3%)", ylab="PC 2 (33.5%)",
          main = "PCA of Conirostrum margaritae 859 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conmar.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             repel = TRUE,
             title = "PCA of Conirostrum margaritae 859 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conmar.pca,
                  #axes = c(1, 2),
                  #xlim=c(-10, 10), ylim=c(-10, 10),
                  #repel = TRUE,     # Avoid text overlapping
                  geom = "point",
                  #habillage=dapc1$grp,
                  #addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2",
                  #ggtheme = theme_minimal(),
                  title = "PCA of Conirostrum margaritae 859 SNPs \naxes 1-2"
                  #geom_point(size = 1)
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conbic_full.pca,
                  geom = "point",
                  title = "PCA of Conirostrum margaritae 859 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(conmar)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
#round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
#assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.003, 
                       thres=.07, lab.jitter=1)

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
