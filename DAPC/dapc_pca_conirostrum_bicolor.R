library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#setwd("/Volumes/Brumfield_Lab_Drive/data/by_species_take2/conirostrum_bicolor/14_formatted_output/adegenet/")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

#conbic_full <- read.structure("conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_random_adegenet_all.stru", n.ind=15, n.loc=1978, onerowperind=FALSE, col.lab=1, col.pop=0, row.marknames=0)
conbic_full <- read.structure("conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", n.ind=15, n.loc=1978, onerowperind=FALSE, col.lab=1, col.pop=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/conirostrum_bicolor_full/")

conbic_full.X <- tab(conbic_full, freq=TRUE, NA.method = "mean")
#conbic_full.X <- scaleGen(conbic_full, NA.method = "mean")
sum(is.na(conbic_full$tab)) #amount of missing data
conbic_full.X[1:2,1:2]
conbic_full.pca <- dudi.pca(conbic_full.X, scale=FALSE, scannf = FALSE, nf = 4)
conbic_full.pca$eig[1]/sum(conbic_full.pca$eig) * 100
conbic_full.pca$eig[2]/sum(conbic_full.pca$eig) * 100

conbic_full.pca$li



#use these to save as high quality figures
tiff("Plot.tiff", width = 10, height = 10, units = 'in', 
     res = 500, compression = 'none')
png("Plot.png", width = 6, height = 6, units = 'in', 
    res = 800)
bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
       type="tifflzw", res=800)

# run the desired format above then the plot, then this dev.off command
dev.off()
par(mfrow = c(1,1))


# pca with colored dots
colorplot(conbic_full.pca$li, conbic_full.pca$li, add.plot=FALSE, cex=3, xlim=c(-10, 10), ylim=c(-10, 10), xlab="PC 1 (28.1%)", ylab="PC 2 (11.0%)")
abline(v=0,h=0,col="grey", lty=2)
title("PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2")

# pca with labeled points
fviz_pca_ind(conbic_full.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             repel = TRUE,
             title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
)

# pca without labeled points - groups
x <- fviz_pca_ind(conbic_full.pca,
             #axes = c(1, 2),
             #xlim=c(-10, 10), ylim=c(-10, 10),
             #repel = TRUE,     # Avoid text overlapping
             geom = "point",
             habillage=dapc1$grp,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2",
             #ggtheme = theme_minimal(),
             title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
             #geom_point(size = 1)
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)


# pca without labeled points - black
x <- fviz_pca_ind(conbic_full.pca,
                  geom = "point",
                  title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
                  ) + geom_point(size = 3)

x + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
x + theme(
  plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
  ,panel.border = element_blank()
)

x + geom_point(aes(colour = factor(dapc1$grp)), size = 3) + theme_bw() +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
       #theme(legend.title=element_blank()) +
       scale_color_discrete(name="",
                breaks=c("1", "2"),
                labels=c("Coastal", "Amazon"))

conbic_full.grp <- find.clusters(conbic_full, n.pca = 100)

xval <- xvalDapc(conbic_full.X, conbic_full.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

# save these to a text file
indNames(conbic_full)
xval$DAPC$assign
xval[4]
xval[6]
#xval[2:6]

dapc1 <- dapc(conbic_full.X, conbic_full.grp$grp)

# group assignment probabilities
round(xval$DAPC$posterior,6)
summary(xval$DAPC)
assignplot(xval$DAPC)

temp1 <- optim.a.score(dapc1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.003, 
                       thres=.07, lab.jitter=1)

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)








# for conbic, not conbic_full
conbic <- read.structure("conirostrum_bicolor_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1978, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/conirostrum_bicolor/")

conbic.X <- tab(conbic, freq=TRUE, NA.method = "mean")
#conbic.X <- scaleGen(conbic, NA.method = "mean")
sum(is.na(conbic$tab)) #amount of missing data
conbic.X[1:2,1:2]
conbic.pca <- dudi.pca(conbic.X, scale=FALSE)
conbic.pca$eig[1]/sum(conbic.pca$eig) * 100
conbic.pca$eig[2]/sum(conbic.pca$eig) * 100

#s.label(conbic.pca$li, xlim=c(-10, 10), ylim=c(-10, 10), clabel=.6, cpoint=0, boxes = FALSE, grid=FALSE)

# pca with colors
tiff("conirostrum_bicolor_color.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(conbic.pca$li, conbic.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-10, 10), ylim=c(-10, 10), 
          xlab="PC 1 (35.3%)", ylab="PC 2 (33.5%)",
          main = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("conirostrum_bicolor_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conbic.pca,
             xlim=c(-10, 10), ylim=c(-10, 10),
             repel = TRUE,
             title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("conirostrum_bicolor_groups.tiff", width = 10, height = 10, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conbic.pca,
             #axes = c(1, 2),
             xlim=c(-10, 10), ylim=c(-10, 10),
             #repel = TRUE,     # Avoid text overlapping
             geom = "point",
             #habillage=dapc1$grp,
             #addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2",
             #ggtheme = theme_minimal(),
             title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
             #geom_point(size = 1)
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("conirostrum_bicolor_black.tiff", width = 10, height = 10, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(conbic.pca,
             geom = "point",
             xlim=c(-10, 10), ylim=c(-10, 10),
             title = "PCA of Conirostrum bicolor 1978 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

#conbic.pca$li
conbic.grp <- find.clusters(conbic, n.pca = 100)
dapc1 <- dapc(conbic.X, conbic.grp$grp)
#scatter(dapc1)

xval <- xvalDapc(conbic.X, conbic.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

#conbic.pca <- dudi.pca(conbic.X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
#barplot(conbic.pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

#s.class(conbic.pca$li, pop(conbic))

# save these to a text file
indNames(conbic)
xval$DAPC$assign
xval[4]
xval[6]

# group assignment probabilities
round(xval$DAPC$posterior,6)
summary(xval$DAPC)
assignplot(xval$DAPC)


#dapc2 <- dapc(conbic.X, conbic.grp$grp, n.da=100, n.pca=10)
#temp <- a.score(dapc1)
temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score

#contrib <- loadingplot(dapc1$var.contr, 
                       #thres=.07, lab.jitter=1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.003, 
                       thres=.07, lab.jitter=1)

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign