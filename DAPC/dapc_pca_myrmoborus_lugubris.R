library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# myrmoborus_lugubris
# Myrmoborus lugubris
# myrmoborus_lugubris
# 1539 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

myrmoborus_lugubris <- read.structure("myrmoborus_lugubris_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=8, n.loc=1539, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/myrmoborus_lugubris")
myrmoborus_lugubris.X <- tab(myrmoborus_lugubris, freq=TRUE, NA.method = "mean")
myrmoborus_lugubris.X <- scaleGen(myrmoborus_lugubris, NA.method = "mean")
sum(is.na(myrmoborus_lugubris$tab)) #amount of missing data
myrmoborus_lugubris.X[1:2,1:2]
myrmoborus_lugubris.pca <- dudi.pca(myrmoborus_lugubris.X, scale=FALSE, scannf = FALSE, nf = 4)
myrmoborus_lugubris.pca$eig[1]/sum(myrmoborus_lugubris.pca$eig) * 100
myrmoborus_lugubris.pca$eig[2]/sum(myrmoborus_lugubris.pca$eig) * 100

myrmoborus_lugubris.pca$li

myrmoborus_lugubris.grp <- find.clusters(myrmoborus_lugubris, n.pca = 100)

xval <- xvalDapc(myrmoborus_lugubris.X, myrmoborus_lugubris.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca
dapc1 <- dapc(myrmoborus_lugubris.X, myrmoborus_lugubris.grp$grp)

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
tiff("myrmoborus_lugubris_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(myrmoborus_lugubris.pca$li, myrmoborus_lugubris.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 15), ylim=c(-15, 15), 
          xlab="PC 1 (29.3%)", ylab="PC 2 (17.8%)",
          main = "PCA of Myrmoborus lugubris 1539 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("myrmoborus_lugubris_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmoborus_lugubris.pca,
      	     xlim=c(-15, 15), ylim=c(-15, 15), 
             repel = TRUE,
             title = "PCA of Myrmoborus lugubris 1539 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("myrmoborus_lugubris_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmoborus_lugubris.pca,
             xlim=c(-15, 15), ylim=c(-15, 15),
             geom = "point",
             title = "PCA of Myrmoborus lugubris 1539 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("myrmoborus_lugubris_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrmoborus_lugubris.pca,
             geom = "point",
             xlim=c(-15, 15), ylim=c(-15, 15), 
             title = "PCA of Myrmoborus lugubris 1539 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(myrmoborus_lugubris)
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
