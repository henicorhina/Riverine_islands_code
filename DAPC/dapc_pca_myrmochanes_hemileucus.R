library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace with your species:
# myrhem 
# Myrmochanes hemileucus
# myrmochanes_hemileucus
# 898 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12),

myrhem <- read.structure("myrmochanes_hemileucus_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=898, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/myrmochanes_hemileucus")
myrhem.X <- tab(myrhem, freq=TRUE, NA.method = "mean")
myrhem.X <- scaleGen(myrhem, NA.method = "mean", center=FALSE)
sum(is.na(myrhem$tab)) #amount of missing data
myrhem.X[1:2,1:2]
myrhem.pca <- dudi.pca(myrhem.X, scale=FALSE, scannf = FALSE, nf = 4)
myrhem.pca$eig[1]/sum(myrhem.pca$eig) * 100
myrhem.pca$eig[2]/sum(myrhem.pca$eig) * 100

myrhem.pca$li

myrhem.grp <- find.clusters(myrhem, n.pca = 100)

xval <- xvalDapc(myrhem.X, myrhem.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(myrhem.X, myrhem.grp$grp)

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
tiff("myrmochanes_hemileucus_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(myrhem.pca$li, myrhem.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (22.3%)", ylab="PC 2 (17.8%)",
          main = "PCA of Myrmochanes hemileucus 898 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("myrmochanes_hemileucus_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrhem.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             repel = TRUE,
             title = "PCA of Myrmochanes hemileucus 898 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("myrmochanes_hemileucus_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrhem.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Myrmochanes hemileucus 898 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("myrmochanes_hemileucus_black.tiff", width = 10, height = 10, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(myrhem.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12),
             title = "PCA of Myrmochanes hemileucus 898 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(myrhem)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.0025, 
                       thres=.07, lab.jitter=1)

contrib$var.values
compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
