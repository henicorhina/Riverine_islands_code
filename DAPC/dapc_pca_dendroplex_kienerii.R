library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace denkie with your species
# replace Dendroplex kienerii
# replace dendroplex_kienerii
# replace 1306 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12),

denkie <- read.structure("dendroplex_kienerii_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.str", 
                         n.ind=7, n.loc=1306, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/data/1_analysis/dapc_pca/dendroplex_kienerii")
denkie.X <- tab(denkie, freq=TRUE, NA.method = "mean")
denkie.X <- scaleGen(denkie, NA.method = "mean")
denkie.X <- na.omit(denkie)
sum(is.na(denkie$tab)) #amount of missing data
denkie.X[1:2,1:2]
denkie.pca <- dudi.pca(denkie.X, scale=FALSE, scannf = FALSE, nf = 4)
denkie.pca$eig[1]/sum(denkie.pca$eig) * 100
denkie.pca$eig[2]/sum(denkie.pca$eig) * 100

denkie.pca$li

denkie.grp <- find.clusters(denkie, n.pca = 100)

xval <- xvalDapc(denkie.X, denkie.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(denkie.X, denkie.grp$grp)

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
tiff("dendroplex_kienerii_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
colorplot(denkie.pca$li, denkie.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-12, 12), ylim=c(-12, 12), 
          xlab="PC 1 (25.1%)", ylab="PC 2 (19.4%)",
          main = "PCA of Dendroplex kienerii 1306 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("dendroplex_kienerii_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(denkie.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             repel = TRUE,
             title = "PCA of Dendroplex kienerii 1306 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("dendroplex_kienerii_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(denkie.pca,
             xlim=c(-12, 12), ylim=c(-12, 12),
             geom = "point",
             title = "PCA of Dendroplex kienerii 1306 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("dendroplex_kienerii_black.tiff", width = 10, height = 10, units = 'in', 
     res = 800, compression = 'none')
fviz_pca_ind(denkie.pca,
             geom = "point",
             xlim=c(-12, 12), ylim=c(-12, 12),
             title = "PCA of Dendroplex kienerii 1306 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(denkie)
xval$DAPC$assign
xval[4]
xval[6]
# group assignment probabilities
round(xval$DAPC$posterior,6)
round(dapc1$posterior,6)


#summary(xval$DAPC)
#assignplot(xval$DAPC)
assignplot(dapc1)

contrib <- loadingplot(dapc1$var.contr, threshold = 0.0025, 
                       thres=.07, lab.jitter=1)
contrib$var.values

compoplot(xval$DAPC, show.lab = TRUE,  posi=list(x=12,y=-.01), cleg=.7)
#xval[2:6]
#xval$DAPC$assign
