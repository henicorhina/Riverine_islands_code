library("adegenet") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Campephilus_melanoleucos
# Campephilus melanoleucos
# Campephilus_melanoleucos
# 1764 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Campephilus_melanoleucos <- read.structure("Campephilus_melanoleucos_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=1764, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
#dir.create("Campephilus_melanoleucos_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Campephilus_melanoleucos_AmazonOnly")
Campephilus_melanoleucos.X <- tab(Campephilus_melanoleucos, freq=TRUE, NA.method = "mean")
Campephilus_melanoleucos.X <- scaleGen(Campephilus_melanoleucos, NA.method = "mean")
sum(is.na(Campephilus_melanoleucos$tab)) #amount of missing data
Campephilus_melanoleucos.X[1:2,1:2]
Campephilus_melanoleucos.pca <- dudi.pca(Campephilus_melanoleucos.X, scale=FALSE, scannf = FALSE, nf = 4)
Campephilus_melanoleucos.pca$eig[1]/sum(Campephilus_melanoleucos.pca$eig) * 100
Campephilus_melanoleucos.pca$eig[2]/sum(Campephilus_melanoleucos.pca$eig) * 100

Campephilus_melanoleucos.pca$li

Campephilus_melanoleucos.grp <- find.clusters(Campephilus_melanoleucos.X, n.pca = 100)

xval <- xvalDapc(Campephilus_melanoleucos.X, Campephilus_melanoleucos.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca
dapc1 <- dapc(Campephilus_melanoleucos.X, Campephilus_melanoleucos.grp$grp)

temp1 <- optim.a.score(dapc1)

#use these to save as high quality figures
# tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
#      res = 800, compression = 'none')
# png("Plot.png", width = 6, height = 6, units = 'in', 
#     res = 800)
# bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
#        type="tifflzw", res=800)
# # run the desired format above then the plot, then this dev.off command
# dev.off()
# par(mfrow = c(1,1))


# pca with colored dots
tiff("Campephilus_melanoleucos_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Campephilus_melanoleucos_colors.png", width = 6, height = 6, units = 'in',
    res = 800)
colorplot(Campephilus_melanoleucos.pca$li, Campephilus_melanoleucos.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-20, 75), ylim=c(-50, 20), 
          xlab="PC 1 (16.3%)", ylab="PC 2 (12.0%)",
          main = "PCA of Campephilus melanoleucos 1764 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Campephilus_melanoleucos_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Campephilus_melanoleucos_labeled.png", width = 12, height = 12, units = 'in',
    res = 800)
fviz_pca_ind(Campephilus_melanoleucos.pca,
             xlim=c(-20, 75), ylim=c(-50, 20), 
             repel = TRUE,
             title = "PCA of Campephilus melanoleucos 1764 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Campephilus_melanoleucos_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Campephilus_melanoleucos_groups.png", width = 6, height = 6, units = 'in',
    res = 800)
fviz_pca_ind(Campephilus_melanoleucos.pca,
             xlim=c(-20, 75), ylim=c(-50, 20), 
             geom = "point",
             title = "PCA of Campephilus melanoleucos 1764 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Campephilus_melanoleucos_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Campephilus_melanoleucos_black.png", width = 6, height = 6, units = 'in',
    res = 800)
fviz_pca_ind(Campephilus_melanoleucos.pca,
             geom = "point",
             xlim=c(-20, 75), ylim=c(-50, 20), 
             title = "PCA of Campephilus melanoleucos 1764 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Campephilus_melanoleucos)
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
