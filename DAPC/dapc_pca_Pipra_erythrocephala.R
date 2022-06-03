library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Pipra_erythrocephala
# Pipra_erythrocephala
# Pipra_erythrocephala
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Pipra_erythrocephala <- read.structure("Pipra_erythrocephala_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=11, n.loc=2008, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Pipra_erythrocephala_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Pipra_erythrocephala_AmazonOnly")
Pipra_erythrocephala.X <- tab(Pipra_erythrocephala, freq=TRUE, NA.method = "mean")
Pipra_erythrocephala.X <- scaleGen(Pipra_erythrocephala, NA.method = "mean")
sum(is.na(Pipra_erythrocephala$tab)) #amount of missing data
Pipra_erythrocephala.X[1:2,1:2]
Pipra_erythrocephala.pca <- dudi.pca(Pipra_erythrocephala.X, scale=FALSE, scannf = FALSE, nf = 4)
Pipra_erythrocephala.pca$eig[1]/sum(Pipra_erythrocephala.pca$eig) * 100
Pipra_erythrocephala.pca$eig[2]/sum(Pipra_erythrocephala.pca$eig) * 100

Pipra_erythrocephala.pca$li

Pipra_erythrocephala.grp <- find.clusters(Pipra_erythrocephala, n.pca = 100)

xval <- xvalDapc(Pipra_erythrocephala.X, Pipra_erythrocephala.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Pipra_erythrocephala.X, Pipra_erythrocephala.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
# tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
#      res = 800, compression = 'none')
# png("Plot .png", width = 6, height = 6, units = 'in', 
#     res = 800)
# bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
#        type="tifflzw", res=800)
# # run the desired format above then the plot, then this dev.off command
# dev.off()
# par(mfrow = c(1,1))


# pca with colored dots
tiff("Pipra_erythrocephala_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Pipra_erythrocephala_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Pipra_erythrocephala.pca$li, Pipra_erythrocephala.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-45, 20), ylim=c(-35, 35), 
          xlab="PC 1 (13.3%)", ylab="PC 2 (13.1%)",
          main = "PCA of Pipra_erythrocephala 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Pipra_erythrocephala_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Pipra_erythrocephala_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Pipra_erythrocephala.pca,
             xlim=c(-45, 20), ylim=c(-35, 35), 
             repel = TRUE,
             title = "PCA of Pipra_erythrocephala 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Pipra_erythrocephala_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Pipra_erythrocephala_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Pipra_erythrocephala.pca,
             xlim=c(-45, 20), ylim=c(-35, 35), 
             geom = "point",
             title = "PCA of Pipra_erythrocephala 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Pipra_erythrocephala_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Pipra_erythrocephala_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Pipra_erythrocephala.pca,
             geom = "point",
             xlim=c(-45, 20), ylim=c(-35, 35), 
             title = "PCA of Pipra_erythrocephala 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Pipra_erythrocephala)
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
