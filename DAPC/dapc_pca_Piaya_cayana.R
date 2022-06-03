library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Piaya_cayana
# Piaya cayana
# Piaya_cayana
# 1018 with your number of snps
# xlim=c(-15, 15), ylim=c(-15, 15), 

Piaya_cayana <- read.structure("Piaya_cayana_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=10, n.loc=1988, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Piaya_cayana_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Piaya_cayana_AmazonOnly")
Piaya_cayana.X <- tab(Piaya_cayana, freq=TRUE, NA.method = "mean")
Piaya_cayana.X <- scaleGen(Piaya_cayana, NA.method = "mean")
sum(is.na(Piaya_cayana$tab)) #amount of missing data
Piaya_cayana.X[1:2,1:2]
Piaya_cayana.pca <- dudi.pca(Piaya_cayana.X, scale=FALSE, scannf = FALSE, nf = 4)
Piaya_cayana.pca$eig[1]/sum(Piaya_cayana.pca$eig) * 100
Piaya_cayana.pca$eig[2]/sum(Piaya_cayana.pca$eig) * 100

Piaya_cayana.pca$li

Piaya_cayana.grp <- find.clusters(Piaya_cayana, n.pca = 100)

xval <- xvalDapc(Piaya_cayana.X, Piaya_cayana.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca

dapc1 <- dapc(Piaya_cayana.X, Piaya_cayana.grp$grp)

temp1 <- optim.a.score(dapc1)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score


#use these to save as high quality figures
tiff("Plot.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Plot .png", width = 6, height = 6, units = 'in', 
    res = 800)
bitmap("Plot.tiff", height = 6, width = 6, units = 'in', 
       type="tifflzw", res=800)
# run the desired format above then the plot, then this dev.off command
dev.off()
par(mfrow = c(1,1))


# pca with colored dots
tiff("Piaya_cayana_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Piaya_cayana_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Piaya_cayana.pca$li, Piaya_cayana.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-60, 25), ylim=c(-60, 25), 
          xlab="PC 1 (13.1%)", ylab="PC 2 (12.9%)",
          main = "PCA of Piaya cayana 1018 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Piaya_cayana_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Piaya_cayana_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Piaya_cayana.pca,
             xlim=c(-60, 25), ylim=c(-60, 25), 
             repel = TRUE,
             title = "PCA of Piaya cayana 1018 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Piaya_cayana_groups_forced_to_2.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Piaya_cayana_groups_forced_to_2.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Piaya_cayana.pca,
             xlim=c(-60, 25), ylim=c(-60, 25), 
             geom = "point",
             title = "PCA of Piaya cayana 1018 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Piaya_cayana_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Piaya_cayana_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Piaya_cayana.pca,
             geom = "point",
             xlim=c(-60, 25), ylim=c(-60, 25), 
             title = "PCA of Piaya cayana 1018 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Piaya_cayana)
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
