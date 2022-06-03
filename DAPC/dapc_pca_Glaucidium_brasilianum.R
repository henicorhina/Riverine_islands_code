library("adegenet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ade4", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("factoextra", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/structure_files_seqcap_pop_random_snp/")

# replace all these with info for your species
# Glaucidium_brasilianum
# Glaucidium brasilianum
# Glaucidium_brasilianum
# 1225 with your number of snps
# xlim=c(-12, 12), ylim=c(-12, 12), 

Glaucidium_brasilianum <- read.structure("Glaucidium_brasilianum_SNPs_phased_rmIndels_75_QC_DP_random_seqcap_pop_structure.AmazonOnly.str", 
                         n.ind=8, n.loc=1225, onerowperind=FALSE, 
                         col.lab=1, col.pop=0, col.others=0, row.marknames=0)

setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/")
dir.create("Glaucidium_brasilianum_AmazonOnly")
setwd("/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dapc_pca/Glaucidium_brasilianum_AmazonOnly")
Glaucidium_brasilianum.X <- tab(Glaucidium_brasilianum, freq=TRUE, NA.method = "mean")
Glaucidium_brasilianum.X <- scaleGen(Glaucidium_brasilianum, NA.method = "mean")
sum(is.na(Glaucidium_brasilianum$tab)) #amount of missing data
Glaucidium_brasilianum.X[1:2,1:2]
Glaucidium_brasilianum.pca <- dudi.pca(Glaucidium_brasilianum.X, scale=FALSE, scannf = FALSE, nf = 4)
Glaucidium_brasilianum.pca$eig[1]/sum(Glaucidium_brasilianum.pca$eig) * 100
Glaucidium_brasilianum.pca$eig[2]/sum(Glaucidium_brasilianum.pca$eig) * 100

Glaucidium_brasilianum.pca$li

Glaucidium_brasilianum.grp <- find.clusters(Glaucidium_brasilianum, n.pca = 100)

xval <- xvalDapc(Glaucidium_brasilianum.X, Glaucidium_brasilianum.grp$grp, n.pca.max = 20, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$DAPC$n.pca
dapc1 <- dapc(Glaucidium_brasilianum.X, Glaucidium_brasilianum.grp$grp)

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
tiff("Glaucidium_brasilianum_colors.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Glaucidium_brasilianum_colors.png", width = 6, height = 6, units = 'in', 
    res = 800)
colorplot(Glaucidium_brasilianum.pca$li, Glaucidium_brasilianum.pca$li, add.plot=FALSE, cex=2, 
          xlim=c(-15, 50), ylim=c(-50, 15),
          xlab="PC 1 (19.9%)", ylab="PC 2 (17.6)",
          main = "PCA of Glaucidium brasilianum 1225 SNPs \naxes 1-2") 
abline(v=0,h=0,col="grey", lty=2)
dev.off()


# pca with labeled points
tiff("Glaucidium_brasilianum_labeled.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Glaucidium_brasilianum_labeled.png", width = 12, height = 12, units = 'in', 
    res = 800)
fviz_pca_ind(Glaucidium_brasilianum.pca,
             xlim=c(-15, 50), ylim=c(-50, 15),
             repel = TRUE,
             title = "PCA of Glaucidium brasilianum 1225 SNPs \naxes 1-2"
)
dev.off()

# pca without labeled points - groups
tiff("Glaucidium_brasilianum_groups.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Glaucidium_brasilianum_groups.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Glaucidium_brasilianum.pca,
             xlim=c(-15, 50), ylim=c(-50, 15),
             geom = "point",
             title = "PCA of Glaucidium brasilianum 1225 SNPs \naxes 1-2"
) + geom_point(aes(colour = factor(dapc1$grp)), size = 3)
dev.off()


# pca without labeled points - black
tiff("Glaucidium_brasilianum_black.tiff", width = 6, height = 6, units = 'in', 
     res = 800, compression = 'none')
png("Glaucidium_brasilianum_black.png", width = 6, height = 6, units = 'in', 
    res = 800)
fviz_pca_ind(Glaucidium_brasilianum.pca,
             geom = "point",
             xlim=c(-15, 50), ylim=c(-50, 15),
             title = "PCA of Glaucidium brasilianum 1225 SNPs \naxes 1-2"
) + geom_point(size = 3)
dev.off()

# save these to a text file
indNames(Glaucidium_brasilianum)
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
