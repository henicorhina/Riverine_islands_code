x <- read.table("/Volumes/Brumfield_Lab_Drive/data/1_analysis/heterozygosity_stats/allspecies.txt", header=TRUE)
head(x)
hist(x[ , c('F')])
F <- x[ , c('F')]
summary(F)
outliers <- mean(F) - (3*sd(F))
F_filtered <- F < outliers
