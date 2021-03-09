library("phytools")
library(ape)
library(nlme)
library(geiger)
library(caper)

tree_for_caper <- read.tree("/Volumes/Backup_Plus/data/all_species_tree/all_species_tree_pruned_rooted.tre")
setwd("/Users/home/Documents/Research/River_Islands/results/")

#for all data:
tree <- read.tree("/Volumes/Backup_Plus/data/all_species_tree/all_species_tree_pruned.tre")

rr.interactive<-reroot(tree,interactive=TRUE)
plotTree(rr.interactive)

all_for_caper = read.csv("results_for_pgls_take1_caper.csv", header = TRUE) #,row.names=1)
name.check(rr.interactive,all_for_caper)
head(all_for_caper)
comp.data<-comparative.data(rr.interactive, all_for_caper, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
plot(rr.interactive)
comp.data$dropped

#PHYLOGENY VARIABLES
all_for_caper_p = read.csv("results_for_pgls_take1_caper_phylogeny.csv", header = TRUE)
comp.data_p<-comparative.data(rr.interactive, all_for_caper_p, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)


modelo4<-pgls(structure_prob_cutoff~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(BAPS~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(dapc_laptop~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_groups~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(SNPs_per_bp~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_UCE_gene_tree_length~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(av_mtDNA_tree_depth_corr~all, data=comp.data)
modelo4<-pgls(trans_av_mtDNA_tree_depth_corr~all, data=comp.data)
modelo4<-pgls(trans2_av_mtDNA_tree_depth_corr~all, data=comp.data)
summary(modelo4)

modelo4<-pgls(D~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(theta_per_locus~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(seg_sites~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(nuc_div~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(heterozygosity~all, data=comp.data)
summary(modelo4)

modelo5<-pgls(subtending_branch~habitat_control, data=comp.data_p)
summary(modelo5)

modelo6<-pgls(stem_length~habitat_control, data=comp.data_p)
summary(modelo6)



#square root transformed data

modelo4<-pgls(trans_structure_prob_cutoff~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(trans_BAPS~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(trans_dapc_laptop~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(trans_Av_groups~habitat_control, data=comp.data)
summary(modelo4)

modelo4<-pgls(trans_HW_index~habitat_control, data=comp.data_kipp)
summary(modelo4)

all_for_caper$trans_HW_index <- sqrt(all_for_caper$HW_index)
all_for_caper_p$trans_stem_length <- sqrt((all_for_caper_p$stem_length)+1)
all_for_caper_dxy$trans_dxy_Oscar_seqcap_pop <- sqrt(all_for_caper_dxy$dxy_Oscar_seqcap_pop)

all_for_caper_kipp<-all_for_caper[complete.cases(all_for_caper[ , 31]),]
comp.data_kipp<-comparative.data(rr.interactive, all_for_caper_kipp, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
head(all_for_caper_kipp)
length(all_for_caper)
all_for_caper_kipp$trans_HW_index




#for vs floodplain

all_for_caper = read.csv("results_for_pgls_take1_caper_vs_floodplain.csv", header = TRUE)
name.check(rr.interactive,all_for_caper)
head(all_for_caper)
comp.data<-comparative.data(rr.interactive, all_for_caper, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
plot(rr.interactive)





#PHYLOGENY VARIABLES
all_for_caper_p = read.csv("results_for_pgls_take1_caper_phylogeny_vs_floodplain.csv", header = TRUE)
comp.data_p<-comparative.data(rr.interactive, all_for_caper_p, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)


modelo4<-pgls(structure_prob_cutoff~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_structure_prob_cutoff~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(BAPS~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_BAPS~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(dapc_laptop~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_dapc_laptop~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_groups~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_Av_groups~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(SNPs_per_bp~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_UCE_gene_tree_length~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(D~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(theta_per_locus~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(seg_sites~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(nuc_div~habitat, data=comp.data)
summary(modelo4)

modelo5<-pgls(subtending_branch~habitat_control, data=comp.data_p)
summary(modelo5)

modelo6<-pgls(stem_length~habitat_control, data=comp.data_p)
summary(modelo6)







#for vs upland

all_for_caper = read.csv("results_for_pgls_take1_caper_vs_upland.csv", header = TRUE)
name.check(rr.interactive,all_for_caper)
head(all_for_caper)
comp.data<-comparative.data(rr.interactive, all_for_caper, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
plot(rr.interactive)
comp.data$dropped

#PHYLOGENY VARIABLES
all_for_caper_p = read.csv("results_for_pgls_take1_caper_phylogeny_vs_upland.csv", header = TRUE)
comp.data_p<-comparative.data(rr.interactive, all_for_caper_p, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
comp.data_p$dropped
head(all_for_caper_p)


modelo4<-pgls(structure_prob_cutoff~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_structure_prob_cutoff~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(BAPS~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_BAPS~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(dapc_laptop~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_dapc_laptop~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_groups~habitat, data=comp.data)
summary(modelo4)
modelo4<-pgls(trans_Av_groups~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(SNPs_per_bp~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(Av_UCE_gene_tree_length~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(D~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(theta_per_locus~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(seg_sites~habitat, data=comp.data)
summary(modelo4)

modelo4<-pgls(nuc_div~habitat, data=comp.data)
summary(modelo4)

modelo5<-pgls(subtending_branch~habitat_control, data=comp.data_p)
summary(modelo5)

modelo6<-pgls(stem_length~habitat_control, data=comp.data_p)
summary(modelo6)





#KIPP'S INDEX

#for all data:
all_for_caper_kipp = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_kipps/results_for_pgls_take1_caper.csv", header = TRUE)
all_for_caper_kipp_f = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_kipps/results_for_pgls_take1_caper_vs_floodplain.csv", header = TRUE)
all_for_caper_kipp_u = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_kipps/results_for_pgls_take1_caper_vs_upland.csv", header = TRUE)

comp.data_kipp<-comparative.data(rr.interactive, all_for_caper_kipp, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
comp.data_kipp_f<-comparative.data(rr.interactive, all_for_caper_kipp_f, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
comp.data_kipp_u<-comparative.data(rr.interactive, all_for_caper_kipp_u, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)


name.check(rr.interactive,all_for_caper_kipp)
head(all_for_caper_kipp)
plot(rr.interactive)
#all_for_caper_f<-all_for_caper_kipp[!(all_for_caper_kipp$habitat=="u"),]
#all_for_caper_u<-all_for_caper[!(all_for_caper$habitat=="f"),]
#head(all_for_caper_f)

modelo4<-pgls(HW_index~habitat_control, data=comp.data_kipp)
summary(modelo4)

modelo4<-pgls(HW_index~habitat_control, data=comp.data_kipp_f)
summary(modelo4)

modelo4<-pgls(HW_index~habitat_control, data=comp.data_kipp_u)
summary(modelo4)


modelo4<-pgls(trans_HW_index~habitat_control, data=comp.data_kipp)
summary(modelo4)

modelo4<-pgls(trans_HW_index~habitat_control, data=comp.data_kipp_f)
summary(modelo4)

modelo4<-pgls(trans_HW_index~habitat_control, data=comp.data_kipp_u)
summary(modelo4)



#Dxy Fst

all_for_caper_dxy = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_DxyFst/results_for_pgls_take1_caper_both.csv", header = TRUE)
all_for_caper_dxy_f = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_DxyFst/results_for_pgls_take1_caper_f.csv", header = TRUE)
all_for_caper_dxy_u = read.csv("/Users/home/Documents/Research/River_Islands/results/results_spreadsheets_caper_DxyFst/results_for_pgls_take1_caper_u.csv", header = TRUE)

name.check(rr.interactive,all_for_caper_dxy)
head(all_for_caper_dxy)
comp.data_dxy<-comparative.data(rr.interactive, all_for_caper_dxy, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
plot(rr.interactive)
#all_for_caper_f<-all_for_caper_dxy[!(all_for_caper_dxy$habitat=="upland_forest"),]
#all_for_caper_u<-all_for_caper_dxy[!(all_for_caper_dxy$habitat=="floodplain_forest"),]
#all_for_caper_f$habitat_key <- all_for_caper_f$habitat[all_for_caper_f$habitat=="floodplain_forest"] <- 1


head(all_for_caper_f)
comp.data_dxy_f<-comparative.data(rr.interactive, all_for_caper_dxy_f, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)
comp.data_dxy_u<-comparative.data(rr.interactive, all_for_caper_dxy_u, names.col="species", vcv=TRUE, vcv.dim=2, na.omit=TRUE, warn.dropped=TRUE)

modelo4<-pgls(dxy_Oscar_seqcap_pop~habitat, data=comp.data_dxy)
summary(modelo4)

modelo4<-pgls(dxy_Oscar_seqcap_pop~habitat, data=comp.data_dxy_f)
summary(modelo4)

modelo4<-pgls(dxy_Oscar_seqcap_pop~habitat, data=comp.data_dxy_u)
summary(modelo4)

