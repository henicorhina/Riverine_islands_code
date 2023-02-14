#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(caper)


setwd("/Volumes/Brumfield_Lab_Drive/River_islands")


key <- read.csv("Habitat_key_all_samples.csv")
ri_inds <- read.csv("all_sample_individuals_NoOutgroups.csv")
het <- read.csv("1_analysis/heterozygosity_stats/allspecies_heterozygosity.csv")
island.Tree <- read.tree(file = "all_species_tree/incomplete_taxon_set/mafft-raxml-nexus-edge-trimmed-95percent/RAxML_bipartitions.all_species_95percent_final.newick.phy")

ri_inds$keep <- "Y"
ri_inds <- ri_inds[,c(1,3)]
island.Tree <- drop.tip(island.Tree, "Thamnophilus_cryptoleucus2")

key <- dplyr::left_join(key, ri_inds, by = c("species_sample" = "sample" ))

key$keep <- if_else(key$habitat != "island", "Y", key$keep)

key <- key %>% 
  dplyr::filter(non_amazonian != "Y") %>%
  dplyr::filter(keep == "Y")
key <- key[,1:3]

het.t <- dplyr::left_join(key, het, by = c("species_sample" = "INDV"))
colnames(het.t)[7] <- "inbreeding_coef"

# remove extreme outliers that may have sequencing issues
het.t <- het.t %>%
  dplyr::filter(inbreeding_coef > -0.5)

write.csv(het.t, file = "1_analysis/heterozygosity_stats/allspecies_heterozygosity.formatted.v2.csv", row.names = FALSE)
het.t <- read.csv("1_analysis/heterozygosity_stats/allspecies_heterozygosity.formatted.v2.data_sources_added.csv")



ggplot(het.t, aes(x=habitat, y=inbreeding_coef)) + 
  geom_boxplot(aes(fill=habitat)) + 
  geom_point(aes(group=habitat)) +
  # geom_point(position=position_dodge(width=0.75),aes(group=habitat)) +
  scale_x_discrete(limits=c("island", "floodplain", "upland")) + 
  ggtitle("inbreeding coefficient") +
  labs(y = "inbreeding coefficient")



# make a copy for t-test
# merge floodplain and upland
het.t2 <- het.t

# bin by exact data source for each sample
# peer with new
het.t2$data_source.2 <- if_else(het.t2$data_source != "old", "new", het.t2$data_source)
t.test(dplyr::filter(het.t2, data_source.2 == "new")$inbreeding_coef, dplyr::filter(het.t2, data_source.2 == "old")$inbreeding_coef)

# re-position peer samples to lump with old
het.t2$data_source.2 <- if_else(het.t2$data_source != "new", "old", het.t2$data_source) 
t.test(dplyr::filter(het.t2, data_source.2 == "new")$inbreeding_coef, dplyr::filter(het.t2, data_source.2 == "old")$inbreeding_coef)

ggplot(het.t2, aes(x=data_source.2, y=inbreeding_coef)) + 
  geom_boxplot() + 
  geom_point() +
  ggtitle("inbreeding coefficient") +
  labs(y = "inbreeding coefficient")


# bin by exact data source for each sample
# lump peer samples with old
# and plot heterozygosity
het.t2$data_source.2 <- if_else(het.t2$data_source != "new", "old", het.t2$data_source) 
het.t2 <- het.t2 %>% dplyr::mutate(O.HET = (N_SITES - O.HOM.) / N_SITES) 
t.test(dplyr::filter(het.t2, data_source.2 == "new")$O.HET, dplyr::filter(het.t2, data_source.2 == "old")$O.HET)


# t = 10.123, df = 114.44, p-value < 2.2e-16
# this is the one to use


ggplot(het.t2, aes(x=data_source.2, y=O.HET)) + 
  geom_boxplot() + 
  geom_point() +
  ggtitle("observed heterozygosity") +
  labs(y = "heterozygosity")

ggplot(het.t2, aes(x=habitat, y=O.HET)) + 
  geom_boxplot(aes(fill=habitat)) + 
  geom_point(aes(group=habitat)) +
  # geom_point(position=position_dodge(width=0.75),aes(group=habitat)) +
  scale_x_discrete(limits=c("island", "floodplain", "upland")) + 
  ggtitle("observed heterozygosity") +
  labs(y = "heterozygosity")


# plot only old data sources, and separate island from the rest
# so, do island actually have higher heterozygosity? for a given data source
het.t3 <- het.t2 %>% dplyr::filter(data_source.2 == "old")
het.t3$habitat.binned <- if_else(het.t3$habitat != "island", "floodplain_upland", het.t3$habitat) 
t.test(dplyr::filter(het.t3, habitat.binned == "island")$O.HET, dplyr::filter(het.t3, habitat.binned == "floodplain_upland")$O.HET)

# still significant
# t = 7.0378, df = 54.573, p-value = 3.397e-09

ggplot(het.t3, aes(x=habitat, y=O.HET)) + 
  geom_boxplot(aes(fill=habitat)) + 
  geom_point(aes(group=habitat)) +
  scale_x_discrete(limits=c("island", "floodplain", "upland")) + 
  ggtitle("observed heterozygosity") +
  labs(y = "heterozygosity")


# this was binning by habitat
# het.t2$habitat <- if_else(het.t2$habitat != "island", "floodplain_upland", het.t2$habitat)
# t.test(dplyr::filter(het.t2, habitat == "island")$inbreeding_coef, dplyr::filter(het.t2, habitat == "floodplain_upland")$inbreeding_coef)

one.way <- aov(inbreeding_coef ~ habitat, data = het.t)
summary(one.way)



sp_het <- het.t %>%
  group_by(species) %>%
  summarise(mean = mean(inbreeding_coef), n = n())
sp_het <- sp_het %>%
  left_join(sp_het, het.t, by = "species")

pgls.char <- function(df.t, char) {
  # df.t = dataframe
  # char = column number of the character in question
  
  
  colnames(df.t)[char] <- "char"
  temp <- df.t %>%
    dplyr::select(habitat, species, char) %>% 
    dplyr::filter(!is.na(char)) %>%
    droplevels()
  
  temp$habitat <- gsub("island", "0", temp$habitat)
  temp$habitat <- gsub("floodplain", "1", temp$habitat)
  temp$habitat <- gsub("upland", "2", temp$habitat)
  temp$habitat <- as.numeric(temp$habitat)
  
  island.Tree.t <- keep.tip(island.Tree, as.character(temp$species))
  
  hab.dat <- comparative.data(island.Tree.t, temp, species)
  
  if (char == 25) { # absolute value for log-transform Tajima's D. all are negative
    res.pgls <- pgls(log(abs(char)) ~ habitat, hab.dat)
  } else {
    # res.pgls <- pgls(log(char) ~ habitat, hab.dat)
    res.pgls <- pgls(habitat ~ log(char), hab.dat)
  }
  # summary(res.pgls)
  # coef(res.pgls)
  
  
  # plot(log(temp$char) ~ temp$habitat)
  # abline(a = coef(res.pgls)[1], b = coef(res.pgls)[2])
  
  
  res_list <- list("pgls" = res.pgls, 
                   "estimate" = summary(res.pgls)[5]$coefficients[2,1], 
                   "p" = summary(res.pgls)[5]$coefficients[2,4], 
                   "t" = summary(res.pgls)[5]$coefficients[2,3], 
                   "aicc" = res.pgls$aicc,
                   "slope" = coef(res.pgls)[2],
                   "int" = coef(res.pgls)[1],
                   "sterr" = res.pgls$sterr)
  
  
  return(res_list)
}

# summarize across species for species-average inbreeding coefficient values
# and (N_SITES - O.HOM.) / N_SITES (per-site heterozygosity)

het_observed.nsites <- het.t %>% group_by(species) %>% dplyr::summarise(heterozygosity.nsites = mean(N_SITES - O.HOM.))
het_observed <- het.t %>% group_by(species) %>% dplyr::summarise(heterozygosity = mean((N_SITES - O.HOM.) / N_SITES))
inbreeding <- het.t %>% group_by(species) %>% dplyr::summarise(inbreeding = mean(inbreeding_coef))
N_SITES <- het.t %>% group_by(species) %>% dplyr::summarise(N_SITES = mean(N_SITES))

het.all <- left_join(het_observed.nsites, inbreeding, by = "species")
het.all <- left_join(het_observed, het.all, by = "species")
het.all <- left_join(het.all, N_SITES, by = "species")
write.csv(het.all, file = "1_analysis/heterozygosity_stats/allspecies_heterozygosity.formatted.csv", row.names = FALSE)

df.het <- read.csv("3_results/trait.database.formatted.final.csv")
df.het <- left_join(df.het, het_observed, by = "species")
df.het <- left_join(df.het, inbreeding, by = "species")
df.het <- left_join(df.het, N_SITES, by = "species")

df.het <- df.het %>% dplyr::select(-structure_prob_cutoff, -Dxy_Mike_and_Oscar, -Fst_Mike_and_Oscar, -log.total_bp, -log.theta)
df.het <- df.het[,c(1:2,23,3:22,24:59)]
write.csv(df.het, file = "4_River_Island_tables_figures/Supplemental_data_table.csv", row.names = FALSE)

df.t <- df.het %>% filter(species != "Elaenia_pelzelni")

res.1 <- pgls.char(df.t, char=57)
res.1
res.2 <- pgls.char(df.t, char=58)
res.2

df.t %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = heterozygosity)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nHeterozygosity") 

ggsave("3_results/1_plots/heterozygosity.v1.pdf", width = 6,
       height = 6,
       units = c("in"))

df.t %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = inbreeding)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nInbreeding coefficient") 

ggsave("3_results/1_plots/inbreeding.v1.pdf", width = 6,
       height = 6,
       units = c("in"))


N_SITES.p <- df.t %>%
  mutate(habitat = fct_relevel(habitat, "island", "floodplain", "upland")) %>%
  ggplot( aes(x = habitat, y = log10(N_SITES))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  theme(plot.margin = unit(c(10,10,10,10), "pt"), legend.position = "none") +
  xlab("\n") +
  ylab("\n\nNumber of sites (log)") 


