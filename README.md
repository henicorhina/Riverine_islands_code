# River Island code pipelines

Configuration files and code used in assembling UCE data for comparative river island phylogeography project. Most pipelines based on the Phyluce pipeline, with some modifications from seq_cap_pop

Remaining files are for data manipulation and PGLS / ANOVA analyses

-----------------

Folders:

`DAPC`: scripts used for Discriminant Analysis of Principle Components (DAPC) for each species

`seqcap_pop`: scripts used for calling SNPs using the seqcap_pop pipeline: https://github.com/mgharvey/seqcap_pop

-----------------

Scripts:

`PopGenome.R`: Fst and Dxy calculations with the R package PopGenome

`PopGenome.DAPC_assignments.R`: Fst and Dxy calculations with the R package PopGenome, but with samples assigned to genetic clusters from DAPC

`calculate_average_dxy_fst.py`: calculate Fst and Dxy metrics from alignments, averaging across all pairwise comparisons

`calculate_average_dxy_fst_all_species.sh`: run the `calculate_average_dxy_fst.py` script for all species

`gene_trees-all_species.sh`: estimate UCE gene trees in RAxML

`calculate_mtdna_tree_depth.R` and `calculate_uce_gene_tree_depth.R`: calculate average branch lengths of mitochondrial and UCE gene trees

`dendropy_popgenstats.py`: calculate population genetic metrics in DendroPy

`distruct_structure_all_species.sh`: run Distruct for visualizing STRUCTURE results

`genepop_ibd.R`: calculate isolation-by-distance slopes with the R package genepop and process results: `genepop_ibd_process_results.R`

`heterozygosity_outliers.R`: detect extreme heterozygosity outliers for removal from dataset

`pgls_caper.R`: run PGLS analysis in R package caper (old)

`phylANOVA.v4.R`: format input data and run phylogenetic ANOVA and PGLS analyses. Most analyses and plotting are included in this file.
