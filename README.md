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



