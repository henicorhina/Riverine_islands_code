#!/bin/sh

source activate phyluce
# conirostrum_bicolor conirostrum_margaritae cranioleuca_vulpecula dendroplex_kienerii elaenia_pelzelni furnarius_minor furnarius_minor_full knipolegus_orenocensis leucippus_chlorocercus mazaria_propinqua myrmoborus_lugubris myrmochanes_hemileucus myrmotherula_assimilis myrmotherula_klagesi ochthornis_littoralis serpophaga_hypoleuca stigmatura_napensis stigmatura_napensis_full thamnophilus Campephilus_melanoleucos Campephilus_rubricollis Cantorchilus_leucotis Celeus_flavus Celeus_grammicus Crypturellus_undulatus Crypturellus_variegatus Formicarius_analis Formicarius_colma Glaucidium_brasilianum Glaucidium_hardyi Hylophylax_naevia Hylophylax_punctulata Megascops_choliba Megascops_watsonii Monasa_morphoeus Monasa_nigrifrons Myrmeciza_fortis Myrmeciza_hyperythra Myrmoborus_leucophrys Myrmoborus_myotherinus Phaethornis_bourcieri Phaethornis_hispidus Pheugopedius_coraya Piaya_cayana Piaya_melanogaster Pipra_erythrocephala Pipra_filicauda Saltator_coerulescens Saltator_grossus Schiffornis_major Schiffornis_turdina Synallaxis_gujanensis Synallaxis_rutilans Tachyphonus_cristatus Tachyphonus_luctuosus Trogon_rufus Xiphorhynchus_elegans Xiphorhynchus_obsoletus
for species in Trogon_collaris ; do 
	cd /Volumes/Brumfield_Lab_Drive/data/2_phasing/2_complete-taxon-set/1_cleaned_alignments-phylip-no_outgroups/${species}-phased-mafft-phylip-untrimmed-complete-clean/

	for file in *.phylip-relaxed;
	do
		raxmlHPC-PTHREADS-SSE3 \
		    -m GTRGAMMA \
		    -p 19877 \
		    -n ${file}_best.tre \
		    -s ${file} \
		    -T 8

	done;

	mv *.tre /Volumes/Brumfield_Lab_Drive/data/1_analysis/UCE_gene_trees/${species}/ ;
	mv *.reduced /Volumes/Brumfield_Lab_Drive/data/1_analysis/UCE_gene_trees/${species}/ ;

done

source deactivate
