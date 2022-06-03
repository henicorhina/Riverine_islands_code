#!/bin/sh

for species in Trogon_rufus Pipra_erythrocephala Crypturellus_variegatus Cantorchilus_leucotis Glaucidium_brasilianum Glaucidium_hardyi Crypturellus_undulatus Celeus_grammicus Campephilus_melanoleucos Campephilus_rubricollis Celeus_flavus conirostrum_bicolor_full conirostrum_bicolor conirostrum_margaritae cranioleuca_vulpecula dendroplex_kienerii elaenia_pelzelni furnarius_minor furnarius_minor_full myrmochanes_hemileucus knipolegus_orenocensis leucippus_chlorocercus serpophaga_hypoleuca  mazaria_propinqua myrmoborus_lugubris stigmatura_napensis stigmatura_napensis_full thamnophilus_cryptoleucus  myrmotherula_assimilis myrmotherula_klagesi ochthornis_littoralis thamnophilus_nigrocinereus thamnophilus Trogon_collaris Tachyphonus_luctuosus Schiffornis_turdina Saltator_grossus Saltator_coerulescens Piaya_cayana Phaethornis_bourcieri Myrmoborus_myotherinus Formicarius_analis Hylophylax_naevia Pheugopedius_coraya Hylophylax_punctulata Pipra_filicauda Schiffornis_major Synallaxis_rutilans Formicarius_colma Monasa_morphoeus Monasa_nigrifrons Myrmeciza_fortis Myrmeciza_hyperythra Myrmoborus_leucophrys Phaethornis_hispidus Piaya_melanogaster Synallaxis_gujanensis Tachyphonus_cristatus Xiphorhynchus_elegans Xiphorhynchus_obsoletus Megascops_choliba Megascops_watsonii ; do 
 
	cd /Volumes/Brumfield_Lab_Drive/data/1_analysis/distruct_plots/${species}/

	for K in {1..6}; do	
		
		python /Users/mharvey/src/distruct2.3/distruct23/distruct2.3.py -K ${K} \
			--input=ClumppIndFile.output \
			--output=${species}_K${K}_by_individual.pdf \
			--title="${species} K=${K}" \
			--popfile=${species}_popfile_by_individual \
			--poporder=${species}_poporder_by_individual
		
		cp ${species}_K${K}_by_individual.pdf ../all_species/${species}_K${K}_by_individual.pdf
		
	done;
done

