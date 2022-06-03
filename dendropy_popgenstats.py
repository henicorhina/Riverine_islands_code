#!/usr/bin/env python


"""
Uses dendropy to calculate basic pop gen summary stats for each of some varzea/terra firme birds

Oscar Johnson
7 February 2022
"""
import pandas as pd
import dendropy
from dendropy.calculate import popgenstat

res = pd.DataFrame(columns = ['pairwise_differences', 'num_segregating_sites', 'wattersons_theta', 'tajimas_d', 'nucleotide_diversity'],
                   index = ['Campephilus_melanoleucos', 'Cantorchilus_leucotis', 'Formicarius_analis', 'Formicarius_colma', 'Glaucidium_brasilianum', 'Megascops_choliba', 'Monasa_nigrifrons', 'Pheugopedius_coraya', 'Piaya_cayana', 'Pipra_erythrocephala', 'Saltator_grossus', 'Schiffornis_turdina', 'Saltator_coerulescens', 'Tachyphonus_cristatus', 'Tachyphonus_luctuosus', 'Trogon_collaris', 'Trogon_rufus'])

# Campephilus_melanoleucos
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Campephilus_melanoleucos-phased-mafft-concatenated/Campephilus_melanoleucos-phased-mafft-concatenated.phylip", schema="phylip")
print("Campephilus_melanoleucos")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
	#.startswith is critical, because the phased alleles are stored as different individuals, with the same prefix
     if t.label.startswith('Campephilus_melanoleucos_LSUMNS2802'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Campephilus_melanoleucos'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
#save data after every species, in case something goes wrong
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')

#Cantorchilus_leucotis
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Cantorchilus_leucotis-phased-mafft-concatenated/Cantorchilus_leucotis-phased-mafft-concatenated.phylip", schema="phylip")
print("Cantorchilus_leucotis")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Cantorchilus_leucotis_LSUMNS10899'):
         p1.append(seqs[t])
     elif t.label.startswith('Cantorchilus_leucotis_LSUMNS46117'):
         p1.append(seqs[t])
     elif t.label.startswith('Cantorchilus_leucotis_LSUMNS9529'):
         p1.append(seqs[t])
     elif t.label.startswith('Cantorchilus_leucotis_LSUMNS7343'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Cantorchilus_leucotis'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Formicarius_analis
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Formicarius_analis-phased-mafft-concatenated/Formicarius_analis-phased-mafft-concatenated.phylip", schema="phylip")
print("Formicarius_analis")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Formicarius_analis_AMNH11914'):
         p1.append(seqs[t])
     elif t.label.startswith('Formicarius_analis_MPEG6767'):
         p1.append(seqs[t])
     elif t.label.startswith('Formicarius_analis_INPA1349'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Formicarius_analis'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')

#Formicarius_colma
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Formicarius_colma-phased-mafft-concatenated/Formicarius_colma-phased-mafft-concatenated.phylip", schema="phylip")
print("Formicarius_colma")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Formicarius_colma_MPEG9519'):
         p1.append(seqs[t])
     elif t.label.startswith('Formicarius_colma_MPEG12536'):
         p1.append(seqs[t])
     elif t.label.startswith('Formicarius_colma_LSUMNS13068'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Formicarius_colma'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')

#Glaucidium_brasilianum
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Glaucidium_brasilianum-phased-mafft-concatenated/Glaucidium_brasilianum-phased-mafft-concatenated.phylip", schema="phylip")
print("Glaucidium_brasilianum")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Glaucidium_brasilianum_AMNH2937'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Glaucidium_brasilianum'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')



#Megascops_choliba
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Megascops_choliba-phased-mafft-concatenated/Megascops_choliba-phased-mafft-concatenated.phylip", schema="phylip")
print("Megascops_choliba")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Megascops_choliba_LSUMNS42284'):
         p1.append(seqs[t])
     elif t.label.startswith('Megascops_choliba_MPEG13970'):
         p1.append(seqs[t])
     elif t.label.startswith('Megascops_choliba_AMNH4811'):
         p1.append(seqs[t])
     elif t.label.startswith('Megascops_choliba_LSUMNS7420'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Megascops_choliba'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Monasa_nigrifrons
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Monasa_nigrifrons-phased-mafft-concatenated/Monasa_nigrifrons-phased-mafft-concatenated.phylip", schema="phylip")
print("Monasa_nigrifrons")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Monasa_nigrifrons_MPEG10345'):
         p1.append(seqs[t])
     elif t.label.startswith('Monasa_nigrifrons_MPEG17327'):
         p1.append(seqs[t])
     elif t.label.startswith('Monasa_nigrifrons_MPEG2312'):
         p1.append(seqs[t])
     elif t.label.startswith('Monasa_nigrifrons_USNMB07074'):
         p1.append(seqs[t])
     elif t.label.startswith('Monasa_nigrifrons_MPEG11481'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Monasa_nigrifrons'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Pheugopedius_coraya
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Pheugopedius_coraya-phased-mafft-concatenated/Pheugopedius_coraya-phased-mafft-concatenated.phylip", schema="phylip")
print("Pheugopedius_coraya")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Pheugopedius_coraya_KU16858'):
         p1.append(seqs[t])
     elif t.label.startswith('Pheugopedius_coraya_MPEG12256'):
         p1.append(seqs[t])
     elif t.label.startswith('Pheugopedius_coraya_LSUMNS4133'):
         p1.append(seqs[t])
     elif t.label.startswith('Pheugopedius_coraya_AMNH4240'):
         p1.append(seqs[t])
     elif t.label.startswith('Pheugopedius_coraya_MPEG9534'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Pheugopedius_coraya'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Piaya_cayana
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Piaya_cayana-phased-mafft-concatenated/Piaya_cayana-phased-mafft-concatenated.phylip", schema="phylip")
print("Piaya_cayana")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Piaya_cayana_USNMB13940'):
         p1.append(seqs[t])
     elif t.label.startswith('Piaya_cayana_LSUMNS44417'):
         p1.append(seqs[t])
     elif t.label.startswith('Piaya_cayana_KU16728'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Piaya_cayana'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Pipra_erythrocephala
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Pipra_erythrocephala-phased-mafft-concatenated/Pipra_erythrocephala-phased-mafft-concatenated.phylip", schema="phylip")
print("Pipra_erythrocephala")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Pipra_chloromeros_KU18533'):
         p1.append(seqs[t])
     elif t.label.startswith('Pipra_chloromeros_LSUMNS106768'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Pipra_erythrocephala'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')



#Saltator_grossus
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Saltator_grossus-phased-mafft-concatenated/Saltator_grossus-phased-mafft-concatenated.phylip", schema="phylip")
print("Saltator_grossus")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Saltator_grossus_AMNH12692'):
         p1.append(seqs[t])
     elif t.label.startswith('Saltator_grossus_MPEG7434'):
         p1.append(seqs[t])
     elif t.label.startswith('Saltator_grossus_USNMB09368'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Saltator_grossus'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Schiffornis_turdina
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Schiffornis_turdina-phased-mafft-concatenated/Schiffornis_turdina-phased-mafft-concatenated.phylip", schema="phylip")
print("Schiffornis_turdina")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Schiffornis_turdina_LSUMNS27903'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Schiffornis_turdina'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Saltator_coerulescens
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Saltator_coerulescens-phased-mafft-concatenated/Saltator_coerulescens-phased-mafft-concatenated.phylip", schema="phylip")
print("Saltator_coerulescens")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Saltator_coerulescens_INPA2123'):
         p1.append(seqs[t])
     elif t.label.startswith('Saltator_coerulescens_USNMB11140'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Saltator_coerulescens'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')

#Tachyphonus_cristatus
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Tachyphonus_cristatus-phased-mafft-concatenated/Tachyphonus_cristatus-phased-mafft-concatenated.phylip", schema="phylip")
print("Tachyphonus_cristatus")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Tachyphonus_cristatus_LSUMNS2693'):
         p1.append(seqs[t])
     elif t.label.startswith('Tachyphonus_cristatus_AMNH2985'):
         p1.append(seqs[t])
     elif t.label.startswith('Tachyphonus_cristatus_INPA1636'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Tachyphonus_cristatus'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')



#Tachyphonus_luctuosus
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Tachyphonus_luctuosus-phased-mafft-concatenated/Tachyphonus_luctuosus-phased-mafft-concatenated.phylip", schema="phylip")
print("Tachyphonus_luctuosus")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Tachyphonus_luctuosus_MVZ169546'):
         p1.append(seqs[t])
     elif t.label.startswith('Tachyphonus_luctuosus_MPEG4450'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Tachyphonus_luctuosus'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Trogon_collaris
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Trogon_collaris-phased-mafft-concatenated/Trogon_collaris-phased-mafft-concatenated.phylip", schema="phylip")
print("Trogon_collaris")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Trogon_collaris_MPEG8909'):
         p1.append(seqs[t])
     elif t.label.startswith('Trogon_collaris_USNMB22128'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Trogon_collaris'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')


#Trogon_rufus
seqs = dendropy.DnaCharacterMatrix.get(path="/Volumes/Brumfield_Lab_Drive/River_islands/2_phasing/4_alignments/Trogon_rufus-phased-mafft-concatenated/Trogon_rufus-phased-mafft-concatenated.phylip", schema="phylip")
print("Trogon_rufus")
# calculate wattersons theta and tajimas d for all groups
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
     if t.label.startswith('Trogon_rufus_FMNH456559'):
         p1.append(seqs[t])
     elif t.label.startswith('Trogon_rufus_INPA1668'):
         p1.append(seqs[t])
     elif t.label.startswith('Trogon_rufus_INPA1170'):
         p1.append(seqs[t])
     elif t.label.startswith('Trogon_rufus_LSUMNS4256'):
         p1.append(seqs[t])
     elif t.label.startswith('Trogon_rufus_LSUMNS27570'):
         p1.append(seqs[t])
     else:
         p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
nuc = popgenstat.nucleotide_diversity(seqs)
res.loc['Trogon_rufus'] = [pp.average_number_of_pairwise_differences, pp.num_segregating_sites, pp.wattersons_theta, pp.tajimas_d, nuc]
res.to_csv(path_or_buf='/Volumes/Brumfield_Lab_Drive/River_islands/3_results/dendropy/AmazonOnly_reanalysis_dendropy.csv')
