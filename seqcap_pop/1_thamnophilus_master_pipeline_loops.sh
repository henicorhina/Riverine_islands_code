#!/bin/sh

cd /Volumes/Brumfield_Lab_Drive/data/by_species_take2/thamnophilus/

mkdir 5_mapping 5_mapping/bam 5_mapping/sam 6_picard 6_picard/bam_cleaned 6_picard/mark_pcr_duplicates 
mkdir 6_picard/mark_pcr_metrics 6_picard/read_groups 7_merge-bams 8_GATK
mkdir 9_SNP-tables 10_sequences 11_fasta-parts 12_raw-alignments 13_processed-phylip 14_formatted_output
mkdir 14_formatted_output/adegenet 14_formatted_output/gphocs 14_formatted_output/structure 14_formatted_output/faststructure 14_formatted_output/genepop
cd /Volumes/Brumfield_Lab_Drive/data/by_species_take2/

# bwa index

bwa index -a is /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta

for SAMPLE in thamnophilus_cryptoleucus_139_LSU25431 thamnophilus_cryptoleucus_141_LSU74103 thamnophilus_cryptoleucus_142_LSU7285 thamnophilus_cryptoleucus_144_LSU93318 thamnophilus_cryptoleucus_145_MPEG_T22677 thamnophilus_cryptoleucus_146_MPEG_T23595 thamnophilus_cryptoleucus_147_MPEG_T22714 thamnophilus_nigrocinereus_151_INPA_A10843 thamnophilus_nigrocinereus_154_INPA_A2165 thamnophilus_nigrocinereus_155_INPA_A3093 thamnophilus_nigrocinereus_158_INPA_A295 thamnophilus_nigrocinereus_165_MPEG_T77261 thamnophilus_nigrocinereus_173_ICN39301 thamnophilus_nigrocinereus_174_LSU20234 thamnophilus_nigrocinereus_175_MZUSP92841 thamnophilus_nigrocinereus_176_MZUSP93302 thamnophilus_nigrocinereus_177_MZUSP_OI001 thamnophilus_cryptoleucus_244_MPEG_T22585
do

	# step 4: map reads to contigs

	/Users/mharvey/src/bwa-0.7.17/bwa mem /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
	thamnophilus/2_clean_reads/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ1.fastq.gz \
	thamnophilus/2_clean_reads/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ2.fastq.gz \
	> thamnophilus/5_mapping/sam/${SAMPLE}-aln-pe.sam

	# step 5: convert sam to bam
	
	samtools view -bS thamnophilus/5_mapping/sam/${SAMPLE}-aln-pe.sam \
	> thamnophilus/5_mapping/bam/${SAMPLE}-aln-pe.bam

	# step 6: clean bam

	java -jar ~/anaconda/jar/CleanSam.jar \
	I=thamnophilus/5_mapping/bam/${SAMPLE}-aln-pe.bam \
	O=thamnophilus/6_picard/bam_cleaned/${SAMPLE}-aln-pe_CL.bam \
	VALIDATION_STRINGENCY=SILENT


	# step 7: add read groups

	java -Xmx2g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar \
	I=thamnophilus/6_picard/bam_cleaned/${SAMPLE}-aln-pe_CL.bam \
	O=thamnophilus/6_picard/read_groups/${SAMPLE}-aln_RG.bam \
	SORT_ORDER=coordinate \
	RGPL=illumina \
	RGPU=TestXX \
	RGLB=Lib1 \
	RGID=${SAMPLE} \
	RGSM=${SAMPLE} \
	VALIDATION_STRINGENCY=LENIENT

	# step 8: mark PCR duplicate reads

	java -Xmx2g -jar ~/anaconda/jar/MarkDuplicates.jar \
	I=thamnophilus/6_picard/read_groups/${SAMPLE}-aln_RG.bam \
	O=thamnophilus/6_picard/mark_pcr_duplicates/${SAMPLE}-aln_MD.bam \
	METRICS_FILE=thamnophilus/6_picard/mark_pcr_metrics/${SAMPLE}.metrics \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=false

done

# step 9: merge bams

java -Xmx2g -jar ~/anaconda/jar/MergeSamFiles.jar \
SO=coordinate \
AS=true \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_139_LSU25431-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_141_LSU74103-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_142_LSU7285-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_144_LSU93318-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_145_MPEG_T22677-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_146_MPEG_T23595-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_147_MPEG_T22714-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_cryptoleucus_244_MPEG_T22585-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_151_INPA_A10843-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_154_INPA_A2165-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_155_INPA_A3093-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_158_INPA_A295-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_165_MPEG_T77261-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_173_ICN39301-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_174_LSU20234-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_175_MZUSP92841-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_176_MZUSP93302-aln_MD.bam \
I=thamnophilus/6_picard/mark_pcr_duplicates/thamnophilus_nigrocinereus_177_MZUSP_OI001-aln_MD.bam \
O=thamnophilus/7_merge-bams/thamnophilus.bam 


# step 10: index bams


samtools index thamnophilus/7_merge-bams/thamnophilus.bam 


# steps 11 and 12: create and index dictionary

java -Xmx2g -jar ~/anaconda/pkgs/picard-1.106-0/jar/CreateSequenceDictionary.jar \
R=/Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
O=/Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.dict

samtools faidx /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta


# step 13: call indels

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/7_merge-bams/thamnophilus.bam   \
--minReadsAtLocus 7 \
-o thamnophilus/8_GATK/thamnophilus.intervals



# step 14: realign indels

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/7_merge-bams/thamnophilus.bam  \
-targetIntervals thamnophilus/8_GATK/thamnophilus.intervals \
-LOD 3.0 \
-o thamnophilus/8_GATK/thamnophilus_RI.bam


# step 15: call SNPs

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/8_GATK/thamnophilus_RI.bam \
-gt_mode DISCOVERY \
-o thamnophilus/8_GATK/thamnophilus_raw_SNPs.vcf \
-ploidy 2 \
-rf BadCigar


# step 16: annotate SNPs

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/8_GATK/thamnophilus_RI.bam \
-G StandardAnnotation \
-V:variant,VCF thamnophilus/8_GATK/thamnophilus_raw_SNPs.vcf \
-XA SnpEff \
-o thamnophilus/8_GATK/thamnophilus_SNPs_annotated.vcf \
-rf BadCigar      



# step 17: annotate indels

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/8_GATK/thamnophilus_RI.bam \
-gt_mode DISCOVERY \
-glm INDEL \
-o thamnophilus/8_GATK/thamnophilus_SNPs_indels.vcf \
-rf BadCigar         



# step 18: mask indels

java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-V thamnophilus/8_GATK/thamnophilus_raw_SNPs.vcf \
--mask thamnophilus/8_GATK/thamnophilus_SNPs_indels.vcf \
--maskExtension 5 \
--maskName InDel \
--clusterWindowSize 10 \
--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filterName "BadValidation" \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
--filterExpression "QD < 5.0" \
--filterName "LowVQCBD" \
-o thamnophilus/8_GATK/thamnophilus_SNPs_no_indels.vcf  \
-rf BadCigar


# step 19: Restrict to high-quality SNPs (bash)

cat thamnophilus/8_GATK/thamnophilus_SNPs_no_indels.vcf | grep 'PASS\|^#' > thamnophilus/8_GATK/thamnophilus_SNPs_pass-only.vcf 


# step 20: Read-backed phasing (GATK)


java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
-I thamnophilus/8_GATK/thamnophilus_RI.bam \
--variant thamnophilus/8_GATK/thamnophilus_SNPs_pass-only.vcf \
-L thamnophilus/8_GATK/thamnophilus_SNPs_pass-only.vcf \
-o thamnophilus/8_GATK/thamnophilus_SNPs_phased.vcf \
--phaseQualityThresh 20.0 \
-rf BadCigar


for SAMPLE in thamnophilus_cryptoleucus_139_LSU25431 thamnophilus_cryptoleucus_141_LSU74103 thamnophilus_cryptoleucus_142_LSU7285 thamnophilus_cryptoleucus_144_LSU93318 thamnophilus_cryptoleucus_145_MPEG_T22677 thamnophilus_cryptoleucus_146_MPEG_T23595 thamnophilus_cryptoleucus_147_MPEG_T22714 thamnophilus_nigrocinereus_151_INPA_A10843 thamnophilus_nigrocinereus_154_INPA_A2165 thamnophilus_nigrocinereus_155_INPA_A3093 thamnophilus_nigrocinereus_158_INPA_A295 thamnophilus_nigrocinereus_165_MPEG_T77261 thamnophilus_nigrocinereus_173_ICN39301 thamnophilus_nigrocinereus_174_LSU20234 thamnophilus_nigrocinereus_175_MZUSP92841 thamnophilus_nigrocinereus_176_MZUSP93302 thamnophilus_nigrocinereus_177_MZUSP_OI001 thamnophilus_cryptoleucus_244_MPEG_T22585
do

# steps 21-23

	java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
	--variant thamnophilus/8_GATK/thamnophilus_SNPs_phased.vcf \
	-o thamnophilus/8_GATK/${SAMPLE}_SNPs.vcf \
	-sn ${SAMPLE} \
	-rf BadCigar


	java -Xmx2g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
	-T VariantsToTable \
	-R /Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
	-V thamnophilus/8_GATK/${SAMPLE}_SNPs.vcf \
	-F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
	-o thamnophilus/9_SNP-tables/${SAMPLE}_SNPs_phased-table.txt \
	-rf BadCigar


	python /Users/mharvey/seqcap_pop-master/bin/add_phased_snps_to_seqs_filter.py \
	/Volumes/Brumfield_Lab_Drive/data/2_itero/mafft-nexus-edge-trimmed-exploded-uncontaminated/thamnophilus_cryptoleucus_144_LSU93318.fasta \
	thamnophilus/9_SNP-tables/${SAMPLE}_SNPs_phased-table.txt \
	thamnophilus/10_sequences/${SAMPLE}_sequences.txt \
	1


done

: '


# step 24 

python /Users/mharvey/seqcap_pop-master/bin/collate_sample_fastas_GATK.py \
	thamnophilus/10_sequences/ \
	thamnophilus/11_fasta-parts/ \
	sequences.txt

# step 25: Align the sequences

python /Users/mharvey/seqcap_pop-master/bin/run_mafft.py \
	thamnophilus/11_fasta-parts/ \
	thamnophilus/12_raw-alignments/


# step 26: Process the alignments

python /Users/mharvey/seqcap_pop-master/bin/process_mafft_alignments_GATK.py \
	thamnophilus/12_raw-alignments/ \
	thamnophilus/13_processed-phylip/
	
# gphocs

python /Users/mharvey/seqcap_pop-master/bin/gphocs_from_phy.py \
	thamnophilus/13_processed-phylip/ \
	thamnophilus/14_formatted_output/gphocs/thamnophilus_GPHOCS.txt


# structure

python /Users/mharvey/seqcap_pop-master/bin/structure_from_vcf.py \
	thamnophilus/8_GATK/conirostrum_margaritae_SNPs_phased.vcf \
	thamnophilus/14_formatted_output/structure/conirostrum_margaritae_STRUCTURE.txt
	
	
# faststructure

python /Users/mharvey/seqcap_pop-master/bin/faststructure_from_vcf.py \
	thamnophilus/8_GATK/conirostrum_margaritae_SNPs_phased.vcf \
	thamnophilus/14_formatted_output/faststructure/ \
	thamnophilus \
	--random


# adegenet

python /Users/mharvey/seqcap_pop-master/bin/adegenet_from_vcf.py \
	thamnophilus/8_GATK/conirostrum_margaritae_SNPs_phased.vcf \
	thamnophilus/14_formatted_output/adegenet/ \
	thamnophilus \
	--random


# genepop

python /Users/mharvey/seqcap_pop-master/bin/genepop_from_vcf.py \
	thamnophilus/8_GATK/conirostrum_margaritae_SNPs_phased.vcf \
	thamnophilus/14_formatted_output/genepop/ \
	thamnophilus \
	--random

'

