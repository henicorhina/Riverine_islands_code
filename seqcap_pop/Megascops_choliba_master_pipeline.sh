#!/bin/sh

cd /Volumes/Brumfield_Lab_Drive/mike_data/6_seqcap_pop/Megascops_choliba/

mkdir 5_mapping 5_mapping/bam 5_mapping/sam 6_picard 6_picard/bam_cleaned 6_picard/mark_pcr_duplicates 
mkdir 6_picard/mark_pcr_metrics 6_picard/read_groups 7_merge-bams 8_GATK
mkdir 9_SNP-tables 10_sequences 11_fasta-parts 12_raw-alignments 13_processed-phylip 14_formatted_output
mkdir 14_formatted_output/adegenet 14_formatted_output/gphocs 14_formatted_output/structure 14_formatted_output/faststructure 14_formatted_output/genepop


# in this document you will need to replace the following:
# sample names in steps 4-8, step 9, and steps 21-23
# reference individual "Megascops_choliba_MPEG16711.fasta"
# all instances of "Megascops_choliba" with the name of your taxon

# paths to GATK, picard, and seqcap_pop scripts
# the maximum java ram permitted on your system (e.g. 'java -Xmx16g')

# bwa index

bwa index -a is /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta

# steps 4-8

for SAMPLE in Megascops_choliba_AMNH4811 Megascops_choliba_FMNH392673 Megascops_choliba_LSUMNS7420 Megascops_choliba_LSUMNS9593 Megascops_choliba_LSUMNS15318 Megascops_choliba_LSUMNS38217 Megascops_choliba_LSUMNS42284 Megascops_choliba_MPEG6986 Megascops_choliba_MPEG13970 Megascops_choliba_MPEG16509 Megascops_choliba_MPEG16711 Megascops_choliba_MPEG17491 Megascops_choliba_UWBM116721 
do
	# step 4: map reads to contigs

	/Users/mharvey/src/bwa-0.7.17/bwa mem /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
	/Volumes/Brumfield_Lab_Drive/mike_data/2_clean_reads/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ1.fastq.gz \
	/Volumes/Brumfield_Lab_Drive/mike_data/2_clean_reads/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ2.fastq.gz \
	> 5_mapping/sam/${SAMPLE}-aln-pe.sam

	# step 5: convert sam to bam
	
	samtools view -bS 5_mapping/sam/${SAMPLE}-aln-pe.sam \
	> 5_mapping/bam/${SAMPLE}-aln-pe.bam

	# step 6: clean bam

	java -jar ~/anaconda/jar/CleanSam.jar \
	I=5_mapping/bam/${SAMPLE}-aln-pe.bam \
	O=6_picard/bam_cleaned/${SAMPLE}-aln-pe_CL.bam \
	VALIDATION_STRINGENCY=SILENT


	# step 7: add read groups

	java -Xmx16g -jar ~/anaconda/jar/AddOrReplaceReadGroups.jar \
	I=6_picard/bam_cleaned/${SAMPLE}-aln-pe_CL.bam \
	O=6_picard/read_groups/${SAMPLE}-aln_RG.bam \
	SORT_ORDER=coordinate \
	RGPL=illumina \
	RGPU=TestXX \
	RGLB=Lib1 \
	RGID=${SAMPLE} \
	RGSM=${SAMPLE} \
	VALIDATION_STRINGENCY=LENIENT

	# step 8: mark PCR duplicate reads

	java -Xmx16g -jar ~/anaconda/jar/MarkDuplicates.jar \
	I=6_picard/read_groups/${SAMPLE}-aln_RG.bam \
	O=6_picard/mark_pcr_duplicates/${SAMPLE}-aln_MD.bam \
	METRICS_FILE=6_picard/mark_pcr_metrics/${SAMPLE}.metrics \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=false

done

# step 9: merge bams

java -Xmx16g -jar ~/anaconda/jar/MergeSamFiles.jar \
SO=coordinate \
AS=true \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_AMNH4811-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_FMNH392673-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_LSUMNS7420-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_LSUMNS9593-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_LSUMNS15318-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_LSUMNS38217-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_LSUMNS42284-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_MPEG6986-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_MPEG13970-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_MPEG16509-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_MPEG16711-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_MPEG17491-aln_MD.bam \
I=6_picard/mark_pcr_duplicates/Megascops_choliba_UWBM116721-aln_MD.bam \
O=7_merge-bams/Megascops_choliba.bam 


# step 10: index bams


samtools index 7_merge-bams/Megascops_choliba.bam 


# steps 11 and 12: create and index dictionary

java -Xmx16g -jar ~/anaconda/pkgs/picard-1.106-0/jar/CreateSequenceDictionary.jar \
R=/Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
O=/Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.dict

samtools faidx /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta


# step 13: call indels

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 7_merge-bams/Megascops_choliba.bam   \
--minReadsAtLocus 7 \
-o 8_GATK/Megascops_choliba.intervals



# step 14: realign indels

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 7_merge-bams/Megascops_choliba.bam  \
-targetIntervals 8_GATK/Megascops_choliba.intervals \
-LOD 3.0 \
-o 8_GATK/Megascops_choliba_RI.bam


# step 15: call SNPs

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 8_GATK/Megascops_choliba_RI.bam \
-gt_mode DISCOVERY \
-o 8_GATK/Megascops_choliba_raw_SNPs.vcf \
-ploidy 2 \
-rf BadCigar


# step 16: annotate SNPs

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 8_GATK/Megascops_choliba_RI.bam \
-G StandardAnnotation \
-V:variant,VCF 8_GATK/Megascops_choliba_raw_SNPs.vcf \
-XA SnpEff \
-o 8_GATK/Megascops_choliba_SNPs_annotated.vcf \
-rf BadCigar      



# step 17: annotate indels

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 8_GATK/Megascops_choliba_RI.bam \
-gt_mode DISCOVERY \
-glm INDEL \
-o 8_GATK/Megascops_choliba_SNPs_indels.vcf \
-rf BadCigar         



# step 18: mask indels

java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-V 8_GATK/Megascops_choliba_raw_SNPs.vcf \
--mask 8_GATK/Megascops_choliba_SNPs_indels.vcf \
--maskExtension 5 \
--maskName InDel \
--clusterWindowSize 10 \
--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filterName "BadValidation" \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
--filterExpression "QD < 5.0" \
--filterName "LowVQCBD" \
-o 8_GATK/Megascops_choliba_SNPs_no_indels.vcf  \
-rf BadCigar


# step 19: Restrict to high-quality SNPs (bash)

cat 8_GATK/Megascops_choliba_SNPs_no_indels.vcf | grep 'PASS\|^#' > 8_GATK/Megascops_choliba_SNPs_pass-only.vcf 


# step 20: Read-backed phasing (GATK)


java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
-I 8_GATK/Megascops_choliba_RI.bam \
--variant 8_GATK/Megascops_choliba_SNPs_pass-only.vcf \
-L 8_GATK/Megascops_choliba_SNPs_pass-only.vcf \
-o 8_GATK/Megascops_choliba_SNPs_phased.vcf \
--phaseQualityThresh 20.0 \
-rf BadCigar



# steps 21-23

for SAMPLE in Megascops_choliba_AMNH4811 Megascops_choliba_FMNH392673 Megascops_choliba_LSUMNS7420 Megascops_choliba_LSUMNS9593 Megascops_choliba_LSUMNS15318 Megascops_choliba_LSUMNS38217 Megascops_choliba_LSUMNS42284 Megascops_choliba_MPEG6986 Megascops_choliba_MPEG13970 Megascops_choliba_MPEG16509 Megascops_choliba_MPEG16711 Megascops_choliba_MPEG17491 Megascops_choliba_UWBM116721 
do

	java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
	--variant 8_GATK/Megascops_choliba_SNPs_phased.vcf \
	-o 8_GATK/${SAMPLE}_SNPs.vcf \
	-sn ${SAMPLE} \
	-rf BadCigar


	java -Xmx16g -jar ~/anaconda/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
	-T VariantsToTable \
	-R /Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
	-V 8_GATK/${SAMPLE}_SNPs.vcf \
	-F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
	-o 9_SNP-tables/${SAMPLE}_SNPs_phased-table.txt \
	-rf BadCigar


	python /Users/mharvey/seqcap_pop-master/bin/add_phased_snps_to_seqs_filter.py \
	/Volumes/Brumfield_Lab_Drive/mike_data/mafft-nexus-edge-trimmed-exploded-uncontaminated/Megascops_choliba_MPEG16711.fasta \
	9_SNP-tables/${SAMPLE}_SNPs_phased-table.txt \
	10_sequences/${SAMPLE}_sequences.txt \
	1


done

: '

# step 24 

python /Users/mharvey/seqcap_pop-master/bin/collate_sample_fastas_GATK.py \
	10_sequences/ \
	11_fasta-parts/ \
	sequences.txt

# step 25: Align the sequences

python /Users/mharvey/seqcap_pop-master/bin/run_mafft.py \
	11_fasta-parts/ \
	12_raw-alignments/


# step 26: Process the alignments

python /Users/mharvey/seqcap_pop-master/bin/process_mafft_alignments_GATK.py \
	12_raw-alignments/ \
	13_processed-phylip/
	
# gphocs

python /Users/mharvey/seqcap_pop-master/bin/gphocs_from_phy.py \
	13_processed-phylip/ \
	14_formatted_output/gphocs/Megascops_choliba_GPHOCS.txt


# structure

python /Users/mharvey/seqcap_pop-master/bin/structure_from_vcf.py \
	8_GATK/Megascops_choliba_SNPs_phased.vcf \
	14_formatted_output/structure/Megascops_choliba_STRUCTURE.txt
	
	
# faststructure

python /Users/mharvey/seqcap_pop-master/bin/faststructure_from_vcf.py \
	8_GATK/Megascops_choliba_SNPs_phased.vcf \
	14_formatted_output/faststructure/ \
	Megascops_choliba \
	--random


# adegenet

python /Users/mharvey/seqcap_pop-master/bin/adegenet_from_vcf.py \
	8_GATK/Megascops_choliba_SNPs_phased.vcf \
	14_formatted_output/adegenet/ \
	Megascops_choliba \
	--random


# genepop

python /Users/mharvey/seqcap_pop-master/bin/genepop_from_vcf.py \
	8_GATK/Megascops_choliba_SNPs_phased.vcf \
	14_formatted_output/genepop/ \
	Megascops_choliba \
	--random


'

source activate phyluce


cd 8_GATK/

# remove indels, restrict to biallelic loci, min mean depth, min quality score
vcftools --vcf Megascops_choliba_SNPs_phased.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minDP 5.5 --minQ 30 --recode --out Megascops_choliba_SNPs_phased_rmIndels_QC_DP

# allow 25% missing data
vcftools --vcf Megascops_choliba_SNPs_phased_rmIndels_QC_DP.recode.vcf --max-missing 0.75 --recode --out Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP

# export mean depth per site
vcftools --vcf Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP.recode.vcf --site-mean-depth --out Megascops_choliba_SNPs_phased

# format files to get one snp per locus
cat Megascops_choliba_SNPs_phased.ldepth.mean | awk '{print $1,$2}' > Megascops_choliba_SNPs_phased_ok.txt

egrep -o "\w+.\w+\s" Megascops_choliba_SNPs_phased_ok.txt | sort -u | parallel 'egrep "{}" Megascops_choliba_SNPs_phased_ok.txt | sort -r | head -1 >> Megascops_choliba_SNPs_phased_random_snps.temp'


# format files
sort -u Megascops_choliba_SNPs_phased_random_snps.temp > Megascops_choliba_SNPs_phased_random_snps.txt

# extract one snp per locus
vcftools --vcf Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP.recode.vcf --positions Megascops_choliba_SNPs_phased_random_snps.txt --recode --out Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP_random

cp Megascops_choliba_SNPs_phased_rmIndels_75_QC_DP_random.recode.vcf /Volumes/Brumfield_Lab_Drive/data/1_analysis/vcf_files_final/

source deactivate

