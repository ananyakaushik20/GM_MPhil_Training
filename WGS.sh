#! bin/bash

 tar -xvf GM04_AdvBioinfo_NGS_Mar2020.tar

 fastqc lane1/s-7-1.fastq
 fastqc lane1/s-7-2.fastq

 fastqc lane2/s-7-1.fastq
 fastqc lane2/s-7-2.fastq

##trimming adapter sequences

#lane 1
trimmomatic PE -phred33 -threads 1 -trimlog \ 
lane1/trimm_logfile lane1/s-7-1.fastq lane1/s-7-2.fastq \ 
lane1/s-7-1.paired.fastq lane1/s-7-1.unpaired.fastq \ 
lane1/s-7-2.paired.fastq lane1/s-7-2.unpaired.fastq \ 
ILLUMINACLIP:primers_adapters.fa:2:30:10 MINLEN:36

#lane2
trimmomatic PE -phred33 -threads 1 -trimlog \ 
lane2/trimm_logfile lane2/s-7-1.fastq lane2/s-7-2.fastq \ 
lane2/s-7-1.paired.fastq lane2/s-7-1.unpaired.fastq \ 
lane2/s-7-2.paired.fastq lane2/s-7-2.unpaired.fastq \ 
ILLUMINACLIP:primers_adapters.fa:2:30:10 MINLEN:36

##trimming low quality bases

#lane 1
trimmomatic PE -phred33 -threads 1 -trimlog \ 
lane1/trimm_logfile2 \
lane1/s-7-1.paired.fastq lane1/s-7-2.paired.fastq \ 
lane1/s-7-1.trim.paired.fastq lane1/s-7-1.trim.unpaired.fastq \ 
lane1/s-7-2.trim.paired.fastq lane1/s-7-2.trim.unpaired.fastq \ 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#lane 2
trimmomatic PE -phred33 -threads 1 -trimlog \ 
lane2/trimm_logfile2 \
lane2/s-7-1.paired.fastq lane2/s-7-2.paired.fastq \ 
lane2/s-7-1.trim.paired.fastq lane2/s-7-1.trim.unpaired.fastq \ 
lane2/s-7-2.trim.paired.fastq lane2/s-7-2.trim.unpaired.fastq \ 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## indexing

samtools faidx Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa
bwa index -a is Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa
picard CreateSequenceDictionary \ -R Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa \ -O Saccharomyces_cerevisiae.EF4.68.dna.toplevel.dict

#alignment
bwamem-R\ '@RG\tID:1\tLB:library\tPL:Illumina\tPU:lane1\tSM:yeast' \ Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa \ lane1/s-7-1.trim.paired.fastq lane1/s-7-2.trim.paired.fastq | \ samtools view -b - | samtools sort - -o lane1_sorted.bam &&

samtools index lane1_sorted.bam

bwamem-R\ '@RG\tID:2\tLB:library\tPL:Illumina\tPU:lane2\tSM:yeast' \ Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa \ lane2/s-7-1.trim.paired.fastq lane2/s-7-2.trim.paired.fastq | \ samtools view -b - | samtools sort - -o lane2_sorted.bam &&

samtools index lane2_sorted.bam

#check BAM files
samtools flagstat lane1_sorted.bam
samtools flagstat lane2_sorted.bam

#merge BAM files
picard MergeSamFiles INPUT=lane1_sorted.bam INPUT=lane2_sorted.bam OUTPUT=library.bam

samtools index library.bam

### REFINEMENT ###
#Remove Duplicates
picard MarkDuplicates INPUT=library.bam OUTPUT=library_final.bam METRICS_FILE=dupl_metrics.txt
samtools index library_final.bam

##BAM QC

qualimap bamqc -bam library_final.bam -outdir qualimap_report
samtools depth -a -r Mito:1-85779 library_final.bam > mito_coverage

##BAM visualisation

samtools view -bh -o mito.bam library_final.bam Mito
samtools index mito.bam

##GATK BASE RECALIBRATION
gatk BaseRecalibrator -R Yeast.fasta -I library_final.bam - knownSites latest_dbsnp.vcf -o recal_data.table
gatk ApplyBQSR -R Yeast.fa -I library_final.bam --bqsr-recal- file recal_data.table -o library_BQSR.bam

### GATK Haplotype Caller ###

#Calling raw variants with no filters on chromosome 1
gatk HaplotypeCaller -R Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa -I library_final.bam -L I -O gatk_variants_raw_I.vcf

#Calling raw variants with no filters on chromosome 1 using a filter of minimum base quality 20 and minimum mapping quality 50
gatk HaplotypeCaller -R Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa -I library_final.bam -L I -mbq 20 --minimum-mapping-quality 50 -O gatk_variants_raw_I_bq20_mq50.vcf

##Variant selection
gatk SelectVariants -R Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa --variant gatk_variants_raw_I.vcf -O gatk_variants_raw_I_SNP.vcf --select- type SNP

##FreeBayes
#Call raw variants using a filter of minimum base quality 20 and minimum mapping quality 50.
freebayes -q 20 -m 50 -u -f Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa library_final.bam -r I > freebayes_variants_raw_I_bq20_mq50.vcf

gatk SelectVariants -R Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa --variant freebayes_variants_raw_I_bq20_mq50.vcf -O freebayes_variants_raw_I_bq20_mq50_SNP.vcf
--select-type SNP

##comparison of variant calling
vcftools --vcf gatk_variants_raw_I_bq20_mq50_SNP.vcf --diff freebayes_variants_raw_I_bq20_mq50_SNP.vcf --diff-site --out compare


###Variant Filtering###
vcftools --vcf gatk_variants_raw_I_bq20_mq50_SNP.vcf --minDP 3 - -minQ 20 --out gatk_variants_I_flt2 --recode --recode-INFO-all
vcftools --vcf gatk_variants_I_flt2.recode.vcf --max-missing 1 - -out gatk_variants_I_flt2_nomissing --recode --recode-INFO-all

grep -v “#” <your vcf file> | awk '{if($6>100) print}' > gatk_QUAL_100_variants
grep -v ”#” gatk_QUAL_100_variants | wc -l

