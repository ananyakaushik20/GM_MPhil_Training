#! bin/bash

# download reference files
wget -P /Users/ananya/WGS/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /Users/ananya/WGS/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx /Users/ananya/WGS/supporting_files/hg38/hg38.fa

# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/Users/ananya/WGS/supporting_files/hg38/hg38.fa O=~/Desktop/demo/supporting_files/hg38/hg38.dict

# download known sites files for BQSR from GATK resource bundle
wget -P /Users/ananya/WGS/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /Users/ananya/WGS/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# directories
ref="/Users/ananya/WGS/supporting_files/hg38.fa"
known_sites="/Users/ananya/WGS/supporting_files/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/Users/ananya/WGS/aligned_reads"
reads="/Users/ananya/WGS/reads"
results="/Users/ananya/WGS/results"
data="/Users/ananya/WGS/data"

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
#bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
