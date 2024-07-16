#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ExtractTar {
    input:
    path tar_file from 'GM04_AdvBioinfo_NGS_Mar2020.tar'

    output:
    path 'GM04_AdvBioinfo_NGS_Mar2020' into extracted_files

    script:
    """
    tar -xvf ${tar_file}
    """
}

process FastQC {
    input:
    path fastq_files from extracted_files.collect()
    
    output:
    path '*_fastqc.zip', path '*_fastqc.html'

    script:
    """
    fastqc ${fastq_files.join(' ')}
    """
}

process TrimAdapters {
    input:
    path fastq1 from extracted_files.filter { it.name.endsWith('s-7-1.fastq') }
    path fastq2 from extracted_files.filter { it.name.endsWith('s-7-2.fastq') }
    path adapters from 'primers_adapters.fa'

    output:
    path 'lane*paired.fastq', path 'lane*unpaired.fastq'

    script:
    """
    trimmomatic PE -phred33 -threads 1 -trimlog trimm_logfile \\
        ${fastq1} ${fastq2} \\
        ${fastq1.baseName}.paired.fastq ${fastq1.baseName}.unpaired.fastq \\
        ${fastq2.baseName}.paired.fastq ${fastq2.baseName}.unpaired.fastq \\
        ILLUMINACLIP:${adapters}:2:30:10 MINLEN:36
    """
}

process TrimLowQuality {
    input:
    path fastq1 from TrimAdapters.out.filter { it.name.endsWith('1.paired.fastq') }
    path fastq2 from TrimAdapters.out.filter { it.name.endsWith('2.paired.fastq') }

    output:
    path 'lane*trim.paired.fastq', path 'lane*trim.unpaired.fastq'

    script:
    """
    trimmomatic PE -phred33 -threads 1 -trimlog trimm_logfile2 \\
        ${fastq1} ${fastq2} \\
        ${fastq1.baseName}.trim.paired.fastq ${fastq1.baseName}.trim.unpaired.fastq \\
        ${fastq2.baseName}.trim.paired.fastq ${fastq2.baseName}.trim.unpaired.fastq \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process IndexReference {
    input:
    path reference from 'Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa'

    output:
    path 'Saccharomyces_cerevisiae.EF4.68.dna.toplevel.*'

    script:
    """
    samtools faidx ${reference}
    bwa index -a is ${reference}
    picard CreateSequenceDictionary \\
        -R ${reference} \\
        -O ${reference.baseName}.dict
    """
}

process Align {
    input:
    path reference from IndexReference.out
    path fastq1 from TrimLowQuality.out.filter { it.name.contains('1.trim.paired.fastq') }
    path fastq2 from TrimLowQuality.out.filter { it.name.contains('2.trim.paired.fastq') }
    each lane_id, lane_name from Channel.of(['1', 'lane1'], ['2', 'lane2'])

    output:
    path "lane${lane_id}_sorted.bam"

    script:
    """
    bwa mem -R '@RG\\tID:${lane_id}\\tLB:library\\tPL:Illumina\\tPU:${lane_name}\\tSM:yeast' \\
        ${reference.baseName}.fa ${fastq1} ${fastq2} | \\
        samtools view -b - | \\
        samtools sort - -o lane${lane_id}_sorted.bam && \\
    samtools index lane${lane_id}_sorted.bam
    """
}

process MergeBam {
    input:
    path bam_files from Align.out.collect()

    output:
    path 'library.bam'

    script:
    """
    picard MergeSamFiles \\
        ${bam_files.collect { "INPUT=" + it }.join(' ')} \\
        OUTPUT=library.bam
    samtools index library.bam
    """
}

process MarkDuplicates {
    input:
    path bam from MergeBam.out

    output:
    path 'library_final.bam'

    script:
    """
    picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=library_final.bam \\
        METRICS_FILE=dupl_metrics.txt
    samtools index library_final.bam
    """
}

process BamQC {
    input:
    path bam from MarkDuplicates.out

    output:
    path 'qualimap_report', path 'mito_coverage'

    script:
    """
    qualimap bamqc -bam ${bam} -outdir qualimap_report
    samtools depth -a -r Mito:1-85779 ${bam} > mito_coverage
    """
}

process MitoBam {
    input:
    path bam from MarkDuplicates.out

    output:
    path 'mito.bam'

    script:
    """
    samtools view -bh -o mito.bam ${bam} Mito
    samtools index mito.bam
    """
}

process BaseRecalibration {
    input:
    path reference from 'Yeast.fasta'
    path bam from MarkDuplicates.out
    path known_sites from 'latest_dbsnp.vcf'

    output:
    path 'recal_data.table', path 'library_BQSR.bam'

    script:
    """
    gatk BaseRecalibrator -R ${reference} -I ${bam} -knownSites ${known_sites} -O recal_data.table
    gatk ApplyBQSR -R ${reference} -I ${bam} --bqsr-recal-file recal_data.table -O library_BQSR.bam
    """
}

process HaplotypeCaller {
    input:
    path reference from 'Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa'
    path bam from BaseRecalibration.out.filter { it.name == 'library_BQSR.bam' }

    output:
    path 'gatk_variants_raw_I.vcf', path 'gatk_variants_raw_I_bq20_mq50.vcf'

    script:
    """
    gatk HaplotypeCaller -R ${reference} -I ${bam} -L I -O gatk_variants_raw_I.vcf
    gatk HaplotypeCaller -R ${reference} -I ${bam} -L I -mbq 20 --minimum-mapping-quality 50 -O gatk_variants_raw_I_bq20_mq50.vcf
    """
}

process SelectVariants {
    input:
    path reference from 'Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa'
    path vcf from HaplotypeCaller.out.filter { it.name == 'gatk_variants_raw_I.vcf' }

    output:
    path 'gatk_variants_raw_I_SNP.vcf'

    script:
    """
    gatk SelectVariants -R ${reference} --variant ${vcf} -O gatk_variants_raw_I_SNP.vcf --select-type SNP
    """
}

process FreeBayes {
    input:
    path reference from 'Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa'
    path bam from BaseRecalibration.out.filter { it.name == 'library_BQSR.bam' }

    output:
    path 'freebayes_variants_raw_I_bq20_mq50.vcf'

    script:
    """
    freebayes -q 20 -m 50 -u -f ${reference} ${bam} -r I > freebayes_variants_raw_I_bq20_mq50.vcf
    """
}

process CompareVariants {
    input:
    path vcf1 from SelectVariants.out
    path vcf2 from FreeBayes.out

    output:
    path 'compare.diff'

    script:
    """
    vcftools --vcf ${vcf1} --diff ${vcf2} --diff-site --out compare
    """
}

process FilterVariants {
    input:
    path vcf from SelectVariants.out

    output:
    path 'gatk_variants_I_flt2.recode.vcf', path 'gatk_variants_I_flt2_nomissing.recode.vcf'

    script:
    """
    vcftools --vcf ${vcf} --minDP 3 --minQ 20 --out gatk_variants_I_flt2 --recode --recode-INFO-all
    vcftools --vcf gatk_variants_I_flt2.recode.vcf --max-missing 1 --out gatk_variants_I_flt2_nomissing --recode --recode-INFO-all
    """
}

process QualFilterVariants {
    input:
    path vcf from FilterVariants.out.filter { it.name == 'gatk_variants_I_flt2_nomissing.recode.vcf' }

    output:
    path 'gatk_QUAL_100_variants'

    script:
    """
    grep -v "#" ${vcf} | awk '{if(\$6>100) print}' > gatk_QUAL_100_variants
    grep -v "#" gatk_QUAL_100_variants | wc -l
    """
}

workflow {
    tar_file = 'GM04_AdvBioinfo_NGS_Mar2020.tar'

    extracted_files = ExtractTar(tar_file)
    fastqc_files = FastQC(extracted_files)
    trimmed_adapters = TrimAdapters(extracted_files)
    trimmed_low_quality = TrimLowQuality(trimmed_adapters)
    indexed_reference = IndexReference('Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa')

    alignment = Align(indexed_reference, trimmed_low_quality, ['1', 'lane1'], ['2', 'lane2'])
    merged_bam = MergeBam(alignment)
    marked_duplicates = MarkDuplicates(merged_bam)
    bam_qc = BamQC(marked_duplicates)
    mito_bam = MitoBam(marked_duplicates)
    base_recalibration = BaseRecalibration('Yeast.fasta', marked_duplicates, 'latest_dbsnp.vcf')
    haplotype_calling = HaplotypeCaller('Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa', base_recalibration)

    selected_variants = SelectVariants('Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa', haplotype_calling)
    freebayes_variants = FreeBayes('Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa', base_recalibration)

    compare_variants = CompareVariants(selected_variants, freebayes_variants)
    filtered_variants = FilterVariants(selected_variants)
    qual_filtered_variants = QualFilterVariants(filtered_variants)
}
