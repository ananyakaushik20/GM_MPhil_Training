#!/usr/bin/env nextflow

params.input_tar = "GM04_AdvBioinfo_NGS_Mar2020.tar"
params.reference = "Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa"
params.adapters = "primers_adapters.fa"
params.dbsnp = "latest_dbsnp.vcf"

process extractTar {
    input:
    path input_tar

    output:
    path "lane1", emit: lane1
    path "lane2", emit: lane2

    """
    tar -xvf ${input_tar}
    """
}

process runFastQC {
    input:
    path lane

    output:
    path "*_fastqc.zip"

    """
    fastqc ${lane}/s-7-1.fastq
    fastqc ${lane}/s-7-2.fastq
    """
}

process trimAdapters {
    input:
    path lane
    path adapters

    output:
    tuple path("${lane}/s-7-1.paired.fastq"), path("${lane}/s-7-2.paired.fastq")

    """
    trimmomatic PE -phred33 -threads 1 -trimlog ${lane}/trimm_logfile \
    ${lane}/s-7-1.fastq ${lane}/s-7-2.fastq \
    ${lane}/s-7-1.paired.fastq ${lane}/s-7-1.unpaired.fastq \
    ${lane}/s-7-2.paired.fastq ${lane}/s-7-2.unpaired.fastq \
    ILLUMINACLIP:${adapters}:2:30:10 MINLEN:36
    """
}

process trimLowQuality {
    input:
    tuple path(r1), path(r2)

    output:
    tuple path("${r1.simpleName}.trim.paired.fastq"), path("${r2.simpleName}.trim.paired.fastq")

    """
    trimmomatic PE -phred33 -threads 1 -trimlog trimm_logfile2 \
    ${r1} ${r2} \
    ${r1.simpleName}.trim.paired.fastq ${r1.simpleName}.trim.unpaired.fastq \
    ${r2.simpleName}.trim.paired.fastq ${r2.simpleName}.trim.unpaired.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process indexReference {
    input:
    path reference

    output:
    path "${reference}.*"

    """
    samtools faidx ${reference}
    bwa index -a is ${reference}
    picard CreateSequenceDictionary -R ${reference} -O ${reference.baseName}.dict
    """
}

process alignReads {
    input:
    tuple path(r1), path(r2)
    path reference
    val lane_id

    output:
    path "${lane_id}_sorted.bam"

    """
    bwa mem -R '@RG\\tID:${lane_id}\\tLB:library\\tPL:Illumina\\tPU:lane${lane_id}\\tSM:yeast' \
    ${reference} ${r1} ${r2} | \
    samtools view -b - | samtools sort - -o ${lane_id}_sorted.bam
    samtools index ${lane_id}_sorted.bam
    """
}

process mergeBams {
    input:
    path bams

    output:
    path "library.bam"

    """
    picard MergeSamFiles ${bams.collect { "INPUT=$it" }.join(' ')} OUTPUT=library.bam
    samtools index library.bam
    """
}

process markDuplicates {
    input:
    path bam

    output:
    path "library_final.bam"

    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=library_final.bam METRICS_FILE=dupl_metrics.txt
    samtools index library_final.bam
    """
}

process baseRecalibration {
    input:
    path bam
    path reference
    path dbsnp

    output:
    path "library_BQSR.bam"

    """
    gatk BaseRecalibrator -R ${reference} -I ${bam} --known-sites ${dbsnp} -O recal_data.table
    gatk ApplyBQSR -R ${reference} -I ${bam} --bqsr-recal-file recal_data.table -O library_BQSR.bam
    """
}

process callVariants {
    input:
    path bam
    path reference

    output:
    path "gatk_variants_raw_I.vcf"

    """
    gatk HaplotypeCaller -R ${reference} -I ${bam} -L I -O gatk_variants_raw_I.vcf
    """
}

workflow {
    extractTar(params.input_tar)
    runFastQC(extractTar.out.lane1.mix(extractTar.out.lane2))
    trimAdapters(extractTar.out.lane1.mix(extractTar.out.lane2), params.adapters)
    trimLowQuality(trimAdapters.out)
    indexReference(params.reference)
    alignReads(trimLowQuality.out, params.reference, Channel.from(1..2))
    mergeBams(alignReads.out.collect())
    markDuplicates(mergeBams.out)
    baseRecalibration(markDuplicates.out, params.reference, params.dbsnp)
    callVariants(baseRecalibration.out, params.reference)
}
