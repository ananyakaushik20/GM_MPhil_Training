rule all:
    input:
        "reference_genome.fasta.bwt",
        "fastqc_out/reads_fastqc.html",
        "aligned_sorted.bam.bai",
        "variants.vcf"


rule fastqc:
    input:
        fastq="reads.fastq"
    output:
        "fastqc_out/reads_fastqc.html",
        "fastqc_out/reads_fastqc.zip"
    shell:
        "fastqc -o fastqc_out {input.fastq}"

rule bwa_index:
    input:
        ref="reference_genome.fasta"
    output:
        "reference_genome.fasta.bwt"
    shell:
        "bwa index {input.ref}"
        
rule bwa_mem:
    input:
        ref="reference_genome.fasta",
        fastq="reads.fastq"
    output:
        "aligned.sam"
    shell:
        "bwa mem {input.ref} {input.fastq} > {output}"

rule samtools_view:
    input:
        "aligned.sam"
    output:
        "aligned.bam"
    shell:
        "samtools view -bS {input} > {output}"

rule samtools_sort:
    input:
        "aligned.bam"
    output:
        "aligned_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index:
    input:
        "aligned_sorted.bam"
    output:
        "aligned_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_mpileup:
    input:
        ref="reference_genome.fasta",
        bam="aligned_sorted.bam"
    output:
        "variants.vcf"
    shell:
        "samtools mpileup -f {input.ref} {input.bam} > {output}"
