rule all: 
    input:
        "lane1/s-7-1_fastqc.html", 
        "lane1/s-7-2_fastqc.html",
        "lane2/s-7-1_fastqc.html", 
        "lane2/s-7-2_fastqc.html",
        "lane1/s-7-1.paired.fastq", 
        "lane1/s-7-1.unpaired.fastq", 
        "lane1/s-7-2.paired.fastq", 
        "lane1/s-7-2.unpaired.fastq",
        "lane2/s-7-1.paired.fastq", 
        "lane2/s-7-1.unpaired.fastq", 
        "lane2/s-7-2.paired.fastq", 
        "lane2/s-7-2.unpaired.fastq",

rule fastqc_a_file: 
    input:
        "{filename}.fastq" 
    output:
        "{filename}_fastqc.html", 
        "{filename}_fastqc.zip"
    shell:
        "fastqc {input}"

rule trim_reads: 
    input:
        "{lane}/{filename}-1.fastq",
        "{lane}/{filename}-2.fastq"
    output: 
        "{lane}/{filename}-1.paired.fastq",
        "{lane}/{filename}-1.unpaired.fastq",
        "{lane}/{filename}-2.paired.fastq", 
        "{lane}/{filename}-2.unpaired.fastq"
    shell:
        "trimmomatic PE {input} {output} LEADING:2 TRAILING:2\
        SLIDINGWINDOW:4:15\
        MINLEN:25"

rule bwa_mem: 
    input:
        "Yeast.fa", 
        "lane1/{filename}-1.fastq", 
        "lane1/{filename}-2.fastq"
    output: 
        "lane1.sam"
    shell:
        "bwa mem {input} > {output}"