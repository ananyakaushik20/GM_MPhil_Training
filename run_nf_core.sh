#!/bin/bash

#SBATCH -p long
#SBATCH --job-name=mouse_wt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak2470@cam.ac.uk 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=10:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 3.14.0 \
-resume \
-profile singularity \
--input design.csv \
--outdir /Users/ananya/Desktop/RNAseq/data \
--reads /Users/ananya/Desktop/RNAseq/data/*{_R1,_R2}.fastq.gz \
--fasta /Users/ananya/Desktop/RNAseq/data/GRCm38.p6.genome.fa \
--gtf /Users/ananya/Desktop/RNAseq/data/gencode.vM25.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--email ak2470@cam.ac.uk \
-c nextflow.config
