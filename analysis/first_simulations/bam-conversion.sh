#!/bin/bash

zcat 0_human/human_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 0_human/human_from_fa.bam &
zcat 1_macaque/macaque_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 1_macaque/macaque_from_fa.bam &
zcat 2_dog/dog_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 2_dog/dog_from_fa.bam &
zcat 3_bear/bear_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 3_bear/bear_from_fa.bam &
zcat 4_bison/bison_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 4_bison/bison_from_fa.bam &
zcat 5_mouse/mouse_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 5_mouse/mouse_from_fa.bam &
zcat 6_chicken/chicken_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 6_chicken/chicken_from_fa.bam &
zcat 7_maize/maize_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 7_maize/maize_from_fa.bam &
zcat 8_drosophila/drosophila_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 8_drosophila/drosophila_from_fa.bam &
zcat 9_yeast/yeast_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > 9_yeast/yeast_from_fa.bam &
zcat X_anthrax/anthrax_d.fa.gz | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > X_anthrax/anthrax_from_fa.bam &
