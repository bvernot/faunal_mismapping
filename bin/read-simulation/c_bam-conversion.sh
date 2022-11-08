#!/bin/bash

fastaDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/2_mutated-seqs"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/3_bams"

for speciesDir in "$fastaDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        # input
        damageFasta=$speciesDir"/"$species"_mutated.fa.gz"

        # output
        mkdir -p $resultsDir"/"$withNumber
        outputBam=$resultsDir"/"$withNumber"/"$species"_simulated-and-mutated.bam"

        zcat $damageFasta | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > $outputBam &
    fi
done
wait
