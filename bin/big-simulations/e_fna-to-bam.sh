#!/bin/bash
fastaDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/3_merged-fastas"
bamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/4_bams"

for speciesDir in $fastaDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    speciesBamDir=$bamDir"/"$species
    mkdir -p $speciesBamDir

    for numberDir in $speciesDir"/"*; do
        name=`echo $numberDir | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'`
        resultsFilename=$speciesBamDir"/"$name".bam"

        zcat $numberDir | fasta_to_fastq | sed 's/\^/]/g' | fastqtobam qualitymax=70 > $resultsFilename &
        
    done
    wait
done