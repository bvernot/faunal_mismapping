#!/bin/bash
genomeRef="/mnt/solexa/Genomes/hg19_evan/bwa-0.4.9"

bamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/4_bams"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/5_mapped"

for speciesDir in $bamDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    speciesResults=$resultsDir"/"$species
    mkdir -p $speciesResults

    for queryBam in $speciesDir"/"*; do
        name=`echo $queryBam | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}'`
        jobName=`echo $queryBam | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'`"-to-human"

        finalDir=$speciesResults"/"$name
        mkdir -p $finalDir

        outputBam=$finalDir"/"$jobName".bam"

        /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeRef -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned

    done
done