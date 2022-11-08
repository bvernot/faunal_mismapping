#!/bin/bash

bamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings"

for speciesDir in "$bamDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        # input
        queryBam=$bamDir"/"$withNumber"/"$species"_simulated.bam"

        ## REF
        jobName=$species"-to-human_REF"
        $mapDir=$resultsDir"/0_ref/"$withNumber
        mkdir -p $mapDir
        $outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeRef -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned


        ## THIRD
        jobName=$species"-to-human_THIRD"
        $mapDir=$resultsDir"/1_third/"$withNumber
        mkdir -p $mapDir
        $outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeThird -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned


        ## N
        jobName=$species"-to-human_N"
        $mapDir=$resultsDir"/2_N/"$withNumber
        mkdir -p $mapDir
        $outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $gneomeN -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned
    fi
done
wait
