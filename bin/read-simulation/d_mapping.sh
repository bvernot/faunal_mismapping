#!/bin/bash

bamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/3_bams"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings"

genomeRef="/mnt/solexa/Genomes/hg19_evan/bwa-0.4.9"
genomeThird="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/whole_genome.third.fa"
genomeN="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/whole_genome.N.fa"

for speciesDir in "$bamDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        # input
        queryBam=$bamDir"/"$withNumber"/"$species"_simulated-and-mutated.bam"

        ## REF
        jobName=$species"-to-human_REF"
        mapDir=$resultsDir"/0_ref/"$withNumber
        mkdir -p $mapDir
        outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeRef -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned


        ## THIRD
        jobName=$species"-to-human_THIRD"
        mapDir=$resultsDir"/1_third/"$withNumber
        mkdir -p $mapDir
        outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeThird -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned


        ## N
        jobName=$species"-to-human_N"
        mapDir=$resultsDir"/2_N/"$withNumber
        mkdir -p $mapDir
        outputBam=$mapDir"/"$jobName".bam"

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $genomeN -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned
    fi
done
wait
