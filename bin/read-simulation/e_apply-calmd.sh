#!/bin/bash

bamBaseDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings"
resultsBaseDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/5_calmd-and-indexed"

genomeDir="/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"

for speciesDir in "$bamBaseDir"/0_ref/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        ## REF ##
        bamFile=$bamBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF.bam"
        bamResult=$resultsBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF-calmd.bam"
        log=$resultsBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF-calmd.log"
        mkdir -p $resultsBaseDir"/0_ref/"$withNumber
        
        samtools calmd -b --reference $genomeDir $bamFile > $bamResult 2> $log &


        ## THIRD ##
        bamFile=$bamBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD.bam"
        bamResult=$resultsBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD-calmd.bam"
        log=$resultsBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD-calmd.log"
        mkdir -p $resultsBaseDir"/1_third/"$withNumber

        samtools calmd -b --reference $genomeDir $bamFile 1> $bamResult 2> $log &


        ## N ##
        bamFile=$bamBaseDir"/2_N/"$withNumber"/"$species"-to-human_N.bam"
        bamResult=$resultsBaseDir"/2_N/"$withNumber"/"$species"-to-human_N-calmd.bam"
        log=$resultsBaseDir"/2_N/"$withNumber"/"$species"-to-human_N-calmd.log"
        mkdir -p $resultsBaseDir"/2_N/"$withNumber

        samtools calmd -b --reference $genomeDir $bamFile 1> $bamResult 2> $log &
    fi
done
wait
