#!/bin/bash

bamBaseDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/5_calmd-and-indexed"
resultsBaseDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/6_parsed-bams"

for speciesDir in "$bamBaseDir"/0_ref/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        echo $species

        ## REF ##
        bamFile=$bamBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF-calmd.bam"
        tsvFile=$resultsBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF.tsv"
        mkdir -p $resultsBaseDir"/0_ref/"$withNumber
        
        python3 x-g_bam-parser.py $bamFile $tsvFile &


        ## THIRD ##
        bamFile=$bamBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD-calmd.bam"
        tsvFile=$resultsBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD.tsv"
        mkdir -p $resultsBaseDir"/1_third/"$withNumber

        python3 x-g_bam-parser.py $bamFile $tsvFile &
        

        ## N ##
        bamFile=$bamBaseDir"/2_N/"$withNumber"/"$species"-to-human_N-calmd.bam"
        tsvFile=$resultsBaseDir"/2_N/"$withNumber"/"$species"-to-human_N.tsv"
        mkdir -p $resultsBaseDir"/2_N/"$withNumber

        python3 x-g_bam-parser.py $bamFile $tsvFile &
    fi
done
wait
