#!/bin/bash

bamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/3_bams"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings"

for speciesDir in "$bamDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        # input
        queryBam=$bamDir"/"$withNumber"/"$species"_simulated.bam"

        ## REF


        ## THIRD


        ## N

        # /mnt/solexa/bin/mappr-cli -z $jobName -a -g $gneomeN -f $outputBam $queryBam -n 0.01 -o 2 -l 16500 --only-aligned
    fi
done
wait
