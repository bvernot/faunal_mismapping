#!/bin/bash

bamBaseDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/5_calmd-and-indexed"

for speciesDir in "$bamBaseDir"/0_ref/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        echo $species

        ## REF ##
        bamFile=$bamBaseDir"/0_ref/"$withNumber"/"$species"-to-human_REF-calmd.bam"
        samtools index $bamFile

        ## THIRD ##
        bamFile=$bamBaseDir"/1_third/"$withNumber"/"$species"-to-human_THIRD-calmd.bam"
        samtools index $bamFile
        
        ## N ##
        bamFile=$bamBaseDir"/2_N/"$withNumber"/"$species"-to-human_N-calmd.bam"
        samtools index $bamFile
    fi
done
wait
