#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings/0_ref"
thirdgenome="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/whole_genome.third.fa"


for species in 0_human 2_dog; do 
    #echo $species 
    directoryname="$basedirectory/$species"
    queryfilename=$directoryname/*"-to-human_REF.bam"
    ls $queryfilename
    outfile="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/processing_comparisons/$species-to-human_REF-to-human_THIRD.bam"

    bwa bam2bam --only-aligned -g $thirdgenome -n 0.01 -o 2 -l 16500 -t 30 -f $outfile $queryfilename
done

exit
