#!/bin/bash
unmergedDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/5_mapped"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/6_merged-bams"

for speciesDir in $unmergedDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    mkdir -p $resultsDir"/"$species

    ## first ind
    f1=$speciesDir"/1st/1st_"$species"-to-human.bam"
    f2=$speciesDir"/2st/2st_"$species"-to-human.bam"
    samtools merge -f $resultsDir"/"$species"/"$species"_1.bam" $f1 $f2 &

    ## second ind
    f3=$speciesDir"/3st/3st_"$species"-to-human.bam"
    f4=$speciesDir"/4th/4th_"$species"-to-human.bam"
    samtools merge -f $resultsDir"/"$species"/"$species"_2.bam" $f3 $f4 &

    ## third ind
    f5=$speciesDir"/5th/5th_"$species"-to-human.bam"
    f6=$speciesDir"/6th/6th_"$species"-to-human.bam"
    samtools merge -f $resultsDir"/"$species"/"$species"_3.bam" $f5 $f6 &

done
wait