#!/bin/bash
dataDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/1_simulations"

for speciesDir in $dataDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    for numberDir in $speciesDir"/"*; do
        echo "#" $numberDir

        for unnededFile in species_a.fa.gz species.b.fa.gz species.c.fa.gz species.e.fa.gz species_s1.fq.gz species_s2.fq.gz; do
            fullFilePath=$numberDir"/"$unnededFile

            rm $fullFilePath
        
        done
    done
done