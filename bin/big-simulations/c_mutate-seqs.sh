#!/bin/bash
dataDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/1_simulations"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/2_mutated"

for speciesDir in $dataDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    mutatedSpeciesDir=$resultsDir"/"$species
    mkdir -p $mutatedSpeciesDir

    for numberDir in $speciesDir"/"*; do
        origFile=$numberDir"/species_d.fa.gz"
        echo "#" $origFile

        number=`echo $origFile | cut -d "/" -f 10`

        mutatedFile=$mutatedSpeciesDir"/"$species"-"$number"_mutated.fna.gz"

        python2 "/mnt/expressions/robin_warner/3_faunal-mismapping/bin/read-simulation/x-b_mutate.py" $origFile 0.002 | gzip -c > $mutatedFile &
        
    done
    wait
done