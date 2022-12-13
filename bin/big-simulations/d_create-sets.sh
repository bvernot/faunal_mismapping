#!/bin/bash
dataDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/2_mutated"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/3_merged-fastas"

for speciesDir in $dataDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`

    speciesResults=$resultsDir"/"$species
    mkdir -p $speciesResults

    firstResult=$speciesResults"/1st_"$species".fna.gz"
    secondResult=$speciesResults"/2st_"$species".fna.gz"
    thirdResult=$speciesResults"/3st_"$species".fna.gz"

    zcat $speciesDir"/"* | head -n 200000000 | gzip -c > $firstResult &
    zcat $speciesDir"/"* | head -n 400000000 | tail -n 200000000 | gzip -c > $secondResult &
    zcat $speciesDir"/"* | tail -n 200000000 | gzip -c > $thirdResult &

done
wait