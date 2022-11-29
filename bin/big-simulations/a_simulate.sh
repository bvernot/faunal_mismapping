#!/bin/bash
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/1_simulations"
sizeDistribution="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/size_distributions.tsv"

# genome directory
dogGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/2_dog"
bearGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/3_bear"
sheepGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/0_genomes/0_sheep"

for genomeDir in $dogGenomeDir $bearGenomeDir $sheepGenomeDir; do
    species=`echo $genomeDir | awk -F "_" '{print $NF}'`
    echo $species

    speciesDir=$resultsDir"/"$species
    mkdir -p $speciesDir

    for i in {1..20}; do
        finalDir=$speciesDir"/"$i
        mkdir -p $finalDir

        stdout=$finalDir"/stdout.log"
        stderr=$finalDir"/stderr.log"

        gargammel -f $sizeDistribution -n 30000000 -o $finalDir"/"species -damage 0.03,0.4,0.01,0.3 $genomeDir 1> $stdout 2> $stderr &

    done
    wait
done