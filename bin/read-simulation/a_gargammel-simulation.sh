#!/bin/bash

dataDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/1_simulations"
sizeDistribution=$dataDir"/size_distribution.tsv"

for speciesDir in "$dataDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        # input
        speciesGenome=$speciesDir"/endo/"$species"_genome.fna"

        # output
        outputDir=$resultsDir"/"$withNumber"/"$species
        mkdir -p $resultsDir"/"$withNumber
        stdout=$resultsDir"/"$withNumber"/stdout.log"
        stderr=$resultsDir"/"$withNumber"/stderr.log"

        gargammel -f $sizeDistribution -n 10000000 -o $outputDir -damage 0.03,0.4,0.01,0.3 $speciesDir 1> $stdout 2> $stderr &
    fi
done
wait
