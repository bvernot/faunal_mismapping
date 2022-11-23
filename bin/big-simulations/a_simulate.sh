#!/bin/bash

resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/1_simulations"
sizeDistribution="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/size_distributions.tsv"

# genome directory
dogGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/2_dog"
bearGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/0_genomes/3_bear"
sheepGenomeDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/0_genomes/0_sheep"

# gargammel usage
echo gargammel -f $sizeDistribution -n 10000000 -o $outputDir -damage 0.03,0.4,0.01,0.3 $speciesDir 1> $stdout 2> $stderr

for genomeDir in $dogGenomeDir $bearGenomeDir $sheepGenomeDir; do
    species=`echo $genomeDir | awk -F "_" '{print $NF}'`
    echo $genomeDir $species


done