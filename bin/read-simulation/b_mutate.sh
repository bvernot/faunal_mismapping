#!/bin/bash

simuDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/1_simulations"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/2_mutated-seqs"

for speciesDir in "$simuDir"/*; do 
    if [ -d "$speciesDir" ]; then
        species=${speciesDir##*_}
        withNumber=${speciesDir##*/}

        mkdir -p $resultsDir"/"$withNumber

        # input
        inputFasta=$simuDir"/"$withNumber"/"$species"_d.fa.gz"
        
        # output
        outputFasta=$resultsDir"/"$withNumber"/"$species"_mutated.fa.gz"

        python2 "x-b_mutate.py" $inputFasta 0.002 | gzip -c > $outputFasta &
    fi
done
wait
