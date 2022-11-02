#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
sizedistribution=$basedirectory"/size_distribution.tsv"

for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        fastafilename=$directoryname"/endo/"$species"_genome.fna"
        results=$directoryname/$species
        stdout=$directoryname"/stdout.log"
        stderr=$directoryname"/stderr.log"

        gargammel -f $sizedistribution -n 10000 -o $results -damage 0.03,0.4,0.01,0.3 $directoryname 1> $stdout 2> $stderr &
    fi
done

